################################################################################
##
## This library is free software; you can redistribute it and/or
## modify it under the terms of the GNU Lesser General Public
## License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## 
## This library is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## Lesser General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public
## License along with this library; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
##
## (C) Copyrights Dr. Michel F. Sanner and TSRI 2016
##
################################################################################

#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2016
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/PmvApp/cartoonCmds.py,v 1.1.4.2 2017/09/25 21:16:20 annao Exp $
# 
# $Id: cartoonCmds.py,v 1.1.4.2 2017/09/25 21:16:20 annao Exp $
#
import numpy

from math import sqrt

from DejaVu2.Geom import Geom
from DejaVu2.Cylinders import Cylinders
from DejaVu2.Spheres import Spheres
from DejaVu2.IndexedPolygons import IndexedPolygons
from DejaVu2.Shapes import Rectangle2D, Circle2D, Ellipse2D

from MolKit2.selection import Selection

from AppFramework.App import DeleteObjectEvent

from PmvApp.displayCmds import DisplayCommand
from PmvApp.Pmv import MVCommand  #MVAtomICOM
from PmvApp.splineFrames import SplineFrames
from PmvApp.extrude import Extruder
from PmvApp.cubicInterpolate import ResampleControlPoints, GetSmoothNormals, GetFrames

from MolKit2.selection import SelectionSet

class ComputeCartoon(MVCommand):
    """The ComputeCartoon command computes a graphical representation of molecule's backbones. This represention relies on a 2D path representing a plane following a spline through Calpha positions of amino acids and P atoms for nucleic acids. Different shapes can be extruded along this path for varous secondary structure elements. The quality controls the number of points used for each amino or nucleic acid to discretize the spline.The calculation is carried out separately for each chain in the molecule.

    Synopsis:
        None <- ComputeCartoon(atoms, helixShape=None, sheetShape=None, coilShape=None, shapeNABB=None, quality=3)

    Required Arguments :
        a set of atoms --- the command will compute the cartoon for all molecules with an atoms in the set.

    Optional Arguments:

        quality --- the number of points used to represent a single aminio or nucleic acid backbone

        helixShape --- DejaVu2.Shape 2D object extruded for helical segments.
                       A value of 'None' defaults to Ellipse2D(0.6, 0.1, quality=quality)

        sheetShape --- DejaVu2.Shape 2D object extruded for beta-sheets segments.
                       A value of 'None' defaults to shapeSheet = Rectangle2D(width=1.2, height=0.2)

        coilShape  --- DejaVu2.Shape 2D object extruded for segements that are neither helical
                       nor beta-sheets.
                       A value of 'None' defaults to Circle2D(radius=0.1, quality=quality)
        shapeNABB  --- DejaVu2.Shape 2D object extruded for nucleic acids backbone
                       A value of 'None' defaults to Circle2D(radius=0.3, quality=quality)

    Package : PmvApp
    Module  : cartoonCmds
    Class   : ComputeCartoonPath
    Command name : computeCartoonPath
    """
    _argNames = ['helixShape', 'sheetShape', 'coilShape', 'shapeNABB', 'quality', 'scale']

    def checkArguments(self, nodes, **kw):
        for name in kw.keys():
            if name not in self._argNames:
                raise RuntimeError("%s: unrecognized keyword argument '%s', valid names are %s"%(
                    self.name, name, str(self._argNames)))
        return (nodes,), kw

    def onAddCmdToApp(self):
        self._extruder = Extruder()
        app = self.app()
        app.eventHandler.registerListener(DeleteObjectEvent,
                                          self.handleDeleteObject)
        self._splineFrames = SplineFrames()
        
    def handleDeleteObject(self, event):
        mol = event.object 
        initialized = self.initializedFor.get(mol, False)
        if initialized:
            del mol._cartoonGeoms

    def isInitialized(self, mol):
        return self.initializedFor.get(mol, False)
    
    def initialize(self,mol):
        # call self.initializedFor is needed
        isInitialized = self.initializedFor.get(mol, False)
        if not isInitialized:
            self.initializeMolForCmd(mol) 

    def initializeMolForCmd(self, mol):
        """
        """
        self.initializedFor[mol] = True
        geomC = mol.geomContainer
        mol._cartoonGeoms = {}
        chids = numpy.unique(mol._ag.getChids())
        #if mol._multi=="conformations":

        if mol._colors.get('cartoon', None) is None:
            # per atom properties used for lines
            # initialize with line colors
            mol._colors['cartoon'] = mol._colors['lines'].copy()

        if mol._ag.getData('colorsIndices_cartoon') is None:
            mol._ag.setData('colorsIndices_cartoon',
                            mol._ag.getData('colorsIndices_lines').copy())

        mol._cartoonData = [None]*mol._ag._coords.shape[0]
        if not geomC.geoms.has_key('cartoon'):
            master = Geom('cartoon', inheritLineWidth=False)
            geomC.addGeom( master, parent=geomC.masterGeom, redo=0, makeAtomSet=False,
                           displayCmd=self.app().displayCartoon, undisplayCmd=self.app().undisplayCartoon)
            for chid in chids:
                geom = IndexedPolygons("chain_%s"%chid, inheritMaterial=False, inheritLineWidth=False,
                                       )
                geomC.addGeom( geom, parent=master, redo=0,
                               displayCmd=self.app().displayCartoon, undisplayCmd=self.app().undisplayCartoon)
                mol._cartoonGeoms[chid] = geom

    def pickedVerticesToAtoms(self, geom, vertInd):
        mol = geom.mol()
        geomC = mol.geomContainer
        confNum = mol._ag.getACSIndex()
        cartoonDataDict = mol._cartoonData[confNum]
        for chnum, chid in enumerate(numpy.unique(mol._ag.getChids())):
            g = cartoonDataDict[chid][0]
            if g == geom:
                v2At = cartoonDataDict[chid][4]
                atIndices = numpy.unique(v2At[list(vertInd)])
                res = numpy.unique(Selection(mol._ag, atIndices, "").getResindices())
                resStr = str(res.tolist()).replace(',',' ')[1:-1]
                atoms = mol._ag.select('resindex %s'%resStr)
                atmInds = atoms.getIndices()
                return atmInds

    def gaps(self, coords, cutOff):
        v = coords[1:]-coords[:-1]
        d = numpy.sum(v*v, axis=1)
        mask = d >= cutOff
        gapIndices = (numpy.where(mask==True)[0])
        gapDistances = d[gapIndices]
        return (gapIndices+1).tolist(), gapDistances

    def doit(self, molSel,  helixShape=None, sheetShape=None, coilShape=None, shapeNABB=None, quality=3, scale=1.):
        """
        compute the spline path for the entire molecules and extrude 2D shapes
        """
        mol = molSel.getAtomGroup().getMolecule()
        molSel = molSel.select("not deleted")
        isInitialized = self.initializedFor.get(mol, False)
        if not isInitialized:
            self.initializeMolForCmd(mol)

        # assign secondary structure if needed
        if mol._ag.getSecstrs() is None:
            mol.assignSecondaryStructureWithPross()            

        if helixShape is None:
            helixShape = Ellipse2D(0.6*scale, 0.1, quality=quality)
        if sheetShape is None:
            sheetShape = Rectangle2D(width=1.2*scale, height=0.2*scale, vertDup=1)
        if coilShape is None:
            coilShape = Circle2D(radius=0.1*scale, quality=quality)
        if shapeNABB is None:
            shapeNABB = Circle2D(radius=0.3*scale, quality=quality)
        if mol._renderingProp.get('cartoon', None) is None:
                mol._renderingProp['cartoon'] = {
                    }
        mol._renderingProp['cartoon']['helixShape'] = helixShape
        mol._renderingProp['cartoon']['sheetShape'] = sheetShape
        mol._renderingProp['cartoon']['coilShape'] = coilShape
        mol._renderingProp['cartoon']['shapeNABB'] = shapeNABB
        mol._renderingProp['cartoon']['quality'] = quality
        mol._renderingProp['cartoon']['scale'] = scale
        geomC = mol.geomContainer
        sf = self._splineFrames
        cartoonDataDict = {}
        confNum = mol._ag.getACSIndex()
        #for chid in numpy.unique(mol._ag.getChids()):
        for chain in mol._ag.getHierView():
            chid = chain.getChid()
            if len(chain) < 4:
                cartoonDataDict[chid] = None
                continue
            # calculate 2D sheet
            allcaAtoms = chain.select("ca")
            if allcaAtoms is None:
                allpAtoms = chain.select("nucleotide and element P and not deleted")
                if allpAtoms is None:
                    cartoonDataDict[chid] = None

                else:
                    # build cartoon for Nucleic Acids
                    # look for gaps
                    gapIndices, gd = self.gaps(allpAtoms.getCoords(), 64.) # 8.0**2 P-P
                    gapIndices.insert(0, 0)
                    gapIndices.append(len(allpAtoms._indices))
                    for n, gindex in enumerate(gapIndices[1:-1]):
                        ind = allpAtoms._indices[gindex]
                        data = mol._ag._data
                        print 'GAP in chain %s after residue %s%d (distance %.2f)'%(
                            chid, data['resname'][ind], data['resnum'][ind], sqrt(gd[n]))

                    verts = []
                    faces = []
                    normals = []
                    vert2At = []
                    faces4At= {}
                    vc = 0
                    fc = 0
                    cylCoords =[]
                    cylcc = 0
                    facesDict = {}
                    allCtrlAtoms = mol.emptySelection()

                    #allNAResnums = chain.select('nucleotide and not deleted').getResnums()
                    #withP = numpy.unique(chain.select("nucleotide and element P and not deleted").getResnums())
                    #allNA = numpy.unique(chain.select("nucleotide").getResnums())
                    #noP = 
                    allNAResnums = numpy.unique(chain.select('nucleotide and not deleted').getResnums()).tolist()

                    for gindex, start in enumerate(gapIndices[:-1]):
                        end = gapIndices[gindex+1]

                        patoms = Selection(mol._ag, allpAtoms._indices[start:end], '')
                        ## Assumption: we only cartoon consecutive segments of NA containing a P atom
                        ##             but check if there is 1 NA before and after without a P and use O5'
                        ##             at the begining and O3' at the end in these cases
                        ctrl, ctrlAtoms = mol.getChainCtrlPointsForNA(patoms, chain, allNAResnums)
                        # ctrl has 2 more points than ctrlAtoms for spline interpolation
                        #print chid, ctrlAtoms.getNames()
                        # add tot he set of control atoms
                        allCtrlAtoms = allCtrlAtoms | ctrlAtoms
                        nbres = len(numpy.unique(ctrlAtoms.getResindices()))
                        resIndices = ctrlAtoms.getResindices()
                        atIndices = ctrlAtoms.getIndices()
                        
                        resNums = ctrlAtoms.getResnums()
                        smoothPoints, ids = ResampleControlPoints(ctrl, quality, atIndices)
                        pp1, norm, binormals = GetFrames(smoothPoints)
                        sstype = ['C']*len(smoothPoints)
                        v, f, n, v2A, f4A = self._extruder(
                            smoothPoints, norm, binormals, quality, ids, pp1[0],
                            sstype, coilShape, coilShape, shapeNABB, fc)
                        verts.extend(v)
                        faces.extend((numpy.array(f)+vc).tolist())
                        normals.extend(n)
                        vert2At.extend(v2A)
                        faces4At.update(f4A)
                        vc += len(v)
                        fc += len(f)

                        ##
                        ## add cylinder and sphere geoms for ladder steps

                        # compute middle of P-P as base for cylinder
                        mcoords = (0.5*(ctrl[1:-2]+ctrl[2:-1])).tolist()
                        
                        #import pdb; pdb.set_trace()
                        for ii, resIndex in enumerate(numpy.unique(resIndices)):
                            res = chain[resNums[ii]]
                            resname = res.getResname()
                            cylCoords.append(mcoords[ii])
                            if resname=='A' or  resname=='G' or resname=='DA' or resname=='DG':
                                cylCoords.append(res.select("nucleotide and name N1").getCoords()[0])
                            else:
                                cylCoords.append(res.select("nucleotide and name N3").getCoords()[0])
                            facesDict[resIndex] = ( 2*ii+cylcc, 2*ii+cylcc+1)
                        cylcc = len(cylCoords)
                        
                    geom = mol._cartoonGeoms[chid]
                    geom._nuclotide = True # used by display to decide if ladder bar have
                                           # to be displayed
                    geom._allAtoms = Selection(mol._ag, allCtrlAtoms._indices, "p") # used by display
                    geom.Set(vertices=verts, vnormals=normals, faces=[], visible=0)
                    cartoonDataDict[chid] = (geom, numpy.array(verts), numpy.array(faces),
                                             numpy.array(normals), numpy.array(vert2At),
                                             faces4At)
                    geomC.geomPickToAtoms["chain_%s"%chid] = self.pickedVerticesToAtoms
                    geomC.geomPickToBonds["chain_%s"%chid] = None

                    cyl = Cylinders('steps_%s'%chid, vertices=cylCoords, faces=[], radii=0.2)
                    cyl._facesDict = facesDict
                    self.app().gui().viewer.AddObject(cyl, parent=geom)

                    sph = Spheres('stepsCap_%s'%chid, vertices=cylCoords, faces=[], radii=(0.2,))
                    self.app().gui().viewer.AddObject(sph, parent=geom)

                    
            else: # protein cartoon for this chain
                allcaoAtoms = chain.select("protein and name CA O and not deleted")

                # look for gaps
                #gapIndices = self.gaps(allcaAtoms.getCoords(), 16.81) # 4.1**2 CA-CA
                gapIndices, gd = self.gaps(allcaAtoms.getCoords(), 18.5) # 4.1**2 CA-CA
                gapIndices.insert(0, 0)
                gapIndices.append(len(allcaAtoms._indices))
                for n, gindex in enumerate(gapIndices[1:-1]):
                    ind = allcaAtoms._indices[gindex]
                    data = mol._ag._data
                    print 'GAP in chain %s after residue %s%d (distance %.2f)'%(
                        chid, data['resname'][ind], data['resnum'][ind], sqrt(gd[n]))

                # handle each segment
                verts = []
                faces = []
                normals = []
                vert2At = []
                faces4At= {}
                vc = 0
                fc = 0
                for gindex, start in enumerate(gapIndices[:-1]):
                    end = gapIndices[gindex+1]
                    
                    caoAtoms = Selection(mol._ag, allcaoAtoms._indices[2*start:2*end], '')
                    caAtoms = Selection(mol._ag, allcaAtoms._indices[start:end], '')
                    caoCoords = caoAtoms.getCoords()
                    nbres = len(caoCoords)/2
                    sstype = caAtoms.getSecstrs()
                    resIndices = caAtoms.getResindices()
                    atIndices = caAtoms.getIndices()
                    smooth = sf.calculate(caoCoords, sstype=sstype, quality=quality,
                                          ident=atIndices)

                    # extrude shapes along 2D sheet
                    v, f, n, v2A, f4A = self._extruder(
                        sf.path, sf.pathNV, sf.pathBNV, quality, sf.identity, sf._ca1V,
                        sstype, helixShape, sheetShape, coilShape, fc)
                    verts.extend(v)
                    faces.extend((numpy.array(f)+vc).tolist())
                    normals.extend(n)
                    vert2At.extend(v2A)
                    faces4At.update(f4A)
                    vc += len(v)
                    fc += len(f)
                    
                geom = mol._cartoonGeoms[chid]
                geom._nuclotide = False
                geom._allAtoms = Selection(mol._ag, allcaAtoms._indices, "ca")
                geom.Set(vertices=verts, vnormals=normals, faces=[], visible=0)
                cartoonDataDict[chid] = (geom, numpy.array(verts), numpy.array(faces),
                                         numpy.array(normals), numpy.array(vert2At),
                                         faces4At)
                geomC.geomPickToAtoms["chain_%s"%chid] = self.pickedVerticesToAtoms
                geomC.geomPickToBonds["chain_%s"%chid] = None
                    
                
        #if mol._multi=='conformations':
        mol._cartoonData[confNum] = cartoonDataDict
    
from .displayCmds import DisplayCommand
from PmvApp.Pmv import RefreshDisplayEvent

class DisplayCartoon(DisplayCommand):
    """The displayCartoon command displays the cartoon representation computed by computeCartoon for a user specified of atoms.

    Synopsis:
        None <- displayCartoon(atoms, frontRendering='solid', linewidth=2, pointwidth=2)

    arguments:
        atoms    --- set of atoms

        # global parameters: (i.e. affecting all cpk spheres of a molecule)

        frontRendering --- can be 'solid', 'mesh', or 'points'. Default:'solid' 
        backRendering --- can be 'solid', 'mesh', or 'points' or 'sameAsFront'
        backfaceColor --- a color for back faces (numpy array [[r,g,b]]) or the string 'sameAsFront'

        culling  --- face culling; can be 'front', 'back', 'none'
        linewidth --- integer > 1 specifying the with of lines for drawing spheres as a mesh

        pointwidth --- integer > 1 specifying the size of poitns for drawing spheres as points

    Package : PmvApp
    Module  : cartoonCmds
    Class   : DisplayCartoon
    Command name : displayCartoon
    """

    _argNames = ['frontRendering', 'backRendering', 'culling', 'backfaceColor', 'linewidth', 'pointwidth']
    _cullingModes = ['front', 'back', 'none']

    def onAddCmdToApp(self):
        self.app().eventHandler.registerListener(RefreshDisplayEvent, self.refreshDisplay_cb)

    def __init__(self):
        DisplayCommand.__init__(self)
        self._atomSetName = 'cartoon' # name of the atom set in GeomContainer.geoms
                                      # indicating atoms for which this representation is shown

    def initializeMolForCmd(self, mol):
        if self.initializedFor.get(mol, False):
            return
        self.initializedFor[mol] = True
        if self.app().undisplayCartoon.isLoader():
            self.app().undisplayCartoon.loadCommand()
        self.app().undisplayCartoon.initializedFor[mol] = True
        if not mol._renderingProp.has_key('cartoon'):
            mol._renderingProp['cartoon'] = {}
        prop = mol._renderingProp['cartoon']
        prop['linewidth'] = 2 # line width
        prop['pointwidth'] = 2 # point size
        prop['frontRendering'] = 'fill'
        prop['backRendering'] = 'fill'
        prop['backfaceColor'] = [ 1.,  1.,  1., 1.]
        prop['culling'] = 'back'
        
    def updateModelGlobals(self, mol, propDict=None, frontRendering=None, backRendering=None, culling=None, backfaceColor=None, linewidth=None, pointwidth=None, **kw):
        self.initialize(mol)
        if propDict is None:
            propDict = mol._renderingProp['cartoon']
        if frontRendering is not None:
            assert frontRendering in self._renderToDejaVu.keys()
            propDict['frontRendering'] = self._renderToDejaVu[frontRendering]
        if backRendering is not None:
            assert backRendering in self._renderToDejaVu.keys()+['sameAsFront']
            propDict['backRendering'] = self._renderToDejaVu[backRendering]
        if culling is not None:
            assert culling in self._cullingModes
            propDict['culling'] = culling
        if  backfaceColor is not None:
            assert backfaceColor=='sameAsFront' or len(backfaceColor) in [3,4]
            propDict['backfaceColor'] = backfaceColor
        if linewidth is not None:
            assert isinstance(linewidth, int) and linewidth >=1
            propDict['linewidth'] = linewidth
        if pointwidth is not None:
            assert isinstance(pointwidth, int) and pointwidth >=1
            propDict['pointwidth'] = pointwidth

    def updateModel(self, sel, **kw):
        mol = sel.getAtomGroup().getMolecule()
        self.initialize(mol)
        gc = mol.geomContainer

        if mol._multi == 'molecules':
            if mol._renderingBits & mol._CARTOON: # if set
                mol._renderingBits -= mol._CARTOON # remove the bit

        ## CHECKME will this work if atoms are deleted ??
        if len(sel)==len(mol._ag):
            if not self.negate:
                if mol._multi == 'molecules':
                    mol._renderingBits += mol._CARTOON

        cartoonData = mol._cartoonData[mol._ag.getACSIndex()]
        confNum = mol._ag.getACSIndex()
        for chain in mol._ag.getHierView():
            chid = chain.getChid()
            data = cartoonData[chid]
            if data is None:
                continue
            geom, v, f, n, v2a, f4a = data
            
            if len(sel)==len(mol._ag):
                if self.negate: #if negate, remove current atoms from displayed _set
                    _set = mol.emptySelection()
                else:
                    _set = geom._allAtoms
            else:
                _set = gc.atoms["chain_%s"%chid]
                if self.negate: #if negate, remove current atoms from displayed _set
                    _set = _set - sel
                else:
                    _set = _set | (sel & geom._allAtoms)
            mol.geomContainer.setAtomsForGeom("chain_%s"%chid, _set)

    def doit(self, molSel, **kw):
        sel = molSel.select('not deleted')
        mol = sel.getAtomGroup().getMolecule()

        self.initialize(mol)
        
        self.updateModelGlobals(mol, **kw)
        self.updateModel(sel, **kw)

        self.refreshDisplay(mol)

    def refreshDisplay_cb(self, event):
        self.refreshDisplay(event.molecule)
    
    def refreshDisplay(self, mol):
        if not self.isInitialized(mol):
            return

        gc = mol.geomContainer

        cartoonData = mol._cartoonData[mol._ag.getACSIndex()]
        confNum = mol._ag.getACSIndex()

        # get a list of all atom colors for lines
        colorInds = mol._ag.getData('colorsIndices_cartoon')
        colors = mol._colors['cartoon'][colorInds.tolist()]
        globProp = mol._renderingProp['cartoon']
        bfcolor = globProp['backfaceColor']
        for chain in mol._ag.getHierView():
            chid = chain.getChid()
            data = cartoonData[chid]
            if data is None:
                continue
            geom, v, f, n, v2a, f4a = data
            faceIndices = []
            
            _set = gc.atoms["chain_%s"%chid]

            for index in _set.getIndices():
                faceIndices.extend(f4a[index])

            # identify colors for faces
            mat = []
            for index in _set.getIndices():
                mat.extend( [colors[index]]*len(f4a[index]) )

            if len(faceIndices)>0:
                if bfcolor=='sameAsFront':
                    geom.frontAndBack=True
                    polyFace = 'front'
                else:
                    geom.frontAndBack=False
                    geom.Set(materials=[bfcolor], polyFace='back')
                    polyFace = 'front'
                geom.Set(visible=True, faces=f[faceIndices],
                         materials=mat, lineWidth=globProp['linewidth'],
                         frontPolyMode=globProp['frontRendering'],
                         backPolyMode=globProp['backRendering'],
                         culling=globProp['culling'], polyFace=polyFace)
            else:
                geom.Set(visible=False)

            # handle DNA ladder
            if geom._nuclotide is True:
                ## geom has 2 children geometries (cylinders and spheres)
                ## we need to set faces for them to only show ladder steps for
                ## atoms in _set
                cyl = geom.children[0]
                if isinstance(cyl, Cylinders):
                    sph = geom.children[1]
                else:
                    sph = cyl
                    cyl = geom.children[1]
                facesD = cyl._facesDict
                facesCyl = []
                facesSph = []
                for index in _set.getResindices():
                    facesCyl.append(facesD[index])
                    facesSph.extend(facesD[index])
                if bfcolor=='sameAsFront':
                    cyl.frontAndBack=True
                    sph.frontAndBack=True
                    polyFace = 'front'
                else:
                    cyl.frontAndBack=False
                    cyl.Set(materials=[bfcolor], polyFace='back')
                    sph.frontAndBack=False
                    sph.Set(materials=[bfcolor], polyFace='back')
                    polyFace = 'front'    
                cyl.Set(faces=facesCyl,
                        lineWidth=globProp['linewidth'],
                        frontPolyMode=globProp['frontRendering'],
                        backPolyMode=globProp['backRendering'],
                        culling=globProp['culling'],
                        polyFace=polyFace)
                sph.Set(faces=[[x] for x in numpy.unique(facesSph)],
                        lineWidth=globProp['linewidth'],
                        frontPolyMode=globProp['frontRendering'],
                        backPolyMode=globProp['backRendering'],
                        culling=globProp['culling'], polyFace=polyFace)


class UndisplayCartoon(DisplayCartoon):
    """The undisplayCartoon command removes the cartoon repesentation for the specified set ot atoms

    Synopsis:
        None <- undisplayCartoon(atoms)

    arguments:
        atoms  --- set of atoms

    Package : PmvApp
    Module  : cartoonCmds
    Class   : UndisplayCartoon
    Command name : undisplayCartoon
    """

    def __init__(self):
        DisplayCartoon.__init__(self)
        self.negate = True

    ## def initializeMolForCmd(self, mol):
    ##     if not self.initializedFor[mol]:
    ##         return
            #raise RuntimeError("undisplayCartoon called before displayCartoon for molecule %s"%mol.name)
    
    def onAddCmdToApp(self):
        DisplayCommand.onAddCmdToApp(self)
        if not self.app().commands.has_key('displayCartoon'):
            self.app().lazyLoad('msmsCmds', commands=['displayCartoon'], package='PmvApp')
            self.displayCartoon.loadCommand()

commandClassFromName = {
    'computeCartoon' : [ComputeCartoon,  None],
    'displayCartoon': [DisplayCartoon, None],
    'undisplayCartoon': [UndisplayCartoon, None]
}

