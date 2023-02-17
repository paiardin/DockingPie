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
# Copyright: M. Sanner TSRI 2014
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/PmvApp/msmsCmds.py,v 1.7.4.1 2017/07/13 20:55:28 annao Exp $
# 
# $Id: msmsCmds.py,v 1.7.4.1 2017/07/13 20:55:28 annao Exp $
#

import numpy
import os, sys
from MolKit2.molecule import Atom
from MolKit2.selection import Selection, SelectionSet

from DejaVu2.IndexedPolygons import IndexedPolygons

from PmvApp.displayCmds import DisplayIndexedGeomCommand
from PmvApp.Pmv import MVCommand  #MVAtomICOM
#from PmvApp.msmsParser import MSMSParser

from PmvApp.Pmv import DeleteAtomsEvent, DeleteObjectEvent, RefreshDisplayEvent
from PmvApp import Pmv

class ComputeMSMS(MVCommand): #, MVAtomICOM):
    """The computeMSMS command will compute a triangulated solvent excluded surface for the given atoms. The surface will always be a closed surface enclosing the set given set of atoms. If the molecule already has a surface with the same name, the command will recompute the surface and replace it unless none of the parameters has changed.

    Synopsis:
        None <- computeMSMS(atoms, surfName=None, probeRadius=1.5, density=3.0,
                            radii=None, hdset=None, hdensity=6.0, cavities=True)

    Required Arguments :
        atoms    --- set of atoms

    Optional Arguments:
    surfName --- name for the surface. If this naem is not given it will
                 default to the molecule's name followed by '_surface'.
                 If this surface name has been used before for this moleculeal
                 the newly computed surface will replace the previous one.

    probeRadius --- probe radius (default: 1.5)

    density  --- density of vertices used to triangle the surface (3.0)

    cavities --- when False only the outter surface is calculated. Default False

    radii    --- If specified this argument has to be a list of floating points of the same
                 length as atoms, providing a radius for each atom. When this argument is not
                 specified the atoms van der Waals radii are used.
              
    hdset    --- A Selection instance or selection name for which
                 high density triangualtion will be generated

    hdensity --- vertex density for high density part of the surface

    Package : PmvApp
    Module  : msmsCmds
    Class   : ComputeMSMS
    Command name : computeMSMS
    """

    ## mol.geomContainer.geoms are called msms_%s_surface mol._basename
    ## mol.atoms
    _argNames = ['surfName', 'probeRadius', 'density', 'radii',
                 'hdset', 'hdensity', 'closedSurface', 'cavities']

    def checkArguments(self, nodes, **kw):
        for name in kw.keys():
            if name not in self._argNames:
                raise RuntimeError("%s: unrecognized keyword argument '%s', valid names are %s"%(
                    self.name, name, str(self._argNames)))
        return (nodes,), kw

    def isInitialized(self, mol):
        return self.initializedFor.get(mol, False)
    
    def initialize(self,mol):
        # call self.initializedFor is needed
        isInitialized = self.initializedFor.get(mol, False)
        if not isInitialized:
            self.initializeMolForCmd(mol) 

    def __init__(self):
        MVCommand.__init__(self)
        self.MSMSdisabled = False


    #def checkDependencies(self, vf):
    #    import mslib

    def onAddCmdToApp(self):
        MVCommand.onAddCmdToApp(self)
        app = self.app()
        app.eventHandler.registerListener(DeleteAtomsEvent,
                                          self.handleDeleteAtoms)
        
        if not  hasattr( app, 'numOfSelectedVerticesToSelectTriangle'):
            app.numOfSelectedVerticesToSelectTriangle = 1

        ## app.userpref.add('Number of selected vertices per triangle', 'choice', (1,2,3),
        ##     doc="""specifies how many vertices in a triangle nee to be selected for thsi triangle to be selected when selecting triangles for a set of atoms""")

        app.userpref.add('Color Boundaries for MSMS', 'blurred', ('sharp','blurred'),
            doc="""specifies whether color boundaries are sharp ('yes') or blurred ('no')""")
        app.userpref.add('Compute cavities by default', 'no', ('yes','no'),
            doc="""specifies whether to compute cavity surfaces by default ('yes', 'no')""")

    def handleDeleteAtoms(self, event):
        mol = event.deletedAtoms.getAtomGroup().getMolecule()
        if not self.initializedFor.get(mol, False):
            return
        for name in mol._msmsData['msms'].keys():
            atoms = mol._msmsData['atoms'][name] & event.deletedAtoms
            if len(atoms) == 0: # all atoms for thsi surface were deleted
                del mol._msmsData['atoms'][name]
                del mol._msmsData['msms'][name]
                del mol._msmsData['params'][name]
                
            elif len(mol._msmsData['atoms'][name] & event.deletedAtoms)>0:
                # if some deleted atoms contribute to this surface
                # recompute the surface
                self(mol._msmsData['atoms'][name], surfName=name)

        self.app().displayMSMS.refreshDisplay(mol)

    def handleDeleteObject(self, event):
        if event.objectType=='Molecule':
            mol = event.object
            if hasattr(mol, '_msmsData'):
                mol._msmsData['msms'] = {}
                mol._msmsData['atoms'] = {}
                mol._msmsData['indexMaps'] = {}
        
    def initializeMolForCmd(self, mol):
        """
        """
        if self.initializedFor.get(mol, False):
            return
        self.initializedFor[mol] = True
        
        if mol._colors.get('msms', None) is None:
            # per atom properties used for lines
            mol._colors['msms'] = {}

        if not hasattr(mol, "_msmsData"):
            mol._msmsData = {
                'msms': {}, # stores surfName:[msmsCObj] of len mol._ag.numCoordsets()
                'atoms': {}, # stores surfName:atomSet used to compute this surface
                'params': {}, # stores surfName: {'probeRadius:val, density:val, etc..}
                'indexMaps': {}, # holds surfName: indexMap where indexMap[i] = j
                                 # means atom i in molecule is at index j in set
                                 # used to compute this surface
                }

    def recomputeMSMS(self, mol, molSel, surfName, probeRadius, density,
                      radii, hdset, hdensity, cavities):
        params = mol._msmsData['params'][surfName]
        recompute = False
        kw = {'probeRadius':probeRadius,
              'density':density,
              'radii':radii,
              'hdset':hdset,
              'hdensity':hdensity,
              'cavities':cavities}

        if probeRadius is None:
            kw['probeRadius'] = params['probeRadius']
        elif params['probeRadius'] !=  probeRadius:
            kw['probeRadius'] = probeRadius
            recompute = True
            
        if density is None:
            kw['density'] = params['density']
        elif params['density'] !=  density:
            kw['density'] = density
            recompute = True

        if cavities is None:
            kw['cavities'] = params['cavities']
        elif params['cavities'] !=  cavities:
            kw['cavities'] = cavities
            recompute = True

        if hdensity is None:
            kw['hdensity'] = params['hdensity']
        elif params['hdensity'] !=  hdensity:
            kw['hdensity'] = hdensity
            recompute = True

        if radii is None:
            kw['radii'] = params['radii']
        elif params['radii'] is None:
            kw['radii'] = radii
            recompute = True
        elif numpy.sum(radii-params['radii'])!=0:
            kw['radii'] = radii
            recompute = True
            
        if hdset is None:
            kw['hdset'] = params['hdset']
        elif params['hdset'] is None:
            kw['hdset'] = hdset
            recompute = True
        elif len(hdset-params['hdset'])>0:
            kw['hdset'] = hdset
            recompute = True

        if len(molSel-mol._msmsData['atoms'][surfName])>0:
            recompute = True

        return recompute, kw

    def doit(self, molSel, surfName=None, probeRadius=None, density=None,
             radii=None, hdset=None, hdensity=None, cavities=None, **kw):

        # check for mslib availability (one time)
        try:
            if not self.MSMSdisabled:
                from mslib import MSMS
            else:
                return
        except RuntimeError:
            self.MSMSdisabled = True
            return

        sel = molSel.select('not deleted')
        mol = sel.getAtomGroup().getMolecule()
        self.initialize(mol)
        # set default values for unspecified parameters
        surfName = "%s_surface"%mol._basename if surfName is None else surfName
        # if there already is a surface with that name
        # 1 - check if we need to recompute
        # 2 - if a calcualtion is needed, get values of unspecified parameters
        #     from the previous surface
        if mol._msmsData['msms'].has_key(surfName):
            if mol._msmsData['msms'][surfName][mol._ag.getACSIndex()] is not None:
                recompute, kw = self.recomputeMSMS(mol, molSel, surfName, probeRadius, 
                                                   density, radii, hdset, hdensity, cavities)
            if not recompute:
                if len(sel)!=len(mol._msmsData['atoms'][surfName]) or \
                       len(sel & mol._msmsData['atoms'][surfName]) < len(sel):
                    recompute = True # atoms changed
                else:
                    #print 'NO RECOMPUTE'
                    return
            probeRadius = kw["probeRadius"]
            density = kw["density"]
            radii = kw["radii"]
            hdset = kw["hdset"]
            hdensity = kw["hdensity"]
            cavities = kw["cavities"]
        else: # set default values for omitted arguments
            if probeRadius is None: probeRadius = 1.5
            if density is None:
                if len(mol._ag)>100:
                    density = 3.0
                else:
                    density = 5.0
            if hdensity is None: hdensity = 8.0
            if cavities is None:
                cavities = mol.app().userpref['Compute cavities by default']['value']=='yes'
            
        app = self.app()
        # get the set of molecules and the set of atoms per molecule in the
        # current selection

        # build flags for HDset
        if hdset is not None:
            assert sel._ag == hdset._ag
            hdflags = numpy.zeros(len(mol._ag), "i")
            hdflags[hdset.getIndices()] = 1

        #if surfName is None:
        #    surfName = mol._basename + '_surface'

        # build a map that yields the index of an atom in sel for every atom in mol
        # i.e. indexMap[i] = j means atom i in molecule is at index j in sel
        indexMap = numpy.ones(len(mol._ag), dtype='i')*(-1)
        indices = sel.getIndices()
        indexMap[indices] = range(len(sel))
        mol._msmsData['indexMaps'][surfName] = indexMap
        surf = [1]*len(sel)
        if hdset is not None:
            hd = hdflags[indices].tolist()
        else:
            hd = [0]*len(sel)

        if radii is None:
            atmRadii = sel.getRadii().copy()
        else:
            assert len(radii)==len(sel) and isinstance(radii[0], float)
            atmRadii = radii
        atmCoords = sel.getCoords().tolist()

        # build an MSMS object and compute the surface
        #print 'MSMS', hdset, hd
        srf = MSMS(coords=atmCoords, radii=atmRadii, surfflags=surf, hdflags=hd)
        srf.compute(probe_radius=probeRadius, density=density,
                    allComponents=cavities, hdensity=hdensity)
        srf.compute_numeric_area_vol()
        
        #import pdb; pdb.set_trace()
        if len(atmCoords)==len(mol._ag):
            # computing surface for the entire molecule, add a flag for
            # solvent accessibility to atoms participating to the surface
            #dum, vi, dum = srf.getTriangles()
            #_vi = vi[:,1]-1
            _vi = numpy.zeros( (0,), 'i')
            for nc in range(srf.sesr.nb):
                dum, vi, dum = srf.getTriangles(component=nc)
                _vi = numpy.concatenate((_vi, vi[:,1]-1))
            sa = numpy.zeros((len(mol._ag),), 'bool')
            saAtomIndices = numpy.unique(_vi)
            sa[saAtomIndices] = 1
            mol._ag.setFlags('solventAccessible', sa)
        #s = srf.sesr.fst;
        #while s:
        #    print s.n_ses_volume, s.a_ses_volume, s.n_ses_area, s.a_ses_area 
        #    s = s.nxt
        #import pdb; pdb.set_trace()
        
        # Create a new geometry if needed
        geomC = mol.geomContainer
        #if this is a new surface
        if not mol._msmsData['msms'].has_key(surfName):
            # added default for renderinProps
            if mol._renderingProp.get('msms', None) is None:
                mol._renderingProp['msms'] = {
                    } # dict will be added for each surface of that molecule
                
                mol._renderingProp['msms'][surfName] = {'linewidth':2, # line width
                                                        'pointwidth':2, # point size
                                                        'frontRendering':'fill', #
                                                        'backRendering':'line', #
                                                        'quality':None,
                                                        'culling':'none',
                                                        'shading': 'smooth',
                                                        'backfaceColor': [ 1.,  1.,  1., 1.],
                                                        'components': 'all', # of list of component indices
                                                        }
                mol._renderingProp['msms'][surfName]['compVisible'] = [True]*srf.sesr.nb
                mol._renderingProp['msms'][surfName]['flip'] = [False]*srf.sesr.nb

            # create a  geometry for it
            g = IndexedPolygons('msms_'+surfName, pickableVertices=1, protected=True,
                                inheritMaterial=False, inheritLineWidth=0, 
                                inheritFrontPolyMode=False)
            if self.app().userpref['Color Boundaries for MSMS']['value'] == 'blurred':
                g.Set(inheritSharpColorBoundaries=False, sharpColorBoundaries=False,)
                mol._renderingProp['msms'][surfName]['sharpColorBoundaries'] = False

            else:
                g.Set(inheritSharpColorBoundaries=False, sharpColorBoundaries=True)
                mol._renderingProp['msms'][surfName]['sharpColorBoundaries'] = True
            g._pmvName = surfName
            g._isMSMS = True
            geomC.addGeom(g,  displayCmd=self.app().displayMSMS,
                          undisplayCmd=self.app().undisplayMSMS)
            self.managedGeometries.append(g)
            geomC.geomPickToAtoms[g.name] = self.pickedVerticesToAtoms
            #geomC.geomPickToBonds[g.name] = None
            geomC.setAtomsForGeom(surfName, sel)
            # add color fields and dict
            mol._colors['msms'][surfName] = mol._colors['lines'].copy()
            mol._ag.setData('colorsIndices_msms_%s'%surfName,
                            mol._ag.getData('colorsIndices_lines').copy())
        else:
            mol._renderingProp['msms'][surfName]['compVisible'] = [True]*srf.sesr.nb
            mol._renderingProp['msms'][surfName]['flip'] = [False]*srf.sesr.nb
            
        # if the surface was recompute for this conformation all surfaces of
        # other conformations are now invalid so we reset the verctor [srf]
        mol._msmsData['msms'][surfName] = [None]*mol._ag.numCoordsets()
        mol._msmsData['msms'][surfName][mol._ag.getACSIndex()] = srf
        mol._msmsData['atoms'][surfName] = sel.copy()
        mol._msmsData['params'][surfName] = {
            'probeRadius':probeRadius, 'density': density,'hdensity':hdensity,
            'cavities':cavities}
        if hdset is not None:
            mol._msmsData['params'][surfName]['hdset'] = hdset._indices.copy(),
        else:
            mol._msmsData['params'][surfName]['hdset'] = None
        if radii is not None:
            mol._msmsData['params'][surfName]['radii'] = radii
        else:
            mol._msmsData['params'][surfName]['radii'] = None
            
    def pickedVerticesToBonds(self, geom, parts, vertex):
        return None

    def pickedVerticesToAtoms(self, geom, vertInd):
        """Function called to convert picked vertices into atoms"""

        # this function gets called when a picking or drag select event has
        # happened. It gets called with a geometry and the list of vertex
        # indices of that geometry that have been selected.
        # This function is in charge of turning these indices into an AtomSet
        surfName = geom._pmvName
        mol = geom.mol()
        geomC = mol.geomContainer
        coordSetIndex = mol._ag.getACSIndex()
        surf = mol._msmsData['msms'][surfName][coordSetIndex]
        al = geomC.atoms[surfName].select('not deleted')
        surfinds = mol._msmsData['indexMaps'][surfName]
        indices = al.getIndices()
        atomindices = surfinds[indices]
        #print 'MINI2', min(atomindices)
        #dum1, vi, dum2 = surf.getTriangles(atomindices, keepOriginalIndices=1)
        allAtInds = mol._msmsData['atoms'][surfName].getIndices()
        return allAtInds[ numpy.unique(geomC.geoms['msms_'+surfName]._vertexToAtomInd[list(vertInd)]) ]
        #return [allAtInds[vi[i][1]-1] for i in vertInd]

class DisplayMSMS(DisplayIndexedGeomCommand):
    """The displayMSMS command displays molecular surfaces computed using MSMS.

    Synopsis:
        None <- displayLines(atoms, surfNames=None, 
                             frontRendering='solid', backRendering='mesh',
                             linewidth=2, pointwidth=2, quality=6,
                             shading='smooth', culling='back', backfaceColor=(1,1,1,1))

    arguments:
        atoms    --- set of atoms

        surfNames --- a list of names of surfaces on which to operate. Multiple surfaces
                     can be computed for a molecule. They are distiguished by their names

        components --- 'all' or a list of integers providing indices of closed surfaces to be displayed

        # global parameters: (i.e. affecting the entire molecular surface)

        frontRendering -- can be 'solid', 'mesh', or 'points'. Default:'solid'(rendering mode of the front and back polygons) 

        backRendering --- can be 'solid', 'mesh', or 'points' or 'sameAsFront'
        
        linewidth  --- integer > 1 specifying the with of lines for drawing spheres as a mesh

        pointwidth --- integer > 1 specifying the size of poitns for drawing spheres as points

        quality    --- integer specifying the number of certices per square inch used to triangulate the surface.
        culling  --- face culling; can be 'front', 'back', 'none'

        shading  --- shading technique; can be 'flat' or 'smooth'.

    Package : PmvApp
    Module  : msmsCmds
    Class   : DisplayMSMS
    Command name : displayMSMS
    """

    _argNames = ['surfNames', 'components', 'frontRendering',
                 'backRendering', 'linewidth', 'pointwidth',
                 'quality', 'backfaceColor', 'shading', 'culling',
                 'sharpColorBoundaries']

    def __init__(self):
        DisplayIndexedGeomCommand.__init__(self)
        self._atomSetName = None # FIXME

    def initializeMolForCmd(self, mol):
        if self.initializedFor.get(mol, False):
            return
        self.initializedFor[mol] = True
        self.app().undisplayMSMS.initializedFor[mol] = True
        if self.app().undisplayMSMS.isLoader():
            self.app().undisplayMSMS.loadCommand()
            
    def doit(self, molSel, **kw):
        sel = molSel.select('not deleted')
        mol = sel.getAtomGroup().getMolecule()

        self.initialize(mol)
        
        indices = sel.getIndices()

        surfNames = kw.get('surfNames', None)
        if surfNames is None:
            if not hasattr(mol, '_msmsData'):
                return # no surface for thsi molecule
            surfNames = mol._msmsData['msms'].keys()

        for name in surfNames:
            self.updateModelGlobals(mol, name, mol._renderingProp['msms'][name], **kw)
            self.updateModel(sel, name, **kw)

        self.refreshDisplay(mol, names=surfNames)
        self.app().labelAtoms.refreshDisplay(mol) # to make the labels be in front of MSMS
        self.app().labelResidues.refreshDisplay(mol) # to make the labels be in front of MSMS
        
    def updateModelGlobals(self, mol, name, propDict, components=None,
                           frontRendering=None, backRendering=None,
                           linewidth=None, pointwidth=None, quality=None,
                           shading=None, culling=None, backfaceColor=None, sharpColorBoundaries=None, **kw):
        self.initialize(mol)

        DisplayIndexedGeomCommand.updateModelGlobals(self, mol, propDict,
                           frontRendering=frontRendering, backRendering=backRendering,
                           linewidth=linewidth, pointwidth=pointwidth, quality=quality,
                           shading=shading, culling=culling, backfaceColor=backfaceColor,
                           sharpColorBoundaries=sharpColorBoundaries)

        if components is not None:
            assert components=='all' or isinstance(components, list)
            if isinstance(components, list):
                srf = mol._msmsData['msms'][name][mol._ag.getACSIndex()]
                for v in components:
                    assert isinstance(v, int) and v < srf.rsr.nb
            mol._renderingProp['msms'][name]['components'] = components

    def updateModel(self, sel, name, **kw):
        mol = sel.getAtomGroup().getMolecule()
        self.initialize(mol)

        # restrict sel to what was used to compute this surface
        sel = sel & mol._msmsData['atoms'][name]

        if mol._multi == 'molecules':
            if mol._renderingBits & mol._SURFACE: # if set
                mol._renderingBits -= mol._SURFACE # remove the bit

        gc = mol.geomContainer    
        _set = gc.atoms['msms_'+name]
        if len(sel) == len(mol._ag):
            if self.negate:
                _set = mol.emptySelection()
            else:
                _set = sel
                if mol._multi == 'molecules':
                    mol._renderingBits += mol._SURFACE # set the bit
        else:
            ##if negate, remove current atms from displayed _set
            if self.negate:
                _set = _set - sel
            else:
                _set = _set | sel
        # at this point _set is the set of atoms that will still be displayed as MSMS 
        gc.setAtomsForGeom('msms_'+name, _set)

    def refreshDisplay(self, mol, names=None):
        if not self.isInitialized(mol):
            return

        if names is None:            
            names = mol._msmsData['msms'].keys()

        for name in names:
            gname = 'msms_'+name
            geom = mol.geomContainer.geoms[gname]
            atoms = mol.geomContainer.atoms.get(gname, mol.emptySelection())
            if len(atoms)==0: # nothing is displayed as CPK anymore for this molecule
                geom.Set(visible=0, tagModified=False) # hide the geom
            else:
                atoms = atoms.select('not deleted')
                selected = self.app().activeSelection & SelectionSet([atoms])

                # get the msms surface object for that molecule
                surfinds = mol._msmsData['indexMaps'][name]
                indices = atoms.getIndices()
                atomindices = surfinds[indices]
                if len(atomindices) == 0:
                    g.Set(visible=0, tagModified=False)
                else:
                    globProp = mol._renderingProp['msms'][name]
                    srf = mol._msmsData['msms'][name][mol._ag.getACSIndex()]
                    vf = []
                    f = []
                    vi = []
                    if globProp['components']=='all':
                        compIndices = range(srf.rsr.nb)
                    else:
                        compIndices = globProp['components']
                    if len(globProp['compVisible']) < len(compIndices) and len(globProp['compVisible']) == 1:
                        vis = globProp['compVisible'][0]
                        globProp['compVisible'] = [vis]*len(compIndices)
                    nbvInComp = []
                    for compNum in compIndices:
                        if globProp['compVisible'] and not globProp['compVisible'][compNum]: continue
                        vfc, vic, fc = srf.getTriangles(atomindices, component=compNum,
                                                     selnum=1, keepOriginalIndices=1)
                        nbvInComp.append(len(vfc))
                        nf = fc+len(vf) #offset the faces
                        if globProp['flip'][compNum]:
                            # inverse face vertices order
                            tmp = nf[:,1].copy()
                            nf[:,1] = nf[:,2]
                            nf[:,2] = tmp
                            vfc[:, 3:6] = -vfc[:, 3:6] # flip normals
                        f.extend(nf.tolist())
                        vf.extend(vfc.tolist())
                        vi.extend(vic.tolist())
                    if not len(vf):
                        geom.Set(visible=0, tagModified=False)
                        return
                    f = numpy.array(f)
                    vf = numpy.array(vf)
                    vi = numpy.array(vi)
                    geom._vertexToAtomInd = vi[:,1]-1
                    highlight = numpy.zeros( (len(vf),), 'i')
                    geom._nbvInComp = nbvInComp
                    if selected:
                        selAtomindices = surfinds[selected[0].getIndices()]
                        fs = []
                        for cn, compNum in enumerate(compIndices):
                            if globProp['compVisible'] and not globProp['compVisible'][compNum]: continue
                            dum, dum, fsc = srf.getTriangles(selAtomindices, component=compNum,
                                                             selnum=3, keepOriginalIndices=1)
                            if len(fsc):
                                if cn>0:
                                    fsc = fsc+nbvInComp[cn-1]
                                fs.extend(fsc.tolist())
                        fs = numpy.array(fs)
                        if len(fs):
                            highlight[numpy.unique(fs.flatten())] = 1
                    colorInds = mol._ag.getData('colorsIndices_msms_%s'%name)
                    colors = mol._colors['msms'][name][colorInds]
                    opacity = mol._ag.getData('opacity_msms_%s'%name)
                    # turn atom colors into surface vertices colors
                    col = colors[vi[:,1]-1]
                    op = None
                    if opacity is not None:
                        if len(numpy.unique(opacity)) == 1:
                            op = opacity[0]
                        else:
                            op = opacity[vi[:,1]-1]
                    bfcolor = globProp['backfaceColor']
                    if bfcolor=='sameAsFront':
                        geom.frontAndBack=True
                        polyFace = 'front'
                    else:
                        geom.frontAndBack=False
                        geom.Set(materials=[bfcolor], polyFace='back')
                        polyFace = 'front'
                    
                    geom.Set( vertices=vf[:,:3], vnormals=vf[:,3:6],
                              faces=f[:,:3], materials=col, visible=1,
                              frontPolyMode=globProp['frontRendering'],
                              backPolyMode=globProp['backRendering'],
                              culling=globProp['culling'],
                              shading=globProp['shading'],
                              lineWidth=globProp['linewidth'],
                              tagModified=False, polyFace=polyFace,
                              highlight=highlight, opacity=op, transparent="implicit",
                              sharpColorBoundaries=globProp['sharpColorBoundaries'])
                #if g.transparent:
                #    opac = mol.geomContainer.getGeomOpacity(gname)
                #    g.Set(opacity=opac, redo=0, tagModified=False)
                # update texture coordinate if needed
                #if g.texture and g.texture.enabled and g.texture.auto==0:
                #    mol.geomContainer.updateTexCoords[sName](mol)

        # make sure new surfaces reflect the selection
        #self.app().eventHandler.dispatchEvent(RefreshSelectionEvent())
                

class UndisplayMSMS(DisplayMSMS):
    """The undisplayMSMS command removes the molecular surface repesentation for the specified set ot atoms

    Synopsis:
        None <- undisplayMSMS(atoms)

    arguments:
        atoms  --- set of atoms

    Package : PmvApp
    Module  : msmsCmds
    Class   : UndisplayMSMS
    Command name : undisplayMSMS
    """

    def __init__(self):
        DisplayMSMS.__init__(self)
        self.negate = True

    ## def initializeMolForCmd(self, mol):
    ##     if not self.initializedFor[mol]:
    ##         return
            #raise RuntimeError("undisplayMSMS called before displayMSMS for molecule %s"%mol.name)
    
    def onAddCmdToApp(self):
        DisplayMSMS.onAddCmdToApp(self)
        if not self.app().commands.has_key('displayMSMS'):
            self.app().lazyLoad('msmsCmds', commands=['displayMSMS'], package='PmvApp')
            self.displayMSMS.loadCommand()
    

#### The following commands need to be updated ( prody MolKit)  

## class ComputeSESAndSASArea(MVCommand):
##     """Compute Solvent Excluded Surface and Solvent Accessible Surface Areas\n
##     Package : PmvApp\n
##     Module  : msmsCmds\n
##     Class   : ComputeSESAndSASArea\n
##     Command name : computeSESAndSASArea\n
##     \nSynopsis :\n
##     None--->mv.computeSESAndSASArea(mol)\n    
##     \nDescription:\n
##     Computes Solvent Excluded Surface and Solvent Accessible Surface Areas. Stores numeric values per Atom,
##     Residue, Chain, and Molecule in ses_area and sas_area attributes.
##     """
##     def onAddCmdToApp(self):
##         MVCommand.onAddCmdToApp(self)
##         if not self.app().commands.has_key('computeMSMS'):
##             self.app().lazyLoad('msmsCmds', commands=['ComputeMSMS',], package='PmvApp')


##     def doit(self, mol):
##         try:
##             allrads = mol.defaultRadii()
##             allChains = mol.chains
##             allResidues = mol.chains.residues
##             allAtoms = mol.allAtoms
##             import mslib
##             # compute the surface
##             srf = mslib.MSMS(coords=allAtoms.coords, radii = allrads)
##             srf.compute()
##             srf.compute_ses_area()        
##             # get surface areas per atom
##             ses_areas = []
##             sas_areas = []
##             for i in xrange(srf.nbat):
##                 atm = srf.get_atm(i)
##                 ses_areas.append(atm.get_ses_area(0))
##                 sas_areas.append(atm.get_sas_area(0))
##             # get surface areas to each atom
##             allAtoms.ses_area = ses_areas
##             allAtoms.sas_area = sas_areas
##             # sum up ses areas over resdiues
##             for r in allResidues:
##                 r.ses_area = numpy.sum(r.atoms.ses_area)        
##                 r.sas_area = numpy.sum(r.atoms.sas_area)

##             mol.ses_area = 0
##             mol.sas_area = 0            
##             for chain in allChains:
##                 chain.ses_area = 0
##                 chain.sas_area = 0
##                 for residue in chain.residues:
##                     chain.ses_area += numpy.sum(residue.ses_area)
##                     chain.sas_area += numpy.sum(residue.sas_area)
##                 mol.ses_area += chain.ses_area 
##                 mol.sas_area += chain.sas_area
##             self.app()._executionReport.addSuccess('%s: computed surface for molecule %s successfully'%(self.name, mol.name), obj=allAtoms)
##         except:
##             msg = '%s: Error while computing surface for molecule %s'%(self.name, mol.name)
##             self.app().errorMsg(sys.exc_info(), msg, obj=mol)            

            
##     def checkArguments(self, molecule):
##         """
##     Computes Solvent Excluded Surface and Solvent Accessible Surface Areas. Stores numeric values per Atom,
##     Residue, Chain, and Molecule in ses_area and sas_area attributes.               
##         """
##         nodes=self.app().expandNodes(molecule)

##         return (nodes), kw


## class SaveMSMS(MVCommand):
##     """The SaveMSMS command allows the user to save a chosen MSMS surface (tri-angulated solvant excluded surface) in two files: .vert(vertex coordinates) and .face (vertex indices of the faces) .\n
##     Package : PmvApp\n
##     Module  : msmsCmds\n
##     Class   : SaveMSMS\n
##     Command name : saveMSMS\n
##     Description:\n
##     If the component number is 0, files called filename.vert and filename.face
##     are created.
##     For other components, the component number is appended in the file name,
##     for example for the component number 3 the files are called
##     filename_3.vert and filename_3.face.

##     The face file contains three header lines followed by one triangle per
##     line. The first header line provides a comment and the filename of the
##     sphere set.
##     The second header line holds comments about the content of the third line.
##     The third header line provides the number of triangles, the number of
##     spheres in the set, the triangulation density and the probe sphere radius.
##     The first three numbers are (1 based) vertex indices. The next field
##     can be: 1 for a triangle in a toric reentrant face, 2 for a triangle in
##     a spheric reentrant face and 3 for a triangle in a contact face.
##     The last number on the line is the (1 based) face number in the
##     analytical description of the solvent excluded surface. These values
##     are written in the following format ''%6d %6d %6d %2d %6d''.

##     The vertex file contains three header lines (similar to the header
##     in the .face file) followed by one vertex per line and provides the
##     coordinates (x,y,z) and the normals (nx,ny,nz) followed by the number of
##     the face (in the analytical description of the solvent excluded surface)
##     to which the vertex belongs.
##     The vertices of the analytical surface have a value 0 in that field and
##     the vertices lying on edges of this surface have nega tive values.
##     The next field holds the (1 based) index of the closest sphere.
##     The next field is 1 for vertices which belong to toric reentrant faces
##     (including vertices of the analytical surface), 2 for vertices inside
##     reentrant faces and 3 for vertices inside contact faces.
##     Finally, if atom names were present in the input file, the name of the
##     closest atom is written for each vertex. These values are written in
##     the following format
##     ''%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %7d %7d %2d %s''.\n

##     \nSynopsis:\n
##     None <- saveMSMS(filename, mol, surfacename, withHeader=1, component=0,
##                      format='MS_TSES_ASCII', **kw)\n
##     filename : name of the output file\n
##     mol      : molecule associated with the surface\n
##     surfacename : name of the surface to save\n
##     withHeader  : flag to either write the headers or not\n
##     component   : specifies which component of the surface to write out\n
##     format      : specifies in which format to save the surface. 
##                   It can be one of the following ,\n
##                   'MS_TSES_ASCII' Triangulated surface in ASCII format\n
##                   'MS_ASES_ASCII' Analytical surface in ASCII format. This is
##                   actually a discrete representation of the analytical model.\n
##                   'MS_TSES_ASCII_AVS' Triangulated surface in ASCII with
##                   AVS header\n
##                   'MS_ASES_ASCII_AVS'  Analytical surface in ASCII format
##                   with AVS header\n

##     """

##     def onAddCmdToApp(self):
##         MVCommand.onAddCmdToApp(self)
##         self.formats = ['MS_TSES_ASCII', 
##                         'MS_ASES_ASCII', 
##                         'MS_TSES_ASCII_AVS',
##                         'MS_ASES_ASCII_AVS'
##                         ]
##         # MS_TSES_ASCII : Triangulated surface in ASCII format
##         # MS_ASES_ASCII : Analytical surface in ASCII format, which is
##         #                 a discrete representation of the analytical model
##         # MS_TSES_ASCII_AVS : Triangulated surface in ASCII with AVS header
##         # MS_ASES_ASCII_AVS : Analytical surface in ASCII format with AVS
##         #                     header

    
##     def doit(self, filename, molName, surfName, withHeader=True, component=0,
##              format="MS_TSES_ASCII"):
##         try:
##             mol = self.app().getMolFromName(molName)
##             if mol is None:
##                 raise RuntimeError, "saveMSMS: no molecule %s found"% molName
##             gc = mol.geomContainer
##             if not hasattr(gc, "msms") or not gc.msms.has_key(surfName):
##                 raise RuntimeError, "saveMSMS: no molecule surface %s"% surfName
##             msmsSurf = gc.msms[surfName][0]
##             msmsAtms = gc.msms[surfName]
##             from mslib import msms
##             if not format in self.formats:
##                 format = "MS_TSES_ASCII"
##             format = getattr(msms, format)
##             if component is None : component = 0
##             elif not component in range( msmsSurf.rsr.nb ):
##                 raise RuntimeError, "%s error: %s is an invalid component"%(self.name, component)
##             msmsSurf.write_triangulation(filename, no_header=not withHeader,
##                                          component=component, format=format)
##             self.app()._executionReport.addSuccess('saved surface for molecule %s successfully'%
##                     molName, obj=mol)
##         except:
##             msg = 'Error while saving surface for molecule %s'%molName
##             self.app().errorMsg(sys.exc_info(), msg, obj=mol)


##     def checkArguments(self, filename, molName, surfName, withHeader=True,
##                  component=None, format="MS_TSES_ASCII"):
##         """None <--- mv.saveMSMS(filename, mol, surface, withHeader=True,component=None, format='MS_TSES_ASCII', **kw)\n
##         Required Arguments:\n
##         filename --- path to the output file without an extension two files will be created filename.face and a filename.vert\n 
##         mol      --- Protein associated to the surface\n
##         surface  --- surface name\n

##         Optional Arguments:\n
##         withHeader --- True Boolean flag to specify whether or not to write the headers in the .face and the .vert files\n
##         component  --- msms component to save by default None\n
##         format     --- format in which the surface will be saved. It can be,\n        
##         MS_TSES_ASCII: Triangulated surface in ASCII format.\n
##         MS_ASES_ASCII: Analytical surface in ASCII format.This is a discrete representation of the analytical model.MS_TSES_ASCII_AVS: Triangulated surface in ASCII with AVS header\n
##         MS_ASES_ASCII_AVS: Analytical surface in ASCII format with AVS header\n
        
##         """
##         assert isinstance(filename, str)
##         assert isinstance( molName, str)
##         assert isinstance(surfName, str)
##         assert withHeader in [True, False, 1, 0]
##         assert format in self.formats
##         kw = {}
##         kw['withHeader'] = withHeader
##         kw['component'] = component
##         kw['format'] = format
##         return (filename, molName, surfName), kw
        

## class ReadMSMS(MVCommand):
##     """Command reads .face and .vert file, creates the msms surface and links it to the selection if can\n
##     Package : PmvApp\n
##     Module  : msmsCmds\n
##     Class   : ReadMSMS\n
##     Command name : readMSMS\n
##     nSynopsis :\n
##     None--->self.readMSMS(vertFilename, faceFilename, molName=None)\n
##     \nRequired Arguments :\n
##     vertFilename---name of the .vert file\n 
##     faceFilename---name of the .face file\n    
##     """

##     def __init__(self):
##         MVCommand.__init__(self)
##         self.msmsFromFile = {}
    
    
##     def doit(self, vertFilename, faceFilename, molName=None):
##         try:
##             vertFName = os.path.split(vertFilename)[1]
##             faceFName = os.path.split(faceFilename)[1]
##             vertName = os.path.splitext(vertFName)[0]
##             faceName = os.path.splitext(faceFName)[0]
##             assert vertName == faceName
##             msmsParser = MSMSParser()
##             self.msmsFromFile[vertName] = msmsParser
##             msmsParser.parse(vertFilename, faceFilename)
##             self.surf  = IndexedPolygons(vertName+'_msms', visible=1, 
##                                     pickableVertices=1, protected=True,)
##             if self.app().userpref['Color Boundaries for MSMS']['value'] == 'blurred':
##                 self.surf.Set(inheritSharpColorBoundaries=False, sharpColorBoundaries=False )
##             else:
##                 self.surf.Set(inheritSharpColorBoundaries=False, sharpColorBoundaries=True )
##             #self.surf.RenderMode(GL.GL_FILL, face=GL.GL_FRONT, redo=0)
##             #self.surf.Set(frontPolyMode=GL.GL_FILL, redo=0)
##             self.surf.Set(vertices=msmsParser.vertices, faces=msmsParser.faces,
##                           vnormals=msmsParser.normals, tagModified=False)
##             # The AppGUI should register an AddGeometryEvent listener
##             # that implements a method to add geometry to the application
##             # self.app().GUI.VIEWER.AddObject(self.surf)

##             from AppFramework.App import AddGeometryEvent
##             event = AddGeometryEvent(self.surf)
##             self.app().eventHandler.dispatchEvent(event)


##             if not molName is None:
##                 if not hasattr(self.app(), "bindGeomToMolecularFragment"):
##                     self.app().lazyLoad("displayCmds", commands=['bindGeomToMolecularFragment'],
##                                         package="PmvApp")
##                 self.app().bindGeomToMolecularFragment(self.surf, molName)

##                 # highlight selection
##                 surf = self.surf
##                 bindcmd = self.app().bindGeomToMolecularFragment
##                 selMols, selAtms = self.app().getNodesByMolecule(self.app().activeSelection.get(), Atom)
##                 lMolSelectedAtmsDict = dict( zip( selMols, selAtms) )
##                 #print self.name, "doit", surf.mol
##                 if lMolSelectedAtmsDict.has_key(surf.mol()):
##                     lSelectedAtoms = lMolSelectedAtmsDict[surf.mol()]
##                     if len(lSelectedAtoms) > 0:
##                         lAtomVerticesDict = bindcmd.data[surf.fullName]['atomVertices']
##                         highlight = [0] * len(surf.vertexSet.vertices)
##                         for lSelectedAtom in lSelectedAtoms:
##                             lVertexIndices = lAtomVerticesDict.get(lSelectedAtom, [])
##                             for lVertexIndex in lVertexIndices:
##                                 highlight[lVertexIndex] = 1
##                         surf.Set(highlight=highlight)
##             self.app()._executionReport.addSuccess('read surface from files %s, %s successfully'%(vertName, faceName ), obj=surf)
##         except:
##             msg = 'Error while reading surface for molecule %s'%vertFilename
##             self.app().errorMsg(sys.exc_info(), msg, obj=None)

##     def checkArguments(self, vertFilename, faceFilename, molName=None):
##         """None--->mv.readMSMS(vertFilename,faceFilename,molName=None, **kw)
##         """
##         assert os.path.exists(vertFilename)
##         assert os.path.exists(faceFilename)
##         if molName:
##             assert isinstance(molName, str)
##         kw = {}
##         kw['molName'] = molName
##         #kw['redraw'] = 1
##         return (vertFilename, faceFilename), kw
       


commandClassFromName = {
    'displayMSMS' : [DisplayMSMS, None],
    'undisplayMSMS' : [UndisplayMSMS, None],
    'computeMSMS' : [ComputeMSMS, None],
    #'computeSESAndSASArea' : [ComputeSESAndSASArea, None],
    #'readMSMS' : [ReadMSMS, None],
    #'saveMSMS' : [SaveMSMS, None],
    #'computeMSMSApprox' : [ComputeMSMSApprox, None],
    #'identifyBuriedVertices' : [IdentifyBuriedVertices, None],
    #'displayBuriedTriangles' : [DisplayBuriedTriangles, None],
    #'displayIntermolecularBuriedTriangles' : [DisplayIntermolecularBuriedTriangles, None],
    #'assignBuriedAreas' : [AssignBuriedAreas, None],
    
}

def initModule(viewer, gui=True):
    for cmdName, values in commandClassFromName.items():
        cmdClass, guiInstance = values
        viewer.addCommand(cmdClass(), cmdName, guiInstance)


