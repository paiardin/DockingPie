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

########################################################################
#
# Authors: Michel F. SANNER
#
# Copyright: Michel Sanner and TSRI, 2016
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/PmvApp/labelCmds.py,v 1.1.4.1 2017/07/13 20:52:52 annao Exp $
#
# $Id: labelCmds.py,v 1.1.4.1 2017/07/13 20:52:52 annao Exp $
#

"""
This Module implements commands to label the current selection different ways.
For example:
    by properties  of the current selection
"""

import numpy

from mglutil.util.colorUtil import ToHEX, TkColor

from PmvApp.Pmv import MVCommand

from MolKit2.molecule import Molecule, Atom,  MoleculeSet, \
     MOLECULE_RESIDUE_LABEL, MOLECULE_ATOM_LABEL
#from MolKit2.protein import Protein, Residue, Chain, ResidueSet, ChainSet, ProteinSet
from MolKit2.molecule import Residue, ResidueSet, Chain,  ChainSet
from DejaVu2.glfLabels import GlfLabels

from .displayCmds import DisplayCommand

class LabelBase(DisplayCommand):

    """The labelAtoms command displays a label for the specified set of atoms.
the labels can be passed explicitly to the command or extracted from atoms by
specifying an atomic property to be used as a label.

    Synopsis:
        None <- labelAtoms(atoms, labels=None, propName='name', font={}, formatString=None) 

    arguments:
        atoms  --- set of atoms

   Optional Arguments:    

        labels     --- None a list of labels for each atom (trumps propname)
        
        propname   --- attribute to use for label (default 'name')

        font       --- one of the following glfFonts fonts: (default 'arial1')
                       arial1, courier1, crystal1, techno0, techno1, times_new1,
                       aksent1, alpine1, broadway1, chicago1, compact1, cricket1,
                       garamond1, gothic1, penta1, present_script1

        formatString --- format string for label
    """
    
    _argNames = ['label', 'propName', 'font', 'formatString']

    _fonts = ['arial1', 'courier1', 'crystal1', 'techno0', 'techno1', 'times_new1',
              'aksent1', 'alpine1', 'broadway1', 'chicago1', 'compact1', 'cricket1',
              'garamond1', 'gothic1', 'penta1', 'present_script1']
    
    def _initializeMolForCmd(self, mol):
        """Adds the atomLabels geometry and the atomLabels atom group to 
        the object's geometry container

        adds per atom datafields:
            colorsIndices_atomLabels
            atomLabels or residuesLabels
        adds
        molecule._colors[self._atomSetName] = [colors]
        #molecule._renderingProp[self._atomSetName] = [ labels ]
        """
        if mol._colors.get(self._atomSetName, None) is None:
            # per atom properties used for lines
            # initialize with line colors
            mol._colors[self._atomSetName] = numpy.array( [[1.,1.,1.,1.]], 'f')

        if mol._ag.getData('colorsIndices_%s'%self._atomSetName) is None:
            mol._ag.setData('colorsIndices_%s'%self._atomSetName, numpy.zeros(len(mol._ag), 'i'))

        if mol._ag.getData(self._atomSetName) is None:
            mol._ag.setData(self._atomSetName, numpy.array([""]*len(mol._ag)))

        if mol._renderingProp.get(self._atomSetName, None) is None:
            mol._renderingProp[self._atomSetName] = {
                'font':'arial1', 'formatString':None,
                'linewidth':2, # line width
                'pointwidth':2, # line width
                'frontrendering':'fill', # frontPolyMode
                'backrendering': 'line', #backPolyMode
                'culling':'none',# 'front', 'front_and_back', 'none'
                'shading': 'smooth', # 'flat', 'inherit'
                'backfacecolor': numpy.array([[ 0.,  0.,  0.]], 'f'),
                'polyface': 'front' # 'back', 'frontandback', 'none'
                }
        gc = mol.geomContainer
        if not gc.geoms.has_key(self._atomSetName):
            g = GlfLabels(self._atomSetName, fontStyle='solid3d',
                          fontTranslation=(0,0,.1), # so the label is less hidden by the structure
                          fontScales = self.fontSize, # different size for atom, res, chain, mol
                          visible=0, pickable=0,
                          inheritMaterial=False,
                          )
            gc.addGeom(g, parent=gc.masterGeom, redo=0)
            g.vertexSet.vertices.array = gc.allCoords
            #g._hasAllVertices = True
            g.applyStrokes()

    def doit(self, molSel, **kw):
        sel = molSel.select('not deleted')
        mol = sel.getAtomGroup().getMolecule()

        self.initialize(mol)

        #setup undo
        self.updateModelGlobals(mol, **kw)
        self.updateModel(sel, **kw)
        self.refreshDisplay(mol)

    def updateModelGlobals(self, mol, font=None, formatString=None,
                           frontrendering=None, backrendering=None,
                           linewidth=None, pointwidth=None, 
                           shading=None, culling=None, backfacecolor=None,
                           polyface=None, **kw):
        self.initialize(mol)

        if font is not None:
            assert isinstance(font, dict) or font in self._fonts
            mol._renderingProp[self._atomSetName]['font'] = font

        if formatString is not None:
            mol._renderingProp[self._atomSetName]['formatString'] = formatString

        if frontrendering is not None:
            assert frontrendering in self._renderToDejaVu.keys()
            mol._renderingProp[self._atomSetName]['frontrendering'] = self._renderToDejaVu[frontrendering]
        if frontrendering is not None:
            assert frontrendering in self._renderToDejaVu.keys()
            mol._renderingProp[self._atomSetName]['frontrendering'] = self._renderToDejaVu[frontrendering]
        if backrendering is not None:
            assert backrendering in self._renderToDejaVu.keys()
            mol._renderingProp[self._atomSetName]['backrendering'] = self._renderToDejaVu[backrendering]
        if linewidth is not None:
            assert isinstance(linewidth, int) and linewidth >=1
            mol._renderingProp[self._atomSetName]['linewidth'] = linewidth
        if pointwidth is not None:
            assert isinstance(pointwidth, int) and pointwidth >=1
            mol._renderingProp[self._atomSetName]['pointwidth'] = pointwidth
        if shading is not None:
            assert isinstance(shading, str)
            mol._renderingProp[self._atomSetName]['shading'] = shading
        if culling is not None:
            assert isinstance(culling, str)
            mol._renderingProp[self._atomSetName]['culling'] = culling
        if  backfacecolor is not None:
            mol._renderingProp[self._atomSetName]['backfacecolor'] = backfacecolor  
        if polyface is not None:
            assert isinstance(polyFace, str)
            mol._renderingProp[self._atomSetName]['polyface'] = polyface

    def _updateAtomSet(self, mol, atoms):
        # update atoms set for which atomLabels are shown
        gc = mol.geomContainer
        ggeoms = gc.geoms
        _set = gc.atoms[self._atomSetName]
        
        if mol._multi == 'molecules':
            if mol._labelingBits & self._labBits: # if set
                mol._labelingBits -= self._labBits # remove the bit

        if len(atoms) == len(mol._ag):
            if self.negate:
                _set = mol.emptySelection()
            else:
                _set = atoms
                if mol._multi == 'molecules':
                    mol._labelingBits += self._labBits # set the bit
        else:
            ##if negate, remove current atms from displayed _set
            if self.negate:
                _set = _set - atoms
            else:
                _set = _set | atoms

        # at this point _set is what will still be displayed as atomlabels after this cmd
        gc.setAtomsForGeom(self._atomSetName, _set)

class LabelAtoms(LabelBase):

    def __init__(self):
        LabelBase.__init__(self)
        self._atomSetName = "atomLabels"
        self.fontSize = (0.25, 0.25, .01)
        self.placementModes = ['calpha', 'cbeta', "midAtom", 'gravityCenter']
        self._labBits = MOLECULE_ATOM_LABEL

    def initializeMolForCmd(self, mol):
        """Adds the atomLabels geometry and the atomLabels atom group to 
        the object's geometry container

        adds per atom datafields:
            colorsIndices_atomLabels
            atomLabels
        adds
        molecule._colors[self._atomSetName] = [colors]
        """
        if self.initializedFor.get(mol, False):
            return
        self.initializedFor[mol] = True
        self.app().unlabelAtoms.initializedFor[mol] = True

        self._initializeMolForCmd(mol)

    def getHBLabels(self, labels, ltype):
        newLab = numpy.array(['']*len(labels), ltype)
        # in lables array replace 1 with "A"(ACCEPTOR), 2 -"D" (DONOR),
        #  3 - "AD" (BOTH) and 0 - "".
        for i, l in [[1, "A"],[2, "D"],[3, "AD"]]:
            inds = numpy.where(labels==i)
            if len(inds):
                newLab[inds] = l
        return newLab

    def updateModel(self, atoms, labels=None, propName='name', **kw):
        # update the atomLabels display model for the atoms in sel
        #
        # propNames can be:
        #      None     : do not change the current label for these atoms
        #      string   : atomic property name
        # labels can be : (trumps propName)
        #      None     : use the current labels
        #      string   : use this label for all atoms in atoms
        #      [string] : each atom gets its own display label
        #
        mol = atoms.getAtomGroup().getMolecule()
        self.initialize(mol)
        indices = atoms.getIndices()

        if propName is not None and labels is None:
            if propName == "index":
                newPropLabels = indices
            else:
                newPropLabels = atoms.getData(propName)
            if newPropLabels is None:
                msg= "%s: specified nodes do not have property: %s"%(self.name, propName)
                self.app().warningMsg(msg)
                return
            else:
                if propName == "hbType":
                    newPropLabels = self.getHBLabels(newPropLabels, labelGeom.labels.dtype)

            formatString = mol._renderingProp[self._atomSetName]['formatString']
            if formatString is None:
                labels = [str(item) for item in newPropLabels]
            else:
                labels = [format%item for item in newPropLabels]

        # update the labels
        if labels is not None:
            length = max([len(x) for x in labels])
            allLabels = mol._ag.getData(self._atomSetName)
            allLabels = numpy.array(allLabels, '|S%d'%length)
            allLabels[indices] = labels
            mol._ag.setData(self._atomSetName, allLabels)
        self._updateAtomSet(mol, atoms)

    def refreshDisplay(self, mol):
        if not self.isInitialized(mol):
            return
        
        # find the geometry to use for atoms
        gc = mol.geomContainer
        geom = gc.geoms[self._atomSetName]
        atoms = mol.geomContainer.atoms.get(self._atomSetName, mol.emptySelection())

        if len(atoms)==0: # nothing is diplayed as CPK anymore for this molecule
            geom.Set(visible=0, tagModified=False) # hide the geom
        else:
            globProp = mol._renderingProp[self._atomSetName]
            font = globProp['font']
            if isinstance(font, dict):
                GlfLabels.Set(*(labelGeom, 1, 0), **font)
            else:
                geom.Set(font=font+".glf", redo=0)

            # compute per label offset depending on the geometry shown for this atom
            labOff = numpy.zeros( (len(mol._ag), 3), 'f')
            _seta = gc.atoms.get('sb', None)
            if _seta and len(_seta): # S&B
                labOff[_seta.getIndices(),2] = _seta.getData('sb_ballsRadius')+0.1

            if hasattr(mol, '_msmsData'):
                for name in mol._msmsData['msms'].keys():
                    _seta = gc.atoms.get('msms_%s'%name, None)
                    if _seta and len(_seta):
                        labOff[_seta.getIndices(),2] = _seta.getRadii()+0.1

            _seta = gc.atoms.get('cpk', None)
            if _seta and len(_seta):
                gcpk = gc.geoms['cpk']
                labOff[_seta.getIndices(),2] = _seta.getData('cpk_radius')+0.1

            indices = atoms.getIndices()
            labels = mol._ag.getData(self._atomSetName)
            colorInds = mol._ag.getData('colorsIndices_%s'%self._atomSetName)
            colors = mol._colors[self._atomSetName][colorInds.tolist()]
            
            geom.Set(vertices=atoms.getCoords(), #faces=[[x] for x in indices],
                     labels=labels[indices], labelTranslation=labOff[indices],
                     visible=1, tagModified=False, materials=colors[indices],
                     inheritLineWidth=False, lineWidth=globProp['linewidth'],
                     inheritPointWidth=False, pointWidth=globProp['pointwidth'],
                     frontPolyMode=globProp['frontrendering'], culling=globProp['culling'],
                     shading=globProp['shading'], backPolyMode=globProp['backrendering'],
                     polyFace=globProp['polyface'])

class UnlabelAtoms(LabelAtoms):
    """The UnlabelAtoms command removes the atom labels for the specified set ot atoms

    Synopsis:
        None <- unlabelAtoms(atoms)

    arguments:
        atoms  --- set of atoms

    Package : PmvApp
    Module  : labelCmds
    Class   : UnlabelAtoms
    Command : unlabelAtoms
    
    Examples:
    """
    def __init__(self):
        LabelAtoms.__init__(self)
        self.negate = True

    #def initializeMolForCmd(self, mol):
    #    if not self.isInitialized(mol):
    #        return
            #raise RuntimeError("unlabelAtoms called before labelAtoms for moelcule %s"%mol.name)

class LabelResidues(LabelBase):

    _argNames = ['label', 'font', 'formatString', 'placement']

    def __init__(self):
        LabelBase.__init__(self)
        self._atomSetName = "residueLabels"
        self.fontSize = (0.25, 0.25, .01)
        self.placementModes = ['calpha', 'cbeta', "midAtom", 'gravityCenter']
        self._labBits = MOLECULE_RESIDUE_LABEL

    def initializeMolForCmd(self, mol):
        if self.initializedFor.get(mol, False):
            return
        self.initializedFor[mol] = True
        self.app().unlabelResidues.initializedFor[mol] = True

        self._initializeMolForCmd(mol)
        mol._renderingProp[self._atomSetName]['placement'] = 'calpha'
        
    def updateModelGlobals(self, mol, placement=None, **kw):
        if placement is not None:
            assert placement in self.placementModes
            mol._renderingProp[self._atomSetName]['placement'] = placement
        
        LabelBase.updateModelGlobals(self, mol, **kw)
        
    def updateModel(self, atoms, labels=None, **kw):
        mol = atoms.getAtomGroup().getMolecule()
        self.initialize(mol)
        ## gc = mol.geomContainer

        ## if mol._multi == 'molecules':
        ##     if mol._labelingBits & mol._RESIDUE_LABEL:
        ##         mol._labelingBits -= mol._RESIDUE_LABEL

        ## sel = atoms.select('not deleted')
        ## # set rendering bit for entire molecule
        ## if len(sel) == len(mol._ag):
        ##     if self.negate:
        ##         _set = mol.emptySelection()
        ##     else:
        ##         _set = sel
        ##         if mol._multi == 'molecules':
        ##             mol._labelingBits += mol._RESIDUE_LABEL

        ## if self.negate:
        ##     _set = gc.atoms[self._atomSetName] - sel
        ## else:   
        ##     _set = gc.atoms[self._atomSetName] | sel
        ## gc.setAtomsForGeom(self._atomSetName, _set)

        if labels is None:
            buildLabels = True
            labels = []
        else:
            buildLabels = False

        maxlen = 0
        for res in atoms.getHierView().iterResidues():
            if self.negate:
                lab = ''
                labels.extend([lab]*res.numAtoms())
                maxlen = max(maxlen, 1)
            else:
                if buildLabels:
                    lab = '%c:%s%d'%(res.getChid(), res.getResname(), res.getResnum())
                    labels.extend([lab]*res.numAtoms())
                    maxlen = max(maxlen, len(lab))

        allLabels = mol._ag.getData(self._atomSetName)
        allLabels = numpy.array(allLabels, '|S%s'%maxlen)
        allLabels[atoms.getIndices()] = labels
        mol._ag.setData(self._atomSetName, allLabels)
        
        self._updateAtomSet(mol, atoms)

    def getOffset(self, mol, gc, atom):
        off = 0.0
        _seta = gc.atoms.get('sb', None)
        if _seta and len(set(_seta._indices) & set(atom.getIndices())):
            off = atom.getData('sb_ballsRadius')+0.1
        if hasattr(mol, '_msmsData'):
            for name in mol._msmsData['msms'].keys():
                _seta = gc.atoms.get('msms_%s'%name, None)
                if _seta and len(set(_seta._indices) & set(atom.getIndices())):
                    off = atom.getRadii()[0]+0.1
        _seta = gc.atoms.get('cpk', None)
        if _seta and len(set(_seta._indices) & set(atom.getIndices())):
            off = atom.getData('cpk_radius')+0.1

        return (0,0,off)
    
    def refreshDisplay(self, mol):
        if not self.isInitialized(mol):
            return

        # find the geometry to use for atoms
        gc = mol.geomContainer
        geom = gc.geoms[self._atomSetName]
        atoms = mol.geomContainer.atoms.get(self._atomSetName, mol.emptySelection())

        globProp = mol._renderingProp[self._atomSetName]
        font = globProp['font']
        if isinstance(font, dict):
            GlfLabels.Set(*(labelGeom, 1, 0), **font)
        else:
            geom.Set(font=font+".glf", redo=0)

        # build the placement in refreshDisplay so we donlt have to store the offset
        placement = mol._renderingProp[self._atomSetName]['placement']
        # location of the label
        numRes = len(numpy.unique(atoms.getResindices()))
        coords = numpy.zeros((numRes, 3),  "f")
        i = 0
        colorInds = mol._ag.getData('colorsIndices_residueLabels')
        colors = mol._colors[self._atomSetName][colorInds.tolist()]
        col = []
        labels = []
        allLabels = mol._ag.getData(self._atomSetName)
        if len(atoms):
            for res in atoms.getHierView().iterResidues():
                # use the color of the first atom in the residue
                labels.append(allLabels[res.getIndices()[0]])
                col.append(colors[res.getIndices()[0]])
                if placement=='gravityCenter':
                    coords[i] = numpy.sum(res.getCoords(), 0)/len(res)

                elif placement=='midAtom':
                    coords[i] = res.getCoords()[len(res)/2]

                elif placement=='calpha':
                    ca = res.select('calpha or element P')
                    if ca:
                        coords[i] = ca.getCoords()[0] + self.getOffset(mol, gc, ca)
                    else:
                        coords[i] = numpy.sum(res.getCoords(), 0)/len(res)
                elif placement=='cbeta':
                    cb = res.select('name CB and protein')
                    if cb:
                        coords[i] = cb.getCoords()[0] + self.getOffset(mol, gc, cb)
                    else:
                        ca = res.select('calpha')
                        if ca:
                            coords[i] = ca.getCoords()[0] + self.getOffset(mol, gc, ca)
                        else:
                            coords[i] = numpy.sum(res.getCoords(), 0)/len(res)
                i+=1

        if len(labels)==0: # nothing is diplayed as CPK anymore for this molecule
            geom.Set(visible=0, tagModified=False) # hide the geom
        else:
            geom.Set(vertices=coords, labels=labels,
                     visible=1, tagModified=False, materials=col,
                     inheritLineWidth=False, lineWidth=globProp['linewidth'],
                     inheritPointWidth=False, pointWidth=globProp['pointwidth'],
                     frontPolyMode=globProp['frontrendering'], culling=globProp['culling'],
                     shading=globProp['shading'], backPolyMode=globProp['backrendering'],
                     polyFace=globProp['polyface'])

class UnlabelResidues(LabelResidues):
    """The UnlabelResidues command removes the residue labels for the specified set ot atoms

    Synopsis:
        None <- unlabelResidues(atoms)

    arguments:
        atoms  --- set of atoms

    Package : PmvApp
    Module  : labelCmds
    Class   : UnlabelResidues
    Command : unlabelResidues
    
    Examples:
    """
    def __init__(self):
        LabelResidues.__init__(self)
        self.negate = True

    #def initializeMolForCmd(self, mol):
    #    if not self.isInitialized(mol):
    #        LabelResidues.initializeMolForCmd(mol)
            
            #raise RuntimeError("unlabelResidues called before labelResidues for moelcule %s"%mol.name)

commandClassFromName = {
    'labelAtoms' : [LabelAtoms,  None],
    'labelResidues' : [LabelResidues,  None],
    'unlabelAtoms' : [UnlabelAtoms,  None],
    'unlabelResidues' : [UnlabelResidues,  None],
}

def initModule(viewer, gui=True):
    for cmdName, values in commandClassFromName.items():
        cmdClass, guiInstance = values
        if gui:
            viewer.addCommand(cmdClass(), cmdName, guiInstance)
        else:
            viewer.addCommand(cmdClass(), cmdName, None)

