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
# Author: Michel F. SANNER, Anna Omelchenko
#
# Copyright: M. Sanner TSRI 2016
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/PmvApp/displayCmds.py,v 1.14.4.3 2017/10/25 23:17:47 annao Exp $
#
# $Id: displayCmds.py,v 1.14.4.3 2017/10/25 23:17:47 annao Exp $
#

import numpy, sys
from math import sqrt
from time import time

from PmvApp.Pmv import MVCommand, AddAtomsEvent, MoveAtomsEvent
from PmvApp.Pmv import RefreshDisplayEvent, ShowMoleculesEvent

from DejaVu2.Geom import Geom
from DejaVu2.Cylinders import Cylinders
from DejaVu2.IndexedPolylines import IndexedPolylines
from DejaVu2.IndexedGeom import IndexedGeom
from DejaVu2.Spheres import Spheres
from DejaVu2.Points import Points, CrossSet
from DejaVu2.Arcs3D import Fan3D
from DejaVu2.stipples import computeHalfBonds, stippleLines

from MolKit2.molecule import Atom, Molecule, MoleculeSet, Residue, Chain
from MolKit2.selection import Selection, SelectionSet

from opengltk.OpenGL import GL

##############################################################################
##
## DISPLAY  commands base class 
##
##############################################################################

class DisplayCommand(MVCommand):
    """The DisplayCommand class is the base class from which command displaying graphical representations of molecules are derived.

    This class implements the the following methods:

    class variables:
        _argNames: list of string providing the list of reognized keyword arguments for the command's doit() method

    attribute:
        _atomSetName: name fo the atom set created by this command for each molecule in the molecule's GeomContainer.atoms

    Model: the data associated with a molecule get extended the first time a 
           given display command is invoked on a molecule. This happens in the
           initializeMolForCmd(mol) method. In this method one or more DejaVu2
           geometry objects are created and added to the molecule
           GeometryContainer.geoms dictionary. An atom set called
           cmd._atomSetName is also added the the molecule's
           GeometryContainer.atoms dictionary. This atom set holds the list of
           atoms from that mol (as a ProDy selection) for which this
           representation is currently visible.
           The molecule's AtomGroup (mol._ag) is also potentially extended
           with new fields (e.g. colorIndices). These fields are automatically
           save into session through the prody saving mechanisms for atom
           groups.
           The molecule's _renderingProp dictionary is also be extended with
           a dictionnary of global rendering attributes (i.e. affecting the
           rendering of all atoms in that molecule) e.g. linewidth,
           pointwidth, solid vs mesh etc. 

    Methods:

        initialize(mol): initialize the molecule for thsi command if needed.
        
        refreshDisplay_cb(event): a method to call the command's refreshDisplay method 
                                  when a RefreshDisplayEvent is dispatched

        onAddCmdToApp(): this cmethod resisters the callback for RefreshDisplayEvent events
    
        checkArguments(): this method is call before the command's doit() method is invoked
                          and will raise a RuntimeError if a provided keyword argument is 
                          not recognized. The list of know keyword is the _argNames class 
                          variable defined by each DisplayCommand subclass.
                      
    Virtual Methods:

        initializeMolForCmd(mol):
            create the atomset, and DejaVu geometry objects and extend molecule and atom
            group to store information needed by refreshDisplay to render the molecule

        updateModelGlobals(mol, **kw):
            Set the model's global parameters, i.e. affecting all atoms rendered by this
            command for a given molecule
            This method has to initialize the molecule for the command if it has not
            been done yet. It should NOT call refreshDisplay.
        
         updateModel(self, molsel, **kw):
            Update per-atom parameters of the model
            This method has to initialize the molecule for the command if it has not
            been done yet. It should NOT call refreshDisplay.

         doit(self, molsel, **kw):
             call updateModelGlobals(mol) and updateModel(selection), then call
             refreshDisplay(mol) to update the DejaVu geometry and trigger a Redraw
             This method has to initialize the molecule for the command if it has not
             been done yet.

         refreshDisplay(self, mol):
             Update the DejaVu geometry and trigger a Redraw.
             this method should return if the molecule is not initialized for thsi command.
             
    extends each molecules with attributes (called the model) parametrizing the representation. Some parameters are global, meaning they affect
    
    \nPackage : PmvApp
    \nModule  : displayCmds
    \nClass   : DisplayCommand
    """

    _renderToDejaVu = {'solid':'fill', 'mesh':'line', 'points':'point', 'sameAsFront':'sameAsFront'}
    _renderToPmv = {'fill':'solid', 'line':'mesh', 'point':'points', 'sameAsFront':'sameAsFront'}

    def __init__(self):
        MVCommand.__init__(self)
        # set to True or False by sublasses implementing display vs
        # undisplay commands
        self.negate = False

        # this is the name of the atom set in mol.geomContainer
        # for atoms displayed for a given display command
        self._atomSetName = "NoName"
        
    def onAddCmdToApp(self):
        MVCommand.onAddCmdToApp(self)
        evh = self.app().eventHandler
        evh.registerListener(RefreshDisplayEvent, self.refreshDisplay_cb)
        evh.registerListener(AddAtomsEvent, self.AddAtomsEvent_cb)
        evh.registerListener(MoveAtomsEvent, self.MoveAtomsEvent_cb)
        
    def isInitialized(self, mol):
        return self.initializedFor.get(mol, False)
    
    def initialize(self,mol):
        # call self.initializedFor is needed
        isInitialized = self.initializedFor.get(mol, False)
        if not isInitialized:
            self.initializeMolForCmd(mol) 

    def checkArguments(self, nodes, **kw):
        for name in kw.keys():
            if name not in self._argNames:
                raise RuntimeError("%s: unrecognized keyword argument '%s', valid names are %s"%(
                    self.name, name, str(self._argNames)))
        return (nodes,), kw

    def refreshDisplay_cb(self, event):
        self.refreshDisplay(event.molecule)

    def AddAtomsEvent_handler(self, mol, atoms):
        pass

    def MoveAtomsEvent_handler(self, mol, atoms):
        self.refreshDisplay(mol)
        
    def AddAtomsEvent_cb(self, event):
        self.AddAtomsEvent_handler(event.molecule, event.atoms)

    def MoveAtomsEvent_cb(self, event):
        self.MoveAtomsEvent_handler(event.molecule, event.atoms)

    def firstArgStringToMol(self, selection):
        if not isinstance (selection, str):
            return selection
        return self.app().Mols.selectMolKit(selection)

#############################################################################
##
## DISPLAY LINES
##
#############################################################################

class DisplayLines(DisplayCommand):
    """The displayLines command displays atomic bonds as lines.

    Synopsis:
        None <- displayLines(atoms, stippled=False,
                             displayBondOrder=True,
                             linewidth=2
                             aromaticLinewidth=1,
                             doubleBondSep=0.15,
                             tripleBondSep=0.4,
                             stippleLength=0.2,
                             stippleSpace=0.2)

    arguments:
        atoms    --- set of atoms

        stippled --- boolean flag specifying whether or not to display the
                     bond as stippled lines (default = False)

        # global parameters: (i.e. affecting all bonds in the line representation of a molecule)

        displayBondOrder --- boolean. When True, double, triple, and aromatic bonds are displayed

        linewidth    --- integer > 1 specifying the width of lines for drawing spheres as a mesh. (Default:2)
        aromaticLinewidth --- (default:1)
        doubleBondSep --- float specifying the the distance between the 2 lines representing a double bond (default:0.12)
        tripleBondSep --- float specifying the distance between the 2 lines flanking the bond in the representation of a triple bond (default:0.2)
        stippleLength --- float specifying the length of the stipple dash (default:0.2)
        stippleSpace  --- float specifying the length of the space between stipples (default:0.2)

    Global model parameters can be set by calling:

             displayLines.updateModelGlobals(mol, displayBondOrder=None,
                                                  linewidth=None,
                                                  aromaticLinewidth=None,
                                                  doubleBondSep=None,
                                                  tripleBondSep=None,
                                                  stippleLength=None,
                                                  stippleSpace=None,
                                                  sharpColorBoundaries=None)

         This call will not trigger a refresh of the display. To see the effect a call to
         pmv.displayLines.refreshDisplay(mol) is needed

    Model parameters can be set by calling:

              displayLines.updateModel(sel, stippled=False)

         This call will not trigger a refresh of the display. To see the effect a call to
         pmv.displayLines.refreshDisplay(mol) is needed

    the list of provided arguments is checked when displaylines is called. If an name is not recognized
    a RuntimeError exception is Thrown

    Package : PmvApp
    Module  : displayCmds
    Class   : DisplayLines
    Command : displayLines

    Examples:
        mol = Pmv.Mols[0]

        ## display lines for then entire molecule with vdw radii
        pmv.displayLines(mol)

        ## hide all lines for mol
        pmv.undisplayLines(mol)
    """

    _argNames = ['stippled', 'displayBondOrder', 'linewidth',
                 'aromaticLinewidth', 'doubleBondSep', 'tripleBondSep',
                 'stippleLength', 'stippleSpace', 'sharpColorBoundaries']

    def __init__(self):
        DisplayCommand.__init__(self)
        self._atomSetName = "lines"

    def AddAtomsEvent_handler(self, mol, atoms):
        del self.initializedFor[mol]
        lineAtoms = mol.geomContainer.atoms.get('lines', None)
        if lineAtoms:
            self.app().displayLines(atoms)

    ## display and undisplay lines are now builtin
    ## def onAddCmdToApp(self):
    ##     DisplayCommand.onAddCmdToApp(self)
    ##     if not self.app().commands.has_key('undisplayLines'):
    ##         self.app().lazyLoad('displayCmds', commands=['undisplayLines'], package='PmvApp')
    ##     if self.app().undisplayLines.isLoader():
    ##         self.app().undisplayLines.loadCommand()

    def initializeMolForCmd(self, mol):
        """
        Creates the lines, the points, the dotted lines and the bond order
        lines geometries and add these new geometries to the
        geomContainer.geoms of the new molecule.

        adds per atom datafields:
             colorsIndices_lines
             stippled

        ## the displayCmd.doit method modifies the model (i.e. atom properties)
        ## and some attribute in the molecule. The refreshDisplay method will update
        ## the geometry object rendering the molecule
        ##
        ## model attributes for rendering lines:
        ##    per atom:  ag.colorsIndices[geom] where geom can be noBond, singleBonds,
        ##                             'doubleBonds', 'tripleBonds' or 'aromaticBonds'
        ##               ag.stipples (flag)
        ##    per molecule:
        ##               mol._colors[geom]
        ##

        """
        if self.initializedFor.get(mol, False):
            return
        self.initializedFor[mol] = True
        self.app().undisplayLines.initializedFor[mol] = True

        geomC = mol.geomContainer

        if mol._colors.get('lines', None) is None:
            # per atom properties used for lines
            mol._colors['lines'] = numpy.array( [[1.,1.,1.,1.]], 'f')
        l = len(mol._ag)
        if mol._ag.getData('colorsIndices_lines') is None:
            mol._ag.setData('colorsIndices_lines', numpy.zeros(len(mol._ag), 'i'))

        if mol._ag.getFlags('stippled') is None:
            mol._ag.setFlags('stippled', [False]*len(mol._ag))

        if mol._renderingProp.get('lines', None) is None:
            mol._renderingProp['lines'] = {'linewidth':2,
                                           'displayBondOrder':True,
                                           'aromaticLinewidth':1,
                                           'doubleBondSep': 0.12,
                                           'tripleBondSep': 0.2,
                                           'stippleLength':0.2,
                                           'stippleSpace':0.2,
                                           'sharpColorBoundaries':True}
        # lines representation needs 4 different geometries need to create
        # a master geometry Geom.
        if not geomC.geoms.has_key('lines'):
            wire = Geom(self._atomSetName,
                        shape=(0,0),
                        inheritLighting=False,
                        inheritLineWidth=False,
                        inheritSharpColorBoundaries=False,
                        lighting=False,
                        protected=True, pickable=0)
            geomC.addGeom(wire, parent=geomC.masterGeom, redo=0,
                          displayCmd=self, undisplayCmd=self.app().undisplayLines) 
            self.managedGeometries.append(wire)

            #geomC.atomPropToVertices["bonded"] = self.atomPropToVertices
            #geomC.atomPropToVertices["nobnds"] = self.atomPropToVertices
            #geomC.atomPropToVertices["bondorder"] = self.atomPropToVertices

            # nobnds : Non Bonded atoms are represented by Points
            p = Points( "noBond", shape=(0,3), visible=0,
                        inheritMaterial=0, materials=[[1.,1.,1.,1.]],
                        inheritLineWidth=0, lineWidth=2.,
                        protected=True,
                        disableStencil=True,
                        transparent=True
                        )
            p.vertexSet.vertices.array = geomC.allCoords

            # the _hasAllVertices flag is used in Pmv.handleAddAtoms
            # to update geom.vertexSet when atoms are added
            p._hasAllVertices = True

            geomC.addGeom( p, parent=wire, redo=0)#, displayCmd=self,
                          # undisplayCmd=self.app().undisplayLines)
            # adding to managed geom will ensure that the geometry is deleted
            # when the molecule is removed from the App
            self.managedGeometries.append(p)

            # lines: single bonds are represented by IndexedPolyLines
            l = IndexedPolylines("singleBonds", visible=0, inheritMaterial=0,
                                 protected=True, inheritSharpColorBoundaries=False)
            l.vertexSet.vertices.array = geomC.allCoords
            l._hasAllVertices = True
            geomC.addGeom(l, parent=wire, redo=0, makeAtomSet=False,
                          displayCmd=self, undisplayCmd=self.app().undisplayLines)
            self.managedGeometries.append(l)
            # add geom for stippled single bonds 
            ls = IndexedPolylines("singleBondsStippled", visible=0,
                                  inheritMaterial=0,
                                  inheritLineWidth=False,
                                  inheritSharpColorBoundaries=False,
                                  protected=True, pickable=1)
            ls.pickingVertices = geomC.allCoords
            ls._hasAllVertices = True
            geomC.addGeom(ls, parent=wire, redo=0, makeAtomSet=False,
                          displayCmd=self, undisplayCmd=self.app().undisplayLines)
            self.managedGeometries.append(ls)

            # lines: double bonds are represented by IndexedPolyLines
            b = IndexedPolylines('doubleBonds', visible=0,
                                 inheritMaterial=0, protected=True,
                                 inheritSharpColorBoundaries=False,
                                 inheritLineWidth=False)
            geomC.addGeom( b, parent=wire, redo=0, makeAtomSet=False,
                           displayCmd=self, undisplayCmd=self.app().undisplayLines)
            geomC.geomPickToAtoms[b.name] = self.pickedDoubleBondAtoms
            self.managedGeometries.append(b)
            # add geom for stippled single bonds 
            bs = IndexedPolylines('doubleBondsStippled', visible=0,
                                  inheritMaterial=0, protected=True,
                                  inheritLineWidth=False, inheritSharpColorBoundaries=False)
            #bs.pickingVertices = geomC.allCoords
            #bs._hasAllVertices = True
            geomC.addGeom(bs, parent=wire, redo=0, makeAtomSet=False,
                          displayCmd=self, undisplayCmd=self.app().undisplayLines)
            geomC.geomPickToAtoms[bs.name] = self.pickedDoubleBondAtoms
            self.managedGeometries.append(bs)

            # lines: triple bonds are represented by IndexedPolyLines
            b = IndexedPolylines('tripleBonds', visible=0, pickable=0,
                                 inheritLineWidth=False,
                                 inheritSharpColorBoundaries=False,
                                 inheritMaterial=0, protected=True)
            geomC.addGeom( b, parent=wire, redo=0, makeAtomSet=False,
                           displayCmd=self, undisplayCmd=self.app().undisplayLines)
            self.managedGeometries.append(b)
            # add geom for stippled single bonds 
            bs = IndexedPolylines('tripleBondsStippled', visible=0,
                                  inheritMaterial=0, protected=True,
                                  inheritLineWidth=False,
                                  inheritSharpColorBoundaries=False,
                                  pickable=0,)
            geomC.addGeom(bs, parent=wire, redo=0, makeAtomSet=False,
                          displayCmd=self,
                          undisplayCmd=self.app().undisplayLines)
            self.managedGeometries.append(bs)

            # lines: aromatic bonds are represented IndexedPolyLines
            b = Fan3D('aromaticBonds', fan=False, visible=0, pickable=0,
                      inheritMaterial=0, protected=True,
                      inheritStippleLines=False, stippleLines=True,
                      inheritLineWidth=False, lineWidth=1,
                      disableStencil=True, transparent=True,
                      inheritLighting=False, lighting=False,
                      inheritCulling=False, culling='none')
            geomC.addGeom( b, parent=wire, redo=0, makeAtomSet=False,
                           displayCmd=self, undisplayCmd=self.app().undisplayLines)
            self.managedGeometries.append(b)

    def atomPropToVertices(self, geom, atms, propName, propIndex=None):
        """Function called to compute the array of properties"""
        #print "atomPropToVertices:", geom, propName
        if atms is None or len(atms)==0 : return None
        prop = []
        #if propIndex == 'bondorder':
        #    mol = geom.mol()
        #    atms = mol.geomContainer.atoms['bondorder']

        propIndex = 'lines'
        for a in atms:
            d = getattr(a, propName)
            prop.append( d[propIndex] )
        return prop
        
    def doit(self, molSel, **kw):
        ##
        ##  displayLines uses 2 atomsSets 'noBond' and 'lines'
        ##

        sel = molSel.select('not deleted')
        mol = sel.getAtomGroup().getMolecule()
        self.initialize(mol)

        #setup undo
        # FIXME to undo stipple we neet to set sipple array
        self.app().pushUndoCmd(
            self.app().displayLines if self.negate else self.app().undisplayLines,
            (sel,), {'displayBondOrder': mol._renderingProp['lines']['displayBondOrder'],
                     #'stippled':stipple,
                     'linewidth': mol._renderingProp['lines']['linewidth']})
        self.updateModelGlobals(mol, **kw)
        self.updateModel(sel, **kw)
        self.refreshDisplay(mol)

    def updateModelGlobals(self, mol, displayBondOrder=None, linewidth=None,
                           aromaticLinewidth=None,
                           doubleBondSep=None,
                           tripleBondSep=None,
                           stippleLength=None,
                           stippleSpace=None,
                           sharpColorBoundaries=None, **kw):
        self.initialize(mol)

        if linewidth is not None:
            assert isinstance(linewidth, int) and linewidth >=1
            mol._renderingProp['lines']['linewidth'] = linewidth

        if aromaticLinewidth is not None:
            assert isinstance(aromaticLinewidth, int) and aromaticLinewidth >=1
            mol._renderingProp['lines']['aromaticLinewidth'] = aromaticLinewidth

        if doubleBondSep is not None:
            assert isinstance(doubleBondSep, float)
            sep = mol._renderingProp['lines']['doubleBondSep']
            if sep != doubleBondSep:
                mol._renderingProp['lines']['doubleBondSep'] = doubleBondSep
                v = mol.computeDoubleBondsVertices(mol.select().getBonds()[2], doubleBondSep)
                mol.geomContainer.geoms['doubleBonds'].vertexSet.vertices.array = v

        if tripleBondSep is not None:
            assert isinstance(tripleBondSep, float)
            sep = mol._renderingProp['lines']['tripleBondSep']
            if sep != tripleBondSep:
                mol._renderingProp['lines']['tripleBondSep'] = tripleBondSep
                v = mol.computeTripleBondsVertices(mol.select().getBonds()[3], tripleBondSep)
                mol.geomContainer.geoms['tripleBonds'].vertexSet.vertices.array = v

        if stippleLength is not None:
            assert isinstance(stippleLength, float) and stippleLength >0.
            mol._renderingProp['lines']['stippleLength'] = stippleLength

        if stippleSpace is not None:
            assert isinstance(stippleSpace, float) and  stippleSpace>0.
            mol._renderingProp['lines']['stippleSpace'] = stippleSpace

        if displayBondOrder is not None:
            assert isinstance(displayBondOrder, bool)
            mol._renderingProp['lines']['displayBondOrder'] = displayBondOrder

        if sharpColorBoundaries is not None:
            assert sharpColorBoundaries in [1, 0, False, True]
            mol._renderingProp['lines']['sharpColorBoundaries'] = sharpColorBoundaries

    def updateModel(self, molsel, stippled=None, **kw):
        mol = molsel.getAtomGroup().getMolecule()
        self.initialize(mol)
        gc = mol.geomContainer

        if mol._multi == 'molecules':
            if mol._renderingBits & mol._LINES:
                mol._renderingBits -= mol._LINES

        # set stippled flag for the selection
        if stippled is not None:
            mol._ag._flags['stippled'][molsel.getIndices()] = stippled

        sel = molsel.select('not deleted')
        # set rendering bit for entire molecule
        if len(sel) == len(mol.select()):
            if self.negate:
                _set = mol.emptySelection()
                gc.setAtomsForGeom('noBond',  _set)
            else:
                _set = sel
                if mol._multi == 'molecules':
                    mol._renderingBits += mol._LINES
                bonds = _set.getBonds()
                gc.setAtomsForGeom('noBond',  Selection(mol._ag, bonds[0], ''))
            gc.setAtomsForGeom('lines', _set)
        else:
            ## update atom sets
            # first unbonded atoms
            _set = gc.atoms['noBond'] # all atoms displayed as lines
            if self.negate:
                _set = _set - sel
            else:   
                _set = _set | sel
            bonds = _set.getBonds()
            gc.setAtomsForGeom('noBond', Selection(mol._ag, bonds[0], ''))

            # now bonded atoms
            _set = gc.atoms['lines'] # all atoms displayed as lines
            if self.negate:
                _set = _set - sel
            else:   
                _set = _set | sel
            gc.setAtomsForGeom('lines', _set)

    def refreshDisplay(self, mol):
        if not self.isInitialized(mol):
            return

        atomSets = mol.geomContainer.atoms
        atomsNoBond = atomSets.get('noBond', mol.emptySelection())
        gc = mol.geomContainer
        ggeoms = gc.geoms

        # get a list of all atom colors for lines
        colorInds = mol._ag.getData('colorsIndices_lines')
        colors = mol._colors['lines'][colorInds.tolist()]
        opacity = mol._ag.getData('opacity_lines')
        if opacity is not None and len(numpy.unique(opacity))==1:
            opacity = [opacity[0]]
        if len(atomsNoBond)==0:
            ggeoms['noBond'].Set(visible=False, tagModified=False)
        else:
            indices = atomsNoBond.getIndices()
            ggeoms['noBond'].Set(
                faces=[[x] for x in indices],
                pointWidth=mol._renderingProp['lines']['linewidth'],
                visible=1, tagModified=False,
                materials=colors, opacity=opacity, transparent="implicit",
                sharpColorBoundaries=mol._renderingProp['lines']['sharpColorBoundaries']     
)
        
        atomsWithBonds = atomSets.get('lines', mol.emptySelection())

        if len(atomsWithBonds)==0: # no bonds at all
            for i, name in enumerate(['singleBonds', 'doubleBonds',
                                      'tripleBonds', 'aromaticBonds']):
                ggeoms[name].Set(visible=False, tagModified=False)

        else: # we have bonds to display

            atomsBondsSolid = atomsWithBonds.select('not stippled')
            atomsBondsStippled = atomsWithBonds.select('stippled')
            allBonds = atomsWithBonds.getBonds()

            for atoms, solid in zip([atomsBondsSolid, atomsBondsStippled], (True, False)):
                if atoms is None or len(atoms) == 0:
                    if solid:
                        ggeoms['singleBonds'].Set(visible=False, tagModified=False)
                        ggeoms['doubleBonds'].Set(visible=False, tagModified=False)
                        ggeoms['tripleBonds'].Set(visible=False, tagModified=False)
                        ggeoms['aromaticBonds'].Set(visible=False, tagModified=False)
                    else:
                        ggeoms['singleBondsStippled'].Set(visible=False, tagModified=False)
                        ggeoms['doubleBondsStippled'].Set(visible=False, tagModified=False)
                        ggeoms['tripleBondsStippled'].Set(visible=False, tagModified=False)
                    continue

                # bonds[1] has single and aromatic, [2] has double and [3] triple bonds
                # bonds[4] has aromatic bonds only
                if not mol._renderingProp['lines']['displayBondOrder'] or mol._ag._bondOrder is None:
                    # draw all bonds as single bonds
                    singleBonds = allBonds[1]+allBonds[2]+allBonds[3]
                    allBonds[2] = []
                    allBonds[3] = []
                    allBonds[4] = []
                else:
                    singleBonds = allBonds[1]

                # segregate solid bonds and stippled bonds
                # bonds between atom that are solid and stippled will be solid
                sf = mol._ag._flags['stippled']
                if solid:
                    lbonds = []
                    for b in singleBonds:
                        if not sf[b[0]] or not sf[b[1]]:
                            lbonds.append(b)
                else:
                    lbonds = []
                    for b in singleBonds:
                        if sf[b[0]] and sf[b[1]]:
                            lbonds.append(b)
                singleBonds = lbonds
                
                if len(singleBonds) == 0: # no single bonds
                    if solid:
                        ggeoms['singleBonds'].Set(visible=False, tagModified=False)
                    else:
                        ggeoms['singleBondsStippled'].Set(visible=False, tagModified=False)
                else:
                    # process single bonds
                    if len(numpy.unique(colorInds))==1:
                        mat = [mol._colors['lines'][colorInds[0]]]
                    else:
                        mat = colors
                    if solid:
                        ggeoms['singleBonds'].Set(faces=numpy.array(singleBonds, 'int'),
                                                  lineWidth=mol._renderingProp['lines']['linewidth'],
                                                  inheritLineWidth=False,
                                                  materials=mat, visible=True, tagModified=False, opacity=opacity, transparent="implicit", sharpColorBoundaries=mol._renderingProp['lines']['sharpColorBoundaries'])
                    else:
                        vert = mol._ag.getCoords()
                        #faces = numpy.array(singleBonds, 'int')
                        v, f, c = computeHalfBonds(vert, singleBonds, colors)
                        v1, f1, cf1, cv1 = stippleLines(v, f, c,
                                                        segLen=mol._renderingProp['lines']['stippleLength'],
                                                        spaceLen=mol._renderingProp['lines']['stippleSpace'])
                        ggeoms['singleBondsStippled'].Set(
                            visible=True, vertices=v1, faces=f1, materials=cf1,
                            lineWidth=mol._renderingProp['lines']['linewidth'],
                            inheritLineWidth=False, tagModified=False,
                            sharpColorBoundaries=mol._renderingProp['lines']['sharpColorBoundaries'])
                        pfaces = numpy.array(singleBonds)
                        ggeoms['singleBondsStippled'].pickingFaces = pfaces
                    ##
                    ## double bonds
                    ##
                    bondList = allBonds[2] # double bonds
                    # if we do not display bond order this list has been set to emty above
                    doubleBonds = mol._bondOrderData['doubleBonds']
                    if len(bondList)==0 or doubleBonds is None or len(doubleBonds)==0:
                        if solid:
                            ggeoms['doubleBonds'].Set(visible=False, tagModified=False)
                        else:
                            ggeoms['doubleBondsStippled'].Set(visible=False, tagModified=False)
                    else:
                        lbonds = []
                        if solid:
                            for b in bondList:
                                if not sf[b[0]] or not sf[b[1]]:
                                    lbonds.append(b)
                        else:
                            lbonds = []
                            for b in bondList:
                                if sf[b[0]] and sf[b[1]]:
                                    lbonds.append(b)
                        bondList = lbonds
                        faces = []
                        
                        col = numpy.array([[1.,1.,1.,1.]]*len(ggeoms['doubleBonds'].vertexSet.vertices))
                        op = None
                        addOpacity = False
                        if opacity is not None:
                            if len(opacity) == 1 : 
                                op = [opacity[0]]
                            else:
                                op = numpy.ones(len(ggeoms['doubleBonds'].vertexSet.vertices)).astype('f')
                                addOpacity = True
                        # loop over bonds
                        for i1,i2 in bondList:
                            k = "%d %d" %(i1, i2)
                            if doubleBonds.has_key(k):
                                n = doubleBonds[k] * 4
                                faces.append( (n, n+1) )
                                faces.append( (n+2, n+3) )
                                ## col.append(colors[i1])
                                ## col.append(colors[i2])
                                ## col.append(colors[i1])
                                ## col.append(colors[i2])
                                col[n] = colors[i1]
                                col[n+1] = colors[i2]
                                col[n+2] = colors[i1]
                                col[n+3] = colors[i2]
                                if addOpacity:
                                    op[[n,n+1,n+2,n+3]] = [opacity[i1], opacity[i2], opacity[i1], opacity[i2]]
                        if solid:
                            ggeoms['doubleBonds'].Set(faces=faces, materials=col,
                                                      lineWidth=mol._renderingProp['lines']['linewidth'],
                                                      inheritLineWidth=False,
                                                      visible=True, tagModified=False,
                                                      opacity=op, transparent="implicit",
                                                      sharpColorBoundaries=mol._renderingProp['lines']['sharpColorBoundaries'])
                        else:
                            v0 = ggeoms['doubleBonds'].getVertices()
                            v, f, c = computeHalfBonds(v0, faces, col)
                            v1, f1, cf1, cv1 = stippleLines(v, f, c,
                                                        segLen=mol._renderingProp['lines']['stippleLength'],
                                                        spaceLen=mol._renderingProp['lines']['stippleSpace'])
                            ggeoms['doubleBondsStippled'].Set(
                                visible=True, vertices=v1, faces=f1, materials=cf1,
                                lineWidth=mol._renderingProp['lines']['linewidth'],
                                inheritLineWidth=False,
                                sharpColorBoundaries=mol._renderingProp['lines']['sharpColorBoundaries']
                                )
                            ggeoms['doubleBondsStippled'].pickingVertices = v0
                            ggeoms['doubleBondsStippled'].pickingFaces = numpy.array(faces)
                    ##
                    ## triple bonds
                    ##
                    bondList = allBonds[3] # triple bonds
                    tripleBonds = mol._bondOrderData['tripleBonds']
                    if len(bondList)==0 or tripleBonds is None or len(tripleBonds)==0:
                        if solid:
                            ggeoms['tripleBonds'].Set(visible=False, tagModified=False)
                        else:
                            ggeoms['tripleBondsStippled'].Set(visible=False, tagModified=False)
                    else:
                        lbonds = []
                        if solid:
                            for b in bondList:
                                if not sf[b[0]] or not sf[b[1]]:
                                    lbonds.append(b)
                        else:
                            lbonds = []
                            for b in bondList:
                                if sf[b[0]] and sf[b[1]]:
                                    lbonds.append(b)
                        bondList = lbonds
                        faces = []
                        col = numpy.array([[1.,1.,1.,1,]]*len(ggeoms['tripleBonds'].vertexSet.vertices))
                        # loop over bonds
                        op = None
                        if opacity is not None:
                            if len(opacity) == 1 : op = [opacity[0]]
                            else:
                                op = numpy.ones(len(ggeoms['tripleBonds'].vertexSet.vertices)).astype('f')
                        for i1,i2 in bondList:
                            k = "%d %d" %(i1, i2)
                            if tripleBonds.has_key(k):
                                n = tripleBonds[k] * 6
                                faces.append( (n, n+1) )
                                faces.append( (n+2, n+3) )
                                col[[n,n+1,n+2,n+3]] = [colors[i1], colors[i2],
                                                        colors[i1], colors[i2]]
                                if opacity is not None and len(opacity)>1:
                                    op[[n,n+1,n+2,n+3]] = [opacity[i1], opacity[i2], 
                                                           opacity[i1], opacity[i2]]
                        if solid:
                            ggeoms['tripleBonds'].Set(faces=faces, materials=col,
                                                      lineWidth=mol._renderingProp['lines']['linewidth'],
                                                      visible=True, tagModified=False,
                                                      opacity=op, transparent="implicit",
                                                      sharpColorBoundaries=mol._renderingProp['lines']['sharpColorBoundaries'])
                        else:
                            vert = mol._ag.getCoords()
                            v, f, c = computeHalfBonds(ggeoms['tripleBonds'].getVertices(), faces, col)
                            v1, f1, cf1, cv1 = stippleLines(v, f, c,
                                                        segLen=mol._renderingProp['lines']['stippleLength'],
                                                        spaceLen=mol._renderingProp['lines']['stippleSpace'])
                            ggeoms['tripleBondsStippled'].Set(
                                visible=True, vertices=v1, faces=f1, materials=cf1,
                                lineWidth=mol._renderingProp['lines']['linewidth'],
                                inheritLineWidth=False,
                                sharpColorBoundaries=mol._renderingProp['lines']['sharpColorBoundaries'])
                    ##
                    ## aromatic bonds
                    ##
                    bondList = allBonds[4] # aromatic bonds
                    aromaticArcs = mol._bondOrderData['aromaticArcs']
                    if len(bondList)==0 or aromaticArcs is None or len(aromaticArcs)==0:
                        if solid:
                            ggeoms['aromaticBonds'].Set(visible=False, tagModified=False)
                    else:
                        if solid:
                            geom = ggeoms['aromaticBonds']
                            faces = []
                            col = numpy.array([[1.,1.,1.,1.,]]*len(geom.vertexSet.vertices))
                            op = None
                            addOpacity = False
                            if opacity is not None:
                                if len(opacity) == 1 : op = [opacity[0]]
                                else:
                                    op = numpy.ones(len(geom.vertexSet.vertices)).astype('f')
                                    addOpacity = True
                            for index in numpy.unique(bondList):
                                if aromaticArcs.has_key(index):
                                    for ii in aromaticArcs[index]:
                                        faces.append( [ii[0]] )
                                        col[ii[0]] = colors[index]
                                        if addOpacity:
                                            op[ii[0]] = opacity[index]
                            geom.Set(faces=faces, materials=col,
                                     lineWidth=mol._renderingProp['lines']['aromaticLinewidth'],
                                     inheritLineWidth=False,
                                     visible=True, tagModified=False,
                                     opacity=op, transparent="implicit",
                                     sharpColorBoundaries=mol._renderingProp['lines']['sharpColorBoundaries'])

    def pickedDoubleBondAtoms(self, geom, vertInds):
        """
        convert the vertex indices in geom into atom indices
        """
        # get indices of atoms in mol._bondOrderData['doubleBonds'] which is a
        # dict 'i j':doubleBondNum for double bonds from which (i1, j1, i2, j2)
        #  vertices in the doubleBonds geometry were generated
        #     
        #     i1-------------j1
        #     i              j
        #     i2-------------j2
        #
        vertInds = numpy.array(vertInds)
        dbIndices = vertInds/4   # index into mol._bondOrderData['doubleBonds']
        dbAtomIndices = vertInds%4  # offset within mol._bondOrderData['doubleBonds'][index]
        atomIndices = []
        db = geom.mol()._bondOrderData['doubleBonds']
        dbNums = db.values()
        dbKeys = db.keys()
        for dbIndex, offset in zip(dbIndices, dbAtomIndices):
            ind = dbNums.index(dbIndex) # find the double bond index in the dict
            i, j = [int(n) for n in dbKeys[ind].split()] # split key 'i j'
            if offset%2==0:
                atomIndices.append(i) # first atom in double bond
            else:
                atomIndices.append(j) # second atom in double bond
        return numpy.unique(atomIndices)

class UndisplayLines(DisplayLines):
    """The undisplayLine command removes the line repesentation for the specified set ot atoms

    Synopsis:
        None <- undisplayL:ines(atoms)

    arguments:
        atoms  --- set of atoms

    Package : PmvApp
    Module  : displayCmds
    Class   : UndisplayLines
    Command : undisplayLines

    Example:
        mol = Pmv.Mols[0]
        ## undisplay lines entire molecule with vdw radii
        pmv.undisplayLines(mol)
    """

    def __init__(self):
        DisplayLines.__init__(self)
        self.negate = True

    ## def initializeMolForCmd(self, mol):
    ##     if not self.initializedFor[mol]:
    ##         return
            #raise RuntimeError("undisplayLines called before displayLinesLines for molecule %s"%mol.name)

    def onAddCmdToApp(self):
        DisplayLines.onAddCmdToApp(self)
        if not self.app().commands.has_key('displayLines'):
            self.app().lazyLoad('displayCmds', commands=['displayLines'], package='PmvApp')

class DisplayIndexedGeomCommand(DisplayCommand):

    _shadingModes = ['flat', 'smooth']

    _cullingModes = ['front', 'back', 'none']

    def updateModelGlobals(self, mol, propDict, frontRendering=None, backRendering=None,
                           linewidth=None, pointwidth=None, quality=None,
                           shading=None, culling=None, backfaceColor=None, sharpColorBoundaries=None, **kw):
        # update the CPK display model gloabl options
        # frontRendering and backRendering: string, default:'solid', or 'mesh' or 'points' global (i.e. affect all CPK spheres for this molecule)
        #           
        # linewidth: integer > 1 specifying the with of lines for drawing spheres as wires
        # quality: integer specifying the level of tesselation of the spheres. When set to 0
        #          the tesselation will adapt to molecule size (i.e. molecules with few atoms
        #          will use higher quality sphere than large molecules).

        self.initialize(mol)

        if frontRendering is not None:
            assert frontRendering in self._renderToDejaVu.keys()
            propDict['frontRendering'] = self._renderToDejaVu[frontRendering]
        if backRendering is not None:
            assert backRendering in self._renderToDejaVu.keys()+['sameAsFront']
            propDict['backRendering'] = self._renderToDejaVu[backRendering]
        if quality is not None:
            assert isinstance(quality, int) and quality >=0
            propDict['quality'] = quality
        if linewidth is not None:
            assert isinstance(linewidth, int) and linewidth >=1
            propDict['linewidth'] = linewidth
        if pointwidth is not None:
            assert isinstance(pointwidth, int) and pointwidth >=1
            propDict['pointwidth'] = pointwidth
        if shading is not None:
            assert shading in self._shadingModes
            propDict['shading'] = shading
        if culling is not None:
            assert culling in self._cullingModes
            propDict['culling'] = culling
        if  backfaceColor is not None:
            assert backfaceColor=='sameAsFront' or len(backfaceColor) in [3,4]
            propDict['backfaceColor'] = backfaceColor
        if sharpColorBoundaries is not None:
            assert sharpColorBoundaries in [1, 0, False, True]
            propDict['sharpColorBoundaries'] = sharpColorBoundaries

##############################################################################################
##
## DISPLAY CPK
##
##############################################################################################
class DisplayCPK(DisplayIndexedGeomCommand):
    """The displayCPK command displays a set of atoms as spheres.

    Synopsis:
        None <- displayCPK(atoms, radii=None,
                           rendering='solid', linewidth=2, pointwidth=2, quality=0)

    arguments:
        atoms  --- set of atoms

        radii  --- radii used to draw atomic spheres, Valid values are:
                       None         : use the current radius for each atom (initial default vdw radii)
                       float        : all atoms get the same radius for display
                       List(floats) : each atom gets its own display radius

        # global parameters: (i.e. affecting all cpk spheres of a molecule)

        frontRendering -- can be 'solid', 'mesh', or 'points'. Default:'solid'(rendering mode of the front and back polygons) 

        backRendering --- can be 'solid', 'mesh', or 'points' or 'sameAsFront'
        
        linewidth  --- integer > 1 specifying the with of lines for drawing spheres as a mesh

        pointwidth --- integer > 1 specifying the size of poitns for drawing spheres as points

        quality    --- integer specifying the level of tesselation of the spheres. When set to 0
                       the tesselation will adapt to molecule size (i.e. molecules with few atoms
                       will use higher quality sphere than large molecules).

        backfaceColor --- a color for back faces (numpy array [[r,g,b]]) or the string 'sameAsFront'

        culling  --- face culling; can be 'front', 'back', 'none'

        shading  --- shading technique; can be 'flat' or 'smooth'.

    Global model parameters can be set by calling:

             displayCPK.updateModelGlobals(mol, rendering=None, linewidth=None, pointwidth=None, quality=None)

         This call will not trigger a refresh of the display. To see the effect a call to
         pmv.displayCPK.refreshDisplay(mol) is needed
         
    Model parameters can be set by calling:

              displayCPK.updateModel(mol, sel, radii=None)

         This call will not trigger a refresh of the display. To see the effect a call to
         pmv.displayCPK.refreshDisplay(mol) is needed

    the list of provided arguments is checked when displayCPK is called. If an name is not recognized
    a RuntimeError exception is Thrown

    Package : PmvApp
    Module  : displayCmds
    Class   : DisplayCPK
    Command : displayCPK

    Examples:
        mol = Pmv.Mols[0]

        ## display CPK for entire molecule with vdw radii
        pmv.displayCPK(mol)

        ## hide all CPK spheres for mol
        pmv.undisplayCPK(mol)

        ## display CPK for atoms in 2 first amino acids
        atomsIn2FirstAA = mol.select('resindex < 3')
        pmv.displayCPK(atomsIn2FirstAA)

        ## display CPK for atoms in amino acids 4 and 5 with radius 0.5
        atomsIn2FirstAA = mol.select('3 < resindex < 6')
        pmv.displayCPK(atomsIn2FirstAA, radii=.5)

        ## display CPK for atoms in amino acids 6 and 7 with radius 0.5
        ## with a radius proportional to their atomic number 
        atoms = mol.select('5 < resindex < 8')
        # get atomic temperature factors for the atoms
        atomnums = mol._ag.getBetas()[atoms.getIndices()]
        # normalize the temperature factor adn multiply by 2.
        atomnums = 2*atomnums/max(atomnums)
        # display the spheres
        pmv.displayCPK(atoms, radii=atomnums)

        ## undisplay CPK for atoms in 2 first amino acids with radius 1
        atomsIn2FirstAA = mol.select('3 < resindex < 6')
        pmv.displayCPK(atomsIn2FirstAA, radii=.5)
        atSet = mol.geomContainer.atoms['cpk']

        ## make all CPK sphere wire with linewidth 1
        pmv.displayCPK.updateModelGlobals(mol, frontRendering=None, linewidth=None)
        pmv.displayCPK.refreshDisplay(mol)

        # check that a bad keyword argument will cause a Runtime Error
        try:
            pmv.displayCPK(mol, radius=2)
            raise ValueError("failed to raise RuntimeError for bad keyword argument")
        except RuntimeError, e:
            print e

        >>> displayCPK: unrecognized keyword argument 'radius', valid names are ['rendering', 'linewidth', 'pointwidth', 'quality', 'radii']
    """

    _argNames = ['frontRendering', 'backRendering', 'linewidth', 'pointwidth', 'quality', 'radii', 'culling', 'shading', 'backfaceColor']

    def __init__(self):
        DisplayIndexedGeomCommand.__init__(self)
        self._atomSetName = "cpk"

    ## these commands are now builtin
    ## def onAddCmdToApp(self):
    ##     DisplayCommand.onAddCmdToApp(self)

    ##     if not self.app().commands.has_key('undisplayCPK'):
    ##         self.app().lazyLoad('displayCmds', commands=['undisplayCPK'], package='PmvApp')
    ##     if self.app().undisplayCPK.isLoader():
    ##         self.app().undisplayCPK.loadCommand()
            
    ##     if not self.app().commands.has_key('assignAtomsRadii'):
    ##         self.app().lazyLoad('editCmds', commands=['assignAtomsRadii'], package='PmvApp')
 
    def initializeMolForCmd(self, mol):
        """Adds the cpk geometry and the cpk atom group to the object's
        geometry container

        adds per atom datafields:
            cpkRadScale
            cpkRadOff
            colorsIndices_cpk
        adds molecule._renderingProp['cpk'] = {'quality':0}
        """
        if self.initializedFor.get(mol, False):
            return
        self.initializedFor[mol] = True
        self.app().undisplayCPK.initializedFor[mol] = True

        if mol._ag._data.get('cpk_radius', None) is None:
            mol._ag.setData('cpk_radius', mol._ag.getRadii())

        if mol._renderingProp.get('cpk', None) is None:
            mol._renderingProp['cpk'] = {'quality':0, # sphere tesselation
                                         'linewidth':2, # line width
                                         'pointwidth':2, # line width
                                         'frontRendering':'fill', # frontPolyMode
                                         'backRendering': 'fill', #backPolyMode
                                         'culling':'back',# 'front', 'front_and_back', 'none'
                                         'shading': 'smooth', # 'flat', 'inherit'
                                         'backfaceColor': [ 1.,  1.,  1., 1.],
                                         }

        if mol._colors.get('cpk', None) is None:
            # per atom properties used for lines
            # initialize with line colors
            mol._colors['cpk'] = mol._colors['lines'].copy()

        if mol._ag.getData('colorsIndices_cpk') is None:
            mol._ag.setData('colorsIndices_cpk',
                            mol._ag.getData('colorsIndices_lines').copy())

        geomC = mol.geomContainer

        # CPK spheres
        if not geomC.geoms.has_key('cpk'):
            g = Spheres(self._atomSetName, quality=0, visible=0, protected=True,
                        inheritMaterial=False, inheritFrontPolyMode=False)
            geomC.addGeom(g, parent=geomC.masterGeom, redo=0,
                          displayCmd=self, undisplayCmd=self.app().undisplayCPK)
            g.vertexSet.vertices.array = geomC.allCoords
            g._hasAllVertices = True
            self.managedGeometries.append(g)
            geomC.geomPickToBonds['cpk'] = None
            self.managedGeometries.append(g)
        
    def doit(self, molSel, **kw):
        sel = molSel.select('not deleted')
        mol = sel.getAtomGroup().getMolecule()

        self.initialize(mol)

        indices = sel.getIndices()
        cpkRad = mol._ag.getData("cpk_radius")
        
        #setup undo
        self.app().pushUndoCmd(self.app().displayCPK if self.negate else self.app().undisplayCPK,
                               (sel,), {'frontRendering': self._renderToPmv[mol._renderingProp['cpk']['frontRendering']],
                                        'backRendering': self._renderToPmv[mol._renderingProp['cpk']['backRendering']],
                                        'linewidth': mol._renderingProp['cpk']['linewidth'],
                                        'pointwidth': mol._renderingProp['cpk']['pointwidth'],
                                        'quality': mol._renderingProp['cpk']['quality'],
                                        'culling': mol._renderingProp['cpk']['culling'],
                                        'shading': mol._renderingProp['cpk']['shading'],
                                        'backfaceColor':mol._renderingProp['cpk']['backfaceColor'],
                                        'radii': cpkRad[indices],
                                                                                })

        self.updateModelGlobals(mol, mol._renderingProp['cpk'], **kw)
        self.updateModel(sel, **kw)
        self.refreshDisplay(mol)
        self.app().labelAtoms.refreshDisplay(mol) # to make the labels be in front of CPK speheres
        self.app().labelResidues.refreshDisplay(mol) # to make the labels be in front of MSMS

    def updateModel(self, sel, radii=None, **kw):
        # update the CPK display model for the atoms in sel
        # radii can be :
        #      None   : use the vdw radii (i.e. del.getRadii())
        #      float  : all atoms get the same radius for display
        #      [float]: each atom gets its own display radius
        #
        mol = sel.getAtomGroup().getMolecule()
        self.initialize(mol)

        indices = sel.getIndices()
        cpkRad = mol._ag.getData("cpk_radius")

        # handle radii
        if radii is None:
            allRads = mol._ag.getRadii()
            cpkRad[indices] = allRads[indices]
        elif isinstance(radii, float):
            cpkRad[indices] = radii
        else:
            try:
                assert len(radii)==len(sel)
                cpkRad[indices] = radii
            except TypeError:
                raise ValueError("%s bad parameter radius, expected None, float or [float], got %s"%(
                    self.name, str(radii)))
        mol._ag.setData("cpk_radius", cpkRad)
        
        # update atoms set for which CPK spheres are shown
        gc = mol.geomContainer
        ggeoms = gc.geoms
        _set = gc.atoms['cpk']

        if mol._multi == 'molecules':
            if mol._renderingBits & mol._CPK: # if set
                mol._renderingBits -= mol._CPK # remove the bit

        if len(sel) == len(mol.select()):
            if self.negate:
                _set = mol.emptySelection()
            else:
                _set = sel
                if mol._multi == 'molecules':
                    mol._renderingBits += mol._CPK # set the bit
        else:
            ##if negate, remove current atms from displayed _set
            if self.negate:
                _set = _set - sel
            else:
                _set = _set | sel

        # at this point _set is what will still be displayed as CPK after this cmd
        gc.setAtomsForGeom('cpk', _set)
        
    def refreshDisplay(self, mol):
        if not self.isInitialized(mol):
            return
        
        atoms = mol.geomContainer.atoms.get('cpk', mol.emptySelection())
        geom = mol.geomContainer.geoms['cpk']

        if len(atoms)==0: # nothing is diplayed as CPK anymore for this molecule
            geom.Set(visible=0, tagModified=False) # hide the geom
        else:
            cpkRad = mol._ag._data['cpk_radius']
            setInds = atoms.getIndices()
            globProp = mol._renderingProp['cpk']
            colorInds = mol._ag.getData('colorsIndices_cpk')
            colors = mol._colors['cpk'][colorInds]
            opacity = mol._ag.getData('opacity_cpk')
            bfcolor = globProp['backfaceColor']
            if bfcolor=='sameAsFront':
                geom.frontAndBack=True
                polyFace = 'front'
                #polyFace = 'frontandback'
            else:
                geom.frontAndBack=False
                geom.Set(materials=[bfcolor], polyFace='back')
                polyFace = 'front'
                    
            selected = self.app().activeSelection & SelectionSet([atoms])
            if selected:
                highlight = numpy.zeros( (len(geom.vertexSet),), 'i')
                highlight[selected[0].getIndices()] = 1
            else:
                highlight=[]
            geom.Set(faces=[[x] for x in setInds], materials=colors,
                     radii=cpkRad, visible=1, tagModified=False,
                     quality=globProp['quality'],
                     inheritLineWidth=False, lineWidth=globProp['linewidth'],
                     inheritPointWidth=False, pointWidth=globProp['pointwidth'],
                     frontPolyMode=globProp['frontRendering'],
                     culling=globProp['culling'],
                     shading=globProp['shading'],
                     backPolyMode=globProp['backRendering'],
                     polyFace=polyFace, highlight=highlight,
                     opacity=opacity, transparent="implicit")

class UndisplayCPK(DisplayCPK):
    """The undisplayCPK command removes the sphere repesentation for the specified set ot atoms

    Synopsis:
        None <- undisplayCPK(atoms)

    arguments:
        atoms  --- set of atoms

    Package : PmvApp
    Module  : displayCmds
    Class   : UndisplayCPK
    Command : undisplayCPK

    Examples:
        mol = Pmv.Mols[0]

        ## display CPK for entire molecule with vdw radii
        pmv.displayCPK(mol)
        pmv.
    """
    def __init__(self):
        DisplayCPK.__init__(self)
        self.negate = True

    ## def initializeMolForCmd(self, mol):
    ##     if not self.initializedFor[mol]:
    ##         return
            #raise RuntimeError("undisplayCPK called before displayCPK for molecule %s"%mol.name)
        
    def onAddCmdToApp(self):
        DisplayCPK.onAddCmdToApp(self)
        if not self.app().commands.has_key('displayCPK'):
            self.app().lazyLoad('displayCmds', commands=['displayCPK'], package='PmvApp')

##############################################################################################
##
## DISPLAY S&B
##
##############################################################################################

class DisplaySticksAndBalls(DisplayIndexedGeomCommand):
    """
    The displaySticksAndBalls command displays small spheres for
    atom and cylinders for bonds between these atoms.

    Synopsis:
        None <- displaySticksAndBalls(atoms, ballsRadii=None, cylRadii=None, 
                                      displayBondOrder=True,
                                      defaultBallRad=0.3, defaultCylRad=0.2,
                                      aromaticSphRad=0.15, aromaticSphPerAtom=2,
                                      doubleBondSep='auto', tripleBondSep='auto',
                                      frontRendering='solid', backRendering='mesh',
                                      linewidth=2, pointwidth=2, quality=4,
                                      shading='smooth', culling='back', backfaceColor=(1,1,1,1)
                                      )
    arguments:
        atoms      --- set of atoms

        ballsRadii --- radii used to draw atomic spheres, Valid values are:
                           None         : use the vdw radii (i.e. del.getRadii())
                           float        : all atoms get the same radius for display
                           List(floats) : each atom gets its own display radius

        cylRadii    --- radii used to draw atomic spheres, Valid values are:
                           None         : use the vdw radii (i.e. del.getRadii())
                           float        : all atoms get the same radius for display
                           List(floats) : each bond gets its own display radius

        # global parameters: (i.e. affecting all spheres and cylinders in the S&B representation of a molecule)

        displayBondOrder --- boolean. When True, double triple and aromatic bonds are displayed

        doubleBondSep --- separation between the 2 cylinders drawn for a double bond, fefault value 'auto'
                          the double bond cylinders are 1/2 the value of the bond radius of the double bond.
                          the 'auto separation will be 1.0*bondRadius where bondRadius is the cylinder radius
                          assigned to the double bond.

        tripleBondSep --- separartion betweent the central bond and the 2 bonds flanking it. The radius
                          of the central cylinder is 0.5*bondRadius and the 2 flanking ones have radii
                          0.25*bondRadius where bondRadius is the cylinder radius assigned to the tripple bond.
                          'auto' will separate the flanking cylinders by 1.0 bondRadius from the central one. 

        aromaticSphRad --- size of sphere used to draw aromatic circle. (default 0.15)

        aromaticSphPerAtom --- number of spheres use to draw the arc corresponding to an atom in the aromatic circle
                               (default 2)

        frontRendering -- can be 'solid', 'mesh', or 'points'. Default:'solid'(rendering mode of the front and back polygons) 

        backRendering --- can be 'solid', 'mesh', or 'points' or 'sameAsFront'
        
        linewidth  --- integer > 1 specifying the with of lines for drawing spheres as a mesh

        pointwidth --- integer > 1 specifying the size of poitns for drawing spheres as points

        quality    --- integer specifying the level of tesselation of the spheres. When set to 0
                       the tesselation will adapt to molecule size (i.e. molecules with few atoms
                       will use higher quality sphere than large molecules).

        backfaceColor --- a color for back faces (numpy array [[r,g,b]]) or the string 'sameAsFront'

        culling  --- face culling; can be 'front', 'back', 'none'

        shading  --- shading technique; can be 'flat' or 'smooth'.

        defaultBallRad --- radius used for balls when ballsRadii is None (default 0.3)

        defaultCylRad --- radius used for balls when cylRadii is None (default 0.2)

    Global model parameters can be set by calling:

             displaySB.updateModelGlobals(mol, displayBondOrder=None,
                                          defaultBallRad=None, defaultCylRad=None,
                                          aromaticSphRad=None, aromaticSphPerAtom=None,
                                          doubleBondSep=None, tripleBondSep=None,
                                          frontRendering=None, backRendering=None,
                                          linewidth=None, pointwidth=None, quality=None,
                                          shading=None, culling=None, backfaceColor=None,
                                          )

         This call will not trigger a refresh of the display. To see the effect a call to
         pmv.displaySB.refreshDisplay(mol) is needed
         
    Model parameters can be set by calling:

              displaySB.updateModel(mol, sel, ballsRadii=None, cylRadii=None, stippled=False)

         This call will not trigger a refresh of the display. To see the effect a call to
         pmv.displaySB.refreshDisplay(mol) is needed

    the list of provided arguments is checked when displaySB is called. If an name is not recognized
    a RuntimeError exception is Thrown

    Package : PmvApp
    Module  : displayCmds
    Class   : DisplaySticksAndBalls
    Command : displaySB

    Examples:
        mol = Pmv.Mols[0]

        ## display SticksAndBalls for entire molecule with vdw radii
        pmv.displaySB(mol)

        ## hide all Sticks and Balls spheres
        pmv.undisplaySB(mol)

        # specify ball radii and cylinder radii for a subset of atoms
        pmv.displaySB(mol.select('resname PRO'), ballsRadii=0.5, cylRadii=0.1)

        # licorice for all proline residues
        pmv.displaySB(mol.select('resname PRO'), ballsRadii=0.3, cylRadii=0.3)
    """

    _argNames = ['ballsRadii', 'cylRadii', 'displayBondOrder',
                 'aromaticSphRad', 'aromaticSphPerAtom',
                 'doubleBondSep', 'tripleBondSep',
                 'frontRendering', 'backRendering',
                 'linewidth', 'pointwidth', 'quality',
                 'shading', 'culling', 'backfaceColor',
                 'defaultBallRad', 'defaultCylRad', 'sharpColorBoundaries']

    ## MS stippled cylinders are not fucntional
                 
    def __init__(self):
        DisplayIndexedGeomCommand.__init__(self)
        self._atomSetName = "sb"

    ## these commands are now builtin
    ## def onAddCmdToApp(self):
    ##     DisplayCommand.onAddCmdToApp(self)

    ##     if not self.app().commands.has_key('undisplaySB'):
    ##         self.app().lazyLoad('displayCmds', commands=['undisplaySticksAndBalls'], package='PmvApp')
    ##     if self.app().undisplaySB.isLoader():
    ##         self.app().undisplaySB.loadCommand()

    ##     if not self.app().commands.has_key('assignAtomsRadii'):
    ##         self.app().lazyLoad('editCmds', commands=['assignAtomsRadii'], package='PmvApp')
        
    def initializeMolForCmd(self, mol):
        """
        Creates the balls, the sticks, and the stippled balls
        and sticks geometries. Adds these new geometries to the
        geomContainer.geoms of the new molecule.

        adds per atom datafields:
             colorsIndices_sb_balls
             colorsIndices_sb_cyl
             stippled
        """
        if self.initializedFor.get(mol, False):
            return
        self.initializedFor[mol] = True
        self.app().undisplaySB.initializedFor[mol] = True

        patoms = mol._ag
        if mol._ag._data.get('sb_ballsRadius', None) is None:
            mol._ag.setData('sb_ballsRadius', mol._ag.getRadii()*0.3)

        if mol._ag._bondData.get('sb_cylRadius', None) is None:
            mol.addAGBondData('sb_cylRadius',
                              numpy.ones(len(mol._ag._bonds,), 'f')*0.2, 0.2)

        if mol._ag.getFlags('stippled_sb') is None:
            mol._ag.setFlags('stippled_sb', [False]*len(mol._ag))

        if mol._colors.get('sb_balls', None) is None:
            # per atom properties used for lines
            # initialize with line colors
            mol._colors['sb_balls'] = mol._colors['lines'].copy()
            mol._colors['sb_cyl'] = mol._colors['lines'].copy()

        if mol._ag.getData('colorsIndices_sb_balls') is None:
            mol._ag.setData('colorsIndices_sb_balls',
                            mol._ag.getData('colorsIndices_lines').copy())
            mol._ag.setData('colorsIndices_sb_cyl',
                            mol._ag.getData('colorsIndices_lines').copy())

        if mol._renderingProp.get('sb', None) is None:
            mol._renderingProp['sb'] = {
                'displayBondOrder':True,
                'defaultBallRad':0.3,
                'defaultCylRad':0.2,
                'doubleBondSep':'auto',
                'tripleBondSep':'auto',
                'aromaticSphRad': 0.15,
                'aromaticSphPerAtom': 2,
                'frontRendering':'fill', 'backRendering':'line',
                'linewidth':2, 'pointwidth':2, 'quality':0,
                'shading':'smooth', 'culling':'back', 'backfaceColor':(1,1,1,1),
                'sharpColorBoundaries':True,
                'stippleLength':0.2,
                'stippleSpace':0.2
                }
        if not hasattr(mol, '_aromaticSphereCenters'):
            mol._aromaticSphereCenters = {}
            # mol._aromaticSphereCenters hold a list of sphere centers
            # associated to atoms that are aromatic

        geomC = mol.geomContainer

        if not geomC.geoms.has_key("sb"):
            sb = Geom(self._atomSetName, shape=(0,0),
                      inheritLineWidth=False,
                      protected=True, pickable=0, inheritSharpColorBoundaries=False)
            geomC.addGeom(sb, parent=geomC.masterGeom, redo=0,
                          displayCmd=self,
                          undisplayCmd=self.app().undisplaySB) 
            self.managedGeometries.append(sb)

            # Cylinders (sticks)
            g = Cylinders("sticks", visible=0, inheritMaterial=False,
                          protected=True, pickable=0, inheritSharpColorBoundaries=False)
            geomC.addGeom(g, parent=sb, redo=0, makeAtomSet=False,
                          displayCmd=self, undisplayCmd=self.app().undisplaySB)
            g.vertexSet.vertices.array = geomC.allCoords
            g._hasAllVertices = True
            self.managedGeometries.append(g)

            g = Cylinders("sticksStippled", visible=0, inheritMaterial=False,
                          protected=True, pickable=0, inheritSharpColorBoundaries=False)
            geomC.addGeom(g, parent=sb, redo=0, makeAtomSet=False,
                          displayCmd=self,
                          undisplayCmd=self.app().undisplaySB)
            self.managedGeometries.append(g)

            g = Spheres("ballsStippled", visible=0, inheritMaterial=False,
                        protected=True, pickable=0)
            geomC.addGeom(g, parent=sb, redo=0, makeAtomSet=False,
                          displayCmd=self,
                          undisplayCmd=self.app().undisplaySB)
            self.managedGeometries.append(g)

            # Spheres at atomic positions (used for stick and balls)
            g = Spheres( "balls", radii=0.15, quality=4 ,visible=0,
                         inheritMaterial=False, protected=True)
            geomC.addGeom(g, parent=sb, redo=0, makeAtomSet=False,
                          displayCmd=self, undisplayCmd=self.app().undisplaySB)
            g.vertexSet.vertices.array = geomC.allCoords  
            g._hasAllVertices = True
            self.managedGeometries.append(g)
            #geomC.geomPickToBonds['balls'] = None

            # the bo cylinders mostly end inside ball spheres so we don;t draw extra spheres
            #g = Spheres("boBalls", visible=0, inheritMaterial=False, protected=True)
            #geomC.addGeom(g, parent=sb, redo=0, makeAtomSet=False)
            #self.managedGeometries.append(g)

            bonds = mol.select().getBonds()
            # DoubleBond Cylinders (sticks)
            g = Cylinders("doubleBondsSticks", visible=0, pickable=0,
                          inheritMaterial=False, protected=True, inheritSharpColorBoundaries=False)
            geomC.addGeom(g, parent=sb, redo=0, makeAtomSet=False,
                          displayCmd=self, undisplayCmd=self.app().undisplaySB)
            self.managedGeometries.append(g)
            self.updateDoubleBondsVertices(mol, bonds[2])

            g = Cylinders("doubleBondsSticksStippled", visible=0, pickable=0,
                          inheritMaterial=False, protected=True, inheritSharpColorBoundaries=False)
            geomC.addGeom(g, parent=sb, redo=0, makeAtomSet=False,
                          displayCmd=self, undisplayCmd=self.app().undisplaySB)
            self.managedGeometries.append(g)
                    
            # TripleBond Cylinders (sticks)
            g = Cylinders( "tripleBondsSticks", visible=0, pickable=0,
                           inheritMaterial=False, protected=True, inheritSharpColorBoundaries=False)
            geomC.addGeom(g, parent=sb, redo=0, makeAtomSet=False,
                          displayCmd=self, undisplayCmd=self.app().undisplaySB)
            self.managedGeometries.append(g)
            self.updateTripleBondsVertices(mol, bonds[3])

            g = Cylinders("tripleBondsSticksStippled", visible=0, pickable=0,
                          inheritMaterial=False, protected=True, inheritSharpColorBoundaries=False)
            geomC.addGeom(g, parent=sb, redo=0, makeAtomSet=False,
                          displayCmd=self, undisplayCmd=self.app().undisplaySB)
            self.managedGeometries.append(g)

            g = Spheres("aromaticSpheres", radii=0.15, quality=4, visible=0,
                        pickable=0, inheritMaterial=False, protected=True)
            geomC.addGeom(g, parent=sb, redo=0, makeAtomSet=False,
                          displayCmd=self, undisplayCmd=self.app().undisplaySB)
            self.managedGeometries.append(g)
            geomC.geomPickToBonds['aromaticSpheres'] = None
            if len(bonds[4]):
                aromaticArcs = mol._bondOrderData.get('aromaticArcs', None)
                if aromaticArcs is None or len(aromaticArcs)==0:
                    mol.computeAromaticRingsArcs()
                if aromaticArcs is not None:
                    self.updateAromaticArcs(mol, bonds[4])

    def doit(self, molSel, **kw):
        sel = molSel.select('not deleted')
        mol = sel.getAtomGroup().getMolecule()
        self.initialize(mol)

        indices = sel.getIndices()
        brad = mol._ag.getData("sb_ballsRadius")
        crad = mol._ag._bondData['sb_cylRadius']
        bonds = self._bonds = sel.getBonds()
        self._bnums = [mol._ag._bondIndex['%d %d'%(b[0],b[1])] for b in bonds[1]+bonds[2]+bonds[3]]
        #setup undo
        # FIXME to undo stipple we neet to set sipple array
        self.app().pushUndoCmd(self.app().displaySB if self.negate else self.app().undisplaySB,
                               (sel,), {'ballsRadii': brad[indices], 'cylRadii':crad[self._bnums],
                                        'frontRendering': self._renderToPmv[mol._renderingProp['sb']['frontRendering']],
                                        'backRendering': self._renderToPmv[mol._renderingProp['sb']['backRendering']],
                                        'linewidth': mol._renderingProp['sb']['linewidth'],
                                        'pointwidth': mol._renderingProp['sb']['pointwidth'],
                                        'quality': mol._renderingProp['sb']['quality'],
                                        'shading': mol._renderingProp['sb']['shading'],
                                        'culling': mol._renderingProp['sb']['culling'],
                                        'backfaceColor': mol._renderingProp['sb']['backfaceColor'],
                                        'defaultBallRad': mol._renderingProp['sb']['defaultBallRad'],
                                        'defaultCylRad': mol._renderingProp['sb']['defaultCylRad'],
                                        'displayBondOrder': mol._renderingProp['sb']['displayBondOrder'],
                                        'sharpColorBoundaries': mol._renderingProp['sb']['sharpColorBoundaries']
                                        })

        self.updateModelGlobals(mol, mol._renderingProp['sb'], **kw)
        self.updateModel(sel, **kw)
        self.refreshDisplay(mol)
        self._bonds = None
        self._bnums = None
        self.app().labelAtoms.refreshDisplay(mol) # to make the labels be in front of S&B
        self.app().labelResidues.refreshDisplay(mol) # to make the labels be in front of MSMS

    def updateModelGlobals(self, mol, propDict, displayBondOrder=None,
                           defaultBallRad=None, defaultCylRad=None,
                           aromaticSphRad=None, aromaticSphPerAtom=None,
                           doubleBondSep=None, tripleBondSep=None,
                           frontRendering=None, backRendering=None,
                           linewidth=None, pointwidth=None, quality=None,
                           shading=None, culling=None, backfaceColor=None,
                           sharpColorBoundaries=None, **kw):
        
        self.initialize(mol)
        
        DisplayIndexedGeomCommand.updateModelGlobals(self, mol, propDict,
                           frontRendering=frontRendering, backRendering=backRendering,
                           linewidth=linewidth, pointwidth=pointwidth, quality=quality,
                           shading=shading, culling=culling, backfaceColor=backfaceColor,
                           sharpColorBoundaries=sharpColorBoundaries)
                                                     
        if defaultBallRad is not None:
            assert isinstance(defaultBallRad, float)
            propDict['defaultBallRad'] = defaultBallRad

        if defaultCylRad is not None:
            assert isinstance(defaultCylRad, float)
            propDict['defaultCylRad'] = defaultCylRad

        if displayBondOrder is not None:
            assert isinstance(displayBondOrder, bool)
            propDict['displayBondOrder'] = displayBondOrder

        if doubleBondSep is not None:
            if isinstance(doubleBondSep, str):
                assert doubleBondSep=='auto'
            else:
                assert isinstance(doubleBondSep, float)
            
            sep = propDict['doubleBondSep']
            if sep != doubleBondSep:
                propDict['doubleBondSep'] = doubleBondSep
                self.updateDoubleBondsVertices(mol, mol.select().getBonds()[2])

        if tripleBondSep is not None:
            if isinstance(tripleBondSep, str):
                assert tripleBondSep=='auto'
            else:
                assert isinstance(tripleBondSep, float)
            sep = propDict['tripleBondSep']
            if sep != tripleBondSep:
                propDict['tripleBondSep'] = tripleBondSep
                self.updateTripleBondsVertices(mol, mol.select().getBonds()[3])

        if aromaticSphRad is not None:
            assert isinstance(aromaticSphRad, float)
            rad = propDict['aromaticSphRad']
            propDict['aromaticSphRad'] = aromaticSphRad

        if aromaticSphPerAtom is not None:
            assert isinstance(aromaticSphPerAtom, int)
            nb = propDict['aromaticSphPerAtom']
            if nb != aromaticSphPerAtom:
                propDict['aromaticSphPerAtom'] = aromaticSphPerAtom
                self.updateAromaticArcs(mol, mol.select().getBonds()[4])

    def updateDoubleBondsVertices(self, mol, bonds=None):
        if bonds is None:
            bonds = mol.select().getBonds()[2]
        sep = mol._renderingProp['sb']['doubleBondSep']
        cRad = mol._ag._bondData['sb_cylRadius']
        if len(bonds):
            verts = []
            bnum = 0
            for i1,i2 in bonds:
                if sep=='auto':
                    sep = cRad[bnum]
                verts.extend(mol.computeDoubleBond(i1, i2, sep))
                bnum+= 1
            mol.geomContainer.geoms['doubleBondsSticks'].vertexSet.vertices.array = numpy.array(verts, 'f')

    def updateTripleBondsVertices(self, mol, bonds=None):
        if bonds is None:
            bonds = mol.select().getBonds()[3]
        sep = mol._renderingProp['sb']['tripleBondSep']
        cRad = mol._ag._bondData['sb_cylRadius']
        if len(bonds):
            verts = []
            bnum = 0
            for  i1,i2 in bonds:
                if sep=='auto':
                    sep = 2*cRad[bnum]
                verts.extend(mol.computeTripleBond(i1, i2, sep))
                bnum+= 1
            mol.geomContainer.geoms['tripleBondsSticks'].vertexSet.vertices.array = numpy.array(verts, 'f')

    def updateModel(self, sel, ballsRadii=None, cylRadii=None, stippled=None, **kw):
        mol = sel.getAtomGroup().getMolecule()
        self.initialize(mol)

        indices = sel.getIndices()
        bRad = mol._ag.getData('sb_ballsRadius')
        cRad = mol._ag._bondData['sb_cylRadius']
        if self._bonds:
            bonds = self._bonds
        else:
            bonds = sel.getBonds()
         
        # set stippled flag for the selection
        if stippled is not None:
            mol._ag._flags['stippled_sb'][sel.getIndices()] = stippled

        # handle ball radii
        if ballsRadii is None:
            bRad[indices] = mol._renderingProp['sb']['defaultBallRad']
        elif isinstance(ballsRadii, float):
            bRad[indices] = ballsRadii
        else:
            try:
                assert len(ballsRadii)==len(sel)
                bRad[indices] = ballsRadii
            except TypeError:
                raise ValueError("%s bad parameter radius, expected None, float or [float], got %s"%(
                    self.name, str(ballsRadii)))
        mol._ag.setData("sb_ballsRadius", bRad)

        # handle cylinder radii
        if cylRadii is None:
            rad = mol._renderingProp['sb']['defaultCylRad']
            # set default cylinder radius for single bonds (also aromatic and triple)
            for b in bonds[1]:
                cRad[mol._ag._bondIndex['%d %d'%(b[0],b[1])]] = rad
            for b in bonds[2]:
                cRad[mol._ag._bondIndex['%d %d'%(b[0],b[1])]] = rad*.5
            # override cylinder radius for triple bonds
            for b in bonds[3]:
                cRad[mol._ag._bondIndex['%d %d'%(b[0],b[1])]] = rad*0.5
        elif isinstance(cylRadii, float):
            bnums = [mol._ag._bondIndex['%d %d'%(b[0],b[1])] for b in bonds[1]+bonds[2]+bonds[3]]
            cRad[bnums] = cylRadii
            #cRad[indices] = cylRadii
        else:
            try:
                #assert len(ballsRadii)==len(sel)
                #cRad[indices] = cylRadii
                bnums = [mol._ag._bondIndex['%d %d'%(b[0],b[1])] for b in bonds[1]+bonds[2]+bonds[3]]
                cRad[bnums] = cylRadii
            except TypeError:
                raise ValueError("%s bad parameter radius, expected None, float or [float], got %s"%(
                    self.name, str(cylRadii)))

        mol._ag._bondData['sb_cylRadius'] = cRad

        # recompute bo vertices in case cylinder radii changed
        if len(bonds[2]):
            self.updateDoubleBondsVertices(mol)
        if len(bonds[3]):
            self.updateTripleBondsVertices(mol)

        # update atom set
        gc = mol.geomContainer
        _set = gc.atoms['sb']

        if mol._multi == 'molecules':
            if mol._renderingBits & mol._SB:
                mol._renderingBits -= mol._SB
        if len(sel) == len(mol.select()):
            if self.negate:
                _set = mol.emptySelection()
            else: 
                _set = sel
                if mol._multi == 'molecules':
                    mol._renderingBits += mol._SB
        else:
            if self.negate:
                _set = _set - sel
            else:
                _set = _set | sel 

        gc.setAtomsForGeom('sb', _set)

    def refreshDisplay(self, mol):
        if not self.isInitialized(mol):
            return

        atomSets = mol.geomContainer.atoms
        gc = mol.geomContainer
        ggeoms = gc.geoms

        # get a list of all atom colors for sb
        ballColorInds = mol._ag.getData('colorsIndices_sb_balls')
        ballColors = mol._colors['sb_balls'][ballColorInds.tolist()]
        cylColorInds = mol._ag.getData('colorsIndices_sb_cyl')
        cylColors = mol._colors['sb_cyl'][cylColorInds.tolist()]
        opacity = mol._ag.getData('opacity_sb')
        if opacity is not None and len(numpy.unique(opacity))==1:
            opacity = [opacity[0]]
        # show balls
        atomsSB = atomSets.get('sb', mol.emptySelection())
        if len(atomsSB)==0:
            # hide sticks and balls
            ggeoms['balls'].Set(visible=False, tagModified=False)
            ggeoms['sticks'].Set(visible=False, tagModified=False)
            ggeoms['sticksStippled'].Set(visible=False, tagModified=False)
            ggeoms['ballsStippled'].Set(visible=False, tagModified=False)
            ggeoms['doubleBondsSticks'].Set(visible=False, tagModified=False)
            ggeoms['doubleBondsSticksStippled'].Set(visible=False, tagModified=False)
            ggeoms['tripleBondsSticks'].Set(visible=False, tagModified=False)
            ggeoms['tripleBondsSticksStippled'].Set(visible=False, tagModified=False)
            ggeoms['aromaticSpheres'].Set(visible=False, tagModified=False)
        else:
            globProp = mol._renderingProp['sb']
            selected = self.app().activeSelection & SelectionSet([atomsSB])
            if selected:
                highlight = numpy.zeros( (len(ggeoms['balls'].vertexSet),), 'i')
                highlight[selected[0].getIndices()] = 1
            else:
                highlight=[]
            
            bfcolor = globProp['backfaceColor']
            if bfcolor=='sameAsFront':
                ggeoms['balls'].frontAndBack=True
                polyFace = 'front'
            else:
                ggeoms['balls'].frontAndBack = False
                ggeoms['balls'].Set(materials=[bfcolor], polyFace='back')
                polyFace = 'front'

            # show balls
            indices = atomsSB.getIndices()
            ggeoms['balls'].Set(
                faces=[[x] for x in indices],
                radii = mol._ag.getData('sb_ballsRadius'),
                visible=1, tagModified=False, inheritMaterial=False, materials=ballColors, 
                quality=globProp['quality'],
                inheritLineWidth=False, lineWidth=globProp['linewidth'],
                inheritPointWidth=False, pointWidth=globProp['pointwidth'],
                frontPolyMode=globProp['frontRendering'],
                backPolyMode=globProp['backRendering'],
                culling=globProp['culling'],
                shading=globProp['shading'],
                polyFace=polyFace, highlight=highlight,
                opacity=opacity, transparent='implicit')

            # show sticks
            atomsBondsSolid = atomsSB.select('not stippled_sb')
            atomsBondsStippled = atomsSB.select('stippled_sb')
            allBonds = atomsSB.getBonds()

            for atoms, solid in zip([atomsBondsSolid, atomsBondsStippled], (True, False)):
                if atoms is None:
                    if solid:
                        ggeoms['sticks'].Set(visible=False, tagModified=False)
                        ggeoms['doubleBondsSticks'].Set(visible=False, tagModified=False)
                        ggeoms['tripleBondsSticks'].Set(visible=False, tagModified=False)
                    else:
                        ggeoms['sticksStippled'].Set(visible=False, tagModified=False)
                        ggeoms['doubleBondsSticksStippled'].Set(visible=False, tagModified=False)
                        ggeoms['tripleBondsSticksStippled'].Set(visible=False, tagModified=False)
                    continue

                # bonds[1] has single and aromatic, [2] has double and [3] triple bonds
                # bonds[4] has aromatic bonds only
                if not mol._renderingProp['sb']['displayBondOrder'] or mol._ag._bondOrder is None:
                    # draw all bonds as single bonds
                    singleBonds = allBonds[1]+allBonds[2]+allBonds[3]
                    allBonds[2] = []
                    allBonds[3] = []
                    allBonds[4] = []
                else:
                    singleBonds = allBonds[1]

                # segregate solid bonds and stippled bonds
                # bonds between atom that are solid and stippled will be solid
                sf = mol._ag._flags['stippled_sb']
                if solid:
                    lbonds = []
                    for b in singleBonds:
                        if not sf[b[0]] or not sf[b[1]]:
                            lbonds.append(b)
                else:
                    lbonds = []
                    for b in singleBonds:
                        if sf[b[0]] and sf[b[1]]:
                            lbonds.append(b)
                singleBonds = lbonds

                if len(singleBonds) == 0: # no single bonds
                    if solid:
                        ggeoms['sticks'].Set(visible=False, tagModified=False)
                    else:
                        ggeoms['sticksStippled'].Set(visible=False, tagModified=False)
                        ggeoms['ballsStippled'].Set(visible=False, tagModified=False)
                else:
                    # process single bonds
                    cylRad = mol._ag._bondData['sb_cylRadius']
                    if len(numpy.unique(cylColorInds))==1:
                        mat = [mol._colors['sb_cyl'][cylColorInds[0]]]
                    else:
                        mat = cylColors

                    bnums = [mol._ag._bondIndex['%d %d'%(b[0],b[1])] for b in singleBonds]
                    stippledBallsVerts = []
                    stippledBallsCol = []
                    stippledBallsRad = []
                    aromaticSpheresVerts = []
                    aromaticSpheresCols = []
                    op = None
                    addOpacity = False
                    if opacity is not None:
                        op = []
                        if len(opacity) == 1 : 
                            op = [opacity[0]]
                        else: addOpacity = True
                    selected = self.app().activeSelection & SelectionSet(
                        [Selection(mol._ag, numpy.unique(numpy.array(singleBonds).flatten()), "")])
                    if selected:
                        highlight = numpy.zeros( (len(ggeoms['sticks'].vertexSet),), 'i')
                        highlight[selected[0].getIndices()] = 1
                        
                    if solid:
                        
                        bfcolor = globProp['backfaceColor']
                        if bfcolor=='sameAsFront':
                            ggeoms['sticks'].frontAndBack=True
                            polyFace = 'front'
                        else:
                            ggeoms['sticks'].frontAndBack=False
                            ggeoms['sticks'].Set(materials=[bfcolor], polyFace='back')
                            polyFace = 'front'
                    
                        ggeoms['sticks'].Set( faces=numpy.array(singleBonds, 'int'), radii=cylRad[bnums],
                                              materials=cylColors, 
                                              visible=1, tagModified=False,
                                              quality=globProp['quality'],
                                              inheritLineWidth=False, lineWidth=globProp['linewidth'],
                                              inheritPointWidth=False, pointWidth=globProp['pointwidth'],
                                              frontPolyMode=globProp['frontRendering'],
                                              culling=globProp['culling'],
                                              shading=globProp['shading'],
                                              backPolyMode=globProp['backRendering'],
                                              polyFace=polyFace, highlight=highlight,
                                              opacity=opacity, transparent='implicit',
                                              sharpColorBoundaries=globProp['sharpColorBoundaries'])
                    else:
                        vert = mol._ag.getCoords()
                        v, f, c, r = computeHalfBonds(vert, singleBonds, cylColors, cylRad[bnums])
                        v1, f1, cf1, cv1, r1 = stippleLines(v, f, c, r,
                                                            segLen=mol._renderingProp['sb']['stippleLength'],
                                                            spaceLen=mol._renderingProp['sb']['stippleSpace'])
                        ggeoms['sticksStippled'].Set(tagModified=False,
                            visible=True, vertices=v1, faces=f1, materials=cf1, radii=r1,
                            lineWidth=mol._renderingProp['sb']['linewidth'], inheritLineWidth=False,
                            inheritPointWidth=False, pointWidth=globProp['pointwidth'],
                            quality=globProp['quality'], frontPolyMode=globProp['frontRendering'],
                            sharpColorBoundaries=globProp['sharpColorBoundaries'])
                        stippledBallsVerts.extend(v1)
                        stippledBallsCol.extend(cv1)
                        for v in r1:
                            stippledBallsRad.append(v)
                            stippledBallsRad.append(v)
                    ##
                    ## double bonds
                    ##
                    bondList = allBonds[2] # double bonds
                    # if we do not display bond order this list has been set to emty above
                    doubleBonds = mol._bondOrderData['doubleBonds']
                    if len(bondList)==0 or doubleBonds is None or len(doubleBonds)==0:
                        if solid:
                            ggeoms['doubleBondsSticks'].Set(visible=False, tagModified=False)
                        else:
                            ggeoms['doubleBondsSticksStippled'].Set(visible=False, tagModified=False)
                    else:
                        faces = []
                        #col = numpy.array([[1.,1.,1.,1.]]*len(ggeoms['doubleBondsSticks'].vertexSet.vertices))
                        col = []
                        rad = []
                        # selected atoms in double bonds
                        activeSel = self.app().activeSelection & SelectionSet(
                            [Selection(mol._ag, numpy.unique(numpy.array(bondList).flatten()), "")])
                        highlight = numpy.zeros( (len(ggeoms['doubleBondsSticks'].vertexSet),), 'i')
                        # loop over bonds
                        if opacity is not None:
                            if len(opacity) == 1 : 
                                op = [opacity[0]]
                            else:
                                op = numpy.ones(len(ggeoms['doubleBondsSticks'].vertexSet.vertices)).astype('f')
                                addOpacity = True
                        for i1,i2 in bondList:
                            k = "%d %d" %(i1, i2)
                            if doubleBonds.has_key(k):
                                n = doubleBonds[k] * 4
                                faces.append( (n, n+1) )
                                faces.append( (n+2, n+3) )
                                #col[[n,n+1,n+2,n+3]] = [cylColors[i1], cylColors[i2],
                                #                       cylColors[i1], cylColors[i2]]
                                col.extend([cylColors[i1], cylColors[i2],
                                                       cylColors[i1], cylColors[i2]])
                                if addOpacity:
                                    op[[n,n+1,n+2,n+3]] = [opacity[i1], opacity[i2], 
                                                           opacity[i1], opacity[i2]]
                                rad.append( cylRad[mol._ag._bondIndex[k]]*.5 )
                                rad.append( cylRad[mol._ag._bondIndex[k]]*.5 )
                                if activeSel and i1 in activeSel[0]._indices:
                                    highlight[n]= 1
                                    highlight[n+2]= 1
                                if activeSel and i2 in activeSel[0]._indices:
                                    highlight[n+1]= 1
                                    highlight[n+3]= 1

                        if solid:
                            bfcolor = globProp['backfaceColor']
                            if bfcolor=='sameAsFront':
                                ggeoms['doubleBondsSticks'].frontAndBack=True
                                polyFace = 'front'
                            else:
                                ggeoms['doubleBondsSticks'].frontAndBack=False
                                ggeoms['doubleBondsSticks'].Set(materials=[bfcolor], polyFace='back')
                                polyFace = 'front'
                            ggeoms['doubleBondsSticks'].Set(faces=faces, materials=col, radii=rad,
                                              visible=1, tagModified=False,
                                              quality=globProp['quality'],
                                              inheritLineWidth=False, lineWidth=globProp['linewidth'],
                                              inheritPointWidth=False, pointWidth=globProp['pointwidth'],
                                              frontPolyMode=globProp['frontRendering'],
                                              culling=globProp['culling'],
                                              shading=globProp['shading'],
                                              backPolyMode=globProp['backRendering'],
                                              polyFace=polyFace, highlight=highlight,
                                              opacity=opacity, transparent='implicit',
                                              sharpColorBoundaries=globProp['sharpColorBoundaries'])
                        else:
                            # double cylRadius
                            _r = []
                            for r in rad:
                                _r.append(r)
                                _r.append(r)
                            v, f, c, r = computeHalfBonds(ggeoms['doubleBondsSticks'].getVertices(), faces, col, _r)
                            v1, f1, cf1, cv1, r1 = stippleLines(
                                v, f, c, r,
                                segLen=mol._renderingProp['sb']['stippleLength'],
                                spaceLen=mol._renderingProp['sb']['stippleSpace'])
                            ggeoms['doubleBondsSticksStippled'].Set(tagModified=False,
                                visible=True, vertices=v1, faces=f1, materials=cf1, radii=r1,
                                lineWidth=mol._renderingProp['sb']['linewidth'], inheritLineWidth=False,
                                inheritPointWidth=False, pointWidth=globProp['pointwidth'],
                                quality=globProp['quality'], frontPolyMode=globProp['frontRendering'],
                                sharpColorBoundaries=globProp['sharpColorBoundaries'])
                            stippledBallsVerts.extend(v1)
                            stippledBallsCol.extend(cv1)
                            for v in r1:
                                stippledBallsRad.append(v)
                                stippledBallsRad.append(v)
                    ##
                    ## triple bonds
                    ##
                    bondList = allBonds[3] # triple bonds
                    tripleBonds = mol._bondOrderData['tripleBonds']
                    if len(bondList)==0 or tripleBonds is None or len(tripleBonds)==0:
                        if solid:
                            ggeoms['tripleBondsSticks'].Set(visible=False, tagModified=False)
                        else:
                            ggeoms['tripleBondsSticksStippled'].Set(visible=False, tagModified=False)
                    else:
                        faces = []
                        col = numpy.array([[1.,1.,1.,1,]]*len(ggeoms['tripleBondsSticks'].vertexSet.vertices))
                        op = None
                        addOpacity = False
                        if opacity is not None:
                            if len(opacity) == 1 : 
                                op = [opacity[0]]
                            else:
                                op = numpy.ones(len(ggeoms['tripleBondsSticks'].vertexSet.vertices)).astype('f')
                                addOpacity = True
                        rad = []
                        # selected atoms in triple bonds
                        activeSel = self.app().activeSelection & SelectionSet(
                            [Selection(mol._ag, numpy.unique(numpy.array(bondList).flatten()), "")])
                        highlight = numpy.zeros( (len(ggeoms['tripleBondsSticks'].vertexSet),), 'i')
                        for i1,i2 in bondList:
                            k = "%d %d" %(i1, i2)
                            if tripleBonds.has_key(k):
                                n = tripleBonds[k] * 6
                                faces.append( (n, n+1) )
                                faces.append( (n+2, n+3) )
                                col[[n,n+1,n+2,n+3]] = [cylColors[i1], cylColors[i2],
                                                        cylColors[i1], cylColors[i2]]
                                if addOpacity:
                                    op[[n,n+1,n+2,n+3]] = [opacity[i1], opacity[i2], 
                                                           opacity[i1], opacity[i2]]
                                rad.append( cylRad[mol._ag._bondIndex[k]]*0.25)
                                rad.append( cylRad[mol._ag._bondIndex[k]]*0.25)
                                if activeSel and i1 in activeSel[0]._indices:
                                    highlight[n]= 1
                                    highlight[n+2]= 1
                                if activeSel and i2 in activeSel[0]._indices:
                                    highlight[n+1]= 1
                                    highlight[n+3]= 1
                        if solid:
                            bfcolor = globProp['backfaceColor']
                            if bfcolor=='sameAsFront':
                                ggeoms['tripleBondsSticks'].frontAndBack=True
                                polyFace = 'front'
                            else:
                                ggeoms['tripleBondsSticks'].frontAndBack=False
                                ggeoms['tripleBondsSticks'].Set(materials=[bfcolor], polyFace='back')
                                polyFace = 'front'
                            ggeoms['tripleBondsSticks'].Set(faces=faces, materials=col, radii=rad,
                                              visible=1, tagModified=False,
                                              quality=globProp['quality'],
                                              inheritLineWidth=False, lineWidth=globProp['linewidth'],
                                              inheritPointWidth=False, pointWidth=globProp['pointwidth'],
                                              frontPolyMode=globProp['frontRendering'],
                                              culling=globProp['culling'],
                                              shading=globProp['shading'],
                                              backPolyMode=globProp['backRendering'],
                                              polyFace=polyFace, highlight=highlight,
                                              opacity=opacity, transparent='implicit',
                                              sharpColorBoundaries=globProp['sharpColorBoundaries'])
                        else:
                            # double cylRadius
                            _r = []
                            for r in rad:
                                _r.append(r)
                                _r.append(r)
                            v, f, c, r = computeHalfBonds(ggeoms['tripleBonds'].getVertices(), faces, col, _r)
                            v1, f1, cf1, cv1, r1 = stippleLines(
                                v, f, c, r,
                                segLen=mol._renderingProp['lines']['stippleLength'],
                                spaceLen=mol._renderingProp['lines']['stippleSpace'])
                            ggeoms['tripleBondsSticksStippled'].Set(
                                visible=True, vertices=v1, faces=f1, materials=cf1, radii=r1,
                                lineWidth=globProp['linewidth'], inheritLineWidth=False,
                                inheritPointWidth=False, pointWidth=globProp['pointwidth'],
                                quality=globProp['quality'],
                                frontPolyMode=globProp['frontRendering'],
                                culling=globProp['culling'],
                                shading=globProp['shading'],
                                backPolyMode=globProp['backRendering'],
                                polyFace=polyFace,
                                sharpColorBoundaries=globProp['sharpColorBoundaries'])
                            
                            stippledBallsVerts.extend(v1)
                            stippledBallsCol.extend(cv1)
                            for v in r1:
                                stippledBallsRad.append(v)
                                stippledBallsRad.append(v)
                    ##
                    ## aromatic bonds
                    ##
                    bondList = allBonds[4] # aromatic bonds
                    aromaticArcs = mol._bondOrderData['aromaticArcs']
                    if len(bondList)==0 or aromaticArcs is None or len(aromaticArcs)==0:
                        if solid:
                            ggeoms['aromaticSpheres'].Set(visible=False, tagModified=False)
                    else:
                        # loop over atoms displayed as S&B
                        highlight = numpy.zeros( (len(ggeoms['aromaticSpheres'].vertexSet),), 'i')
                        for ind in atoms.getIndices():
                            # mol._aromaticSphereCenters hold a list of sphere centers
                            # associated to atoms that are aromatic
                            v = mol._aromaticSphereCenters.get(ind, None)
                            if v is not None:
                                aromaticSpheresVerts.extend(v)
                                aromaticSpheresCols.extend( [ballColors[ind].tolist()]*len(v))
                                if addOpacity:
                                    op.extend([opacity[ind]]*len(v))

                    ##
                    ## configure the spheres at the end of stippled cylinders
                    ##
                    # done globaly because aromatic spheres are same with stipple or solid
                    if len(aromaticSpheresVerts):
                        bfcolor = globProp['backfaceColor']
                        if bfcolor=='sameAsFront':
                            ggeoms['sticks'].frontAndBack=True
                            polyFace = 'front'
                        else:
                            ggeoms['sticks'].frontAndBack=False
                            ggeoms['sticks'].Set(materials=[bfcolor], polyFace='back')
                            polyFace = 'front'

                        ggeoms['aromaticSpheres'].Set(vertices=aromaticSpheresVerts,
                                 materials=aromaticSpheresCols,
                                 visible=True, tagModified=False,
                                 quality=globProp['quality'],
                                 radii=(globProp['aromaticSphRad'],),                    
                                 inheritLineWidth=False, lineWidth=globProp['linewidth'],
                                 inheritPointWidth=False, pointWidth=globProp['pointwidth'],
                                 frontPolyMode=globProp['frontRendering'],
                                 culling=globProp['culling'],
                                 shading=globProp['shading'],
                                 backPolyMode=globProp['backRendering'],
                                 polyFace=polyFace, highlight=highlight,
                                 opacity=op, transparent='implicit')
                    else:
                        ggeoms['aromaticSpheres'].Set(visible=False, tagModified=False)
                    ##
                    ## configure the spheres at the end of stippled cylinders
                    if not solid:
                        if len(stippledBallsVerts):
                            ggeoms['ballsStippled'].Set(
                                visible=True, vertices=stippledBallsVerts,
                                tagModified=False, materials=stippledBallsCol,
                                radii=stippledBallsRad, sharpColorBoundaries=globProp['sharpColorBoundaries'])
                        else:
                            ggeoms['ballsStippled'].Set(visible=False, tagModified=False)

    def updateAromaticArcs(self, mol, bonds):
        if len(bonds)==0:
            return
        nb = mol._renderingProp['sb']['aromaticSphPerAtom']+1
        rad = mol._renderingProp['sb']['aromaticSphRad']
        aromaticArcs = mol._bondOrderData.get('aromaticArcs', None)
        if aromaticArcs is not None and len(aromaticArcs):
            verts = []
            vals = []
            col = []
            mol._aromaticSphereCenters = {}
            for atomIndex, arcParams in aromaticArcs.items():
                for aromAtomInd, center, normal, startVect, angle, radius in arcParams:
                    verts = self.arcVertices(center, numpy.array(startVect, 'f'),
                                             normal, radius, angle, nb)[:-1]
                    if mol._aromaticSphereCenters.get(atomIndex, None):
                        mol._aromaticSphereCenters[atomIndex].extend(verts)
                    else:
                        mol._aromaticSphereCenters[atomIndex] = verts

    def arcVertices(self, center, vec, normal, radius, angle, nvertices):
        # nvertices is the number of spheres per arc i.e. atom
        from math import pi, cos, sin
        angRad = angle*pi*0.00555555555556
        #nvertices = int(ang/self.degreesPerSegment) + 1
        d = angRad / (nvertices-1) # increment
        #a = 0
        a=d*0.5
        # starting angle
        vx, vy, vz = normal
        vec2 = numpy.zeros(3, 'f')
        vec2[0] = vec[1]*vz - vec[2]*vy
        vec2[1] = vec[2]*vx - vec[0]*vz
        vec2[2] = vec[0]*vy - vec[1]*vx
        points = []
        for j in range(nvertices):
            points.append( center + cos(a)*vec*radius + sin(a)*vec2*radius )
            a = a+d
        return points

class UndisplaySticksAndBalls(DisplaySticksAndBalls):
    """The undisplaySB command removes the sticks and balls repesentation for the specified set ot atoms

    Synopsis:
        None <- undisplaySB(atoms)

    arguments:
        atoms  --- set of atoms

    Package : PmvApp
    Module  : displayCmds
    Class   : UndisplaySB
    Command : undisplaySB

    Example:
        mol = Pmv.Mols[0]
        ## undisplay Sticks and Balls for entire molecule with vdw radii
        pmv.undisplaySB(mol)
    """

    def __init__(self):
        DisplaySticksAndBalls.__init__(self)
        self.negate = True

    ## def initializeMolForCmd(self, mol):
    ##     if not self.initializedFor[mol]:
    ##         return
            #raise RuntimeError("undisplaySB called before displaySB for molecule %s"%mol.name)

    def onAddCmdToApp(self):
        DisplaySticksAndBalls.onAddCmdToApp(self)
        if not self.app().commands.has_key('displaySB'):
            self.app().lazyLoad('displayCmds', commands=['DisplaySticksAndBalls'], package='PmvApp')


class DisplayBackboneTrace(DisplaySticksAndBalls):
    """The displayBackboneTrace command allows the user to display/undisplay the given nodes using the sticks and balls representation, where the bonds are represented by cylinders and the atoms by balls.The radii of the cylinders and the balls, and the quality of the spheres are user defined.The user can chose to display 'Licorice', 'Sticks only' or 'Sticks and Balls'.
    \nPackage : PmvApp
    \nModule  : displayCommands
    \nClass   : DisplayBackboneTrace
    \nCommand : displayBackboneTrace
    \nSynopsis:\n
        None <- displayBackboneTrace(nodes,  only=0, negate=0, noballs=0,
                                      bradii=0.4, bquality=4, cradius=0.2,
                                      absolute=1, **kw)\n
        nodes   : any set of MolKit2.Selection describing molecular components\n
        only    : Boolean flag specifying whether or not to only display
                  the current selection\n
        negate  : Boolean flag specifying whether or not to undisplay
                  the current selection\n
        cradius : specifies the cylinder radius\n
        bradii  : specifies the radius of the balls if displayed.\n
        bquality: specifies the quality of the balls if displayed.\n
        setScale --- when True atm.caballRad, atm.caballScale and atm.cacRad are set;\n
                     if set to False, atm.caballRad, atm.caballScale and atm.cacRad \n
                     are used  instead of bRad, bScale, and cradius. 
        sticksBallsLicorice --- string specifying the type of rendering,
                     can be 'Licorice', 'Sticks only' or 'Sticks and Balls'
        keywords: display BackboneTrace representation\n
    """

    def onAddObjectToViewer(self, obj):
        self.objectState[obj] = {'onAddObjectCalled':True}
        self.objectState[obj] = {'onAddObjectCalled':True}
        defaultValues = self.lastUsedValues['default']
        obj.allAtoms.caballRad = defaultValues['bRad']
        obj.allAtoms.caballScale = defaultValues['bScale']
        obj.allAtoms.cacRad = defaultValues['cradius']

        geomC = obj.geomContainer

        # Cylinders (sticks)
        g = Cylinders( "CAsticks", visible=0, vertices=geomC.allCoords,
                       protected=True)
        g._hasAllVertices = True
        geomC.addGeom(g)
        self.managedGeometries.append(g)

        # Spheres at atomic positions (used for stick and balls)
        g = Spheres( "CAballs", vertices = geomC.allCoords,
                     radii = 0.4, quality = 4 ,visible=0, protected=True)
        g._hasAllVertices = True
        geomC.addGeom(g)
        self.managedGeometries.append(g)
        geomC.geomPickToBonds['CAballs'] = None

        for atm in obj.allAtoms:
            #atm.colors['CAsticks'] = (1.0, 1.0, 1.0)
            atm.opacities['CAsticks'] = 1.0
            #atm.colors['CAballs'] = (1.0, 1.0, 1.0)
            atm.opacities['CAballs'] = 1.0


    def undoCmdBefore(self, nodes, only=False, negate=False,
                        bRad=0.3, bScale=0.0, bquality=0,
                        cradius=0.2, cquality=0, setScale=True, 
                        sticksBallsLicorice='Sticks and Balls',
                        redraw=True):
        
        kw={}
        defaultValues = self.lastUsedValues['default']

        kw['bRad'] = defaultValues['bRad']
        kw['bScale'] = defaultValues['bScale']
        kw['bquality'] = defaultValues['bquality']
        kw['cradius'] = defaultValues['cradius']
        kw['cquality'] = defaultValues['cquality']
        kw['sticksBallsLicorice'] = defaultValues['sticksBallsLicorice']
        
        caballset = AtomSet()
        castickset = AtomSet()
        for mol in self.app().Mols:
            if not self.objectState.has_key(mol):
                self.onAddObjectToViewer(mol)
            caballset =  caballset+mol.geomContainer.atoms['CAballs']
            castickset = castickset+mol.geomContainer.atoms['CAsticks']
           
        if len(caballset)==0:
            # no balls displayed
            if len(castickset) == 0: # negate
                # no stick displayed 
                kw['negate'] = True
                kw['redraw'] = True
                return ( [(self, (nodes,), kw)], self.name)
            else:
                # noballs was on
                kw['negate'] = False
                kw['redraw'] = True
                kw['only'] = True
                return ([(self, (castickset,), kw)], self.name)
        else:

            kw['redraw'] = True
            kw['only'] = True
            return ([(self, (castickset,), kw)], self.name)



    def doit(self, nodes, only=False, negate=False, bRad=0.3,
             bScale=0.0, bquality=0, cradius=0.2, cquality=0, 
             sticksBallsLicorice='Sticks and Balls', setScale=True, redraw=True):

        ########################################################
            
        def drawCABalls(mol, atm, only, noBalls, negate, bRad, bScale,
                      bquality, setScale):
            if setScale:
                atm.caballRad = bRad
                atm.caballScale = bScale
            
            _set = mol.geomContainer.atoms['CAballs']
            ## case noballs:
            if noBalls:
                if only:
                    if len(atm) == len(mol.allAtoms):
                        _set = atm
                    else:
                        _set = atm.union(_set)
                    mol.geomContainer.geoms['CAballs'].Set(visible=0,
                                                         tagModified=False)
                    mol.geomContainer.setAtomsForGeom('CAballs', _set)
                    return
                else:
                    negate = True
            
            ## Then handle only and negate    
            ##if negate, remove current atms from displayed _set
            if len(atm) == len(mol.allAtoms):
                if negate: _set = AtomSet([])
                else: _set = atm
            else:
                if negate:
                    _set = _set - atm
                else:
                    ##if only, replace displayed _set with current atms
                    if only:
                        _set = atm
                    else: 
                        _set = atm.union(_set)
            if len(_set) == 0:
                mol.geomContainer.geoms['CAballs'].Set(visible=0,
                                                     tagModified=False)
                mol.geomContainer.setAtomsForGeom('CAballs', _set)
                return
            
            mol.geomContainer.setAtomsForGeom('CAballs', _set)
            #_set.sort()
            bScale = numpy.array(_set.caballScale)
            bRad = numpy.array(_set.caballRad)
            aRad = numpy.array(_set.radius)
            ballRadii = bRad + bScale * aRad
            ballRadii = ballRadii.tolist()

            b = mol.geomContainer.geoms['CAballs']
            # this assumes that the lines geometry has been added
            bcolors = [x.colors.get('CAballs', x.colors['lines']) for x in _set]
            b.Set(vertices=_set.coords, radii=ballRadii,
                  materials=bcolors, inheritMaterial=False, 
                  visible=1, quality=bquality, tagModified=False)

            # highlight selection
            selMols, selAtms = self.app().getNodesByMolecule(self.app().activeSelection.get())
            lMolSelectedAtmsDict = dict( zip( selMols, selAtms) )
            lSelectedAtoms = lMolSelectedAtmsDict.get(mol, None)
            if lSelectedAtoms is not None:
                    lVertexClosestAtomSet = mol.geomContainer.atoms['CAballs']
                    if len(lVertexClosestAtomSet) > 0:
                        lVertexClosestAtomSetDict = dict(zip(lVertexClosestAtomSet,
                                                             range(len(lVertexClosestAtomSet))))
                        highlight = [0] * len(lVertexClosestAtomSet)
                        for i in range(len(lSelectedAtoms)):
                            lIndex = lVertexClosestAtomSetDict.get(lSelectedAtoms[i], None)
                            if lIndex is not None:
                                highlight[lIndex] = 1
                        b.Set(highlight=highlight)


        def drawCASticks(mol, atm, only, negate, cradius, cquality, setScale):

                atm.sort()
                if setScale:
                    atm.cacRad = cradius
                _set = mol.geomContainer.atoms['CAsticks']
                 
                if negate:
                    _set = mol.geomContainer.atoms['CAsticks']
                    _set = _set - atm
                         
                else:                       
                     if only: 
                            _set = atm                            
                     else: 
                            _set = atm.union(_set)
                            
                if len(_set) == 0:
                    mol.geomContainer.geoms['CAsticks'].Set(visible=0,
                                                      tagModified=False)
                    mol.geomContainer.setAtomsForGeom('CAsticks', _set)
                    return

                mol.geomContainer.setAtomsForGeom('CAsticks', _set)
             
                indices =[]
                                   
                for i in range(0,len(_set)):
                        if i+1 <=len(_set)-1:
                            indices.append((i,i+1))
                                     
                for ch in range(0,len(mol.chains)-1):
                      list =[]
                    
                      if len(mol.chains)>1: 
                        a=mol.chains[ch].getAtoms().get(lambda x: x.name=='CA')
                        a.sort()
                        for l in _set:
                            for k in a:
                                if l==k:            
                                    list.append(l)
                        i =len(list)
                        if _set !=list:
                            if (i-1,i) in indices:
                                ind=indices.index((i-1,i))
                                del indices[ind]

                for ch in range(0,len(mol.chains)):    
                    #finding tranformed coordsd to pass in to FindGap
                    chatoms=_set
                    mats=[]
                    for ca in chatoms:
                        c = self.app().transformedCoordinatesWithInstances(AtomSet([ca]))
                        mat = numpy.array(c[0], 'f')
                        mats.append(mat)
                        
                    from MolKit2.GapFinder import FindGap
                    #calling Find Gap such that chains with residues not connected will 
                    #have an attribute 'hasGap' and CA atoms have "gap" attribute
                    
                    if len(_set)!=0:
                         _set.sort()
                         mygap = FindGap(mats,mol.chains[ch],_set)
                         if hasattr(mol.chains[ch],'hasGap'):
                             for i in range(0,len(chatoms)):
                                 if hasattr(chatoms[i],"gap"):
                                     if i+1<=len(chatoms)-1:
                                         if chatoms[i].gap=='start':
                                             #indi=chatoms.index(chatoms[i])
                                             for at in _set:
                                                 if chatoms[i]==at:
                                                     indi =_set.index(at)
                                             if (indi,indi+1) in indices:
                                                 ind=indices.index((indi,indi+1))    
                                                 del indices[ind] 
                                            
                mol.geomContainer.setAtomsForGeom('CAsticks', _set)
                g = mol.geomContainer.geoms['CAsticks'] 
                if len(indices) == 0:
                    g.Set(visible=0, tagModified=False)
                else:
                    cRad = _set.cacRad
                    scolors = [x.colors.get('CAsticks', x.colors['lines']) for x in _set]
                    g.Set( vertices=_set.coords, faces=indices, radii=cRad,
                           materials=scolors, visible=1, quality=cquality,
                           tagModified=False,  inheritMaterial=False)

                # highlight selection
                selMols, selAtms = self.app().getNodesByMolecule(self.app().activeSelection.get())
                lMolSelectedAtmsDict = dict( zip( selMols, selAtms) )
                lSelectedAtoms = lMolSelectedAtmsDict.get(mol, None)
                if lSelectedAtoms is not None:
                    lVertexClosestAtomSet = mol.geomContainer.atoms['CAsticks']
                    if len(lVertexClosestAtomSet) > 0:
                        lVertexClosestAtomSetDict = dict(zip(lVertexClosestAtomSet,
                                                             range(len(lVertexClosestAtomSet))))
                        highlight = [0] * len(lVertexClosestAtomSet)
                        for i in range(len(lSelectedAtoms)):
                            lIndex = lVertexClosestAtomSetDict.get(lSelectedAtoms[i], None)
                            if lIndex is not None:
                                highlight[lIndex] = 1
                        g.Set(highlight=highlight)


        ########################################################
        
        molecules, atmSets = self.app().getNodesByMolecule(nodes, Atom)
        try:
            radii = molecules.allAtoms.radius
        except:
            self.app().assignAtomsRadii(molecules, 
                                     overwrite=False,
                                     united=False)

        for mol, atom in map(None, molecules, atmSets):
            if not self.objectState.has_key(mol):
                self.onAddObjectToViewer(mol)
            #for chain in mol.chains:
            try:
                if self._failForTesting:
                    a = 1/0.
                atm = atom.get(lambda x:x.name =='CA')
                if len(atm)==0:
                    raise RuntimeError("No CA atoms in %s"%mol.name) 
                if sticksBallsLicorice == 'Licorice':
                    drawCABalls(mol, atm, only, 0, negate, cradius, 0,
                                cquality, setScale)
                elif sticksBallsLicorice == 'Sticks only':
                    drawCABalls(mol, atm, only, 1, negate, bRad, bScale,
                                bquality, setScale)
                else:
                    drawCABalls(mol, atm, only, 0, negate, bRad, bScale,
                           bquality, setScale)
                drawCASticks(mol, atm, only, negate, cradius, cquality, setScale)
                self.app()._executionReport.addSuccess('displayed backbone trace for molecule %s successfully'%
                    mol.name, obj=atm)
            except:
                if self.app().trapExceptions is False:
                    exc_info = sys.exc_info()
                    raise exc_info[1], None, exc_info[2]
                else:
                    msg = 'Error while displaying backbone trace(%s) for molecule %s' %\
                          (sticksBallsLicorice, mol.name)
                    self.app().errorMsg(sys.exc_info(), msg, obj=atm)
            
        if only and len(self.app().Mols) != len(molecules):
            # if only is True we need to undisplay back bone for molecules
            # that are not included in nodes
            mols = self.app().Mols - molecules
            for mol in mols:
                try:
                    negate = True
                    only=0
                    if sticksBallsLicorice == 'Licorice':
                        drawCABalls(mol, mol.allAtoms, only, 0, negate, cradius,
                              0, cquality, setScale)
                    elif sticksBallsLicorice == 'Sticks only':
                        drawCABalls(mol, mol.allAtoms, only, 1, negate, bRad,
                              bScale, bquality, setScale)
                    else:
                        drawCABalls(mol, mol.allAtoms, only, 0, negate, bRad,
                              bScale, bquality, setScale)
                    drawCASticks(mol, mol.allAtoms, only, negate, cradius, cquality, setScale)
                except:
                    if self.app().trapExceptions is False:
                        exc_info = sys.exc_info()
                        raise exc_info[1], None, exc_info[2]
                    else:
                        msg = 'Error while undisplaying backbone trace for molecule %s to enforce only=True option'%mol.name
                        self.app().errorMsg(sys.exc_info(), msg, obj=mol) 
        if self.createEvents:
            event = EditGeomsEvent('trace', [nodes,[only, negate, bRad,bScale, bquality, 
						cradius, cquality,sticksBallsLicorice]])
            self.app().eventHandler.dispatchEvent(event)



class UndisplayBackboneTrace(DisplayCommand):
    """ The UndisplayBackboneTrace command is an interactive command to undisplay part of the molecule when represented as sticks and balls.
    \nPackage : PmvApp
    \nModule  : displayCmds
    \nClass   : UnDisplayBackboneTrace
    \nCommand : undisplayBackboneTrace
    \nSynopsis:\n
        None <- undisplayBackboneTrace(nodes, **kw)\n
        nodes --- any set of MolKit2.Selection describing molecular components\n
        keywords --- undisplay, lines\n 
    """

    def onAddCmdToApp(self):
        DisplayCommand.onAddCmdToApp(self)
        if not self.app().commands.has_key('displayBackboneTrace'):
            self.app().lazyLoad('displayCmds', commands=['displayBackboneTrace'],
                                package='PmvApp')


    def checkArguments(self, nodes, **kw):
        """None <- undisplayBackboneTrace(nodes, **kw)
        nodes: TreeNodeSet holding the current selection"""
        
        kw['negate']=1
        return self.app().displayBackboneTrace.checkArguments(nodes, **kw)
 

    def doit(self, nodes, **kw):
        self.app().displayBackboneTrace(nodes, **kw)


class DisplayBoundGeom(DisplayCommand):
    """Command to display/undisplay geometries that were bound to molecules with
    'bindGeomToMolecularFragment' command. """
    
    def checkDependencies(self, vf):
        from bhtree import bhtreelib

    def checkArguments(self, molSel, geomName, negate=False, nbVert=3, redraw=True):
        """
           molSel : prody selection
           negate : flag when set to 1 undisplay the current selection
           nbVert : number of vertices per triangle needed to select a triangle"""
        assert negate in [True, False, 0, 1]
        assert isinstance(nbVert, int)
        assert isinstance(geomName, str)
        kw = {"negate": negate, "nbVert": nbVert, "redraw": redraw}
        return (molSel, geomName), kw
        
    def doit(self, molSel, geomName, negate=False, nbVert=3, redraw=True):
        mol = molSel.getAtomGroup().getMolecule()
        # get all objects(geometries) that were bound with  self.app().bindGeomToMolecularFragment() method:
        # find geoms that were bound to molecule mol:
        geomC = mol.geomContainer
        geom = geomC.geoms.get(geomName, None)
        if not geom:
            msg = "%s: % geometry container does not contain %s\n" (self.name, mol.name, geomName) 
            raise RuntimeError(msg)
        if not hasattr(geomC, "boundGeom"):
            msg = "%s geometry is not bound to molecule %s\n" (geomName, mol.name)
            raise RuntimeError(msg)
        if not geomC.boundGeom.has_key(geomName):
            msg = "%s geometry is not bound to molecule %s\n" (geomName, mol.name)
            raise RuntimeError(msg)
        # get the set of atoms for the geometry
        if len(molSel) == len(mol.select()):
            if negate:
                _set = mol.emptySelection()
            else:
                _set = mol.select()
        else:
            _set = geomC.atoms[geomName]
            #if negate, remove current atms from displayed _set
            if negate:
                _set = _set - molSel
            else:
                _set = molSel|_set
        geomC.setAtomsForGeom(geomName, _set)
        
        if len(_set) == 0:
            geom.Set(visible=0, tagModified=False)
            return
        # get atom indices
        atomindices = []
        cl_atoms = geomC.boundGeom[geomName]['cl_atoms']
        vs = []
        mol_lookup = geomC.boundGeom[geomName]['mol_lookup']
        for i in _set.getIndices():
            ind = mol_lookup.get(i, None)
            if ind is not None:
                atomindices.append(ind)
        cond = [x in atomindices for x in cl_atoms]
        #a new subset of vertex indices:
        nvs = numpy.nonzero(cond)[0] #returns indices of cond,
                                    #where its elements != 0

        #print "atom inds to display: " , nvs
        faces = geomC.boundGeom[geomName]['fs']
        norms = geomC.boundGeom[geomName]['fns']
        s = faces.shape
        # find a subset of faces and face normals:
        from bhtree import bhtreelib
        from time import time
        t1=time()
        nvs = nvs.astype('i')
        nfs_ind = bhtreelib.findFaceSubset(nvs, faces, nbVert)#indices of
                                                              #subset of faces
        nfs = numpy.take(faces, nfs_ind, axis=0)
        nfns = numpy.take(norms, nfs_ind, axis=0)
        t2=time()
        #print "time to loop: ", t2-t1
        #print "nfs.shape: ", nfs.shape
        if len(atomindices)==0:
            geom.Set(visible=0, tagModified=False)
        else:
            col = self.getVertColors(geomName, mol)
            geom.Set(faces=nfs, fnormals = nfns, materials=col,
                  inheritMaterial=False, visible=1, tagModified=False)
            # update texture coordinate if needed
            if geom.texture and g.texture.enabled and g.texture.auto==0:
                mol.geomContainer.updateTexCoords[o](mol)
        if self.createEvents: # and len(rsetOn)+len(rsetOff):
            event = EditGeomsEvent('boubdGeom_ds',
                                   [molSel,[negate, nbVert]],
                                   setOn=[], setOff=[])
            self.app().eventHandler.dispatchEvent(event)

    def getVertColors(self, geomName, mol):
        """Function called to map atomic properties to the vertices of the
        geometry"""
        geomC = mol.geomContainer
        geom = geomC.geoms[geomName]
        atoms = geomC.atoms[geomName]
        surfAtomInds = geomC.boundGeom[geomName]['atoms'].getIndices()
        #print "coarseMSInds:", surfinds
        # array of colors of all atoms for the msms.
        from DejaVu2.viewerConst import OVERALL, PER_VERTEX
        from opengltk.OpenGL import GL
        dataLabel = geomName+"colors"
        if dataLabel.find("-") != -1:
           dataLabel = dataLabel.replace("-", "")
        if not dataLabel in mol._ag.getDataLabels():
            # use colors of the lines
            lines = geomC.geoms['singleBonds']
            if lines.materials[GL.GL_FRONT].binding[1]==OVERALL:
                col = lines.materials[GL.GL_FRONT].prop[1][:,:3]
                col = numpy.resize( col, (len(mol._ag), 3) )
            elif lines.materials[GL.GL_FRONT].binding[1]==PER_VERTEX:
                col = lines.materials[GL.GL_FRONT].prop[1][:,:3]
            mol._ag.setData(dataLabel, col)
        else:
            col = mol._ag.getData(dataLabel)
        #col - is an array of all molecule atoms colors. 
        surfCol = col[surfAtomInds] # colors of surface atoms
        cl_atoms = geomC.boundGeom[geomName]['cl_atoms']
        return numpy.take(surfCol, cl_atoms, axis=0).astype('f')

    def refreshDisplay(self, mol):
        #import pdb;pdb.set_trace()
        atomSets = mol.geomContainer.atoms
        name = 'coarseMolSurf-%s'%(mol._basename)
        atoms = (atomSets.get(name, mol.emptySelection()))
        atomSets[name] = mol.emptySelection()
        if atoms:
            self(atoms, surfName=name)
    

class UndisplayBoundGeom(DisplayBoundGeom):
    
    def checkArguments(self, nodes, **kw):
        """None <- undisplayBoundGeom(nodes, **kw)
           nodes  : TreeNodeSet holding the current selection (mv.activeSelection.get())
           """
        kw['negate']= 1
        return self.app().displayBoundGeom.checkArguments(nodes, **kw)
    
    def doit(self, nodes, **kw):
        kw['negate']= 1
        return self.app().displayBoundGeom(nodes, **kw)



class ShowMolecules(MVCommand):
    """The showMolecules command allows the user to show or hide chosen molecules. \n
    Package : PmvApp \n
    Module  : displayCmds \n
    Class   : ShowMolecules \n
    Command : showMolecules \n
    Synopsis:\n
        None <--- showMolecules(molName, show = True, **kw)\n
        molName --- list of the string representing the name of the molecules to be hidden or shown\n
        show --- Boolean Flag when True the molecules corresponding to the given names are hidden, when set to 0 they will be shown\n
        keywords --- show hide molecules\n
    Events: generates ShowMoleculesEvent
    """
    def expandArg0(self, obj):
        if isinstance(obj, (list, tuple)):
            res = SelectionSet()
            for o in obj:
                res.extend(MVCommand.expandArg0(self,o))
            return res
        else: return MVCommand.expandArg0(self,obj)

    def onRemoveObjectFromViewer(self, obj):
        if hasattr(obj, 'displayMode'):
            delattr(obj, 'displayMode')

    def initializeMolForCmd(self, patoms):
        self.initializedFor[patoms.getMolecule()] = True

    def doit(self, selection, show):
        mol = selection.getAtomGroup().getMolecule()
        self.app().pushUndoCmd(self, (selection,), {'show': not show})
        mol.geomContainer.geoms['master'].Set(visible = show,
                                              tagModified=False)

        if self.createEvents:
            event = ShowMoleculesEvent( mol, visible = show)
            self.app().eventHandler.dispatchEvent(event)

    def checkArguments(self, selectionSet, show=True):
        """None <- showMolecules(molName, show=True) \n
           molName : list of string name of molecules (mol.name) \n
           show  : flag when set to 1 hide the given molecule
        """
        # can not do "assert molName"  here,
        # it will fail when molName is  <MoleculeSetNoSelection instance> empty
        # (show molecule from "All Molecules" menu in the dashboard)  
        # In this case doit has to return 
        assert show in [True, False, 1, 0]
        kw = {'show' : show}
        return (selectionSet,), kw



commandClassFromName = {
    'displayLines' : [DisplayLines, None],
    'undisplayLines' : [UndisplayLines,  None],
    'displayCPK': [DisplayCPK, None],
    'undisplayCPK' : [UndisplayCPK,  None],
    'displaySB': [DisplaySticksAndBalls, None],
    'undisplaySB' : [UndisplaySticksAndBalls, None],
    
    ## 'displayBackboneTrace' : [DisplayBackboneTrace, None],
    ## 'undisplayBackboneTrace': [UndisplayBackboneTrace, None],
     'displayBoundGeom': [DisplayBoundGeom, None],
    ## 'undisplayBoundGeom': [UndisplayBoundGeom, None],
    'showMolecules' : [ShowMolecules, None],
    }


def initModule(viewer, gui=True):
    for cmdName, values in commandClassFromName.items():
        cmdClass, guiInstance = values
        viewer.addCommand(cmdClass(), cmdName, guiInstance)
