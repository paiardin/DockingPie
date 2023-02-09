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
# Copyright: M. Sanner TSRI 2013
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/PmvApp/Pmv.py,v 1.34.4.3 2017/10/12 00:21:33 annao Exp $
#
# $Id: Pmv.py,v 1.34.4.3 2017/10/12 00:21:33 annao Exp $
#

"""
define a class for a molecular viewer
"""
import os, sys, weakref, tempfile, tarfile, shutil, threading
from time import time

import prody, openbabel, numpy
from numpy import savez, load
from prody.atomic.functions import saveAtoms, loadAtoms
from glob import glob

from AppFramework.App import AppFramework, GeomContainer, AddGeometryEvent, DeleteObjectEvent
from AppFramework.AppCommands import CmdLoader

from PmvApp.group import Group
from MolKit2.molecule import Molecule, MoleculeSet, Atom, Chain, Residue
from MolKit2.selection import Selection as MolKitSelection
from MolKit2.selection import Selection, SelectionSet, selector
from MolKit2.openBabelInterface import OBMolToPrody, ProdyToOBMol

from mglutil.util.callback import CallbackFunction
from mglutil.events import Event

def reduceName(name, length):
    """ turn the string name into 'start of name ... end of name' if length of name
is greater than length
"""
    if len(name) < length: return name
    return name[:length/2]+'...'+name[-length/2:]

def formatName(repStr, maxLength):
    """format the string name so that no line exceeds maxLength characters
    """
    if len(repStr)>maxLength:
        pt = 0
        frepStr = ""
        while pt < len(repStr):
            frepStr += repStr[pt:min(pt+maxLength, len(repStr))]+'\n'
            pt += maxLength
        frepStr = frepStr[:-1]
    else:
        frepStr = repStr
    return frepStr
    
numOfSelectedVerticesToSelectTriangle = 1 # 1 , 2 or 3

class RenameSelectionEvent(Event):
    # event.item item associate with event.selection
    # event.selection selection that was renamed
    # event.setCurrent True if curSelection was renamed and was active
    pass

class RenameGroupEvent(Event):
    # event.item item associate with event.object
    # event.object group that was renamed
    pass
    
class RenameTreeNodeEvent(Event):
    # event.object TreeNodeE that was renamed
    pass
    
class AddGroupEvent(Event):
    # event is created when group is added to Pmv
    # event.group is the group that was added
    pass

class DeleteGroupsEvent(Event):
    # event is created when group is added to Pmv
    # event.group is the group that was deleted
    pass

class ReparentGroupObject(Event):
    # event is created when a group or a molecule is move to a new parent
    # event.newGroup is the group which is the new parent of event.object
    # event.oldGroup is the group which was the old parent of event.object
    # event.object is the object reparented
    pass

class ActiveSelectionChangedEvent(Event):
    # event is created when pmv's active selection object is replaced by a new one
    # event.current is the new Selection object
    # event.previous is the old Selection object
    pass

class DeleteNamedSelectionEvent(Event):
    # event is created when a named selection is deleted
    # event.selection the deleted Selection object
    pass


class BeforeAddMoleculeEvent(Event):
    # event is created at the begging of the addition of an molecule
    #
    # event.name is the name of the molecule added
    # event.object is the molecule being added
    pass

## class AddMoleculeCmdEvent(Event):
##     # event is created each time we call a command to be carried out on
##     #  an molecule that is added
##     #
##     # event.number number of command currently called
##     # event.total total number of commands to be called
##     # event.cmdName name of the command currently called
##     pass


class AfterAddMoleculeEvent(Event):
    # event is created after the molecule has been added
    # 
    # event.object is the molecule that has been added
    pass


class BeforeDeleteMoleculeEvent(Event):
    # event is created before the molecule is deleted
    # 
    # event.object is the molecule that is being deleted
    pass

class AfterDeleteMoleculeEvent(Event):
    # event is created after the molecule is deleted
    # 
    # event.object the molecule that has been deleted
    pass


class EditAtomsEvent(Event):
    # event is created when atom properties are modified
    # event.objects - AtomSet
    pass
    #def __init__(self, fields=[], atoms=None, *args, **kw):
    #    Event.__init__(self, *args, **kw)
    #    self.fields = objects
    #    self.name = name

class DeleteAtomsEvent(Event):
    pass

class AfterDeleteAtomsEvent(Event):
    pass

class AddAtomsEvent(Event):
    pass

class MoveAtomsEvent(Event):
    pass

class ShowMoleculesEvent(Event):
    pass

class AddGeomsEvent(Event):
    pass

class DeleteGeomsEvent(Event):
    pass

class EditGeomsEvent(Event):
     def __init__(self, name=None, objects=[], *args, **kw):
        Event.__init__(self, *args, **kw)
        self.objects = objects
        self.name = name

class BondAndAtomTypeChangeEvent(Event):
    
    def __init__(self, molecule=None, bondType='double',
                 vertices=[], colors=[], **kw):
        Event.__init__(self, (), **kw)
        self.molecule = molecule
        self.bondType=bondType
        self.vertices=vertices
        self.colors=colors

class RefreshDisplayEvent(Event):
    pass

moleculeEvents = {
    'beforeAddObject' : BeforeAddMoleculeEvent,
    'afterAddObject' : None,#AfterAddMoleculeEvent,
    'beforeDeleteObject' : BeforeDeleteMoleculeEvent,
    'afterDeleteObject' : AfterDeleteMoleculeEvent
    }

from MolKit2.selection import Selection as MolKitSelection

class PmvSelection(SelectionSet):
    # class to store a selection in Pmv

    def __init__(self, app, name='selection', itr=[]):
        SelectionSet.__init__(self, itr=itr, name=name)
        self.app = app
        
    ## def copy(self):
    ##     newSelSet = PmvSelection(self.app, self.name)
    ##     for sel in self:
    ##         newSelSet.append(sel.copy())
    ##     return newSelSet

##     def empty(self):
##         if len(self.atoms): return False
##         return True
    
##     def __sub__(self, other):
##         sel = self.copy()
##         if isinstance(obj, Selection):
##             sel.atoms -= other.atoms
##         else:
##             sel.atoms -= other.findType(Atom)
##         sel.update()
##         #print 'asda', self.children, len(self.children)
##         return sel
    
##     def __add__(self, other):
##         sel = self.copy()
##         if isinstance(obj, Selection):
##             sel.atoms += other.atoms
##         else:
##             sel.atoms += other.findType(Atom)
##         sel.update()
##         #print 'asda', self.children, len(self.children)
##         return sel
    
##     def __or__(self, other):
##         sel = self.copy()
##         if isinstance(obj, Selection):
##             sel.atoms = sel.atoms | other.atoms
##         else:
##             sel.atoms = sel.atoms | other.findType(Atom)
##         sel.update()
##         #print 'asda', self.children, len(self.children)
##         return sel
    
##     def __inter__(self, other):
##         sel = self.copy()
##         if isinstance(obj, Selection):
##             sel.atoms = sel.atoms.__inter__(other.atoms)
##         else:
##             sel.atoms = sel.atoms.__inter__(other.findType(Atom))
##         sel.update()
##         #print 'asda', self.children, len(self.children)
##         return sel
    
##     def __xor__(self, other):
##         sel = self.copy()
##         if isinstance(obj, Selection):
##             sel.atoms = sel.atoms.__xor__(other.atoms)
##         else:
##             sel.atoms = sel.atoms.__xor__(other.findType(Atom))
##         sel.update()
##         #print 'asda', self.children, len(self.children)
##         return sel
    
##     def __and__(self, other):
##         sel = self.copy()
##         if isinstance(obj, Selection):
##             sel.atoms = sel.atoms & other.atoms
##         else:
##             sel.atoms = sel.atoms & other.findType(Atom)
##         sel.update()
##         #print 'asda', self.children, len(self.children)
##         return sel
    
##     def get(self, klass=Atom):
##         if klass==Atom:
##             return self.atoms.copy()
##         elif klass==Residue:
##             return self.atoms.parent.uniq()
##         elif klass==Chain:
##             return self.atoms.parent.parent.uniq()
##         elif klass==Molecule or klass==Protein:
##             return self.atoms.top.uniq()
        
##     def set(self, obj):
##         #print 'setting selection', self.name, str(obj)
##         if isinstance(obj, AtomSet):
##             self.atoms = obj
##         else:
##             self.atoms = obj.findType(Atom)
##         self.atoms = self.atoms.uniq()
##         #print 'asda', self.children, len(self.children)
##         self.update()
    
    ## def isSelected(self, obj):
    ##     # returns True if all atoms in obj are in this PmvSelection
    ##     #         False is all atoms in obj are deselected
    ##     #         'partial' else

    ##     #if not isinstance(obj, AtomSet):
    ##     #    atoms = obj.findType(Atom)
    ##     #else:
    ##     #    atoms = obj

    ##     if len(obj)==0: return False
        
    ##     # selection status of first atom
    ##     indices = obj.getIndices()
    ##     sel0 = self.atomsDict.has_key(indices[0])
    ##     for a in indices[1:]:
    ##         sel = self.atomsDict.has_key(a)
    ##         if sel!=sel0:
    ##             return 'partial'
    ##     return sel0

##     def isNotSelected(self, obj, mode='fast'):
##         if not isinstance(obj, AtomSet):
##             atoms = obj.findType(Atom)
##         else:
##             atoms = obj
##         if mode=='fast':
##             for a in atoms:
##                 d = self.atomsDict.get(a, None)
##                 if d is None:
##                     return True
##             return False
##         elif mode=='complete':
##             raise
##         else:
##             raise RuntimeError("Bad mode %s, expected 'fast' or 'complete'"%mode)

##     def isAnySelected(self, obj, mode='fast'):
##         if not isinstance(obj, AtomSet):
##             atoms = obj.findType(Atom)
##         else:
##             atoms = obj
##         if mode=='fast':
##             for a in atoms:
##                 d = self.atomsDict.get(a, None)
##                 if d is not None:
##                     return True
##             return False
##         elif mode=='complete':
##             raise
##         else:
##             raise RuntimeError("Bad mode %s, expected 'fast' or 'complete'"%mode)

##     def add(self, obj):
##         #print 'setting selection', self.name, str(obj)
##         if isinstance(obj, Selection):
##             self.atoms += obj.atoms
##         else:
##             self.atoms += obj.findType(Atom)
##         self.atoms = self.atoms.uniq()
##         #print 'asda', self.children, len(self.children)
##         self.update()

##     def remove(self, obj):
##         #print 'setting selection', self.name, str(obj)
##         if isinstance(obj, Selection):
##             self.atoms -= obj.atoms
##         else:
##             self.atoms -= obj.findType(Atom)
##         #print 'asda', self.children, len(self.children)
##         self.update()

##     def xor(self, obj):
##         #print 'setting selection', self.name, str(obj)
##         if isinstance(obj, Selection):
##             self.atoms = self.atoms ^ obj.atoms
##         else:
##             self.atoms = self.atoms ^ obj.findType(Atom)
##         #print 'asda', self.children, len(self.children)
##         self.update()

##     def inter(self, obj):
##         #print 'setting selection', self.name, str(obj)
##         if isinstance(obj, Selection):
##             self.atoms = self.atoms & obj.atoms
##         else:
##             self.atoms = self.atoms & obj.findType(Atom)
##         #print 'asda', self.children, len(self.children)
##         self.update()
        
## from MolKit2.listSet import ListSet
## class SelectionSet(ListSet):

##     def __init__(self, data=None, elementType=Selection, **kw):
##         ListSet. __init__(self, data=data, elementType=Selection)

        
class MolApp(AppFramework):
    """
    Application that supports lazy loading of commands
    """

    def exit(self):
        AppFramework.exit(self)

    def handleDeleteObject(self, event):
        obj = event.object
        if isinstance(obj, Molecule):
            obj._colors = {}
            for label in obj._ag.getDataLabels(which='user'):
                if label.startswith('colorsIndices_'):
                    obj._ag.setData('colorsIndices_lines', numpy.zeros(len(obj._ag), 'i'))
            # not sure about mol._renderindProp
            # del mol._renderindProp

    def handleAddAtoms(self, event):
        mol = event.molecule
        gc = mol.geomContainer
        # update the coordinate stored in geomContainer
        gc.allCoords = mol._ag._coords[mol._ag._acsi].astype('f')
        # update the _ag in the atoms sets in geomContainer
        for atomset in gc.atoms.values():
            atomset._ag = mol._ag
        coords = mol._ag.getCoords().astype('f')
        # fix geometry vertices
        for geom in gc.geoms.values():
            if hasattr(geom, '_hasAllVertices') and geom._hasAllVertices:
                geom.vertexSet.vertices.array = coords

    def handleMoveAtoms(self, event):
        mol = event.molecule
        gc = mol.geomContainer
        gc.allCoords[:] = mol._ag.getCoords().astype('f')

    #def addCommand(self, command, name):
    #    AppFramework.addCommand(self, command, name)
        #cb = CallbackFunction( self.handleDeleteObject, command)
        #self.eventHandler.registerListener(DeleteObjectEvent, cb)
        
    def __init__(self, title="Molecule Viewer", eventHandler=None):

        AppFramework.__init__(self, name=title, eventHandler=eventHandler)

        self.Mols = MoleculeSet('All Molecules')
        self.lastMolNumer = 0 # used to give molecules unique numbers
        
        #self.sets = Sets()  # store user-defined sets in this dict
        self.grids3D = {} # list of Grid3D objects storing volumetric data
        #self.allAtoms = AtomSet() # set of all atoms (across molecules)

        ##
        ## opened session is a list of folders in which we untarred the session
        ## but has to keep the session open for instance because it contains a
        ## multimolecule object we are indexing into
        ## this list of folder will be deleted upon the termination of the app
        self._openedSessions = []

        ## list of running threads
        self._threads = []
        
        def moleculeValidator(mol):
            return isinstance(mol, Molecule)
          
        self.addObjectType("Molecule", moleculeValidator, moleculeEvents)
        self.objects["Molecule"] = self.Mols

        ## selections
        ##
        self.curSelection = PmvSelection(self, u'Current Selection')
        # active selection is the selection to which pmv.select will add atoms
        # activeSelection always hold a Selection which is either
        # self.curSelection or a named selection
        self.activeSelection = self.curSelection
        self.namedSelections = {}
        #self.selectionLevel = Atom

        ##
        ## register listeners for events
        self.registerlisteners()
        self.eventHandler.registerListener(AddAtomsEvent,
                                           self.handleAddAtoms)
        self.eventHandler.registerListener(MoveAtomsEvent,
                                           self.handleMoveAtoms)
        self.eventHandler.registerListener(DeleteObjectEvent,
                                           self.handleDeleteObject)
 
        ## groups
        ##
        self.groups = {} # name: MoleculeGroup Instance
        
        # this levelColors attribute used to be in ICmdCaller (???)
        # the level variable is used by selection commands
        self.levelColors = {
            'Atom':(1.,1.,0.),
            'Residue':(0.,1.,0.),
            'Chain':(0.,1.,1.),
            'Molecule':(1.,0.,0.),
            'Protein':(1.,0.,0.),
            }
        choices = ['molecules', 'chains', 'conformations']
        choice = 'molecules'
        self.userpref.add('Save session multiMolSize', 50, category="Pmv",
                          doc = """define the threshold (in MB) above which the user is asked to decided whether or not to include multi-molecules into session files""")
        self.userpref.add('Read molecules as', choice, validValues=choices,
                          category="Molecules",
                          doc = """for pdb file with multi MODEL, can be read as 'molecules', 'chains', 'conformations'.chains' is not yet implemented (7/09)""")

        choices = ['caseSensitive', 'caseInsensitive',
                   'caseInsensWithEscapedChars']

        self.userpref.add('String Matching Mode', 'caseSensitive', validValues=choices,
                          doc = """When set to caseSensitive the string match
mode will be case sensitive the other possibility is to be case insensitive or
case insensitive with escaped characters.
""")
        # overwrite firstObject only with firstMoleculeOnly
        self.userpref.add(
            'Center Scene','firstMoleculeOnly',
            validValues=['firstMoleculeOnly','always', 'never', 'ask'],
            category="3DViewer",
            doc="""the value of this preference defines whether\
the sceen is centered and tranformed to fit into the field\
of view when a new molecule is loaded""")
        #self.userpref['Center Scene']['validValues'][0] = 'firstMoleculeOnly'
        #self.userpref.set('Center Scene', 'always')

        self.gui = None # will be set if this gets bound to a GUI

        self.loadBuiltinCmds()

    def cleanup(self):
        """close the application after cleaning up"""
        for folder in self._openedSessions:
            try:
                shutil.rmtree(folder)
            except OSError:
                pass
        self._openedSessions = []

    def checkForThreads(self):
        #check for running threads
        alive = []
        _threads = []
        for thread, descr in self._threads:
            if thread.isAlive():
                alive.append(descr)
                _threads.append( (thread, descr) )
        self._threads = _threads
        return alive

    def close(self):
        alive = self.checkForThreads()
        if len(alive)>0:
            print "the following threads are still running:"
            for d in alive: print "    ",d
        print "Waiting for these threads to finish"
        self.cleanup()
        sys.exit(1)

    def runInThread(self, func, description, args=(), kw={}):
        """run fun(*args, **kw) in a thread and keep track of the threads and descriptions"""
        t1 = threading.Thread(target=func, args=args, kwargs=kw)
        t1.start()
        self._threads.append((t1, description))

    def loadBuiltinCmds(self):
        ## the most basic commands are always loaded
        from .bondsCmds import BuildBondsByDistance
        self.addCommand(BuildBondsByDistance(), 'buildBondsByDistance')

        from .fileCmds import MoleculesReader, Fetch
        self.addCommand(MoleculesReader(), 'readMolecules')
        self.addCommand(Fetch(), 'fetch')

        from .displayCmds import DisplayLines, UndisplayLines, \
             DisplayCPK, UndisplayCPK, \
             DisplaySticksAndBalls, UndisplaySticksAndBalls,\
             ShowMolecules
        self.addCommand(DisplayLines(), 'displayLines')
        self.addCommand(UndisplayLines(), 'undisplayLines')
        self.addCommand(DisplayCPK(), 'displayCPK')
        self.addCommand(UndisplayCPK(), 'undisplayCPK')
        self.addCommand(DisplaySticksAndBalls(), 'displaySB')
        self.addCommand(UndisplaySticksAndBalls(), 'undisplaySB')
        self.addCommand(ShowMolecules(), 'showMolecules')
        #self.addCommand((), '')

        #from .msmsCmds import  DisplayMSMS, UndisplayMSMS
        #self.addCommand(DisplayMSMS(), 'displayMSMS')
        #self.addCommand(UndisplayMSMS(), 'undisplayMSMS')
        self.lazyLoad('msmsCmds', package='PmvApp')
        self.displayMSMS.loadCommand()
        self.undisplayMSMS.loadCommand()

        from .hbCmds import DisplayHB
        self.addCommand(DisplayHB(), 'displayHB')
        from .labelCmds import LabelAtoms, LabelResidues,\
              UnlabelAtoms, UnlabelResidues
        self.addCommand(LabelAtoms(), 'labelAtoms')
        self.addCommand(LabelResidues(), 'labelResidues')
        self.addCommand(UnlabelAtoms(), 'unlabelAtoms')
        self.addCommand(UnlabelResidues(), 'unlabelResidues')
        
        self.lazyLoad('editCmds', package='PmvApp')
        self.assignAtomsRadii.loadCommand()

        from .colorCmds import ColorCommand, ColorByAtomType, ColorByMolecule,\
             ColorByChain, ColorRainbow, CustomColor, MoleculesRainbow,\
             ColorBySecondaryStructure, ColorRainbowChain
        self.addCommand(ColorCommand(), 'color')
        self.addCommand(ColorByAtomType(), 'colorByAtomType')
        self.addCommand(ColorByMolecule(), 'colorByMolecules')
        self.addCommand(ColorByChain(), 'colorByChains')
        self.addCommand(ColorRainbow(), 'colorRainbow')
        self.addCommand(ColorRainbowChain(), 'colorRainbowChain')
        self.addCommand(CustomColor(), 'customColor')
        self.addCommand(MoleculesRainbow(), 'moleculesRainbow')
        self.addCommand(ColorBySecondaryStructure(), 'colorBySS')
        #self.addCommand((), '')

        from .selectionCmds import MVSelectCommand, MVDeSelectCommand,\
             MVClearSelection, MVExpandSelection, MVSelectAround,\
             MVInvertSelection
        self.addCommand(MVSelectCommand(), 'select')
        self.addCommand(MVDeSelectCommand(), 'deselect')
        self.addCommand(MVClearSelection(), 'clearSelection')
        self.addCommand(MVExpandSelection(), 'expandSelection')
        self.addCommand(MVSelectAround(), 'selectAround')
        self.addCommand(MVInvertSelection(), 'invertSelection')
        #self.addCommand((), '')

        from .cartoonCmds import ComputeCartoon, DisplayCartoon, UndisplayCartoon
        self.addCommand(ComputeCartoon(), 'computeCartoon')
        self.addCommand(DisplayCartoon(), 'displayCartoon')
        self.addCommand(UndisplayCartoon(), 'undisplayCartoon')

        from .deleteCmds import DeleteMolecule, DeleteAtoms
        self.addCommand(DeleteMolecule(), 'deleteMolecule')
        self.addCommand(DeleteAtoms(), 'deleteAtoms')
    ##
    ## event listeners
    ##
    def registerlisteners(self):
        ev = self.eventHandler
        ev.registerListener(DeleteAtomsEvent, self.handleDeleteAtomsEvent)

    def handleDeleteAtomsEvent(self, event):
        # remove atoms from the names selections if necessary
        #selection, indMap = event.objects
        selection = event.deletedAtoms
        from MolKit2.selection import SelectionSet
        from PmvApp.selectionCmds import SelectionEvent
        atSet = SelectionSet([selection], "")
        selectedAtoms = self.curSelection.copy()
        if self.curSelection.isAnySelected(atSet):
            if self.curSelection == self.activeSelection:
                d1 = d2 = None
                if hasattr(self.curSelection, '_treeItems'):
                    d1 = self.curSelection._treeItems
                if hasattr(self.activeSelection, '_treeItems'):
                    d2 = self.activeSelection._treeItems
                self.curSelection = self.curSelection - atSet
                self.activeSelection = self.curSelection 
                if d1:
                    self.curSelection._treeItems = d1
                if d2:
                    self.activeSelection._treeItems = d2
            else:
                self.curSelection = self.curSelection - atSet

            for ag, sel in self.curSelection._agDict.items():
                if sel == selection:
                    self.curSelection._agDict.pop(ag)
                    break
            event = SelectionEvent(object=self.activeSelection,
                new=self.curSelection, old=selectedAtoms,
                setOn=self.curSelection-selectedAtoms,
                setOff=selectedAtoms-self.curSelection)
            self.eventHandler.dispatchEvent(event)

        for name, sel in self.namedSelections.items():
            if sel.isAnySelected(atSet):
                sel = sel - atSet
                event = SelectionEvent(object=self.activeSelection,
                    new=self.curSelection, old=selectedAtoms,
                    setOn=sel-selectedAtoms, setOff=selectedAtoms-sel)
                self.eventHandler.dispatchEvent(event)
        
    ##
    ## Renaming objects
    ##
    def rename(self, obj, name, item=None):
        """Rename Pmv objects
None <- pmv.rename(obj, item, name)
obj: object or name of an object (Objects can be selections, groups or molecular fragments
name: usesr defined name for the object
items: is the item associate with the object in the Dashboard, or None

For groups and selections the name will overwrite the old name. For molecular fragments
the name will become the alias for this object
"""
        inUndo = self.undo.inUndo
        inRedo = self.redo.inUndo
        if isinstance(obj, str):
            # turn the string into an object
            if obj == 'Current Selection':
                obj = self.curSelection
            elif obj in self.namedSelections.keys():
                obj = self.namedSelections[obj]
            elif obj in self.groups.keys():
                obj = self.groups[obj]
            else:
                objs = self.expandNodes(obj)
                if len(objs)==0:
                    raise ValueError, "%s is not a selection, group or molecular fragment"%obj
                obj = objs

        if isinstance(obj, SelectionSet):
            # make sure the name is not already used
            if obj.name == name: return
            if name=='Current Selection' or name in self.namedSelections.keys():
                raise ValueError, "%s is alreday use for a selection"%name
            if obj == self.curSelection:
                obj.name = name
                # save it in named selections
                self.namedSelections[name] = obj
                self.curSelection = PmvSelection(self, u'Current Selection')
                event = RenameSelectionEvent(item=item, selection=obj,
                                             setCurrent=(obj==self.activeSelection))
                self.eventHandler.dispatchEvent(event)
                def mkSelCurSel(selection, item=None):
                    self.curSelection = selection
                    old = selection.name
                    del self.namedSelections[selection.name]
                    selection.name = u'Current Selection'
                    if self.undo.inUndo>=0: # we are in undo mode 
                        self.undo._cmdList[0].extend([(self.rename, (obj, old), {'item':item})]  )
                        if self.undo.inUndo==0:
                            self.redo.addUndoCall(*self.undo._cmdList )
                    event = RenameSelectionEvent(item=item, selection=selection,
                                                 setCurrent=(selection==self.activeSelection))
                    self.eventHandler.dispatchEvent(event)
                if inUndo==-1 and inRedo==-1: # not doing Undo or Redo
                    self.undo.addUndoCall( [(mkSelCurSel, (obj,), {'item':item})],
                                           'rename "Current Selection' )
                else:
                    if inRedo>=0: # we are in redo mode 
                        self.redo._cmdList[0].extend([(mkSelCurSel, (obj,), {'item':item})] )
                        if inRedo==0:
                            self.undo.addUndoCall( *self.redo._cmdList )
                    
            else:
                namedSelection = obj
                old = obj.name
                del self.namedSelections[namedSelection.name]
                namedSelection.name = name
                self.namedSelections[namedSelection.name] = namedSelection

                event = RenameSelectionEvent(item=item, selection=obj)
                self.eventHandler.dispatchEvent(event)
                if inUndo==-1 and inRedo==-1: # not doing Undo or Redo
                    self.undo.addUndoCall( [(self.rename, (obj, old), {'item':item})],
                                           'rename %s'%old )
                else:
                    if inUndo>=0: # we are in undo mode 
                        self.undo._cmdList[0].extend([(self.rename, (obj, old), {'item':item})] )
                        if inUndo==0:
                            self.redo.addUndoCall(*self.undo._cmdList )
                    else: # redo mode
                        self.redo._cmdList[0].extend([(self.rename, (obj, old), {'item':item})] )
                        if inRedo==0:
                            self.undo.addUndoCall( *self.redo._cmdList )

        elif isinstance(obj, Group):
            if obj.name == name: return
            if name in self.groups.keys():
                raise ValueError, "%s is alreday use for a group"%name
            old = obj.name
            self.groups[name] = obj
            del self.groups[obj.name]
            obj.name = name
            event = RenameGroupEvent(object=obj, item=item)
            self.eventHandler.dispatchEvent(event)
            if inUndo==-1 and inRedo==-1: # not doing Undo or Redo
                self.undo.addUndoCall( [(self.rename, (obj, old), {'item':item})],
                                       'rename %s'%old )
            else:
                if inUndo>=0: # we are in undo mode 
                    self.undo._cmdList[0].extend([(self.rename, (obj, old), {'item':item})] )
                    if inUndo==0:
                        self.redo.addUndoCall(*self.undo._cmdList )
                else: # redo mode
                    self.redo._cmdList[0].extend([(self.rename, (obj, old), {'item':item})] )
                    if inRedo==0:
                        self.undo.addUndoCall( *self.redo._cmdList ) 

        elif isinstance(obj, (Chain, Residue)):
            old = obj.alias
            obj.alias = name
            event = RenameTreeNodeEvent(object=obj)
            self.eventHandler.dispatchEvent(event)
            if inUndo==-1 and inRedo==-1: # not doing Undo or Redo
                self.undo.addUndoCall( [(self.rename, (obj, old), {'item':item})],
                                       'rename %s'%obj.name )
            else:
                if inUndo>=0: # we are in undo mode 
                    self.undo._cmdList[0].extend([(self.rename, (obj, old), {'item':item})] )
                    if inUndo==0:
                        self.redo.addUndoCall(*self.undo._cmdList )
                else: # redo mode
                    self.redo._cmdList[0].extend([(self.rename, (obj, old), {'item':item})] )
                    if inRedo==0:
                        self.undo.addUndoCall( *self.redo._cmdList )             

        elif isinstance(obj, Molecule):
            undo = []
            undoName = ""
            old = obj.alias
            obj.alias = name
            event = RenameTreeNodeEvent(object=obj)
            self.eventHandler.dispatchEvent(event)
            #undo.append( (self.rename, (ob, old), {'item':item}) )
            #undoName += ob.name+', '
            #undoName = reduceName(undoName[:-2], 60)
            #if self.undo.inUndo==-1 and self.redo.inUndo==-1: # not doing Undo or Redo
            #    self.undo.addUndoCall( undo, 'rename %s'%undoName) 
        else:
            print 'RENAMING %s not implemented'%( obj.__class__)

    ###
    ### Clone
    ###

    def clone(self, obj, name=None):
        """Creates a molecule from speciied object(obj). Adds it to the application.
        obj - is a Molecule, Selection or SelectionSet """
        if name is None:
            ll = [mol.name.startswith("clone") for mol in self.Mols]
            name = "clone%d"% (sum(ll)+1)
        newmol = None
        if isinstance(obj, Molecule):
            newmol =  obj.clone(name)
        elif isinstance(obj, (Selection, SelectionSet)):
            newmol =  obj.toMolecule(name)
        if newmol:
            newmol.defaultRadii()
            self.addMolecule(newmol, group=None)
            self.applyDefaultCommands(newmol, "Molecule")
            event = AfterAddMoleculeEvent(molecule=newmol)
            self.eventHandler.dispatchEvent(event)
    ##
    ## Groups
    ##
    def _deleteGroupTree(self, obj):
        if isinstance(obj, Molecule):
            self.deleteMolecules(MoleculeSet([obj]))
        else:
            for ob in obj.children:
                self._deleteGroupTree(ob)
            del self.groups[obj.name]
        
    def deleteGroups(self, groups):
        for group in groups:
            self._deleteGroupTree(group)
        event = DeleteGroupsEvent(groups=groups)
        self.eventHandler.dispatchEvent(event)

    def createGroup(self, name, group=None):
        assert name not in self.groups.keys()
        grp = Group(name)
        self.groups[name] = grp
        if group:
            if isinstance(group, str):
                parent = self.groups.get(group, None)
            else:
                assert isinstance(group, Group)
                parent = group
            group.add(grp)
        else:
            parent = None
        event = AddGroupEvent(group=grp, parent=group)
        self.eventHandler.dispatchEvent(event)
        return grp
    
    def reparentObject(self, obj, group):
        assert group is None or isinstance(group, Group)
        if hasattr (obj, "_group"):
            oldGroup = obj._group
            oldGroup._children.remove(obj)
        else:
            oldGroup = None
        if group is not None:
            if group.validChild(obj):
                group.add(obj)
                
        event = ReparentGroupObject(newGroup=group, oldGroup=oldGroup,
                                    object=obj)
        self.eventHandler.dispatchEvent(event)
    
    ## def addObjectToGroup(self, group, obj):
    ##     assert isinstance(group, MoleculeGroup)
    ##     assert isinstance(obj, (Protein, MoleculeGroup))
    ##     group.children.append(obj)
    ##     obj._group = group
    ##     event = AddObjectToGroupEvent(group=group, object=obj)
    ##     self.eventHandler.dispatchEvent(event)
        
    ## def removeObjectFromGroup(self, group, obj):
    ##     assert isinstance(group, MoleculeGroup)
    ##     assert isinstance(obj, (Protein, MoleculeGroup))
    ##     group.children.remove(obj)
    ##     del obj._group
    ##     event = RemoveObjectFromGroupEvent(group=group, object=obj)
    ##     self.eventHandler.dispatchEvent(event)
        
    ##
    ## Selections
    ##
    def setActiveSelection(self, sele, sendEvent=True):
        """Set the pmv's active selection which is the selection on which \
        pmv.select will operate
Events: ActiveSelectionChangedEvent(old=oldSelection, new=newSelection)
        """
        assert isinstance(sele, SelectionSet)
        if sele==self.activeSelection: return
        oldSele = self.activeSelection
        self.activeSelection = sele
        if sendEvent:
            event = ActiveSelectionChangedEvent(old=oldSele, new=sele)
            self.eventHandler.dispatchEvent(event)

    def deleteNamedSelection(self, sele):
        del self.namedSelections[sele.name]
        if sele == self.activeSelection:
            self.clearSelection()
            self.activeSelection = self.curSelection
        event = DeleteNamedSelectionEvent(selection=sele)
        self.eventHandler.dispatchEvent(event)
        
    def setSelection(self, obj):
        """Set the content of the active selection to be alist of atoms
provided by obj
        """
        if isinstance(obj, SelectionSet):
            self.activeSelection.atoms = obj.atoms[:]
            self.activeSelection.update()
        else:
            self.activeSelection.set(obj)

    def indexColors(self, colorList):
        # build an indexed list for a color vector
        colIndexList = []
        colDict = {}
        colList = []
        nbKeys = 0
        for c in colorList:
            key = "%5.3f,%5.3f,%5.3f,%5.3f"%tuple(c)
            if not colDict.has_key(key):
                colDict[key] = nbKeys
                colList.append(c)
                colIndexList.append(nbKeys)
                nbKeys += 1
            else:
                colIndexList.append(colDict[key])
        return colList, colIndexList
    
    def getGeoms(self, obj, visibleOnly=True):
        # return a dictionary with available geometries 
        # for this selectionSet obj
        geoms = {}
        moleculesAndSelections = [[x.getAtomGroup().getMolecule(),x] for x in obj]
        for mol, sele in moleculesAndSelections:
            gc = mol.geomContainer
            for name, ats in gc.atoms.items():
                if gc.geoms.has_key(name):
                    if visibleOnly and not gc.geoms[name].visible:
                        continue
                    if gc.displayedAs([name], sele, 'fast'):
                        if name in ['balls', 'sticks']:
                            geoms['sb'] = gc.geoms[name]
                        elif name in ['noBond', 'singleBonds', 'doubleBonds', 'tripleBonds', 'aromaticBonds']:
                            continue
                        #    geoms['lines'] = gc.geoms[name]
                        elif name.startswith("chain_") and gc.geoms[name].parent.name == "cartoon":
                            geoms["cartoon %s"%mol._basename] = gc.geoms[name]
                        else:
                            geoms[name] = gc.geoms[name]
        return geoms

    def isNotSelected(self, obj):
        return not self.activeSelection.allSelected( obj)


    def setSelLev(self, value):
        if value==Protein: value = Molecule
        assert value in [Molecule, Chain, Residue, Atom]
        self.setSelectionLevel(value)
      
    
    def readMolecule(self, filename, addToRecent=True, group=None, header=True):
        """Reads molecule in the following formats: pdb, ent, pdbq, \
         pqr, mol2, cif.
        Adds the molecule object to the application
        group is a MoleculeGroup instance or a name(string) of existing group, if not None , the molecule is addded
        to the specified group. 
        """
        t0 = time()
        try:
            import resource
            oldMemUse = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        except:
            oldMemUse=None

        from MolKit2 import Read
        mol = Read(filename, group=group, header=header, indexingMode='background')
        if addToRecent and hasattr(self,'recentFiles'):
            self.recentFiles.add(os.path.abspath(filename), 'readMolecule')
        if oldMemUse is not None:
            newMemUse = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            print 'molecule loaded in %f (s) using %f (kb)'%(
                time()-t0, newMemUse-oldMemUse)

        self.addMolecule(mol, group=group, filename=filename)
        self.applyDefaultCommands(mol, "Molecule")
        event = AfterAddMoleculeEvent(molecule=mol)
        self.eventHandler.dispatchEvent(event)
        
        return mol

    def fetchFile(self, name, ext, group=None, addToRecent=True, pdbFolder=None):
        from MolKit2 import readMMTF
        mol = None
        if pdbFolder is None:
            from mglutil.util.packageFilePath import getCacheFolder
            pdbFolder = getCacheFolder()
        if not os.path.exists(pdbFolder):
            os.mkdir(pdbFolder)
        msg = ""
        if ext == "mmtf":
            from MolKit2 import readMMTF
            filename = os.path.join(pdbFolder, name+".mmtf")
            if not os.path.exists(filename):
                from mmtf import  get_url, MMTFDecoder
                from mmtf.api.default_api import ungzip_data
                import msgpack
                from urllib2 import Request, urlopen, URLError, HTTPError 
                import gzip
                from StringIO import StringIO
                #from mmtf import fetch
                #decoder = fetch(name)
                _data = None
                try:
                    url = get_url(name)
                    request = Request(url)
                    request.add_header('Accept-encoding', 'gzip')
                    response = urlopen(request)
                    if response.info().get('Content-Encoding') == 'gzip':
                        _data = response.read()

                except HTTPError as e:
                    msg = 'ERROR: "%s" is not a valid pdb id.\nError code %s'%(name, str(e.code))
                    if hasattr(e, 'reason'):
                        msg = msg + "(%s)" % e.reason
                    if self.trapExceptions: 
                        raise e
                except URLError as e:
                    msg = "ERROR: the server could not be reached.\nPlease check your internet connection.\n"
                    if hasattr(e, 'reason'):
                        msg = msg+"%s" % (e.reason)
                    if self.trapExceptions: 
                        raise e
                if _data is not None:
                    decoder = MMTFDecoder()
                    data = ungzip_data(_data)
                    decoder.decode_data(msgpack.unpackb(data.read()))
                    mol = readMMTF(filename, decoder=decoder)
                    buf = gzip.GzipFile(fileobj=StringIO(_data), mode='rb', compresslevel=9)
                    with open(filename, 'wb') as outfile:
                        for line in buf:
                            outfile.write(line)
            else:
                 mol = readMMTF(filename)
        elif ext == "pdb":
            if os.path.exists(os.path.join(pdbFolder, name+".pdb")):
                filename = os.path.join(pdbFolder, name+".pdb")
            else:
                from prody import pathPDBFolder, fetchPDB
                filename = fetchPDB(name, folder=pdbFolder, compressed=False)
            if filename:
                from MolKit2 import Read
                mol = Read(filename)
            else:
                msg = "Could not fetch %s molecule" % (name+".pdb")
                if self.trapExceptions:
                    raise RuntimeError(msg)
        if mol:
            self.addMolecule(mol, group=group, filename=filename) 
            self.applyDefaultCommands(mol, "Molecule")
            event = AfterAddMoleculeEvent(molecule=mol)
            self.eventHandler.dispatchEvent(event)
            if addToRecent and hasattr(self,'recentFiles'):
                self.recentFiles.add(os.path.abspath(filename), 'readMolecule')
        return mol, msg

    def addMolecule(self, newmol, group=None, filename=None):
        """
        Add a molecule to the application 
        """
        assert len(newmol._ag) > 0, "Empty molecule. No atom record found for %s." % newmol.name
        if newmol._indexingThread:
            self._threads.append(newmol._indexingThread)

        #IN ANY CASE: change any special characters in name to '-'
        spChar=['?','*','.','$','#',':','-',',']        
        for item in spChar:
            newmol.name = newmol.name.replace(item, '_')
        if len(self.Mols) > 0:
            if newmol._basename in [mol._basename for mol in self.Mols]:
                newmol._basename='%s_%d'%(newmol._basename,len(self.Mols))

        g = MolGeomContainer(newmol, self)
        #if hasattr(newmol, 'spaceGroup'):
        #    self.lazyLoadCommands ('crystalCommands', package='Pmv')
        newmol.unitedRadii = None # set to None to force initial radii assignment
        self.lastMolNumer += 1
        newmol.number = self.lastMolNumer
        
        # create geometries for line rendering
        self.displayLines.initializeMolForCmd(newmol)

        # if bondOrder are present build newmol._bondOrderData and compute vertices
        # for conds that are not dingle bonds
        if newmol._ag._bondOrder is not None: # e.g. when comes from Mol2 file
            newmol.setBondorder(newmol._ag._bondOrder)

        # add molecule to app
        self.addObject(newmol.name, newmol, "Molecule", g)
        newmol.app = weakref.ref(self)

        if group:
            if isinstance(group, str):
                group = self.groups[group]
            assert isinstance(group, Group)
            group.add(newmol)
        else:
            newmol._group = None

        if filename is not None:
            if not newmol.filename:
                newmol.filename = filename
        return newmol

    def getNodesByMolecule(self, nodes, nodeType=None):
        """ moleculeSet, [nodeSet, nodeSet] <- getNodesByMolecule(nodes, nodeType=None)
        nodes can be either: a string, a TreeNode or a TreeNodeSet.
        This method returns a molecule set and for each molecule a TreeNodeSet
        of the nodes belonging to this molecule.
        'nodeType' enables a desired type of nodes to be returned for each
        molecule
        """

        # special case list of complete molecules to be expanded to atoms
        # this also covers the case where nothing is selected
        #if isinstance(nodes, MoleculeSet) or isinstance(nodes, ProteinSet):
        if isinstance(nodes, MoleculeSet):
            if nodeType is Atom:
                atms = []
                for mol in nodes:
                    atms.append(mol.allAtoms)
                return nodes, atms
            elif (nodeType is Molecule):
                return nodes, nodes
        
        # if it is a string, get a bunch of nodes from the string
        if isinstance(nodes, str):
            nodes = self.expandNodes(nodes)

        assert issubclass(nodes.__class__, TreeNode) or \
               issubclass(nodes.__class__, TreeNodeSet)

        # if nodes is a single TreeNode make it a singleton TreeNodeSet
        if issubclass(nodes.__class__, TreeNode):
            nodes = nodes.setClass([nodes])
            nodes.setStringRepr(nodes.full_name())

        if len(nodes)==0: return MoleculeSet([]), []

        # catch the case when nodes is already a MoleculeSet
        if nodes.elementType == Molecule:
            molecules = nodes
        else: # get the set of molecules
            molecules = nodes.top.uniq()

        # build the set of nodes for each molecule
        nodeSets = []

        # find out the type of the nodes we want to return
        searchType=0
        if nodeType is None:
            Klass = nodes.elementType # class of objects in that set
        else:
            assert issubclass(nodeType, TreeNode)
            Klass = nodeType
            if Klass != nodes.elementType:
                searchType=1

        for mol in molecules:
            # get set of nodes for this molecule
            mol_nodes = nodes.get(lambda x, mol=mol: x.top==mol)

            # get the required types of nodes
            if searchType:
                if Klass == Atom and hasattr(mol_nodes, 'allAtoms'):
                    mol_nodes = mol_nodes.allAtoms
                else:
                    mol_nodes = mol_nodes.findType( Klass ).uniq()

            nodeSets.append( mol_nodes )

        return molecules, nodeSets


    def selectionSetFromObject(self, obj):
        """Obj can be a string or a Selection,  or a SelectionSet.
"""
        if isinstance(obj, str):
            # FIXME: we need a reliable way to find out if the string is
            #  a prody selection string or not
            #import pdb; pdb.set_trace()
            if obj.find('chid') > 0 or obj.find('index') > 0: #it is prody
                return self.Mols.select(obj)
            else:
                return self.Mols.selectMolKit(obj)
                #result = selector(obj, self.Mols)
        elif isinstance(obj, MolKitSelection):
            result = SelectionSet([obj])
        elif isinstance(obj, Molecule):
            result = SelectionSet([obj.select()])
        elif isinstance(obj, MoleculeSet):
            result = []
            for mol in obj:
                result.append(mol.select())
            result = SelectionSet(result)
        elif isinstance(obj, SelectionSet):
            result = obj
        else:
            raise ValueError, 'Could not expand nodes %s\n'%str(obj)
        return result


    def getMolFromName(self, name):
        """
        Return the molecule of a given name, or the list of molecules of given
        names.
    
        @type  name: string or list
        @param name: the name of the molecule, or a list of molecule name
        @rtype:   MolKit2.protein or list
        @return:  the molecule or the list of molecules for the given name(s).
        """
        
        if type(name) is list :
            mol = [x for x in self.Mols if x.name in name]
        else : 
            mols = [x for x in self.Mols if x.name==name]
            if len(mols):
                mol = mols[0]
            else:
                mol = None
        return mol
        
    def transformedCoordinatesWithInstances(self, nodes):
        """ for a nodeset, this function returns transformed coordinates.
This function will use the pickedInstance attribute if found.
"""
        # nodes is a list of atoms, residues, chains, etc. where each member
        # has a pickedInstances attribute which is a list of 2-tuples
        # (object, [i,j,..])
        vt = []
        for node in nodes:
            #find all atoms and their coordinates
            coords = nodes.findType(Atom).coords
            if hasattr(node, 'pickedInstances'):
                # loop over the pickedInstances of this node
                for inst in node.pickedInstances:
                    geom, instance = inst # inst is a tuple (object, [i,j,..])
                    M = geom.GetMatrix(geom.LastParentBeforeRoot(), instance[1:])
                    for pt in coords:
                        ptx = M[0][0]*pt[0]+M[0][1]*pt[1]+M[0][2]*pt[2]+M[0][3]
                        pty = M[1][0]*pt[0]+M[1][1]*pt[1]+M[1][2]*pt[2]+M[1][3]
                        ptz = M[2][0]*pt[0]+M[2][1]*pt[1]+M[2][2]*pt[2]+M[2][3]
                        vt.append( (ptx, pty, ptz) )
            else:
                # no picking ==> no list of instances ==> use [0,0,0,...] 
                g = nodes[0].top.geomContainer.geoms['master']
                M = g.GetMatrix(g.LastParentBeforeRoot())
                for pt in coords:
                    ptx = M[0][0]*pt[0]+M[0][1]*pt[1]+M[0][2]*pt[2]+M[0][3]
                    pty = M[1][0]*pt[0]+M[1][1]*pt[1]+M[1][2]*pt[2]+M[1][3]
                    ptz = M[2][0]*pt[0]+M[2][1]*pt[1]+M[2][2]*pt[2]+M[2][3]
                    vt.append( (ptx, pty, ptz) )
                
        return vt
    
    def getSelectionString(self, sel):
        ag = sel.getAtomGroup()
        mol = sel.getAtomGroup().getMolecule()
        if len(sel) == len(ag):
            return mol.name
        hv = mol._ag.getHierView()
        selStr = ""
        for chain in hv.iterChains():
            chid = chain.getChid()
            inChain = sel.select('chid %s'%chid)
            if inChain is not None:
                if len(inChain)==len(mol.select('chid %s'%chid)):
                    selStr += ' or chid %s'%chid
                else:
                    selStr += ' or index %s'%' '.join([str(x) for x in inChain.getIndices()])

        selStr = selStr[4:] # remove first or
        return "%s:%s" %(mol.name, selStr)

    def addNamedSelection(self, selSet, name, group=None):
        assert isinstance(selSet, SelectionSet)
        self.namedSelections[name] = selSet
        if self.gui:
            if group:
                if isinstance(group, str):
                    parent = self.groups.get(group, None)
                else:
                    assert isinstance(group, Group)
                    parent = group._treeItems.keys()[0]
                group.add(selSet)
            else:
                parent = None
            self.gui().objTree.addObject(
                selSet, parent, name,
                color=self.gui().objTree.colors['namedSelections'])
        
    def saveFullSession(self, name, inlcudeMultiMol):
        """
        save a session file

        None <- mv.saveFullSession(self, name, inlcudeMultiMol)

        create name.psf file. The extension is only added if it is not
        already in name

        includeMultiMol is True of False. When set to True molecules loaded
        from files containing multiple molecules will save the entire source
        file and after restoring the session it will be possible to browse
        the content of the file by switching to other molecules in the file.
        When set to False, only the currently shown molecule will be saved
        and restored.
        """
        ##
        ## we create a file called name.psf. if name is /a/b/c.psf we call
        ## 'c' the basename. The name.psf file is a tar giz'ed directory
        ## containing the directory c.psf_dir in which we store molecules
        ## and a python script called session.py
        ##
        bname = os.path.basename(name)+'_dir'  # i.e. mysession.psf

        # create a temporary folder
        folder = tempfile.mktemp()
        sessionFolder = os.path.join(folder, bname)
        #print 'creating folder', sessionFolder
        os.mkdir(folder)
        os.mkdir(sessionFolder)

        # create tar object in the location required by user
        #print 'mktar', name
        # the compression takes too long and can not be interrupted safely
        # write the session uncompressed
        #tar = tarfile.open(name, "w:gz")
        tar = tarfile.open(name, "w")

        # goto folder 
        cwd = os.getcwd()
        #print 'goto folder', sessionFolder
        os.chdir(sessionFolder)

        ## save molecules
        for mol in self.Mols:
            # write prody atom groups to npz files in that folder
            saveAtoms(mol._ag, "%s.%s"%(mol._basename, mol.curMolIndex()))
            if mol._multi=='molecules' and inlcudeMultiMol:
                os.chdir(cwd)
                shutil.copy(mol.filename,os.path.join(\
                    sessionFolder, os.path.split(mol.filename)[1]+'._multiMol'))
                os.chdir(sessionFolder)

            ## save representation sets 
            mol.geomContainer.atoms.keys()
            setDict = {}
            for k,v in mol.geomContainer.atoms.items():
                setDict[k] = v._indices
            ostream = open('%s.sets.npz'%mol._basename, 'wb')
            savez(ostream, **setDict)
            ostream.close()

            ## save MSMS
            if hasattr(mol, '_msmsData'):
                msDict = {}
                for name in mol._msmsData['msms'].keys():
                    msDict[name] = {'atoms': mol._msmsData['atoms'][name]._indices,
                                    'params': mol._msmsData['params'][name]}
                ostream = open('%s.msms.npz'%mol._basename, 'wb')
                savez(ostream, **msDict)
                ostream.close()

            ## save Cartoon
            #if hasattr(mol, '_cartoonGeoms'):
            #    msDict = {}
            #    for name in mol._msmsData['msms'].keys():
            #        msDict[name] = {'atoms': mol._msmsData['atoms'][name]._indices,
            #                        'params': mol._msmsData['params'][name]}
            #    ostream = open('%s.msms.npz'%mol._basename, 'wb')
            #    savez(ostream, **msDict)
            #    ostream.close()

            propDict = {}

            # save colors and mol._renderingProp['lines']
            propDict['colors'] = mol._colors
            propDict['renderingProp'] = mol._renderingProp

            #import pdb; pdb.set_trace()
            ostream = open('%s.prop.npz'%mol._basename, 'wb')
            savez(ostream, **propDict)
            ostream.close()

        ## save named selections
        with open('selections.npz', 'wb') as f:
            for sel in self.curSelection:
                f.write("CURRENT %s\n"%sel._molecule()._basename)
                numpy.save(f, sel._indices)
            
            for name, selSet in self.namedSelections.items():
                for sel in selSet:
                    f.write("%s %s\n"%(name, sel._molecule()._basename))
                    numpy.save(f, sel._indices)
        
        vstate = self.gui().viewer.getViewerStateDefinitionCode('self.gui().viewer')
        ostate = self.gui().viewer.getObjectsStateDefinitionCode('self.gui().viewer')
        f = open('session.py', 'w')
        f.write("mode='both'\n")
        [f.write(line) for line in vstate]
        [f.write(line) for line in ostate]
        f.close()
        
        # go to parent folder
        os.chdir(folder)
        # add subdirectory with basename
        tar.add(bname)
        tar.close()
        os.chdir(cwd)

        #delete the temporary folder
        shutil.rmtree(folder)
        print "Done saving session %s"%name
        
    def readFullSession(self, name):
        """
        read a session file

        None <- mv.readFullSession(self, name)

        name has to be the full path to a .tar.gz file as created by
        mv.saveSession 
        """
        # suspend redraw to save time and avoid flashing
        if self.gui is not None:
            self.gui().viewer.suspendRedraw = True
        #print 'restoring session', name

        bname = os.path.basename(name)+'_dir'  # i.e. mysession.psf
        import tempfile, tarfile, shutil
        folder = tempfile.mktemp()

        # open tar file and extract into folder called name 
        tar = tarfile.open(name)
        #print 'extract to ',folder
        tar.extractall(path=folder)

        # save current directory
        cwd = os.getcwd()
        sessionFolder = os.path.join(folder, os.listdir(folder)[0])
        #print 'going to', sessionFolder
        # go to session folder
        os.chdir(sessionFolder)
        closeSession = True
        #try:
        if True:
            # load molecules
            from MolKit2 import Molecule, Read
            molnames = glob("*.ag.npz")
            multiMolNames = glob("*._multiMol")
            mmbase = []

            mmformat = []
            mmnames = []
            for name in multiMolNames:
                newname = os.path.splitext(name)[0] # remove ._multiMol
                base, ext = os.path.splitext(newname)
                mmbase.append(base)
                mmformat.append(ext)
                mmnames.append(name)
            
            for molname in glob("*.ag.npz"):
                # remove .confNum.ag.npz
                mname, confIndex, et1, ext2 = os.path.basename(molname).split('.')
                confIndex = int(confIndex)
                nag = loadAtoms(molname)
                if mname in mmbase:
                    ind = mmbase.index(mname)
                    newname = "%s%s"%(mmbase[ind],mmformat[ind])
                    os.rename(mmnames[ind], newname)
                    # read the file with an absolute path so that we can
                    # later save it again
                    mol = Read(os.path.abspath(newname))
                    mol.gotoMolecule(confIndex)
                    mol.setAtomgroup(nag)
                    # keep this session open so we can use multimolecules
                    closeSession = False
                    abspath = os.path.abspath(sessionFolder)
                    if abspath not in self._openedSessions:
                        self._openedSessions.append(abspath)
                else:
                    mol = Molecule(mname, nag)
                    if confIndex>0:
                       mol. _ag.setACSIndex(confIndex)

                self.addMolecule(mol)
                event = AfterAddMoleculeEvent(molecule=mol)
                self.eventHandler.dispatchEvent(event)

                # load MSMS
                from MolKit2.selection import Selection
                if os.path.exists('%s.msms.npz'%mname):
                    setDict = load('%s.msms.npz'%mname)
                    surfD = {}
                    for surfName in setDict.files:
                        d = setDict[surfName].reshape(-1)[0]
                        surfAtoms = Selection(nag, d['atoms'], "")
                        p = d['params']
                        hdset = p.get('hdset', None)
                        if hdset is not None:
                            hdset = Selection(nag, d['hdset'], "")
                        radii = p.get('radii', None)
                        self.computeMSMS(surfAtoms, surfName=surfName,
                                         probeRadius=p['probeRadius'], density=p['density'],
                                         radii=radii, hdset=hdset, hdensity=p['hdensity'])
                
                # load atom sets
                setDict = load('%s.sets.npz'%mname)
                hasCartoon = False
                for k in setDict.files:
                    v = setDict[k]
                    mol.geomContainer.setAtomsForGeom(k, Selection(nag, v, ""))
                    if k.startswith('chain_'):
                        hasCartoon = True

                for cmd in [self.displayLines, self.displayCPK, self.displaySB,
                            self.labelAtoms, self.labelResidues]:
                # displayLines initialized mol in addMolecule
                #for cmd in [self.displayCPK, self.displaySB]:
                    atomSet = mol.geomContainer.atoms.get(cmd._atomSetName, None)
                    if atomSet:
                        # create the data structures
                        cmd.initializeMolForCmd(mol)

                # load molecule properties
                propDict = load('%s.prop.npz'%mname)
                for k in propDict.files:
                    v = propDict[k]
                    if k == 'colors':
                        mol._colors = v.reshape(-1)[0]
                    elif k=='renderingProp':
                        mol._renderingProp = v.reshape(-1)[0]
                    elif k=='cpkProperties':
                        mol._cpkProperties = v.reshape(-1)[0]

                if self.gui is not None:
                    event = RefreshDisplayEvent(molecule=mol)
                    self.eventHandler.dispatchEvent(event)
                
                if hasCartoon:
                    self.computeCartoon(mol)
                    if self.displayCartoon.isLoader():
                        cmd = self.displayCartoon.loadCommand()
                    self.displayCartoon.initializeMolForCmd(mol)
                    self.displayCartoon.refreshDisplay(mol)
                    
            # end loop over molecules
            
            # load selections
            with open('selections.npz') as f:
                cursels = []
                selections = {}
                while True:
                    key = f.readline()
                    if not key: break
                    indices = numpy.load(f)
                    selname, molname = key.split()
                    mol = self.getMolFromName(molname)
                    sel = Selection(mol._ag, indices, "")
                    if selname=='CURRENT':
                        #print 'current', mol.name, indices
                        sel = Selection(mol._ag, indices, "")
                        cursels.append(sel)
                    else:
                        #print selname, mol.name, indices
                        if selections.has_key(selname):
                            selections[selname].append(sel)
                        else:
                            selections[selname] = [sel]

            lastSele = None
            for name, sels in selections.items():
                #print 'creating', name, sels
                if lastSele is not None:
                    self.activeSelection = self.curSelection
                    self.gui().objTree.setCurrentSelection(None)
                self.select(SelectionSet(sels, name=name))
                self.rename(self.curSelection, '%s'%name)
                # TRY USING THIS INSTEAD 
                #self.addNamedSelection(SelectionSet(sels, name=name), name)
                #self.clearSelection()
            #print 'creating cureSel', cursels
            self.activeSelection = self.curSelection
            self.gui().objTree.setCurrentSelection(None)
            self.select(SelectionSet(cursels))

            # execute session.py file inside untar'ed folder
            execfile('./session.py', {'self':self})
            os.chdir(cwd)

            # delete the untar'ed folder
            #print 'removing ', name+'_dir'
            if hasattr(self,'recentFiles'):
                self.recentFiles.add(os.path.abspath(str(name)), "readFullSession")

        ## except:
        ##     print "trapExceptions:", self.trapExceptions
        ##     import sys
        ##     if self.trapExceptions is False:
        ##         exc_info = sys.exc_info()
        ##         raise exc_info[1], None, exc_info[2]
        ##     else:
        ##         #import pdb; pdb.set_trace()
        ##         msg = 'Error in reading session file  %s'%(name)
        ##         self.errorMsg(sys.exc_info(), msg, obj=None)
        ##         self._executionReport.finalize()
        ## finally:
        if closeSession:
            shutil.rmtree(folder)
        # allow viewer to redraw
        if self.gui is not None:
            self.gui().viewer.suspendRedraw = False

    def saveFullSessionOLD(self, name):
        """
        save a session file

        None <- mv.saveFullSession(self, name)

        create name.psf file. The extension is only added if it is not
        already in name
        """
        ##
        ## we create a file called name.psf. if name is /a/b/c.psf we call
        ## 'c' the basename. The name.psf file is a tar giz'ed directory
        ## containing the directory c.psf_dir in which we store molecules
        ## and a python script called session.py
        ##

        bname = os.path.basename(name)+'_dir'  # i.e. mysession.psf

        # create a temporary folder
        import tempfile, tarfile, shutil
        folder = tempfile.mktemp()
        sessionFolder = os.path.join(folder, bname)
        #print 'creating folder', sessionFolder
        os.mkdir(folder)
        os.mkdir(sessionFolder)

        # create tar object in the location required by user
        #print 'mktar', name
        tar = tarfile.open(name, "w:gz")

        # goto folder 
        cwd = os.getcwd()
        #print 'goto folder', sessionFolder
        os.chdir(sessionFolder)

        #try:
        #write all molecules into this file
        # FIXME .. we should use the same writer as the parser i.e. PQR etc
        for mol in self.Mols:
            #self.writeMoleculeToSessionFolder(mol)
            filename = mol.name + ".pdb"
            records = mol.getPDBRecords()
            f = open(filename, "w")
            for rr in records:
                f.write(rr)
            f.close()
        # write the session file in the folder
        threshold = numpy.get_printoptions()["threshold"]
        lines = """"""
        for mol in self.Mols:
            lines += self.getStateCodeForMolecule(mol)

        #lines += self.getStateCodeForMeasureGeoms()

        #lines += self.getStateCodeForSets()
        lines += self.getStateCodeForSelection()

        #if self.commands.has_key('sequenceViewer'):
        #    lines += self.sequenceViewer.getStateCodeForSeqViewer()

        vstate = self.gui().viewer.getViewerStateDefinitionCode('self.gui().viewer')
        ostate = self.gui().viewer.getObjectsStateDefinitionCode('self.gui().viewer')
        
        numpy.set_printoptions(threshold=threshold)
        # write vision networks
        #ed = self.vision.ed
        #for name,net in ed.networks.items():
        #    if len(net.nodes):
        #        net.saveToFile(os.path.join(sessionFolder,
        #                                    net.name+'_pmvnet.py'))

        f = open('session.py', 'w')
        [f.write(line) for line in lines]
        f.write("mode='both'\n")
        [f.write(line) for line in vstate]
        [f.write(line) for line in ostate]

        # add lines to session file to load vision networks
        #f.write('ed = self.vision.ed\n')
        #for name,net in ed.networks.items():
        #    if len(net.nodes):
        #        f.write("ed.loadNetwork('%s')\n"%(net.name+'_pmvnet.py'))
        #if animcode is not None:
        #    [f.write(line) for line in animcode]

        f.close()
        # go to parent folder
        os.chdir(folder)

        # add subdirectory with basename
        tar.add(bname)
        tar.close()

        #finally:
        # restore original directory
        os.chdir(cwd)
        #delete the temporary folder
        print 'removing ', folder
        shutil.rmtree(folder)

    def getStateCodeForMolecule(self, mol):
        ##
        ## read in all molecules
        ##
        lines = """##\n## read in all molecules\n##\n"""
        #import pdb; pdb.set_trace()
        ## generate Pmv commands to load the molecule
        molName = mol.name
        ## check for multimodel file ("modelsAs") :
        ## if hasattr(mol.parser, 'modelsAs') and mol.parser.modelsAs != "molecules":
        ##     # add "modelsAs" keyword to self.readMolecule:
        ##     lines += """mols = self.Mols.get('%s')\nif len(mols)==0:\n\tmols = self.readMolecule('%s', addToRecent=False, modelsAs='%s')\n"""%(molName, molFile, mol.parser.modelsAs)
        ## else:
        lines += """mols = self.Mols.getMolecules('%s')\nif len(mols)==0:\n\tmols = self.readMolecule('%s', addToRecent=False)\n"""%(molName, molName+".pdb")
        lines += """mol = mols[0]\n"""
        
        if mol._ag._bondOrder is not None:
            lines += """if mol._ag._bondOrder is None:\n"""
            lines += """    mol.typeAtomsAndBonds()\n"""
        #lines += "currentAtoms = mols.getAtomGroup()\n"
        currentAtoms = mol.getAtomGroup()

        def numarr2str(val):
            import numpy
            threshold = numpy.get_printoptions()["threshold"]
            if val.size > threshold:
                numpy.set_printoptions(threshold=val.size)
            if len(val.shape) == 1:
                valstr = "array(%s ,'%s')" %(numpy.array2string(val, precision=4, separator =",") , val.dtype.char)
            else:
                n = val.shape[1]
                if val.dtype in [numpy.float, numpy.float32, numpy.float64, numpy.double]:
                    formatStr = "[" + "{0[%s]:.3f}, " *n % tuple(range(n))  + "],"
                else:
                    formatStr = "[" + "{0[%s]}, " *n % tuple(range(n))  + "],"
                valstr = "array(["
                for v in val:
                    valstr += formatStr.format(v)
                valstr += "], '%s')" % val.dtype.char
            #numpy.set_printoptions(threshold=threshold)
            return valstr

        for key in mol._ag._data.keys():
            if not key in ['resnum', 'name', 'chain', 'altloc', 'colors', 'occupancy',
                           'element', 'beta', 'icode', 'radius', 'resindex', 'resname',
                           'chindex', 'serial', 'segment', 'atomicNumber', 'segindex', 'numbonds']:
                datastr = "numpy.%s\n"%(numarr2str(mol._ag._data[key]))
                if key in ['hetatm', 'pdbter', 'deleted']:
                    lines += "mol._ag.Flags('%s', %s)\n" % (key, datastr)
                else:
                    lines += "mol._ag.setData('%s', %s)\n" % (key, datastr)

        ##
        ## create geometry
        ##
        lines += """self.displayLines(mol, negate=True)\n"""
        lines += """##\n## create geometries\n##\n"""
        ## generate Pmv commands to create geometry that are currently displayed
        gc = mol.geomContainer
        gca = gc.atoms
        gcg = gc.geoms
        ##
        ## Lines
        ##
        if gca.has_key("singleBonds") and len(gca["singleBonds"]):
            colors = self.color.getAtomColors(gc.geoms['singleBonds'])
            lines += "col1 = numpy.%s\n"%(numarr2str(colors))
            lw = gc.geoms["singleBonds"].lineWidth
            lines += """self.displayLines('%s', displayBO=True, lineWidth=%d)\n"""%(self.getSelectionString(gca["singleBonds"]), lw)
            lines += """mol.geomContainer.geoms['singleBonds'].Set(materials=col1, redo=1, inheritMaterial=False, tagModified=False, transparent='implicit')\n"""
            if mol._doubleBonds is not None and len(mol._doubleBonds):
                bonds = mol.select().getBonds()
                colDB = self.color.getAtomColors(gc.geoms["doubleBonds"])
                for i1,i2 in bonds[2]:
                    k = "%d %d" %(i1, i2)
                    if mol._doubleBonds.has_key(k):
                        n = mol._doubleBonds[k] * 4
                        colDB[n] = colors[i1]
                        colDB[n+1]= colors[i2]
                        colDB[n+2] = colors[i1]
                        colDB[n+3] = colors[i2]
                lines += "col2 = numpy.%s\n"%(numarr2str(colDB))
                lines += """mol.geomContainer.geoms['doubleBonds'].Set(materials=col2, redo=1, inheritMaterial=False, tagModified=False, transparent='implicit')\n"""
            if mol._tripleBonds is not None and len(mol._tripleBonds):
                colTB = self.color.getAtomColors(gc.geoms["tripleBonds"])
                # loop over bonds
                for i1,i2 in bonds[3]:
                    k = "%d %d" %(i1, i2)
                    if mol._tripleBonds.has_key(k):
                        n = mol._tripleBonds[k] * 6
                        colTB[[n, n+1, n+2, n+3, n+4, n+5]] = \
                                  [colors[i1], colors[i2],
                                   colors[i1], colors[i2],
                                   colors[i1], colors[i2]]
                lines += "col3 = numpy.%s\n"%(numarr2str(colTB))
                lines += """mol.geomContainer.geoms['tripleBonds'].Set(materials=col3, redo=1, inheritMaterial=False, tagModified=False, transparent='implicit')\n"""
            if mol._aromaticArcs is not None and len(mol._aromaticArcs):
                colAB = self.color.getAtomColors(gc.geoms['aromaticBonds'])
                for index in mol._aromaticArcs.keys():
                    for ii in mol._aromaticArcs[index]:
                        colAB[ii[0]] = colors[index]
                lines += "col4 = numpy.%s\n"%(numarr2str(colAB))
                lines += """mol.geomContainer.geoms['aromaticBonds'].Set(materials=col4, redo=1, inheritMaterial=False, tagModified=False, transparent='implicit')\n"""
        ##
        ## molecular surfaces
        ##
        if hasattr(gc, "msms"):
            for name, srfa in gc.msms.items():
                srf = srfa[0]
                atoms = gca[name]
                if len(gca[name]): # This surface is displayed
                    lines += """self.computeMSMS('%s', hdensity=%f, hdset=%s, density=%f, pRadius=%f, perMol=%d, noHetatm=%d, display=False, surfName='%s')\n"""%(
                        self.getSelectionString(atoms), srf.hdensity, srf.hdset, srf.density,
                        srf.probeRadius, srf.perMol, srf.noHetatm, name)

                else:
                    lines += """self.computeMSMS('%s', hdensity=%f, hdset=%s, density=%f, pRadius=%f, perMol=%d, noHetatm=%d, display=False, surfName='%s')\n"""%(
                        mol.name, srf.hdensity, srf.hdset, srf.density,
                        srf.probeRadius, srf.perMol, srf.noHetatm, name)
                col = self.color.getAtomColors(gc.geoms[name])
                lines += "col = numpy.%s\n"%(numarr2str(col))
                lines += "mol.geomContainer.geoms['%s'].Set(materials=col, redo=1, inheritMaterial=False, tagModified=False, transparent='implicit')\n"%(name,)

        if hasattr(gc, "boundGeom"):
            for name in gc.boundGeom.keys():
                atoms = gc.boundGeom[name]['atoms']
                opts = mol._coarseMolSurfParams[name]
                if len(atoms) == len(mol._ag):
                    lines += """self.computeCoarseMolecularSurface(mol, surfName='%s', gridSize=%d, padding=%f, resolution=%f, bind_surface_to_molecule=True, isovalue=%f)\n"""%(name, opts['gridSize'], opts['padding'], opts['resolution'], opts['isovalue'])
                else:
                    lines += """self.computeCoarseMolecularSurface('%s', surfName='%s', gridSize=%d, padding=%f, resolution=%f, bind_surface_to_molecule=True, isovalue=%f)\n"""%(
                        self.getSelectionString(atoms), name, opts['gridSize'], opts['padding'], opts['resolution'], opts['isovalue'])
                col = self.color.getAtomColors(gc.geoms[name])
                lines += "col = numpy.%s\n"%(numarr2str(col))
                lines += "mol.geomContainer.geoms['%s'].Set(materials=col, redo=1, inheritMaterial=False, tagModified=False, transparent='implicit')\n"%(name,)
                if len(gca[name]) == 0:
                    lines += """self.displayBoundGeom(mol, '%s', negate=True)\n""" % name
                elif len(gca[name]) ==  len(mol._ag):
                    lines += """self.displayBoundGeom(mol, '%s', negate=False)\n"""% name
                else:
                    lines += """self.displayBoundGeom('%s', '%s')\n"""% (self.getSelectionString(gca[name]), name)
                    
        ## Ribbons (not implemented)
        ##

        ##
        ## beaded Ribbons (not implemented)
        ##
        ## vi = self.GUI.VIEWER
        ## master = vi.FindObjectByName('root|%s|beadedRibbon'%(mol.name,))
        ## if master:
        ##     # code to restore mol.strandVarValue
        ##     d = {}
        ##     for k,v in mol.strandVar.items():
        ##         d[k] = v.get()
        ##     lines += "mols[0].strandVarValue = %s\n"%str(d)
        ##     cmds.append( (setattr, (mol, "strandVarValue", d), {}) )
        ##     lines += "self.beadedRibbons('%s', **%s)\n"%(
        ##         mol.name, str(mol.beadedRibbonParams))
        ##     cmds.append((self.beadedRibbons, (mol.name,),  mol.beadedRibbonParams))
            
        ## display geometry
        ##
        ##
        
        ##
        ## Sticks and Balls
        ##
        data={}
        if gca.has_key('sticks') and gca['sticks']:
            # sticks radius: cRad
            if numpy.allclose(min(mol._cRad), max(mol._cRad)):
                lines += "mol._cRad = numpy.ones( (%i,), 'f') * %f\n"%(len(mol._ag), mol._cRad[0])
            else:
                lines += "mol._cRad = numpy.%s\n"%(numarr2str(mol._cRad))
            # ball radius: ballRad
            if numpy.allclose(min(mol._ballRad), max(mol._ballRad)):
                lines += "mol._ballRad = numpy.ones( (%i,), 'f') * %f\n"%(len(mol._ag), mol._ballRad[0])
            else:
                lines += "mol._ballRad = numpy.%s\n"%(numarr2str(mol._ballRad))
            if numpy.allclose(min(mol._ballScale), max(mol._ballScale)):
                lines += "mol._ballScale = numpy.ones( (%i,), 'f') * %f\n"%(len(mol._ag), mol._ballScale[0])
            else:
                lines += "mol._ballScale = numpy.%s\n"%(numarr2str(mol._ballScale))

            sticksQuality = gc.geoms['sticks'].quality
            ballsQuality = gc.geoms['balls'].quality
            lines += """self.displaySticksAndBalls('%s', cquality=%d, sticksBallsLicorice='Sticks and Balls', bquality=%d, setScale=False)\n"""%(self.getSelectionString(gca['sticks']), sticksQuality, ballsQuality, )
            col1 = self.color.getAtomColors(gc.geoms['sticks'])
            col2 = self.color.getAtomColors(gc.geoms['balls'])
            lines += "stcol = numpy.%s\n"%(numarr2str(col1))
            lines += """mol.geomContainer.geoms['sticks'].Set(materials=stcol, redo=1, inheritMaterial=False, tagModified=False, transparent='implicit')\n"""
            lines += "bcol = numpy.%s\n"%(numarr2str(col2))
            lines += """mol.geomContainer.geoms['balls'].Set(materials=bcol, redo=1, inheritMaterial=False, tagModified=False, transparent='implicit')\n"""
        ##
        ## CPK and Balls
        ##
        if gca.has_key('cpk') and gca['cpk']:
            # save atm.cpkScale unless == 1.0
            # and atm.cpkRad unless == 0.0
            if numpy.allclose(min(mol._cpkRad), max(mol._cpkRad)):
                if mol._cpkRad[0] != 0.0:
                    lines += "mol._cpkRad = numpy.ones( (%i,), 'f') * %f\n"%(len(mol._ag), mol._cpkRad[0])
            else:
                lines += "mol._cpkRad = numpy.%s\n"%(numarr2str(mol._cpkRad))

            if numpy.allclose(min(mol._cpkScale), max(mol._cpkScale)):
                if mol._cpkScale[0] != 1.0:
                    lines += "mol._cpkScale = numpy.ones( (%i,), 'f') * %f\n"%(len(mol._ag), mol._cpkScale[0])
            else:
                lines += "mol._cpkScale = numpy.%s\n"%(numarr2str(mol._cpkScale))

            # display CPK with setScale=False
            quality = gc.geoms['cpk'].quality
            atms = gca['cpk']
            reprstring = self.getSelectionString(atms)
            lines += """self.displayCPK('%s', quality=%d, setScale=False)\n"""%(reprstring, quality)
            col = self.color.getAtomColors(gc.geoms['cpk'])
            lines += "col = numpy.%s\n"%(numarr2str(col))
            lines += """mol.geomContainer.geoms['cpk'].Set(materials=col, redo=1, inheritMaterial=False, tagModified=False, transparent='implicit')\n"""

        ##
        ## Ribbons (not implemented)
        ##
        ## ribbonResidues = ResidueSet([], stringRepr='abc')
        ## for k, v in gca.items():
        ##     if k[:4] in ['Heli', 'Stra', 'Turn', 'Coil'] and gca[k]:
        ##         ribbonResidues += gca[k]
        ## if ribbonResidues:
        ##     lines += """self.displayExtrudedSS('%s')\n"""%ribbonResidues.buildRepr()

        ##
        ## Labels
        ##
        for gname, cmd in [('atomLabels','labelAtoms'), ('residueLabels', 'labelResidues')]:
            if gca.has_key(gname) and gca[gname]:
                labelGeom = gc.geoms[gname]
                state = labelGeom.getState()
                labTran =  state.get("labelTranslation", None)
                if labTran is not None:
                    if len(numpy.unique(labTran) == 1):
                        state.pop("labelTranslation")
                verts = labelGeom.vertexSet.vertices.array
                lines += "from numpy import *\n"
                lines += "verts = numpy.%s\n" % numarr2str(verts)
                lines += "state = %s\n"%str(str(state))
                lines += "state['vertices']=verts\n"
                col = self.color.getAtomColors(labelGeom)
                lines += "col = numpy.%s\n"%(numarr2str(col))
                selStr = self.getSelectionString(gca[gname])
                if gname == 'residuelabels':
                    mode = labelGeom.params['mode']
                    lines += "self.%s('%s', mode ='%s')\n" %(cmd, selStr, mode)
                else:
                    lines += "self.%s('%s')\n" %(cmd, selStr)
                lines += "labelGeom = self.gui().viewer.FindObjectByName('%s')\n"%(labelGeom.fullName)
                lines += "if labelGeom:\n"
                lines += "    labelGeom.Set(**state)\n\n"
                lines += "    mol.geomContainer.geoms['%s'].Set(materials=col, redo=1, inheritMaterial=False, tagModified=False, transparent='implicit')\n"%(gname,)
        if hasattr(gc, "msms"):   
            for name in gc.msms.keys():
                if gca[name]:
                    lines += """self.displayMSMS('%s', negate=False, only=False, surfName='%s', nbVert=1)\n"""%(self.getSelectionString(gca[name]), name)

        ##                
        ## check if molecule is visible:
        ##    
        if not gc.geoms['master'].visible:
            lines += "self.showMolecules(['%s'], show=0)\n" % (mol.name)
        return lines

    def getStateCodeForSelection(self):
        ##
        ## restore selection
        ##
        lines = """##n\## restore selection\n##\n"""
        lastSele = None
        if len(self.namedSelections):
            for sname, sele in self.namedSelections.items():
                if lastSele is not None:
                    lines += "self.activeSelection = self.curSelection\n"
                    lines += "self.gui().objTree.setCurrentSelection(None)\n"
                for item in sele:
                    lines += "self.select('%s')\n"%self.getSelectionString(item)
                lines += "self.rename(self.curSelection, '%s')\n"%sname
                lastSele = sele

        if len(self.curSelection):
            if lastSele is not None:
                lines += "self.activeSelection = self.curSelection\n"
                lines += "self.gui().objTree.setCurrentSelection(None)\n"
            for item in self.curSelection:
                lines += "self.select('%s')\n"%self.getSelectionString(item)

        if self.gui().activeSelection is not None:
            if self.activeSelection != self.curSelection:
                lines += "self.activeSelection = self.curSelection\n"
                lines += "self.gui().objTree.setCurrentSelection(None)\n"
                for sname, sele in self.namedSelections.items():
                   if sele == self.activeSelection:
                       lines += "self.setActiveSelection(self.namedSelections['%s'])\n"% sname
                       break
        else:
            lines += "self.activeSelection = self.curSelection\n"
            lines += "self.gui().objTree.setCurrentSelection(None)\n"
            #lines += "self.activeSelection = []\n"
        return lines

    def readFullSessionOLD(self, name):
        """
        read a session file

        None <- mv.readFullSession(self, name)

        name has to be the full path to a .tar.gz file as created by
        mv.saveSession 
        """
        # suspend redraw to save time and avoid flashing
        if self.gui is not None:
            self.gui().viewer.suspendRedraw = True
        #print 'restoring session', name

        bname = os.path.basename(name)+'_dir'  # i.e. mysession.psf
        import tempfile, tarfile, shutil
        folder = tempfile.mktemp()
        #sessionFolder = os.path.join(folder, bname)
        #sessionFolder = os.path.join(folder, "*.psf_dir")
        #print 'creating folder', folder
        #os.mkdir(folder)

        # open tar file and extract into folder called name 
        tar = tarfile.open(name)
        #print 'extract to ',folder
        tar.extractall(path=folder)

        # save current directory
        cwd = os.getcwd()
        sessionFolder = os.path.join(folder, os.listdir(folder)[0])
        #print 'going to', sessionFolder
        # go to session folder
        os.chdir(sessionFolder)

        try:
            # execute session.py file inside untar'ed folder
            execfile('./session.py')
            os.chdir(cwd)
            
            # delete the untar'ed folder
            #print 'removing ', name+'_dir'
            self.openedSessions.append(folder)
            shutil.rmtree(folder)
            if hasattr(self,'recentFiles'):
                self.recentFiles.add(os.path.abspath(str(name)), "readFullSession")
        except:
            shutil.rmtree(folder)
            #print "trapExceptions:", self.trapExceptions
            import sys
            if self.trapExceptions is False:
                exc_info = sys.exc_info()
                raise exc_info[1], None, exc_info[2]
            else:
                #import pdb; pdb.set_trace()
                msg = 'Error in reading session file  %s'%(name)
                self.errorMsg(sys.exc_info(), msg, obj=None)
                self._executionReport.finalize()
        finally:
            # allow viewer to redraw
            if self.gui is not None:
                self.gui().viewer.suspendRedraw = False


    def bindGeomToMolecularFragment(self, geom, sel, cutoff_from=3.5,
                                    cutoff_to=40.0, instanceMatrices=None):
        """Bind an IndexedGeom geometry (for example a coarse molecular surface) to
        a molecule."""
        mol = sel.getAtomGroup().getMolecule()
        cl_atoms = []
        geomC = mol.geomContainer
        reparent = True
        # check if geom is already under the molecule's hierarchy of
        # DejaVu2 Objects
        obj = geom
        while obj.parent:
            if geom.parent==geomC.masterGeom:
                reparent=False
                break
            obj = obj.parent
        try:
            gname = geom.name
            vs = geom.vertexSet.vertices.array
            from DejaVu2.IndexedGeom import IndexedGeom
            if isinstance(geom, IndexedGeom):
                fs = geom.faceSet.faces.array
                fns = geom.getFNormals()
            else:
                fs = fns = None
            #print "in doit: lenvs:", len(vs)
            from mglutil.bhtfunctions import findClosestAtoms
            atom_coords = sel.getCoords().astype("f")
            cl_atoms = findClosestAtoms(vs, atom_coords, cutoff_from, cutoff_to, instanceMatrices)
            if len(cl_atoms):
                if geom.parent==None:   #add geometry to Viewer:
                    geomC.addGeom(geom, parent=geomC.masterGeom, redo=0)
                    reparent = False
                if reparent: # need to reparent geometry
                    vi= geomC.VIEWER
                    if vi:
                        vi.ReparentObject(geom, geomC.masterGeom)
                    else:
                        oldparent = geom.parent
                        oldparent.children.remove(geom)
                        geomC.masterGeom.children.append(geom)
                        geom.parent = geomC.masterGeom
                    geom.fullName = geom.parent.fullName+'|'+gname
                geomC.geomPickToBonds[gname] = None
                geom.userName = gname
                geomC.geoms[gname] = geom
                import weakref
                geom.mol = weakref.ref(mol)
                inds = numpy.unique(cl_atoms)
                #u_atoms = sel.select("index " +  "%s "*len(inds) %tuple(inds))
                #geomC.setAtomsForGeom(gname, u_atoms)
                geomC.setAtomsForGeom(gname, mol.emptySelection())
                #self.cl_atoms[gname] = cl_atoms
                #create a unique identifier for all atoms
                #ids = [x.getName()+str(i) for i, x in enumerate(sel)]
                ids = sel.getIndices()
                mol_lookup = dict(zip(ids, range(len(ids))))
                
                # highlight selection (done in PmvGUI.py)
                ## lAtomVerticesDict = {}
                ## for lIndex in range(len(cl_atoms)):
                ##     if not lAtomVerticesDict.has_key(cl_atoms[lIndex]):
                ##         lAtomVerticesDict[cl_atoms[lIndex]]=[]
                ##     lAtomVerticesDict[cl_atoms[lIndex]].append(lIndex)
                if not hasattr (geomC, "boundGeom"):
                    geomC.boundGeom = {}
                geomC.boundGeom[geom.name]={'cl_atoms':cl_atoms, 'fs':fs,
                                                'fns':fns, 'atoms':sel,
                                                'mol_lookup':mol_lookup, 
                                       #'atomVertices':lAtomVerticesDict
                                            }
                geom._boundGeomType = True
            else:
                raise RuntimeError("%s: no close atoms found for geometry: %s" % (self.name, geom.name))
        except:
            msg = 'Error in binding geometry %s to molecule %s'%(geom.name, mol.name)
            import sys
            if self.trapExceptions is False:
                exc_info = sys.exc_info()
                raise exc_info[1], None, exc_info[2]
            else:
                self.errorMsg(sys.exc_info(), msg, obj=atoms)                       

        return cl_atoms

    def saveSelectionAsMolecule(self, obj, filename):
        name, ext = os.path.splitext(os.path.basename(filename))
        if isinstance(obj, Molecule):
            mol = obj
        else:
            mol = obj.toMolecule(name)
        if ext == ".pdb":
            records = mol.getPDBRecords()
            f = open(filename, "w")
            for rr in records:
                f.write(rr)
            f.close()
        elif ext == ".mol2":
            from MolKit2.openBabelInterface import ProdyToOBMol
            import openbabel
            self.obmol = ProdyToOBMol(mol.select())
            obconv = openbabel.OBConversion()
            obconv.SetOutFormat("mol2")
            molStr = obconv.WriteString(self.obmol)
            f = open(filename, "w")
            f.write(molStr)
            f.close()

    def focusScene(self, obj=None, viewer=None):
        #self.viewer.stopAutoRedraw()
        if not viewer:
            if self.gui is not None:
                viewer = self.gui().viewer
        viewer.suspendRedraw =True
        if obj is not None:
            if isinstance(obj, SelectionSet):
                sel = obj
            elif isinstance(obj, Selection):
                sel = [obj]
            else:
                sel = [obj.select()]
        else: # obj is None
            if self.gui is not None:
                items = self.gui().objTree.selectedItems()
                sel = []
                if len(self.curSelection):
                    sel.extend(self.curSelection)
                if len(items):
                    for item in items:
                        if isinstance(item._pmvObject, Molecule):
                            mol = item._pmvObject
                            if mol.geomContainer.masterGeom.visible:
                                sel.append(mol.select()) 
                        elif isinstance(item._pmvObject, (Chain, Residue, Atom)):
                            selection = item._pmvObject.select()
                            mol = selection.getAtomGroup().getMolecule()
                            if mol.geomContainer.masterGeom.visible:
                                sel.append(selection) 
                        elif isinstance(item._pmvObject, SelectionSet):
                            for selection in item._pmvObject:
                                mol = selection.getAtomGroup().getMolecule()
                                if mol.geomContainer.masterGeom.visible:
                                    sel.append(item._pmvObject.select())
                if not len(sel):
                    sel = [x.select() for x in self.Mols if x.geomContainer.masterGeom.visible]
            else:
                if not len(sel):
                    sel = [x.select() for x in self.Mols if x.geomContainer.masterGeom.visible]
        #print "SELECTION:", sel
        if len(sel):
            visible = []
            root = viewer.rootObject
            for ch in root.children:
                if ch.visible:
                    visible.append(ch)
                    ch.Set(visible=False)
            coords = sel[0].getCoords()
            for it in sel[1:]: 
                coords = numpy.concatenate([coords, it.getCoords()])
            from DejaVu2.Points import Points
            p = Points( "selatoms", vertices=coords.astype("f"),
                                      faces=[range(len(coords))], visible=1,
                                      inheritMaterial=1, protected=False, disableStencil=True,
                                      transparent=True)
                        
            viewer.AddObject(p)
            rootObj = viewer.rootObject

            from DejaVu2.moveGeom import MoveGeom
            mg = MoveGeom(rootObj)
            transf1 = mg.getTransformation(rootObj)
            rot = transf1[0]
            ca1 = mg.getCameraAttr()
            #print 'BEFORE', ca1
            if viewer.currentCamera.ssao: # ambient occlusion is on
                ssaoNear = 2.0
            else:
                ssaoNear = None
            viewer.Reset_cb()
            viewer.Normalize_cb()
            viewer.Center_cb()
            
            viewer.RemoveObject(p)
            transf2 = mg.getTransformation(rootObj)
            transf2[0] = rot
            ca2 = mg.getCameraAttr()
            #print 'AFTER', ca1, ca2
            if ssaoNear: # ambient occlusion is on
                ca2['near'] = ssaoNear
                
            for geom in visible:
                geom.Set(visible=True)
            viewer.suspendRedraw = False
            mg.interpolate(.5, rootObj, transf=[{rootObj:transf1}, {rootObj:transf2}],
                           camattr=[ca1, ca2])
            #rootObj.Set(rotation=rot)
            #viewer.startAutoRedraw()
            

class MolGeomContainer(GeomContainer):
    """
    Class to hold geometries used to represent molecules in a viewer.
    An instance of this class called geomContainer is added to
    each loaded Molecule.
    """
    def __init__(self, mol, app):
        """constructor of the geometry container"""

        GeomContainer.__init__(self, app)
        self.mol = mol
        mol.geomContainer = self

        ## Dictionary of AtomSets used to hold atoms for each geometry in the container
        self.atoms = {}
        ## Dictionary of geoms: {atoms:None} used to find out if an atom has a certain representation
        self.atomGeoms = {}

        ## Dictionary {geometry name: function}. Function is used to expand an atomic
        ## property to the corresponding vertices in a geometry.
        ## The function has to accept 4 arguments: a geometry name,
        ## a list of atoms,  the name of the property and an optional argument--
        ## propIndex (default is None)-- specifying the index of the property
        ## when needed.
        self.atomPropToVertices = {}

        
        ##  Dictionary {geometry name: function}. Function is used to map a vertex to an atom.
        ## If no function is registered, default mapping is used (1vertex to 1atom).
        ## If None is registered -- this geometry cannot represent atoms
        self.geomPickToAtoms = {}

        ## Dictionary {geometry name: function}. Function is used to map vertices to bonds.
        ## If no function is registered, default mapping is used (1vertex to 1bond).
        ## If None is registered: this geometry cannot represent bonds.
        self.geomPickToBonds = {}

        ## this set of coordinates should really be shared by all geometries
        self.allCoords = mol._ag._coords[mol._ag._acsi].astype('f')

        # master Geometry
        from DejaVu2.Geom import Geom
        g = self.masterGeom = Geom(mol._basename, shape=(0,0), 
                                   pickable=0, protected=True)
        self.masterGeom.isScalable = 0
        self.geoms['master'] = self.masterGeom
        self.masterGeom.replace = True

        event = AddGeometryEvent(g, parent=None, redo=False)
        self.app.eventHandler.dispatchEvent(event)
        # This dictioanary will contain commands used to display/undisplay
        #geometries. {geom.name:[displayCmd, undisplayCmd]}
        self.geomCmds = {}

    def displayedAs(self, geomNames, atoms, mode='fast'):
        # mode = 'fast' means return as soon as we find at least 1 atom
        # mode = 'complete' means after checking all atoms
        #t0 = time()
        l = len(atoms)
        if mode == 'complete':
            counts = []
            count = 0
            for geomName in geomNames:
                d = self.atomGeoms.get(geomName, None)
                if d is None: continue
                for a in atoms:
                    if d.has_key(a): count +=1
                counts.append(count)
            #print 'DISPLAYED as', geomName, l, counts, time()-t0
            return counts
        
        elif mode == 'fast':
            indices = atoms.getIndices()
            for geomName in geomNames:
                d = self.atomGeoms.get(geomName, None)
                if d is None: continue
                for a in indices:
                    if d.has_key(a):
                        #print 'DISPLAYED as', geomName, l, count, count/float(l), time()-t0
                        return 1
            return 0
            
    def setAtomsForGeom(self, geomName, atoms):
        self.atoms[geomName] = atoms
        self.atomGeoms[geomName] = {}.fromkeys(atoms.getIndices(), None)

        
    def addGeom(self, geom, parent=None, redo=False, makeAtomSet=True,
                displayCmd=None, undisplayCmd=None):
        """Add geometry to to geomContainer, create atom set and set pointer
        from geom to molecule"""
        GeomContainer.addGeom(self, geom, parent, redo)
        geom.mol = weakref.ref(self.mol)  #need for backtracking picking
        if makeAtomSet:
            if not self.atoms.has_key(geom.name):
                self.setAtomsForGeom(geom.name, self.mol.emptySelection())
        self.geomCmds[geom.name] = [displayCmd, undisplayCmd]

    ## def getGeomColor(self, geomName):
    ##     """Build a list of colors for a geometry from the atom's colors"""
    ##     if self.atomPropToVertices.has_key(geomName):
    ##         func = self.atomPropToVertices[geomName]
    ##         geom = self.geoms[geomName]
    ##         atms = self.atoms[geomName]
    ##         col = func(geom, atms, 'colors', propIndex=geomName)

    ##     else:
    ##         if geomName in self.atoms.keys():
    ##             col = [x.colors.get(geomName, x.colors['lines']) for x in self.atoms[geomName]]
    ##         else:
    ##             return

    ##     if col is not None:
    ##         colarray = numpy.array(col, 'f')
    ##         diff = colarray - colarray[0]
    ##         maxi = numpy.maximum.reduce(numpy.fabs(diff.ravel()))
    ##         if maxi==0:
    ##             return [colarray[0].tolist()]
    ##         else:
    ##             return col


    ## def updateColors(self, geomName=[], updateOpacity=0):
    ##     for name in geomName:
    ##         if geomName=='master': continue
    ##         if geomName=='selectionSpheres': continue
    ##         if self.atoms.has_key(name) and len(self.atoms[name])==0: continue 
    ##         col = self.getGeomColor(name)

    ##         if updateOpacity:
    ##             self.geoms[name].Set(materials=col, redo=1,
    ##                                  tagModified=False, transparent='implicit')
    ##             opac = self.getGeomOpacity(name)
    ##         else: opac = None
            
    ##         if col is not None and opac is not None:
    ##             self.geoms[name].Set(materials=col, opacity=opac, redo=1,
    ##                                  tagModified=False, transparent='implicit')
    ##         elif col is not None:
    ##             self.geoms[name].Set(materials=col, redo=1, tagModified=False, transparent='implicit')
    ##         elif opac is not None:
    ##             self.geoms[name].Set(opacity=opac, redo=1, tagModified=False, transparent='implicit')


    ## def getGeomOpacity(self, geomName):
    ##     if self.atomPropToVertices.has_key(geomName):
    ##         func = self.atomPropToVertices[geomName]
    ##         geom = self.geoms[geomName]
    ##         atms = self.atoms[geomName]
    ##         col = func(geom, atms, 'opacities', propIndex = geomName)
    ##     else:
    ##         if geomName in self.atoms.keys():
    ##             col = [x.opacities[geomName] for x in self.atoms[geomName]]
    ##         else:
    ##             return
    ##     if col is not None:
    ##         colarray = numpy.array(col, 'f')
    ##         diff = colarray - colarray[0]
    ##         maxi = numpy.maximum.reduce(numpy.fabs(diff.ravel()))
    ##         if maxi==0:
    ##             return colarray[0]
    ##         else:
    ##             return col


    ## def updateOpacity(self, geomName=[]):
    ##     for name in geomName:
    ##         if geomName=='master': continue
    ##         if geomName=='selectionSpheres': continue
    ##         if len(self.atoms[name])==0: continue
    ##         col = self.getGeomColor(name)
    ##         if col:
    ##             col = numpy.array(col, 'f')
    ##             self.geoms[name].Set(materials=col, redo=1, tagModified=False)


from AppFramework.AppCommands import AppCommand

class MVCommand(AppCommand):

    def onAddCmdToApp(self):
        evh = self.app().eventHandler
        evh.registerListener(DeleteObjectEvent, self.handleDeleteObject)

    def handleDeleteObject(self, event):
        try:
            del self.initializedFor[event.object]
        except KeyError:
            pass

    def __init__(self):
        AppCommand.__init__(self)
        self.firstArgIsSelectionSet = True
        self.iterateOverFirstArgument = True
        
    def expandArg0(self, obj):
        # method to be overridden by subclass
        if self.firstArgIsSelectionSet:
            return self.app().selectionSetFromObject(obj)
        else:
            return obj
    
    def getName(self, args, kw, full=False):
        name = self.name + ' '
        if full:
            for arg in args: name += repr(arg)
            name += repr(kw)
        else:
            if len(args):
                if isinstance(args[0], list):
                    for a in args[0]:
                        if isinstance(a, Selection):
                            name += '{0} atoms from {1}'.format(len(a), a.getAtomGroup().getTitle())
                        else:
                            name += repr(a)
                else:
                    name += repr(args[0])
        return name
    
    def _strArg(self, arg):
        """
        Method to turn a command argument into a string for logging purposes
        Add support for TreeNodes and TreeNodeSets
        """
        #if type(arg)==types.InstanceType:
        from mglutil.util.misc import isInstance
        if isInstance(arg) is True:
            if issubclass(arg.__class__, TreeNode):
                return '"' + arg.full_name() + '", ', None
                
            if issubclass(arg.__class__, TreeNodeSet):
                stringRepr = arg.getStringRepr()
                if stringRepr:
                    return '"' + stringRepr + '", ', None
                else:
                    name = ""
                    mols, elems = self.app().getNodesByMolecule(arg)
                    for elem in elems:
                        name = name + elem.full_name() +';'
                    return '"' + name + '", ', None

        return AppCommand._strArg(self, arg)
