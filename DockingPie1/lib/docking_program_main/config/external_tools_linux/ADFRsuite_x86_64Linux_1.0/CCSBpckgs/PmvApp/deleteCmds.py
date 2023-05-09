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
# Author: Ruth Huey, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2014
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/PmvApp/deleteCmds.py,v 1.8.4.1 2017/07/13 20:55:28 annao Exp $
#
# $Id: deleteCmds.py,v 1.8.4.1 2017/07/13 20:55:28 annao Exp $
#

"""
This Module implements commands to delete items from the MoleculeViewer:
for examples:
    Delete Molecule
"""
from PmvApp.Pmv import Event, AfterDeleteAtomsEvent, DeleteAtomsEvent, \
     BeforeDeleteMoleculeEvent, AfterDeleteMoleculeEvent, RefreshDisplayEvent
from PmvApp.Pmv import BeforeDeleteMoleculeEvent, AfterDeleteMoleculeEvent
from PmvApp.selectionCmds import SelectionEvent
from PmvApp.Pmv import MVCommand

from AppFramework.App import DeleteObjectEvent
from MolKit2.selection import Selection, SelectionSet
from MolKit2.molecule import Molecule, MoleculeSet
from PmvApp.selectionCmds import SelectionEvent



class DeleteMolecule(MVCommand):
    """Command to delete molecule from the MoleculeViewer \n
    Package : PmvApp \n
    Module  : deleteCmds \n
    Class   : DeleteMolecule \n
    Command : deleteMolecule \n
    Synopsis:\n
        None<---deleteMolecule(molSel, **kw) \n
    """

    def checkArguments(self, molSel):
        """None <- deleteMolecules(nodes, **kw) \n
        nodes: TreeNodeSet holding the current selection.
        """
        assert molSel
        return (molSel,), {} 
   
    def doit(self, molSel):
        mol = molSel.getAtomGroup().getMolecule()
        # manage selection
        if self.app().activeSelection.nbAtoms():
            curSel = self.app().activeSelection
            for sele in curSel:
                if sele.getAtomGroup() == mol._ag:
                    self.app().select(sele, negate=True, setupUndo=False)
                    # we do not need to create undo for this command
                    self.app().clearTopUndoCmds() # this is to prevent adding the contents
                    # of app()._topUndoCmds list to the undo stack by the top command (deleteCommand)
                    break
            #molSel = SelectionSet([mol.select()])
            #molSel = SelectionSet([molSel])
            #newcursel = curSel - molSel
        for name, selSet in self.app().namedSelections.items():
            newselSet = selSet - SelectionSet([molSel])
            self.app().namedSelections[name] = newselSet
        #
        self.app().removeObject(mol, "Molecule")
        event = DeleteObjectEvent(object=mol, objectType='Molecule')
        self.app().eventHandler.dispatchEvent(event)
        self.updateUndoRedoCmdStack(mol, self.app().undo)
        self.updateUndoRedoCmdStack(mol, self.app().redo)
        mol._treeItems = {}
        del mol

    def updateUndoRedoCmdStack(self, molecule, command):
        # cmd is either undo or redo
        # remove all entries form the undo or redo  command stack that contain the deleted molecule
        remove = []
        for i, cmdlist in enumerate(command.cmdStack):
            # loop over undo/redo items in the stack. Each item is a tuple containing a list of
            # commands to execute an undo/redo action and a string describing the action:
            # ( [ (cmd, (args,) {kw}), (cmd, (args,) {kw}), ... ], "undo/redo name" )
            #import pdb; pdb.set_trace()
            newName = cmdlist[1].split(" ")[0]
            ncmds = len(cmdlist[0])
            for j in xrange(len(cmdlist[0])-1, -1, -1): # loop backwards over the list of commands in undo/redo action 
                cmd, args, kw = cmdlist[0][j]
                if len(args):
                    newargs0 = [] # this list will contain updated atom selections
                    newargs = []  # this list will contain other arguments that are not atomsets or molecules
                    for arg in args:
                        if isinstance(arg, list) and len(arg) and isinstance(arg[0], Molecule):
                            mollist = []
                            for mol in arg:
                                if mol != molecule: mollist.append(mol)
                            if len(mollist):
                                newargs0.append(MoleculeSet("molset", mollist))
                        elif isinstance(arg, Selection):
                            if arg.getAtomGroup() != molecule._ag: newargs0.append(arg)
                        elif isinstance(arg, Molecule):
                            if arg == molecule: newargs0.append(arg)
                        else: newargs.append(arg)

                    if len(newargs0):
                        cmdlist[0][j] = (cmd, tuple(newargs0+newargs), kw)
                        name = self.getName(newargs0, kw)
                        ind = name.find(" ") # "deleteMolecule " + "...."
                        if ind and len(name)>ind+1:
                            newName = newName+name[ind:]
                    else:
                        #remove this command from the list
                        del cmdlist[0][j]
                    
            if len(cmdlist[0]) == 0:
                remove.append(i)
            else:
                if ncmds != len(cmdlist[0]):
                   command.cmdStack[i] = (cmdlist[0], newName)
        n = 0
        for i in remove:
            command.cmdStack.pop(i-n)
            n = n+1

        from AppFramework.notOptionalCommands import AfterUndoEvent
        event =  AfterUndoEvent(objects=command.cmdStack, command=command)
        self.app().eventHandler.dispatchEvent(event)
        
     
  

class DeleteAtoms(MVCommand):
    """ Command to remove an AtomSet from the MoleculeViewer  \n
    Package : PmvApp \n
    Module  : deleteCmds \n
    Class   : DeleteAtoms \n
    Command : deleteAtoms \n
    Synopsis:\n
        None<---deleteAtoms(atoms)\n
    Required Arguments:\n
        atoms --- AtomSet to be deleted.\n
    """
    
    def __init__(self):
        MVCommand.__init__(self)


    def onAddCmdToApp(self):
        if not self.app().commands.has_key('deleteMolecule'):
            self.app().lazyLoad("deleteCmds", commands=["deleteMolecule"], package="PmvApp")
        self.app().lazyLoad("selectionCmds", commands=[
            "select", "clearSelection"], package="PmvApp")


    def checkArguments(self, selection):
        """None <- deleteAtoms(atoms) \n
        atoms: AtomSet to be deleted."""
        assert selection
        return (selection,), {}

    def doit(self, selection):
        """ Function to delete all the references to each atom  of a
        AtomSet."""

        # Remove the atoms of the molecule you are deleting from the
        # the AtomSet self.app().allAtoms
        mol = selection.getAtomGroup().getMolecule()
        app = self.app()
        if len(selection)==len(mol._ag):
            app.deleteMolecule(mol)
            return

        # delete the atoms
        mol.deleteAtoms(selection)
        event = DeleteAtomsEvent(deletedAtoms=selection)
        self.app().eventHandler.dispatchEvent(event)
        event = RefreshDisplayEvent(molecule=mol)
        self.app().eventHandler.dispatchEvent(event)
        
    def updateUndoRedoCmdStack(self, oldAg,  mol, indMap, command):
        remove = []
        for i, cmdlist in enumerate(command.cmdStack):
            # loop over undo/redo items in the stack. Each item is a tuple containing a list of
            # commands to execute an undo/redo action and a string describing the action:
            # ( [ (cmd, (args,) {kw}), (cmd, (args,) {kw}), ... ], "undo/redo name" )
            #import pdb; pdb.set_trace()
            cmdName = cmdlist[1].split(" ")[0]# name of the undo command without the string representing its arguments
            newName = ""
            ncmds = len(cmdlist[0])
            for j in xrange(len(cmdlist[0])-1, -1, -1): # loop backwards over the list of commands in undo/redo action 
                cmd, args, kw = cmdlist[0][j]
                if len(args):
                    newargs0 = [] # this list will contain new atom selections
                    newargs = []  # this list will contain arguments that are not Molecules ar atom selections
                    for arg in args:
                        if isinstance(arg, Selection):
                            if arg.getAtomGroup() == oldAg:
                                inds = arg.getIndices()
                                import numpy
                                _newinds = numpy.where(indMap[inds]> -1)[0]
                                if len (_newinds):
                                    newinds = indMap[inds][_newinds]
                                    # create  a selection string:
                                    stformat = "index " + "%s " *len(newinds)
                                    selstr = stformat%tuple(newinds)
                                    newAt = mol.select(selstr)
                                    newargs0.append(SelectionSet([newAt], ""))
                        elif isinstance(arg, list) and isinstance(arg[0], Molecule):
                            newargs0.append(arg)
                        elif isinstance(arg, Molecule):
                            newargs0.append(arg)
                        else: newargs.append(arg)
                    if len(newargs0):
                        cmdlist[0][j] = (cmd, tuple(newargs0+newargs), kw)
                        # create string representing selection. Can not use cmd.getName()
                        # since cmd can be a commandLoader
                        if newName == "":
                            name = self.getName(newargs0, kw) #this returns "deleteAtoms " + "...."
                            ind = name.find(" ") # get name after the " "
                            if ind and len(name)>ind+1:
                                newName = name[ind:]
                    else:
                        #remove this command from the list
                        del cmdlist[0][j]
            if len(cmdlist[0]) == 0:
                remove.append(i)
            else:
                if newName != "":
                    command.cmdStack[i] = (cmdlist[0], cmdName + " " + newName)
        n = 0
        for i in remove:
            command.cmdStack.pop(i-n)
            n = n+1

        from AppFramework.notOptionalCommands import AfterUndoEvent
        event =  AfterUndoEvent(objects=command.cmdStack, command=command)
        self.app().eventHandler.dispatchEvent(event)
        
     
                    

class DeleteCurrentSelection(DeleteAtoms):
    """ Command to remove an AtomSet from the MoleculeViewer \n
    Package : PmvApp \n
    Module  : deleteCmds \n
    Class   : DeleteCurrentSelection \n
    Command : deleteCurrentSelection \n
    Synopsis:\n
        None<---deleteCurrentSelection()\n
    Required Arguments:\n
        None
    """
    
    def checkArguments(self):
        """None <- deleteCurrentSelection()
        """
        return (), {}


    def doit(self):
        """ Function to delete all the references to each atom of a
        the currentSelection."""

        atoms = old = self.app().activeSelection.get()[:] # make copy of selection
        ats = self.app().expandNodes(atoms)
        if not len(ats):
            return
        ats = ats.findType(Atom)

        self.app().clearSelection(log=0)
        DeleteAtoms.doit(self, ats)



class DeleteHydrogens(DeleteAtoms):
    """ Command to remove hydrogen atoms from the MoleculeViewer  \n
    Package : PmvApp
    Module  : deleteCmds  \n
    Class   : DeleteHydrogens  \n
    Command : deleteHydrogens  \n
    Synopsis:\n
        None<---deleteHydrogens(atoms)\n
    Required Arguments:\n 
        atoms --- Hydrogens found in this atom set are deleted.\n
    """

    def checkArguments(self, atoms):
        """None <- deleteHydrogents(atoms) \n
        atoms: set of atoms, from which hydrogens are to be deleted."""
        if isinstance(atoms, str):
            self.nodeLogString = "'"+atoms+"'"
        
        ats = self.app().expandNodes(atoms)
        assert ats
        return (ats,), {}


    def doit(self, ats):
        hatoms = ats.get(lambda x: x.element=='H')
        if not len(hatoms):
            self.app().warningMsg("No hydrogens to delete.")
            return
        DeleteAtoms.doit(self, hatoms)


commandClassFromName = {
    'deleteMolecule' : [DeleteMolecule,  None],
    'deleteAtoms' : [DeleteAtoms,  None ],
    #'deleteMolecules' : [DeleteMolecules,  None],
    #'deleteAllMolecules' : [DeleteAllMolecules, None],
    
    #'deleteCurrentSelection' : [DeleteCurrentSelection, None ],
    #'deleteHydrogens' : [DeleteHydrogens,  None],
    #'restoreMol' : [RestoreMolecule,  None],
}


def initModule(viewer, gui=True):
    for cmdName, values in commandClassFromName.items():
        cmdClass, guiInstance = values
        if gui:
            viewer.addCommand(cmdClass(), cmdName, guiInstance)
        else:
            viewer.addCommand(cmdClass(), cmdName, None)

