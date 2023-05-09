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
# Author: Michel F. SANNER, Sophie COON, Anna Omelchenko
#
# Copyright: M. Sanner TSRI 2014
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/PmvApp/bondsCmds.py,v 1.5.4.1 2017/07/13 20:55:28 annao Exp $
#
# $Id: bondsCmds.py,v 1.5.4.1 2017/07/13 20:55:28 annao Exp $
#
#
import numpy
from opengltk.OpenGL import GL

from PmvApp.Pmv import EditAtomsEvent
from PmvApp.Pmv import MVCommand  #, MVAtomICOM, MVBondICOM
from DejaVu2.Geom import Geom
from DejaVu2.Spheres import Spheres

from DejaVu2.IndexedPolylines import IndexedPolylines
from MolKit2.molecule import Atom, Chain

var=0
"""
Package: PmvApp
Module : bondsCmds
This module provides a set of commands:
- BuildBondsByDistance (buildBondsByDistance) to compute the bonds between
  atoms of the given nodes and connect the residues when necessary.
- AddBondsCommands (addBonds) to create a bond instance between two given atoms.
    -no gui
- AddBondsGUICCommands (addBondsGC) to create a bond between two picked
  atoms. If a group of atoms is selected by dragging, it will
  buildBondsByDistance between them.  
  This is an interactive command (ICOM) and the GUI version of addBonds.
- RemoveBondsCommands (removeBonds) to delete an existing bonds between the
  two given atoms.
    -no gui
- RemoveBondsGUICCommands (removeBondsGC) to delete the picked bond. If a
  groups of atoms is selected by dragging, it will remove all the bonds
  between them.  This is an interactive command (ICOM) , provides the GUI for removeBonds command.

The menubuttons for BondsCommands are located under the 'Bonds' entry in the 'Edit' menu.
"""

class BuildBondsByDistance(MVCommand): #, MVAtomICOM):
    """This command creates the bonds between atoms of the given set of nodes and connect the residues when necessary. The cut_off distance is based on the bond order radius of the two atoms to be bound. The command will use the bhtree functionality to find the pairs of atoms to be found when available.This command can be applied on a current selection or used interactively when bound to the mouse event.  \n
   Package : PmvApp  \n
   Module  : bondsCmds  \n
   Class   : BuildBondsByDistance  \n
   Command : buildBondsByDistance  \n
   Synopsis:\n
        None<---buildBondsByDistance( nodes, display=True)\n
   Required Arguments:\n     
        nodes --- any set for MolKit2.Selection describing molecular components\n
        display --- when set to 1 the displayLines command is called and th
                 bonds are displayed as lines.\n
   Optional Arguments:\n        
        kw --- any of the additional optional keywords supported by commands
               (see ViewerFramework.VFCommands.Command documentation).\n
    Required Commands:\n
      displayLines 
    Examples:\n
      mol = mv.Mols[ 0 ]\n
      mv.buildBondsByDistance(mol, display=1)\n
    """

    def __init__(self):
        MVCommand.__init__(self)
        #MVAtomICOM.__init__(self)
        #self.flag = self.flag | self.objArgOnly
        

    def undoCmdAfter(self, *args, **kw):
        if len(self.createdBonds):
            return ([(self.app().removeBonds, (self.createdBonds,), {})],
                    self.app().removeBonds.name)
        

    ## def onAddCmdToApp(self):
    ##     #self.app().lazyLoad('displayCmds', commands=['displayLines'], package='PmvApp')

    ##     # try to import bhtree package and set the command flag hasBhtree
    ##     # to the proper value.
    ##     try:
    ##         import bhtree
    ##         self.hasBhtree = 1
    ##     except:
    ##         self.hasBhtree = 0


    def doit(self, nodes, display=False, **kw):
        self.createdBonds = None
        # Get the molecules and the atomSets per molecules in the set of nodes.
        #molecules, atomSets = self.app().getNodesByMolecule(nodes, Atom)
        # Loop over the molecules and the sets of atoms
        bonds = []
        for sel in nodes:
            mol = sel.getAtomGroup().getMolecule()
            bonds = mol.buildBondsByDistance()

        if display:
            # when we call it from the GUI we want the result to be displayed
            self.app().displayLines(nodes, setupUndo=0)
        self.createdBonds = bonds
        return bonds


    def checkArguments(self, nodes, display=False, redraw=1):
        """None <- buildBondsByDistance(nodes, display=1) \n
        nodes   : TreeNodeSet holding the current selection  \n
        display : Boolean flag if set to True the displayLines  \n
                  commands are called. \n
        The buildBondsByDistance commands connects atoms and residues in the \n
        given set of nodes based on the covalent radius of each atoms to be
        connected. \n
        """
        #if isinstance(nodes, str):
        #    self.nodeLogString = "'"+nodes+"'"
        kw = {}
        kw['redraw']=redraw
        kw['display']=display
        nodes = self.app().expandNodes(nodes)

        return (nodes,) , kw


class AddBondsGUICommand(MVCommand): #, MVAtomICOM):
    """
    The AddBondGUICommand provides an interactive way of creating bonds between two given atoms by picking on them. To use this command you need first to load it into PMV. Then you can find the entry 'addBonds' under the Edit menu. To add bonds  you just need to pick on the 2 atoms you want to bind. If you drag select  a bunch of atoms, the command will buildBondsByDistance between them.This command is undoable.\n
   Package : PmvApp \n
   Module  : bondsCmds \n
   Class   : AddBondsGUICommand \n
   Command : addBondsGC \n
   Synopsis:\n
        None<-addBondsGC(atoms)\n
    Required Arguments:\n    
        atoms  : atom(s)\n
    """
    
    def __init__(self):
        MVCommand.__init__(self)
        #MVAtomICOM.__init__(self)
        self.atomList = AtomSet([])
        self.labelStrs = []
        self.createdBonds = [] # this list will store created bonds used to create
                               # negation of command

            
    def onRemoveObjectFromViewer(self, obj):
        removeAts = AtomSet([])
        for at in self.atomList:
            if at in obj.allAtoms:
                removeAts.append(at)
        self.atomList = self.atomList - removeAts
        removeAts = AtomSet([])
        self.update()

       
    def onAddCmdToApp(self):
        #if not self.app().commands.has_key('setICOM'):
        #    self.app().lazyLoadCommands('interactiveCommands', ['setICOM'], 'Pmv')

        self.app().lazyLoad('bondsCmds', commands=[
            'addBonds', 'removeBondsGC'], package='PmvApp')

        self.masterGeom = Geom('addBondsGeom',shape=(0,0), 
                               pickable=0, protected=True)
        self.masterGeom.isScalable = 0
        self.spheres = Spheres(name='addBondsSpheres', shape=(0,3),
                               inheritMaterial=0,
                               radii=0.2, quality=15,
                               materials = ((1.,1.,0.),), protected=True) 
        if not self.app().commands.has_key('labelByExpression'):
            self.app().lazyLoad('labelCmds', commands=['labelByExpression',],
                                    package='PmvApp')
        from AppFramework.App import AddGeometryEvent


        # The AppGUI should register addGeom() method with
        # AddGeometryEvent.
        #This method will do the following:
        #miscGeom = self.app().GUI.miscGeom
        #self.app().GUI.VIEWER.AddObject(self.masterGeom, parent=miscGeom)
        #self.app().GUI.VIEWER.AddObject(self.spheres, parent=self.masterGeom)
        #There should be a way to specify
        # self.app().GUI.miscGeom parent here (???)
                
        miscGeom = None
        event1 = AddGeometryEvent(self.masterGeom, parent=miscGeom)
        self.app().eventHandler.dispatchEvent(event1)
        event2 = AddGeometryEvent(self.spheres, parent=self.masterGeom)
        self.app().eventHandler.dispatchEvent(event2)

        

    def checkArguments(self, atoms, **kw):
        """None<-addBondsGC(atoms)
           \natoms  : atom(s)"""
        if isinstance(atoms, str):
            self.nodeLogString = "'"+atoms+"'"
        ats = self.app().expandNodes(atoms)
        assert isinstance(ats, AtomSet)
        assert len(ats)
        return (ats,), kw


    def doit(self, ats):
        if len(ats)>2:
            if len(self.atomList):
                atSet = ats + self.atomList
            else: atSet = ats
            parent = atSet[0].parent
            #parent.buildBondsByDistanceOnAtoms(atSet)
            self.app().buildBondsByDistance(atSet)
            self.update(True)
            self.atomList = AtomSet([])
            self.app().displayLines(atSet, setupUndo=0)
        else:
            lenAts = len(self.atomList)
            last = None
            if lenAts:
                last = self.atomList[-1]
                top = self.atomList[0].top
            for at in ats:
                #check for repeats of same atom
                if lenAts and at==last:
                    continue
                #lenAts = len(self.atomList)
                #if lenAts and at==self.atomList[-1]:
                #    continue
                if lenAts and at.top!=self.atomList[-1].top:
                    msg = "intermolecular bond to %s disallowed"%(at.full_name())
                    self.app().warningMsg(msg)
                self.atomList.append(at)
                lenAts = len(self.atomList)
            self.update(True)
            #if only have one atom, there is nothing else to do
            if lenAts<2: return
            #now build bonds between pairs of atoms
            atSet = self.atomList
            if lenAts%2!=0:
                atSet = atSet[:-1]
                #all pairs of atoms will be bonded
                #so keep only the last one
                self.atomList = atSet[-1:]
                lenAts = lenAts -1
            else:
                self.app().labelByExpression(self.atomList, negate=1, setupUndo=0)
                self.atomList = AtomSet([])
            for i in range(0, lenAts, 2):
                at1 = atSet[i]
                at2 = atSet[i+1]
                self.app().addBonds( [(at1, at2)], bondOrder=[1], origin='UserDefined',log=0)
        self.update(True)


    def applyTransformation(self, pt, mat):
        pth = [pt[0], pt[1], pt[2], 1.0]
        return numpy.dot(mat, pth)[:3]


    def getTransformedCoords(self, atom):
        if not atom.top.geomContainer:
            return atom.coords
        g = atom.top.geomContainer.geoms['master']
        c = self.applyTransformation(atom.coords, g.GetMatrix(g))
        return  c.astype('f')


    def update(self, event=None):
        if not len(self.atomList):
            self.spheres.Set(vertices=[], tagModified=False)
            self.app().labelByExpression(self.atomList, negate=1, setupUndo=False)
            from AppFramework.App import RedrawEvent
            event = RedrawEvent(self.app())
            self.app().eventHandler.dispatchEvent(event)
            #self.app().GUI.VIEWER.Redraw()
            return
        self.lineVertices=[]
        #each time have to recalculate lineVertices
        for at in self.atomList:
            c1 = self.getTransformedCoords(at)
            self.lineVertices.append(tuple(c1))
            
        if event:
            self.spheres.Set(vertices=self.lineVertices, tagModified=False)
            self.app().labelByExpression(self.atomList,
                                      function = 'lambda x: x.full_name()',
                                      lambdaFunc = 1,
                                      textcolor = 'yellow',
                                      format = '', negate = 0,
                                      location = 'Last', log = 0,
                                      font = 'arial1.glf', only = 1,
                                      setupUndo=False )
        #setting spheres doesn't trigger redraw so do it explicitly
        # The AppGUI should register a listener for the redraw event. The method will do:
        # self.app().GUI.VIEWER.Redraw()
        from AppFramework.App import RedrawEvent
        event = RedrawEvent(self.app())
        self.app().eventHandler.dispatchEvent(event)
        
        
    ## def startICOM(self):
    ##     self.app().setIcomLevel( Atom)

    ## def stopICOM(self):
    ##     if len(self.atomList)!=0:
    ##         self.app().labelByExpression(self.atomList, negate=1, setupUndo=False )
    ##     del self.atomList[:]
    ##     self.labelStrs = []
    ##     self.spheres.Set(vertices=[], tagModified=False)
    ##     self.app().GUI.VIEWER.Redraw()
    ##     self.save = None


class AddBondsCommand(MVCommand):
    """
    The AddBondsCommand is a command to connect two given atoms.This command will not allow the creation inter-molecular bonds.This command doesn't have a GUI interface and can only be called through the pyShell.
  Package : PmvApp  \n
   Module  : bondsCommands  \n
   Class   : AddBondsCommand  \n
   Command : addBonds  \n
   Synopsis:\n
        None<-addBonds(atom1,atom2)\n
   Required Arguments:\n    
     atom1  : first atom\n
     atom2  : second atom \n
   Optional Arguments:\n  
     bondOrder : Integer specifying the bondOrder of the bond that is going to be created between atom1 and atom2.The bondOrder by default is 1.\n
     origin  : string describing how bond was specified \n 
    """

    def onAddCmdToApp(self):
        if not self.app().commands.has_key('removeBonds'):
            self.app().lazyLoad('bondsCmds', commands=['removeBonds'], package='PmvApp')


    def undoCmdAfter(self, *args, **kw):
        if self.createdBonds:
            return ([(self.app().removeBonds, (self.createdBonds,), {})],
                    self.app().removeBonds.name)
        

    def checkArguments(self, bondAtomPairs, bondOrder=None, origin='UserDefined'):
        """None<-addBonds(bondAtomPairs)
           \nbondAtomPairs  : list of pair of atoms forming bonds
           \nbondOrder : list of integer specifying the bondOrder of the bonds that is going to be created between. The bondOrder by default is 1.
           \norigin  : string describing how bond was specified"""
        assert isinstance(bondAtomPairs, (list, tuple))
        assert len(bondAtomPairs)
        if bondOrder is None:
            bondOrder = [1]*len(bondAtomPairs)
        else:
            assert len(bondOrder)==len(bondAtomPairs)

        return (bondAtomPairs,), {'bondOrder':bondOrder,'origin':origin}


    def doit(self, bondAtomPairs, bondOrder, origin):
        bonds = self.createdBonds = []
        for atoms, bo in zip(bondAtomPairs, bondOrder):
            ats = self.app().expandNodes(atoms[0])
            if not len(ats):
                self.app().warningMsg("No atoms found to add bonds") 
                return 'ERROR'
            atom1 = ats[0]
            ats = self.app().expandNodes(atoms[1])
            if not len(ats):
                self.app().warningMsg("No atoms found to add bonds")
                return 'ERROR'
            atom2 = ats[0]
            # check for inter molecular bonds
            if atom1.top!=atom2.top:
               msg = "intermolecular bond between %s and %s disallowed"%(atom1.full_name(), atom2.full_name())
               self.app().warningMsg(msg)
               continue
            # do not duplicate existing bonds
            bnds = AtomSet([atom1, atom2]).bonds[0]
            if len(bnds): 
               msg = "bond between %s and %s already exists"%(atom1.full_name(), atom2.full_name())
               self.app().warningMsg(msg)
               continue
            bonds.append( Bond( atom1, atom2, origin=origin, bondOrder=bo) )

            event = EditAtomsEvent(fields=['coords'], atoms=AtomSet([atom1, atom2]))
            self.app().eventHandler.dispatchEvent(event)
        return bonds
    

class RemoveBondsGUICommand(MVCommand): #, MVBondICOM):
    """ The RemoveBondsGUICommands provides an interactive way of deleting picked bonds. To use this command you need to first load it in the application. Once loaded you will find an entry called 'delete bonds' under the bonds entry in the Edit menu. You can then pick on the bonds you wish to delete. This command is undoable.
   Package : PmvApp  \n
   Module  : bondsCmds  \n
   Class   : RemoveBondsGUICommand  \n
   Command : removeBondsGC  \n
   Synopsis:\n
        None<---removeBondsGC(bonds)\n
   Required Arguments:\n    
        bonds  : bond(s)\n
    """

    def onAddCmdToApp(self):
        if not hasattr(self.app(), 'addBondsGC'):
            self.app().lazyLoad("bondsCmds", commands=["addBondsGC"], package="PmvApp")

        if not self.app().commands.has_key('removeBonds'):
            self.app().lazyLoad('bondsCmds', commands=['removeBonds'], package='PmvApp')


    def __init__(self):
        MVCommand.__init__(self)
        #MVBondICOM.__init__(self)
        self.pickLevel = 'parts'
        self.undoBondList = []


    def getObjects(self, pick):
        for o, val in pick.hits.items(): #loop over geometries
            primInd = map(lambda x: x[0], val)
            g = o.mol.geomContainer
            if g.geomPickToBonds.has_key(o.name):
                func = g.geomPickToBonds[o.name]
                if func: return func(o, primInd)
            else:
                l = []
                bonds = g.atoms[o.name].bonds[0]
                if not len(bonds): return BondSet()
                for i in range(len(primInd)):
                    l.append(bonds[int(primInd[i])])
                return BondSet(l)


    def dismiss(self):
        #self.app().setICOM(self.save, modifier="Shift_L", setupUndo=False)
        self.save = None

    
    def checkArguments(self, bonds, **kw):
        """None <- removeBondsGC(bonds, **kw)
           \nbonds: bonds
           """
        if bonds :
            assert len(bonds)
        return (bonds,), kw


    def doit(self, bonds):
        if not bonds: return
        global var
        var=1 
        ats = AtomSet([])
        for bond in bonds:
            ats.append(bond.atom1)
            ats.append(bond.atom2)
            self.app().removeBonds( [bond] )
        var=0
        # The AppGUI should register a listener for the redraw event.
        # The method will do:
        # self.app().GUI.VIEWER.Redraw()
        from AppFramework.App import RedrawEvent
        event = RedrawEvent(self.app())
        self.app().eventHandler.dispatchEvent(event)

        

    def undoCmdBefore(self, bonds):
        cmds =([], self.app().addBondsGC.name)
        for b in bonds:
            nameStr = AtomSet([b.atom1, b.atom2]).full_name()
            cmds[0].append((self.app().addBondsGC, (nameStr,),{}) )
        return cmds

    

class RemoveBondsCommand(MVCommand):
    """ 
    The RemoveBondsCommand allows the user to remove the bond existing between two given atoms (atom1 and atom2). This command doesn't have a gui and therefore can only be called through the pyShell. \n
   Package : PmvApp \n
   Module  : bondsCommands \n
   Class   : RemoveBondsCommand \n
   Command : removeBonds \n
   Synopsis:\n
        None<---removeBonds(bondList)\n
    Required Arguments:\n    
        bondList  : list of bonds to remove\n
    """

    def onAddCmdToApp(self):
        if not hasattr(self.app(), 'addBonds'):
            self.app().lazyLoad("bondsCmds", commands=["addBonds"], package="PmvApp")


    def undoCmdBefore(self, bondList):
        # The undo is to recreate the bond you just deleted.when
        # removebondsCommand called
        bondAtomPairs = []
        for bond in bondList:
            bondAtomPairs.append( (bond.atom1, bond.atom2) )
            
        return ( [(self.app().addBonds, (bondAtomPairs,), {})], self.app().addBonds.name)


    def checkArguments(self, bondList, **kw):
        """None<-removeBonds(bondList) \n
        abondList  : list of bonds to remove"""

        assert len(bondList)
        return (bondList,), kw
        

    def doit(self, bondList):
        #Have to find the bond first
        
        for theBond in bondList:
            atom1 = theBond.atom1
            atom2 = theBond.atom2
            #remove this bond
            atom2.bonds.remove(theBond)
            #this may not be possible
            atom1.bonds.remove(theBond)
##             atom1.parent.hasBonds=0
##             if atom2.parent!=atom1.parent:
##                 atom2.parent.hasBonds = 0
##                 if atom1.parent.parent:
##                     atom1.parent.parent.hasBonds = 0
##                 if atom2.parent.parent:
##                     atom2.parent.parent.hasBonds = 0
##                 break

##         if not theBond:
##             from warnings import warn
##             warn('bond not found %s-%s'%(atom1, atom2))
##             return 'ERROR'
##         else:
        event = EditAtomsEvent('coords', AtomSet([atom1, atom2]))
        self.app().eventHandler.dispatchEvent(event)

commandClassFromName = {
    'buildBondsByDistance' : [BuildBondsByDistance, None],
    'addBonds' : [AddBondsCommand, None],
    'addBondsGC' : [AddBondsGUICommand, None],
    'removeBonds' : [RemoveBondsCommand, None],
    'removeBondsGC' : [RemoveBondsGUICommand, None],
}


def initModule(viewer):
    for cmdName, values in commandClassFromName.items():
        cmdClass, guiInstance = values
        viewer.addCommand(cmdClass(), cmdName, guiInstance)

