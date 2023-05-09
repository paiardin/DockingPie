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
# Author: Michel F. SANNER and Ruth HUEY
#
# Copyright: M. Sanner TSRI 2014
#
#############################################################################


# $Header: /mnt/raid/services/cvs/PmvApp/selectionCmds.py,v 1.10.4.1 2017/07/13 20:55:28 annao Exp $
#
# $Id: selectionCmds.py,v 1.10.4.1 2017/07/13 20:55:28 annao Exp $
#
from mglutil.events import Event
from Pmv import formatName
from MolKit2 .selection import SelectionSet

class SelectionEvent(Event):
    pass

class RefreshSelectionEvent(Event):
    pass

import time

import numpy

from PmvApp.Pmv import MVCommand, Selection #, MVAtomICOM
from MolKit2.molecule import Atom, Molecule, MoleculeSet, Residue, Chain
#from MolKit2.protein import Residue, Chain, ResidueSet, ChainSet, Protein, ProteinSet
#from MolKit2.stringSelector import StringSelector


from DejaVu2.Spheres import Spheres
from DejaVu2.Cylinders import Cylinders
from DejaVu2.IndexedPolygons import IndexedPolygons

from PmvApp import Pmv
if not hasattr( Pmv, 'numOfSelectedVerticesToSelectTriangle') is False:
    Pmv.numOfSelectedVerticesToSelectTriangle = 1


levels = ['Molecule', 'Chain', 'Residue', 'Atom']
class MVSetSelectionLevel(MVCommand):

    def __init__(self):
        MVCommand.__init__(self)
        self.levelDict={"Molecule":Molecule,
                        "Protein":Molecule,
                        "Chain":Chain,
                        "Residue":Residue,
                        "Atom":Atom}
        self.levelColors = {
            'Atom':'yellow',
            'Residue':'green',
            'Chain':'cyan',
            'Molecule':'red',
            'Protein':'red',
            }
        self.levelVar = "Molecule"


    def onAddCmdToApp(self):
        self.app().levelVar = self.levelVar

        
    def undoCmdBefore(self, Klass, KlassSet=None):
        if not (self.app().selectionLevel in [Molecule, Protein] and \
           Klass in [Molecule, Protein]) and \
           self.app().selectionLevel!= Klass:
            if not self.app().activeSelection.empty():
                # if there is a selection that will expand we add undoing
                # the selection expansion to the undo stack
                
                return ( [ (self.app().select, (self.app().activeSelection.atoms.copy(),),
                            {'only':1} ),
                           (self, (self.app().selectionLevel,),
                            {'KlassSet':None,} ) ],
                         self.name )
            else:
                return ([(self, (self.app().selectionLevel,),
                          {'KlassSet':None})], self.name )
            
    
    def doit(self, Klass, KlassSet=None):
        if isinstance(Klass, str):
            if Klass in self.levelDict.keys():
                Klass = self.levelDict[Klass]
            else:
                msg = Klass + "string does not map to a valid level"
                self.warningMsg(msg)
                return "ERROR"
        if Klass is Protein:
            Klass = Molecule
        
        ## if not self.app().activeSelection.empty():
        ##     self.app().setSelection(self.app().activeSelection.get(Klass).uniq())
        ## else:
        ##     if Klass==Molecule: self.app().setSelection(MoleculeSet([]))
        ##     elif Klass==Chain: self.app().setSelection(ChainSet([]))
        ##     elif Klass==Residue: self.app().setSelection(ResidueSet([]))
        ##     elif Klass==Atom: self.app().setSelection(AtomSet([]))

        self.app().selectionLevel = Klass
        event = selectionFeedbackEvent(objects=self.app().activeSelection.get(), highlightSelection=False)
        self.app().eventHandler.dispatchEvent(event)
        self.levelVar = (self.app().selectionLevel.__name__)
            

    def checkArguments(self, Klass, KlassSet=None, **kw):
        """selectionLevel <- setSelectionLevel(Klass, KlassSet=None, **kw)
        set the current selection level and promotes current selection to
        this level
        """
        kw['KlassSet'] = KlassSet
        return (Klass,), kw



class MVSelectCommand(MVCommand): #, MVAtomICOM):
    """Class for modfying the current selection in a molecule Viewer. \n
    selection --- a TreeNodeSet holding the current selection. Modfied by SubClasses implementing a specific selection operation \n
    level --- level at which selection occurs \n
    Package : PmvApp \n
    Module  : selectionCmds \n
    Class   : MVSelectCommand \n
    Command : select \n
    Synopsis:\n
    currentSelection <- select(nodes, negate=False, only=False, klass=None) \n
    Required Arguments:\n
        nodes --- can be a string, a TreeNode or a treeNodeSet\n
    Optional Arguments:\n    
        negate --- is 1 nodes are removed from current selection \n
        only  --- is  1 the current selection is replaced by the nodes \n
        klass ----is omitted class of objects in nodes is used, else \n
                  nodes is converted to a set of objects of type Klass and\n
                  current selection level is set to klass
    """

    def onAddCmdToApp(self):
        if not hasattr (self.app(), 'setSelectionLevel'):
            self.app().lazyLoad('selectionCmds', commands=['setSelectionLevel'], package='PmvApp')

    def checkArguments(self, molSel, negate=False, xor=False, intersect=False, only=False):
        """
        Required Arguments:\n
            molSel --- can be a string, a molecule or molecule set or selection or selectionSet
        Optional Arguments:\n
            negate --- if True, nodes are removed from current selection \n
            intersect --- if True, intersects nolSel with the currently active selection
            xor --- if True, selected elements of specified nodes will be \n
                    deselected and unselected will be selected
        """
        assert negate in [1,0,True, False]
        assert xor in [1,0,True, False]
        assert intersect in [1,0,True, False]
        assert only in [1,0,True, False]
        kw = {}
        kw['negate']= negate
        kw['xor']= xor
        kw['intersect']= intersect
        kw['only']= only
        return (molSel,), kw 

    def expandArg0(self, obj):
        # method to be overridden by subclass
        return [self.app().selectionSetFromObject(obj)]
    

    def doit(self, inSelectionSet, negate=False, xor=False, intersect=False, only=False):
        """add/remove nodes to the current selection.
        """
        activeSelectionSet = self.app().activeSelection

        if only:
            self.app().pushUndoCmd( self.app().select, (activeSelectionSet.copy(),), {'only':True})
            setOn = inSelectionSet - activeSelectionSet
            setOff = activeSelectionSet - inSelectionSet
            activeSelectionSet.replace(inSelectionSet)
        else:
            setOn = []
            setOff = []
            for inSelection in inSelectionSet:
                curSel = activeSelectionSet.selectionFor(inSelection)
                if curSel is not None:
                    ind = activeSelectionSet.index(curSel)

                if negate: # nodes will be unselected
                    if curSel is not None:
                        setOff = curSel
                        curSel = curSel - inSelection
                        if len(curSel):
                            activeSelectionSet[ind] = curSel
                            activeSelectionSet._agDict[curSel.getAtomGroup()] = curSel
                        else:
                            activeSelectionSet.remove(activeSelectionSet[ind])
                            del  activeSelectionSet._agDict[curSel.getAtomGroup()]
                        self.app().pushUndoCmd( self.app().select, (inSelection,), {'negate':False})

                ## elif xor: # selected elements of node will be deselected and
                ##     if curSel is not None:
                ##         setOff = activeSelectionSet & curSel
                ##         setOn = curSel
                ##         exor = (curSel | inSelection) - (curSel & inSelection)
                ##         ind = activeSelectionSet.index(curSel)
                ##         curSel = curSel ^ inSelection
                ##         activeSelectionSet[ind] = curSel
                ##         activeSelectionSet._agDict[curSel.getAtomGroup()] = curSel
                ##         self.app().pushUndoCmd( self.app().select, (exor,), {})

                elif intersect:
                    if curSel is not None:
                        union = curSel | inSelection
                        setOff = inSelection ^ curSel
                        setOn = inSelection & curSel
                        curSel = curSel & inSelection
                        if len(curSel):
                            activeSelectionSet[ind] = curSel
                            activeSelectionSet._agDict[curSel.getAtomGroup()] = curSel
                        else:
                            activeSelectionSet.remove(activeSelectionSet[ind])
                            del  activeSelectionSet._agDict[curSel.getAtomGroup()]
                        self.app().pushUndoCmd( self.app().select, (union,), {})

                else: # or
                    if curSel is not None:
                        setOn = curSel - inSelection
                    else:
                        setOn = inSelection

                    activeSelectionSet.append(inSelection)
                    self.app().pushUndoCmd( self.app().select, (inSelection,), {'negate':True})

        event = SelectionEvent(object=None, setOn=setOn, setOff=setOff)
        self.app().eventHandler.dispatchEvent(event)


    def startICOM(self):
        #print "in select.startICOM: mv.selectionLevel=", self.app().selectionLevel
        #print "mv.selection=", self.app().selection
        #print "in select.startICOM: mv.ICmdCaller.level.value=", self.app().ICmdCaller.level.value
        self.app().setSelectionLevel( self.app().ICmdCaller.level.value)

    
## !! The following code should go into the gui part of the command !!
    ## def updateSelectionFeedback(self):
    ##     self.updateInfoBar()
    ##     self.updateSelectionIcons()
        

    ## def updateInfoBar(self):
    ##     if self.app().hasGui:
    ##         if self.app().selectionLevel:
    ##             #msg = '\t Current selection: %d %s(s)' % ( len(self.selection),
    ##             #                          self.app().ICmdCaller.level.value.__name__)
    ##             msg = '%d %s(s)' % ( len(self.app().selection),
    ##                                  self.app().selectionLevel.__name__ )
    ##         else:
    ##             msg = '0 %s(s)'% self.app().selectionLevel__name__
    ##         self.app().GUI.pickLabel.configure( text=msg )

    def clear(self):
        setOff = self.app().activeSelection.get().copy()
        app = self.app()
        app.activeSelection.clear()
        event = SelectionEvent(app.activeSelection, new=[], old=setOff, setOn=[], setOff=setOff)
        app.eventHandler.dispatchEvent(event)

class MVExpandSelection(MVCommand):
    """Expands selection within a specified distance around current selection\n
    Package : PmvApp \n
    Module  : selectionCmds \n
    Class   : MVExpandSelection \n
    Command : expandSelection \n
    Synopsis:\n   
         None <- expandSelection(selection, distance, molFrag) \n
    Arguments:\n      
     """

    def __init__(self):
        MVCommand.__init__(self)
        self.iterateOverFirstArgument = False

    def checkArguments(self, selection, distance, molFrag):
        """selection - current selection , to be expanded within specified distance;
           molFrag - molecular fragment. Selection expands only to atoms that belong
                     to this fragment.
        """
        assert isinstance(distance, (int, float))
        assert distance > 0
        assert isinstance(selection, SelectionSet)
        assert len(selection)
        assert isinstance(molFrag, SelectionSet)
        #return (centerList, distance, selection)
        return (selection, distance, molFrag), {}

    def getAtoms(self, centerList, distance, molFrag):
        selSet = SelectionSet([])
        from prody.measure.contacts import findNeighbors
        #import pdb; pdb.set_trace()
        for sele in molFrag:
            # get a list of neighbors that are within *distance* of each other and the
            # distance between them
            res = numpy.array(findNeighbors(centerList, distance, sele))
            if not len(res):
                continue
            selinds = sele.getIndices()
            atinds = selinds[numpy.unique(res[:, 1]).astype("i")]
            selStr = "index  %s" % " ".join([str(item) for item in atinds])
            mol = sele._ag.getMolecule()
            selSet.append(mol.select(selStr))
        return selSet

    def doit(self, selection, distance, molFrag):
        """centerList, distance, selList"""
        #import pdb; pdb.set_trace()
        centers = selection[0].getCoords()
        if len(selection) > 1:
            import numpy
            for sele in selection[1:]:
                centers = numpy.concatenate((centers, sele.getCoords() ))
        selSet = self.getAtoms(centers, distance, molFrag)
        if len(selSet)>0:
            selection.extend(selSet)
            self.app().clearSelection()
            self.app().select(selection)    

            
class MVSelectAround(MVExpandSelection):
    """Select around  de-selects the current selection and selects the atoms
    within the specified distance of the current selection\n
    Package : PmvApp \n
    Module  : selectionCmds \n
    Class   : MVSelectAround \n
    Command : selectAround \n
    Synopsis:\n   
         None <- selectAround(selection, distance, molFrag) \n
    Arguments:\n      
     """
    def doit(self, selection, distance, molFrag):
        #import pdb; pdb.set_trace()
        centers = selection[0].getCoords()
        if len(selection) > 1:
            import numpy
            for sele in selection[1:]:
                centers = numpy.concatenate((centers, sele.getCoords() ))
        selSet = self.getAtoms(centers, distance, molFrag)
        if len(selSet)>0:
            self.app().clearSelection()
            self.app().select(selSet-selection)
            #self.app().select(selection, negate=True)

class MVInvertSelection(MVCommand):
    """Inverts current selection within all molecules or moleuels participating
    to the selection if invertLevel is specified or a within a user specified
    set is subset is specified. \n
    Package : PmvApp \n
    Module  : selectionCmds \n
    Class   : MVInvertSelection \n
    Command : invertSelection \n
    Synopsis:\n   
         None <- invertSelection(invertLevel=None, subset=None) \n
         Arguments:\n      
         invertLevel --- 'all' or 'molecule'  \n
         subset      --- set of atoms, residues, chains or molecules\n
     """

    def __init__(self):
        MVCommand.__init__(self)

    def checkArguments(self, selection, subset=None):
        """selection - current selection 
        subset      --- None or SelectionSet; if it is not None the inversion is made
                        in the given subset
        """
        assert isinstance(selection, SelectionSet)
        assert len(selection)
        if subset is not None:
            assert isinstance(subset, SelectionSet)
            assert len(subset)
        return (selection, subset), {}

    def doit(self, selection, subset=None):
        old=selection.copy()
        if subset is not None:
            for sele in subset:
                allobjects = SelectionSet()
                mol = sele._ag.getMolecule()
                allobjects.extend(mol.select())
        else:
            mol = selection._ag.getMolecule()
            allobjects = mol.select()

        self.app().deselect(old & allobjects, setupUndo=False)
        self.app().select(allobjects - old, setupUndo=False )


class MVAddSelectCommand(MVSelectCommand):
    """This Command adds to the current selection \n
    Package : PmvApp \n
    Module  : selectionCmds \n
    Class   : MVAddSelectCommand \n
    Command : addselect \n
    Synopsis:\n 
        currentSelection <--- addselect(nodes, klass=None) \n
        Required Arguments:\n       
        nodes --- TreeNodeSet holding the current selection \n
        klass --- if specified nodes are converted to a set of that type\n
                  and selection level is set to klass.
    """
    def checkArguments(self, nodes, klass=None):
        if isinstance (nodes, str):
           self.nodeLogString = "'"+nodes+"'"
        nodes = self.app().expandNodes(nodes)
        assert klass in [None, Atom, Residue, Chain, Molecule]
        kw = {}
        kw ['klass'] = klass
        kw['negate'] = 0
        return (nodes,), kw 

        
class MVXorSelectCommand(MVSelectCommand):
    """This Command xors the current selection \n
    Package : PmvApp \n
    Module  : selectionCmds \n
    Class   : MVXorSelectCommand \n
    Command : xorselect \n
    Synopsis:\n 
        currentSelection <--- xorselect(nodes, klass=None) \n
    Required Arguments:\n       
        nodes --- TreeNodeSet holding the current selection \n
        klass --- if specified, nodes are converted to a set of that type \n
                  and selection level is set to klass.
    """
    def checkArguments(self, nodes, klass=None, **kw):
        """ nodes --- TreeNodeSet holding the current selection \n
            klass --- if specified, nodes are converted to a set of that type\n
                      and selection level is set to klass.
        """
        assert klass in [None, Atom, Residue, Chain, Molecule]
        kw = {}
        kw['negate']= 0
        kw['only']= 1
        kw['klass'] = klass
        if len(self.app().activeSelection):
            current = self.app().activeSelection.findType(nodes.setClass)
            new_sel = current ^ nodes
        else:
            new_sel = nodes
        return (new_sel), kw

class MVDeSelectCommand(MVSelectCommand):
    """This Command deselects the current selection \n
    Package : PmvApp \n
    Module  : selectionCmds \n
    Class   : MVDeSelectCommand \n
    Command : deselect \n
    Synopsis:\n 
        currentSelection <--- deselect(nodes, klass=None) \n
    Required Arguments:\n       
        nodes --- TreeNodeSet holding the current selection \n
        klass --- if specified nodes are converted to a set of that type, \n
                  and selection level is set to klass.
    """
    def checkArguments(self, nodes, only=False, negate=False):
        """ nodes --- TreeNodeSet holding the current selection \n
        klass --- if specified nodes are converted to a set of that type otherwise
                  and selection level is set to Klass.
        """
        kw = {}
        kw['negate']= 1
        kw['only'] = only
        return (nodes,), kw 

class MVClearSelection(MVCommand):
    """ MVClearSelect implements method to clear the current selection. \n
    Package : PmvApp \n
    Module  : selectionCmds \n
    Class   : MVClearSelection \n
    Command : clearSelection \n
    Synopsis:\n
    None <- clearSelection(**kw)
   """

    def __init__(self):
       MVCommand.__init__(self)
       self.iterateOverFirstArgument = False
        
    def onAddCmdToApp(self):
        if not hasattr (self.app(), 'select'):
            self.app().lazyLoad('selectionCmds', commands=['select'], package='PmvApp')

    def doit(self):
        app = self.app()
        selection = app.activeSelection
        if len(selection)==0: return

        oldsel= selection.copy()
        self.app().pushUndoCmd( self.app().select, (selection.copy(),), {})
        self.app().activeSelection.clear()

        if self.createEvents:
            event = SelectionEvent(object=None, setOn=[], setOff=oldsel)
            app.eventHandler.dispatchEvent(event)
        
    def checkArguments(self):
        """None <- clearSelection(**kw)"""
        return (), {}



##########################################################################
##########################################################################
#
#  Support for static sets
#
##########################################################################
##########################################################################


sets__ = {} # name: (set,comments)

class MVSaveSetCommand(MVCommand):

    """Save a selection under a user specified name \n
    Package : PmvApp \n
    Module  : selectionCmds \n
    Class   : MVSaveSetCommand \n
    Command : saveSet \n
    Synopsis:\n
         None <- saveSet(nodes, name, comments='No description') \n
    Required Arguments:\n       
         nodes --- TreeNodeSet holding the current selection \n
         name --- name under which the selected set will be saved \n
    Optional Arguments:\n           
         comments --- description of the saved set, default is 'No description' \n
         addToDashboard --whether to add set name to dashboard historylist, default is do not add
    """

    def onAddCmdToViewer(self):
        from PmvApp.moleculeViewer import DeleteAtomsEvent
        self.app().eventHandler.registerListener(DeleteAtomsEvent,
                                 self.handleDeleteAtomsEvents)

    def handleDeleteAtomsEvents(self, event):
        atoms = event.objects
        for setName, _set in self.app().sets.items():
            if isinstance(_set[0], Atom):
                newSet = _set - atoms
                if len(newSet):
                    self.app().sets[setName] = _set - atoms
                else:
                    del self.app().sets[setName]

                    
    def doit(self, nodes, name, comments='No description', addToDashboard=False):
        if nodes.stringRepr:
            nodes.setStringRepr(nodes.extendStringRepr(nodes.stringRepr))        
        #sets__[name] = (nodes.__class__( nodes.data), comments)
        nodes.comments = comments
        self.app().sets.add(name, nodes)
        if addToDashboard:
            n = nodes[0].isBelow(Molecule)
            kstr = n*':' + name
            if  self.app().commands.has_key('dashboard'):
                self.app().dashboard.tree.selectorEntry.insert('end', kstr)
        # if vision is started add the set to the PmvApp library
        if hasattr (self.app(), "visionAPI") and self.app().visionAPI:
            fullname = nodes.full_name()
            self.app().visionAPI.add(nodes, name, kw={
                'set':nodes,
                'selString': fullname,
                'constrkw':{
                    'set':'masterNet.editor.vf.sets["%s"]'%(name),
                    'selString': 'selString',
                    }
                }
                                  )
        return name


    def checkArguments(self, nodes, name, comments='No description'):
        """ nodes --- TreeNodeSet holding the current selection \n
            name --- name of the saved set  \n
            comments --- description of the saved set, default is 'No description'
        """
        if isinstance(nodes, str):
            self.nodeLogString = "'" + nodes +"'"
        nodes = self.app().expandNodes(nodes)
        assert nodes
        assert isinstance (name , str)
        if name in self.app().sets.keys():
            newname = name + '_' + str(len(self.app().sets.keys()))
            self.app().warningMsg('set name %s already used\nrenaming it %s' % (
                name, newname ) )
            name = newname
        assert isinstance (comments, str)
        kw = {'comments': comments}
        return (nodes, name), kw 

        


class MVCreateSetIfNeeded(MVCommand):

    """create a set, but only if it does not already exist \n
    Package : PmvApp \n
    Module  : selectionCmds \n
    Class   : MVCreateSetIfNeeded \n
    Command : createSetIfNeeded \n\
    Synopsis:\n       
         None <- createSetIfNeeded(nodes, name, comments='No description') \n
    Required Arguments:\n        
        nodes --- TreeNodeSet holding the current selection \n
        name --- name under which the selected set will be saved \n
    Optional Argumnets:\n    
        comments --- description of the saved set, default is 'No description'
    """

    def doit(self, nodes, name, comments='No description'):
        if name in self.app().sets.keys():
            return
        #sets__[name] = (nodes.__class__( nodes.data), comments)
        #sets__[name] = (nodes, comments)
        nodes.comments = comments
        self.app().sets.add(name, nodes) 
        if hasattr (self.app(), 'visionAPI') and self.app().visionAPI:
            fullname = nodes.full_name()
            self.app().visionAPI.add(nodes, name, kw={
                'set':nodes,
                'selString': fullname,
                'name': name,
                'constrkw':{
                    'set':'masterNet.editor.vf.select("%s")'%fullname,
                    'selString': fullname,
                    'name':name}
                }
                                  )
        return self.app().sets[name]


    def checkArguments(self, nodes, name, comments='No description', **kw):
        """ None <- createSetIfNeeded(nodes, name, comments='No description',**kw)
        \nnodes --- TreeNodeSet holding the current selection
        \nname --- name under which the selected set will be saved
        \ncomments --- description of the saved set, default is 'No description'
        """
        
        assert isinstance (name, str)
        if isinstnce (nodes, str):
            self.nodeLogString = "'" + nodes +"'"
        nodes = self.app().expandNodes(nodes)
        assert nodes
        kw = {'comments':comments}
        return (nodes, name), kw 
        


class MVSelectSetCommand(MVCommand):
    """This Command is used to select the saved set. \n
    Package : PmvApp \n
    Module  : selectionCmds \n
    Class   : MVSelectSetCommand \n
    Command : selectSet \n
    Synopsis:\n
    None <- selectSet(setNames, only=False, negate=False) \n
    Required Arguments:\n
        setNameStr---name \n
    Optional Argumnets:\n
        only=negate=0:selected set is added to the selection \n
        only=0,negate=1: selected set is removed from the selection \n
        only=1,negate=0: selected set is the only thing selected"""

    def onRemoveObjectFromViewer(self, mol):
        if not self.app().sets.keys(): return
        for k in self.app().sets.keys():
            nodes = self.app().sets[k]
            if mol in nodes.top.uniq():
                if len(nodes.top.uniq())==1:
                    del self.app().sets[k]
                else:
                    molNodes = nodes.get(lambda x, mol=mol: x.top==mol)
                    newSet = nodes - molNodes
                    newComments = self.app().sets[k].comments + ' the elements of this set belonging to %s have been deleted with the molecule'%mol.name
                    del self.app().sets[k]
                    newSet.comments = newComments
                    self.app().sets.add(k, newSet)
                    del molNodes
        del nodes


    def undoCmdBefore(self, name, negate=False, only=False):
        # we overwrite undoCmdBefore enven though this command implements
        # the negate kw because the only kw would not be handled properly
        
        # create command to select the current selection
        # WARNING we cannot use getSelection here since it would
        # return everything if nothing is selected
        select = self.app().select
        if len(self.app().activeSelection)==0:
            return ( [(self.app().clearSelection, (), {})], self.app().clearSelection.name)
        else:
            return ( [(self.app().select, (self.app().activeSelection.atoms.copy(),), {})],
                              self.app().select.name )


    def doit(self, setNames, negate=False, only=False ):
        """None <- selectSet(setNames, only=False, negate=False) \n
        setNames can be a string or a list of strings \n
        only=negate=0:selected set is added to the selection \n
        only=0, negate=1: selected set is removed from the selection \n
        only=1, negate=0: selected set is the only thing selected"""

       
            
        if only:
            self.app().clearSelection()
        
        nodes = []
        set_names = self.app().sets.keys()
        for name in setNames:
            if name not in set_names:
                msgStr='\''+name+'\'' + ' not a previously defined set!'
                self.app().warningMsg(msgStr)
            else:
                #what if nodes and specified sets are not the same level?
                nodes = nodes + self.app().sets[name]
            #try:
            #    nodes = nodes + sets__[name][0]
            #except KeyError:
            #    msgStr='\''+name+'\'' + ' not a previously defined set!'
            #    self.app().warningMsg(msgStr)
        #if called with non-existent set, don't try to select it
        if nodes:
            self.app().select( nodes, negate=negate)




    def checkArguments(self, setNames, negate=False, only=False):
        """None <- selectSet(setNames, only=0, negate=0) \n
        setNames can be a string or a list of strings \n
        only=negate=0:selected set is added to the selection \n
        only=0, negate=1: selected set is removed from the selection \n
        only=1, negate=0: selected set is the only thing selected"""
        
        if isinstance(setNames, str):
            setNames = [setNames]
        else:
            if not isinstance (setNames, (list, tuple)):
                raise RuntimeError("%s, setNames should be a string or a list of strings" % self.name)
            for name in setNames:
                assert isinstance(name, str)
        assert negate in [1, 0 , True , False]
        assert only in [1, 0 , True , False]
        kw = {}
        kw['negate'] = negate
        kw['only'] = only
        return (setNames,), kw



class MVSelectFromStringCommand(MVCommand):
    """ molStr,chainStr,residueStr,atomStr Select items by typed-in strings: one for each level. No entry corresponds to select everything at this level. Strings are parsed as would-be regular expressions, replacing * by .*.... \n
    Package : PmvApp \n
    Module  : selectionCmds \n
    Class   : MVSelectFromStringCommand \n
    Command : selectFromString \n
    Synopsis:\n
        None <- selectFromString(nodes=None, mols='', chains='', res='', atoms='', deselect=False, silent=False) \n
    Arguments:\n
       mols, chains, res and atoms are strings. *?-, are supported \n
       silent --- when true, do not warn user about errors \n
       negate --- when True, the command does a deselect operation
    """

    def __init__(self):
        MVCommand.__init__(self)


    def undoCmdBefore(self, *args, **kw):
        select = self.app().select
        oldselection = self.app().activeSelection.atoms.copy()
        cmds = ( [(self.app().clearSelection ,(), {})], self.name)
        if len(oldselection):
            cmds[0].append( (self.app().select ,(oldselection,), {}) )
        return cmds


    def onAddCmdToApp(self):
        if not hasattr (self.app(), 'setSelectionLevel'):
            self.app().lazyLoad('selectionCmds', commands=['setSelectionLevel'], package='PmvApp')


    def doit(self, mols='', chains='', res='', atoms='', negate=False,
             silent=False, xor=False, intersect=False,  nodes=None):
        from MolKit2.tree import TreeNode
        app = self.app()
        if nodes is None:
            nodes = app.Mols
        elif isinstance(nodes, TreeNode):
            nodes = nodes.setClass([nodes])
        userpref = app.userpref['String Matching Mode']['value']
        caseSensitive = True
        escapeCharacters=True
        if 'caseInsens' in userpref:
            caseSensitive = False
        if 'WithEscapedChars' in userpref:
            escapeCharacters=True
        #if userpref == 'caseInsensWithEscapedChars': 
        #    newitem = 'cIWEC'
        #elif userpref == 'caseInsensitive': 
        #    newitem = 'cI'
        #else: newitem = 'cS'

        selectionList = [mols, chains, res, atoms]
        #need to format this correctly
        from string import join
        if atoms!='':
            selectionString = join(selectionList, ':')
        elif res!='':
            selectionString = join(selectionList[:-1], ':')
        elif chains!='':
            selectionString = join(selectionList[:-2], ':')
        else:
            selectionString = mols
        selitem, selMsg = StringSelector().select(nodes.top.uniq(),
                            selectionString, app.sets, caseSensitive,
                            escapeCharacters)

        if not selitem:
            if not silent:
                msg = 'No new selection:\n'+ selMsg
                app.warningMsg(msg)
            return

        if len(selMsg):
            app.warningMsg(selMsg)

        if selitem and len(selitem) > 0:
            lev = selitem[0].__class__
            if lev == Protein: lev = Molecule
            if app.selectionLevel != lev:
                app.setSelectionLevel(selitem[0].__class__)
            #    if not silent:
            #        t = "Current selection level is %s\nThis set holds %s objects. Do \
    #you want to set the selection level to %s ?" % (app.selectionLevel,
    #     lev, lev)
    #                d = SimpleDialog(app.GUI.ROOT, text=t,
    #                                 buttons=["Yes","No"],
    #                                 default=0, title="change selection level")
    #                ok=d.go()
    #                if ok==0: #answer was yes
    #                    app.setSelectionLevel(selitem[0].__class__, busyIdle=0, log=0)
    #                    #app.setIcomLevel( selitem[0].__class__, log = 0,
    #                           #KlassSet = selitem[0].setClass)
    #                else:
    #                    return 'ERROR' # prevent logging
    #            else:
    #                app.setSelectionLevel(selitem[0].__class__, busyIdle=0, log=0)
            app.select( selitem , negate=negate, xor=xor,
                            intersect=intersect, setupUndo=False)
            #if not negate:
            #    app.select( selitem , negate=negate)
            #else:
            #    #can't just deselect because SFString inherits from select
            #    old = app.getSelection()[:]
            #    if len(old)==0: return
            #    else:
            #        app.clearSelection()
            #        if len(selitem)!=len(old): 
            #            app.select(old-selitem)
        return app.activeSelection.get()


    def checkArguments(self,  mols='', chains='', res='', atoms='',
                 negate=False, silent=True, xor=False, intersect=False, nodes=None):
        """mols, chains, res and atoms are strings. *?-, are supported \n
        silent: when true, do not warn user about errors \n
        negate: when True, the command does a deselect operation \n
        nodes : nodes that will expand to molecules to be used to select from
        """
        #print "in selectFromString.__call__ with xor=", xor
            
        assert isinstance(mols, str)
        assert isinstance(chains, str)
        assert isinstance(res, str)
        assert isinstance(atoms, str)
        kw = {}
        kw['mols'] = mols
        kw['chains'] = chains
        kw['res'] = res
        kw['atoms'] = atoms
        kw['negate'] = negate
        kw['silent'] = silent
        kw['xor'] = xor
        kw['intersect'] = intersect
        kw['nodes'] = nodes
        return (), kw 


class MVDirectSelectCommand(MVSelectFromStringCommand):
    """This Command allows you to directly select from moleculelist,chainlist,setslist \n
    Package : PmvApp \n
    Module  : selectionCmds \n
    Class   : MVDirectSelectCommand \n
    Command : directSelect \n
    Synopsis:\n
        None <- directSelect(nameStr, **kw) \n
    Required Arguments:\n
        nameStr: name of selection from the string list
    """


    def findItemFromName(self,itemList,someName):
        #should add error detection at end
        for i in itemList:
            if i.name==someName:
                return i


    def findItemFromID(self,itemList,someID):
        #keys are "1hvr:A"
        for i in itemList:
            newStr = i.parent.name+":"+i.id
            if newStr == someID:
                return i


    def doit(self, nameStr):
        chNameList =[]
        app = self.app()
        ch = app.Mols.findType(Chain)
        if ch:
            chParents = ch.parent.name
            chIds=ch.id
            for i in range(len(chIds)):
                newItem=chParents[i]+':'+chIds[i]
                chNameList.append(newItem)
        if nameStr in app.Mols.findType(Molecule).name:
            app.setSelectionLevel(Protein)
            newMol=self.findItemFromName(app.Mols, nameStr)
            app.select(newMol)
        elif nameStr in chNameList:
            app.setSelectionLevel(Chain)
            chains = app.Mols.findType(Chain)
            newChain=self.findItemFromID(chains, nameStr)
            app.select(newChain)
        elif nameStr in app.sets.keys():
            newSet = app.sets[nameStr]
            newNode0 = newSet[0]
            lev=newNode0.__class__
            app.setSelectionLevel(lev)
            app.select(newSet)
        else:
            msg=nameStr + " not a Molecule, Chain or Set name"
            app.warningMsg(msg)


    def checkArguments(self, nameStr):
        """None <- directSelect(nameStr)
        \nnameStr --- name of selection from the string list"""
        assert isinstance(nameStr, str)
        return (nameStr,), {}


class MVSelectSphericalRegion(MVCommand):
    """This Command selects nodes from a specified base within specified spherical region(s) \n
    Package : PmvApp \n
    Module  : selectionCmds \n
    Class   : MVSelectSphericalRegion \n
    Command : selectInSphere \n
    Synopsis:\n
        None <- selectInSphere(centerList, radius, selList=['all']) \n
    Required Arguments:\n       
        centerList --- specifies the centers of the selection spheres \n
        radius --- radius of the selection spheres \n
        selList --- specifies selection base (atoms to test), \n
                possible selList values are ['all'], a list of sets or a list
                of molecule names.
    """
    def onAddCmdToViewer(self):
        self.radius = 5.
        self.centerList = None
        self.atar = []

        from DejaVu2.Spheres import Spheres
        from DejaVu2.Geom import Geom
        
        self.masterGeom = Geom("selectInSphereGeoms", 
                               shape=(0,0), protected=True)

        self.selSph = Spheres(name='SelSphReg_selSph', 
                              materials=((1.,1.,0),), shape=(0,3), radii=5.0, 
                              quality=15, inheritMaterial=0, protected=True,
                              inheritFrontPolyMode=0, frontPolyMode='line')
        from opengltk.OpenGL import GL
        self.selSph.frontPolyMode = GL.GL_LINE
        #self.selSph.frontPolyMode = GL.GL_POINT
        self.selSph.Set(visible=1, tagModified=False)
        self.selSph.pickable = 0
        from DejaVu2.Points import CrossSet
        self.cenCross = CrossSet('SelSphReg_cenCross', 
                                 inheritMaterial=0, materials=((1.,0.2,0),),
                                 offset=0.5, lineWidth=2, protected=True)
        self.cenCross.Set(visible=0, tagModified=False)
        self.cenCross.pickable = 0
        miscGeom = None
        event1 = AddGeometryEvent(self.masterGeom, parent=miscGeom)
        self.app().eventHandler.dispatchEvent(event1)
        # The AppGUI should register a listener for AddGeometryEvent
        # that implements  a method to add a geometry to
        # the application:
        #
            ## miscGeom = self.app().GUI.miscGeom
            ## self.app().GUI.VIEWER.AddObject(self.masterGeom, parent=miscGeom)
            ## self.app().GUI.VIEWER.AddObject(self.selSph, redo=0,
            ##                     parent=self.masterGeom)
            ## self.app().GUI.VIEWER.AddObject(self.cenCross, redo=0,
            ##
        from AppFramework.App import AddGeometryEvent
        miscGeom = None
        event1 = AddGeometryEvent(self.masterGeom, parent=miscGeom)
        self.app().eventHandler.dispatchEvent(event1)
        event2 = AddGeometryEvent(self.selSph, redo=0,
                                  parent=self.masterGeom)
        self.app().eventHandler.dispatchEvent(event2)
        event3 = AddGeometryEvent(self.cenCross, redo=0,
                                  parent=self.masterGeom)
        self.app().eventHandler.dispatchEvent(event3)


    def getTransformedCoords(self, atom):
        if not atom.top.geomContainer:
            return atom.coords
        g = atom.top.geomContainer.geoms['master']
        coords = atom.coords
        pth = [coords[0], coords[1], coords[2], 1.0]
        c = numpy.dot(g.GetMatrix(g), pth)[:3]
        return c.astype('f')


    def undoCmdBefore(self,  *args, **kw):
        select = self.app().select
        oldselection = self.app().activeSelection.atoms.copy()
        return ([(self.app().select, (oldselection,),
                  {'only':True})], self.name)


    def checkArguments(self, centerList, radius, selList=['all']):
        """centerList --- specifies the centers of the selection spheres \n
        radius --- radius of the selection spheres \n
        selList --- specifies selection base (atoms to test)\n
                possible selList values are ['all'], a list of sets or a list\n
                of molecule names. 
        """
        assert isinstance(centerList, (list, tuple, numpy.ndarray))
        assert len(centerList)
        if len(centerList) == 3 and isinstance (centerList[0], (int, float)):
            centerList = [centerList]
        assert isinstance(radius, (int, float))
        assert isinstance(selList, (list, tuple))
        if selList == 'all':
            selList = ['all']
        assert isinstance(selList, (list, tuple))
        assert len(selList)
        return (centerList, radius), {'selList':selList}


    def getAtoms(self, centerList, radius, selList):
        #base_nodes is an AtomSet constructed from selList
        app = self.app()
        base_nodes = AtomSet([])
        if selList[0]=='all':
            base_nodes = app.allAtoms
        elif selList[0] in app.sets.keys():
            #elif selList[0] in sets__.keys():
            #need to get all atoms in all sets specified
            for item in selList:
                newnodes = app.sets[item].findType(Atom)
                #newnodes = sets__[item][0].findType(Atom)
                base_nodes = base_nodes+newnodes
        else:
            for item in selList:
                mol = app.Mols.NodesFromName(item)[0]
                newnodes = mol.findType(Atom)
                base_nodes = base_nodes+newnodes
        if len(base_nodes)==0:
            t = '1:no base for selection specified: selList=', selList
            app.warningMsg(t)
            return base_nodes
        atar = numpy.array(base_nodes.data)

        #use transformed coords:
        tcoords = []
        for at in base_nodes:
            tcoords.append(self.getTransformedCoords(at))
        coords = numpy.array(tcoords, 'f')

        ats = []
        for center in centerList:
            d = coords - numpy.array(center, 'f')
            d = d * d
            d = numpy.sum(d,1)
            atindex = numpy.nonzero(numpy.less(d, radius * radius)) [0]
            newats = numpy.take(atar, atindex, axis=0)
            if len(newats)>0: 
                ats = ats + newats.tolist()
        return AtomSet(ats)


    def doit(self, centerList, radius, selList=['all']):
        """centerList, radius, selList"""
        #this cmd logs centerList as a list of arrays of coords
        #this cannot be replayed without from numpy.oldnumeric import array
        #however, in the tests self.app() has no logAllFile
        #hence this ugliness:  -rh 8/05
        if hasattr(self.app(), 'logAllFile'):
            self.app().logAllFile.write("from numpy import array, float32\n")
        ats = self.getAtoms(centerList, radius, selList)

        if len(ats)>0:
            #atSet = AtomSet(list(ats))
            ###FIX ME! this assumes 4 level hierarchy
            lev = ats[0].__class__
            if lev==Protein: 
                lev = Molecule
            vflev = self.app().selectionLevel
            if vflev!=lev:
                self.app().setSelectionLevel(lev)
            self.app().clearSelection()
            self.app().select(ats)



    def drawSphere(self):
        #callback to update command's geoms 
        ##7/21: now self.centerList is a List of centers
        from AppFramework.App import RedrawEvent
        if self.centerList and self.radius:
            self.selSph.Set(vertices=self.centerList, radii=self.radius,
                            visible=1, tagModified=False)
            event = RedrawEvent(self.app())
            self.app().eventHandler.dispatchEvent(event)
            #self.app().GUI.VIEWER.Redraw()
        else:
            self.selSph.Set(visible=0, tagModified=False)
            self.cenCross.Set(visible=0, tagModified=False)
            event = RedrawEvent(self.app())
            self.app().eventHandler.dispatchEvent(event)
            #self.app().GUI.VIEWER.Redraw()




class SelectNoWaterHeteroAtomsCommand(MVSelectFromStringCommand, MVSelectCommand):
    """This class provides a command to select all hetero atoms that are not
    in a water molecule. \n
    Package:PmvApp \n
    Module :selectionCmds \n
    Class:selectNonWaterHeteroAtomsCommand \n
    Command:selectHeteroAtoms \n
    Synopsis:\n
        None <--- selectHetereoAtoms(nodes, negate=False, only=False, xor=False, intersect=False, **kw)\n
        this command will select hetero atoms that are not water combine\n
        them with the current selection. As a consequence, the current\n
        selection will be changed to an atom set. The selected hetero atoms\n
        are returned as an AtomSet object.\n
    Required Arguments:\n   
        nodes --- TreeNodeSet from which to select hetero atoms \n
        example:\n
      >>> SelectNoWaterHeteroAtomsCommand()\n
      >>> SelectNoWaterHeteroAtomsCommand(molecule)\n
    """
    def __init__(self):
        MVSelectFromStringCommand.__init__(self)
        self.flag = self.flag | self.objArgOnly


    def checkArguments(self, nodes=None, negate=False, only=False,
                        xor=False, intersect=False):
        """
        Required Arguments:\n
            nodes --- can be a string, a TreeNode or a treeNodeSet \n
        Optional Arguments:\n
            negate --- is 1 nodes are removed from current selection \n
            only  --- is  1 the current selection is replaced by the nodes \n
            xor --- when True hetero atmos are xor'ed with current selection \n
            intersect -- when True hetero atoms are intersected with current selection
        """
        if nodes is None: nodes = self.app().Mols
        if isinstance(nodes, str):
           self.nodeLogString = "'"+nodes+"'"
        nodes = self.app().expandNodes(nodes)
        assert nodes
        kw = {}
        kw['negate']= negate
        kw['only']= only
        kw['xor']= xor
        kw['intersect']= intersect
        return (nodes,), kw 


    def doit(self, nodes, negate=False, only=False, xor=False, intersect=False):
        """
        select/deselect hetero atoms
        """
        app = self.app()
        allAtoms = nodes.findType(Atom)

        het = AtomSet( [a for a in allAtoms if a.hetatm and
                        a.parent.type!='WAT' and a.parent.type!='HOH'] )

        if len(het)==0:
            return AtomSet( [])
            
        app.setSelectionLevel(Atom, callListener=0)
        if only:
            app.setSelection(het)
                
        curSelection = app.activeSelection
        app.setSelection(curSelection.findType(Atom))
        if negate:
            app.setSelection(app.activeSelection - het)
        elif xor:
            app.setSelection(app.activeSelection ^ het)
        elif intersect:
            app.setSelection(app.activeSelection & het)
        else:
            app.setSelection(app.activeSelection | het)

        # the gui part of the command should register a listener for
        # the following event:
        event = selectionFeedbackEvent(objects=self.app().activeSelection, highlightSelection=True)
        self.app().eventHandler.dispatchEvent(event)
        
        #self.updateSelectionFeedback()
        #self.highlightSelection()
        
        return het



commandClassFromName = {
    'select' : [MVSelectCommand,  None],
    'deselect' : [MVDeSelectCommand,  None],
    'clearSelection' : [MVClearSelection, None],
    'expandSelection' : [MVExpandSelection, None],
    'selectAround' : [MVSelectAround, None],
    'saveSet' : [MVSaveSetCommand, None],
    'createSetIfNeeded' : [MVCreateSetIfNeeded, None],
    'invertSelection' : [MVInvertSelection, None],
    'selectSet' : [MVSelectSetCommand, None],
    'selectFromString' : [MVSelectFromStringCommand,  None],
    'directSelect' : [MVDirectSelectCommand,  None],
    'selectInSphere' : [MVSelectSphericalRegion,  None],
    'selectHeteroAtoms' : [SelectNoWaterHeteroAtomsCommand,  None],
    'setSelectionLevel' : [MVSetSelectionLevel, None]
}



def initModule(viewer, gui=True):
    for cmdName, values in commandClassFromName.items():
        cmdClass, guiInstance = values
        viewer.addCommand(cmdClass(), cmdName, guiInstance)









    





    

