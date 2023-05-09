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
# Date: 2014 Authors: Michel Sanner
#
#    sanner@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Michel Sanner and TSRI
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/PmvApp/GUI/Qt/dashboard.py,v 1.25.4.2 2017/11/01 20:34:29 annao Exp $
#
# $Id: dashboard.py,v 1.25.4.2 2017/11/01 20:34:29 annao Exp $
#

## Selections:
##   when nothing is selected, selecting somethign add current selection to to the tree
##   only one selection can be active at any time.
##   selected items are added to the current selection
##   the active selection has a yellow background in the Tree widget
##   clicking on the active selection (i.e. with yellow background) make the active selection None
##   There is no range or complex combinations of selections (i.e. multiple selections can not be
##     made blue and operated on as one can do with molecular fragments and groups
##

import weakref, os, numpy
from PySide import QtCore, QtGui

from MolKit2.molecule import Atom, Molecule, MultiMolecule, MoleculeSet, Residue, Chain#, MoleculeGroup, TreeObject
from MolKit2.selection import Selection, SelectionSet
from PmvApp.group import Group

from mglutil.util.callback import CallbackFunction
from mglutil.util.packageFilePath import findFilePath
PMVICONPATH = findFilePath('Icons', 'PmvApp.GUI')

from PmvGUI import GridGroup
from PmvApp.Pmv import PmvSelection, AfterAddMoleculeEvent, EditAtomsEvent, RefreshDisplayEvent
from PmvApp.group import Group

class ResTreeWidgetItem(QtGui.QTreeWidgetItem):

    def __lt__(self, other):
        column = self.treeWidget().sortColumn()
        key1 = self.text(column)
        key2 = other.text(column)
        return int(key1[3:]) < int(key2[3:])

class Dashboard(QtGui.QTreeWidget):

    def __init__(self, pmvGUI, parent=None):

        self.pmvGUI = pmvGUI
        #self.objToTreeitem = {}
        #self.treeitemToObj = {}
        QtGui.QTreeWidget.__init__(self, parent)

        self.setColumnCount(1)
        self.setHeaderLabels(['objects', ])
        self.currentItemChanged.connect(self.onSetCurrentItem)
        self.itemExpanded.connect(self.onItemExpanded)
        self.itemDoubleClicked.connect(self.showHide)
        ## self.itemClicked.connect(self.onItemClick)
        ## self.itemActivated.connect(self.onItemActivated)
        self.setAlternatingRowColors(True)
        self.setAutoScroll(True)
        
        self.setAcceptDrops(True)
        self.setDragEnabled(True)
        #self.setDropIndicatorShown(True) #is default

        # this will move the actual node
        self.setDragDropMode(QtGui.QAbstractItemView.InternalMove)

        self.currentSeleItem = None # will be a TreeWidgetItem when a selection is active
        from PmvApp.selectionCmds import SelectionEvent
        self.pmvGUI.pmv.eventHandler.registerListener(
            SelectionEvent, self.selectionEventHandler,
            before=self.pmvGUI.viewer.selectionEventHandler)

        from PmvApp.deleteCmds import DeleteObjectEvent
        self.pmvGUI.pmv.eventHandler.registerListener(
            DeleteObjectEvent, self.deleteObjectEventHandler)

        from PmvApp.Pmv import Event, DeleteAtomsEvent, AddAtomsEvent
        self.pmvGUI.pmv.eventHandler.registerListener(
            DeleteAtomsEvent, self.deleteAtomEventHandler)
        self.pmvGUI.pmv.eventHandler.registerListener(
            AddAtomsEvent, self.addAtomEventHandler)

        #root = self.invisibleRootItem()
        #root.setFlags(root.flags() & ~QtCore.Qt.ItemIsDropEnabled)
        self.colors = {
            'molecules':QtGui.QColor(255, 255, 255),
            'currentSelection':QtGui.QColor(253, 253, 150),
            'namedSelections':QtGui.QColor(119, 158, 203),
            'groups':QtGui.QColor(255, 105, 97),
            'grids':QtGui.QColor(255, 179, 71),
            'molecules':QtGui.QColor(255, 255, 255),
            'currentSelection':QtGui.QColor(253, 253, 150),
            'namedSelections':QtGui.QColor(255, 255, 255),
            'groups':QtGui.QColor(255, 255, 255),
            'grids':QtGui.QColor(255, 255, 255),
            'white':QtGui.QColor(255, 255, 255),
            }
        self.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        #selModel = self.selectionModel()
        #.setSelectionFlag(QtGui.QAbstractItemView.ToggleCurrent)
        self.itemExpanded.connect(self.onItemExpanded)
        #selModel.selectionChanged.connect(self.onSelectionChange)
        self.itemSelectionChanged.connect(self.onItemSelectionChange)
        self.setStyleSheet("QTreeWidget { background: lightGray}")

        # setup timer for playing models and multi molecules
        self._playing = False # toggle model playing using space bar on multimol
        self._modelPlayingTimer = None

    def multiGoto(self, mol, idx, item):
        """go to a given index of a multiMolecule or MultiConformation molecule"""
        # check if we are within bounds
        idxOrig = mol.curMolIndex()
        maxi = mol.numMols()
        idx = min(idx, maxi-1)
        idx = max(0, idx)
        if idx==idxOrig: return
        pmv = self.pmvGUI.pmv
        
        if mol._multi=='molecules':
            # remember index in dashboard so we can put molecule back
            atIndex = self.indexOfTopLevelItem(item)
            wasSelected = item.isSelected()
            # replace mol by new molecule
            mol.gotoMolecule(idx)

            pmv.deleteMolecule(mol, setupUndo=False)
            pmv.addMolecule(mol)
            pmv.applyDefaultCommands(mol, "Molecule")
            event = AfterAddMoleculeEvent(molecule=mol, atIndex=atIndex, select=wasSelected)
            pmv.eventHandler.dispatchEvent(event)

            if mol._renderingBits & mol._CPK and not hasattr(pmv.displayCPK, 'loadCommand'):
                pmv.displayCPK(mol, setupUndo=False)
            if mol._renderingBits & mol._SB:
                pmv.displaySB(mol, setupUndo=False)
            if mol._renderingBits & mol._SURFACE:
                surfName = '%s_surface'%mol._basename
                pmv.computeMSMS.guiCall(mol, surfName=surfName)
                pmv.displayMSMS.guiCall(mol, surfNames=[surfName])
            if mol._renderingBits & mol._CARTOON:
                pmv.computeCartoon(mol)
                pmv.displayCartoon(mol)
            if mol._labelingBits & mol._ATOM_LABEL:
                pmv.labelAtoms(mol)
            if mol._labelingBits & mol._RESIDUE_LABEL:
                pmv.labelResidues(mol)
                    
            if pmv.gui():
                # avoid skipping structure when scrolling fast
                pmv.gui().viewer.OneRedraw()

        else:
            if isinstance(mol, Group):
                mols = [x for x in mol._children if
                        hasattr(x, '_multi') and x._multi=='conformations']
                txt  = '%s %d/%d'%(mol.name.split(" ")[0], idx+1, maxi)
            elif mol._multi=='conformations':
                mol.name = '%s %d/%d'%(mol._basename, idx+1, maxi)
                txt = mol.name
                mols = [mol]
            else:
                return
            item.setText(0, txt)
            mol._currentIndex = idx
            for mol in mols:
                mol._ag.setACSIndex(idx)
                mol.geomContainer.allCoords[:] = mol._ag.getCoords()
                mol.name = '%s %d/%d'%(mol._basename, idx+1, maxi)
                item = mol._treeItems.keys()[0]
                item.setText(0, mol.name)

                pmv.gui().viewer.suspendRedraw = True
                #import pdb; pdb.set_trace()
                geomC = mol.geomContainer
                if hasattr(mol, "_msmsData"):
                    for name in mol._msmsData['msms'].keys():
                        if mol._msmsData['msms'][name][idx] is None:                        
                            pmv.computeMSMS(mol, surfName=name,
                                            **mol._msmsData['params'][name])

                if hasattr(mol, "_cartoonData"):
                    cartoonData = mol._cartoonData[idx]
                    if cartoonData is None:
                        pmv.computeCartoon(mol)
                    else:
                        for chid, data in cartoonData.items():
                            geom, v, f, n, v2a, f4a = data
                            geom.Set(vertices=v, vnormals=n, faces=[], visible=0)

                event = RefreshDisplayEvent(molecule=mol)
                pmv.eventHandler.dispatchEvent(event)
                pmv.gui().viewer.selectionEventHandler()

                # FIXME handle label, SS surfaces etc ..
                pmv.gui().viewer.suspendRedraw = False
                pmv.gui().viewer.OneRedraw()
            
    def wheelEvent(self, event):
        #print 'WHEEL', event.x(), event.y()
        item = self.itemAt(event.x(), event.y())
        ## if hasattr(item._pmvObject, 'gotoIndex'):
        ##     if item.isSelected():
        ##         items = self.selectedItems()
        ##     else:
        ##         items = [item]
        ##     for it in [it for it in items if hasattr(it._pmvObject,
        ##                                              'gotoIndex')]:
        ##         obj = item._pmvObject
        ##         if event.delta()>0:
        ##             self.multiGoto(obj, obj.curMolIndex()+1, item)
        ##         else:
        ##             self.multiGoto(obj, obj.curMolIndex()-1, item)

        if item is not None and hasattr(item._pmvObject, '_multi'):
            if item.isSelected():
                items = self.selectedItems()
            else:
                items = [item]

            for item in [it for it in items if hasattr(it._pmvObject, '_multi')]:
                mol = item._pmvObject
                if event.delta()>0:
                    self.multiGoto(mol, mol.curMolIndex()+1, item)
                else:
                    self.multiGoto(mol, mol.curMolIndex()-1, item)
        else:
            QtGui.QTreeWidget.wheelEvent(self, event)
            
    def mousePressEvent(self, e):
        ## only left mouse button press are handled as mouse press in the tree
        self._MouseBUtton = e.buttons()
        #print 'mousePress', int(e.buttons())
        if e.buttons()==QtCore.Qt.LeftButton:
            #item = self.itemAt(e.pos())
            #self.setCurrentItem(item)
            #item.setExpanded(not item.isExpanded())
            QtGui.QTreeWidget.mousePressEvent(self, e)

    def onItemSelectionChange(self):
        #
        # Here we enforce rules for extended selections
        # 1 - objects have to be of the same class
        # 2 - the roots can not be a selection and something else
        
        selectedItems = self.selectedItems()
        ## when first item is selected we remember class of selected object
        ## and the root of the item
        if len(selectedItems)==1:
            self._selectionClass = selectedItems[0]._pmvObject.__class__
            self._selectionRootObj = self.getRootItem(selectedItems[0])._pmvObject
            #print 'FIRST PICK', self._selectionClass, self._selectionRootObj.name
        else:
            if not hasattr(self, '_selectionClass'): # we got here with multiple selection before
                                                     # could set the class
                # deselect everythign and return
                selModel = self.selectionModel()
                for index in selModel.selection().indexes():
                    selModel.select(index, selModel.Deselect)
                return
            
        ## for all elements in the selection we verify that they are of the proper type
        ## and under the right root
        #print 'selected:'
        selModel = self.selectionModel()
        for index in selModel.selection().indexes():
            item = self.itemFromIndex(index)
            rootObj = self.getRootItem(item)._pmvObject

            # objects need to be of the same class
            if not isinstance(item._pmvObject, self._selectionClass):
                selModel.select(index, selModel.Deselect)
            elif isinstance(self._selectionRootObj, PmvSelection):
                # if objects are f the same class they also nee to have compatible roots
                if not isinstance(rootObj, PmvSelection):
                    selModel.select(index, selModel.Deselect)
            elif not isinstance(self._selectionRootObj, PmvSelection):
                if isinstance(rootObj, PmvSelection):
                    selModel.select(index, selModel.Deselect)
            #else:
            #    print 'keeping2', item.text(0)

    ##
    ## delete atoms
    def deleteObjectEventHandler(self, event):
        if event.objectType=='Molecule':
            mol = event.object
            self.invisibleRootItem().removeChild(mol._treeItems.keys()[0])

    def _removeItems(self, item, key):
        parent = item.parent()
        if parent is None: return
        parent.removeChild(item)
        del item._pmvObject._treeItems[key]
        if parent.childCount()==0:
            self. _removeItems(parent, key)

    def walkTreeDelete(self, node, deleted):
        sel = node._pmvObject._atomSubset
        #print node.text(0)
        #import pdb
        #pdb.set_trace()
        if len(sel & deleted):
            for n in xrange(node.childCount(), 0, -1):
                child = node.child(n-1)
                obj = child._pmvObject
                if isinstance(obj, Atom):
                    if len(obj.selSet[0] & deleted & deleted):
                        for k, item in obj._treeItems.items():
                            #import pdb
                            #pdb.set_trace()
                            #print 'removing', item().text(0)
                            self._removeItems(item(), k)
                            
                elif not hasattr(child, 'dummyChild'):
                    #import pdb
                    #pdb.set_trace()
                    self.walkTreeDelete(child, deleted)

    def deleteAtomEventHandler(self, event):
        deleted = event.deletedAtoms
        mol = deleted.getAtomGroup().getMolecule()
        for root in mol._treeItems.keys():
            if not hasattr(root, 'dummyChild'):
                self.walkTreeDelete(root, deleted)

    ##
    ## add atoms
    def addAtomEventHandler(self, event):
        #import pdb
        #pdb.set_trace()
        added = event.atoms
        mol = added.getAtomGroup().getMolecule()
        # collapse the molecule tree
        for root in mol._treeItems.keys():
            root._pmvObject._atomSubset |= added
            if not hasattr(root, 'dummyChild'):
                for n in xrange(root.childCount(), 0, -1):
                    child = root.child(n-1)
                    root.removeChild(child)
                root.dummyChild = QtGui.QTreeWidgetItem(root)
                root.setExpanded(False)
    ## def removeItem(self, items):
    ##     for k, item in items.items():
    ##         print 'remove %s under %s'%(item.text(0), self.getRootItem(item).text(0))

    ## def removeEmptyBranches(self, item, selection, root):
    ##     #print 'Entering', item.text(0), item, item.childCount()
    ##     #import pdb
    ##     #pdb.set_trace()
    ##     for i in range(item.childCount()):
    ##         child = item.child(i)
    ##         if hasattr(child, '_pmvObject'): # dummyChild does not have it
    ##             selectedAtoms = selection.inter(child._pmvObject.findType(Atom))
    ##             if len(selectedAtoms)==0:
    ##                 #try:
    ##                 del child._pmvObject._treeItems[root]
    ##                 #except KeyError:
    ##                 #    import pdb
    ##                 #    pdb.set_trace()
    ##                 #print 'deleting2', child.text(0), child, child._pmvObject.name
    ##                 child.parent().removeChild(child)
                        
    ##             else:
    ##                 self.removeEmptyBranches(child, selection, root)

    ## def onItemClick(self, item, column):
    ##     print 'mouse press intercepted', item, column
        
    ## def onItemActivated(self, item, column):
    ##     import pdb
    ##     pdb.set_trace()
    ##     print 'mouse activaed intercepted', item, column
        

    def selectionEventHandler(self, event):
        self._selectionEventHandler(event.object, event.setOn, event.setOff)

    def _selectionEventHandler(self, selection, setOn, setOff):
        #print 'SELECTIONEVENT', selection
        app = self.pmvGUI.pmv
        item = None
        #print 'DDDD', app.activeSelection, app.curSelection
        if app.curSelection.nbAtoms()==0:
            # the selection is empty and was shown in the tree => we remove it
            
            # if "Current selection" in the dashboard is not active (has white background),
            # self.pmvGUI.activeSelection is set to None. In this case we can not use
            # self.pmvGUI.activeSelection to remove the "Current selection" items
            # from the dashboard tree.
            
            #if hasattr(self.pmvGUI.activeSelection,'_treeItems') and \
            if hasattr(app.activeSelection,'_treeItems') and \
                   app.activeSelection is app.curSelection:
                #item = self.pmvGUI.activeSelection._treeItems.keys()[0]
                item = app.activeSelection._treeItems.keys()[0]
                self.invisibleRootItem().removeChild(item)
                #del self.pmvGUI.activeSelection._treeItems
                del app.activeSelection._treeItems
                self.pmvGUI.activeSelection = None # remember that no selection is active
                return
        
        ##
        ## selection in PMV is empty
        if app.activeSelection.nbAtoms()==0:
            if self.pmvGUI.activeSelection is None: # no selection is active in the dashboard
                return
            ##
            ##  named selection is active
            if app.activeSelection != app.curSelection:
                # we clear the sub tree to reflect empty selection
                # but keep this selection as active with top entry in dashboard
                item = app.activeSelection._treeItems.keys()[0]
                for i in range(item.childCount()):
                    item.removeChild(item.child(i))
                    item.dummyChild = QtGui.QTreeWidgetItem(item)
                    item.setExpanded(False)
                return
        ##
        ## selection is NOT empty
        else:
            ##
            ## not selection is active, i.e. the selection is in  app.curSelection
            if self.pmvGUI.activeSelection is None:
                # make curSelection the currently active selection
                app.activeSelection = app.curSelection
                self.pmvGUI.activeSelection = app.curSelection
                # add current selection to the tree if needed
                if hasattr(app.curSelection, '_treeItems'):
                    self.setCurrentSelection(app.curSelection._treeItems.keys()[0])
                else:
                    self.addObject(app.curSelection, None, 'Current Selection',
                                   color=self.colors['currentSelection'])
                return
            else:
                if hasattr(app.activeSelection, '_treeItems'):
                    item = app.activeSelection._treeItems.keys()[0]
            
        if item:
            if self.isItemExpanded(item):
                # delete subtree rooted at item
                for n in range(item.childCount()):
                    child = item.child(n)
                    item.removeChild(item.child(n))
                item.dummyChild = QtGui.QTreeWidgetItem(item)
                return

                # NOT sure if I could get this older veersion to work so for now
                # we collapse he selection
                
                ## # reset _treeItems for subtree
                
                ## # remove nodes of deselected atoms
                ## for atom in setOff:
                ##     if hasattr(atom, '_treeItem') and \
                ##            atom._treeItem.has_key(root):
                ##         #print 'FAGA', atom.name, atom._treeItems.keys()
                ##         item = atom._treeItems[root]
                ##         parent = item.parent()
                ##         parent.removeChild(item)
                ##         #print 'deleting1', item.text(0), item
                ##         del atom._treeItems[root]
                ##         while parent.childCount()==0:
                ##             grandParent = parent.parent()
                ##             print 'FUGU removing', parent.text(0)
                ##             grandParent.removeChild(parent)
                ##             del parent._pmvObject._treeItems[root]
                ##             parent = grandParent
                            
                ## if len(setOff):
                ##     self.removeEmptyBranches(item, self.pmvGUI.activeSelection, root)
                
                ## # add nodes of selected atoms
                ## parents = []
                ## for atom in setOn:
                ##     root = item
                ##     # find first ancestor that is shown in tree
                ##     if hasattr(atom, '_treeItems') and atom._treeItems.has_key(root):
                ##         #print 'FUGU12'
                ##         continue # already selected
                ##     obj = atom
                ##     while not hasattr(obj.parent, '_treeItems') or \
                ##               not obj.parent._treeItems.has_key(root):
                ##         if obj.parent is None:
                ##             break
                ##         else:
                ##             obj = obj.parent
                ##     #print 'FAGA', atom, obj, obj.parent, obj.parent._treeItems[root]._pmvObject.name
                ##     #print 'FAGAO', obj.name
                ##     if obj.parent is None:
                ##         parent = root
                ##     else:
                ##         parent = obj.parent._treeItems[root]
                ##     if parent.isExpanded():
                ##         newItem = self.addObject(obj, parent, obj.name.replace(' ', '_'))
                ##     parents.append(parent)

                ## # sort residues
                ## for parent in parents:
                ##     parent.sortChildren(0, QtCore.Qt.AscendingOrder)
                    
        
    def getObjectsForTreeItem(self, item):
        # gets all the objects in the subtree rooted at item
        # obj = self.treeitemToObj[item]
        obj = item._pmvObject
        root = self.getRootItem(item)

        if isinstance(obj, SelectionSet):
            # for selections return the intersection of the selection
            # with the atoms of the node corresponding to item
            return obj

        elif isinstance(root._pmvObject, SelectionSet):
            # for selections return the intersection of the selection
            # with the atoms of the node corresponding to item
            return root._pmvObject

        elif isinstance(obj, Group):
            # for groups return all molecules in subtree
            selections = SelectionSet()
            for n in range(item.childCount()):
                child = item.child(n)
                if isinstance(child._pmvObject, Group):
                    for s in self.getObjectsForTreeItem(child):
                        selections.append(s)
                elif isinstance(child._pmvObject, SelectionSet):
                    selections += child._pmvObject
                else:
                    selections.append( child._pmvObject.select() )
            return selections

        else:
            #print 'getObjectsForTreeItem', obj.__class__
            # for other nodes, i.e. Molecule, chains, Residues and Atoms
            #klass = obj.setClass
            #return klass([obj])
            return [obj.select()]
            #return obj

    def getIcon(self, obj):
        if isinstance(obj, Atom):
            icon = os.path.join(PMVICONPATH, "atom.png")
        elif isinstance(obj, Residue):
            icon = os.path.join(PMVICONPATH, "sidechain.png")
        elif isinstance(obj, Chain):
            icon = os.path.join(PMVICONPATH, "chain.png")
        elif isinstance(obj, Molecule):
            icon = os.path.join(PMVICONPATH, "molecule.png")
        elif isinstance(obj, Group):
            icon = os.path.join(PMVICONPATH, "group.png")
        elif isinstance(obj, SelectionSet):
            icon = os.path.join(PMVICONPATH, "selection.png")
        else:
            icon = os.path.join(PMVICONPATH, "python.gif")

        return QtGui.QIcon(icon)

    def getColor(self, obj):
        if isinstance(obj, Molecule):
            return self.colors['molecules']
        elif isinstance(obj, PmvSelection):
            if obj.name==u'Current Selection':
                return self.colors['currentSelection']
            else:
                return self.colors['namedSelections']
        elif isinstance(obj, Group):
            return self.colors['groups']
        else:
            return None
                
    def addObject(self, obj, parent, name, color=None, atIndex=None):
        # obj is a pmv object such as molecule or group
        # parent is a tree item
        # atIndex is the index in the children of parent
        assert hasattr(obj, 'nbChildren') and callable(obj.nbChildren)
        assert hasattr(obj, 'getChildrenAndNames') and callable(obj.getChildrenAndNames)
        
        if parent is None: # add a root object (not draggable)
            if atIndex is not None:
                root = item = QtGui.QTreeWidgetItem(None)
                self.insertTopLevelItem(atIndex, item)
            else:
                root = item = QtGui.QTreeWidgetItem(self)
            #item.setFlags(item.flags() & ~QtCore.Qt.ItemIsDragEnabled)
            #self.objToTreeitem[root] = {}
        else:
            if isinstance(obj, Residue):
                item = ResTreeWidgetItem(parent)
            else:
                item = QtGui.QTreeWidgetItem(parent)
            root = self.getRootItem(item)

            if isinstance(obj, Group):
                pass

            elif isinstance(obj, Molecule):

                # molecules inside selections cannot be dragged
                if isinstance(root._pmvObject, PmvSelection):
                    item.setFlags(item.flags() & ~QtCore.Qt.ItemIsDragEnabled)

                # disallow dropping on Proteins
                item.setFlags(item.flags() & ~QtCore.Qt.ItemIsDropEnabled)

            else: # chains, residues and atoms are not dragable or droppable
                item.setFlags(item.flags() & ~(QtCore.Qt.ItemIsDragEnabled |
                                               QtCore.Qt.ItemIsDropEnabled))

        if color:
            item.setBackground(0, color)

        if obj.nbChildren() > 0:
            # create a dummy child and add it so that the item is expandable
            # that way we can make lazy expansion
            item.dummyChild = QtGui.QTreeWidgetItem(item)

        # add the icon
        icon = self.getIcon(obj)
        if icon:
            item.setIcon(0, icon)

        #if hasattr(obj, 'alias') and obj.alias:
        #    name = '%s (%s)'%(obj.alias, name)

        item.setText(0, self.tr(name))
        if not hasattr(obj, '_treeItems'):
            obj._treeItems = {}
        obj._treeItems[root] = weakref.ref(item)
        #obj._treeItem = weakref.ref(item)

        #print '%%%%%%', obj.name, id(obj), obj.__class__, obj._treeItems
        item._pmvObject = obj
        #print 'ADDING1', item.text(0), item, name #, obj.name, obj.__class__, item, parent
        return item

    def getRootItem(self, item):
        root = item
        while root.parent() is not None:
            root = root.parent()
        return root                
        
    def onItemExpanded (self, item):
        #print 'Expanded', item, item._pmvObject
        #print 'Expand XXXXXX', item, obj, hasattr(item, 'dummyChild')

        root = self.getRootItem(item)
        if hasattr(item, 'dummyChild'):
            
            item.removeChild(item.dummyChild)
            del item.dummyChild
            obj = item._pmvObject
            if isinstance(root._pmvObject, SelectionSet):
                rootSelection = root._pmvObject
            else:
                rootSelection = None
            children, names = item._pmvObject.getChildrenAndNames(rootSelection)
            for child, name in zip(children, names):
                self.addObject(child, item, name)
                #self.addObject(child, item, name.replace(' ', '_'))
                
            ## #print 'Expanding', obj.name 
            ## #import pdb
            ## #pdb.set_trace()
            ## root = self.getRootItem(item)
            ## #if isinstance(root._pmvObject, PmvSelection):
            ## if isinstance(obj, PmvSelection):
            ##     if len(obj)==0:
            ##         item.setExpanded(False)
            ##         return

            ##     children = [c.getAtomGroup().getMolecule() for c in obj]
            ##     names = [c.name for c in children]

            ##     for c, name in zip(children, names):
            ##        self.addObject(c, item, name.replace(' ', '_'))

            ##     #for c in obj.children:
            ##     #    selectedChildren = selection.inter(c.findType(Atom))
            ##     #    if len(selectedChildren):
            ##     #        newItem = self.addObject(c, item, c.name.replace(' ', '_'))
            ## else:
            ##     item.removeChild(item.dummyChild)
            ##     del item.dummyChild
            ##     children = []

            ##     if isinstance(obj, Molecule):
            ##         if isinstance(root._pmvObject, PmvSelection):
            ##             for x in obj.getAtomGroup().iterChains():
            ##                 chain = Chain(x)
            ##                 if len(root._pmvObject & SelectionSet([chain.select()])):
            ##                     children.append(chain)
            ##         else:
            ##             children =[Chain(x) for x in obj.getAtomGroup().iterChains()]

            ##     elif isinstance(obj, Chain):
            ##         if isinstance(root._pmvObject, PmvSelection):
            ##             for x in obj._atomSubset.iterResidues():
            ##                 res = Residue(x)
            ##                 if len(root._pmvObject & SelectionSet([res.select()])):
            ##                     children.append(res)
            ##         else:
            ##             children = [Residue(x) for x in obj._atomSubset.iterResidues()]

            ##     elif isinstance(obj, Residue):
            ##         if isinstance(root._pmvObject, PmvSelection):
            ##             for x in obj._atomSubset.iterAtoms():
            ##                 atom = Atom(x)
            ##                 if len(root._pmvObject & SelectionSet([atom.select()])):
            ##                     children.append(atom)
            ##         else:
            ##             children = [Atom(x) for x in obj._atomSubset.iterAtoms()]
            ##     else:
            ##         raise ValueError("bad obbject for dashboard tree")
                
            ##     for child in children:
            ##        self.addObject(child, item, child.name.replace(' ', '_'))
        
    def startDrag(self, action):
        # save the object we are dragging
        self._draggedItems = self.selectedItems()#currentItem()
        QtGui.QTreeWidget.startDrag(self, action)

    def dragMoveEvent(self, event):
        ## to constrain dragging we have to monitor the motion and interactively decide
        ## whether the item we are on is suitable to be dropped on based on the type of
        ## the currently dragged object.
        ## in order to show the forbidden cursor icon we need to set the flags of all
        ## parent to non droppable.
        item = self.itemAt(event.pos()) # get item we are on now

        # Disallow dragging 
        if hasattr(item, '_pmvObject'): # the invisible root does not have it
            root = self.getRootItem(item)
            if isinstance(item._pmvObject, Molecule) or \
                   isinstance(root._pmvObject, PmvSelection):
                event.ignore()
                return

        if item is None:
            QtGui.QTreeWidget.dragMoveEvent(self, event)
            return

        obj = item._pmvObject # get object we are over
        items = []
        if isinstance(self._draggedItems[0]._pmvObject, Molecule):
            if isinstance(obj, (Molecule, GridGroup)): # make not droppable
                _obj = item
                for _obj in item, item.parent():
                    if _obj is not None:
                        flags = _obj.flags()
                        _obj.setFlags(_obj.flags() & ~QtCore.Qt.ItemIsDropEnabled)
                        items.append( (_obj, flags) )
                        #print 'resetting', _obj._pmvObject.name
                        _obj = _obj.parent()

        if len(items)==0:
            QtGui.QTreeWidget.dragMoveEvent(self, event)

        for item, flags in items:
            item.setFlags(flags)
        

    def dropEvent(self, e):
        #import pdb
        #pdb.set_trace()
        sourceItems = self._draggedItems
        destItem = self.itemAt(e.pos())

        for sourceItem in sourceItems:
            oldRoot = self.getRootItem(sourceItem)
            if destItem is not None:
                newRoot = self.getRootItem(destItem)
            try:
                #print 'dragged %s onto %s'%(sourceItem.text(0), destItem.text(0))
                if isinstance(destItem._pmvObject, Molecule):
                    #print 'REFUSE, return'
                    return
            except AttributeError:
                pass

            if destItem is None:
                self.pmvGUI.pmv.reparentObject(sourceItem._pmvObject, None)
            else:
                self.pmvGUI.pmv.reparentObject(sourceItem._pmvObject, destItem._pmvObject)
        #QtGui.QTreeWidget.dropEvent(self, e)



    def unsolicitedPick(self, pick):
        """treat an unsollicited picking event"""
        print 'FFFFF unsolicitedPick'

    def showHide(self, item, column, expandCollapse=True):
        if column==0:
            if isinstance (item._pmvObject, Molecule):
                mol = item._pmvObject
                visible = mol.geomContainer.masterGeom.visible
                self.pmvGUI.pmv.showMolecules(mol, show = not visible)
                if hasattr(mol, '_multi') and hasattr(self.pmvGUI.pmv, "displayHB"):
                    self.pmvGUI.pmv.displayHB.refreshDisplay(mol=mol)
            elif isinstance (item._pmvObject, Group):
                group = item._pmvObject
                for mol in group._children:
                    if hasattr(mol, '_multi'):
                        visible = mol.geomContainer.masterGeom.visible
                        self.pmvGUI.pmv.showMolecules(mol, show = not visible)
                        if hasattr(self.pmvGUI.pmv, "displayHB"):
                            self.pmvGUI.pmv.displayHB.refreshDisplay(mol=mol)
            else:
                visible = False
                if list(item.foreground(0).color().getRgb()[:3]) == [0,0,0]:
                    visible = True
                mols_selections = [[x.getAtomGroup().getMolecule(), x] for x in self.getObjectsForTreeItem(item)]
                for mol, sel in mols_selections:
                    # find all available geometries:
                    gc = mol.geomContainer
                    cmds = {}
                    for gname, ats in gc.atoms.items():
                        if gc.geoms.has_key(gname):
                            if gc.geomCmds.has_key(gname):
                                if visible:
                                    # need to undisplay geom
                                    cmd = gc.geomCmds[gname][1] #undisplayCmd
                                else:
                                    cmd = gc.geomCmds[gname][0] #displayCmd
                                if cmd: cmds[cmd] = sel
                    for cmd, sel in cmds.items():
                        cmd(sel)
            if visible: # we are hidding mol
                item.setForeground(0, QtGui.QColor(125, 125, 125)) # make text grey
            else:
                item.setForeground(0, QtGui.QColor(0,0,0)) # make text black
            # redo expand/collapse that double click triggered
            if expandCollapse:
                if item.isExpanded():
                    self.collapseItem(item)
                else:
                    self.expandItem(item)
                item.setSelected(False)
            
    def setCurrentSelection(self, item):
        app = self.pmvGUI.pmv
        if item:
            # make current GUI active selection background white
            if hasattr(self.pmvGUI.activeSelection, '_treeItems'):
                self.pmvGUI.activeSelection._treeItems.keys()[0].setBackground(0, self.colors['white'])
            # make background of new active selection yellow
            item.setBackground(0, self.colors['currentSelection'])

            # set the new selection to be active for the GUI
            self.pmvGUI.activeSelection = item._pmvObject

        else: # item is None means we want no selection to be active in the GUI
            if hasattr(self.pmvGUI.activeSelection, '_treeItems'):
                self.pmvGUI.activeSelection._treeItems.keys()[0].setBackground(0, self.colors['white'])
            self.pmvGUI.activeSelection = None

        self.pmvGUI.viewer.updateSelectionIcons()
        self.pmvGUI.viewer.highlightSelection()

    def onSetCurrentItem(self, current, previous):
        # called when an item is left click on in the dashboard
        #print 'clicked on', current.text(0), previous.text(0)
        if current and isinstance(current._pmvObject, SelectionSet):
            # clicked on a selection in the dashboard
            if current._pmvObject is self.pmvGUI.activeSelection:
                ## clicked on the the GUI's active selection selection ==>
                ##   hide the current selection and make pmv.curSelection the active selection
                self.pmvGUI.pmv.setActiveSelection(self.pmvGUI.pmv.curSelection, sendEvent=False)
                self.setCurrentSelection(None)
            else:
                self.pmvGUI.pmv.setActiveSelection(current._pmvObject, sendEvent=False)
                self.setCurrentSelection(current)

            # set the previous item to the the current tree item to avoid
            # blue background on active selection
            self.setCurrentItem(previous)
        else:
            msg = "toggle with keys: select l:lines, b:S&B, c:CPK, r:ribbon, m:surface"
            self.pmvGUI.statusBar().showMessage(self.tr(msg))
            
    def enterEvent(self, e):
        self.setFocus()
        self.pmvGUI.statusBar().showMessage(self.tr("Ready"))

    def leaveEvent(self, e):
        self.clearFocus()
        self.pmvGUI.statusBar().showMessage(self.tr("Ready"))

    def keyPressEvent(self, e):
        #print "key press", e.key(), e.text(), (e.modifiers() & QtCore.Qt.ControlModifier)== QtCore.Qt.ControlModifier

        # get the current item
        #curItem = self.currentItem()

        # get item under cursor
        #print 'KEY', self.mapFromGlobal(QtGui.QCursor.pos())
        wpos = self.mapFromGlobal(QtGui.QCursor.pos())
        item = self.itemAt(wpos.x()-2, wpos.y()-22) # not sure why the position has these offsets
        #import pdb; pdb.set_trace()
        if item is None:
            QtGui.QTreeWidget.keyPressEvent(self, e)
            return
        if item.isSelected():
            # get all items that are selected (i.e.) highlighted in the tree
            items = self.selectedItems()
        else:
            items = [item]

        app = self.pmvGUI.pmv
        selections = SelectionSet()
        for item in items:
            for sel in self.getObjectsForTreeItem(item):
                selections.append( sel )

        for sel in selections:
            mol = sel.getAtomGroup().getMolecule()
            gc = mol.geomContainer
            if hasattr(self, '_alt'):
                if e.text()=='c':
                    self._altC = True
                elif e.text()=='d':
                    self._altD = True
                del self._alt
            elif hasattr(self, '_altC'): # AltC color commands
                if e.text()=='a':
                    app.colorByAtomType(selections, geomsToColor=['all'])
                elif e.text()=='m':
                    app.colorByMolecules(sel, geomsToColor=['all'])
                elif e.text()=='p':
                    app.colorAtomsUsingDG(sel, geomsToColor=['all'])
                elif e.text()=='r':
                    app.colorRainbow(sel, geomsToColor=['all'])
                elif e.text()=='c':
                    app.colorRainbowByChain(sel, geomsToColor=['all'])
                del self._altC
            elif hasattr(self, '_altD'):
                if e.text()=='l':
                    print 'Sequence alt + D + L'
                del self._altD
            else:
                if e.text()=='s':
                    if app.isNotSelected(sel):            
                        app.select(sel, negate=0)
                    else:
                        app.select(sel, negate=1)
                elif e.text()=='l':
                    negate = gc.displayedAs(['noBonds', 'lines'],sel, 'fast')
                    app.undisplayLines.guiCall(sel) if negate else app.displayLines.guiCall(sel)
                elif e.text()=='b':
                    negate = gc.displayedAs(['sb'], sel, 'fast')
                    app.undisplaySB.guiCall(sel) if negate else app.displaySB.guiCall(sel)
                elif e.text()=='c':
                    negate = gc.displayedAs(['cpk'], sel, 'fast')
                    #print 'Display CPK', len(sel), negate
                    app.undisplayCPK.guiCall(sel) if negate else app.displayCPK.guiCall(sel)
                elif e.text()=='r':
                    #import pdb; pdb.set_trace()
                    g = gc.geoms.get('cartoon', None)
                    if g is None: # not cartoon geometry means no SS was computed so far
                        app.computeCartoon(sel)
                        negate = False
                    else: # a cartoon was Computed
                        if max([c.visible for c in g.children])>0:
                            negate = True
                        else:
                            negate = False
                    app.undisplayCartoon.guiCall(sel) if negate else app.displayCartoon.guiCall(sel)
                elif e.text()=='f':
                    self.pmvGUI.focusScene(sel)
                elif e.text()=='m':
                    ## this code is partly duplicated in PmvGUI. displayMSMSFor
                    ## shoudl eb re-concilde
                    for item in items:
                        if isinstance(item._pmvObject, Molecule):
                            if not hasattr(mol, '_msmsData'):
                                app.computeMSMS(mol)
                                app.displayMSMS(mol)
                            else:
                                for name in mol._msmsData['msms'].keys():
                                    # if not surface compute it and display it 
                                    if mol._msmsData['msms'][name][mol._ag.getACSIndex()] is None:
                                        app.computeMSMS(mol)
                                        app.displayMSMS(mol)
                                    else:
                                        if mol.geomContainer.geoms['msms_'+name].visible:
                                            app.undisplayMSMS(mol)
                                        else:
                                            app.displayMSMS(mol)
                        else:# isinstance(item._pmvObject, Selection):
                            #mol = item._pmvObject.getAtomGroup().getMolecule()
                            print 'NOT YET'
                            
                    #objItems = [ (x._pmvObject, x) for x in items ]
                    #atoms = SelectionSet([ x._pmvObject.select() for x in items ])
                    #import pdb; pdb.set_trace()
                    #self.pmvGUI.viewer.displayMSMSfor(objItems, atoms)
                elif (e.modifiers() & QtCore.Qt.AltModifier) == QtCore.Qt.AltModifier:
                    self._alt = True

                ## multi conformations
                elif e.key()==QtCore.Qt.Key_Left: # NEXT
                    if mol._multi:
                        self.multiGoto(mol, mol.curMolIndex()-1, mol._treeItems.keys()[0])
                elif e.key()==QtCore.Qt.Key_Right: # PREVIOUS
                    if mol._multi:
                        self.multiGoto(mol, mol.curMolIndex()+1, mol._treeItems.keys()[0])
                elif e.key()==QtCore.Qt.Key_Home: # FIRST
                    if mol._multi:
                        self.multiGoto(mol, 0, mol._treeItems.keys()[0])
                elif e.key()==QtCore.Qt.Key_End: # LAST
                    if mol._multi:
                        self.multiGoto(mol, mol.numMols()-1, mol._treeItems.keys()[0])
                elif e.text()==' ': # toggle play
                    multiMols = [it._pmvObject for it in items if isinstance(it._pmvObject, Molecule) and it._pmvObject._multi]
                    if len(multiMols) ==0 : return

                    self._playing = not self._playing
                    if self._playing: # start playing timer
                        self._modelPlayingTimer = QtCore.QTimer(self)
                        cb = CallbackFunction(self.playNextModel, multiMols)
                        self._modelPlayingTimer.timeout.connect(cb)
                        self._modelPlayingTimer.start(500)
                        return
                    else: # stop playing timer
                        self._modelPlayingTimer.stop()
                        return
                    
        ## elif e.text()=='r':
        ##     negate = gc.displayedAs(geoms, atoms, 'fast')
        ##     app.displaySecondaryStructure(atms, negate=negate)
        ## elif e.text()=='m':
        ##     negate = gc.displayedAs(['cpk'], atoms, 'fast')
        ##     app.displayCPK(atms, negate=negate)
    def playNextModel(self, molecules):
        allAtend = True
        for mol in molecules:
            idx = mol.curMolIndex()+1
            if idx < mol.numMols():
                allAtend = False
                self.multiGoto(mol, mol.curMolIndex()+1, mol._treeItems.keys()[0])
        if allAtend:
            self._modelPlayingTimer.stop()
            
    def getLevel(self, item):
        # return the level
        # levels: 0: atom
        #         1: residue
        #         2: chain
        #         3: molecule
        #         4: selection
        #         5: group

        if isinstance(item._pmvObject, Group):
            level = 5
        elif isinstance(item._pmvObject, Molecule):
            level = 4
        elif isinstance(item._pmvObject, SelectionSet):
            level = 3
        elif isinstance(item._pmvObject, Chain):
            level = 2
        elif isinstance(item._pmvObject, Residue):
            level = 1
        elif isinstance(item._pmvObject, Atom):
            level = 0
        return level
    
    def contextMenuEvent(self, event, item=None):
        # find the item on which we clicked
        if item is None:
            item = self.itemAt(event.pos())

        if item is None: return
        level = self.getLevel(item)
        # find the asscoiate Pmv object
        obj = item._pmvObject

        # if the item we clicked on is selected we want the operation to
        # occur for everything selected in the dashboard
        selectedItems = self.selectedItems()

        # check is clicked item is selected
        # for some reason item in selectedItems fails so we have to check manually
        itemInSel = False
        for it in selectedItems:
            if it is item:
                itemInSel = True
                break
        ## if the item is selected
        ## replace obj by a set of object for all selected items
        if itemInSel:
            items = selectedItems
        else:
            items = [item]

        objItems = [(it._pmvObject, it) for it in items]
        ## objItems = []
        ## for it in items:
        ##     if isinstance(it._pmvObject, SelectionSet):
        ##         objItems.append( (it._pmvObject, it) )
        ##     else:
        ##         objItems.append( (it._pmvObject.selSet, it) )

        return self.pmvGUI.viewer.buildMenuForObject(objItems, level, parent=self, parentMenu=None)

        ## if itemInSel:
        ##     obj = [it._pmvObject for it in selectedItems]
        ##     obj = obj[0].setClass(obj)
            
        ## else:
        ##     obj = obj.setClass([item._pmvObject])

        ## at this point obj is always a set of potentially one object
        ## print 'OBJECTS', [o.name for o in obj]

        ## # for nodes under a selection, right click should not makes us loose the selection
        ## # so we force the current item to be the root
        ## root = self.getRootItem(item)
        ## #if isinstance(root._pmvObject, PmvSelection):
        ## #    self.setCurrentItem(root)

        ## if isinstance(obj, Group):
        ##     self.pmvGUI.viewer.buildMenuForObject(obj, parent=self, parentMenu=None,
        ##                                    target=self.getObjectsForTreeItem(item),
        ##                                    rootItem=root)
        ## else:
        ##     self.pmvGUI.viewer.buildMenuForObject(obj, parent=self, parentMenu=None, rootItem=root)
            


    #def onItemPressed (self, item, column):
    #    print 'pressed', item, self.treeitemToObj.get(item, 'dummy')

    #def onExpandItem(self, item):
    #    print 'onExpandItem', item, self.treeitemToObj.get(item, 'dummy')
            

     ## def __init__(self, parent=None):

    ##     QtGui.QTreeWidget.__init__(self, parent)

    ##     self.pmv = None
    ##     self.rootRow = None
        
    ##     tree = self.tree = QtGui.QTreeWidget(parent)
    ##     self.headerLabels = ['geometry', 'info']
    ##     self.nbCol = len(self.headerLabels)
    ##     tree.setColumnCount(self.nbCol)
    ##     tree.setHeaderLabels(self.headerLabels)
            
    ##     mainLayout = QtGui.QVBoxLayout()   
    ##     mainLayout.addWidget(tree)

    ##     self.setLayout(mainLayout)
        
    ##     tree.show()

    ## def bind(self, pmv):
    ##     self.pmv = pmv
    ##     tree = self.tree
    ##     row = QtGui.QTreeWidgetItem(tree)
    ##     row.setText(0, 'All')
    ##     row.buttons = []
    ##     for i in range(1,self.nbCol):
    ##         b1 = QtGui.QCheckBox(tree)
    ##         self.connect(b1, QtCore.SIGNAL("stateChanged(int)"), self.box_cb)
    ##         b1.object = pmv.Molecules
    ##         b1.col = i
    ##         row.buttons.append(b1)
    ##         tree.setItemWidget(row, i, b1)
    ##     self.rootRow = row

    ##     for mol in pmv.Molecules:
    ##         self.addMolecule(mol)

        
    ## def addMolecule(self, mol):
    ##     def addRow(object, parent):
    ##         row = QtGui.QTreeWidgetItem(parent)
    ##         row.setText(0, object.name)
    ##         row.buttons = []
    ##         for i in range(1,self.nbCol):
    ##             b1 = QtGui.QCheckBox()
    ##             b1.object = object
    ##             b1.col = i
    ##             self.connect(b1, QtCore.SIGNAL("stateChanged(int)"), self.box_cb)
    ##             row.buttons.append(b1)
    ##             self.tree.setItemWidget(row, i, b1)
    ##         for child in object.children:
    ##             if isinstance(child, Atom):
    ##                 return
    ##             addRow(child, row)

    ##     addRow(mol, self.rootRow)


    ## def box_cb(self, i):
    ##     button = self.sender()
    ##     frag = button.object
    ##     mol = frag.top
    ##     geomC = mol.geomContainer
    ##     if button.col==3: # cpk
    ##         atmset = geomC.atoms['cpk']
    ##         if i:
    ##             atms = atmset + frag.findType(Atom)
    ##         else:
    ##             atms = atmset - frag.findType(Atom)

    ##         geom = geomC.geoms['cpk']
    ##         if len(atms)==0:
    ##             geom.Set(visible=0, tagModified=False)
    ##         else:
    ##             try:
    ##                 rads = atms.radius
    ##             except AttributeError:
    ##                 rads = mol.defaultRadii()
    ##             geom.Set(visible=1, vertices=atms.coords, radii=atms.radius)

    ##         geomC.atoms['cpk'] = atms
            
    ##         vi = geom.viewer
    ##         vi.Redraw()
    ##         vi.cameras[0].updateGL()
    ##     #print i, frag, self.headerLabels[button.col]
        

if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)

    from MolKit import Read
    from time import time
    t1 = time()
    mol = Read(sys.argv[1])[0]
    print 'read mol in ', time()-t1

    class Pmv:
        def __init__(self):
            self.Molecules = []

    pmv = Pmv()
    pmv.Molecules.append(mol)
    
    t1 = time()
    db = Dashboard(pmv)
    print 'built tree in ', time()-t1

