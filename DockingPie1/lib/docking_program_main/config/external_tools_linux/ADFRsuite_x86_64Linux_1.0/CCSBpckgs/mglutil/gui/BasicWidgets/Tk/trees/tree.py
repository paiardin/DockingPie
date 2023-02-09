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
##
## Author Michel Sanner  Copyright TSRI (C) 2006 
##
########################################################################
#
# $Header: /mnt/raid/services/cvs/python/packages/share1.5/mglutil/gui/BasicWidgets/Tk/trees/tree.py,v 1.59.6.1 2017/07/26 22:35:42 annao Exp $
#
# $Id: tree.py,v 1.59.6.1 2017/07/26 22:35:42 annao Exp $
#
"""Efficient TreeWidget based on Tkinter and Pmw

The speed comes from only drawing visibles parts of the tree.
"""

import os, sys, types, Pmw
from PIL import Image, ImageTk
from weakref import ref
from Pmw import ScrolledCanvas
from Tkinter import Tk, PhotoImage, Label, ALL, Menu, IntVar
from mglutil.util.packageFilePath import findFilePath
from mglutil.gui import widgetsOnBackWindowsCanGrabFocus
from mglutil.events import Event, EventHandler
from mglutil.util.callback import CallbackFunction

## TODO
## - add size subdirectory to IconManager
## - cannot select 2 ranges, might want to add drag selection

##
## BUGS
##

class IconsManager:
    """The IconsManager object simplifies creating PhotoImage icons from
images in the directory.  The directory can be set by specifying a path,
a path within a Python package, or a Pyton package.
Images are loaded upon requests and a reference to the icon is saved to
avoid its garbage collection."""


    def __init__(self, path=[], packageName=''):
        """Create an IconManager object.
obj <- IconManager(path='', packageName='')
If packageName is specified, the path is assumed to be relative to the location
of the package."""
        if packageName:
            mod = __import__(packageName)
            components = packageName.split('.')
            for comp in components[1:]:
                mod = getattr(mod, comp)
            packagePath = os.path.abspath(mod.__path__[0])
        else:
            packagePath=''
            
        self.directory = os.path.join( packagePath, os.path.join(*path))
        
        assert os.path.exists(self.directory)
        assert os.path.isdir(self.directory)
        self.icons = {}


    def get(self, iconName, master=None, directory=None):
        """ return an Tk PhotoImage for a given icon name and saves a reference
"""
        if directory is None:
            directory = self.directory
        icon = self.icons.get(iconName, None)
        if icon is None:
            filename = os.path.join(directory, iconName)
            image = Image.open(filename)
            icon = ImageTk.PhotoImage(master=master, image=image)
            self.icons[iconName] = icon
        return icon



class KeySelectable:
    """Adds the ability to use keystrokes to quickly select items in a list.
root has to be a widget supporting .bind .after
"""
    
    def __init__(self, KeyRootTk):
        self.KeyRootTk = KeyRootTk
        self.afterID = None
        self.matchString = ''
        self.lastMatchString = ''
	KeyRootTk.bind('<KeyPress>', self.key_cb)
	KeyRootTk.bind('<KeyRelease>', self.keyUp_cb)
        self.isControl = False
        self.isShift = False
        self.isAlt = False
        self.ctrlModCallback = None
        self.shiftModCallback = None
        self.altModCallback = None
        

    def timeOut(self, event=None):
        """resets self.matchCharIndex to 0, called after a short period of
        time if no new character has been typed"""
        #print 'timeout'
        self.lastMatchString = self.matchString
        self.matchString = ''
        self.matchCharIndex = 1
        self.afterID = None
        

    def keyUp_cb(self, event=None):
        if event.keysym=='Control_L' or event.keysym=='Control_R':
            self.isControl = False
        elif event.keysym=='Shift_L' or event.keysym=='Shift_R':
            self.isShift = False
        elif event.keysym=='Alt_L' or event.keysym=='Alt_R':
            self.isAlt = False

        
    def key_cb(self, event=None):
        # use key strokes to select entry in listbox
        # strokes placed within 500 miliseconds are concatenated
        #print self.matchCharIndex, '|', self.matchString, '|', event.keysym
        if event.keysym=='Control_L' or event.keysym=='Control_R':
            self.isControl = True
            return
        elif event.keysym=='Shift_L' or event.keysym=='Shift_R':
            self.isShift = True
            return
        elif event.keysym=='Alt_L' or event.keysym=='Alt_R':
            self.isAlt = True
            return

        if self.isControl:
            if self.ctrlModCallback:
                self.ctrlModCallback(event)
            return
        elif self.isShift:
            if self.shiftModCallback:
                self.shiftModCallback(event)
            return
        elif self.isAlt:
            if self.altModCallback:
                self.altModCallback(event)
            return
            
        if event.keysym=='Return':
            str = self.lastMatchString
        else:
            str = self.matchString + event.keysym
        #print str
        item = self.match(str)
        if item:
            self.selectItem(item)
            if self.afterID is not None:
                self.KeyRootTk.after_cancel(self.afterID)
            self.afterID = self.KeyRootTk.after(1000, self.timeOut)
            self.matchString = str

    # SUBCLASS THIS
    def match(self, name):
        """has to return None if no match or an object that matches"""
        print 'jumping to ', name
        return None


    def selectItem(self, item):
        """do what has to be done to show what matches the typed string"""
        print 'selecting item', item


class SelectEvent(Event):
    def __init__(self, object):
        Event.__init__(self)
        self.object = object

class DeselectEvent(Event):
    def __init__(self, object):
        Event.__init__(self)
        self.object = object

class ClearSelectionEvent(Event):
    def __init__(self, object):
        Event.__init__(self)
        self.object = object


class Tree(ScrolledCanvas, KeySelectable, EventHandler):
    """Tree widget for Tk and Pmw"""


    def __init__(self, master, root, iconsManager=None, selectionMode='single',
                 idleRedraw=True, nodeHeight=15, headerHeight=0,
                 balloon=None, **kw):
        """Tree( master, root, idleRedraw=True, kw = {}, **opts)
- Master can be a Tk frame in which the tree is displayed
- root has to be a Node object
- iconsManageer has to be an instance of an IconManager object or None.
If None is passed, an IconManager with the default icons directory will be
created.  If an IconManager is passed we expect to find ...
- selection mode cane be 'single' or 'multiple'. In multiple mode, the shift
modifier is use to select/deselect ranges (ranges can only be defined on
sibling nodes) and the Control modifier is used to add/remove to the current
selection.  Node are selected by click withthe left mouse button on the node's
label or icon.
- Set idleRedraw to False to disable background drawing.
- nodeHeight is the heigh of a node int he Tree in pixels
- kw can contain any keywork allowable for a Pmw.ScrolledCanvas
"""
        assert isinstance(root, Node)
        assert selectionMode in ['single', 'multiple']

        EventHandler.__init__(self)

        self.master = master
        
        self.selectionMode = selectionMode
        self.lastPickedNode = None
        self.balloon = balloon
        
        self.idleRedraw  = idleRedraw  # when set to true redraw operation
                                       # occur when CPU is idle
        
        if iconsManager is None:
            iconsManager = IconsManager(
                ['gui','BasicWidgets','Tk','trees','Icons'], 'mglutil')

        assert isinstance(iconsManager, IconsManager)
        self.iconsManager = iconsManager
        fs = self.fontSize = 10
        self.font = "Arial %d"%fs
        
        # cache the icons use to expand and collapse since we will
        # use them for each node
        mod = __import__('mglutil')
        packagePath = mod.__path__[0]
        path = ['gui','BasicWidgets','Tk','trees','Icons']
        directory = os.path.join( packagePath, os.path.join(*path))
        #self.collapsedIcon  = self.iconsManager.get("1rightarrow.png",
        #                                            directory=directory)
        #self.expandedIcon = self.iconsManager.get("1downarrow.png",
        #                                          directory=directory)
        self.collapsedIcon  = self.iconsManager.get("plus.png",
                                                    directory=directory)
        self.expandedIcon = self.iconsManager.get("minus.png",
                                                  directory=directory)
        self.iconHalfWidth = self.collapsedIcon.width()/2 
        self.selectedNodes = []
        self.selectionHistory = [] # used to undo selection operations
        self.nodeHeight = nodeHeight
        self.headerHeight = headerHeight
        
        self.firstVisibleNodeLineNum = 0 # height of the first visible node
                                # i.e. how many nodes drown above, 0 for root

        self.firstVisibleNode = root # first node visible in window
        
        self.displayedNodes = [] # list of nodes having a graphical represenation
        self.pending_draw = None
        self.suspendRedraw = False

        self.root = root
        root.currenty = 0
        self.maxNodeNum = 0
        
        root.tree = ref(self)

        root.childNumber = 0
        self.nbLines = 1 # height of the tree (i.e. how many lines in canvas)

        # build the GUI
        if not kw.has_key('vscrollmode'):
            kw['vscrollmode'] ='dynamic'
        if not kw.has_key('hscrollmode'):
            kw['hscrollmode'] ='dynamic'
        if not kw.has_key('hull_width'):
            kw['hull_width'] = 550
        if not kw.has_key('hull_height'):
            kw['hull_height'] = 100
        if not kw.has_key('background') and not kw.has_key('bg') and \
               not kw.has_key('canvas_bg'):
            kw['canvas_bg']='white'
        if not kw.has_key('usehullsize'):
            kw['usehullsize'] = 1
            
        self.nbLinesPerPage = (kw['hull_height']-self.headerHeight) / nodeHeight

        if not kw.has_key('yscrollincrement'):
            kw['canvas_yscrollincrement'] = nodeHeight

        kw['horizscrollbar_command'] = self.xview
        kw['vertscrollbar_command'] = self.yview
        kw['borderframe'] = 5
        kw['canvas_highlightthickness'] = 0
        
        #kw['canvasmargin'] = 10
        #kw['hscrollmode'] = 'none'
        
        apply( ScrolledCanvas.__init__, (self, master), kw)

        canvas = self.canvas = self.component('canvas')
        canvas.master.configure(borderwidth=0)
        
        KeySelectable.__init__(self, canvas)
        self.ctrlModCallback = self.handleControlKey

        canvas.bind("<Configure>", self.configure_cb)
        canvas.bind("<Key-Prior>", self.pageUp)
        canvas.bind("<Key-Next>", self.pageDown)
        canvas.bind("<Key-Up>", self.lineUp)
        canvas.bind("<Key-Down>", self.lineDown)
        if os.name == 'nt': #sys.platform == 'win32':
            canvas.bind("<MouseWheel>", self.lineUpDown)
        else:
            canvas.bind("<Button-4>", self.lineUp)
            canvas.bind("<Button-5>", self.lineDown)
        canvas.bind("<Enter>", self.enter_cb)

        self.isControl = False
        self.isShift = False
        self.isAlt = False

        if widgetsOnBackWindowsCanGrabFocus is False:
            lActiveWindow = canvas.focus_get()
            if    lActiveWindow is not None \
              and ( lActiveWindow.winfo_toplevel() != canvas.winfo_toplevel() ):
                return

        canvas.focus_set()


    def redrawHeader(self, *args):
        pass
    
    def yview(self, *args):
        # dragging scroll bar triggers('moveto', '0.024242424242424242')
        # using scrolbar arrows triggers ('scroll', '1', 'units')
        # clicking in scroll bar background ('scroll', '-1', 'pages')
        #
        ## callback for vertscrollbar_command
        # args can be ('scroll', number, type) where type is page or unit
        # or ('moveto', '0.2854982')
        #print 'YVIEW', args
        if args[0] == "scroll":
            # ('scroll', '1', 'pages') when click on back of bar
            # ('scroll', '1', 'units') when click on bar arrow
            #print 'SCROLL', args
            if args[1]=='1' and self.displayedNodes[-1].nextNode() is None:
                    return
            self.yview_scroll(args[1], args[2])
        else:
            # ('moveto', '0.350939') when dragging bar
            percent = float(args[1])
            line = min( self.nbLines-self.nbLinesPerPage,
                        int((self.nbLines-1) * float(percent)))
            if line > self.firstVisibleNodeLineNum:
                if self.displayedNodes[-1].nextNode() is None:
                    return
            if self.firstVisibleNodeLineNum != line:
                v = self.scrollView(line - self.firstVisibleNodeLineNum)
                if v:
                    # compute percentage of motion
                    total = self.nbLines+(self.headerHeight/float(self.nodeHeight))
                    percent = self.firstVisibleNodeLineNum/total
                    ScrolledCanvas.yview_moveto(self, percent)
                    self.redraw()
        
            #self.redrawHeader(args)

    def pageUp(self, event):
        # triggered by page up on keyboard
        self.yview_scroll(-1, "pages")
        return "break"


    def pageDown(self, event):
        # triggered by page down on keyboard
        self.yview_scroll(1, "pages")
        return "break"


    def lineUp(self, event):
        # triggered by up arrow on keyboard
        # print 'LINE UP'
        self.yview_scroll(-1, "units")
        return "break"


    def lineDown(self, event):
        # triggered by down arrow on keyboard
        #print 'LINE DOWN'
        self.yview_scroll(1, "units")
        return "break"


    def lineUpDown(self, event):
        #print 'LINE UPDOWN'
        if event.delta < 0:
            return self.lineDown(event)
        else:
            return self.lineUp(event)


    def yview_scroll(self, *args):
        # mouse wheel or arrow keys up and down args = (1, 'units'), (-1, 'units')
        # page up and down keys (1, 'pages')(13, 'units')
        #                       (-1, 'pages')(-13, 'units')
        # end key args = "jumping to  End"
        # home key args = "jumping to  Home"
        height = self.winfo_height() - self.headerHeight
        #print 'yview_scroll', args

        #height = int(self.canvas['height']) - self.headerHeight
        #print 'GGGG', height, self.nbLines, self.nodeHeight, self.nbLines * self.nodeHeight
        if self.nbLines * self.nodeHeight <= height:
            return

        if args[1] == "pages":
            # FIXME we scroll the whole thing one page
            linesPerPage = int(args[0]) * (height / self.nodeHeight)
            self.yview_scroll(linesPerPage, "units")
        else: # it is a unit scroll of args[0] units
            if args[0]=='1' and self.displayedNodes[-1].nextNode() is None:
                return
            line = min(self.nbLines-self.nbLinesPerPage,
                       self.firstVisibleNodeLineNum+int(args[0]))
            v = self.scrollView(line - self.firstVisibleNodeLineNum)
            if v:
                ScrolledCanvas.yview_scroll(self, *args)
                self.redraw()

        #self.redrawHeader(args)


    def scrollView(self, deltalines):
        """move the visible part of the tree up of down by deltalines lines
        """
        #print 'SCROLL VIEW', deltalines
        if deltalines > 0:
            # compute the index of the line for which the bottom of the tree
            # is drawn at the bottom of the visible window
            last = self.nbLines - self.nbLinesPerPage + 1
            # clamp deltalines so we do no go too far down
            deltalines = min(last-self.firstVisibleNodeLineNum+1, deltalines)
            #print 'GUGU', last, self.firstVisibleNodeLineNum, deltalines
            j = 0
            for i in range(deltalines):
                node = self.firstVisibleNode.nextNode()
                if node is None:
                    break
                self.firstVisibleNodeLineNum += 1
                self.firstVisibleNode = node
                j+=1
            return j
        else:
            j = 0
            for i in range(-deltalines):
                node = self.firstVisibleNode.previousNode()
                if node is None:
                    break
                self.firstVisibleNodeLineNum -= 1
                self.firstVisibleNode = node
                j+=1
            return j
        

    def showNode(self, node):
        self.firstVisibleNodeLineNum = node.childNumber
        self.firstVisibleNode = node
        self.redraw()

        
    def enter_cb(self, event=None):
        if widgetsOnBackWindowsCanGrabFocus is False:
            lActiveWindow = self.focus_get()
            if    lActiveWindow is not None \
              and ( lActiveWindow.winfo_toplevel() != self.winfo_toplevel() ):
                return

        self.canvas.focus_set()


    def configure_cb(self, event=None):
        nl = (self.winfo_height()-self.headerHeight) / self.nodeHeight
        self.nbLinesPerPage = nl
        last = self.nbLines - nl
        if self.firstVisibleNodeLineNum > last:
            self.scrollView(last - self.firstVisibleNodeLineNum)
        self.redraw()


    def destroy(self):
        if self.root:
            self.root.destroy()
        ScrolledCanvas.destroy(self)

    ##
    ##  Drawing
    ##
    def undisplay(self):
        for node in self.displayedNodes:
            node.deleteNodeIcons()

    def reparentNodes(self, tree):
        n = self.root
        while n:
            n.tree = ref(tree)
            n = n.nextNode()

    def redraw(self, force=False):
        """post a redraw or actually redraw depending on self.idleRedraaw
"""
        if self.suspendRedraw: return
        if self.root:
            if self.idleRedraw:
                if self.pending_draw:
                    self.after_cancel(self.pending_draw)
                cb = CallbackFunction(self.reallyRedraw, force=force)
                self.pending_draw = self.after(10, cb)
            else:
                self.reallyRedraw(force)


    def reallyRedraw(self, force=False):
        """actually redraw
        if force is True we for redrawing each line
        else we move existing lines
        """
        self.pending_draw = None
        canvas = self.canvas

        if force:
            # destroy representation of visible nodes
            for node in self.displayedNodes:
                #print 'destroying', node.object.name
                node.deleteNodeIcons()
        else:
            # make a copy of previousely displayed nodes
            prevDisplayNodes = {}.fromkeys(self.displayedNodes[:])

        # build list of visible nodes
        nodeHeight = self.nodeHeight
        nodes = []
        node = self.firstVisibleNode
        for i in range(self.nbLinesPerPage):
            if not force:
                if prevDisplayNodes.has_key(node):
                    del prevDisplayNodes[node]
            nodes.append(node)
            node = node.nextNode()
            if node is None:
                break

        # destroy representation of previously visible nodes
        # that are not longer visible
        if not force:
            for node in prevDisplayNodes.keys():
                #print 'destroying', node.object.name
                node.deleteNodeIcons()

##  this approach does not destoy some canvas items of node that move
## for instance after an expand

##         # destroy representation of nodes no longer visible
##         for node in self.displayedNodes:
##             if nodesd.has_key(node):
##                 print 'deleteIcon', node.object
##                 node.deleteNodeIcons()

        # Draw visible nodes
        ypos = self.firstVisibleNodeLineNum * nodeHeight + self.headerHeight
        for node in nodes:
            if node.labelTkid is None or node.currenty is None:
                node.redraw(ypos)
                #print 'DRAWING', node.object.name, ypos, node.currenty, node.uniqID
            else:
                dy = ypos - node.currenty
                #print 'MOVING', node.object.name, dy, self.firstVisibleNodeLineNum, ypos, node.currenty
                if dy:
                    canvas.move(node.nodeTkTag, 0, dy)
                    node.currenty += dy
                    
            #if node.currenty is None:
            #    node.currenty = 0
            #print 'redraw', node.needsRedraw, node.object.name
            #node.redraw(ypos)
            ypos += nodeHeight
        
        self.displayedNodes = nodes
      
        # Update the scroll-bar size
        if self.root:
            canvas = self.canvas
            x1, y1, x2, y2 = canvas.bbox(ALL)
            canvas.configure( scrollregion=(
                0, 0, x2, self.nbLines*nodeHeight+self.headerHeight))


    def updateTreeHeight(self):
        """ """
        self.nbLines = self.root.countSubtreeLines()


    def clearSelection(self, history=True):
        """Deselect all selected Nodes."""
        canvas = self.canvas
        if history:
            self.selectionHistory.append(self.selectedNodes[:])
        selectedNodes = self.selectedNodes[:]
        for node in selectedNodes:
            node.deselect(history=False)


    def undoSelect(self):
        """ """
        if len(self.selectionHistory):
            self.clearSelection(history=False)
            self.selectedNodes = self.selectionHistory.pop()
            for n in self.selectedNodes:
                n.select(only=False, history=False)


    def handleControlKey(self, event):
        """ """
        if event.keysym in ['z', 'Z']:
            self.undoSelect()

        elif event.keysym in ['l', 'L']:
            self.redraw()

##     def clearTree(self):
##         if self.pending_draw:
##             self.after_cancel(self.pending_draw)
##             self.pending_draw = None
      
##         if self.root:
##             self.root.destroy()
##             for n in self.displayedNodes:
##                 n.deleteNodeIcons()


##     def refresh(self):
##         """Refresf the tree. Notice that it is speeder to update only some Nodes (with Node.update() or Node.updatetree() and then tree.redraw()), if you know which Nodes have changed."""
##         if self.root:
##             self.root.refresh()
##             self.redraw()



class Node:
    """Base class for Nodes in a Tree"""

    def __init__(self, object, parent):
        """Create a Node.
object is the Python object represented by this node in the tree.
Object is expected to provide an interable sequence in its .children attribute.
Parent can be either the parent's Node, or the tree (for the root Node).
"""
        self.object = object
        self.isExpanded = False
        self.hasExpandIcon = True
        self.isSelected = False
        self.hasBeenExpanded = False
        self.children = []
        self.childNumber = 0 # index of this node in the list of children of its parent
        
        self.iconWidth = 0
        self.font = None
        
        self.currenty = None     # y value on canvas when node is visible
                                 # if None the node is not drawn
        self.needsRedraw = False
        
        self.parent = parent
        if parent is None:
            #print 'FOFOFOFOFOFOFOFOFO', object
            self.generation = 0   # how many ancestors
            self.tree = None      # weakref to tree
            self.uniqID = 0
            self.balloon = None
        else:
            self.generation = parent.generation + 1
            self.tree = parent.tree
            self.tree().maxNodeNum += 1
            self.balloon = self.tree().balloon
            
            # canvas ids
            self.uniqID = self.tree().maxNodeNum

        self.nodeTkTag = 'node_'+str(self.uniqID)
        self.labelTkid = None  # canvas id of node name
        self.iconTkid = None  # canvas id if node's icon
        self.expandIconTkid = None # canvas id of expand/collapse icon
        self.selectionBoxId = None
        self.canvasIDs = [] # list of other ids
        
    ##
    ## to be overridden
    ##
    def isExpandable(self):
        """Returns true if this node has children"""
        return 0


    def getChildren(self):
        """
        return children for object associated with this node.
        By default we return object.children. Override this method to
        selectively show children
        """
        return self.object.children

    
    def createChildrenNodes(self):
        """Create node for all children of self.object"""
        children = []
        for child in self.getChildren():
            if hasattr(child, 'treeNodeClass'):
                klass = child.treeNodeClass
            else:
                klass = self.__class__
            children.append(klass(child, self))
        return children


    ##
    ##  recursive traversals
    ##
    def countSubtreeLines(self):
        """Recusrsively traverse the tree and count the number of lines needed
to draw this subtree on a canvas"""
        nbLines = 1
        if self.isExpanded:
            for child in self.children:
                nbLines = nbLines + child.countSubtreeLines()
        return nbLines


    def lastVisibleChild(self):
        """Return the last visible child of this node"""
        if self.hasBeenExpanded and self.isExpanded:
            return self.children[-1].lastVisibleChild()
        else:
            return self


    def findNextNode(self, node):
        """recursively traverse the tree to find the first node visible
below self"""
        # if this is not the last child, return the next one
        if node.childNumber + 1 < len(self.children):
            return self.children[node.childNumber + 1]

        if self.parent: # we need to walk the tree to find te next node
            return self.parent.findNextNode(self)
        else: # root node
            return None


    ##
    ## tree navigation methods
    ##
    def previousNode(self):
        """return the node right above this one in the tree"""
        if self.parent is None: # root node
            return None
        if self.childNumber==0: # node is first child, return parent
            return self.parent
        else: # return last visible child of sibling before self
            return self.parent.children[self.childNumber-1].lastVisibleChild()


    def nextNode(self):
        """return the node right below this one in the tree"""
        if self.isExpanded and self.children:
            return self.children[0]
        else:
            if self.parent:
                return self.parent.findNextNode(self)
            else:
                return None # Root Node


    ##
    ##  expanding and collapsing the tree
    ##
    def expand(self, event = None):
        """Expand this node and show its children"""
        if self.isExpanded:
            return

        tree = self.tree()
        if not self.hasBeenExpanded:
            self.children = self.createChildrenNodes()
            for i,c in enumerate(self.children):
                c.childNumber = i
                tree.objectToNode[c.object] = c
            self.hasBeenExpanded = True

        if len(self.children):
            self.needsRedraw = True
            self.isExpanded = True
            tree.updateTreeHeight()
            self.deleteNodeIcons()
            tree.redraw()


    def collapse(self, event = None):
        """Collapse this node, i.e hid its children"""
        if not self.isExpanded:
            return
        
        self.isExpanded = False
        self.needsRedraw = True
        tree = self.tree()
        tree.updateTreeHeight()
        self.deleteNodeIcons()
        tree.redraw()


    def doubleLabel1(self, event=None):
        """Callback for double click on label with button 1"""
        self.toggleExpansion(event)


    def toggleExpansion(self, event=None):
        """Toggles expanded/collapsed state of this node"""
        if self.isExpanded:
            self.collapse()
        else:
            self.expand()


    ##
    ##  Drawing
    ##
    def getIcon(self):
        """return node's icons"""
        iconsManager = self.tree().iconsManager
        if hasattr(self.object, 'treeIconName'):
            icon = iconsManager.get(self.object.treeIconName,
                                    self.tree().master)
        elif self.isExpandable():
            if self.isExpanded:
                icon = iconsManager.get("expandedIcon.pgm")
            else:
                icon = iconsManager.get("collapsedIcon.pgm")
        else:
            icon = iconsManager.get("python.pgm")

        if icon:
            self.iconWidth = icon.width()
        else:
            self.iconWidth = 0
        return icon


    def drawExpandCollapseIcon(self, x, y):
        """Draw the icon used to expand and collapse nodes"""
        tree = self.tree()
        canvas = tree.canvas

        self.nodeStartX = x

        # delete old icon
        if self.expandIconTkid:
            canvas.delete(self.expandIconTkid)

        # draw new one
        if self.isExpandable() and self.hasExpandIcon:
            if self.isExpanded:
                icon = tree.expandedIcon
            else:
                icon = tree.collapsedIcon

            tkid = self.expandIconTkid = canvas.create_image(
                x + tree.iconHalfWidth, y, image=icon, tags=(self.nodeTkTag,))
            canvas.tag_bind(tkid, "<Button-1>", self.toggleExpansion)
            iconWidth = icon.width()
            return iconWidth+tree.iconHalfWidth
        else:
            iconWidth = 0
            self.expandIconTkid = None
            return 0


    def drawNodeIcon(self, x, y):
        """Draw the node's icon"""
        tree = self.tree()
        canvas = tree.canvas

        if self.iconTkid:
            canvas.delete(self.iconTkid)

        icon = self.getIcon()
        if icon:
            self.iconTkid = canvas.create_image(
                x, y, image=icon, anchor="nw",
                tags=(self.nodeTkTag,))
            canvas.tag_bind(self.iconTkid, "<Double-1>", self.toggleExpansion)
            canvas.tag_bind(self.iconTkid, "<Shift-Button-1>",
                            self.selectRange_cb)
            canvas.tag_bind(self.iconTkid, "<Control-Button-1>",
                            self.modifySelection_cb)
            canvas.tag_bind(self.iconTkid, "<Button-1>", self.toggleSelection)
            if os.name != 'nt': #sys.platform != 'win32':
                canvas.tag_bind(self.iconTkid, "<Button-4>", self.tree().lineUp)
                canvas.tag_bind(self.iconTkid, "<Button-5>", self.tree().lineDown)
        else:
            self.iconTkid = None

        return self.iconWidth

    def fitText(self, x, y, tree, canvas, text, font=None):
        # added so that NodeWithButton can change label to fit tree width
        return text
        

    def drawNodeLabel(self, x, y):
        """Draw the node's label"""
        tree = self.tree()
        canvas = tree.canvas

        origtext = self.object.name
        if not origtext:
            self.needsRedraw = True
            return 0

        if self.font:
            font = self.font
        else:
            font = self.tree().font

        text = self.fitText(x, y, tree, canvas, origtext, font)
        if self.labelTkid:
            curText = canvas.itemconfig(self.labelTkid)['text'][4]
        else:
            curText = ''
        if curText==text:
            canvas.coords(self.labelTkid, x, y)
            return 0
        
        if self.labelTkid:
            canvas.delete(self.labelTkid)

        ## if len(text) < len(origtext):
        ##     balloon = Pmw.Balloon(canvas)
        ## else:
        ##     balloon = None

        labelTkid = canvas.create_text(
            x, y, anchor='nw', text=text, font=font, tags=(self.nodeTkTag,))

        canvas.tag_bind(labelTkid, "<Button-1>", self.toggleSelection)
        canvas.tag_bind(labelTkid, "<Shift-Button-1>", self.selectRange_cb)
        canvas.tag_bind(labelTkid, "<Control-Button-1>", self.modifySelection_cb)
        canvas.tag_bind(labelTkid, "<Button-3>", self.button3OnLabel)
        canvas.tag_bind(labelTkid, "<Double-Button-1>", self.doubleLabel1)
        canvas.tag_bind(labelTkid, "<Motion>", tree.move_cb)
        if os.name == 'nt': #sys.platform == 'win32':
            #this does not work:
            #canvas.tag_bind(labelTkid, "<MouseWheel>", tree.lineUpDown)
            canvas.tag_bind(labelTkid, "<Button-4>", tree.lineUpDown)
        else:
            canvas.tag_bind(labelTkid, "<Button-4>", tree.lineUp)
            canvas.tag_bind(labelTkid, "<Button-5>", tree.lineDown)
        self.labelTkid = labelTkid

        if self.balloon:
            self.balloon.tagbind(canvas, labelTkid, origtext)
            
        bb = canvas.bbox(self.labelTkid)
        return bb[2]-bb[0]


    def drawNodeCustomization(self, x, y):
        """Draw additional things on the canvas"""
        return 0
  

    def nodeRedraw(self):
        # this functions performs a redraw of a node in the same place
        ypos = self.currenty
        self.deleteNodeIcons()
        self.redraw(ypos)

        
    def redraw(self, y):
        """Redraw this node at position y on the canvas"""
        # this function actually draws the node
        
        #print 'AAAA Redraw', self.object.name
        #if self.currenty is None:
        #    return

        tree = self.tree()
        if y != self.currenty or self.needsRedraw:
            x = 2 + max(0, self.generation) * tree.nodeHeight
            if not hasattr(self, 'nodeTkTag'):
                self.nodeTkTag = 'node_'+str(self.uniqID)
            x += self.drawExpandCollapseIcon(x, y+tree.iconHalfWidth)
            x += self.drawNodeIcon(x, y)
            x += self.drawNodeLabel(x, y)
            x += self.drawNodeCustomization(x, y)
            self.currenty = y
            if self.isSelected:
                if self.selectionBoxId:
                    tree.canvas.delete(self.selectionBoxId)
                self.drawSelectionBox()
        self.needsRedraw = False


    def deleteNodeIcons(self):
        """Delect canvas items representing this node"""

        canvas = self.tree().canvas
        canvas.delete(self.nodeTkTag)
        self.expandIconTkid = None
        self.labelTkid = None
        self.selectionBoxId = None
            
        ## if self.expandIconTkid:
        ##     canvas.delete(self.expandIconTkid)
        ##     self.expandIconTkid = None

        ## if self.labelTkid:
        ##     canvas.delete(self.labelTkid)
        ##     self.labelTkid = None

        ## if self.selectionBoxId:
        ##     canvas.delete(self.selectionBoxId)
        ##     self.selectionBoxId = None
        ##     self.labelTkid = None

        #for id in self.canvasIDs:
        #    canvas.delete(id)
        self.canvasIDs = []

        self.currenty = None


    def destroy(self):
        for child in self.children:
            child.destroy()


    ##
    ## selection
    ##
    def drawSelectionBox(self):
        """draw a yellow box"""
        tree = self.tree()
        canvas = tree.canvas
        w = tree.winfo_width()
        y = self.currenty
        id = canvas.create_rectangle(0, y-2, w, y+20, outline='yellow',
                                     fill='yellow', tags=(self.nodeTkTag,)) 
        canvas.lower(id)
        self.selectionBoxId = id


    def button3OnLabel(self, event=None):
        # override in subclass
        return

    
    def toggleSelection(self, event=None, only=True):
        tree = self.tree()
        tree.lastPickedNode = self
        if self.isSelected:
            self.deselect()
        else:
            self.select(only=only)


    def select(self, only=True, history=True):
        if self.isSelected:
            return

        tree = self.tree()
        if history:
            tree.selectionHistory.append(tree.selectedNodes[:])
        if only:
            tree.clearSelection(history=False)
        self.isSelected = True
        tree.selectedNodes.append(self)
        if self in tree.displayedNodes:
            self.drawSelectionBox()
        tree.dispatchEvent( SelectEvent(self))


    def deselect(self, history=True):
        if not self.isSelected:
            return
        tree = self.tree()
        if history:
            tree.selectionHistory.append(tree.selectedNodes[:])
        tree.selectedNodes.remove(self)
        if self in tree.displayedNodes:
            tree.canvas.delete(self.selectionBoxId)

        self.isSelected = False
        self.selectionBoxId = None
        tree.dispatchEvent( DeselectEvent(self))


    def selectRange_cb(self, event=None):
        tree = self.tree()
        if tree.selectionMode=='single':
            return

        if len(tree.selectedNodes)==0:
            return

        last = tree.lastPickedNode
        all = self.parent.children
        try:
            i1 = all.index(self)
            i2 = all.index(last)
        except ValueError:
            raise ValueError, "range only work over siblings"

        if i1 > i2:
            tmp=i1; i1=i2; i2=tmp

        tree.selectionHistory.append(tree.selectedNodes[:])

        #print last, i1, i2+1
        for node in all[i1:i2+1]:
            if node==last:
                continue
            if node.isSelected:
                node.deselect(history=False)
            else:
                node.select(only=False, history=False)

    
    def modifySelection_cb(self, event=None):
        tree = self.tree()
        tree.lastPickedNode = self
        if tree.selectionMode=='single':
            return
        self.toggleSelection(event=event, only=False)


    def undoSelect(self, event=None):
        tree.lastPickedNode = None
        self.tree().undoSelect()


    def refresh(self):
        self.needsRedraw = True
        self.tree().canvas.itemconfig(self.labelTkid, text=unicode(self))


    def refreshSubTree(self):
        self.refresh()
        for child in self.children:
            child.refreshSubTree()
        self.tree().redraw()
        

    def refreshChildren(self, redraw=True):
        ## bizarre
        if self.hasBeenExpanded:
            if self.isExpanded:
                tree = self.tree()
                # save current children
                oldchildren = self.children
                # delete all icons for current children
                for child in self.children:
                    child.deleteNodeIcons()
                # create Nodes for new children while saving existing ones
                oldChildObj = {}
                for c in oldchildren: # dict of obj:existing nodes
                    oldChildObj[c.object]= c
                i = 0
                newchildren = [] # build the list and create new children
                for childobj in self.getChildren():
                    try:
                        c = oldChildObj[childobj]
                    except KeyError:
                        if hasattr(childobj, 'treeNodeClass'):
                            klass = childobj.treeNodeClass
                        elif  hasattr(childobj, 'molIteratorNodeClass'):
                            klass = childobj.molIteratorNodeClass
                        else:
                            klass = self.__class__
                        c = klass(childobj, self)
                        tree.objectToNode[c.object] = c

                    newchildren.append(c)
                    c.childNumber = i
                    i = i + 1
                self.children = newchildren

                if len(self.children)==0:
                    self.isExpanded = 0 # Cannot be expanded if no children
                if redraw:
                    self.tree().redraw()
            else:
                for child in self.children:
                    child.destroy()
                self.children = []
                self.hasBeenExpanded = False
                if redraw and self.currenty:
                    self.redraw(self.currenty)
        else:
            if redraw and self.currenty:
                self.redraw(self.currenty)


##     def selectTree(self):
##         if not self.isSelected: self.select()
##          for child in self.children:
##              child.selectTree()


##     def deselectTree(self):
##         self.deselect()
##         for child in self.children:
##             child.deselectTree()

## #root = Tk()
## #iconsmanager = IconsManager('Icons/32x32/', 'Pmv')
## #ico = iconmanager.get('ss.png', root)

if __name__=='__main__':

    from mglutil.util.callback import CallbackFunction

    from MolKit.molecule import Atom, Molecule
    from MolKit.protein import Residue, Chain

    class ObjectTree(Tree):
        """Each node in the tree has an object associated in the node's .object
    attribute.  The objects are expected to have a .parent and a .children
    attribute describing the hierarchy."""

        def __init__(self, master, root, iconsManager=None, idleRedraw=True,
                     nodeHeight=15, **kw):
            Tree.__init__(self, master, root, iconsManager=None, idleRedraw=True,
                          nodeHeight=nodeHeight, **kw)
            self.objectToNode = {}  # key is object, values if Node instance
            self.objectToNode[root.object] = root

            self.nbCol = 10
            self.menuEntries = []
            for i in range(self.nbCol):
                self.menuEntries.append([])
            self.menuEntries[0] = [
                'displayLines', 'display S&B', 'display CPK',
                'display second. struct.', 'display molecular surface',
                'display second and S&B'
                'undisplayLines', 'undisplay S&B', 'undisplay CPK',
                'undisplay second. struct.', 'undisplay molecular surface',
                'undisplay second and S&B'
                ]
            self.menuEntries[1] = [
                'label Atoms', 'label residues', 'label chains', 'label molecules'
                ]
            self.menuEntries[2] = [
                'color by atom types', 'color by molecule', 'color by chains',
                'color by residue (RASMOL)', 'color by residue (SHAPELY)',
                'color by DG', 'color by instance', 'color by second. struct.'
                ]

        def createNodeForObject(self, object):
            """given an object if the corresponding Node is not yet created
    for its creation by expanding all its ancestors"""
            try:
                return self.objectToNode[object]
            except KeyError:
                p = object.parent
                ancestors = [p]
                while p and not self.objectToNode.get(p, None):
                    ancestors.append(p)
                    p = p.parent
                ancestors.append(p)
                ancestors.reverse()
                for p in ancestors:
                    self.objectToNode[p].expand()
                return self.objectToNode[object]


    class ObjectNode(Node):
        """ """

        def __repr__(self):
            return self.object.name

        def getIcon(self):
            """return node's icons"""
            iconsManager = self.tree().iconsManager
            object = self.object
            
            if isinstance(object, Atom):
                icon = iconsManager.get("atom.png", self.tree().master)
            elif isinstance(object, Residue):
                icon = iconsManager.get("residue.png", self.tree().master)
            elif isinstance(object, Chain):
                icon = iconsManager.get("chain.png", self.tree().master)
            elif isinstance(object, Molecule):
                icon = iconsManager.get("ms.png", self.tree().master)
            else:
                icon = None
                
            if icon:
                self.iconWidth = icon.width()
            else:
                self.iconWidth = 0
            return icon



        def isExpandable(self):
            """Returns true if this node has children"""
            return len(self.getChildren())


        def createChildrenNodes(self):
            """Create node for all children of self.object"""
            children = []
            for child in self.getChildren():
                children.append(ObjectNode(child, self))
            return children


        def drawNodeCustomization(self, x, y):
            """Add things to the rigth side of the tree
    """
            tree = self.tree()
            canvas = tree.canvas
            level = 2
            col = ['gray75', 'red', 'cyan', 'green', 'yellow']

            nbButtons = 10
    ##         if not hasattr(self, 'chkbt'):
    ##             self.chkbtVar = []
    ##             self.chkbt = []
    ##             for i in range(nbButtons):
    ##                 v = IntVar()
    ##                 self.chkbtVar.append(v)
    ##                 cb = CallbackFunction(self.buttonClick, i )
    ##                 button = Checkbutton(
    ##                     canvas, variable=v, command=cb, padx=0, pady=0,
    ##                     background=col[level-1], anchor='nw')
    ##                 self.chkbt.append(button)

            x = 150
            fill = ['white', 'green']
            if not hasattr(self, 'chkbt'):
                nbButtons = tree.nbCol
                self.chkbt = [0]*nbButtons
                self.chkbtid = [None]*nbButtons

                # create a pull down menu and allocate variables for each column
                self.menu = Menu(canvas, title='Choose', tearoff=False)
                self.menuVars = []  # list of [column][menu entries] IntVar
                for i, entries in enumerate(tree.menuEntries):
                    l = []
                    for j in range(len(entries)):
                        l.append(IntVar())
                    self.menuVars.append(l) # a list for each column

            for i, val in enumerate(self.chkbt):
                xo = x+i*35
                cid = canvas.create_oval(xo, y, xo+15,y+15, fill=fill[val])
                cb = CallbackFunction(self.buttonClick, i )
                canvas.tag_bind(cid, "<Button-1>", cb)
                cb = CallbackFunction(self.shiftButtonClick, i )
                canvas.tag_bind(cid, "<Shift-Button-1>", cb)
                self.chkbtid[i] = cid
                self.canvasIDs.append(cid)
    ##             button = self.chkbt[i]
    ##             cid = canvas.create_window( x+i*35, y, window=button,
    ##                                         width=20, height=15)
    ##             self.canvasIDs.append(cid)

                # add secondary structure glyph
                molFrag = self.object
                if isinstance(molFrag, Residue):
                    if hasattr(molFrag, 'secondarystructure'):
                        ssname = molFrag.secondarystructure.name
                        if ssname[:6]=='Strand': color = '#FFF700'
                        elif ssname[:4]=='Coil': color = 'grey45'
                        elif ssname[:5]=='Helix': color = '#FF198C'
                        elif ssname[:4]=='Turn': color = 'blue'

                        cid = canvas.create_rectangle(
                            130, y-10, 140, y+10, outline=color,fill=color)
                        self.canvasIDs.append(cid)

                        cid = canvas.create_text(
                            152, y, text=ssname, anchor='nw',
                            font=self.tree.font)
                        self.canvasIDs.append(cid)
            return x+i*35
    ##             func = tree.buttonValFunc[i]
    ##             if func:
    ##                 func(self)

        def menu_cb(self, column, what, menuEntryIndex, event=None):
            print 'FFFF', column, what, menuEntryIndex


        def buttonClick(self, column, event=None):
            # get called for each checkbutton
            tree = self.tree()

            if self.chkbt[column]==0:
                tree.canvas.itemconfigure(self.chkbtid[column], fill='green')
                self.chkbt[column] = 1
            else:
                tree.canvas.itemconfigure(self.chkbtid[column], fill='white')
                self.chkbt[column] = 0

            for i,v in enumerate(self.menuVars[column]):
                if v.get():
                    print 'DDD', tree.menuEntries[column][i], self.chkbt[column]


        def shiftButtonClick(self, column, event=None):
            # get called for each checkbutton
            tree = self.tree()
            menu = self.menu

            if menu.index(0)==0: # there is something in the menu, remove it
                menu.delete(0, 100)

            for i, menuEntry in enumerate(tree.menuEntries[column]):
                v = self.menuVars[column][i]
                cb = CallbackFunction(self.menu_cb, column, menuEntry, i)
                menu.add_checkbutton(label=menuEntry, variable=v, command=cb)
            menu.post(event.x_root, event.y_root)

    ##         if self.chkbt[column]==0:
    ##             tree.canvas.itemconfigure(self.chkbtid[column], fill='green')
    ##             self.chkbt[column] = 1
    ##         else:
    ##             tree.canvas.itemconfigure(self.chkbtid[column], fill='white')
    ##             self.chkbt[column] = 0

    ##   def getIcon(self):
    ##       """return node's icons"""
    ##       return None

    from MolKit import Read
    from MolKit.molecule import MolecularSystem
    syst = MolecularSystem ('world')
    mols = Read('../dev23/1crn.pdb')
    syst.adopt(mols[0])
    
    #mols = Read('../dev23/2plv.pdb')
    #syst.adopt(mols[0])
    #mols = Read('../dev23/1gav.pdb')

    root = Tk()
    rootnode = ObjectNode(syst, None)
    tree = ObjectTree(root, rootnode, selectionMode='multiple')
    tree.pack(expand=1, fill="both")
