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

## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
# Revision: Guillaume Vareille
#
#############################################################################
#
# $Header: /mnt/raid/services/cvs/DejaVu2/Tk/ViewerGUI.py,v 1.1.1.1.4.1 2017/07/13 22:22:07 annao Exp $
#
# $Id: ViewerGUI.py,v 1.1.1.1.4.1 2017/07/13 22:22:07 annao Exp $
#
 
import types, numpy, Tkinter, Pmw, string, os
import numpy.oldnumeric as Numeric
from Tkinter import *

from opengltk.OpenGL import GL
from Slider import Slider
from DejaVu2.ColorChooser import ColorChooser
from DejaVu2 import viewerConst
from DejaVu2 import jitter
from DejaVu2 import colorTool
from DejaVu2 import PropertyEditor

from mglutil.util.callback import CallBackFunction
from mglutil.util.colorUtil import TkColor
from mglutil.gui.BasicWidgets.Tk.fileBrowsers import FileOpenBrowser, \
     FileSaveBrowser
from mglutil.gui.BasicWidgets.Tk.RelabelingCascadeMenu import RelabelingCascadeMenu
import tkMessageBox
from mglutil.gui.BasicWidgets.Tk.thumbwheel import ThumbWheel
from DejaVu2.Transformable import Transformable
from DejaVu2.IndexedGeom import IndexedGeom
from DejaVu2.Geom import Geom
from DejaVu2.SelectionGUI import SelectionGUI
import DejaVu2
from DejaVu2.Spheres import Spheres


class TwoSlidersGUI(Frame):
    """Class for a gui with two sliders
"""
    def mouseDown(self, event):
        # remember where the mouse went down
        self.draw.itemconfig(CURRENT, fill="red")
        self.draw.addtag('selected', 'withtag', CURRENT)
        self.lastx = event.x


    def mouseUp(self, event):
        # remember where the mouse went down
        self.draw.itemconfig("selected", fill="blue")
        self.draw.dtag("selected")


    def mouseMove(self, event):
        # whatever the mouse is over gets tagged as CURRENT for free by tk.
         t = event.widget.gettags(CURRENT)
         deltax = event.x - self.lastx
         if 'near' in t:
            newx = self.nearx+deltax
            if newx > self.left-6 and newx < self.farx:
                self.nearx = self.nearx + deltax
                self.near = self.cst1 * (self.nearx-self.left)
                self.draw.move(CURRENT, deltax, 0 )
                self.lastx = self.nearx
                if self.near < self.far:
                    self.viewer.currentCamera.Set(near=self.near)
                self.viewer.Redraw()
         elif 'far' in t:
            newx = self.farx+deltax
            if newx > self.nearx :#and newx < 191:
                self.farx = self.farx + deltax
                self.far = self.cst1 * (self.farx-self.left)
                self.draw.move(CURRENT, deltax, 0 )
                self.lastx = self.farx
                if self.far > self.near:
                    self.viewer.currentCamera.Set(far=self.far)
                self.viewer.Redraw()
         elif 'start' in t:
            newx = self.startx+deltax
            if newx > self.left-6 and newx < self.endx:
                self.startx = self.startx + deltax
                self.start = self.cst1 * (self.startx-self.left)
                self.draw.move(CURRENT, deltax, 0 )
                self.lastx = self.startx
                if self.start < self.end:
                    self.viewer.currentCamera.fog.Set(start=self.start)
                self.viewer.Redraw()
         elif 'end' in t:
            newx = self.endx+deltax
            if newx > self.startx :#and newx < 191:
                self.endx = self.endx + deltax
                self.end = self.cst1 * (self.endx-self.left)
                self.draw.move(CURRENT, deltax, 0 )
                self.lastx = self.endx
                if self.end > self.start:
                    self.viewer.currentCamera.fog.Set(end=self.end)
                self.viewer.Redraw()


    def createWidgets(self):
        self.draw = Canvas(self, width=300, height=30)
        self.draw.pack(expand=1, fill='x')
        self.left = 40
        self.draw.create_line( 40, 15, 9999, 15, width=2)

        fnt='-*-helvetica-medium-r-narrow-*-*-120-*-*-*-*-*-*'
        self.draw.create_text( 16, 7, text='clipZ', font=fnt )
        self.draw.create_text( 15, 22, text='Fog', font=fnt)

        self.nearC = self.draw.create_polygon( 1, 1, 10, 1, 5, 10, 1, 1,
                                               fill='blue', tag='near')
        self.farC = self.draw.create_polygon( 1, 1, 10, 1, 5, 10, 1, 1,
                                              fill='blue', tag='far')
        self.startC = self.draw.create_polygon( 1, 10, 5, 1, 10, 10, 1, 10,
                                                fill='blue', tag='start')
        self.endC = self.draw.create_polygon( 1, 10, 5, 1, 10, 10, 1, 10,
                                              fill='blue', tag='end')

        Widget.bind(self.draw, "<1>", self.mouseDown)
        Widget.bind(self.draw, "<B1-Motion>", self.mouseMove)
        Widget.bind(self.draw, "<ButtonRelease-1>", self.mouseUp)


    def ComputeCst(self):
        self.cst = 155.0 / (self.max-self.min)
        self.cst1 = 1.0 / self.cst


    def GetX(self, val):
        """Compute the slider x offset in pixels for a given value"""
        return int( (self.left-5) + (val*self.cst) )

    def Set(self, near, far, start, end):
        """Set the cursors
"""
        if not near is None:
            if hasattr(self, 'nearx'):
                lPrev = self.nearx
                lYpos = 0
            else:
                lPrev = 0
                lYpos = 5
            self.nearx = self.GetX(near)
            self.near = near
            self.draw.move( self.nearC, self.nearx - lPrev, lYpos)
        if not far is None:
            if hasattr(self, 'farx'):
                lPrev = self.farx
                lYpos = 0
            else:
                lPrev = 0
                lYpos = 5
            self.farx = self.GetX(far)
            self.far = far
            self.draw.move( self.farC, self.farx - lPrev, lYpos)
        if not start is None:
            if hasattr(self, 'startx'):
                lPrev = self.startx
                lYpos = 0
            else:
                lPrev = 0
                lYpos = 15
            self.startx = self.GetX(start)
            self.start = start
            self.draw.move( self.startC, self.startx - lPrev, lYpos)
        if not end is None:
            if hasattr(self, 'endx'):
                lPrev = self.endx
                lYpos = 0
            else:
                lPrev = 0
                lYpos = 15
            self.endx = self.GetX(end)
            self.end = end
            self.draw.move( self.endC, self.endx - lPrev, lYpos)


    def __init__(self, viewer, master=None, min=0.0, max=100.0):
        Frame.__init__(self, master)
        Pack.config(self, expand=1, fill='x')
        self.max=max
        self.min=min
        self.viewer = viewer
        self.ComputeCst()
        self.createWidgets()
        
        if self.viewer:
            c = self.viewer.currentCamera
            f = self.viewer.currentCamera.fog
            self.Set(c.near, c.far, f.start, f.end)



class ViewerGUI:

    def __init__(self, viewer, maxLights, maxClip, name='Viewer',
                 nogui=0, master=None):

        if master is None:
            master = Toplevel()
            master.geometry('+80+180')
            master.title(name)
            
        self.root = master
        f = CallBackFunction(self.dialog)
        if isinstance(self.root,Toplevel):
            self.root.protocol('WM_DELETE_WINDOW', f)
        self.shown = True

        fnt = '-*-helvetica--r-narrow-*-*-140-*-*-*-*-*-*'
        self.viewer = viewer

        #self.addOcclusionCamera()

        self.top = Frame(self.root)
        # create menu bar
        self.mBar = Frame(self.top, relief=RAISED, borderwidth=2)
        self.mBar.pack(side=TOP, fill=X, expand=0)
    
        # FILE
        # create menu button
        self.menuFile = Menubutton(self.mBar, text='File', underline=0 )
        self.menuFile.pack(side=LEFT, padx="1m")
        # create pull down menuand add entries
        self.menuFile.menu = Menu(self.menuFile)

        # File Browser obj
        self.openFileBrowser = FileOpenBrowser(#lastDir='.',
                filetypes=[('all', '*'), ('py', '*.py')], title='Choose File') 
        self.saveFileBrowser = FileSaveBrowser(#lastDir='.',
                filetypes=[('all', '*'), ('py', '*.py')], title='Choose File')
    
        items = ["save Viewer's State", "Save Object's state",
                 "save Viewer and Objects states", "restore state",
                 'Load Transformation', 'Save Transformation',
                 'Save Representations', 'Load Representations']
        cmds = [self.saveViewerState_cb,
                self.saveObjectsStates_cb,
                self.saveViewerAndObjectsStates_cb,
                self.restoreState,
                self.loadTransformCurrentGeom, self.saveTransformCurrentGeom,
                self.saveReprToFile, self.loadReprFromFile]
        if not nogui:
            items.append('Exit')
            cmds.append(self.Exit_cb)
        for m in range(len(items)):
            self.menuFile.menu.add_command(label=items[m], underline=0,
                                               command = cmds[m] )
                 
        # attach pull down menu to button
        self.menuFile['menu'] = self.menuFile.menu
    
        
        # EDIT
        # create menu button
        self.menuEdit = Menubutton(self.mBar, text='Edit', underline=0 )
        self.menuEdit.pack(side=LEFT, padx="1m")
        # create pull down menuand add entries
        self.menuEdit.menu = Menu(self.menuEdit)
    
        items = ('Apply Transformation',
                     )
        cmds = (self.applyTransformation_cb,
                    )
    
        for m in range(len(items)):
            self.menuEdit.menu.add_command(label=items[m], underline=0,
                                               command = cmds[m] )
        # attach pull down menu to button
        self.menuEdit['menu'] = self.menuEdit.menu

    
        # PREF
        self.menuPref = Menubutton(self.mBar, text='Preferences', underline=0)
        self.menuPref.pack(side=LEFT)
        self.menuPref.menu = Menu(self.menuPref)
        items = (##'WYSIWIG Color Editor',
                     'Transf. Root Only',
                     'Show Picked Vertex',
                     'Display Value in the Object List',
                     'Display Quick Keys panel',
                )
        cmds = (##self.Wysiwyg_cb,
                    self.MoveRootOnly_cb,
                    self.showPickedVertex_cb,
                    self.displayValueInObjList_cb,
                    self.showHideQuickKeys_cb,
               )
        ##self.wysiwyg = IntVar()
        self.moveRootOnly = IntVar()
        self.showPickedVertex = IntVar()
        self.displayValueInObjList = IntVar()
        self.showHideQuickKeysVar = IntVar()

        vars = (##self.wysiwyg,
                    self.moveRootOnly, 
                    self.showPickedVertex, \
                    self.displayValueInObjList, self.showHideQuickKeysVar,
               )
        for m in range(len(items)):
            self.menuPref.menu.add_checkbutton(label=items[m],
                               var = vars[m],
                               command = cmds[m],
                               onvalue=viewerConst.YES,
                               offvalue=viewerConst.NO )
        self.menuPref.menu.add_command(label='Add Quick Key',
                                       command=self.addQuickKey_cb)

        self.menuPref['menu'] = self.menuPref.menu
    
        # HELP
        self.menuHelp = Menubutton(self.mBar, text='Help', underline=0)
        self.menuHelp.pack(side=RIGHT)
        self.menuHelp.menu = Menu(self.menuHelp)
        self.menuHelp['menu'] = self.menuHelp.menu

        self.QuickKeysFrame = Tkinter.Toplevel(width=50, height=100)
        self.QuickKeysFrame.withdraw()
        self.QuickKeysFrame.title('Quick Keys')
        #self.QuickKeysFrame.forget()

        #
        # Radio button for transformation binding
        #
        self.Xform = viewer.Xform

        import Pmw
        self.sf = Pmw.ScrolledFrame(self.root,horizflex='expand',vertflex='expand')
        self.sframe = self.sf.interior()
        frame1 = Tkinter.Frame(self.sframe)

        self.XformRadio = Pmw.Group(frame1, tag_text="Mouse transforms:")
        w = self.XformRadio.interior()
        lFrameTransformButtons = Frame(w)
        Radiobutton(lFrameTransformButtons, text="Object",
                variable=self.Xform, command = self.TObject,
                value="Object", anchor=W,
                indicatoron=0, selectcolor="#c67171",
                activebackground="#28686b").pack(side=LEFT)
        Radiobutton(lFrameTransformButtons, text="Camera",
                    variable=self.Xform, command = self.TCamera,
                    value="Camera", anchor=W,
                    indicatoron=0, selectcolor="#c67171",
                    activebackground="#28686b").pack(side=LEFT)
        Radiobutton(lFrameTransformButtons, text="Clip",
                    variable=self.Xform, command = self.TClip,
                    value="Clip", anchor=W,
                    indicatoron=0, selectcolor="#c67171",
                    activebackground="#28686b").pack(side=LEFT)
        Radiobutton(lFrameTransformButtons, text="Light",
                    variable=self.Xform, command = self.TLight,
                    value="Light", anchor=W,
                    indicatoron=0, selectcolor="#c67171",
                    activebackground="#28686b").pack(side=LEFT)
        Radiobutton(lFrameTransformButtons, text="Texture",
                    variable=self.Xform, command = self.TMap,
                    value="Texture", anchor=W,
                    indicatoron=0, selectcolor="#c67171",
                activebackground="#28686b").pack(side=LEFT)
        Radiobutton(lFrameTransformButtons, text="Scissor",
                    variable=self.Xform, command = self.Scissor,
                    value="Scissor", anchor=W,
                    indicatoron=0, selectcolor="#c67171",
                    activebackground="#28686b").pack(side=LEFT)
        lFrameTransformButtons.pack(side='top')

        #set mouse transformInfo
        lFrameTransformInfo = Frame(w)
        Label(lFrameTransformInfo, text='left').grid(row=0, column=0, sticky='w')
        Label(lFrameTransformInfo, text='middle').grid(row=0, column=1, sticky='w')
        Label(lFrameTransformInfo, text='right').grid(row=0, column=2, sticky='w')
        Label(lFrameTransformInfo, text=' ').grid(row=0, column=3, sticky='w')
        Label(lFrameTransformInfo, text='wheel').grid(row=0, column=4, sticky='w')
        Label(lFrameTransformInfo, text='zoom',
                                   relief='sunken',
                                   borderwidth=1,
                                   anchor='w').grid(row=1, column=4, sticky='w')
        self.transformInfo = {}
        for lButtonIndex in (1, 2, 3):
            self.transformInfo[lButtonIndex] = Tkinter.Label(
                                        lFrameTransformInfo, 
                                        width=10,
                                        relief='sunken',
                                        borderwidth=1,
                                        anchor='w')
            self.transformInfo[lButtonIndex].grid(row=1, column=lButtonIndex-1, sticky='w')
        lFrameTransformInfo.pack(side='top')
        self.fillTransformInfo_cb(self.Xform.get())
        #self.bindModifersToTransformInfo(master)
        for lCamera in self.viewer.cameras:
            self.bindModifersToTransformInfo(lCamera.master)

        # root only checkbutton
#        self.moveRootOnly = IntVar()
        Checkbutton(w,
                    text='mouse transforms apply to "root" object only',
                    variable=self.moveRootOnly,
                    command=self.MoveRootOnly_cb).pack(side=TOP)       

##         if DejaVu2.enableScenario:
##             addScenarioButton = True
##             if hasattr(self.viewer, "addScenarioButton"):
##                 addScenarioButton = self.viewer.addScenarioButton
##             if addScenarioButton:
##                 from mglutil.util.packageFilePath import findFilePath
##                 icondir = findFilePath('32x32', 'scenario.Icons')
##                 PI = Tkinter.PhotoImage
##                 self.animatorIcon = PI(file=os.path.join(icondir, "scenario.gif"))
##                 icon = self.animatorIcon
##                 self.scenarioBt = Tkinter.Button(
##                     frame1, image=icon, width=icon.width(), height=icon.height(),
##                     command=self.startScenario)
##                 self.scenarioBt.pack(side='right', padx=2, pady=2, ipadx=2, ipady=2)
        
        frame1.pack(fill=X, pady=3)
        #
        # Object browser
        #
        self.ObjectList = Frame(self.sframe, borderwidth=2)
        scrollbar = Scrollbar(self.ObjectList, orient=VERTICAL)
    
        #
        # New Tree Object browser ( Class TreeView )
        #
        self.ObjectListFrame = Frame(self.sframe, borderwidth=2, relief='ridge')
        from mglutil.gui.BasicWidgets.Tk.TreeWidget.tree import TreeView
        self.tvolist = TreeView(master=self.ObjectListFrame,displayValue=True)
        self.tvolist.setAction(event='select',
                               function=self.setCurrentObjectFromTree_cb)
        self.tvolist.double1PickNode_Func_cb = self.toggleVisibility_cb
        
        #
        # Reset Normalize Center buttons
        #
        w = self.RNC = Frame(self.sframe, borderwidth=2)
        w1 = Frame(w, borderwidth=2)
        self.resetB = Button(w1, text='Reset')
        self.resetB.pack(side=LEFT, fill=Y, expand=1)
        self.normalizeB = Button(w1, text='Norm.')
        self.normalizeB.pack(side=LEFT, fill=Y, expand=1)
        self.centerB = Button(w1, text='Center')
        self.centerB.pack(side=LEFT, fill=Y, expand=1)
        self.deleteB = Button(w1, text='Delete')
        self.deleteB.pack(side=LEFT, fill=Y, expand=1)
        self.OwnGuiB = Button(w1, text='Settings')
        self.OwnGuiB.pack(side=LEFT, fill=Y, expand=1)
        w1.pack(fill=X)

        #
        # RADIO BUTTONS FOR PROPERTIES SUB PANELS
        #
        self.panel = StringVar()
        self.panel.set('Object')

#        self.PropRadio = Pmw.Group(self.sframe, tag_text='Show properties panel for:')
#        w = self.PropRadio.interior()
#
#        Radiobutton(w, text='Object', var=self.panel,
#                    command = self.PObject, value = 'Object',
#                    indicatoron=0, selectcolor='#c67171',
#                    activebackground='#28686b').pack(side=LEFT)
#        Radiobutton(w, text='Camera', var=self.panel,
#                    command = self.PCamera, value = 'Camera',
#                    indicatoron=0, selectcolor='#c67171',
#                    activebackground='#28686b').pack(side=LEFT)
#        Radiobutton(w, text='Clip', var=self.panel,
#                    command = self.PClip, value = 'Clip',
#                    indicatoron=0, selectcolor='#c67171',
#                    activebackground='#28686b').pack(side=LEFT)
#        Radiobutton(w, text='Light', var=self.panel,
#                    command = self.PLight, value = 'Light',
#                    indicatoron=0, selectcolor='#c67171',
#                    activebackground='#28686b').pack(side=LEFT)

        def propertyRaiseCommand(arg):
            if arg == 'Object':
                self.PObject()
            elif arg == 'Camera':
                self.PCamera()
            elif arg == 'Clip':
                self.PClip()
            elif arg == 'Light':
                self.PLight()
            elif arg == 'Bookmarks':
                self.PViews()
            self.propertyNoteBook.setnaturalsize(pageNames=(arg,))
            
        self.propertyNoteBook = Pmw.NoteBook(self.sframe, raisecommand=propertyRaiseCommand)
        lPropertyPage = {}
        lPropertyPage['Object'] = self.propertyNoteBook.add('Object')
        lPropertyPage['Camera'] = self.propertyNoteBook.add('Camera')
        lPropertyPage['Clip'] = self.propertyNoteBook.add('Clip')
        lPropertyPage['Light'] = self.propertyNoteBook.add('Light')
        lPropertyPage['Bookmarks'] = self.propertyNoteBook.add('Bookmarks')

    
        ##################################
        # OBJECT PROPERTIES FRAME
        ##################################
        self.ObjProp = Frame(lPropertyPage['Object'], relief=RIDGE, borderwidth=2)
    
        # INHERIT PROPERTIES
        self.inheritF = Frame(self.ObjProp)
        
        # create menu button fo inheriting properties
        self.inheritMenu_b = Menubutton(self.inheritF,
                                   text='Current geom properties',
                                   #underline=0,
                                   relief='raise',
                                   direction='right',
                                   )
        self.inheritMenu_b.indices = []

        def toggleMenu(event):
            if event.widget.menu.winfo_ismapped(): event.widget.menu.grab_set()
        self.inheritMenu_b.bind('<Button-1>', toggleMenu, '+')

        self.inheritMenu_b.pack(fill='x')
        
        self.inheritMenu_b.menu = Menu(self.inheritMenu_b)

        majorProp = [
                        'protected',
                        'visible',
                      ]

        moreProp = [
                        #'antialiased',
                        'disableTexture',
                        'faceNormals',
                        'invertNormals',
                        'scissor',
                        'vertexNormals',
                        'transparent',
                      ]

        self.cascadeProp = [
                        'lighting',
                        'sharpColorBoundaries',
                        'stippleLines',
                        'stipplePolygons',
                        ]

        inheritProp = [
                        'inheritMaterial',
                        #'inheritStippleLines',
                        #'inheritStipplePolygons',
                        'inheritXform',
                      ]

        inheritProp2 = [
                        'inheritPointWidth',
                        'inheritLineWidth',
                       ]

        opacityProp = [ 'depthMask',
                      ]

        self.inheritVar = {}
        for prop in majorProp:
            var = Tkinter.IntVar()
            self.inheritVar[prop] = var
            f = CallBackFunction(self.setInherit_cb, prop, var)
            self.inheritMenu_b.menu.add_checkbutton(label=prop,
                                               command=f,
                                               variable=var,)
        
        self.inheritMenu_b.menu.add_separator()
        for prop in inheritProp:
            var = Tkinter.IntVar()
            self.inheritVar[prop] = var
            f = CallBackFunction(self.setInherit_cb, prop, var)
            self.inheritMenu_b.menu.add_checkbutton(label=prop,
                                               command=f,
                                               variable=var,)
            self.inheritMenu_b.indices.append(self.inheritMenu_b.menu.index(prop))
        
        for prop in inheritProp2:
            var = Tkinter.IntVar()
            self.inheritVar[prop] = var
        
        self.relabelingCascadeMenus = {}
        for prop in self.cascadeProp:
                lCascadeMenu = RelabelingCascadeMenu(label=prop,
                                                     variable=Tkinter.IntVar(),
                                                     master=self.inheritMenu_b.menu,
                                                     tearoff=0)
                self.relabelingCascadeMenus[lCascadeMenu.baseLabel] = lCascadeMenu
                self.inheritMenu_b.menu.add_cascade(label=lCascadeMenu.baseLabel,
                                               menu=lCascadeMenu)
                lFuncCall = CallBackFunction(self.setPropMode,
                                             lCascadeMenu.baseLabel,
                                             lCascadeMenu.cascadeVariable)
                lCascadeMenu.add_radiobutton(
                                     label='on',
                                     variable=lCascadeMenu.cascadeVariable,
                                     command=lFuncCall,
                                     value=1,
                                     underline=0)                                                       
                lCascadeMenu.add_radiobutton(
                                     label='off',
                                     variable=lCascadeMenu.cascadeVariable,
                                     command=lFuncCall,
                                     value=0,
                                     underline=0)                                                       
                lCascadeMenu.add_radiobutton(
                                     label='inherit',
                                     variable=lCascadeMenu.cascadeVariable,
                                     command=lFuncCall,
                                     value=viewerConst.INHERIT,
                                     underline=0)
                #lCascadeMenu.setWithoutCallbackFunction('inherit')

        # shading
        lCascadeMenu = RelabelingCascadeMenu(label='shading',
                                             variable=Tkinter.IntVar(),
                                             master=self.inheritMenu_b.menu,
                                             tearoff=0)
        self.relabelingCascadeMenus[lCascadeMenu.baseLabel] = lCascadeMenu
        self.inheritMenu_b.menu.add_cascade(label=lCascadeMenu.baseLabel,
                                       menu=lCascadeMenu)
        lFuncCall = CallBackFunction(self.setPropMode,
                                     lCascadeMenu.baseLabel,
                                     lCascadeMenu.cascadeVariable)
        for lab, mode in viewerConst.SHADINGS.items():
            lCascadeMenu.add_radiobutton(label=lab,
                                         variable=lCascadeMenu.cascadeVariable,
                                         command=lFuncCall,
                                         value=mode,
                                         )

        # culling
        lCascadeMenu = RelabelingCascadeMenu(label='culling',
                                             variable=Tkinter.IntVar(),
                                             master=self.inheritMenu_b.menu,
                                             tearoff=0)
        self.relabelingCascadeMenus[lCascadeMenu.baseLabel] = lCascadeMenu
        self.inheritMenu_b.menu.add_cascade(label=lCascadeMenu.baseLabel,
                                       menu=lCascadeMenu)
        lFuncCall = CallBackFunction(self.setPropMode,
                                     lCascadeMenu.baseLabel,
                                     lCascadeMenu.cascadeVariable)
        for lab, mode in zip(viewerConst.CULLINGS_keys, viewerConst.CULLINGS_values):
            lCascadeMenu.add_radiobutton(label=lab,
                                         variable=lCascadeMenu.cascadeVariable,
                                         command=lFuncCall,
                                         value=mode,
                                         )

        # frontPolyMode
        self.frontPolyMode = Tkinter.IntVar()
        lCascadeMenu = RelabelingCascadeMenu(label='frontPolyMode',
                                             variable=self.frontPolyMode,
                                             master=self.inheritMenu_b.menu,
                                             tearoff=0)
        self.relabelingCascadeMenus[lCascadeMenu.baseLabel] = lCascadeMenu
        self.inheritMenu_b.menu.add_cascade(label=lCascadeMenu.baseLabel,
                                       menu=lCascadeMenu)
        for lab, mode in zip(viewerConst.Front_POLYGON_MODES_keys,
                             viewerConst.Front_POLYGON_MODES_values):
            if not self.viewer.hasOffsetExt:
                continue
            lCascadeMenu.add_radiobutton(label=lab,
                                         variable=lCascadeMenu.cascadeVariable,
                                         command=self.SetFrontPolyMode,
                                         value=mode,
                                         )

        # backPolyMode
        lCascadeMenu = RelabelingCascadeMenu(label='backPolyMode',
                                             variable=Tkinter.IntVar(),
                                             master=self.inheritMenu_b.menu,
                                             tearoff=0)
        self.relabelingCascadeMenus[lCascadeMenu.baseLabel] = lCascadeMenu
        self.inheritMenu_b.menu.add_cascade(label=lCascadeMenu.baseLabel,
                                            menu=lCascadeMenu)
        for lab, mode in zip(viewerConst.Back_POLYGON_MODES_keys,
                             viewerConst.Back_POLYGON_MODES_values):
            if not self.viewer.hasOffsetExt:
                continue
            lCascadeMenu.add_radiobutton(label=lab,
                                         variable=lCascadeMenu.cascadeVariable,
                                         command=self.SetBackPolyMode,
                                         value=mode,
                                         )

        # lineWidth
        lCascadeMenu = RelabelingCascadeMenu(label='lineWidth',
                                             variable=Tkinter.IntVar(),
                                             master=self.inheritMenu_b.menu,
                                             tearoff=0)
        self.relabelingCascadeMenus[lCascadeMenu.baseLabel] = lCascadeMenu
        self.inheritMenu_b.menu.add_cascade(label=lCascadeMenu.baseLabel,
                                            menu=lCascadeMenu)
        def lineWidthFunc():
            self.setPropMode(
                'lineWidth',
                self.relabelingCascadeMenus['lineWidth'].cascadeVariable)
            value = self.relabelingCascadeMenus['lineWidth'].cascadeVariable.get()
            if value != viewerConst.INHERIT:
                self.lw.Set( value, update=False )
                if self.inheritVar['inheritLineWidth'].get() != 0:
                    self.inheritVar['inheritLineWidth'].set(0)
            else:
                self.inheritVar['inheritLineWidth'].set(1)
        for lIndex in range(1, 11):
            lCascadeMenu.add_radiobutton(label=' '+str(lIndex),
                                         variable=lCascadeMenu.cascadeVariable,
                                         command=lineWidthFunc,
                                         value=lIndex,
                                         )
        lCascadeMenu.add_radiobutton(
                                 label='inherit',
                                 variable=lCascadeMenu.cascadeVariable,
                                 command=lineWidthFunc,
                                 value=viewerConst.INHERIT,
                                 )

        # pointWidth
        lCascadeMenu = RelabelingCascadeMenu(label='pointWidth',
                                             variable=Tkinter.IntVar(),
                                             master=self.inheritMenu_b.menu,
                                             tearoff=0)
        self.relabelingCascadeMenus[lCascadeMenu.baseLabel] = lCascadeMenu
        self.inheritMenu_b.menu.add_cascade(label=lCascadeMenu.baseLabel,
                                            menu=lCascadeMenu)
        def pointWidthFunc():
            self.setPropMode(
                'pointWidth',
                self.relabelingCascadeMenus['pointWidth'].cascadeVariable)
            value = self.relabelingCascadeMenus['pointWidth'].cascadeVariable.get()
            if value != viewerConst.INHERIT:
                self.pw.Set( value, update=False )
                self.inheritVar['inheritPointWidth'].set(0)
            else:
                self.inheritVar['inheritPointWidth'].set(1)
        for lIndex in range(1, 11):
            lCascadeMenu.add_radiobutton(label=' '+str(lIndex),
                                         variable=lCascadeMenu.cascadeVariable,
                                         command=pointWidthFunc,
                                         value=lIndex,
                                         )
        lCascadeMenu.add_radiobutton(
                                 label='inherit',
                                 variable=lCascadeMenu.cascadeVariable,
                                 command=pointWidthFunc,
                                 value=viewerConst.INHERIT,
                                 )

        self.inheritMenu_b.menu.add_separator()
        for prop in moreProp:
            var = Tkinter.IntVar()
            self.inheritVar[prop] = var
            f = CallBackFunction(self.setInherit_cb, prop, var)
            self.inheritMenu_b.menu.add_checkbutton(label=prop,
                                               command=f,
                                               variable=var,)

        for prop in opacityProp:
            var = Tkinter.IntVar()
            self.inheritVar[prop] = var
            f = CallBackFunction(self.setInherit_cb, prop, var)
            self.inheritMenu_b.menu.add_checkbutton(label=prop,
                                               command=f,
                                               variable=var,)

        self.opacityRelabelingCascadeMenus = {}
        self.srcBFtk = IntVar()
        lCascadeMenu = RelabelingCascadeMenu(label='srcBF',
                                             variable=self.srcBFtk,
                                             master=self.inheritMenu_b.menu,
                                             tearoff=0)
        self.opacityRelabelingCascadeMenus[lCascadeMenu.baseLabel] = lCascadeMenu
        self.inheritMenu_b.menu.add_cascade(label=lCascadeMenu.baseLabel,
                                       menu=lCascadeMenu)
        for lab, mode in viewerConst.srcBFnames.items():
            lCascadeMenu.add_radiobutton(
                             label=lab,
                             variable=lCascadeMenu.cascadeVariable,
                             command=self.SetBF,
                             value=mode,
                             underline=0)                                                       

        self.dstBFtk = IntVar()
        lCascadeMenu = RelabelingCascadeMenu(label='dstBF',
                                             variable=self.dstBFtk,
                                             master=self.inheritMenu_b.menu,
                                             tearoff=0)
        self.opacityRelabelingCascadeMenus[lCascadeMenu.baseLabel] = lCascadeMenu
        self.inheritMenu_b.menu.add_cascade(label=lCascadeMenu.baseLabel,
                                       menu=lCascadeMenu)
        dstBFnames = [ 'GL.GL_ZERO', 'GL.GL_ONE',
                       'GL.GL_SRC_COLOR', 'GL.GL_ONE_MINUS_SRC_COLOR',
                       'GL.GL_SRC_ALPHA', 'GL.GL_ONE_MINUS_SRC_ALPHA',
                       'GL.GL_DST_ALPHA', 'GL.GL_ONE_MINUS_DST_ALPHA' ]
        for dstBFname in dstBFnames:
            lCascadeMenu.add_radiobutton(
                             label=dstBFname,
                             variable=lCascadeMenu.cascadeVariable,
                             command=self.SetBF,
                             value=eval(dstBFname),
                             underline=0)                                                       


        self.inheritMenu_b['menu'] = self.inheritMenu_b.menu

        otherProp = [ 
                      'inheritBackPolyMode',
                      'inheritCulling',
                      'inheritFrontPolyMode',
                      'inheritLighting',
                      'inheritShading',
                      'inheritSharpColorBoundaries',
                      'inheritStippleLines',
                      'inheritStipplePolygons',                      
                    ]

        self.cascadeProp += ['shading', 'culling', 'lineWidth', 'pointWidth']
        
        lChidrenInheritProp1 = inheritProp \
                             + self.cascadeProp
        lChidrenInheritProp2 = moreProp \
                             + opacityProp 
                             #+ inheritProp2#+ otherProp +#majorProp
        lChidrenInheritProp = lChidrenInheritProp1 + lChidrenInheritProp2
        
        self.cInheritMenu_b = Menubutton(self.inheritF,
                                   text='Propagate property',
                                   relief='raise',
                                   #direction='right',
                                   )
        self.cInheritMenu_b.bind('<Button-1>', toggleMenu, '+')
        self.cInheritMenu_b.pack(fill='x')
        self.cInheritMenu_b.menu = Menu(self.cInheritMenu_b, tearoff=0)
        self.cInheritMenu_b['menu'] = self.cInheritMenu_b.menu

        lCascade = ( Tkinter.Menu(self.cInheritMenu_b.menu, tearoff=1),
                     Tkinter.Menu(self.cInheritMenu_b.menu, tearoff=1) )
        self.cInheritMenu_b.menu.add_cascade(label='to direct children only',
                                             menu=lCascade[0])
        self.cInheritMenu_b.menu.add_cascade(label='to children recursively',
                                             menu=lCascade[1])
        for lRecursive in (0,1):
            for prop in inheritProp:
                f = CallBackFunction(self.setChildProp_cb, prop, lRecursive)
                lCascade[lRecursive].add_command(label=prop, command=f)
            for prop in self.cascadeProp:
                f = CallBackFunction(self.setChildInheritProp_cb, prop, lRecursive)
                lCascade[lRecursive].add_command(label=prop, command=f)
            lCascade[lRecursive].add_separator()
            for prop in lChidrenInheritProp2:
                f = CallBackFunction(self.setChildProp_cb, prop, lRecursive)
                lCascade[lRecursive].add_command(label=prop, command=f)

        # spin settings
        self.spinMenuButton = Button(
                self.inheritF,
                text='Spin settings',
                #underline=0,
                relief='raise',
                command=self.viewer.currentCamera.trackball.toggleSpinGui,
              )
        self.spinMenuButton.pack(fill='x')

        self.inheritF.pack(fill='x')

        # outline mesh properties
        self.outlineMeshButtonInitialized = False
        self.outlineMeshButton = Button(
                self.inheritF,
                text='Outline-Mesh Properties',
                relief='raise',
                command=self.outlineMeshProp_cb,
              )
        self.outlineMeshButton.pack(fill='x')
        
        # strokes properties
        self.strokesButtonInitialized = False
        self.strokesButton = Button(
            self.inheritF,
            text='Strokes',
            relief='raise',
            command=self.strokesProp_cb,
            )
        self.strokesButton.pack(fill='x')
        
        self.objectInheritButtons = {}

        # EDIT OBJECT MATERIAL
        self.objMatEdTk = IntVar()
        self.objMatEdTk.set(0)
        f = Frame(self.ObjProp)
        Label(f, text='Material:').grid(row=0, column=0, sticky='w')
        self.MatEdB1 = Button(f,
                              text='Front',
                              command=self.ObjMatEdFront_cb,
                              )
        self.MatEdB1.grid(row=0, column=2)
        self.MatEdB2 = Button(f,
                              text='Back',
                              command=self.ObjMatEdBack_cb,
                              )
        self.MatEdB2.grid(row=0, column=3)
        #self.inheritVar['inheritMaterial'] = Tkinter.IntVar()
        lFunc = CallBackFunction(self.setInherit_cb,
                                 'inheritMaterial',
                                 self.inheritVar['inheritMaterial'])
        lCheckbutton = Checkbutton(f,
                    text='inherit', 
                    variable=self.inheritVar['inheritMaterial'],
                    command=lFunc,
                    #indicatoron=0,
                    )
        lCheckbutton.grid(row=0, column=4)
        self.objectInheritButtons['inheritMaterial'] = lCheckbutton
        
        # LINE WIDTH
        Label(f, text='Line width:').grid(row=1, column=0, sticky='w')
        self.lw = Slider(f,
                 minval=1, maxval = 10, immediate=1,
                 incr=1, labelformat = '%d', cursortype='int',
                 )
        self.lw.frame.grid(row=1, column=1, columnspan=3)
        def inheritLineWidthFunc():
            if self.inheritVar['inheritLineWidth'].get():
                self.relabelingCascadeMenus['lineWidth'].setWithCallbackFunction(
                    viewerConst.INHERIT)
            else:
                self.relabelingCascadeMenus['lineWidth'].setWithCallbackFunction(
                    self.lw.Get() )
        lCheckbutton = Checkbutton(f,
                    text='inherit', 
                    variable=self.inheritVar['inheritLineWidth'],
                    command=inheritLineWidthFunc,
                    )
        lCheckbutton.grid(row=1, column=4)
        self.objectInheritButtons['inheritLineWidth'] = lCheckbutton
        def CurrentObjectLineWidth(val):
            self.relabelingCascadeMenus['lineWidth'].setWithCallbackFunction(val)
        self.lw.AddCallback(CurrentObjectLineWidth)

        # POINT WIDTH
        Label(f, text='Point width:').grid(row=2, column=0, sticky='w')
        self.pw = Slider(f,
                 minval=1, maxval = 10, immediate=1,
                 incr=1, labelformat = '%d', cursortype='int',
                 )
        self.pw.frame.grid(row=2, column=1, columnspan=3)
        def inheritPointWidthFunc():
            if self.inheritVar['inheritPointWidth'].get():
                self.relabelingCascadeMenus['pointWidth'].setWithCallbackFunction(
                    viewerConst.INHERIT)
            else:
                self.relabelingCascadeMenus['pointWidth'].setWithCallbackFunction(
                    self.pw.Get() )
        lCheckbutton = Checkbutton(f,
                    text='inherit', 
                    variable=self.inheritVar['inheritPointWidth'],
                    command=inheritPointWidthFunc,
                    )
        lCheckbutton.grid(row=2, column=4)
        self.objectInheritButtons['inheritPointWidth'] = lCheckbutton
        def CurrentObjectPointWidth(val):
            self.relabelingCascadeMenus['pointWidth'].setWithCallbackFunction(val)
        self.pw.AddCallback(CurrentObjectPointWidth)

        self.polyModeMenus = []
        # FRONT AND BACK POLYGONS MODE
        #self.polyMode = Frame(self.ObjProp)
        Label(f, text='Polygon mode:').grid(row=3, column=0, sticky='w')
        #self.frontPolyMode = IntVar()
        #self.backPolyMode = IntVar()
        Radio_b1 = Menubutton(f, text='Front',
                                  #underline=0, 
                                  relief='raised',
                                  padx=5,
                                  )
        Radio_b1.grid(row=3, column=1)
        Radio_b1.menu = Menu(Radio_b1)
        self.polyModeMenus.append(Radio_b1)
    
        Radio_b2 = Menubutton(f, text='Back',
                                  #underline=0,
                                  relief='raised',
                                  padx=5,
                                  )
        Radio_b2.grid(row=3, column=2)
        Radio_b2.menu = Menu(Radio_b2)
        self.polyModeMenus.append(Radio_b2)

        labels = viewerConst.Front_POLYGON_MODES_keys[:]
        modes = viewerConst.Front_POLYGON_MODES_values[:]
        for label,mode in zip(labels,modes):
            Radio_b1.menu.add_radiobutton(
                  label=label,
                  var=self.relabelingCascadeMenus['frontPolyMode'].cascadeVariable,
                  value=mode,
                  command=self.relabelingCascadeMenus['frontPolyMode'].setWithCallbackFunction)
        Radio_b1['menu'] = Radio_b1.menu
        
        labels = viewerConst.Back_POLYGON_MODES_keys[:]
        modes = viewerConst.Back_POLYGON_MODES_values[:]
        for label,mode in zip(labels,modes):
            Radio_b2.menu.add_radiobutton(
                  label=label,
                  var=self.relabelingCascadeMenus['backPolyMode'].cascadeVariable,
                  value=mode,
                  command=self.relabelingCascadeMenus['backPolyMode'].setWithCallbackFunction)
        Radio_b2['menu'] = Radio_b2.menu

        # Culling
        Radio_b3 = Menubutton(f, text='culling', relief='raised')
        Radio_b3.grid(row=3, column=3)
        Radio_b3.menu = Menu(Radio_b3)
        for lab, mode in zip(viewerConst.CULLINGS_keys, viewerConst.CULLINGS_values):
            Radio_b3.menu.add_radiobutton(
                  label=lab,
                  var=self.relabelingCascadeMenus['culling'].cascadeVariable,
                  value=mode,
                  command=self.relabelingCascadeMenus['culling'].setWithCallbackFunction,
                  )
        Radio_b3['menu'] = Radio_b3.menu
        self.polyModeMenus.append(Radio_b3)
        

        #self.zsortFrame = Frame(self.ObjProp)
        b = Label(f, text='Transparency order:')
        b.grid(row=4, column=0, columnspan=2, sticky='w')
        # sort polygons along z direction
        b = Button(f,
                   text='Zsort',
                   command=self.zSortPoly)
        b.grid(row=4, column=2)
    
        # sort polygons along -z direction
        b = Button(f,
                   text='-Zsort',
                   command=self.minuszSortPoly)
        b.grid(row=4, column=3)
            
        f.pack(fill='x', expand=1)
    


    #         b = Button(self.ObjProp, text='Read Object', command=self.read)
    #         b.pack(fill=X, expand=1)
    #         b = Button(self.ObjProp, text='Write Object', command=self.write)
    #         b.pack(fill=X, expand=1)
    #         b = Button(self.ObjProp, text='Edit Prop.', command=self.ObjectProp)
    #         b.pack(fill=X, expand=1)

        ##################################
        # CAMERA PROPERTIES PANEL
        ##################################
        self.CameraProp = Frame(lPropertyPage['Camera'], relief=RIDGE, borderwidth=2)
    
        self.drawSceneBB = IntVar()
        Radiobutton_button = Menubutton(self.CameraProp,
                        text='Bounding Box',
                        #underline=0,
                                                relief='raised')       
        Radiobutton_button.bind('<Button-1>', toggleMenu, '+')
        Radiobutton_button.pack(fill=X, padx=5)
        # the primary pulldown
        bbmodes = ('NO', 'ONLY', 'WITHOBJECT')
        Radiobutton_button.menu = Menu(Radiobutton_button)
        for v in (0,1,2):
            Radiobutton_button.menu.add_radiobutton(label=bbmodes[v],
                          var=self.drawSceneBB,
                          value = viewerConst.BB_MODES[v],
                          command = self.DrawSceneBB_cb)

        Radiobutton_button['menu'] = Radiobutton_button.menu

        # COLOR CHOOSER

        def setSelectionContourColor_cb(color):
            DejaVu2.selectionContourColor = color
            self.selectionGUI.contourColorButton.configure(background=TkColor(DejaVu2.selectionContourColor))
            self.viewer.Redraw()

        lTargetDict = {
                       'ambient light':   
                           (self.viewer.lightModel.ambient[:3],
                            'RGB',
                            self.viewer.LMColor
                           ),
                       'background':
                           (self.viewer.currentCamera.backgroundColor[:3],
                            'RGB',
                            self.viewer.CurrentCameraBackgroundColor
                           ),
                       'selection contour':
                           (DejaVu2.selectionContourColor,
                            'RGB',
                            setSelectionContourColor_cb
                           ),
                      }
        for i in range(maxClip):
            lTargetText = 'clip %d'%(i+1)
            lTargetDict[lTargetText] = (self.viewer.clipP[i].color[:3],
                                        'RGB',
                                        self.viewer.clipP[i].setColor,
                                       )
        for i in range(maxLights):
            ii = i + 1
            lTargetText = 'light %d - ambient'%ii
            lTargetDict[lTargetText] = (self.viewer.lights[i].ambient[:3],
                                        'RGB',
                                        self.viewer.lights[i].setAmbient,
                                       )
            lTargetText = 'light %d - diffuse'%ii
            lTargetDict[lTargetText] = (self.viewer.lights[i].diffuse[:3],
                                        'RGB',
                                        self.viewer.lights[i].setDiffuse,
                                       )
            lTargetText = 'light %d - specular'%ii
            lTargetDict[lTargetText] = (self.viewer.lights[i].specular[:3],
                                        'RGB',
                                        self.viewer.lights[i].setSpecular,
                                       )
        lTopColorChooser = Tkinter.Toplevel() 
        self.colorChooser = ColorChooser(master=lTopColorChooser, 
                                         targetDict=lTargetDict,
                                         targetKey='background',
                                        )
        lTopColorChooser.withdraw()

        # BACKGROUND COLOR
        lFunc = CallBackFunction(self.colorChooser.showColorChooser,'background')
        b = Button(self.CameraProp,
                   text='Background Color',
                   command=lFunc,
                  )
        b.pack(fill=X, padx=5)

        # AUTODEPTHCUE
        b = Button(self.CameraProp, text='Auto Depthcue',
                   command=self.AutoDepthcue).pack(fill=X, padx=5)

        # add a button for video recorder
        # check if the camera is recordable:
        try:
            from DejaVu2.Camera import RecordableCamera
            isrecordable = True
        except:
            isrecordable = False
        if isrecordable:
            camera = self.viewer.currentCamera
            if isinstance(camera, RecordableCamera):
                b =  Button(self.CameraProp, text='Video Recorder',
                            command = self.CameraVideoRecorder_cb)
                b.pack(fill=X, padx=5)
                if not hasattr(camera, "videoRecorder"):
                    camera.videoRecorder = None
    
        # PROJECTION TYPE
        self.projType = IntVar()
        Radiobutton_button = Menubutton(self.CameraProp,
                        text='Projection',
                        underline=0, relief='raised')
        Radiobutton_button.bind('<Button-1>', toggleMenu, '+')
        Radiobutton_button.pack(fill=X, padx=5)
        # the primary pulldown
        Radiobutton_button.menu = Menu(Radiobutton_button)
        for v,s in {0:'Perspective', 1:'Orthographic'}.items():
            Radiobutton_button.menu.add_radiobutton(label=s,
                                var=self.projType,
                                value = v,
                            command = self.SetProjection)
        self.projType.set(self.viewer.currentCamera.projectionType)
        Radiobutton_button['menu'] = Radiobutton_button.menu
    
        # ANTIALIASED
        self.nbJitter = IntVar()
        self.nbJitter.set(self.viewer.currentCamera.antiAliased)
        Radiobutton_button = Menubutton(self.CameraProp,
                        text='Scene Antialiasing',
                        underline=0, relief='raised')
        Radiobutton_button.bind('<Button-1>', toggleMenu, '+')
        Radiobutton_button.pack(fill=X, padx=5)
        # the primary pulldown
        Radiobutton_button.menu = Menu(Radiobutton_button)
        for v in jitter.jitterList:
            Radiobutton_button.menu.add_radiobutton(label=str(v),
                                var=self.nbJitter,
                                value = v,
                                command = self.SetJitter)
    
        Radiobutton_button['menu'] = Radiobutton_button.menu

        #NPR OUTLINES
        Checkbutton_button = Menubutton(self.CameraProp,
                        text='Cartoon Outlines',
                        underline=0, relief='raised')
        Checkbutton_button.bind('<Button-1>', toggleMenu, '+')
        Checkbutton_button.pack(fill=X, padx=5)
        
        self.contourTk = viewer.contourTk
        Checkbutton_button.menu=Menu(Checkbutton_button)
                      
        # CONTOUR
        Checkbutton_button.menu.add_checkbutton(label='On',
                        var = self.contourTk,
                        command=self.toggleOutline)
        #SET NPR PARAMETERS
        Checkbutton_button.menu.add_command(
                        label='Edit Parameters',
                        command=self.showCurveTool)

        Checkbutton_button['menu'] = Checkbutton_button.menu

        #GraphToolWidget
        from mglutil.gui.BasicWidgets.Tk.graphtool import GraphApp
        self.GraphToolpanel=Tkinter.Toplevel()
        self.GraphToolpanel.title("GraphTool")
        self.GraphToolpanel.withdraw()
        self.GraphToolpanel.protocol('WM_DELETE_WINDOW', self.GraphToolpanel.withdraw)
        f1 = Tkinter.Frame(self.GraphToolpanel)
        f1.pack(side='top', fill='both', expand=1)
        self.curvetool=GraphApp(f1,callback=self.continuousRamp)
        self.curvetool.defaultcurve_cb() 
        self.d1scalewheel=self.curvetool.d1scalewheel
        self.d1scalewheel.set(self.viewer.currentCamera.d1scale)
        self.d1scalewheel.callbacks.AddCallback(self.setNPR_cb)


#===============================================================================
# Screen Space Ambient Occlusion
#===============================================================================
        Checkbutton_button = Menubutton(self.CameraProp,
                        text='Screen Space Ambient Occlusion',
                        underline=0, relief='raised')
        Checkbutton_button.bind('<Button-1>', toggleMenu, '+')
        Checkbutton_button.pack(fill=X, padx=5)
        
        self.ssaoTk = BooleanVar()
        self.ssaoTk.set(self.viewer.currentCamera.ssao)
        Checkbutton_button.menu=Menu(Checkbutton_button)
                      
        # TOGGLE
        Checkbutton_button.menu.add_checkbutton(label='On',
                        var = self.ssaoTk,
                        command=self.toggleSSAO)
        #SET SSAO PARAMETERS
        Checkbutton_button.menu.add_command(
                        label='Edit Parameters',
                        command=self.showSSAOTool)
        Checkbutton_button['menu'] = Checkbutton_button.menu


        # occlusion cam GUI
        ## b = Button(self.CameraProp,
        ##            text='Ambient Occlusion',
        ##            underline=0,
        ##            command=self.showOcclusionCamGUI)
        ## b.pack(fill=X, padx=5)

        # SelectionGui
        self.selectionGUI = SelectionGUI(viewerGui=self)
        b = Button(self.CameraProp,
                   text='Selection Settings',
                   underline=0,
                   command=self.selectionGUI.toggle)
        b.pack(fill=X, padx=5)

        # DEPTHCUE
        self.depthcued = IntVar()
        b = Checkbutton(self.CameraProp, text='Depthcueing',
                variable=self.depthcued,
                command=self.viewer.ToggleDepth,
                onvalue=viewerConst.YES,
                offvalue=viewerConst.NO).pack(fill=X, padx=5)
        self.depthcued.set(self.viewer.depthcued)

        # THUMBNAIL
        self.drawThumbnail = IntVar()
        b = Checkbutton(self.CameraProp, text='Thumbnail',
                        variable=self.drawThumbnail,
                        command=self.ToggleThumbnail,
                        onvalue=viewerConst.YES,
                        offvalue=viewerConst.NO).pack(fill=X, padx=5)

        ##SET OVER ALL LIGHTING
        b=Checkbutton(self.CameraProp,text='Overall Lighting',
                      variable = self.viewer.OverAllLightingIsOn,
                      command = self.viewer.deleteOpenglListAndCallRedrawAndCallDisableGlLighting )
        b.pack(fill=X, padx=5)



        # NEAR/FAR DEPTCUEING SLIDER
        self.NearFarFog = TwoSlidersGUI(viewer, self.CameraProp)
    
        #################################
        # Light properties frame
        ##################################
    
        self.LightProp = Frame(lPropertyPage['Light'], relief=RIDGE, borderwidth=2)
    
        ###################################
        #LIGHT MODEL
        lmframe = Frame(self.LightProp, relief=RIDGE, borderwidth=2)
    
        self.localViewer = IntVar()
        self.localViewer.set(self.viewer.lightModel.localViewer is True)
        b = Checkbutton(lmframe, text='Local Viewer', width=10,
                variable=self.localViewer,
                command=self.ToggleLocalViewer,
                onvalue=viewerConst.YES,
                offvalue=viewerConst.NO).grid(row=0, column=0)
    
        self.twoSide = IntVar()
        self.twoSide.set(self.viewer.lightModel.twoSide is True)
        b = Checkbutton(lmframe, text='Two Side', width=7,
                variable=self.twoSide,
                command=self.ToggleTwoSide,
                onvalue=viewerConst.YES,
                offvalue=viewerConst.NO).grid(row=0, column=1)

        lFunc = CallBackFunction(self.colorChooser.showColorChooser,'light 1 - diffuse')
        b = Button(lmframe,
                   text='Light Colors',
                   command=lFunc,
                  )
        b.grid(row=1, column=0, columnspan=4)

        #self.lightModelColor = IntVar()   
    
        ###################################
        #LIGHT Sources
    
        lsframe = Frame(self.LightProp, relief=RIDGE, borderwidth=2)

        self.CurrentLight = IntVar()
        self.CurrentLight.set( 1 )
        f1 = Frame(lsframe, borderwidth=0)
        self.lightSourceOnOffCheckButtons = []
        b = Checkbutton(f1, 
                        text='1 \'key\'',
                        variable=self.CurrentLight,
                        onvalue=1,
                        offvalue=0,
                        command=self.LightSelect)
        self.lightSourceOnOffCheckButtons.append(b)
        b.grid(row=0, column=0)
        b = Checkbutton(f1, 
                        text='2 \'fill\'',
                        variable=self.CurrentLight,
                        onvalue=2,
                        offvalue=0,
                        command=self.LightSelect)
        self.lightSourceOnOffCheckButtons.append(b)
        b.grid(row=0, column=1)
        b = Checkbutton(f1, 
                        text='3 \'reflective\'',
                        variable=self.CurrentLight,
                        onvalue=3,
                        offvalue=0,
                        command=self.LightSelect)
        self.lightSourceOnOffCheckButtons.append(b)
        b.grid(row=0, column=2)
        f12 = Frame(lsframe, borderwidth=0)
        for i in range(3, maxLights):
            b = Checkbutton(f12,
                            text=str(i+1),
                            variable=self.CurrentLight,
                            onvalue=i+1,
                            offvalue=0,
                            command=self.LightSelect)
            self.lightSourceOnOffCheckButtons.append(b)
            b.grid(row=1, column=i-3)

        self.lightOnOff = IntVar()
        self.lightOnOff.set(1)
        f2 = Frame(lsframe, borderwidth=0)
        self.lightOnOffButton = Checkbutton(f2, text='Light On',
                                                variable = self.lightOnOff,
                                                command = self.LightOnOff,
                                                width=14,
                                                onvalue=viewerConst.YES,
                                                offvalue=viewerConst.NO)
            
        self.lightOnOffButton.pack()
    
        self.showLight = IntVar()
        b = Checkbutton(f2, text='Show Lights', command = self.LightShow,
                variable=self.showLight, width=14,
                onvalue=viewerConst.YES,
                offvalue=viewerConst.NO)
        b.pack()
    
        self.lightColor = IntVar()

    #        f3 = Frame(lsframe, borderwidth=2)
    #         self.LightType = StringVar()
    #         self.LightType.set('Directional')
    #         Radiobutton(f3, text="Directional",
    #                     command = self.LightDirectional, value = 'Directional',
    #                     variable = self.LightType,
    #                     indicatoron=0, selectcolor="#c67171",
    #                     activebackground="#28686b").pack(fill=X)
    #         Radiobutton(f3, text="Positional",
    #                     command = self.LightPositional, value = 'Positional',
    #                     variable = self.LightType,
    #                     indicatoron=0, selectcolor="#c67171",
    #                     activebackground="#28686b").pack(fill=X)
    #         Radiobutton(f3, text="Spot",
    #                     command = self.LightSpot, value = 'Spot',
    #                     variable = self.LightType,
    #                     indicatoron=0, selectcolor="#c67171",
    #                     activebackground="#28686b").pack(fill=X)

        f1.pack(fill=X, pady = 3, expand = 1)
        f12.pack(fill=X, pady = 3, expand = 1)
        f2.pack(fill=X, pady = 3, expand = 1)
        lmframe.pack(fill=X, pady = 3, expand = 1)
        lsframe.pack(fill=X, pady = 3, expand = 1)
    #        f3.pack(fill=X, expand = 1)

        #################################
        # CLIPPING PLANES PANEL
        ##################################
    
        self.ClipProp = Frame(lPropertyPage['Clip'], relief=RIDGE, borderwidth=2)
    
        self.CurrentClip = IntVar()
        self.CurrentClip.set( 1 )
        f1 = Frame(self.ClipProp, borderwidth=0)
    
        lab = ('on', 'side', 'clip\nchildren', 'display', 
               'capGL', 'capMesh', 'current')
        lenlab = len(lab)
        lenlabm1 = lenlab - 1
        for i in range(lenlab):
            Label(f1, text=lab[i]).grid(row=0, column=i+1)

        self.clipvar = [ ]
        self.clipw = [ ]
        funcs = (self.ClipOnOff, self.ClipSide,
                 self.ClipInherit, self.ClipVisible,
                 self.ClipCap, self.ClipCapMesh)

        for i in range(1, maxClip+1):
            l = []
            self.clipvar.append(l)
            lw = []
            self.clipw.append(lw)
            w = Label(f1, text=str(i))
            w.grid(row=i+1, column=0)
            for j in range(lenlabm1):
                v = IntVar()
                l.append(v)
                w = Checkbutton(f1, variable=v, onvalue=i,
                        offvalue=0, state=DISABLED,
                        command=funcs[j])
                w.grid(row=i+1, column=j+1)
                lw.append(w)
            Radiobutton(f1, variable=self.CurrentClip, value=i,
                        command=self.ClipSelect).grid(row=i+1, column=lenlab)

        for i in range(lenlabm1):
            self.clipw[0][i].configure(state=NORMAL)

        f2 = Frame(self.ClipProp, borderwidth=0)

        lFunc = CallBackFunction(self.colorChooser.showColorChooser,'clip 1')
        b = Button(f2,
                   text='Clip plane colors',
                   command=lFunc,
                  )
        b.pack()

#        self.clipColor = IntVar()
#        self.editClipCol = Checkbutton(f2, text='Edit ClipPlane Color',
#                    width=16,
#                variable=self.clipColor,
#                command = self.ClipColorButton_cb,
#                onvalue=viewerConst.YES,
#                offvalue=viewerConst.NO)
#        self.editClipCol.pack()

        f1.pack(fill=X, pady = 3, expand = 1)
        f2.pack(fill=X, pady = 3, expand = 1)

        self.top.pack(fill=BOTH)
        self.XformRadio.pack(fill='x', padx=2, pady=2, ipadx=2, ipady=2)
        #self.ObjectList.pack(fill=BOTH, expand=1,padx=2, pady=2)
        self.ObjectListFrame.pack(fill=BOTH, expand=1,padx=2, pady=2)
        self.RNC.pack(padx=2, pady=2, ipadx=2, ipady=2)
    

        self.ObjProp.pack(fill='x', ipadx=2, ipady=2, padx=2, pady=2)
        self.propertyNoteBook.setnaturalsize(pageNames=lPropertyPage.keys())
        self.propertyNoteBook.pack(padx=2, pady=2, ipadx=2, ipady=2, fill='x')
        self.CurrentPropPanel = self.ObjProp
            #################################
            # OBJECT MATERIAL PANNEL
        ##################################
            # NEW TOPLEVEL
    ##          self.colorRoot = Toplevel()
    ##          self.colorRoot.title("Object Material")
    ##          self.colorTop=Frame(self.colorRoot)

        # MATERIAL EDITOR
#        self.materialEditor = PropertyEditor.MaterialEditor(
#                                    self.sframe,
#                                    self.colorChooser,
#                                    )
        # LIGHT COLOR EDITOR
#        self.lightColorEditor = PropertyEditor.LightColorEditor(
#                                        self.sframe,
#                                        self.colorChooser,
#                                        )

        self.sf.interior().pack(fill = 'both', expand = 1)
        self.sf.reposition()
        self.sf.pack(fill = 'both', expand = 1)

        ######################################################################
        ###################### Bookmarks pane ################################
        ######################################################################
        self.BookmarksProp = Frame(lPropertyPage['Bookmarks'], relief=RIDGE,
                                   borderwidth=2)
        self.createBookmarkPane(self.BookmarksProp)


    def strokesProp_cb(self, event=None):
        self.viewer.currentObject.applyStrokes(strength=2) # modified Stefano
        

    def _setOMP(self, attr, obj, val=None):
        if attr=='lighting':
            val = self.lightOutlineMesh.get()
        elif attr=='color':
            obj._setTransparent('implicit')
            
        d = {attr:val}
        #print 'setting outline', obj.fullName, attr, val
                
        obj.outline.Set(**d)
        self.viewer.Redraw()

        
    def outlineMeshProp_cb(self, event=None, geometry=None):
        
        if geometry is None:
            obj = self.viewer.currentObject
        else:
            obj = geometry
            
        if not hasattr(obj, 'outline'):
            return

        if not self.outlineMeshButtonInitialized:
                
            frame = self.outlineMeshPanel = Tkinter.Toplevel()
            frame.geometry("+%d+%d"%frame.winfo_pointerxy())
            frame.protocol('WM_DELETE_WINDOW', frame.withdraw)
          
            lFuncCall = CallBackFunction(self._setOMP, 'factor', obj)
            self.ompfactortw = ThumbWheel(
                frame, width=70, height=16, type=float,
                value=obj.outline.factor, callback=lFuncCall,
                continuous=True, oneTurn=10., wheelPad=2, min=0.1,
                labCfg = {'text': 'factor:', 'side':'top'})
            self.ompfactortw.grid(column=0, row=0)

            lFuncCall = CallBackFunction(self._setOMP, 'unit', obj)
            self.ompunittw = ThumbWheel(
                frame, width=70, height=16, type=float,
                value=obj.outline.unit, callback=lFuncCall,
                continuous=True, oneTurn=10., wheelPad=2, min=0.1,
                labCfg = {'text': 'unit:', 'side':'top'})
            self.ompunittw.grid(column=1, row=0)

            lFuncCall = CallBackFunction(self._setOMP, 'lineWidth', obj)
            self.omplwtw = ThumbWheel(
                frame, width=70, height=16, type=int,
                value=obj.outline.lineWidth, callback=lFuncCall,
                continuous=True, oneTurn=10., wheelPad=2, min=1,
                labCfg = {'text':'line width:', 'side':'top'})
            self.omplwtw.grid(column=2, row=0)

            lFuncCall = CallBackFunction(self._setOMP, 'lighting', obj)
            self.lightOutlineMesh = IntVar()
            self.omplightingb = Checkbutton(
                frame, text='lighting', 
                    variable=self.lightOutlineMesh,
                    command=lFuncCall,
                )
            self.omplightingb.grid(column=3, row=0)

            cc = self.ompcolorchooser = ColorChooser(
                frame, gridCfg={'column':0, 'row':1, 'columnspan':4})
            lFuncCall = CallBackFunction(self._setOMP, 'color', obj)
            cc.Set(obj.outline.color, 'RGB', run=False)
            cc.AddCallback(lFuncCall)

            frame.title('outline-mesh properties for '+obj.fullName)
            self.outlineMeshButtonInitialized = True
        else:
            frame = self.outlineMeshPanel
            frame.deiconify()
            frame.geometry("+%d+%d"%frame.winfo_pointerxy())
            self.updateOMPgui(obj)

    def updateOMPgui(self, obj=None):
        if self.outlineMeshButtonInitialized:
            if obj is None:
                obj = self.viewer.currentObject
            if self.outlineMeshPanel.winfo_ismapped():
                if not hasattr(obj, 'outline'):
                    self.outlineMeshPanel.withdraw()
                    return
                else:
                    self.outlineMeshPanel.title('outline-mesh properties for '+obj.fullName)
                    self.ompfactortw.set(obj.outline.factor, update=0)
                    self.ompunittw.set(obj.outline.unit, update=0)
                    self.omplwtw.set(obj.outline.lineWidth, update=0)
                    self.lightOutlineMesh.set(obj.outline.lighting)
                    self.ompcolorchooser.Set(obj.outline.color, run=False)
            
            
    def createBookmarkPane(self, master):
        self.masterFrame = master

        self.allRepr = 0  # repr counter used to name representations by default
        self.reprMem = {}

	w = self.reprContainer = Pmw.Group(master,
                                           tag_pyclass = Tkinter.Button,
                                           tag_text='Save representation')
        w.configure(tag_command = self.addRepr)
	w.pack(fill = 'both', expand = 1, padx = 6, pady = 6)

        self.viewsBalloon = Pmw.Balloon(self.root)
        self.reprMenu = Tkinter.Menu(self.root, title = "Representation")
        self.reprMenu.add_command(label="Rename repr")
        self.reprMenu.add_command(label="Remove repr")
        self.reprMenu.add_command(label="Dismiss")

        
    def addOrient(self, event=None):
        
        name = 'orient%d'% self.allOrients
        self.allOrients = self.allOrients + 1
        c = self.viewer.currentCamera
        mini, maxi = self.viewer.rootObject.ComputeBB(self.viewer.currentCamera)
        lBox = maxi - mini
        lHalfObject = max(lBox)/2.
        if lHalfObject == 0.:
            lHalfObject = 1.
            
        import math
        from scenario.interpolators import matToQuaternion, quatToMatrix

        dist = lHalfObject / math.tan(c.fovyNeutral/2*math.pi/180.0)
        lookFrom = c.nearDefault+dist+lHalfObject
        
        self.orients[name] = {
            'quat' : matToQuaternion(self.viewer.rootObject.rotation),
            'trans' : self.viewer.rootObject.translation[:],
            'scale' : self.viewer.rootObject.scale[:],
            'pivot' : self.viewer.rootObject.pivot,
            'fovy' : c.fovy,
            'lfrom' : lookFrom,
            'near': c.near,
            'far': c.far,
            'start':c.fog.start,
            'end': c.fog.end
            }
        #positionNames = positions.keys()
        #positionNames.sort()
        #positionList.setlist(positionNames)
        #self.addOrientScenario(name)
        orient = self.orients[name]
        self.orientsDirectors[name] = self.createOrientScenario(name, orient)
        self.addOrientButton(name)
        

    def getButtonIcon(self):
        import Image, ImageChops, ImageTk
        c = self.viewer.currentCamera
        im = c.GrabFrontBuffer()
        def autocrop(im, bgcolor):
            if im.mode != "RGB":
                im = im.convert("RGB")
            bg = Image.new("RGB", im.size, bgcolor)
            diff = ImageChops.difference(im, bg)
            bbox = diff.getbbox()
            if bbox:
                return im.crop(bbox)
            return None # no contents

        imc = autocrop(im,
                tuple([int(round(x*255)) for x in c.backgroundColor[:3]]))
        if imc is None:
            imc = im
        ims = imc.resize((50, 50), Image.ANTIALIAS)
        return ImageTk.PhotoImage(ims), ims


    def saveRepr(self):
        statesMem = {}
        for g in self.viewer.rootObject.AllObjects():
            if not g.visible:
                state = {'visible':0}
            else:
                state = g.getState(full=1)
                del state['rotation']
                del state['translation']
                del state['scale']
                del state['pivot']
                del state['name']
            statesMem[g.fullName] = state
        return statesMem


    def restoreRepr(self, name, button = None):
        mem = self.reprMem[name]
        for g in self.viewer.rootObject.AllObjects():
            if not mem.has_key(g.fullName):
                g.visible = 0
            else:
                #print g.fullName
                tmp = mem[g.fullName].copy()
                mat = tmp.pop('rawMaterialF', None)
                if mat: g.materials[GL.GL_FRONT].Set(**mat)
                mat = tmp.pop('rawMaterialB', None)
                if mat: g.materials[GL.GL_BACK].Set(**mat)
                g.Set( **tmp )
        if button: 
            if not hasattr(button, "photo"):
                #print "updating image for button %s" % name
                self.viewer.master.master.lift()
                self.viewer.OneRedraw()
                photo, ims = self.getButtonIcon()
                button.configure(compound='left', image = photo, text = "",
                                 width=photo.width(), height = photo.height())
                button.photo = photo

                
    def removeRepr(self, name):
        reprB = None
        frame = self.reprContainer.interior()
        for b in frame.grid_slaves():
            if hasattr(b, "name"):
                if b.name == name:
                    reprB = b
                    break
        if reprB:
            reprB.destroy()
            buttons = frame.grid_slaves()
            buttons.reverse()
            for i, b in enumerate(buttons):
                b.grid(row=1+ i/5, column= i%5, sticky='w')
            self.propertyNoteBook.setnaturalsize(pageNames=('Bookmarks',))
            self.reprMem.pop(name)

        
    def renameRepr(self, name):
        from tkSimpleDialog import askstring
        newname = askstring("Rename repr %s"%name, "Enter new name:",
                            initialvalue = name,
                            parent = self.reprContainer.interior())
        if newname != None and newname != name:
            if self.reprMem.has_key(newname):
                from tkMessageBox  import showwarning
                showwarning("Warning", "Representation name %s already exists"%newname, parent = self.root)
                return
            #find cooresponding button, rename it and update the bindings:
            for b in self.reprContainer.interior().grid_slaves():
                if hasattr(b, "name"):
                    if b.name == name:
                       b.name = newname
                       b.configure(command = CallBackFunction( self.restoreRepr, newname))
                       b.bind('<Button-3>', CallBackFunction( self.showReprMenu_cb, newname))
                       self.viewsBalloon.bind(b, newname)
                       repr = self.reprMem.pop(name)
                       self.reprMem[newname] = repr
                       break 

        
    def showReprMenu_cb(self, name, event = None):
        entry = self.reprMenu.index("Rename repr")
        self.reprMenu.entryconfigure(
            entry, command = CallBackFunction(self.renameRepr, name))

        entry = self.reprMenu.index("Remove repr")
        self.reprMenu.entryconfigure(
            entry, command = CallBackFunction(self.removeRepr, name))

        self.reprMenu.post(event.x_root, event.y_root)

        
    def addRepr(self, event=None):
        
        name = 'repr%d'% self.allRepr
        self.allRepr += 1
        self.reprMem[name] = self.saveRepr()
        self.addReprButton(name)
        

    def addReprButton(self, name, icon = True):
        if self.allRepr == 0:
            self.viewsBalloon.bind(self.reprContainer._ring, "Right click on image\nto display its menu")
        master = self.reprContainer.interior()
        if icon == True:
            photo, ims = self.getButtonIcon()
            b = Tkinter.Button(master=master ,compound='left',
                               image=photo)
            b.photo = photo
        else:
            b = Tkinter.Button(master=master , text=name, width = 4, height = 3, wraplength=50)
        self.viewsBalloon.bind(b, name)
        b.name = name
        b.configure(command = CallBackFunction( self.restoreRepr, name, b))
        b.grid(row=1+(len(self.reprMem)-1)/5, column=(len(self.reprMem)-1)%5, sticky='w')
        self.propertyNoteBook.setnaturalsize(pageNames=('Bookmarks',))
        self.sframe.pack()
        b.bind('<Button-3>', CallBackFunction( self.showReprMenu_cb, name))


    def saveReprToFile(self, file = None):
        if len(self.reprMem) == 0:
            from tkMessageBox import showwarning
            showwarning("Warning", "There are no representations in Bookmarks to save", parent=self.root)
            return
        if file is None:
            oldtypes = self.saveFileBrowser.filetypes
            self.saveFileBrowser.filetypes = [('', '*_repr.db'), ('all', '*')]
            file = self.getSaveFileName()
            self.saveFileBrowser.filetypes = oldtypes 
        if file:
            import shelve
            f = shelve.open(file)
            if len(f.keys()):
                f.clear()
            for name in self.reprMem.keys():
                f[name] = self.reprMem[name]
            f.close()


    def loadReprFromFile(self, file = None):
        if file is None:
            oldtypes = self.openFileBrowser.filetypes
            self.openFileBrowser.filetypes = [('', '*_repr.db'), ('all', '*.py')]
            file = self.getLoadFileName()
            self.openFileBrowser.filetypes = oldtypes
        if file:
            import shelve
            newrepr = shelve.open(file)
            action = "ADD"
            if len(self.reprMem):
                # ask the user if he wants to add or replace the representations:
                dialog = Pmw.Dialog(self.root, title = "Representations",
                                    buttons = ('ADD', 'REPLACE', 'Cancel'),
                                    defaultbutton = 'ADD')
                dialog.withdraw()
                label = Tkinter.Label(dialog.interior(),
                    text = "Would you like to ADD file representaions\n or REPLACE existing ones?", pady = 10) 
                label.pack(expand = 1, fill = 'both', padx = 4, pady = 4)
                action = dialog.activate(geometry = '+%d+%d' % (self.root.winfo_x()+100, self.root.winfo_y()+200))
                if action == "Cancel": return
                elif action == "REPLACE":
                    # remove all existing representations
                    for b in self.reprContainer.interior().grid_slaves():
                        b.destroy()
                    self.reprMem = {}
            nrepr = self.allRepr
            self.viewer.master.master.lift()
            #add new repr
            for i, name in enumerate(newrepr.keys()):
                newname = name
                # if saved representation has "custom" name,
                # check if  that name already exists in Representations:
                if not name.startswith("repr"):
                    if self.reprMem.has_key(name):
                        # name exists -> create newname by adding an
                        # integer to the name 
                        count = 0
                        for k in self.reprMem.keys():
                            if k.startswith(name):
                                count = count +1
                        newname = name+"_"+str(count)
                else:
                    newname = 'repr%d'% (nrepr + i)
                self.reprMem[newname] = newrepr[name]
                # ?how can we create an image button of the saved repr?
                # display that view and "take picture" of it? ->
                #self.restoreRepr(newname) #This could be very slow
                self.addReprButton(newname, icon = False)
                self.allRepr = self.allRepr + 1
            newrepr.close()
        

    def bindModifersToTransformInfo(self, master):
        #print "bindModifersToTransformInfo", master
        if master is not None:
            lKeys = ('Alt_L','Alt_R','Control_L','Control_R','Meta_L','Meta_R','Shift_L','Shift_R')
            for lKey in lKeys:
                master.bind('<KeyPress-%s>'%lKey, self.fillTransformInfo_cb)
                master.bind('<KeyRelease-%s>'%lKey, self.fillTransformInfo_cb)
            if hasattr(master, 'master'):
                self.bindModifersToTransformInfo(master.master)


    def fillTransformInfo_cb(self, event=None):
        #print "fillTransformInfo_cb", event, dir(event)

        if isinstance(event, Tkinter.Event) and (event.type == '2'):
            lEvent = event.keysym
        elif type(event) == types.StringType:
            lEvent = event
        else:
            lEvent = 'None'

        if lEvent.startswith('Alt'):
            lModifier = 'Alt'
        elif lEvent.startswith('Shift'):
            lModifier = 'Shift'
        elif lEvent.startswith('Control'):
            lModifier = 'Control'
        elif lEvent.startswith('Meta'):
            lModifier = 'Meta'
        else:
            lModifier = 'None'

        for lButtonIndex in (1, 2, 3):
           lText=self.viewer.currentCamera.mouseButtonActions[self.Xform.get()][lButtonIndex][lModifier]
           if lText == 'None':
               lText = ''
           self.transformInfo[lButtonIndex].configure(text=lText)


    ## def addOcclusionCamera(self):
    ##     from DejaVu2.Camera import OcclusionCamera
    ##     root = Tkinter.Toplevel()
    ##     root.withdraw()
    ##     root.protocol('WM_DELETE_WINDOW', self.hideOcclusionCamGUI)
 
    ##     self.occlusionCamRoot = root
    ##     self.occlusionCamera = OcclusionCamera(root, self.viewer,
    ##                                            width=50, height=50)


    ## def hideOcclusionCamGUI(self, *dummy):
    ##     self.occlusionCamRoot.withdraw()


    ## def showOcclusionCamGUI(self, *dummy):
    ##     self.occlusionCamera.listGeoms()
    ##     if self.occlusionCamRoot.winfo_ismapped() == 0:
    ##         self.occlusionCamRoot.deiconify()
    ##     self.occlusionCamRoot.lift()


##     def startScenario(self, event=None):
##         if not hasattr(self.viewer, 'director'):
##             from scenarioInterface import DejaVu2Scenario
##             from scenario.director import Director
##             director = self.viewer.director = Director()
##             scenario = DejaVu2Scenario(self.viewer, director)
##             director.addScenario("DejaVu2Scenario", scenario)  
##         else:
##             scenario = self.viewer.director.scenarios["DejaVu2Scenario"]
##         scenario.start()


    def dialog(self):
#        t="Do you Wish to Quit?"
#        from SimpleDialog import SimpleDialog
#        d = SimpleDialog(self.root, text=t,        
#                                     buttons=["Quit","Cancel"],
#                                     default=0, title="Quit?")
#        ok=d.go()
        ok = tkMessageBox.askokcancel("Quit?","Do you Wish to Quit?")
        if ok:
            self.quit_cb()
        else:
            return


    def quit_cb(self):
                   
        self.root.master.destroy()
        
        
    def setCurrentObjectFromTree_cb(self, treeNode):
        if treeNode.object:
            self.viewer.SetCurrentObject(treeNode.object)
            
        
    def toggleVisibility_cb(self, treeNode):
        geom = treeNode.object
        if geom:
            geom.Set(visible=not geom.visible)
            
        
    def withdraw(self):
        if self.shown:
            self.root.withdraw()
            self.shown = False
            

    def deiconify(self):
        if not self.shown:
            self.root.deiconify()
            self.shown = True
            

    def Exit(self):
        # kept for backwards compatibility
        self.viewer.Exit()
    

    def Exit_cb(self, event=None):
        self.Exit()

##      def Wysiwyg_cb(self):
##          self.colorChooser.Wysiwyg(self.wysiwyg.get())

    def MoveRootOnly_cb(self):
        self.viewer.TransformRootOnly( self.moveRootOnly.get() )
        

    def showPickedVertex_cb(self):
        self.viewer.showPickedVertex = self.showPickedVertex.get()
        

    def displayValueInObjList_cb(self):
        v = self.displayValueInObjList.get()
        self.tvolist.displayValueInTree(v)
        

    def showHideQuickKeys_cb(self):
        if self.showHideQuickKeysVar.get():
            self.QuickKeysFrame.deiconify()
        else:
            self.QuickKeysFrame.withdraw()
        

    def addQuickKey_cb(self):
        viewer = self.viewer
        camera = viewer.currentCamera
        # find what the trackball is bound to
        xform = self.Xform.get()

        # find the currentObject
        if xform=="Object" or xform=="Texture" or xform=="Scissor":
            obj = viewer.currentObject
        elif xform=="Camera":
            obj = viewer.currentCamera
        elif xform=="Clip":
            obj = viewer.currentClip
        elif xform=="Light":
            obj = viewer.currentLight
        else:
            raise RuntimeError("bad value for xform", xform)

        # find value of rootOnly
        rootOnly= self.moveRootOnly.get()

        cb = CallBackFunction( self.quickKey_cb, xform, obj, rootOnly )
        if xform=="Object":
            if rootOnly or obj==viewer.rootObject:
                label = "Xform Scene"
            else:
                label = "Xform "+obj.fullName
        elif xform=="Clip" or xform=="Light":
            label = "Xform "+xform+str(obj.num)
        elif xform=="Camera":
            label = "Xform Camera"
        elif xform=="Texture" or xform=="Scissor":
            label = "Xform %s %s"%(xform, obj.fullName)

        # create a button and add it to the Quick Keys panel
        button = Tkinter.Button(self.QuickKeysFrame, text=label, command=cb)
        button.pack(side='top', expand=1, fill='y')


    def quickKey_cb(self, xform, obj, rootOnly):
        if xform=="Object":
            self.viewer.SetCurrentObject(obj)
            self.viewer.TransformRootOnly(rootOnly)
            self.TObject()
        elif xform=="Camera":
            self.viewer.SetCurrentCamera(obj)
            self.TCamera()
        elif xform=="Clip":
            self.viewer.SetCurrentClip(obj)
            self.TClip()
        elif xform=="Light":
            self.viewer.SetCurrentLight(obj)
            self.TLight()
        elif xform=="Texture":
            self.viewer.SetCurrentObject(obj)
            self.viewer.SetCurrentTexture(obj)
            self.TMap()
        elif xform=="Scissor":
            self.viewer.SetCurrentObject(obj)
            self.Scissor()
            

    def bindResetButton(self, func):
        self.resetB.configure(command=func)
        

    def bindNormalizeButton(self, func):
        self.normalizeB.configure(command=func)
        

    def bindCenterButton(self, func):
        self.centerB.configure(command=func)


    def bindDeleteButton(self, func):
        self.deleteB.configure(command=func)


    def bindOwnGuiButton(self, func):
        self.OwnGuiB.configure(command=func)


    def enableResetButton(self, val):
        self.resetB.configure(state = val)
        

    def enableNormalizeButton(self, val):
        self.normalizeB.configure(state = val)
        

    def enableCenterButton(self, val):
        self.centerB.configure(state = val)
        

    def enableDeleteButton(self, val):
        self.deleteB.configure(state = val)


    def enableOwnGuiButton(self, val):
        self.OwnGuiB.configure(state = val)


    def countParents(self, object):
        c = 0
        while object.parent:
            c = c+1
            object = object.parent
        return c


    def lstripChar(self, name, char):
        n = string.count(name,'~')
        return n, name[ n : ]


    def objectByName(self, name):
        ol = self.viewer.rootObject.AllObjects()
        n, name = self.lstripChar(name, '~')
        for o in ol:
            if o.name==name: return o
        return None

         

    """


    def objectIndex(self, object):
        # object is a geometry and we find this object's index in the list of
        # names displayed in te widget. If the ibecjt is not shown we
        # return -1
        l = self.olist.get(0, 'end')
        for i in range(len(l)):
            indent, n = self.lstripChar(l[i], '~')
            if n==object.name: break
        if i==len(l): return -1
        else: return i


    def countDecendentsInWidget(self, object):
        # object is a geometry, we count and return the number of 
        # decendents shown in widget
        ind = self.objectIndex(object)
        allNames = self.olist.get(0, 'end')
        nbTild = string.count(allNames[ind],'~')+1
        # count children in widget
        nbChildren = 0
        for i in range(ind+1, len(allNames)):
            nbt = string.count(allNames[i],'~')
            if nbt >= nbTild:
                nbChildren = nbChildren + 1
            else:
                break
        return nbChildren

    # NOT USED 
    def countChildrenInWidget(self, object):
        # object is a geoemtry, we count and return the number of
        # children of object shown in widget
        ind = self.objectIndex(object)
        allNames = self.olist.get(0, 'end')
        nbTild = string.count(allNames[ind],'~')+1
        # count children in widget
        nbChildren = 0
        for i in range(ind+1, len(allNames)):
            nbt = string.count(allNames[i],'~')
            if nbt == nbTild:
                nbChildren = nbChildren + 1
            else:
                break
        return nbChildren


    
    def getObjectIndexInList(self, object):
        # object is a geometry and we find this object's index in the list of
        # names displayed in the widget. If the object is not shown we
        # return -1
        l = list(self.olist.get(0, 'end'))
        ofn = object.fullName.split('|')
        ofn = map(lambda x: ofn.index(x)*'~'+x, ofn)
        findex=0
        for name in ofn:
            en = enumerate(l)
            for ind, lname in en:
                if name==lname:
                    findex = ind + findex
                    l = l[ind:]
                    break
        if findex == len(self.olist.get(0, 'end')): return -1
        else: return findex
        
    
    def expand(self, object):
        # object is a geometry
        if object.isExpandedInObjectList: return
        object.isExpandedInObjectList = 1
        geoms = object.children
        ind = self.objectIndex(object) + 1
        c = self.countParents(object) + 1
        prefix = '~'*c
        for i in range(len(geoms)):
            g = geoms[i]
            if g==object: continue
            if not g.listed: continue
            self.olist.insert(ind, prefix + g.name)
            ind = ind + 1

            
    def collapse(self, object):
        # object is a geometry, we recursively collapse the sub-tree
        object.isExpandedInObjectList = 0

        # delete the names from the bject list widget
        nbChildren = self.countDecendentsInWidget(object)
        ind = self.objectIndex(object) + 1
        for i in range(ind, ind+nbChildren):
            self.olist.delete(ind)
        # toggle isExpandedInObjectList for all descendents
        for child in object.AllObjects():
            if child.listed:
                child.isExpandedInObjectList = 0
                
                
    def getFullName(self, ind):
        # strip the leading ~
        allNames = self.olist.get(0, 'end')
        nbTild = string.count(allNames[ind],'~')
        fullName = allNames[ind][nbTild:]
        for i in range(ind-1, -1, -1):
            nbt, name = self.lstripChar(allNames[i], '~')
            if nbt >= nbTild: continue
            nbTild = nbt
            fullName = name + '|' + fullName
        return fullName
    
    
    def toggleExpansion(self, event):
        # get a 0-based index into list of names
        o = self.olist.nearest(event.y)
            fullName = self.getFullName(o)
        #obj = self.objectByName(self.olist.get(o))
            obj = self.viewer.FindObjectByName(fullName)
        if obj:
                childGeoms = obj.AllObjects()
                if len(childGeoms)==1:  # this geoemtry has no children
                    return
                else: # this geometry has children
                    if obj.isExpandedInObjectList: self.collapse(obj)
                    else: self.expand(obj)

        
    def select(self, event):
        # get a 0-based index into list of names
        o = self.olist.nearest(event.y)
            fullName = self.getFullName(o)
        #obj = self.objectByName(self.olist.get(o))
            obj = self.viewer.FindObjectByName(fullName)
        if obj:
    # SetCurrentObject is called by BindTrackballToObject
    #            self.viewer.SetCurrentObject(obj)
            self.viewer.BindTrackballToObject(obj)
            self.viewer.currentTransfMode = 'Object'
"""
    def renameObject(self, obj, name):
        if not obj.listed:
            return
        treeNode = self.tvolist.objToNode[obj]
        treeNode.rename(name)

        

    def addObject(self, obj, parent):
        """adds object and its children recursively to the treewidget"""
        if not obj.listed:
            return

        if not parent:
            #self.olist.insert('end', obj.name)
            self.tvolist.addNode(parent=None, 
                                 name=obj.name, object=obj)
            #self.olist.select_set(0,0)

        else:
            self.tvolist.addNode(parent=self.tvolist.objToNode[parent],
                                 name=obj.name, object=obj)

        # recursively add the subtree under obj
        for child in obj.children:
            self.addObject( child, obj)

     
      #      if not parent.isExpandedInObjectList:
      #          return
            # return here.. the following lines should be removed 
            
            #i = self.objectIndex(parent)
            # Get the index of the parent in the list to insert the
            # child after it,
        #    i = self.getObjectIndexInList(parent)
        #if i==-1:
        #        return
        #c = self.countParents(obj)
        #    prefix = '~'*c
        #name = prefix + obj.name
            # now we need to skip all children already there
        #    l = self.olist.get(0, 'end')
        #    while 1:
        #        i = i + 1
        #        if i==len(l): break
        #        if self.olist.get(i)[:c]!=prefix: break

            #self.olist.insert(i, name)


    def deleteObject(self, obj):
        if not obj.listed:
            return
        # Get the index of the object in the list to
        # delete it.
        #i = self.getObjectIndexInList(obj)
        #self.olist.delete(i)
        self.tvolist.deleteNode(self.tvolist.objToNode[obj])


    def SetCurrentObject(self, obj):
        if obj.protected is True or obj.protected == 1:
            self.enableDeleteButton(Tkinter.DISABLED)
        else:
            self.enableDeleteButton(Tkinter.NORMAL)        

        if hasattr(obj, 'createOwnGui'):
            self.bindOwnGuiButton(obj.showOwnGui)
            self.enableOwnGuiButton(Tkinter.NORMAL)
        else:
            self.enableOwnGuiButton(Tkinter.DISABLED)

        # MS. the part below should happen even if there is no GUI
        # but cant be done in ViewerBase because self.Xform.get()
        if self.Xform.get() == 'Object':
            if isinstance(obj, Transformable):
                self.viewer.currentCamera.bindAllActions('Object')
            elif isinstance(obj, Insert2d):
                self.viewer.currentCamera.bindAllActions('Insert2d')
        #print "SetCurrentObject"            
        #Update GUI for current object
        
        #if obj.listed:
        #    names = string.split(obj.fullName, '|')
        #    top = self.viewer.topObjects
        #    l = self.olist.get(0, END)
        #    i = 0
        #    for name in l:
        #        if names[0]==l[i]:
        #            break
        #        i = i + 1
        #    for name in names[1:]:
        #        while (i<len(l)):
        #            ind, n = self.lstripChar(l[i], '~')
        #            if n==name: break
        #            i = i + 1

        #    if i>0 and i==len(l):
        #        return
        #    cl = self.olist.curselection()
        #    if len(cl):
        #        self.olist.select_clear(cl)

        #    self.olist.select_set(i,i)
        #    self.olist.see(i)
        #else:
        #    cl = self.olist.curselection()
        #    if len(cl): self.olist.select_clear(cl)
        
        lIsInstanceTransformable = isinstance(obj, Transformable)
        
        # Set the Node corresponding to 'obj' be selected
        p=self.viewer.currentObject
        try:
            node = self.tvolist.objToNode[p]
        except KeyError:
            self.viewer.currentObject = self.viewer.rootObject
            node = self.tvolist.objToNode[self.viewer.currentObject]
            
        if self.tvolist.current_selected != node:
            #print "trying to select", node.name , "over :", \
            #                    self.tvolist.current_selected
            self.tvolist.Select(node)
        
        # update the object properties panel
        for prop, var in self.inheritVar.items():
            if hasattr(obj,prop):
                if (obj is self.viewer.rootObject) \
                  and (prop.startswith('inherit')):
                    var.set(False)
                else:
                    propVal = getattr(obj, prop)
                    var.set(propVal)
            elif prop == 'disableTexture':
                if hasattr(obj, 'texture') and obj.texture:
                    var.set(not obj.texture.enabled)
                else:
                    var.set(False)

#        # update the children properties panel
#        self.initChildPropPanel(self.viewer.currentObject)

        #self.showObject.set(obj.visible)

        mated = self.viewer.materialEditor
        val = self.objMatEdTk.get()
        if    lIsInstanceTransformable \
          and obj.inheritMaterial \
          and obj != self.viewer.rootObject:
            mated.dismiss()
            self.objMatEdTk.set(0)
            self.MatEdB1.configure(state='disabled')
            self.MatEdB2.configure(state='disabled')
        else:
            self.MatEdB1.configure(state='normal')
            self.MatEdB2.configure(state='normal')
            if val==1:
                mated.setObject(obj, GL.GL_FRONT)
            elif val==2:
                mated.setObject(obj, GL.GL_BACK)

            
        for cp in self.viewer.clipP:
            self.UpdateGuiCP(cp, obj)

##         if obj.pickableVertices == viewerConst.YES:
##             self.pickableVertices.set(1)
##         else: self.pickableVertices.set(0)

##         if obj.antialiased == viewerConst.YES:
##             self.antialiasedLines.set(1)
##         else: self.antialiasedLines.set(0)

##         self.scissortk.set( obj.scissor)

#        self.depthMasktk.set( obj.depthMask )
        
        if lIsInstanceTransformable:
            self.lw.Set(obj.lineWidth, update=0)
            self.pw.Set(obj.pointWidth, update=0)
##          self.opacw.Set(obj.materials[GL.GL_FRONT].prop[1][0][3], update=0)

            lCascadeMenu = self.opacityRelabelingCascadeMenus['srcBF']
            lCascadeMenu.setWithoutCallbackFunction( obj.srcBlendFunc )
            lCascadeMenu = self.opacityRelabelingCascadeMenus['dstBF']
            lCascadeMenu.setWithoutCallbackFunction( obj.dstBlendFunc )

            lCascadeMenu = self.relabelingCascadeMenus['frontPolyMode']
            if    (obj is not self.viewer.rootObject) \
              and (obj.inheritFrontPolyMode is True):
                lCascadeMenu.setWithoutCallbackFunction(viewerConst.INHERIT)
            elif obj.drawOutline[0]:
                lCascadeMenu.setWithoutCallbackFunction(viewerConst.OUTLINED)
            else:
                lCascadeMenu.setWithoutCallbackFunction(obj.frontPolyMode)

            lCascadeMenu = self.relabelingCascadeMenus['backPolyMode']
            if obj.frontAndBack:
                lCascadeMenu.setWithoutCallbackFunction(GL.GL_FRONT_AND_BACK)
            elif    (obj is not self.viewer.rootObject) \
              and (obj.inheritBackPolyMode is True):
                lCascadeMenu.setWithoutCallbackFunction(viewerConst.INHERIT)
            elif obj.drawOutline[1]:
                lCascadeMenu.setWithoutCallbackFunction(viewerConst.OUTLINED)
            else:
                lCascadeMenu.setWithoutCallbackFunction(obj.backPolyMode)

            for prop in self.cascadeProp:
                lCascadeMenu = self.relabelingCascadeMenus[prop]
                inheritProp = 'inherit' + prop[0].upper() + prop[1:]
                if    (obj is not self.viewer.rootObject) \
                  and (getattr(obj, inheritProp) is True):
                    lCascadeMenu.setWithoutCallbackFunction(viewerConst.INHERIT)
                else:
                    lCascadeMenu.setWithoutCallbackFunction( getattr(obj, prop) )

        if len(obj.children) == 0:
            self.menuEdit.menu.entryconfig('Apply Transformation',
                state='normal') 
        else:
            self.menuEdit.menu.entryconfig('Apply Transformation', 
                state='disabled')

        if obj == self.viewer.rootObject:
            for lButton in self.objectInheritButtons.values():
                lButton.configure(state='disabled')
            for lIndex in self.inheritMenu_b.indices:
                self.inheritMenu_b.menu.entryconfigure(lIndex, state='disabled')
            for lCascadeMenu in self.relabelingCascadeMenus.values():
                lCascadeMenu.entryconfigure('inherit', state='disabled')
            for lPolyModeMenu in self.polyModeMenus:
                lPolyModeMenu.menu.entryconfigure('inherit', state='disabled')
        else:
            for lButton in self.objectInheritButtons.values():
                lButton.configure(state='normal')
            for lIndex in self.inheritMenu_b.indices:
                self.inheritMenu_b.menu.entryconfigure(lIndex, state='normal')
            for lCascadeMenu in self.relabelingCascadeMenus.values():
                lCascadeMenu.entryconfigure('inherit', state='normal')
            for lPolyModeMenu in self.polyModeMenus:
                lPolyModeMenu.menu.entryconfigure('inherit', state='normal')

        self.updateOMPgui()
        

    def read(self):
        print 'read'


    def write(self):
        print 'write'


    def delete(self):
        print 'delete'


    def ObjectProp(self):
        print 'ObjProp'


## can't get gluProject to work properly !
##          v = obj.vertexSet.vertices.array
##          mini = Numeric.minimum.reduce(v)
##          mini = obj.Project(mini)
##          print 'mini:', mini
##          maxi = Numeric.maximum.reduce(v)
##          maxi = obj.Project(maxi)
##          print 'maxi:', maxi
##          diffX = (maxi[0]-mini[0])*0.25
##          diffY = (maxi[1]-mini[1])*0.25
##          obj.Set( scissor = self.scissortk.get(),
##                   scissorX = mini[0]+diffX, scissorY = mini[0]+diffY,
##                   scissorW = mini[1]+(2*diffX), scissorH = mini[1]+(2*diffY) )

    def zSortPoly(self):
        obj = self.viewer.currentObject
        if isinstance(obj, IndexedGeom): # or isinstance(obj, Spheres):
            obj.sortPoly()
        self.viewer.Redraw()
        

    def minuszSortPoly(self):
        obj = self.viewer.currentObject
        if isinstance(obj, IndexedGeom): # or isinstance(obj, Spheres):
            obj.sortPoly(order=1)
        self.viewer.Redraw()
        

    def SetBF(self):
        obj = self.viewer.currentObject
        if isinstance(obj, IndexedGeom):
            obj.Set( blendFunctions=(self.srcBFtk.get(), self.dstBFtk.get()) )
        self.viewer.Redraw()


    def setChildProp_cb(self, propname, recursive):
        #print "setChildProp_cb", propname
        obj = self.viewer.currentObject
        if propname == 'scissor' :
            obj.SetForChildren( scissor=obj.scissor,
                                scissorX=obj.scissorX,
                                scissorY=obj.scissorY,
                                scissorW=obj.scissorW,
                                scissorH=obj.scissorH )
        elif hasattr(obj, propname):
            propval = getattr(obj, propname)
            apply( obj.SetForChildren, (recursive,), {propname:propval})

        self.viewer.deleteOpenglList()
        self.viewer.Redraw()


    def setChildInheritProp_cb(self, propname, recursive):
        #print "setChildInheritProp_cb", propname
        obj = self.viewer.currentObject
        if isinstance(obj, Geom):
            inheritPropName = 'inherit' + propname[0].upper() + propname[1:]
            #print "inheritPropName", inheritPropName
            inheritPropVal = getattr(obj, inheritPropName)
            apply( obj.SetForChildren, (recursive,), {inheritPropName:inheritPropVal})
            propVal = getattr(obj, propname)
            apply( obj.SetForChildren, (recursive,), {propname:propVal})
            self.viewer.deleteOpenglList()
            self.viewer.Redraw()


    def setInherit_cb(self, propname, variable):
        #print "setInherit_cb", propname, variable.get()
        # 1- Get a handle to the current object.
        obj = self.viewer.currentObject
        # 2- Get the propname and the property value
        propval = variable.get()

        # 3- props that dont need to call obj.set
        if propname == 'disableTexture':
            if hasattr(obj, 'texture') and obj.texture:
                obj.texture.Set(enable=not propval)
        elif propname == 'vertexNormals':
            if propval:
                obj.addVertexNormalsGeom()
            else:
                obj.removeVertexNormalsGeom()
        elif propname == 'faceNormals':
            if propval:
                obj.addFaceNormalsGeom()
            else:
                obj.removeFaceNormalsGeom()
        else:
            # 4- Set the propname of the current object to the given value
            apply(obj.Set, (), {propname:propval})
            if propname=='scissor':
                if propval:
                    self.viewer.currentScissor=obj
                else:
                    self.viewer.currentScissor=None
            # 5- If inheritMaterial needs to do some other stuff.
            if propname == "inheritMaterial":
                if propval: # and obj != self.viewer.rootObject:
                    self.viewer.materialEditor.dismiss()
                    self.objMatEdTk.set(0)
                    self.MatEdB1.configure(state='disabled')
                    self.MatEdB2.configure(state='disabled')
                else:
                    self.MatEdB1.configure(state='normal')
                    self.MatEdB2.configure(state='normal')
            elif propname == 'scissor':
                if propval:
                    c = self.viewer.currentCamera
                    w = c.width/4
                    h = c.height/4
                    obj.Set( scissor = propval,
                             scissorX = w, scissorY = h,
                             scissorW = 2*w,
                             scissorH = 2*h )
                else:
                    obj.Set( scissor = propval)
            elif propname == 'invertNormals':
                if propval:
                    obj.invertNormals = True
                else:
                    obj.invertNormals = False
            
        # 6- Calls the redodisplaylist of the viewer and trigger a redraw
        self.viewer.deleteOpenglList()
        self.viewer.Redraw()


    def ObjMatEdFront_cb(self):
        """start material editor define the diffuse component of the current object"""
        self.objMatEdTk.set(1)
        face = GL.GL_FRONT
        self.viewer.materialEditor.root.title('Front Material Editor')
        mated = self.viewer.materialEditor
        mated.show()
        obj = self.viewer.currentObject
        mated.setObject(obj, face)
        #mated.defineMaterial(obj.materials[face].prop, face)


    def ObjMatEdBack_cb(self):
        """start material editor define the diffuse component of the current object"""
        self.objMatEdTk.set(2)
        face = GL.GL_BACK
        self.viewer.materialEditor.root.title('Back Material Editor')
        mated = self.viewer.materialEditor
        mated.show()
        obj = self.viewer.currentObject
        mated.setObject(obj, face)
        #mated.defineMaterial(obj.materials[face].prop, face)


    def ObjectDelete(self):
        """Remove current Object from viewer and make its parent current"""
        obj = self.viewer.currentObject
        if obj == self.viewer.rootObject: p = obj
        else: 
            p = obj.parent
            self.viewer.RemoveObject(obj)
        self.viewer.SetCurrentObject(p)


    def SetFrontPolyMode(self):
        """Modify the current Object's frontPolyMode"""
        obj = self.viewer.currentObject
        if isinstance(obj, Geom):
            obj.SetFrontPolyMode(
               self.relabelingCascadeMenus['frontPolyMode'].cascadeVariable.get()
              )
            self.viewer.Redraw()


    def SetBackPolyMode(self):
        """Modify the current Object's frontPolyMode"""

        obj = self.viewer.currentObject
        if isinstance(obj, Geom):
            lValue = self.relabelingCascadeMenus['backPolyMode'].cascadeVariable.get()
            if lValue == GL.GL_FRONT_AND_BACK:
                obj.frontAndBack = viewerConst.YES
                obj.SetBackPolyMode( 
                   self.relabelingCascadeMenus['frontPolyMode'].cascadeVariable.get()
                  )
            else:
                obj.frontAndBack = viewerConst.NO
                obj.SetBackPolyMode( lValue )
            self.viewer.Redraw()


    def setPropMode(self, prop, val):
        """Modify the current Object's prop mode to val
"""
        #print "setPropMode", prop, val.get()
        obj = self.viewer.currentObject
        if isinstance(obj, Geom):
            val = val.get()
            kw = {}
            inheritProp = 'inherit' + prop[0].upper() + prop[1:]
            #print "inheritProp", inheritProp
            if val == viewerConst.INHERIT:
               kw[inheritProp] = True
            else:
               kw[inheritProp] = False
               kw[prop] = val
                
            obj.Set(**kw)
            #print prop, getattr(obj, inheritProp), getattr(obj, prop)


    def reset(self, event):
        print 'reset'


    def normalize(self, event):
        print 'normalize'


    def center(self, event):
        print 'center'


    def TObject(self):
        self.viewer.BindTrackballToObject(self.viewer.currentObject)


    def TCamera(self):
        self.viewer.BindTrackballToCamera(self.viewer.currentCamera)


    def TClip(self):
        self.viewer.BindTrackballToClip(self.viewer.currentClip)


    def TLight(self):
        self.viewer.BindTrackballToLight(self.viewer.currentLight)


    def TMap(self):
        self.viewer.BindTrackballToTexture(self.viewer.currentObject)


    def Scissor(self):
        self.viewer.BindTrackballToScissor(self.viewer.currentObject)


    def PObject(self):
        """Make the object panel visible"""
        if hasattr(self, 'CurrentPropPanel') and self.CurrentPropPanel != self.ObjProp:
            self.CurrentPropPanel.forget()
            self.CurrentPropPanel = self.ObjProp
            self.ObjProp.pack(fill=X, expand=1)
            self.sframe.forget()
            self.sframe.pack()
    

    def PCamera(self):
        """Make the camera panel visible"""
        if self.CurrentPropPanel != self.CameraProp:
            self.CurrentPropPanel.forget()
            self.CurrentPropPanel = self.CameraProp
            self.CameraProp.pack(fill=X, ipadx=5, ipady=3, expand=1)
            self.sframe.forget()
            self.sframe.pack()


    def PLight(self):
        """Make the light panel visible"""
        if self.CurrentPropPanel != self.LightProp:
            self.CurrentPropPanel.forget()
            self.CurrentPropPanel = self.LightProp
            self.LightProp.pack(fill=X, expand=1)
            self.sframe.forget()
            self.sframe.pack()


    def PClip(self):
        """Make the clipping planes panel visible"""
        if self.CurrentPropPanel != self.ClipProp:
            self.CurrentPropPanel.forget()
            self.CurrentPropPanel = self.ClipProp
            self.ClipProp.pack(fill=X, expand=1)
            self.sframe.forget()
            self.sframe.pack()


    def PViews(self):
        """Make the views planes panel visible"""
        if self.CurrentPropPanel != self.BookmarksProp:
            self.CurrentPropPanel.forget()
            self.CurrentPropPanel = self.BookmarksProp
            self.BookmarksProp.pack(fill=X, expand=1)
            self.sframe.forget()
            self.sframe.pack()


    def DrawSceneBB_cb(self):
        """Toggles depthcueing"""
        camera = self.viewer.currentCamera
        camera.Set(boundingbox=self.drawSceneBB.get())
        self.viewer.Redraw()


    def ToggleThumbnail(self):
        """Toggles drawing thumbnail"""
        camera = self.viewer.currentCamera
        camera.Set(drawThumbnail=self.drawThumbnail.get())
        self.viewer.Redraw()


    def showCurveTool(self, event=None):
        
        if self.GraphToolpanel.winfo_ismapped()==0:
            self.GraphToolpanel.deiconify()
        self.GraphToolpanel.lift()

    def showSSAOTool(self,event=None):
        #should show some windowswith ssa optio...and pass the option
        #pass
        camera = self.viewer.currentCamera
        if not hasattr(self,"ssao_options_w"):
            self.ssao_options_w = Toplevel(self.root)
            self.ssao_options_w.title("SSAO")
            #should put here preset, quality ? scale ?
            #pull-down menu with preset
            #need an advanced button that will show the advance oiption -> collapsable ?
            self.ssao_options_w_widget ={}
            Label(self.ssao_options_w, text="Screen Space Ambient Occlusion").grid(row=0, column=0, columnspan = 2,sticky='w')
            r=1
            #order?
            for k in camera.SSAO_OPTIONS_ORDER:
                elem = camera.SSAO_OPTIONS[k]
                Label(self.ssao_options_w, text=k).grid(row=r, column=0, sticky='w')
                def updateSSAO_OPTIONS_cb(*args,**kw):
                    self.viewer.Redraw()                
                if k in ["do_noise","only_ssao","use_mask_depth","show_mask_depth",
                         "fog","use_fog","only_depth","mix_depth","negative"]:
                    pw = Tkinter.IntVar()
                    w = Tkinter.Checkbutton(self.ssao_options_w, 
                                            var=pw,command=updateSSAO_OPTIONS_cb)
                    w.grid(row=r, column=1, columnspan=3) 	                
                    pw.set(camera.SSAO_OPTIONS[k][0])
                else :
                    pw = ThumbWheel(self.ssao_options_w,
                showLabel=1, width=60, height=16, type=elem[3], value=elem[0],
                callback=updateSSAO_OPTIONS_cb, continuous=True, oneTurn=10.,min=elem[1],
                wheelPad=2, #labCfg = {'text':'fog end:', 'side':'left'}
                )
                    pw.grid(row=r, column=1, columnspan=3) 
                camera.SSAO_OPTIONS[k].append(pw)
                self.ssao_options_w_widget[k] = pw
                r+=1
            self.ssao_options_w.protocol('WM_DELETE_WINDOW', self.ssao_options_w.withdraw)
        else :
            for k in camera.SSAO_OPTIONS:
                elem = camera.SSAO_OPTIONS[k]
                self.ssao_options_w_widget[k].set(elem[0])
                if len(elem) < 5 :
                    camera.SSAO_OPTIONS[k].append(self.ssao_options_w_widget[k])
        self.ssao_options_w.deiconify()

    def toggleSSAO(self, event=None):
        camera = self.viewer.currentCamera
        camera.ssao = self.ssaoTk.get()
        self.viewer.Redraw()
         
    def showHistogram(self):
        camera = self.viewer.currentCamera
        im=camera.GrabZBuffer()
        from PIL import Image
        #print im.getbands()
        self.curvetool.histvalues=im.histogram()
        self.curvetool.drawHistogram()
        
        
    def continuousRamp(self, event=None):
         
        camera = self.viewer.currentCamera
        camera.Set(d1ramp=self.curvetool.get())
        self.viewer.Redraw()
        if self.GraphToolpanel.winfo_ismapped()==1:
            self.showHistogram() 
        
        
    def setNPR_cb(self,event=None):
        camera = self.viewer.currentCamera
        camera.Set(d1scale= self.d1scalewheel.get())
        self.viewer.Redraw()

    
    def toggleOutline(self, event=None):
        camera = self.viewer.currentCamera
        camera.Set(contours=self.contourTk.get())
        self.viewer.Redraw()
        

    def AutoDepthcue(self):
        """set depthcueing using bounding box
"""
        #print "ViewerGUI.AutoDepthcue"
        vi = self.viewer
        vi.currentCamera.AutoDepthCue()
        vi.Redraw()


    def SetProjection(self):
        """Set jitter values for scene anti aliasing"""

        camera = self.viewer.currentCamera
        camera.Set(projectionType=self.projType.get())
        self.viewer.Redraw()


    def SetJitter(self):
        """Set jitter values for scene anti aliasing"""

        camera = self.viewer.currentCamera
        camera.Set(antialiased=self.nbJitter.get())
        self.viewer.Redraw()


    def CameraVideoRecorder_cb(self):
        """Show video recorder GUI."""
        camera = self.viewer.currentCamera
        
        if camera.videoRecorder:
            st = camera.videoRecorder.form.root.winfo_ismapped()
            if st:
                camera.videoRecorder.close_cb()
            else:
                camera.videoRecorder.buildForm()
        else:
            filename = "out.mpg"
            from DejaVu2.videoRecorder import Recorder
            root = camera.master.master.master
            camera.videoRecorder = Recorder(root,
                                            filetypes=[("MPG", ".mpg")],
                                            fileName = filename, camera=camera)


    def BindTrackballToObject(self, obj, allCameras=None):

	if allCameras:
	    for c in self.cameras:
                self.bindResetButton( self.Reset_cb)
                self.enableNormalizeButton(Tkinter.NORMAL)
                self.enableCenterButton(Tkinter.NORMAL)
	else:
	    c = self.currentCamera
            self.bindResetButton( self.Reset_cb)
            self.enableNormalizeButton(Tkinter.NORMAL)
            self.enableCenterButton(Tkinter.NORMAL)

        self.fillTransformInfo_cb()


    def BindTrackballToClip(self, c, allCameras=None):
	if allCameras:
	    for c in self.cameras:
                self.bindResetButton( self.Reset_cb)
                self.enableNormalizeButton(Tkinter.DISABLED)
                self.enableCenterButton(Tkinter.DISABLED)
	else:
            self.bindResetButton( self.Reset_cb)
            self.enableNormalizeButton(Tkinter.DISABLED)
            self.enableCenterButton(Tkinter.DISABLED)

        self.fillTransformInfo_cb()


    def BindTrackballToCamera(self, c, allCameras=None):
	"""Bind the trackball to the current camera"""

	if allCameras:
	    for c in self.cameras:
                self.bindResetButton( self.Reset_cb )
                self.enableNormalizeButton(Tkinter.DISABLED)
                self.enableCenterButton(Tkinter.DISABLED)
	else:
            self.bindResetButton( self.Reset_cb )
            self.enableNormalizeButton(Tkinter.DISABLED)
            self.enableCenterButton(Tkinter.DISABLED)

        self.fillTransformInfo_cb()


    def BindTrackballToLight(self, light, allCameras=None):
	"""Bind the trackball to the current ligth source"""
	
	if light.positional==viewerConst.NO:  # directional light
	    if allCameras:
		for c in self.cameras:
                    self.bindResetButton( self.Reset_cb)
                    self.enableNormalizeButton(Tkinter.DISABLED)
                    self.enableCenterButton(Tkinter.DISABLED)
	    else:
                self.bindResetButton( self.Reset_cb)
                self.enableNormalizeButton(Tkinter.DISABLED)
                self.enableCenterButton(Tkinter.DISABLED)

        self.fillTransformInfo_cb()


    def BindTrackballToTexture(self, o, allCameras=None):
        """Bind trackball to the texture of the current object
"""

        if allCameras:
            for c in self.cameras:
                self.bindResetButton( self.Reset_cb)
                self.enableNormalizeButton(Tkinter.DISABLED)
                self.enableCenterButton(Tkinter.DISABLED)
        else:
            self.bindResetButton( self.Reset_cb)
            self.enableNormalizeButton(Tkinter.DISABLED)
            self.enableCenterButton(Tkinter.DISABLED)

        self.fillTransformInfo_cb()


    def BindTrackballToScissor(self, o, allCameras=None):
        """Bind trackball to the scissor of the current object
"""
        if allCameras:
            for c in self.cameras:
                c.bindAllActions('Scissor')
                self.bindResetButton( self.Reset_cb )
                self.enableNormalizeButton(Tkinter.DISABLED)
                self.enableCenterButton(Tkinter.DISABLED)
        else:
            self.bindResetButton( self.Reset_cb )
            self.enableNormalizeButton(Tkinter.DISABLED)
            self.enableCenterButton(Tkinter.DISABLED)

        self.fillTransformInfo_cb()


    def SetCurrentCamera(self):
        """Update GUI for current camera
"""
        c = self.viewer.currentCamera
        f = self.viewer.currentCamera.fog
        self.NearFarFog.Set(c.near, c.far, f.start, f.end)
        if f.enabled:
            self.viewer.depthcued = 1
            self.depthcued.set(1)
        else:
            self.viewer.depthcued = 0
            self.depthcued.set(0)
        self.projType.set(self.viewer.currentCamera.projectionType)


    def ToggleLocalViewer(self):
        """turn local viewer mode on and off"""
    
        lm = self.viewer.lightModel
        lm.Set(localViewer = self.localViewer.get())
        self.viewer.Redraw()


    def ToggleTwoSide(self):
        """Turn two side polgon on and off"""
    
        lm = self.viewer.lightModel
        lm.Set(twoSide = self.twoSide.get())
        self.viewer.Redraw()


    def LightSelect(self):
        """Set current Light"""
    
        v = self.viewer
        l = v.lights[self.CurrentLight.get()-1]
        v.SetCurrentLight( l )


    def SetCurrentLight(self, obj):
        """Update GUI for current light"""
    
        v = self.viewer
        l = v.lights[self.CurrentLight.get()-1]
        self.lightOnOff.set(l.enabled)
        self.showLight.set(l.visible)


    def LightOnOff(self):
        """turn current light on and off"""
    
        v = self.viewer
        l = v.lights[self.CurrentLight.get()-1]
        l.Set(enabled=self.lightOnOff.get())
        v.Redraw()


    def LightShow(self):
        """Mke light source visible"""
        self.viewer.currentLight.Set(visible=self.showLight.get())
        self.viewer.Redraw()


    def LightDirectional(self):
        print 'LightDirectional'

    def LightPositional(self):
        print 'LightPositional'

    def LightSpot(self):
        print 'LightSpot'

    ###################################
    #  CLIP
    ###################################

    def UpdateGuiCP(self, cp, obj):
        """Update the Clip panel buttons for cp"""
        
        if isinstance(obj, Transformable) is False:
            return

        v = self.viewer
        cpn = cp.num
        clipI = cp in obj.clipPI
        clip = cp in obj.clipP

        self.clipvar[cpn][0].set( (clipI or clip) ) # on/off
        side = max(0, v.currentObject.clipSide[cpn])
        self.clipvar[cpn][1].set( side ) # side
        self.clipvar[cpn][2].set( clipI ) # inherit
        self.clipvar[cpn][3].set( cp.visible ) # visible
        if v.currentObject.cap:
            self.clipvar[cpn][4].set( 1 ) # capped GL
        else:
            self.clipvar[cpn][4].set( 0 ) # not capped GL
        if v.currentObject.capMesh:
            self.clipvar[cpn][5].set( 1 ) # capped Mesh
        else:
            self.clipvar[cpn][5].set( 0 ) # not capped Mesh


    def SetCurrentClip(self, oldcp):
        """Update GUI for current Cliping plane"""
    
        v = self.viewer
        cp = v.currentClip
        self.UpdateGuiCP(cp, v.currentObject)
        for i in range(6):
            self.clipw[oldcp.num][i].configure(state=DISABLED)
            self.clipw[cp.num][i].configure(state=NORMAL)
        self.CurrentClip.set( cp.num+1 )


    def centerClipPlane(self, obj, cp, inh):
        # compute center of clipped object
        center = None
        if inh: # inherited
            verts = [0,0,0]
            i = 0
            for g in obj.AllObjects():
                ve = g.getVertices()
                if len(ve):
                    verts += numpy.sum(ve, 0)/len(ve)
                    i += 1
            center = verts/i
        else:
            ve = obj.getVertices()
            if len(ve):
                center = numpy.sum(ve, 0)/len(ve)

        # center clipping plane
        if center is not None:
            cp.Set(translation=center)
        

    def ClipOnOff(self):
        """Add this clipping plane to the current object and make it current"""
        v = self.viewer
        cp = v.currentClip
        cpn = cp.num
        on = self.clipvar[cpn][0].get()
        obj = self.viewer.currentObject

        if on:
            if obj.culling!=GL.GL_NONE or obj.inheritCulling:
                obj.Set(culling='none', inheritCulling=False, setCurrent=False)
            if self.clipvar[cpn][1].get(): side=1
            else: side = -1
            if self.clipvar[cpn][2].get(): inh = True
            else: inh = False
            #if self.clipvar[cpn][3].get(): cp.visible=1
            #else: cp.visible = 0
            # show the plane by default
            self.clipvar[cpn][3].set(1)
            cp.visible=1
            obj.AddClipPlane( cp, side, inh )
            self.centerClipPlane(obj, cp, inh)
        else:
            obj.RemoveClipPlane( cp )

        if obj != self.viewer.rootObject:
            self.viewer.objectsNeedingRedo[obj] = None

        v.Redraw()


    def ClipSide(self):
        v = self.viewer
        cp = v.currentClip
        obj = self.viewer.currentObject
        if self.clipvar[cp.num][1].get():
            v.currentObject.clipSide[cp.num] = 1
        else:
            v.currentObject.clipSide[cp.num] = -1

        if obj != self.viewer.rootObject:
            self.viewer.objectsNeedingRedo[obj] = None

        v.Redraw()


    def ClipInherit(self):
        v = self.viewer
        cp = v.currentClip
        obj = self.viewer.currentObject
        inh = self.clipvar[cp.num][2].get()
        if inh:
            if cp in v.currentObject.clipPI: return
            if cp in v.currentObject.clipP:
                obj.clipPI.append(cp)
                obj.clipP.remove(cp)
        else:
            if cp in v.currentObject.clipP: return
            if cp in v.currentObject.clipPI:
                obj.clipP.append(cp)
                obj.clipPI.remove(cp)
        if obj != self.viewer.rootObject:
            self.viewer.objectsNeedingRedo[obj] = None

        self.centerClipPlane(obj, cp, inh)

        v.Redraw()


    def ClipVisible(self):
        v = self.viewer
        cp = v.currentClip
        val = self.clipvar[cp.num][3].get()
        cp.Set(visible= (val!=0))
        obj = self.viewer.currentObject
        if obj != self.viewer.rootObject:
            self.viewer.objectsNeedingRedo[obj] = None
        v.Redraw()


    def ClipCapMesh(self):
        cp = vi.currentClip
        obj = vi.currentObject
        obnOff = self.clipvar[cp.num][5].get()
        cp.ClipCapMesh(obj, onOff)


    def ClipCap(self):
        v = self.viewer
        cp = v.currentClip
        obj = self.viewer.currentObject
        val = self.clipvar[cp.num][4].get()
        if val:
            obj.cap = GL.GL_CLIP_PLANE0 + cp.num
            obj.immediateRendering = 1
        else:
            obj.cap = None
            obj.immediateRendering = 0
        if obj != self.viewer.rootObject:
            self.viewer.objectsNeedingRedo[obj] = None
        v.Redraw()


    def ClipSelect(self):
        """Set current Clipping plane"""
    
        v = self.viewer
        if v.currentClip:
            cpn = v.currentClip.num
            for i in range(6):
                self.clipw[cpn][i].configure(state=DISABLED)
        cp = v.clipP[self.CurrentClip.get()-1]
        v.SetCurrentClip( cp )


    def saveViewerState_cb(self, event=None):
        filename = self.getSaveFileName()
        if filename:
            self.viewer.saveViewerState(filename)


    def saveObjectsStates_cb(self, event=None):
        filename = self.getSaveFileName()
        if filename:
            self.viewer.saveObjectsStates(filename)


    def saveViewerAndObjectsStates_cb(self, event=None):
        filename = self.getSaveFileName()
        if filename:
            self.viewer.saveViewerAndObjectsStates(filename)


    def restoreState(self, event=None):
        filename = self.getLoadFileName()
        if filename:
            execfile(filename, {'vi':self.viewer, 'mode':'both'})
            self.viewer.Redraw()
        
        
    # added by D. Stoffler, TSRI, Jul 2003
    def saveTransformCurrentGeom(self, event=None, askGUI=1):
        geom = self.viewer.currentObject
        if askGUI:
            yesno = tkMessageBox.askokcancel('Saving Transformation',
                      'Saving matrix for geom "%s"?'%geom.name)
            if not yesno:
                return
        filename = self.getSaveFileName()
        if filename:
            self.saveTransformation(filename, geom)
            

    def loadTransformCurrentGeom(self, event=None, askGUI=1):
        geom = self.viewer.currentObject
        if askGUI:
            yesno = tkMessageBox.askokcancel('Restoring Transformation',
                      'Apply matrix on geom "%s"?'%geom.name)
            if not yesno:
                return
        filename = self.getLoadFileName()
        if filename:
            self.loadTransformation(filename, geom)


    def getSaveFileName(self, event=None):
        file = self.saveFileBrowser.get()
        return file

    
    def getLoadFileName(self, event=None):
        file = self.openFileBrowser.get()
        return file
    

    def saveTransformation(self, filename, geom):
        if filename:
            # FIXME won't work with instance matrices
            data = geom.GetMatrix(geom) # in relation to root, if called
            # without arguments, it also adds root transformation
            data = Numeric.reshape(data, (4,4))
            self.lastDir = os.path.split(filename)[0]
            f = open(filename, 'w')
            f.writelines('# File saved by DejaVu2. Do not modify!\n')
            f.writelines('import numpy.oldnumeric as Numeric\n')
            f.writelines('trans = Numeric.array(%s)\n'%list(geom.translation))
            f.writelines('rot = Numeric.array(%s).astype("f")\n'%list(
                geom.rotation))
            f.writelines('scale = Numeric.array(%s)\n'%list(geom.scale))
            f.close()


    def loadTransformation(self, filename, geom):
        if filename:
            mydict = {}
            self.lastDir = os.path.split(filename)[0]
            execfile(filename,  globals(), mydict)

            old = 0
            if geom.viewer.redirectTransformToRoot == 1:
                old = 1
                geom.viewer.TransformRootOnly(0)
            geom.FrameTransform()
            geom.translation = mydict['trans']
            geom.rotation = mydict['rot']
            geom.scale = mydict['scale']
            geom.viewer.deleteOpenglList()
            geom.viewer.Redraw()
            if old == 1:
                geom.viewer.TransformRootOnly(1)
 

    def applyTransformation_cb(self, event=None, askGUI=1):
        """Apply the current transformation to the coordinates of the
        currently selected geometry in the viewer."""
        geom = self.viewer.currentObject
        if askGUI:
            yesno = tkMessageBox.askokcancel('Apply Transformation',
                      'Apply current transformation on "%s"?'%geom.name)
            if not yesno:
                return

        # build Transformation matrix for geom only
        # FIXME won't work with instance matrices
        mat = geom.GetMatrix(root=geom)
        def applyMat(mat, pt):
            ptx = mat[0][0]*pt[0]+mat[0][1]*pt[1]+mat[0][2]*pt[2]+mat[0][3]
            pty = mat[1][0]*pt[0]+mat[1][1]*pt[1]+mat[1][2]*pt[2]+mat[1][3]
            ptz = mat[2][0]*pt[0]+mat[2][1]*pt[1]+mat[2][2]*pt[2]+mat[2][3]
            return (ptx, pty, ptz)

        # compute new coordinates
        newPts = []
        for p in geom.vertexSet.vertices.array:
            newPts.append( applyMat(mat, p) )

        # apply to vertex Set
        geom.Set(vertices=newPts)
        geom.ResetTransformation()
        geom.RedoDisplayList()
