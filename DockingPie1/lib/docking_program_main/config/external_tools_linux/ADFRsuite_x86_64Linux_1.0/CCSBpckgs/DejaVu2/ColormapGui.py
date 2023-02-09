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
# Copyright: M. Sanner TSRI 2000
#
# Revision: Guillaume Vareille
#
#############################################################################
#
# $Header: /mnt/raid/services/cvs/DejaVu2/ColormapGui.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#
# $Id: ColormapGui.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#

import types, numpy, os, sys
import Tkinter, Pmw
import tkFileDialog
from string import strip, split, find

from opengltk.OpenGL import GL
from opengltk.extent import _gllib as gllib
from mglutil.util.misc import deepCopySeq
from mglutil.gui.BasicWidgets.Tk.thumbwheel import ThumbWheel
from mglutil.util.colorUtil import HSL2RGB, RGB2HSL, \
     RGBA2HSLA_list, HSLA2RGBA_list, ToHSV, ToHEX
from mglutil.gui import widgetsOnBackWindowsCanGrabFocus

from DejaVu2 import viewerConst
from DejaVu2.colorMap import ColorMap
from DejaVu2.MaterialEditor import OGLWidget
from DejaVu2.colorMapLegend import ColorMapLegend
from DejaVu2.Geom import Geom
from DejaVu2.colorTool import RGBARamp, RedWhiteBlueRamp, RedWhiteRamp, \
     WhiteBlueRamp, TkColor
from DejaVu2.Legend import drawLegendOnly


class ComboBoxRename(Pmw.ComboBox):

    def setExternalList(self, externalList, reverse=False, numOfBlockedLabels=0):
        self.externalList = externalList
        self.reverse = reverse
        self.numOfBlockedLabels = numOfBlockedLabels


    def _addHistory(self):
        #print "_addHistory"
        lCurSelection = self.curselection()
        if len(lCurSelection) > 0:
            lCurSelection0 = int(lCurSelection[0])
            if self.reverse is True:
                lBoolTest = (len(self.externalList)-1-lCurSelection0) >= self.numOfBlockedLabels
            else:
                lBoolTest = lCurSelection0 >= self.numOfBlockedLabels
            if lBoolTest:
                input = self._entryWidget.get()
                if (input != '') and (input not in self._list.get(0, 'end') ):
                            self.delete(lCurSelection0)
                            self.insert(lCurSelection0, input)
                            if self.reverse is True:
                                self.externalList[len(self.externalList)-1-lCurSelection0] = input
                            else:
                                self.externalList[lCurSelection0] = input
            self.selectitem(lCurSelection0)
        else:
            self.selectitem(0)


class DummyEvent:
    """dummy Tkinter event """
    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y

COLORMAPGUI_MAXIMAL_HEIGHT = 1024

class ColorMapGUI(Tkinter.Frame, ColorMap):
    """The ColorMapGUI is an object providing an interface to expose
and to modify a ColorMap object.

The ColorMapGUI has a local copy of the colormap.ramp, colormap.mini and
colormap.maxi which will be modified through the GUI.
The associated colormap object will only be updated when by the apply_cb
which is called whenever a mouseUp event or update is called if the gui
is in continuous mode or when the Aplly button is pressed when the gui is
not in a continuous mode.

The GUI allows the user to edit each of 4 different color
properties: 'Hue','Sat' or saturation, 'Val' or value 
and 'Opa' or opacity for each entry in the ramp of the
ColorMap. This is implemented via 4 Canvases, one for
each color property. (Drawing on a canvas is done by 
pressing the left mouse button over the canvas and 
holding it down while moving the mouse.)
Drawing at a point 'x','y', on one of the canvases changes 
the value of the corresponding color property in this way: 
'y' specifies the entry in the ColorMap.ramp to be changed.
'x' is first normalized to give a value between 0 and 1.0. 
(by dividing it by the effective canvas width).
A new rgba value is built from the hsva list comprised of
3 previous values and this new normalized 'x' in the
appropriate position. (that is, if the user is drawing on
the 'Hue' canvas, the normalized 'x' value replaces the 'h',
or zeroth entry,etc). 

If the GUI has a viewer, an OpenGLColorMapWidget is constructed
which displays the GUI's current rgbMap. Also, the user can choose
to add a ColorMapLegend to the viewer which displays the current
ColorMap.ramp and which can be labeled.
"""

    staticRamp = RGBARamp(size=128, upperValue=.85).tolist()

    def __init__(self, cmap=None, master=None, continuous=0,  allowRename=True, modifyMinMax=False,
                 viewer=None, xoffset=25, width=200, height=256,
                 geoms=None, name=None, ramp=None, labels=None, mini=None, maxi=None,
                 filename=None, show=True, numOfBlockedLabels=0, **kw):

        #this is the active colormap
        if cmap is None:
            ColorMap.__init__(self, name=name, ramp=ramp, labels=labels,
                              filename=filename, mini=mini, maxi=maxi)
        else:
            ColorMap.__init__(self, name=cmap.name, ramp=cmap.ramp, labels=cmap.labels,
                              mini=cmap.mini, maxi=cmap.maxi)

        self.history = []
        self.legend = None
        if geoms is None:
            geoms = {}
        self.geoms = geoms
        self.viewer = None
        self.currentOnStack = False
        self.numOfBlockedLabels = numOfBlockedLabels
        
        # initialize a bunch of flags.
        self.cmapCurrent = True
        self.cmapCompToReset = None
        self.cmapResetAll = False
        self.nbStepBack = 0
        
        # Initialize the local guiRamp, guiMini and guiMaxi values
        # which will be modified
        self.guiRamp = []
        self.guiMini = None
        self.guiMaxi = None
        if self.labels is not None:        
            self.guiLabels = deepCopySeq(self.labels)
        self.colorbin = []
        self.labelbin = []
        self.lengthRamp = len(self.ramp)

        if height < self.lengthRamp:
            height = self.lengthRamp
            
        self.height = height
        self.linesPerRampValue = self.height / float(self.lengthRamp)
        self.currentHue = []
        self.currentSat = []

        # Create parts of the GUI
        if master is None:
            if viewer:
                theMaster = Tkinter.Toplevel(viewer.master)
            else:
                theMaster = Tkinter.Toplevel()
        else:
            theMaster = master

        if hasattr(theMaster, 'protocol'):
            theMaster.protocol('WM_DELETE_WINDOW', self.dismiss)

        # we hide it if it was asked
        if show is False:
            master2 = theMaster
            while hasattr(master2, 'withdraw') is False and hasattr(master2, 'master') is True:
                master2 = master2.master
            master2.withdraw()
        
        Tkinter.Frame.__init__(self, theMaster)
        Tkinter.Pack.config(self, expand=1, fill='both')
        if hasattr(theMaster, 'title'):
            theMaster.title(self.name)
        # Canvas width (default=200)
        self.width = width

        self.canvasStart = 4

        # left border width (default=25)
        self.xoffset = xoffset
        # Available drawing width default=width-xoffset
        self.xrange = float(width - xoffset)
        
        # DejaVu2 viewer to which colorMapLegend is added 
        self.viewer = None 

        self.continuousUpdate = Tkinter.IntVar()
        self.continuousUpdate.set(continuous)

        #if allowRename is False, name entry is disabled.
        self.allowRename = allowRename
        self.modifyMinMax = modifyMinMax
        
        self.idList = ['Hue','Sat','Val','Opa']
        self.getColorFunc = { 'Hue':self.hueColor,
                              'Sat':self.satColor,
                              'Val':self.valColor,
                              'Opa':self.opaColor
                              }

        #initialize dictionaries
        #NB all the lines end at xpoints corresponding to the
        #the colorproperty value at that point... (except
        #'Hue' is backwards, width-property+xoffset instead of
        #property+xoffset)

        self.rightXVals = {}
        self.lines = {}
        for idStr in self.idList:
            self.rightXVals[idStr] = []
            self.lines[idStr] = []

        self.current = 'Hue' # name of the currently visible prop. canvas
        self.callbacks = []  # list of functions to be called when apply
                             # button is pressed or ramp is modified and
                             # mode is continuous update

        # Calls the createWidgets method to create the 4 canvases
        # and the buttons.
        self.createWidgets(theMaster)
        self.configureGui(ramp=self.ramp, geoms=self.geoms, 
                          legend=self.legend, labels=self.labels,
                          mini=self.mini, maxi=self.maxi, **kw)

        # Call the update to initialize the gui with the given cmap values,
        # without configuring the cmap
        self.update( mini=self.mini, maxi=self.maxi,
                     ramp=deepCopySeq(self.ramp), cfgCmap=False)
        self.cmapCurrent = True

        #set current colorFunc, lines and values to 'Hue'
        self.getColor = self.getColorFunc['Hue']
        self.currentLines = self.lines['Hue']
        self.currentValues = self.rightXVals['Hue']

        Tkinter.Widget.bind(self, "<Enter>", self.enter_cb)
        Tkinter.Widget.bind(self, "<Leave>", self.leave_cb)

        #set-up canvas mouse bindings:
        for idStr in ['Hue', 'Sat', 'Val', 'Opa']:
            Tkinter.Widget.bind(self.canvas[idStr], "<ButtonPress-1>",
                                self.mouseDown)
            Tkinter.Widget.bind(self.canvas[idStr], "<Button1-Motion>",
                                self.mouseMotion)
            Tkinter.Widget.bind(self.canvas[idStr], "<ButtonRelease-1>",
                                self.mouseUp)

            Tkinter.Widget.bind(self.canvas[idStr], "<Motion>",
                                self.updateCurTk)

            Tkinter.Widget.bind(self.canvas[idStr], "<ButtonPress-3>",
                                self.mouseDownRight)
            Tkinter.Widget.bind(self.canvas[idStr], "<Button3-Motion>",
                                self.mouseMotionRight)
            Tkinter.Widget.bind(self.canvas[idStr], "<ButtonRelease-3>",
                                self.mouseUpRight)

            self.canvas[idStr].configure(cursor='cross')
            
        self.straightLine = None

        # create the cml
        if self.legend is None:
            self.createCML()

        self.SetViewer(viewer)

#        # fields created by the CMLwidget
#        self.legend.createOwnGui()
#        self.legend.hideOwnGui()

        if master is not None:
            master2 = theMaster
            if hasattr(master2, 'master') is True and isinstance(master2.master,Tkinter.Widget):
                master2 = master2.master
            if hasattr(master2, 'master') is True and isinstance(master2.master,Tkinter.Widget):
                master2 = master2.master
            #Tkinter.Widget.bind(master2, '<Configure>', self.configure_cb)
            self.inVision = True
        else:
            Tkinter.Widget.bind(self.fullWidgetFrame, '<Configure>', self.configure_cb)
            self.inVision = False


    def configure_cb(self, event=None):
        #print "configure_cb", event, dir (event)
        #print "configure_cb event.height", event.height
        
#        self.width = event.width
#        self.xrange = float(self.width - self.xoffset)
        
        if self.inVision is True:
            height = event.height - 2*self.canvasStart \
                     - self.menuFrame1.winfo_reqheight() \
                     - self.buttonFrame.winfo_reqheight() \
                     - 73
        elif os.name == 'nt': #sys.platform == 'win32':
            height = event.height - 2*self.canvasStart \
                     - self.menuFrame1.winfo_reqheight() \
                     - self.buttonFrame.winfo_reqheight() \
                     - self.frame2.winfo_reqheight() \
                     - 2
        else:
            height = event.height - 2*self.canvasStart \
                     - self.menuFrame1.winfo_reqheight() \
                     - self.buttonFrame.winfo_reqheight() \
                     - self.frame2.winfo_reqheight()
                 
        #print "configure_cb height", height
        if height >= self.lengthRamp and height > 0:
            self.height = height
            #print "self.height 1", self.height
            self.linesPerRampValue = self.height / float(self.lengthRamp)
            self.resizeCanvases(self.height)
            self.drawRampCanvases(drawRamp=True)


    def changeRampLength_cb(self, numOfRampValue):
        #print "changeRampLength_cb", numOfRampValue, len(self.guiRamp), len(self.ramp)
        #import traceback;traceback.print_stack()

        lenguiramp = len(self.guiRamp)
        if numOfRampValue < self.numOfBlockedLabels:
            self.numOfRampValues.set(lenguiramp)
            return
        elif numOfRampValue == lenguiramp:
            return
        elif numOfRampValue < lenguiramp:
            for i in range(numOfRampValue, lenguiramp):
                self.colorbin.append(self.guiRamp.pop())
        else: # numOfRampValue > lenguiramp:
            lNumOfColorsToAdd = numOfRampValue - lenguiramp
            lNumOfColorsToCreate = lNumOfColorsToAdd - len(self.colorbin)
            if lNumOfColorsToCreate > 0:
                for i in range(len(self.colorbin)):
                    self.guiRamp.append(self.colorbin.pop())    
                for i in range(lNumOfColorsToCreate):
                    lIndex = len(self.guiRamp) % len(self.staticRamp)
                    self.guiRamp.append( self.staticRamp[ lIndex ] )  
            else:
                for i in range(lNumOfColorsToAdd):
                    self.guiRamp.append(self.colorbin.pop())    

        self.configureGui( 
                          ramp=self.guiRamp,
                          mini=self.mini, 
                          maxi=self.maxi,
                          updateGui=True,
                          theGivenRampIsTheGuiRamp=True,
                          #guiRamp=self.guiRamp
                          )
        self.cmapCurrent=False


    def adjustGuiLabels(self, newNumOfGuiLabels):
        #print "adjustGuiLabels", newNumOfGuiLabels

        lenguiLabels = len(self.guiLabels)
        if newNumOfGuiLabels == lenguiLabels:
            return
        elif newNumOfGuiLabels < lenguiLabels:
            for i in range(newNumOfGuiLabels, lenguiLabels):
                self.labelbin.append(self.guiLabels.pop())
        else: 
            lNumOfLabelsToAdd = newNumOfGuiLabels - lenguiLabels
            lNumOfLabelsToCreate = lNumOfLabelsToAdd - len(self.labelbin)
            if lNumOfLabelsToCreate > 0:
                for i in range(len(self.labelbin)):
                    self.guiLabels.append(self.labelbin.pop())
                for i in range(lNumOfLabelsToCreate):
                    self.guiLabels.append( str(len(self.guiLabels) ) )
            else:
                for i in range(lNumOfLabelsToAdd):
                    self.guiLabels.append(self.labelbin.pop())
        self.configureGui(labels=self.guiLabels)


    def createCML(self):
        """create the ColorMap Legend
"""
        if self.legend:
            return
        #width defaults to len(ramp); height defaults to 1
        #mini, maxi and interp set from self
        cmlOptions = {'width':10, 'height':1,
                      'name':self.name, 'ramp':self.ramp}

        if self.mini is not None:
            cmlOptions['mini'] = self.mini

        if self.maxi is not None:
            cmlOptions['maxi'] = self.maxi

        self.legend = apply(ColorMapLegend, (self,), cmlOptions)
        

#    def resetComp(self, comp):
#        """resets on component to first entry in history list"""
#        compNum = ['Hue','Sat','Val','Opa'].index(comp)
#        # Get the first entry in the history list
#        hist = self.history[0]
#        # Remove all the entries in the history list
#        self.history = self.history[:1]
#        
#        # Transform the history RGBramp into a HSL ramp
#        histHSL = map(lambda x: list(RGBA2HSLA_list(x)), hist)
#        # Transform the current ramp in a HSL ramp
#        hsl = map(lambda x: list(x), self.asHSL())
#        for i in range(len(hsl)):
#            hsl[i][compNum] = histHSL[i][compNum]
#        ramp = map(lambda x: list(HSLA2RGBA_list(x)), hsl)
#        self.configure(ramp=ramp, updateGui=False)


    def reset(self):
        """return to first entry in history list"""
        ramp = self.history[0]
        self.history = []
        self.configure(ramp=deepCopySeq(ramp), updateGui=False)


    def pushRamp(self):
        """
        The pushRamp method appends the current ramp stored in cmap.ramp at the
        end of the history list.
        The attribute self.currentOnStack is set to True.
        """
        # append the current ramp to the history if the ramp is not None and
        # if the ramp has not been pushed onto the history stack yet.
        if self.ramp is None or self.currentOnStack:
            return
        self.history.append(deepCopySeq(self.ramp))
        self.currentOnStack = True


    def popRamp(self, index=-1):
        """
        optional arguments
        index -- default -1, specifies the index of the ramp to pop.
                 It is a negative index
                 as we are popping ramps from the end of the history list.
                 
        The popRamp method will save the ramp stored at the given index
        in cmap.ramp then remove all the entry from the ramp (included)
        till the end of the history list.
        After a popRamp the attribute self.currentOnStack is False.
        """
        # Need at least one entry in the history
        if not len(self.history): return

        # Always keep the first entry of the history 
        if len(self.history)==1:
            if not self.currentOnStack:
                self.ramp = self.history[0]
                self.currentOnStack = True
            return
        # PopRamp removes entry from the end of the history list
        # which is why the index has to be negative.
        assert (index < 0)

        # 1- set the actual ramp to be the entry of history
        # corresponding to the given index
        # 2- history = history[:index]
        # the popped history is not on the stack any longer except
        # when it is the first entry

        # Get the index from the beginning of the history list
        newind = len(self.history)+index
        pushRamp = False
        if newind <= 0:
            newind = 0
            # Always keep the first entry in the history list
            pushRamp = True
            
        elif newind >= len(self.history):
            newind = -1
            pushRamp = False

        ramp = self.history[newind]
        self.history = self.history[:newind]
        # We do not want to push the ramp we are popping onto the history stack
        # the current ramp will not be on the history stack any longer.
        self.configure(ramp=deepCopySeq(ramp), updateGui=False,
                       pushRamp=pushRamp)


    def SetViewer(self, viewer):
        """to give a viewer reference to the ColorMapGUI even
after it's creation
"""
        if viewer :
            if self.viewer is not None:
                
                self.viewer.RemoveObject(self.legend)

            self.viewer = viewer
            self.createCML()

            self.legend.name = self.viewer.ensureUniqueName(self.legend.name)
            self.update(cmapName=self.legend.name)

            self.legend.replace = False
            #self.legend.protected = True
            self.viewer.AddObject(self.legend, redo=0)   

            self.showLegendVar = Tkinter.IntVar()
            self.showLegendButton = Tkinter.Checkbutton(
                self.buttonFrame, text='show legend',
                variable=self.showLegendVar, command=self.showLegendButton_cb)       
            self.showLegendButton.grid(row=5, column=0, columnspan=2, sticky='w')
                        

    def createWidgets(self, master):
        """create Tkinter widgets: 4 canvas and buttons in 2 frames
"""
        #print "createWidgets"
        self.fullWidgetFrame = Tkinter.Frame(master)
        self.fullWidgetFrame.pack(side='top', expand=1, fill='both')

        # create menu frame
        self.menuFrame1 = Tkinter.Frame(self.fullWidgetFrame, relief='raised', borderwidth=3)
        self.menuFrame1.pack(side='top', expand=1, fill='x')
        
        filebutton = Tkinter.Menubutton(self.menuFrame1, text='File')
        filebutton.pack(side='left')
        filemenu = Tkinter.Menu(filebutton, {})
        filemenu.add_command(label='Read', command=self.read_cb)
        filemenu.add_command(label='Write', command=self.write_cb)
        filebutton['menu'] = filemenu

        editbutton = Tkinter.Menubutton(self.menuFrame1, text='Edit')
        editbutton.pack(side='left', anchor='w')
        editmenu = Tkinter.Menu(editbutton, {})
        editmenu.add_command(label='Reset to first in history', command=self.resetAll_cb)
        editmenu.add_command(label='Step back in history loop', command=self.stepBack_cb)
##         editmenu.add_checkbutton(label='Preview', variable=self.preview)
        editmenu.add_checkbutton(label='Continuous',
                                 variable=self.continuousUpdate)

        editmenu.add_command(label='Edit Legend', command=self.editLegend_cb)
                                                                                     
        editbutton['menu'] = editmenu

        self.numOfRampValues = ThumbWheel(
                                    self.menuFrame1,
                                    labCfg={'text':'Num of colors:', 'side':'left'},
                                    showLabel=1, 
                                    width=80,
                                    height=16,
                                    min=1, 
                                    max=COLORMAPGUI_MAXIMAL_HEIGHT,
                                    type=int,
                                    value=self.lengthRamp,
                                    callback=self.changeRampLength_cb,
                                    continuous=True,
                                    oneTurn=60,
                                    wheelPad=2
                                    )
        self.numOfRampValues.pack(side='right', anchor='e')

        # create canvas frames
        self.canvasesFrame = Tkinter.Frame(self.fullWidgetFrame)

        # create a frame for rproperties canvases
        self.propCanvasFrame = Tkinter.Frame(self.canvasesFrame)
        
        # create 4 canvases
        self.canvas = {}
        for idStr in self.idList:
            self.canvas[idStr] = Tkinter.Canvas(
                self.propCanvasFrame, relief='sunken', borderwidth=3,
                width=self.width, height=self.height)
        #pack Hue to start with
        self.canvas['Hue'].pack(side='left', expand=1, fill='both')
        self.propCanvasFrame.pack(side='left', expand=1, fill='both')

        self.ogl_cmw = OGLColorMapWidget( self.canvasesFrame, self,
                                          width=19, height=self.height,
                                        )
        
        self.canvasesFrame.pack(side = 'top', expand=1, fill='both')

        # create a frame for buttons
        self.buttonFrame = Tkinter.Frame(self.fullWidgetFrame, relief='ridge',
                                         borderwidth=3)

        # create radio buttons to switch between canvas
        self.currentCanvasVar = Tkinter.StringVar()
        self.buttonHue = Tkinter.Radiobutton(
            self.buttonFrame, text='Hue', value = 'Hue', width=8,
            indicatoron = 0, variable = self.currentCanvasVar,
            command=self.button_cb)
        self.buttonHue.grid(row=0, column=0, sticky='we')#pack(side='left')

        self.buttonSat = Tkinter.Radiobutton(
            self.buttonFrame, text='Sat.', value = 'Sat', width=8,
            indicatoron = 0, variable = self.currentCanvasVar,
            command=self.button_cb)
        self.buttonSat.grid(row=0, column=1, sticky='we')#pack(side='left')

        self.buttonVal = Tkinter.Radiobutton(
            self.buttonFrame, text='Lum.', value = 'Val', width=8,
            indicatoron = 0, variable = self.currentCanvasVar,
            command=self.button_cb)
        self.buttonVal.grid(row=0, column=2, sticky='we')#pack(side='left')

        self.buttonOpa = Tkinter.Radiobutton(
            self.buttonFrame, text='Opa.', value = 'Opa', width=8,
            indicatoron = 0, variable = self.currentCanvasVar,
            command=self.button_cb)
        self.buttonOpa.grid(row=0, column=3, sticky='we')#pack(side='left')
        
        self.currentCanvas = self.canvas['Hue']
        self.currentCanvasVar.set('Hue')

        # name
        self.nameVar = Tkinter.StringVar()
        if self.name is not None:
            self.nameVar.set(self.name)
        a = Tkinter.Label(self.buttonFrame, text='name')
        a.grid(row=2, column=0, sticky='e')
        if self.allowRename is True:        
            self.nameEntry = Tkinter.Entry(
                self.buttonFrame, textvariable=self.nameVar, width=18)
            self.nameEntry.bind('<Return>', self.rename_cb)
            self.nameEntry.bind('<Leave>', self.rename_cb)
        else:
            self.nameEntry = Tkinter.Label(
                self.buttonFrame, textvariable=self.nameVar, width=18,
                relief='groove', padx=2, pady=1)
        self.nameEntry.grid(row=2, column=1, columnspan=3, sticky='w')

        # min max
        self.maxTk = Tkinter.StringVar()
        if self.guiMaxi is not None:
            self.maxTk.set(('%8.2f'%self.guiMaxi))
        a = Tkinter.Label(self.buttonFrame, text='max')
        a.grid(row=3, column=2, sticky='e')
        if self.modifyMinMax:
            self.maxEntry = Tkinter.Entry(
                self.buttonFrame, textvariable=self.maxTk, width=8)
            self.maxEntry.bind('<Return>', self.max_cb)
            self.maxEntry.bind('<Leave>', self.max_cb)
        else:
            self.maxEntry = Tkinter.Label(
                self.buttonFrame, textvariable=self.maxTk, width=8,
                relief='groove', justify='left', padx=2, pady=1)
        self.maxEntry.grid(row=3, column=3,sticky='w')

        self.curyTk = Tkinter.StringVar()
        a = Tkinter.Label(self.buttonFrame, text='y ')
        a.grid(row=4, column=2, sticky='e')
        self.curYLabel = Tkinter.Label(
            self.buttonFrame, textvariable=self.curyTk, width=8,
            relief='groove', justify='left', padx=2, pady=1)
        self.curYLabel.grid(row=4, column=3, sticky='w')

        self.minTk = Tkinter.StringVar()
        if self.guiMini is not None:
            self.minTk.set(('%8.2f'%self.guiMini))
        a = Tkinter.Label(self.buttonFrame, text='min')
        a.grid(row=5, column=2, sticky='e')
        if self.modifyMinMax:
            self.minEntry = Tkinter.Entry(
                self.buttonFrame, textvariable=self.minTk, width=8)
            self.minEntry.bind('<Return>', self.min_cb)
            self.minEntry.bind('<Leave>', self.min_cb)
        else:
            self.minEntry = Tkinter.Label(
                self.buttonFrame, textvariable=self.minTk, width=8,
                relief='groove', justify='right', padx=2, pady=1)
        self.minEntry.grid(row=5, column=3, sticky='w')

        self.curxTk = Tkinter.StringVar()
        a = Tkinter.Label(self.buttonFrame, text='x')
        a.grid(row=4, column=1, sticky='w')
        self.curXLabel = Tkinter.Label(
            self.buttonFrame, textvariable=self.curxTk, width=8,
            relief='groove', justify='left', padx=2, pady=1)
        self.curXLabel.grid(row=4, column=0, sticky='w')

        if self.labels is not None:
            #self.labelsInComboBox = Tkinter.StringVar()
            self.labelsComboBox = ComboBoxRename(
                                       self.buttonFrame, 
                                       label_text='label',
                                       labelpos='w',
                                       history=1)
            self.labelsComboBox.setExternalList(
                                    self.guiLabels, 
                                    reverse=True, 
                                    numOfBlockedLabels=self.numOfBlockedLabels)
            self.labelsComboBox.grid(row=6, column=0, columnspan=4, sticky='w')

        #Apply and Dismiss go here
        self.frame2 = Tkinter.Frame(self.fullWidgetFrame, relief='raise', borderwidth=3)
        f2 = self.frame2
        self.apply = Tkinter.Button(f2, text='Apply', command=self.apply_cb)
        self.apply.pack(side='left', expand=1, fill='x')
        self.dismiss = Tkinter.Button(f2, text='Dismiss', command=self.dismiss)
        self.dismiss.pack(side='left', expand=1, fill='x')

        f2.pack(side='bottom', expand=1, fill='x')
        self.frame2.pack(side='bottom', expand=1, fill=Tkinter.X)
        self.buttonFrame.pack(side='bottom', expand=1, fill='x')


    def update(self, ramp=None, mini=None, maxi=None,
               cmapName=None, drawRamp=True, cfgCmap=True, 
               theGivenRampIsTheGuiRamp=False,
               #guiRamp=None
               ):

        # Update the name of cmg if relevant
        if hasattr(self.master, 'title'):
            self.master.title(cmapName)
            
        # Update the maxi and mini values of cmg
        if maxi is not None and mini is not None:
            if mini <= maxi:
                self.guiMaxi = maxi
                self.guiMini = mini
                self.minTk.set(('%8.2f'%self.guiMini))
                self.maxTk.set(('%8.2f'%self.guiMaxi))
            else:
                self.guiMaxi = None
                self.guiMini = None
                self.minTk.set('')
                self.maxTk.set('')
        elif maxi is not None and mini is None:
            if self.guiMini <= maxi:
                self.guiMaxi = maxi
                self.maxTk.set(('%8.2f'%self.guiMaxi))
            else:
                self.guiMaxi = None
                self.maxTk.set('')
        elif mini is not None and maxi is None:
            if mini <= self.guiMaxi:
                self.guiMini = mini
                self.minTk.set(('%8.2f'%self.guiMini))
            else:
                self.guiMini = None
                self.minTk.set('')

        # Update the cmg ramp with the one given
        if not ramp is None:
            if theGivenRampIsTheGuiRamp is False:
                self.cmapCurrent = False
                ramp = deepCopySeq(self.checkRamp(ramp))
                self.lengthRamp = len(ramp)
                if self.lengthRamp > self.height :
                    self.height = self.lengthRamp
                    #print "self.height 2", self.height
                self.numOfRampValues.set(self.lengthRamp)            
                self.linesPerRampValue = self.height / float(self.lengthRamp)
                self.resizeCanvases(self.height)
                self.guiRamp = ramp
                if self.labels is not None:
                    self.adjustGuiLabels(self.lengthRamp)
                self.drawRampCanvases(drawRamp=drawRamp)
            else:
                #print "guiRamp"
                self.lengthRamp = len(ramp)
                if self.lengthRamp > self.height :
                    self.height = self.lengthRamp
                    #print "self.height 3", self.height
                #self.numOfRampValues.set(self.lengthRamp)            
                self.linesPerRampValue = self.height / float(self.lengthRamp)
                self.resizeCanvases(self.height)
                #self.guiRamp = guiRamp
                if self.labels is not None:
                    self.adjustGuiLabels(self.lengthRamp)
                self.drawRampCanvases(drawRamp=drawRamp)

        # Update the OpenGL ramp
        self.ogl_cmw.ramp = deepCopySeq(self.guiRamp)
        self.ogl_cmw.tkRedraw()

        if self.continuousUpdate.get():
            self.configureCmap(cfgCmap=cfgCmap)


    def configure(self, name=None, ramp=None, geoms=None, legend=None, labels=None,
                  mini='not passed', maxi='not passed', viewer=None, updateGui=True,
                  pushRamp=True, **kw):
        #print "ColorMapGUI.configure", mini, maxi
        ColorMap.configure(self, name=name, ramp=ramp, labels=labels, mini=mini, maxi=maxi)

        if ramp is not None:
            ramp = self.ramp
        if labels is not None:
            labels = self.labels # because it was just set in the configure

        self.configureGui(ramp=ramp, labels=labels,
                          geoms=geoms, legend=legend,
                          viewer=viewer, updateGui=updateGui,
                          pushRamp=pushRamp, **kw)

   
    def configureGui(self, name=None, ramp=None, labels=None, 
                     geoms=None, legend=None,
                     viewer=None, updateGui=True,
                     pushRamp=True, 
                     #guiRamp=None, 
                     theGivenRampIsTheGuiRamp=False, 
                     **kw):
        """Configure the colormapGui with the given values.
"""        
        if (ramp is not None) and (theGivenRampIsTheGuiRamp is False):
            # The ramp is new put has not been pushed onto the history
            # stack yet
            self.currentOnStack = False
            # When pushRamp is True this new ramp is pushed on the stack
            if pushRamp:
                self.pushRamp()

        if labels is not None and self.labels is not None:
            #print "labels", labels
            #import traceback;traceback.print_stack()
            self.labelsComboBox.delete(0,'end')
            for label in labels:
                self.labelsComboBox.insert(0, str(label) )

        if geoms is not None and len(geoms):
            self.geoms.update(geoms)
        if viewer is not None and viewer != self.viewer:
            self.SetViewer(viewer)
            if self.name != self.legend.name:
                self.nameVar.set(self.legend.name)
                self.rename_cb()
        # if a legend is specified then set self.legend
        if legend is not None:
            self.legend = legend

        # Then update the legend with the given values
        if self.legend is not None:
            cmlOptions = {'ramp':self.ramp}

            # Need to update the legendValues as well
            cmlOptions['mini'] = self.mini
            cmlOptions['maxi'] = self.maxi
                
            if hasattr(self, 'showLegendVar'):
                cmlOptions['visible'] = self.showLegendVar.get()
            apply(self.legend.Set, (), cmlOptions)

        if updateGui:
            self.update(ramp=ramp, mini=self.mini, maxi=self.maxi, 
                        cfgCmap=False, 
                        theGivenRampIsTheGuiRamp=theGivenRampIsTheGuiRamp,
                        )
        

    def Map(self, values, mini='not passed', maxi='not passed'):

        col = ColorMap.Map(self, values, mini=mini, maxi=maxi)

        if self.legend is not None:
            cmlOptions = {'ramp':self.ramp}

            # Need to update the legendValues with the value that where just used
            cmlOptions['mini'] = self.lastMini     
            cmlOptions['maxi'] = self.lastMaxi

            if hasattr(self, 'showLegendVar'):
                cmlOptions['visible'] = self.showLegendVar.get()
            apply(self.legend.Set, (), cmlOptions)
        
        return col


    def configureCmap(self, cfgCmap=True):
        """ This method configures the associated cmap with the GUI new values.
"""
        if self.cmapCompToReset is not None:
            self.resetComp(self.cmapCompToReset)
            self.cmapCompToReset = None
            cfgCmap = False
            self.cmapCurrent = True

        if self.cmapResetAll is True:
            self.reset()
            self.cmapResetAll = False
            cfgCmap = False
            self.cmapCurrent = True

        if self.nbStepBack != 0:
            if self.nbStepBack == len(self.history):
                cfgCmap = False
                self.cmapCurrent = True
            self.popRamp(-self.nbStepBack)
            self.nbStepBack = 0
            
        if cfgCmap:
#            if self.allowRename: 
#                name = self.nametk.get()
#            else: 
#                name=None
                
            # donnot update the cmap if already current...
            if self.cmapCurrent:
                ramp=None
            else:
                ramp = deepCopySeq(self.guiRamp)
            self.configure(ramp=ramp, mini=self.guiMini, maxi=self.guiMaxi)

            if self.legend is not None:
                if self.legend.viewer:
                    self.tk.call(self.legend.viewer.currentCamera._w, 'makecurrent')
                    self.legend.RedoDisplayList()
                    self.legend.viewer.Redraw()

            self.cmapCurrent = True

        if self.labels is not None:
            self.labels = deepCopySeq(self.guiLabels)
            #print "self.labels", self.labels

        # Then will call the callbacks with the new cmap values.
        self.callCallbacks()


    def drawRampCanvases(self, event=None, drawRamp=True):
        """draw all ramp canvases
"""
        #print "drawRampCanvases"
        # update Ramp

        self.currentHue = range(self.lengthRamp)
        self.currentSat = range(self.lengthRamp)

        if drawRamp:
            for v in self.idList:
                self.deleteCanvasLines(v)
                self.setRightXVals(v)
            var = self.currentCanvasVar.get()
            self.currentValues = self.rightXVals[var]
            self.drawHue()
            self.drawSaturation()
            self.drawValue()
            self.drawOpacity()


    #################################################################
    ### MOUSE CALLBACKS
    #################################################################
    def mouseUp(self, event=None):
        self.cmapCurrent=False
        self.update(drawRamp=False)


    def mouseDown(self, event):
        j = self.canvasStart
        # canvas x and y take the screen coords from the event and translate
        # them into the coordinate system of the canvas object
        Y = min(j+self.height-1, event.y)
        x = min(self.width, event.x)
        x = max(self.xoffset, x)
        if Y < j:
            Y = j
        y = int ( (Y-j) / self.linesPerRampValue)
        if y < 0:
            y = 0
        self.startx = x
        self.startY = Y
        self.starty = y
        self.updateCurTk(event)
        c = self.currentCanvas
        lineIndex = self.lengthRamp-1-y
        line = self.currentLines[lineIndex]
        col, graybitmap = self.getColor(x, lineIndex) #y
        newline = c.create_rectangle( j, 
                                      y*self.linesPerRampValue+j,
                                      1+x,#+j, 
                                      (1+y)*self.linesPerRampValue+j, 
                                      outline='',
                                      fill=col, 
                                      stipple=graybitmap)
        self.currentLines[lineIndex] = newline
        #update self.guiRamp
        self.updateRGBMap(x, lineIndex)
        #update rightXVals
        self.currentValues[lineIndex] = x
        c.delete(line)


    def mouseMotion(self, event):
        j = self.canvasStart
        # canvas x and y take the screen coords from the event and translate
        # them into the coordinate system of the canvas object
        # x,y are float
        #x = self.canvasHue.canvasx(event.x)
        #y = self.canvasHue.canvasy(event.y)
        # event.x, event.y are same as x,y but int
        Y = min(j+self.height-1, event.y)
        x = min(self.width, event.x)
        x = max(self.xoffset, x)
        if Y < j:
            Y = j
        y = int ( (Y-j) / self.linesPerRampValue)
        if y < 0:
            y = 0
        elif self.startx is None:
            self.mouseDown(event)
        c = self.currentCanvas
        if self.startY == Y:
            lineIndex = self.lengthRamp-1-y
            #print "lineIndex 0 ...", lineIndex 
            line = self.currentLines[lineIndex]
            col, graybitmap = self.getColor(x, lineIndex)
            newline = c.create_rectangle( j,
                                          y*self.linesPerRampValue+j,
                                          1+x,#+j,
                                          (1+y)*self.linesPerRampValue+j,
                                          outline='',
                                          fill=col,
                                          stipple=graybitmap)
            self.currentLines[lineIndex] = newline
            c.delete(line)
            self.currentValues[lineIndex] = x
            self.updateRGBMap(x, lineIndex)
        else: # we need to interpolate for all y's between self.starty and y
            dx = x-self.startx
            dy = Y-self.startY
            rat = float(dx)/float(dy)
            if Y > self.startY:
                for Yl in range(self.startY, Y+1):
                    yl = int ( (Yl-j) / self.linesPerRampValue)
                    ddx = int(rat*(Yl-self.startY)) + self.startx
                    lineIndex = self.lengthRamp-1-yl
                    #print "lineIndex 1 ...", lineIndex 
                    line = self.currentLines[lineIndex]
                    col, graybitmap = self.getColor(ddx, lineIndex)
                    newline = c.create_rectangle( j,
                                                  yl*self.linesPerRampValue+j,
                                                  1+ddx,#+j,
                                                  (yl+1)*self.linesPerRampValue+j,
                                                  outline='',
                                                  fill=col,
                                                  stipple=graybitmap)
                    self.currentLines[lineIndex] = newline
                    c.delete(line)
                    self.currentValues[lineIndex] = ddx
                    self.updateRGBMap(ddx, lineIndex)
            else:
                for Yl in range(self.startY, Y-1, -1):
                    yl = int ( (Yl-j) / self.linesPerRampValue)
                    ddx = int(rat*(Yl-self.startY)) + self.startx
                    lineIndex = self.lengthRamp-1-yl
                    #print "lineIndex 2 ...", lineIndex 
                    line = self.currentLines[lineIndex]
                    col, graybitmap = self.getColor(ddx, lineIndex)
                    newline = c.create_rectangle( j,
                                                  yl*self.linesPerRampValue+j,
                                                  1+ddx,#+j,
                                                  (1+yl)*self.linesPerRampValue+j,
                                                  outline='',
                                                  fill=col,
                                                  stipple=graybitmap)
                    self.currentLines[lineIndex] = newline
                    c.delete(line)
                    self.currentValues[lineIndex] = ddx
                    self.updateRGBMap(ddx, lineIndex)
            self.startY = Y
            self.startx = x
            # this flushes the output, making sure that 
            # the rectangle makes it to the screen 
            # before the next event is handled
        self.update_idletasks()
        self.updateCurTk(event)


    def mouseDownRight(self, event):
        j = self.canvasStart
        Y = min(j+self.height-1, event.y)
        x = min(self.width, event.x)
        x = max(self.xoffset, x)
        if Y < j:
            Y = j
        y = int ( (Y-j) / self.linesPerRampValue)
        if y < 0:
            y = 0
        self.startx = x
        self.startY = Y
        self.starty = y


    def mouseMotionRight(self, event):
        j = self.canvasStart
        Y = min(j+self.height-1, event.y)
        x = min(self.width, event.x)
        x = max(self.xoffset, x)
        if Y < j:
            Y = j
        y = int ( (Y-j) / self.linesPerRampValue)
        if y < 0:
            y = 0
        elif self.startx is None:
            self.mouseDown(event)
        c = self.currentCanvas
        if self.straightLine is not None:
            c.delete(self.straightLine[0])
            c.delete(self.straightLine[1])

        dy = Y-self.startY
        if abs(dy) < (self.linesPerRampValue/2):
            return

        straightLineBlack = c.create_line( 
                       self.startx,#+j, 
                       (self.starty+.5)*self.linesPerRampValue+j, 
                       x,#+j,
                       (y+.5)*self.linesPerRampValue+j,
                       fill='black',
                       width=3)
        straightLineWhite = c.create_line( 
                       self.startx,#+j, 
                       (self.starty+.5)*self.linesPerRampValue+j, 
                       x,#+j,
                       (y+.5)*self.linesPerRampValue+j,
                       fill='white',
                       width=1)
        self.straightLine = [straightLineBlack, straightLineWhite]
        self.update_idletasks()
        self.updateCurTk(event)


    def mouseUpRight(self, event=None):
        j = self.canvasStart
        c = self.currentCanvas
        if self.straightLine is not None:
            c.delete(self.straightLine[0])
            c.delete(self.straightLine[1])
        Y = min(j+self.height-1, event.y)
        x = min(self.width, event.x)
        x = max(self.xoffset, x)
        if Y < j:
            Y = j
        y = int ( (Y-j) / self.linesPerRampValue)
        if y < 0:
            y = 0
        elif self.startx is None:
            self.mouseDown(event)
        dx = x-self.startx
        dY = Y-self.startY
        dy = y-self.starty
        if abs(dY) < (self.linesPerRampValue/2):
            return
        rat = float(dx)/float(dy)
        if y > self.starty:
            for yl in range(self.starty, y+1):
                ddx = int(rat*(yl-self.starty)) + self.startx
                lineIndex = self.lengthRamp-1-yl
                line = self.currentLines[lineIndex]
                col, graybitmap = self.getColor(ddx, lineIndex)
                newline = c.create_rectangle( j,
                                              yl*self.linesPerRampValue+j,
                                              1+ddx,#+j,
                                              (yl+1)*self.linesPerRampValue+j,
                                              outline='',
                                              fill=col,
                                              stipple=graybitmap)
                self.currentLines[lineIndex] = newline
                c.delete(line)
                self.currentValues[lineIndex] = ddx
                self.updateRGBMap(ddx, lineIndex)
        else:
            for yl in range(self.starty, y-1, -1):
                ddx = int(rat*(yl-self.starty)) + self.startx
                lineIndex = self.lengthRamp-1-yl
                line = self.currentLines[lineIndex]
                col, graybitmap = self.getColor(ddx, lineIndex)
                newline = c.create_rectangle( j,
                                              yl*self.linesPerRampValue+j,
                                              1+ddx,#+j,
                                              (yl+1)*self.linesPerRampValue+j,
                                              outline='',
                                              fill=col,
                                              stipple=graybitmap)
                self.currentLines[lineIndex] = newline
                c.delete(line)
                self.currentValues[lineIndex] = ddx
                self.updateRGBMap(ddx, lineIndex)
        self.starty = y
        self.startY = Y
        self.startx = x
        # this flushes the output, making sure that 
        # the rectangle makes it to the screen 
        # before the next event is handled
        self.update_idletasks()
        self.updateCurTk(event)
        self.cmapCurrent=False
        self.update(drawRamp=False)


    ###########################################################################
    ##   CALLBACK FUNCTIONS
    ###########################################################################
    def apply_cb(self, event=None, cfgCmap=True):
        """
"""
        #print "ColorMapGUI apply_cb"
        # Need to set self.guiMini, self.guiMaxi,from, to, from fill and to fill
        # to be set to the proper values.
        
        # Method which reconfigures the colormap associated to the gui
        self.configureCmap()

        
    def fileOpenAsk(self, idir=None, ifile=None, types=None,
                    title='Open'):
        if types==None: types = [ ('All files', '*') ]
        file = tkFileDialog.askopenfilename( filetypes=types,
                                             initialdir=idir,
                                             initialfile=ifile,
                                             title=title)
        if file=='': file = None
        return file
    

    def fileSaveAsk(self, idir=None, ifile=None, types = None,
                title='Save'):
        if types==None: types = [ ('All files', '*') ]
        file = tkFileDialog.asksaveasfilename( filetypes=types,
                                               initialdir=idir,
                                               initialfile=ifile,
                                               title=title)
        if file=='': file = None
        return file


    def button_cb(self, event=None):
        """call back function for the buttons allowing to toggle between
        different canvases.
        This function hides the currentCanvas and shows the canvas
        corresponding to the active radio button.
        In addition it sets
             self.currentCanvas : used to hide it next time we come in
             self.currentLines : list of Canvas Line objects (one per canvas)
             self.currentValues : list of numerical rightXVals (one per canvas)
             self.getColor =  function to be called to represent a color as a
                              Tk string
        """
        var = self.currentCanvasVar.get()

        if var == 'Hue':
            self.drawHue()
        elif var == 'Sat':
            self.drawSaturation()
        elif var == 'Val':
            self.drawValue()
        elif var == 'Opa':
            self.drawOpacity()

        newCanvas = self.canvas[var]
        self.currentCanvas.forget()
        newCanvas.pack(side='top')
        self.currentCanvas = newCanvas
        self.currentLines = self.lines[var]
        self.currentValues = self.rightXVals[var]
        self.getColor = self.getColorFunc[var]
        self.current = var


    def min_cb(self, event=None):
        minVal = self.minTk.get()
        try:
            minVal = float(minVal)
            if self.guiMaxi is not None and minVal >= self.guiMaxi:
                raise ValueError
            self.update(mini=minVal)
        except:
            self.minTk.set('')
            self.guiMini = None
            return  


    def max_cb(self, event=None):
        maxVal = self.maxTk.get()
        try:
            maxVal = float(maxVal)
            if self.guiMini is not None and maxVal <= self.guiMini:
                raise ValueError
            self.update(maxi=maxVal)
        except:
            self.maxTk.set('')
            self.guiMaxi = None
            return  


    def rename_cb(self, event=None):
        self.name = self.nameVar.get()
        self.update(cmapName=self.name)
        self.legend.Set(name=self.name)
        if self.legend.ownGui is not None:
            self.legend.ownGui.title(self.legend.name)
        if self.viewer is not None:
            self.viewer.Redraw()


    def quit(self):
        #print "quit"
        self.master.destroy()


    def dismiss(self):
        #print "dismiss"
        if self.master.winfo_ismapped():
            self.master.withdraw()


    def resetAll_cb(self, event=None):
        self.cmapResetAll = True
        self.update(ramp = deepCopySeq(self.history[0]), cfgCmap=False)


#    def reset_cb(self, event=None):
#        """reset only the currently visible property canvas"""
#        var = self.currentCanvasVar.get()
#        self.deleteCanvasLines(var)
#        ramp = deepCopySeq(self.history[0])
#        
#        self.lengthRamp = len(ramp)
#        if self.lengthRamp > self.height :
#            self.height = self.lengthRamp
#        self.numOfRampValues.set( self.lengthRamp )
#        
#        hsva = map(lambda x, conv=RGBA2HSLA_list: conv(x), ramp)
#        ind = self.idList.index(var)
#        #format?...list?....array??
#        newVals = numpy.array(hsva)[:,ind]
#        self.guiRamp = self.buildRGBMap(ind, newVals)
#        
#        self.cmapCompToReset = var
#        if var == 'Hue':
#            self.setRightXVals('Hue')
#            self.drawHue()
#        elif var == 'Sat':
#            self.setRightXVals('Sat')
#            self.drawSaturation()
#        elif var == 'Val':
#            self.setRightXVals('Val')
#            self.drawValue()
#        elif var == 'Opa':
#            self.setRightXVals('Opa')
#            self.drawOpacity()
#        if self.continuousUpdate.get():
#            self.configureCmap()
#        #self.mouseUp(None) # to call callbacks
#        self.currentValues = self.rightXVals[var]
#        self.ogl_cmw.ramp = self.guiRamp[:]
#        self.ogl_cmw.tkRedraw()
        

    def stepBack_cb(self):
        #update self.ramp to last ramp on history list
        if len(self.history)==1 and self.cmapCurrent is True:
            return
        
        if self.cmapCurrent and self.nbStepBack == 0:
            index = self.nbStepBack + 2
        else:
            index = self.nbStepBack + 1

        self.nbStepBack = index
        newind = len(self.history)-index
        cfgCmap = True

        if newind < 0:
            # tried to go back too far
            newind = 0
            self.nbStepBack = 0

        if newind == 0:
            # only valid for continuous mode.
            # want to go back to the first entry of the history which
            # stays so we don't need to update the cmap
            cfgCmap=False
        self.update(ramp=self.history[newind], cfgCmap=cfgCmap)


    def read_cb(self):
        fileTypes = [("ColorMap",'*_map.py'), ("any file",'*.*')]
        fileBrowserTitle = "Read Color Map"
        fileName  = self.fileOpenAsk(types=fileTypes,
                                     title=fileBrowserTitle)
        if not fileName:
            return

        self.read(fileName)


    def read(self, filename):

        ColorMap.read(self, filename)
        self.configureGui(ramp=self.ramp, mini=self.mini, maxi=self.maxi)
        
        if self.legend.viewer:
            #self.legend.RedoDisplayList()
            self.legend.viewer.Redraw()
            #self.configureCmap()

    def enter_cb(self, event=None):
        if widgetsOnBackWindowsCanGrabFocus is False:
            lActiveWindow = self.focus_get()
            if    lActiveWindow is not None \
              and ( lActiveWindow.winfo_toplevel() != self.winfo_toplevel() ):
                return

        self.focus_set()


    def leave_cb(self, event=None):
        pass


    ############################################################
    ### UTILITY FUNCTIONS
    ############################################################
    def makeEvent(self, x, y):
        xval = x*self.xrange + self.xoffset
        return DummyEvent( xval, y)


    def hueColor(self, val, rampIndex):
        """ TkColorString <- hueColor(val)
val is an integer between self.xoffset and self.width.
(default:val is an integer between 25 and 200.)
returns color to be use to draw hue lines in hue canvas
"""
        #print "hueColor", h
        if val > self.width:
            h = 1
        else:
            h = (float(self.width-val)/self.xrange) #* .999 #*.6666666667
#        if h == 0:
#            h = 1
        self.currentHue[rampIndex] = h
        rgb = HSL2RGB(h, 1., .5)
        graybitmap = ''
        return TkColor(rgb), graybitmap


    def satColor(self, val, rampIndex):
        """ TkColorString <- satColor(val)
val is an integer between self.xoffset and self.width.
(default:val is an integer between 25 and 200.)
returns color to be use to draw saturation lines in saturation canvas
"""
        #print "satColor", val, rampIndex
        if val < self.xoffset:
            x = 0
        else:
            x = float(val-self.xoffset)/self.xrange
        self.currentSat[rampIndex] = x
        rgb = HSL2RGB(self.currentHue[rampIndex], x, .5)
        rgb.append(1.)
        graybitmap = ''
        return TkColor(rgb), graybitmap


    def valColor(self, val, rampIndex):
        """ TkColorString <- valColor(val)
val is an integer between 25 and 200
returns color to be used to draw value and opacity lines in canvas
"""
        #print "valColor", val, rampIndex
        if val < self.xoffset:
            x = 0
        else: 
            x = float(val-self.xoffset)/self.xrange
        rgb = HSL2RGB(self.currentHue[rampIndex], self.currentSat[rampIndex], x)
        rgb.append(1.)
        graybitmap = ''
        return TkColor(rgb), graybitmap


    def opaColor(self, val, rampIndex=None):
        #print "opaColor"       
        if val < self.xoffset:
            x = 0
        else:
            x = float(val-self.xoffset)/self.xrange

        if x >= 1:
            graybitmap = ''
            col = 'black' 
        elif x > .75:
            graybitmap = 'gray75'
            col = 'black'
        elif x > .50:
            graybitmap = 'gray50'
            col = 'black' 
        elif x > .25:
            graybitmap = 'gray25'
            col = 'black'
        elif x > 0:
            graybitmap = 'gray12'
            col = 'black' 
        else:
            graybitmap = ''
            col=''
        return col, graybitmap


    def opaColor_Unused(self, val, rampIndex=None):
        # the color of the dot changes 
        # we don't use it because the limit line drawn is less visible

        x = float(val- self.xoffset)/self.xrange
        background = self.canvas['Opa'].configure()['background']
        backgroundHSV = ToHSV(background[4], mode='HEX', flag255=0)
        bgHSV = [backgroundHSV[0], backgroundHSV[1], backgroundHSV[2] * (1-x)]
        col = ToHEX(bgHSV, mode='HSV')
        if x >= 1:
            graybitmap = ''
        elif x > .75:
            graybitmap = 'gray75'
        elif x > .50:
            graybitmap = 'gray50'
        elif x > .25:
            graybitmap = 'gray25'
        elif x > 0:
            graybitmap = 'gray12'
        else:
            graybitmap = ''
        return col, graybitmap


#    def createCMLWidget(self):
#        self.legend_gui = Tkinter.Toplevel()
#        self.legend_gui.title(self.legend.name)
#        self.legend_gui.protocol('WM_DELETE_WINDOW', self.legend_gui.withdraw )
#
#        frame1 = Tkinter.Frame(self.legend_gui)
#        frame1.pack(side='top')
#
#        #unit
#        self.unitsEnt = Pmw.EntryField(frame1, 
#                                       label_text='Units  ',
#                                       labelpos='w',
#                                       command=self.updateCML_cb)
#        self.unitsEnt.pack(side='top', fill='x')
#
#        #glf vector font
#        self.glfFont = Tkinter.StringVar()
#        self.glfFont.set('chicago1.glf')
#        self.glfFontCB = Pmw.ComboBox(frame1, label_text='Font    ',
#                                   labelpos='w',
#                                   entryfield_value=self.glfFont.get(),
#                                   scrolledlist_items=ColorMapLegend.glfVectorFontList,
#                                   selectioncommand=self.updateCML_cb)
#        self.glfFontCB.pack(side='top', fill='x')
#
#        #fontScale
#        self.fontScale = ThumbWheel(frame1,
#                                    labCfg={'text':'font scale            ', 'side':'left'},
#                                    showLabel=1, 
#                                    width=90,
#                                    height=14,
#                                    min=0, 
#                                    max=200,
#                                    type=int, 
#                                    value=self.legend.fontScale,
#                                    callback=self.updateCML_cb,
#                                    continuous=True,
#                                    oneTurn=10,
#                                    wheelPad=0)
#        self.fontScale.pack(side='top')
#
#        #label
#        self.labelValsEnt = Pmw.EntryField(
#                                frame1, 
#                                label_text='numpy labels    ',
#                                labelpos='w',
#                                command=self.updateCML_cb)
#        self.labelValsEnt.component('entry').config(width=6)
#        self.labelValsEnt.pack(side='top', fill='x')
#
#        #numOfLabel
#        self.numOfLabelsCtr = ThumbWheel(frame1,
#                                    labCfg={'text':'Automatic labels', 'side':'left'},
#                                    showLabel=1, 
#                                    width=90,
#                                    height=14,
#                                    min=0, 
#                                    max=200,
#                                    type=int, 
#                                    value=5,
#                                    callback=self.updateCML_cb,
#                                    continuous=True,
#                                    oneTurn=20,
#                                    wheelPad=0)
#        self.numOfLabelsCtr.pack(side='top')
#
#        # Interpolate
#        self.interpVar = Tkinter.IntVar()
#        self.interpVar.set(0)
#        self.checkBoxFrame = Tkinter.Checkbutton(
#                                frame1, 
#                                text='Interpolate',
#                                variable=self.interpVar, 
#                                command=self.updateCML_cb)
#        self.checkBoxFrame.pack(side='top')
#
#        # frame
#        self.frameVar = Tkinter.IntVar()
#        self.frameVar.set(1)
#        self.checkBoxFrame = Tkinter.Checkbutton(
#                                frame1, 
#                                text='Frame',
#                                variable=self.frameVar, 
#                                command=self.updateCML_cb)
#        self.checkBoxFrame.pack(side='top')
#
#        # invert labels color
#        self.invertLabelsColorVar = Tkinter.IntVar()
#        self.invertLabelsColorVar.set(0)
#        self.checkBoxinvertLabelsColor = Tkinter.Checkbutton(
#                                frame1, 
#                                text='Invert labels color',
#                                variable=self.invertLabelsColorVar, 
#                                command=self.updateCML_cb)
#        #self.checkBoxFrame.pack(side='top')
#        self.checkBoxinvertLabelsColor.pack(side='top')
#
#        # colormapguiwidget:
#        self.launchColormapWidget = Tkinter.Button(
#                                        frame1, 
#                                        text="Show colormap settings",
#                                        command=self.showColormapSettings_cb )
#        self.launchColormapWidget.pack(side='top', fill='x')


    def showColormapSettings_cb(self, event=None):
        #print "showColormapSettings_cb"
        master = self.master
        while hasattr(master, 'deiconify') is False and hasattr(master, 'master') is True:
            master = master.master
        if hasattr(master, 'deiconify'):
            if master.winfo_ismapped() == 0:
                master.deiconify()
            master.lift()
            #else: master.withdraw()


    def showLegendButton_cb(self, event=None):
        #print "showLegendButton_cb"
        if not self.viewer: 
            return
        if not self.legend:
            self.createCML()
        visible = self.showLegendVar.get()
        self.legend.Set(visible=visible)
        #self.legend.setWithOwnGui()

        if visible == 0:
            if self.legend.viewer.currentObject is self.legend:
                self.legend.viewer.SetCurrentObject(self.legend.viewer.rootObject)

        #do you have to do a redraw to see it?
        self.viewer.Redraw()


    def showLegend(self):
        self.showLegendVar.set(1)
        self.showLegendButton_cb()


    def hideLegend(self):
        self.showLegendVar.set(0)
        self.showLegendButton_cb()


    def editLegend_cb(self, event=None):
        #print "editLegend_cb"
        if not self.legend: 
            self.createCML()
        self.legend.showOwnGui()


    def resizeCanvases(self, height=None, width=None):
        #print "resizeCanvases"
        if height is not None:
            for c in self.canvas.values():
                c.configure(height=height)
            self.ogl_cmw.configure(height=height)
            self.ogl_cmw.height = height

        if width is not None:
            for c in self.canvas.values():
                c.configure(width=width)
            self.ogl_cmw.configure(width=width)
            self.ogl_cmw.width = width


    def buildRGBMap(self, ind, newVals):
        """ ind is 0 for Hue, 1 for Sat etc
newVals is list of new rightXValues for that column
"""
        #print "buildRGBMap"
        rgbs = []
        hsva = map(list, map(RGBA2HSLA_list, self.guiRamp))
        for i in range(len(self.guiRamp)):
            hsva[i][ind] = newVals[i]
            rgbs.append(list(HSLA2RGBA_list(hsva[i])))
        return rgbs
            

    def buildRightXVals(self):
        """ idList = ['Hue','Sat','Opa','Val']
"""
        for id in self.idList:
            self.setRightXVals(id)

       
    def setRightXVals(self, idStr):
        #each canvas line is drawn at constant y value
        #between x1=0 on the left and x2 which is current
        #ramp at y (which ranges between 0. and 1.) times 175
        #NB: 175 is the effective width of the canvas
        #idList = ['Hue','Sat','Opa','Val']

        #print "setRightXVals"

        if idStr not in self.idList:
            return
        ind = self.idList.index(idStr)
        #apparently these are 200-this value and
        #for entire Hue ramp delta x = .66401
        hsvaRamp = map(list, map(RGBA2HSLA_list, self.guiRamp))
        ramp_vals = numpy.array(hsvaRamp)[:,ind]
        norm_n = ramp_vals
        if idStr=='Hue':
            # self.width is 200 which is 25 border + 175 map
            # 25 and 175 are the default values but can be changed
            val = self.width - (norm_n * self.xrange)
        else:
            val = norm_n * self.xrange + self.xoffset
        #self.rightXVals[idStr] = map(int, val.tolist())
        self.rightXVals[idStr] = val.tolist()


    def write_cb(self):
        fileTypes = [("ColorMap",'*_map.py'), ("any file",'*.*')]
        fileBrowserTitle = "Write Color Map"
        fileName  = self.fileSaveAsk(types=fileTypes,
                                     title=fileBrowserTitle)
        if not fileName:
            return
        self.write(fileName)
        

    def deleteCanvasLines(self, name):
        """delete all lines in ramp canvas"""
        assert name in self.idList
        #assert name in ['Hue', 'Sat', 'Val', 'Opa' ]
        l = self.lines[name]
        c = self.canvas[name]
        for i in range(len(l)):
            c.delete(l[i])
        self.lines[name] = []
        self.rightXVals[name] = []
        self.curxTk.set("0.0")
        self.curyTk.set("0.0")
        if self.current==name:
            self.currentValues = self.rightXVals[name]
            self.currentLines = self.lines[name]


    def drawHue(self):
        j = self.canvasStart
        #print "drawHue"
        #this is for initializing and resets
        l = self.lines['Hue']
        c = self.canvas['Hue']
        v = self.rightXVals['Hue']

        for i in range(len(l)):
            c.delete(l[i])
        self.lines['Hue'] = []
        l = self.lines['Hue']
        
        for i in range(self.lengthRamp):
            #nb val is the LENGTH of the line to draw
            val = v[i]
            col, graybitmap = self.hueColor(val, i)
            l.append( c.create_rectangle( j,
                                          (self.lengthRamp-1-i)*self.linesPerRampValue+j,
                                          1+val,
                                          ((self.lengthRamp-1-i)+1)*self.linesPerRampValue+j,
                                          fill=col,
                                          outline='') )

        if self.current=='Hue':
            self.currentLines = l


    def drawSaturation(self):
        #print "drawSaturation"     
        c = self.canvas['Sat']
        l = self.lines['Sat']
        v = self.rightXVals['Sat']
        j = self.canvasStart
        
        for i in range(len(l)):
            c.delete(l[i])
        self.lines['Sat'] = []
        l = self.lines['Sat']
        
        for i in range(self.lengthRamp):
            val = v[i]
            col, graybitmap = self.satColor(val, i)
            l.append( c.create_rectangle( j,
                                          (self.lengthRamp-1-i)*self.linesPerRampValue+j,
                                          1+val,
                                          ((self.lengthRamp-1-i)+1)*self.linesPerRampValue+j,
                                          fill=col,
                                          outline='') )
        if self.current=='Sat':
            self.currentLines = l


    def drawValue(self):
        c = self.canvas['Val']
        l = self.lines['Val']
        v = self.rightXVals['Val']
        j = self.canvasStart
        
        for i in range(len(l)):
            c.delete(l[i])
        self.lines['Val'] = []
        l = self.lines['Val']
        
        for i in range(self.lengthRamp):
            val = v[i]
            col, graybitmap = self.valColor(val, i)
            #print "val i col", val, i, col
            l.append( c.create_rectangle( j,
                                          (self.lengthRamp-1-i)*self.linesPerRampValue+j,
                                          1+val,
                                          ((self.lengthRamp-1-i)+1)*self.linesPerRampValue+j,
                                          fill=col,
                                          outline='') )
        if self.current=='Val':
            self.currentLines = l


    def drawOpacity(self):
        c = self.canvas['Opa']
        l = self.lines['Opa']
        v = self.rightXVals['Opa']
        j = self.canvasStart

        for i in range(len(l)):
            c.delete(l[i])
        self.lines['Opa'] = []
        l = self.lines['Opa']
                
        for i in range(self.lengthRamp):
            val = v[i]
            col, graybitmap = self.opaColor(val)
            l.append( c.create_rectangle( j,
                                          (self.lengthRamp-1-i)*self.linesPerRampValue+j,
                                          1+val,
                                          ((self.lengthRamp-1-i)+1)*self.linesPerRampValue+j,
                                          outline='',
                                          fill=col,
                                          stipple=graybitmap) )
        if self.current=='Opa':
            self.currentLines = l


    def addCallback(self, function):
        assert callable(function)
        self.callbacks.append( function )


    def callCallbacks(self):
        for f in self.callbacks:
            f( self )
        

    def valueAtY(self, y):
        # compute the value for q given line in canvas
        y = y - self.canvasStart
        y = max(y, 0)
        y = min(y, self.height-1)
        if self.guiMini is not None and self.guiMaxi is not None:
            range = self.guiMaxi-self.guiMini
            return (float(self.height-y)/(self.height-1))*range + self.guiMini
        else:
            return int((self.height-y)/self.linesPerRampValue)


    def indexAtY(self, y):
        # compute the value for q given line in canvas
        j = self.canvasStart
        Y = min(j+self.height-1, y)
        if Y < j:
            Y = j
        y = int ( (Y-j) / self.linesPerRampValue)
        return y


    def updateCurTk(self, event):
        x = event.x# - self.canvasStart
        x = max(x, self.xoffset)
        x = min(x, self.width)
        val = (x-self.xoffset)/self.xrange
        self.curxTk.set("%4.2f"%val)
        lVal = self.valueAtY(event.y)
        if lVal is not None:
            self.curyTk.set("%4.2f"%lVal)
        else:
            self.curyTk.set('')
        if self.labels is not None:            
            self.labelsComboBox.selectitem(self.indexAtY(event.y))


    def updateOGLwidget(self):
        self.ogl_cmw.tkRedraw()

        
    def updateRGBMap(self, x, y):

        cur_val = list(RGBA2HSLA_list(self.guiRamp[y]))
        cur_val[0] = self.currentHue[y]
        cur_val[1] = self.currentSat[y]

        if self.current=='Hue':
            h = ((self.width-x)/self.xrange) * .999 #*.6666666667
            cur_val[0] = round(h, 3)
        else:
            #use keyList = ['Sat', 'Val', 'Opa']
            keyList = self.idList[1:]
            ind = keyList.index(self.current) + 1
            h = (x-self.xoffset)/self.xrange
            cur_val[ind] = round(h, 3)
        self.guiRamp[y] = list(HSLA2RGBA_list(cur_val))



class OGLColorMapWidget(OGLWidget):
    """This object provides a OpenGL window that displays an RGBA color map"""

    def __init__(self, master, cmap, title='OGLColorMapWidget', width=None,
                 height=None, cnf=None, **kw):

        if width is None:
            kw['width'] = 19
        else:
            kw['width'] = width
        if height is None:
            kw['height'] = 180
        else:
            kw['height'] = height

        self.callback = None

        assert isinstance(cmap, ColorMap)
        self.ramp = deepCopySeq(cmap.ramp)
        self.step = 1 # set to larger values to accelerate redraw
        
        tf = self.topFrame = Tkinter.Frame(master, borderwidth=1)
        self.frame = Tkinter.Frame(tf, borderwidth=3, relief='sunken')

        if cnf is None:
            cnf = {}
        apply( OGLWidget.__init__, (self, self.frame, title, cnf, 0), kw )

        self.frame.pack(expand=1,fill='both')
        self.topFrame.pack(expand=1,fill='both')
        self.pack(side='left',expand=1,fill='both')
        #self.initTexture()
        
 
    def tkRedraw(self, *dummy):
        """guillaume: probably useless, the code is almost identical
in the overriden function OGLWidget.tkRedraw
"""
        self.update_idletasks()
        self.tk.call(self._w, 'makecurrent')
        self.initProjection()
        GL.glPushMatrix()
        apply( self.redraw, dummy ) #guillaume: self.redraw() in OGLWidget.tkRedraw
        GL.glFlush()
        GL.glPopMatrix()
        self.tk.call(self._w, 'swapbuffers')

       
    def redraw(self, line=None):
        #print "redraw"
        drawLegendOnly(
                    fullWidth=self.width,
                    fullHeight=self.height,
                    ramp=self.ramp,
                    verticalLegend=True,
                    roomLeftToLegend=0,
                    roomBelowLegend=0,
                    legendShortSide=self.width,
                    legendLongSide=self.height,
                    interpolate=False,
                    selected=False,
                    )
        return
        



if __name__ == '__main__':
    #import pdb

    #test = ColorMapGUI()
    #def cb(ramp, min, max):
        #print len(ramp), min, max
    #test.addCallback(cb)

    ##cm = ColorMap('cmap', filename='/home/rhuey/python/dev/rgb.map')
    #cmg = ColorMapGUI(cm, allowRename=1)
    
    l = {}
    g = {}
    execfile('Tests/rgb256_map.py', g,l)
    cm = None
    #for name, object in l.items():
    #    if isinstance(object, ColorMap):
    #        cm = object
    #        break
    cm = l['cm']

    from DejaVu2 import Viewer
    vi = Viewer()
    cmg = ColorMapGUI(cm, allowRename=1, viewer=vi )

#test.mainloop()
