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
# Date: April 2003 Authors: Michel Sanner
#
#    vareille@scripps.edu
#    sanner@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Michel Sanner and TSRI
#
# revision: Guillaume Vareille
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/DejaVu2/colorMapLegend.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#
# $Id: colorMapLegend.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#

import os
import types
import Tkinter
import Pmw
from weakref import ref 
from copy import deepcopy
import string

from mglutil.gui.BasicWidgets.Tk.thumbwheel import ThumbWheel
from viewerFns import checkKeywords
from opengltk.OpenGL import GL
from DejaVu2.IndexedGeom import IndexedGeom
from colorTool import RGBRamp, resetMaterialMemory
from DejaVu2.Insert2d import Insert2d
import viewerConst
from pyglf import glf


class ColorMapLegend(Insert2d):
    """Class for displaying a colormap legend.
    arguments for the constructor or Set method are:
       ramp: numpy array of shape Nx3 or Nx4 (default is RBGRamp())
       height: height of colormap legend
       width:  width of colormap legend
       interp: 1 or 0 to turn color interpolation on and off resp.
       mini:   minimum values (i.e which corresponds to the first color)
       maxi:   maximum values (i.e which corresponds to the last color)
    
    If interp is set to 1 QUAD_STRIPS are used and colors are interpolated,
    else one QUAD is drawn for each entry in the colormap.
    If the color provide alpha values, a chechered background is drawn.
    """

    keywords = Insert2d.keywords + [
        'ramp',
        'height',
        'width',
        'interp',       # 1: for color interpolation, or 0
        'mini',         # 
        'maxi',         #
        'labelValues',  # floating point numbers to be written below cml
        'glfFont',
        'fontScale',
        'numOfLabels',
        'unitsString',
        'interp',
        'visibleFrame',
        'invertLabelsColor',
        ]    

    glfVectorFontList = [
        'arial1.glf',
        'courier1.glf',
        'crystal1.glf',
        'techno0.glf',
        'techno1.glf',
        'times_new1.glf',
        'aksent1.glf',
        'alpine1.glf',
        'broadway1.glf',
        'chicago1.glf',
        'compact1.glf',
        'cricket1.glf',
        'garamond1.glf',
        'gothic1.glf',
        'penta1.glf',
        'present_script1.glf'
    ]


    def __init__(self, colormapgui, name='color map legend', check=1, **kw):

        # GLF FONTS Initialisations
        glf.glfInit()
        glf.glfEnable(glf.GLF_CONSOLE_MESSAGES)
        lGlfModulePath = os.path.split(glf.__file__)[-2]
        lPathToFonts = lGlfModulePath+os.sep+'fonts'+os.sep
        self.glfVectorFontLoadIdDict = {}
        for font in self.glfVectorFontList:
            self.glfVectorFontLoadIdDict[font] = glf.glfLoadFont(lPathToFonts+font)   
        self.fontScale = 8
        self.glfFontID = 0

        # other initialisations
        self.colormapguiRef = ref(colormapgui) 
        self.resizeSpotRadius = 5 
        self.resizeSpot = None
        self.verticalLegend = True

        self.mini = None
        self.maxi = None
        
        kw['ramp'] = RGBRamp()
        kw['height'] = 1.
        kw['width'] = len(kw['ramp'])
        kw['interp'] = 1
#        kw['mini'] = None
#        kw['maxi'] = None
        kw['labelValues'] = []
        kw['numOfLabels'] = 5
        kw['unitsString'] = 'law'
        kw['invertLabelsColor'] = False
        kw['visibleFrame'] = True

        #kw['protected'] = True
        kw['immediateRendering'] = False
        kw['visible'] = False
        kw['transparent'] = True

        kw['size'] = [12, 120] # override the default value in Insert2d

        self.clickPosFromLegendBottomLeft = [0, 0]

        apply( Insert2d.__init__, (self, name, check), kw)

        # Insert2d initialisations 
        # (need to be done after the Insert2d.__init__ otherwise it overrides)
        self.needsRedoDpyListOnResize = True
        self.initialPosition = [1, 350] # override the default value in Insert2d
        self.coord2d = deepcopy(self.initialPosition) # 2d coordinates in pixels from top left
        self.animatable = False
        self.hiddenInCamera = {} # key:Camera object, value dummy


    def Set(self, check=1, redo=1, updateOwnGui=True, **kw):
        """set data for this object: add faces (polygon or lines) to this object
check=1 : verify that all the keywords present can be handle by this func 
redo=1 : append self to viewer.objectsNeedingRedo
updateOwnGui=True : allow to update owngui at the end this func
"""
        #print "colorMapLegend.Set"
        redoFlags = apply( Insert2d.Set, (self, check, 0), kw)
        
        ramp=kw.get('ramp')
        if ramp is not None:
            assert len(ramp) > 0
            assert len(ramp[0])==3 or len(ramp[0])==4
            self.ramp = ramp
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        height = kw.get('height')
        if height:
            assert height > 0.0
            self.height = height
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        width = kw.get('width')
        if width:
            assert width > 0.0
            self.width = width
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        mini = kw.get('mini')
        if mini is not None:
            #assert isinstance(mini, types.FloatType)
            self.mini = mini
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        maxi = kw.get('maxi')
        if maxi is not None:
            #assert isinstance(maxi, types.FloatType)
            self.maxi = maxi
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        labelValues = kw.get('labelValues')
        if labelValues is not None:
            # for v in labelValues:
            #    assert isinstance(v, types.FloatType)
            self.labelValues = labelValues
            redoFlags |= self._redoFlags['updateOwnGuiFlag']
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        glfFont = kw.get('glfFont')
        if glfFont is not None:
            if glfFont in self.glfVectorFontLoadIdDict.keys():
                self.glfFontID = self.glfVectorFontLoadIdDict[glfFont]
                redoFlags |= self._redoFlags['updateOwnGuiFlag']
                redoFlags |= self._redoFlags['redoDisplayListFlag']

        fontScale = kw.get('fontScale')
        if fontScale is not None:
            assert isinstance(fontScale, types.IntType)
            self.fontScale = fontScale
            redoFlags |= self._redoFlags['updateOwnGuiFlag']
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        numOfLabels = kw.get('numOfLabels')
        if numOfLabels is not None:
            assert isinstance(numOfLabels, types.IntType)
            self.numOfLabels = numOfLabels
            redoFlags |= self._redoFlags['updateOwnGuiFlag']
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        unitsString = kw.get('unitsString')
        if unitsString is not None:
            self.unitsString = unitsString
            redoFlags |= self._redoFlags['updateOwnGuiFlag']
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        interp = kw.get('interp')
        if interp is not None:
            assert interp in (0, 1)
            self.interp = interp
            redoFlags |= self._redoFlags['updateOwnGuiFlag']
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        visibleFrame = kw.get('visibleFrame')
        if visibleFrame is not None:
            self.visibleFrame = visibleFrame
            redoFlags |= self._redoFlags['updateOwnGuiFlag']
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        invertLabelsColor = kw.get('invertLabelsColor')
        if invertLabelsColor is not None:
            self.invertLabelsColor = invertLabelsColor
            redoFlags |= self._redoFlags['updateOwnGuiFlag']
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        return self.redoNow(redo, updateOwnGui, redoFlags)


    def Draw(self):
        #print "colorMapLegend.Draw", self

        if self.viewer.tileRender:
            tile = self.getTile()
            #print "tile", tile
        else:
            tile = None

        fullWidth = self.viewer.currentCamera.width
        fullHeight = self.viewer.currentCamera.height

        if self.invertLabelsColor is False:
            backgroundColor = (
                           self.viewer.currentCamera.backgroundColor[0],
                           self.viewer.currentCamera.backgroundColor[1],
                           self.viewer.currentCamera.backgroundColor[2],
                           .5)
        else:
            backgroundColor = (
                           1-self.viewer.currentCamera.backgroundColor[0],
                           1-self.viewer.currentCamera.backgroundColor[1],
                           1-self.viewer.currentCamera.backgroundColor[2],
                           .5)

        from DejaVu2.Legend import drawSelfOrientedLegend
        self.polygonContour , self.resizeSpot , self.verticalLegend = drawSelfOrientedLegend( 
                fullWidth=fullWidth,
                fullHeight=fullHeight,
                tile=tile,
                ramp=self.ramp,
                mini=self.mini,
                maxi=self.maxi,
                name=self.name,
                unit=self.unitsString,
                labelValues=self.labelValues,
                roomLeftToLegend=self.coord2d[0],
                roomBelowLegend=fullHeight-self.coord2d[1],
                legendShortSide=self.size[0],
                legendLongSide=self.size[1],
                significantDigits=3,
                backgroundColor=backgroundColor,
                interpolate=self.interp,
                frame=self.visibleFrame,
                selected=(self.viewer.currentObject == self),
                numOfLabels=self.numOfLabels,
                resizeSpotRadius=self.resizeSpotRadius,
                fontScale=self.fontScale,
                glfFontID=self.glfFontID,
        )

        return 1


    def pickDraw(self):
        """called by the picking process to operate the selection
"""
        #print "colorMapLegend.pickDraw", self
        # we draw just flat quad of the insert2d
        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glPushMatrix()       
        #GL.glLoadIdentity()
        GL.glLoadMatrixf(self.viewer.currentCamera.pickMatrix) 
        GL.glOrtho(0, float(self.viewer.currentCamera.width),
                   0, float(self.viewer.currentCamera.height), -1, 1)
        GL.glMatrixMode(GL.GL_MODELVIEW)
        GL.glPushMatrix()
        GL.glLoadIdentity()
        GL.glPolygonMode(GL.GL_FRONT, GL.GL_FILL)
        #GL.glColor3f(1,0,0)

        if self.resizeSpot is not None:
            GL.glPushName(1)
            GL.glBegin(GL.GL_QUADS)
            GL.glVertex2f(float(self.resizeSpot[0]+self.resizeSpotRadius),
                          float(self.resizeSpot[1]-self.resizeSpotRadius))
            GL.glVertex2f(float(self.resizeSpot[0]+self.resizeSpotRadius),
                          float(self.resizeSpot[1]+self.resizeSpotRadius))
            GL.glVertex2f(float(self.resizeSpot[0]-self.resizeSpotRadius),
                          float(self.resizeSpot[1]+self.resizeSpotRadius))
            GL.glVertex2f(float(self.resizeSpot[0]-self.resizeSpotRadius),
                          float(self.resizeSpot[1]-self.resizeSpotRadius))
            GL.glEnd()
            GL.glPopName()

        GL.glPushName(0)
        GL.glBegin(GL.GL_QUADS)
        GL.glVertex2fv(self.polygonContour[0])
        GL.glVertex2fv(self.polygonContour[1])
        GL.glVertex2fv(self.polygonContour[2])
        GL.glVertex2fv(self.polygonContour[3])
        GL.glEnd()
        GL.glPopName()

        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glPopMatrix()
        GL.glMatrixMode(GL.GL_MODELVIEW)
        GL.glPopMatrix()


    def setPosition(self, event, redo=1):
        """the trackball transmit the translation info
"""
        #print "colorMapLegend.setPosition", event.x, event.y
        self.coord2d[0] = event.x - self.clickPosFromLegendBottomLeft[0]
        self.coord2d[1] = event.y - self.clickPosFromLegendBottomLeft[1]

        if self.coord2d[0] < 0:
            self.coord2d[0] = 0
        if self.coord2d[1] < 0:
            self.coord2d[1] = 0

        if self.coord2d[0] > self.viewer.currentCamera.width:
            self.coord2d[0] = self.viewer.currentCamera.width
        if self.coord2d[1] > self.viewer.currentCamera.height:
            self.coord2d[1] = self.viewer.currentCamera.height

        self.viewer.objectsNeedingRedo[self] = None


    def setSize(self, event, redo=1):
        """override the Insert2d function
"""
        #print "colorMapLegend.setSize", self
        if self.verticalLegend is True:
            self.size[0] = event.x - self.coord2d[0]        
            self.size[1] = self.coord2d[1] - event.y

            if self.size[0] > self.viewer.currentCamera.width:
                self.size[0] = self.viewer.currentCamera.width
            if self.size[1] > self.viewer.currentCamera.height:
                self.size[1] = self.viewer.currentCamera.height
        else:
            self.size[1] = event.x - self.coord2d[0]        
            self.size[0] = self.coord2d[1] - event.y

            if self.size[1] > self.viewer.currentCamera.width:
                self.size[1] = self.viewer.currentCamera.width
            if self.size[0] > self.viewer.currentCamera.height:
                self.size[0] = self.viewer.currentCamera.height

        if self.size[0] < 1:
            self.size[0] = 1
        if self.size[1] < 1:
            self.size[1] = 1

        if self.needsRedoDpyListOnResize and self.viewer:
            self.viewer.objectsNeedingRedo[self] = None


    def ResetPosition(self):
        self.coord2d = deepcopy(self.initialPosition)
        if self.viewer:
            self.viewer.objectsNeedingRedo[self] = None


    def respondToDoubleClick(self, event):
        """
"""
        self.showOwnGui()

        if self.needsRedoDpyListOnResize and self.viewer:
            self.viewer.objectsNeedingRedo[self] = None


    def processHit_cb(self, pick):
        #print "colorMapLegend.processHit_cb", self
        #print "pick",pick
        #print "pick.event",dir(pick)
        #print "pick.type",pick.type
        #print "pick.event",dir(pick.event)
        #print "pick.event",pick.event
        #print "pick.event.type",pick.event.type
        #print "pick.event.state",pick.event.state
        #print "pick.event.time",pick.event.time
        #print "pick.hits",pick.hits

        if ( len(pick.hits) == 1) and  pick.hits.has_key(self):
            if self.viewer.currentObject != self:
                    # if the only hit is the legend, 
                    # it becomes the current object
                    self.viewer.SetCurrentObject(self)
                    self.isMoving = True
            elif pick.event.time - self.lastPickEventTime < 200: #double click
                self.viewer.SetCurrentObject(self.viewer.rootObject)
                self.respondToDoubleClick(pick.event)
            elif pick.hits[self][0][0] == 1:
                # the click in inside the resize button
                #print "resize"
                self.isMoving = False
            elif pick.hits[self][0][0] == 0:
                # the click in inside the legend but outside 
                # the resize button
                self.isMoving = True
                self.clickPosFromLegendBottomLeft = [pick.event.x - self.coord2d[0],
                                                     pick.event.y - self.coord2d[1]]
                #print "self.clickPosFromLegendBottomLeft", self.clickPosFromLegendBottomLeft

            if self.viewer:
                self.viewer.objectsNeedingRedo[self] = None

        elif self.viewer.currentObject == self:
            #print "the insert2d is selected, but picking is outside"
            self.isMoving = None
            self.viewer.SetCurrentObject(self.viewer.rootObject)
            if self.needsRedoDpyListOnResize and self.viewer:
                self.viewer.objectsNeedingRedo[self] = None

        self.lastPickEventTime = pick.event.time


    def createOwnGui(self):
        self.ownGui = Tkinter.Toplevel()
        self.ownGui.title(self.name)
        self.ownGui.protocol('WM_DELETE_WINDOW', self.ownGui.withdraw )

        frame1 = Tkinter.Frame(self.ownGui)
        frame1.pack(side='top')

        #unit
        self.unitsEnt = Pmw.EntryField(frame1, 
                                       label_text='Units  ',
                                       labelpos='w',
                                       value=self.unitsString,
                                       command=self.setWithOwnGui)
        self.unitsEnt.pack(side='top', fill='x')

        #glf vector font
        self.glfFont = Tkinter.StringVar()
        self.glfFont.set('chicago1.glf')
        self.glfFontCB = Pmw.ComboBox(frame1, label_text='Font    ',
                                   labelpos='w',
                                   entryfield_value=self.glfFont.get(),
                                   scrolledlist_items=self.glfVectorFontList,
                                   selectioncommand=self.setWithOwnGui)
        self.glfFontCB.pack(side='top', fill='x')

        #fontScale
        self.fontScaleThumb = ThumbWheel(frame1,
                                    labCfg={'text':'font scale            ', 'side':'left'},
                                    showLabel=1, 
                                    width=90,
                                    height=14,
                                    min=0, 
                                    max=200,
                                    type=int, 
                                    value=self.fontScale,
                                    callback=self.setWithOwnGui,
                                    continuous=True,
                                    oneTurn=10,
                                    wheelPad=0)
        self.fontScaleThumb.pack(side='top')

        #label
        lLabelValuesString = ''
        for lLabelValue in self.labelValues:
            lLabelValuesString += str(lLabelValue) + ' '
        self.labelValsEnt = Pmw.EntryField(
                                frame1, 
                                label_text='Numeric labels    ',
                                labelpos='w',
                                value=lLabelValuesString,
                                command=self.setWithOwnGui)
        self.labelValsEnt.component('entry').config(width=6)
        self.labelValsEnt.pack(side='top', fill='x')

        #numOfLabel
        self.numOfLabelsCtr = ThumbWheel(frame1,
                                    labCfg={'text':'Automatic labels', 'side':'left'},
                                    showLabel=1, 
                                    width=90,
                                    height=14,
                                    min=0, 
                                    max=200,
                                    type=int, 
                                    value=self.numOfLabels,
                                    callback=self.setWithOwnGui,
                                    continuous=True,
                                    oneTurn=20,
                                    wheelPad=0)
        self.numOfLabelsCtr.pack(side='top')

        # Interpolate
        self.interpVar = Tkinter.IntVar()
        self.interpVar.set(0)
        self.checkBoxFrame = Tkinter.Checkbutton(
                                frame1, 
                                text='Interpolate',
                                variable=self.interpVar, 
                                command=self.setWithOwnGui)
        self.checkBoxFrame.pack(side='top')

        # frame
        self.frameVar = Tkinter.IntVar()
        self.frameVar.set(1)
        self.checkBoxFrame = Tkinter.Checkbutton(
                                frame1, 
                                text='Frame',
                                variable=self.frameVar, 
                                command=self.setWithOwnGui)
        self.checkBoxFrame.pack(side='top')

        # invert labels color
        self.invertLabelsColorVar = Tkinter.IntVar()
        self.invertLabelsColorVar.set(0)
        self.checkBoxinvertLabelsColor = Tkinter.Checkbutton(
                                frame1, 
                                text='Invert labels color',
                                variable=self.invertLabelsColorVar, 
                                command=self.setWithOwnGui)
        #self.checkBoxFrame.pack(side='top')
        self.checkBoxinvertLabelsColor.pack(side='top')

        # colormapguiwidget:
        self.launchColormapWidget = Tkinter.Button(
                                        frame1, 
                                        text="Show colormap settings",
                                        command=self.colormapguiRef().showColormapSettings_cb 
                                        )
        self.launchColormapWidget.pack(side='top', fill='x')


    def setWithOwnGui(self, event=None):
        #print "setWithOwnGui"

        glfFont = self.glfFontCB.get()
        fontScale = int(self.fontScaleThumb.get())
        labelValues = map(float, string.split(self.labelValsEnt.get()))
        unitsString = self.unitsEnt.get()
        numOfLabels = int(self.numOfLabelsCtr.get())

        if self.interpVar.get() == 1:
            interp = True
        else:
            interp = False

        if self.frameVar.get() == 1:
            visibleFrame = True
        else:
            visibleFrame = False

        if self.invertLabelsColorVar.get() == 1:
            invertLabelsColor = True
        else:
            invertLabelsColor = False

        self.Set(
                glfFont=glfFont,
                fontScale=fontScale,
                labelValues=labelValues, 
                numOfLabels=numOfLabels,
                unitsString=unitsString,
                interp=interp,
                visibleFrame=visibleFrame,
                invertLabelsColor=invertLabelsColor,
                updateOwnGui=False)
        self.viewer.Redraw()


    def updateOwnGui(self):
        if self.ownGui is None:
            return
        self.ownGui.title(self.name)
        self.glfFontCB.selectitem(self.glfFont.get())
        self.fontScaleThumb.set(self.fontScale)
        lLabelValuesString = ''
        for lLabelValue in self.labelValues:
            lLabelValuesString += str(lLabelValue) + ' '
        self.labelValsEnt.setentry(lLabelValuesString)
        self.unitsEnt.setentry(self.unitsString)
        self.numOfLabelsCtr.set(self.numOfLabels)
        self.interpVar.set(self.interp)
        self.invertLabelsColorVar.set(self.visibleFrame)
        self.invertLabelsColorVar.set(self.invertLabelsColor)
