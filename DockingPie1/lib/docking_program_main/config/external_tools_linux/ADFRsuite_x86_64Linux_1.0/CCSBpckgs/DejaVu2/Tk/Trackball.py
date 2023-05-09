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
# $Header: /mnt/raid/services/cvs/DejaVu2/Tk/Trackball.py,v 1.1.1.1.4.1 2017/07/13 22:22:07 annao Exp $
#
# $Id: Trackball.py,v 1.1.1.1.4.1 2017/07/13 22:22:07 annao Exp $
#

import Tkinter
import Pmw
import math
from time import time
from copy import deepcopy
from weakref import ref 

from mglutil.gui.BasicWidgets.Tk.thumbwheel import ThumbWheel
from mglutil.math.transformation import UnitQuaternion
from mglutil.gui.BasicWidgets.Tk.vector3DGUI import vectorGUI
from mglutil.util.misc import ensureFontCase

from opengltk.extent.utillib import glTrackball

import DejaVu2
from DejaVu2.Transformable import Transformable
from DejaVu2.viewerConst import FPRECISION, YES, NO
from DejaVu2.viewerFns import getkw
from DejaVu2.cursors import cursorsDict

# remove this line when we go to version 1.46
if hasattr( DejaVu2, 'defaultSpinningMode') is False: DejaVu2.defaultSpinningMode = 0


class Trackball:

    def __init__(self, camera, size=0.8, rscale=2.0, tscale=0.05,
                 sscale=0.01, renorm=97):
        #print "Trackball.__init___", DejaVu2.defaultSpinningMode
        #import traceback;traceback.print_stack()

        for mod in camera.mouseButtonModifiers:
            if mod=='None': mod=''
            setattr(self, mod+'B1motion', self.NoFunc)
            setattr(self, mod+'B2motion', self.NoFunc)
            setattr(self, mod+'B3motion', self.NoFunc)

        self.reportFrameRate = 0
        
        # create a trackball
        self.tb = glTrackball(size, rscale, renorm)
        self.vectorXY = [0.0 ,0.0]
        self.vectorZ = 0.0
        self.transScale = tscale
        self.scaleFactor = 0.0
        self.scaleScale = sscale

        # save the Opengl widget for this trackball
        self.camera = ref(camera)

        # Current coordinates of the mouse.
        self.xmouse = 0
        self.ymouse = 0
        self.xmouseDown = 0
        self.ymouseDown = 0
        self.xmouseUp = 0
        self.ymouseUp = 0

        # setup the spin
        self.cancelId = 0 # used to stop the spinning
        self.spinRotationMat = UnitQuaternion( data=(0.,-1.,0.,1.) ).getRotMatrix(shape=(16,))

        self.spinVar = self.camera().viewer.spinVar

        self.lastSpinVar = 0
        self.resetSpinValues(updateSpinGui=False)
        self.createSpinGui()
        self.spinGui.withdraw()

        # Basic bindings for the virtual trackball
        self.ResetBindings()


    def set(self, otherTrackball):
        self.spinRotationMat = otherTrackball.spinRotationMat
        self.spinVar.set( otherTrackball.spinVar.get() )
        self.lastSpinVar = otherTrackball.lastSpinVar
        self.spinAxis = otherTrackball.spinAxis
        self.spinAnglePerFrame = otherTrackball.spinAnglePerFrame
        self.rockAngle = otherTrackball.rockAngle
        self.numOfFrameRock = otherTrackball.numOfFrameRock
        self.rockCount = otherTrackball.rockCount
        self.updateSpinGui()
        if otherTrackball.cancelId != 0:
            self.spinCycle()


    def resetSpinValues(self, updateSpinGui=True):
        self.spinAxis = [0., 1., 0.]
        self.spinAnglePerFrame = .5
        self.rockAngle = 30
        self.numOfFrameRock = abs(self.rockAngle / self.spinAnglePerFrame)
        self.rockCount = self.numOfFrameRock / 2
        if updateSpinGui is True:
            self.updateSpinGui()
        if self.cancelId != 0:
            self.camera().restoreAA(None)
            self.camera().after_cancel(self.cancelId)
            self.cancelId = 0


    def ResetBindings(self):
	"""binds default callbacks to events"""

        cam = self.camera()
        evm = cam.eventManager

        # Any-Button-1 acts the same as Button-1 when bound to widget;
        # therefore, Modifier-ButtonPress-1 overwrites it; must specify for
        # modifier
        for mod in cam.mouseButtonModifiers:
            if mod=='None': mod=''
            else: mod=mod+'-'

            evm.AddCallback('<'+mod+'ButtonPress-1>', self.RecordMouse)
            evm.AddCallback('<'+mod+'ButtonPress-2>', self.RecordMouse)
            evm.AddCallback('<'+mod+'ButtonPress-3>', self.RecordMouse)
	                        
            evm.AddCallback('<'+mod+'ButtonPress-1>', self.MouseButtonDown)
            evm.AddCallback('<'+mod+'ButtonPress-2>', self.MouseButtonDown)
            evm.AddCallback('<'+mod+'ButtonPress-3>', self.MouseButtonDown)

            evm.AddCallback('<'+mod+'ButtonRelease-1>', self.MouseButtonUp)
            evm.AddCallback('<'+mod+'ButtonRelease-2>', self.MouseButtonUp)
            evm.AddCallback('<'+mod+'ButtonRelease-3>', self.MouseButtonUp)

            evm.AddCallback('<'+mod+'ButtonRelease-1>', self.setDefaultCursor)
            evm.AddCallback('<'+mod+'ButtonRelease-2>', self.setDefaultCursor)
            evm.AddCallback('<'+mod+'ButtonRelease-3>', self.setDefaultCursor)

        for mod in cam.mouseButtonModifiers:
            if mod=='None': mod=mod1=''
            else: mod1=mod+'-'
            evm.AddCallback('<'+mod1+'B1-Motion>',
                            eval('self.'+mod+'B1motion_cb'))
            evm.AddCallback('<'+mod1+'B2-Motion>',
                            eval('self.'+mod+'B2motion_cb'))
            evm.AddCallback('<'+mod1+'B3-Motion>',
                            eval('self.'+mod+'B3motion_cb'))

        # we put this after mouse up so we leave more time for the mouse
        # to move and launch the spin
        lButton, lModifier= cam.findButton('rotation', 'Object')
        self.rotationButton = lButton
        self.rotationModifier = lModifier
        if lModifier == 'None':
            lModifier = ''
        else:
            lModifier += '-'
        evm.AddCallback('<'+lModifier+'ButtonRelease-'+str(lButton)+'>', self.launchTrackbalSpin_cb)


    def setDefaultCursor(self, event):
        self.camera().configure(cursor=cursorsDict['default'])

        
    def trackMotion_cb(self, event):
        self.xrootBmotionPrevious = self.xrootBmotion
        self.yrootBmotionPrevious = self.yrootBmotion
        self.xrootBmotion = event.x_root
        self.yrootBmotion = event.y_root


    def B1motion_cb(self, event):
        """Default callback for rotation, computes rotation matrix"""
        self.computeTransformation(event)
        self.B1motion(event, self.tb.mat, self.vectorXY, self.vectorZ)
        self.RecordMouse(event)


    def B2motion_cb(self, event):
        """Default callback for rotation, computes rotation matrix"""
        #print "B2motion_cb"
        self.computeTransformation(event)
        self.B2motion(event, self.tb.mat, self.vectorXY, self.vectorZ)
        self.RecordMouse(event)


    def B3motion_cb(self, event):
        """Default callback for rotation, computes rotation matrix"""
        self.computeTransformation(event)
        self.B3motion(event, self.tb.mat, self.vectorXY, self.vectorZ)
        self.RecordMouse(event)


    def ShiftB1motion_cb(self, event):
        """Default callback for rotation, computes rotation matrix"""
        self.computeTransformation(event)
        self.ShiftB1motion(event, self.tb.mat, self.vectorXY, self.vectorZ)
        self.RecordMouse(event)


    def ShiftB2motion_cb(self, event):
	"""Default callback for rotation, computes rotation matrix"""
        self.computeTransformation(event)
	self.ShiftB2motion(event, self.tb.mat, self.vectorXY, self.vectorZ)
	self.RecordMouse(event)


    def ShiftB3motion_cb(self, event):
	"""Default callback for rotation, computes rotation matrix"""
        self.computeTransformation(event)
	self.ShiftB3motion(event, self.tb.mat, self.vectorXY, self.vectorZ)
	self.RecordMouse(event)


    def ControlB1motion_cb(self, event):
	"""Default callback for rotation, computes rotation matrix"""
        self.computeTransformation(event)
	self.ControlB1motion(event, self.tb.mat, self.vectorXY, self.vectorZ)
	self.RecordMouse(event)


    def ControlB2motion_cb(self, event):
	"""Default callback for rotation, computes rotation matrix"""
        self.computeTransformation(event)
	self.ControlB2motion(event, self.tb.mat, self.vectorXY, self.vectorZ)
	self.RecordMouse(event)


    def ControlB3motion_cb(self, event):
	"""Default callback for rotation, computes rotation matrix"""
        self.computeTransformation(event)
	self.ControlB3motion(event, self.tb.mat, self.vectorXY, self.vectorZ)
	self.RecordMouse(event)


    def AltB1motion_cb(self, event):
	"""Default callback for rotation, computes rotation matrix"""
        self.computeTransformation(event)
	self.AltB1motion(event, self.tb.mat, self.vectorXY, self.vectorZ)
	self.RecordMouse(event)


    def AltB2motion_cb(self, event):
	"""Default callback for rotation, computes rotation matrix"""
        self.computeTransformation(event)
	self.AltB2motion(event, self.tb.mat, self.vectorXY, self.vectorZ)
	self.RecordMouse(event)


    def AltB3motion_cb(self, event):
	"""Default callback for rotation, computes rotation matrix"""
        self.computeTransformation(event)
        self.AltB3motion(event, self.tb.mat, self.vectorXY, self.vectorZ)
	self.RecordMouse(event)


    def MetaB1motion_cb(self, event):
	"""Default callback for rotation, computes rotation matrix"""
        self.computeTransformation(event)
	self.MetaB1motion(event, self.tb.mat, self.vectorXY, self.vectorZ)
	self.RecordMouse(event)


    def MetaB2motion_cb(self, event):
	"""Default callback for rotation, computes rotation matrix"""
        self.computeTransformation(event)
	self.MetaB2motion(event, self.tb.mat, self.vectorXY, self.vectorZ)
	self.RecordMouse(event)


    def MetaB3motion_cb(self, event):
	"""Default callback for rotation, computes rotation matrix"""
        self.computeTransformation(event)
        self.MetaB3motion(event, self.tb.mat, self.vectorXY, self.vectorZ)
	self.RecordMouse(event)


    def computeTransformation(self, event):
        """compute the rotation corresponding to the latest mouse motion
and the (X,Y) translation corresponding to the latest mouse motion
and the Z translation corresponding to the latest mouse motion
"""
        #print "computeTransformation", event
        #print "self.xmouse, x", self.xmouse, event.x
        #print "self.ymouse, y", self.ymouse, event.y

        self.trackMotion_cb(event)

	self.tb.update(self.xmouse, self.ymouse, event.x, event.y,
			    self.camera()._width, self.camera()._height, 1 )
        # now self.tb.mat has the rotation matrix
	self.vectorXY[0] = dx = event.x - self.xmouse
	self.vectorXY[1] = dy = self.ymouse - event.y
        
	if abs(dx) > abs(dy): delta = dx
	else: delta = dy
	self.vectorZ = delta
        # new self.vector[:2] has the (X,Y) translation vector
        # and self.vector[2] has the Z translation


    def NoFunc(self, event, matrix, transXY, transZ):
	pass


    def cursorPosition(self):
	"""get new mouse position
"""
	c = self.camera()
	x = c.tk.getint(c.tk.call('winfo', 'pointerx', c._w))
	y = c.tk.getint(c.tk.call('winfo', 'pointery', c._w))
	return x,y


    def MouseButtonDown(self, event):
	"""Do all we have to do every time a mouse button goes down
"""
        #print "MouseButtonDown", event.num
        #print "Focus on MouseButtonDown", self.camera().focus_get()
        #if event.num in [2]:

        # set cursor
        cam = self.camera()
        cam.setCursor(event.num)
        v = cam.viewer

        # find picked pixel
        point, bg = v.get3DPointFromPick(event)
        v.lastPickedPixel = point
        
        self.xrootBmotionPrevious = event.x_root
        self.yrootBmotionPrevious = event.y_root
        self.xrootBmotion = event.x_root
        self.yrootBmotion = event.y_root

        if ((event.num in [1, self.rotationButton]) and (self.cancelId != 0)):
            cam.restoreAA(None)
            cam.after_cancel(self.cancelId)
            self.cancelId = 0
            self.resetSpinValues()

        for c in v.cameras:
            for func in c.onButtonDownCBlist:
                func(event)
                
        v.lastPickedCamera = v.currentCamera
        self.xmouseDown = event.x
        self.ymouseDown = event.y
        v.timeButtonDown = time()
        v.numberOfRedraws = 0
        #print 'Down', self.xmouseDown, self.ymouseDown


    def MouseButtonUp(self, event):
	"""Do all we have to do every time a mouse button is released
"""
        self.camera().configure(cursor=cursorsDict['default'])
        #print "Trackball.MouseButtonUp", event.num, event, time()
        v = self.camera().viewer
        v.lastPickedPixel = None
        if self.cancelId == 0: # if we are spinning, we do as if the button wasn't up
            v.timeButtonUp = time()
            if self.reportFrameRate:
                print v.numberOfRedraws / (v.timeButtonUp - v.timeButtonDown)
            for c in v.cameras:
                for func in c.onButtonUpCBlist:
                    #print "MouseButtonUp func", func
                    func(event)

            v.lastPickedCamera = v.currentCamera
            self.xmouseUp = event.x
            self.ymouseUp = event.y
            #print 'UP', self.xmouseUp, self.ymouseUp
            v.Redraw()
        #print "Focus after MouseButtonUP", self.camera().focus_get()


    def RecordMouse(self, event):
	"""Record the current mouse position
"""
        #print "RecordMouse"
        #v = self.camera().viewer
        #v.trackball = self
        self.xmouse = event.x
        self.ymouse = event.y


    def launchTrackbalSpin_cb(self, event):
        """If motion after button release start auto spin
"""
        #print "Trackball.launchTrackbalSpin_cb", event, time()
        c = self.camera()
        vi = c.viewer
        if vi.Xform.get() == 'Object':
            lSpinVar = self.spinVar.get()
            if lSpinVar in (1, 2, 3):
                x,y = self.cursorPosition()
                lDistanceX = (self.xrootBmotionPrevious - x)
                lDistanceX2 = lDistanceX * lDistanceX
                lDistanceY = (self.yrootBmotionPrevious - y)
                lDistanceY2 = lDistanceY * lDistanceY
                lSquaredDistance = lDistanceX2 + lDistanceY2
                if lSquaredDistance > 81:
                    #it seems that the quaternions in dejavu.i
                    #are not the same as those in transformation.py
                    lUnitQuat = UnitQuaternion( data=(self.tb.quat[3], self.tb.quat[:-1] ) )
                    self.spinAxis, lAngle = lUnitQuat.getAxisAndAngleDegres()
                    self.spinAnglePerFrame = lAngle / 4
                    self.spinGui.anglePerFrameThumbwheel.set(self.spinAnglePerFrame, update=0)
                    self.spinGui.vectorGUI.vector = list(self.spinAxis)
                    self.spinGui.vectorGUI.setEntries()
                    self.spinGui.rockAngleThumbwheel.set(120, update=0)
                    self.setWithSpinGui()
                    self.lastSpinVar = self.spinVar.get()
                    self.rockCount = self.numOfFrameRock / 2
                    self.spinCycle()
                elif self.cancelId != 0:
                    c.restoreAA(None)
                    c.after_cancel(self.cancelId)
                    self.cancelId = 0
                    self.resetSpinValues()


    def makeSpinRotationMat(self, spinAxis, spinAngle):
        """
"""
        #print "Trackball.makeSpinRot"
        #it seems that the quaternions in dejavu.i
        #are not the same as those in transformation.py
        lUnitQuat = UnitQuaternion( data=(spinAxis[0], spinAxis[1],
                                          spinAxis[2], spinAngle ) )
        return lUnitQuat.getRotMatrix(shape=(16,))


    def toggleCycle(self, docb=True):
        #print "toggleCycle", self.cancelId
        lCurrentSpinVar = self.spinVar.get()
        if   (self.cancelId == 0) \
          or (self.lastSpinVar != lCurrentSpinVar):
            if lCurrentSpinVar in [0, 1]:
                self.rockCount = self.numOfFrameRock / 2
            self.spinRotationMat = self.makeSpinRotationMat(
                                        self.spinAxis, 
                                        self.spinAnglePerFrame
                                        )
            self.spinCycle()
        else:
            self.camera().restoreAA(None)
            self.camera().after_cancel(self.cancelId)
            self.cancelId = 0

        if docb and self.camera().viewer.spinCallBack is not None:
            self.camera().viewer.spinCallBack(lCurrentSpinVar)


    def spinCycle(self):
        """Spin the object in continuous mode
"""
        #print "Trackball.doSpin", self.tb.mat
        c = self.camera()
        c.after_cancel(self.cancelId)
        vi = c.viewer
        if vi.Xform.get() == 'Object':
            c.suspendAA(None)
            self.lastSpinVar = self.spinVar.get()
            if self.lastSpinVar == 1:
                vi.rootObject.ConcatRotation(self.spinRotationMat)
                c.Redraw()
                vi.afterRedraw()
                self.cancelId = c.after(10, self.spinCycle)
            elif self.lastSpinVar in [2, 3]:
                self.cancelId = c.after(10, self.bounceCycle)
            else:
                c.restoreAA(None)      
        else:
            c.restoreAA(None)      


    def bounceCycle(self):
        """Bounce the object in continuous mode
"""
        #print "Trackball.doRock", self.tb.mat
        c = self.camera()
        c.after_cancel(self.cancelId)    
        vi = c.viewer
        if vi.Xform.get() == 'Object':
            c.suspendAA(None)
            self.lastSpinVar = self.spinVar.get()
            if self.lastSpinVar == 2:
                if self.rockCount >= self.numOfFrameRock:
                    self.rockCount = 0
                    self.spinAnglePerFrame = - self.spinAnglePerFrame
                    self.spinRotationMat = self.makeSpinRotationMat(
                                                self.spinAxis, 
                                                self.spinAnglePerFrame
                                                )
                else:
                    self.rockCount += 1

                vi.rootObject.ConcatRotation(self.spinRotationMat)
                c.Redraw()
                vi.afterRedraw()
                self.cancelId = c.after(10, self.bounceCycle)
            elif self.lastSpinVar == 3:
                self.cancelId = c.after(10, self.oscillateCycle)
            elif self.lastSpinVar == 1:
                self.rockCount = self.numOfFrameRock / 2                    
                self.cancelId = c.after(10, self.spinCycle)
            else:
                c.restoreAA(None)      
        else:
            c.restoreAA(None)      


    def oscillateCycle(self):
        """Oscillate the object in continuous mode
"""
        #print "Trackball.bounceCycle", self.tb.mat
        c = self.camera()
        c.after_cancel(self.cancelId)
        vi = c.viewer
        vi.master.update_idletasks()
        if vi.Xform.get() == 'Object':
            c.suspendAA(None)
            self.lastSpinVar = self.spinVar.get()
            if self.lastSpinVar == 3:
                lNumOfFrameRock = self.numOfFrameRock
                if self.rockCount >= lNumOfFrameRock:
                    self.rockCount = 0
                    self.spinAnglePerFrame = - self.spinAnglePerFrame
                else:
                    self.rockCount += 1

                lAtenuationFactor = math.sin(3.14159*self.rockCount/lNumOfFrameRock) * 1.57
                self.spinRotationMat = self.makeSpinRotationMat(
                                                self.spinAxis, 
                                                self.spinAnglePerFrame*lAtenuationFactor
                                                )
                vi.rootObject.ConcatRotation(self.spinRotationMat)
                c.Redraw()
                vi.afterRedraw()
                self.cancelId = c.after(10, self.oscillateCycle)
            elif self.lastSpinVar == 2:
                self.spinRotationMat = self.makeSpinRotationMat(
                                                self.spinAxis, 
                                                self.spinAnglePerFrame
                                                )
                self.cancelId = c.after(10, self.bounceCycle)
            elif self.lastSpinVar == 1:
                self.rockCount = self.numOfFrameRock / 2                    
                self.spinRotationMat = self.makeSpinRotationMat(
                                                self.spinAxis, 
                                                self.spinAnglePerFrame
                                                )
                self.cancelId = c.after(10, self.spinCycle)
            else:
                c.restoreAA(None)      
        else:
            c.restoreAA(None)      


    def showSpinGui(self, event=None):
        #print "showSpinGui", self
        if self.spinGui.winfo_ismapped() == 0:
            self.spinGui.deiconify()
        self.spinGui.lift()


    def hideSpinGui(self, event=None):
        #print "hideSpinGui", self
        if self.spinGui.winfo_ismapped() == 1:
            self.spinGui.withdraw()


    def toggleSpinGui(self, event=None):
        #print "toggleSpinGui", self
        if self.spinGui.winfo_ismapped() == 1:
            self.spinGui.withdraw()
        else:
            self.spinGui.deiconify()
            self.spinGui.lift()


    def createSpinGui(self):
        #print "createSpinGui"
        
        self.spinGui = Tkinter.Toplevel()
        self.spinGui.title('Spin Settings')
        self.spinGui.protocol('WM_DELETE_WINDOW', self.spinGui.withdraw )

        mainFrame = Tkinter.Frame(self.spinGui, relief='ridge', borderwidth=3)

        lLabelSpinBold = Tkinter.Label(mainFrame,
                      text="Spinning only occurs in \"Object\" transform mode",
                      font=(ensureFontCase('helvetica'),10,'bold'))
        lLabelSpinBold.pack(side='top')

        lLabelSpin = Tkinter.Label(mainFrame,
                      font=(ensureFontCase('helvetica'),'9'),
                      justify='left',
                      text="""start: release the middle button while the mouse is still moving
stop: left click
"""
                     )
        lLabelSpin.pack(side='top')

        anglesFrame = Tkinter.Frame(mainFrame, relief='ridge', borderwidth=3)
        # off/spin/rock radio buttons
        radioFrame = Tkinter.Frame(anglesFrame)#, relief='ridge', borderwidth=3)

        radioOff = Tkinter.Radiobutton(
            radioFrame,
            text='Off',
            value=0,
            variable=self.spinVar,
            width=8,
            indicatoron=0,
            command=self.resetSpinValues,
            )
        radioOff.grid(row=0, column=0, sticky='we')
        radioSpin = Tkinter.Radiobutton(
            radioFrame,
            text='Spin',
            value=1,
            variable=self.spinVar,
            width=8,
            indicatoron=0,
            command=self.toggleCycle,
            )
        radioSpin.grid(row=0, column=1, sticky='we')
        radioBounce = Tkinter.Radiobutton(
            radioFrame,
            text='Bounce',
            value=2,
            variable=self.spinVar,
            width=8,
            indicatoron=0,
            command=self.toggleCycle,
            )
        radioBounce.grid(row=0, column=2, sticky='we')
        radioOscillate = Tkinter.Radiobutton(
            radioFrame,
            text='Oscillate',
            value=3,
            variable=self.spinVar,
            width=8,
            indicatoron=0,
            command=self.toggleCycle,
            )
        radioOscillate.grid(row=0, column=3, sticky='we')
        radioFrame.pack(side='top')

        self.spinGui.anglePerFrameThumbwheel = ThumbWheel(
                                    anglesFrame,
                                    labCfg={'text':'Angular speed (degrees per frame)', 'side':'left'},
                                    showLabel=1, 
                                    width=90,
                                    height=14,
                                    min=.01, 
                                    max=10.,
                                    type=float, 
                                    value=self.spinAnglePerFrame,
                                    callback=self.setWithSpinGui,
                                    continuous=True,
                                    oneTurn=5,
                                    wheelPad=0,
                                    )
        self.spinGui.anglePerFrameThumbwheel.pack(side='top', anchor='e')
        self.spinGui.rockAngleThumbwheel = ThumbWheel(
                                    anglesFrame,
                                    labCfg={'text':'Rock angle (degrees)', 'side':'left'},
                                    showLabel=1, 
                                    width=90,
                                    height=14,
                                    min=1, 
                                    max=360,
                                    type=int, 
                                    value=self.rockAngle,
                                    callback=self.setWithSpinGui,
                                    continuous=True,
                                    oneTurn=60,
                                    wheelPad=0,
                                    )
        self.spinGui.rockAngleThumbwheel.pack(side='top', anchor='e')
        def shiftRockCenter(val):
            lSpinRotMat = self.makeSpinRotationMat(self.spinAxis, -val)
            c = self.camera()
            c.viewer.rootObject.ConcatRotation(lSpinRotMat)
            c.Redraw()
            vi.afterRedraw()
            shiftRockCenterThumbwheel.set(0, update=0)
        shiftRockCenterThumbwheel = ThumbWheel(
                                    anglesFrame,
                                    labCfg={'text':'Shift rock center', 'side':'left'},
                                    showLabel=0, 
                                    width=90,
                                    height=14,
                                    min=-100, 
                                    max=100,
                                    type=int, 
                                    value=0,
                                    callback=shiftRockCenter,
                                    continuous=True,
                                    oneTurn=90,
                                    wheelPad=0,
                                    )
        shiftRockCenterThumbwheel.pack(side='top', anchor='e')
        anglesFrame.pack(side='top')

        self.spinGui.vectorGUI = vectorGUI(mainFrame,
                                           size=150,
                                           vector=self.spinAxis,
                                           callback=self.setWithSpinGui,
                                           )

        mainFrame.pack(side='top')


    def updateSpinGui(self):
        self.spinGui.rockAngleThumbwheel.set(self.rockAngle, update=0)
        self.spinGui.anglePerFrameThumbwheel.set(self.spinAnglePerFrame, update=0)
        self.spinGui.vectorGUI.vector = self.spinAxis
        self.spinGui.vectorGUI.setEntries()

        if self.camera().viewer.spinCallBack is not None:
            self.camera().viewer.spinCallBack(self.spinVar.get())


    def setWithSpinGui(self, event=None):
        """
"""
        #print "setWithSpinGui", event
        self.rockAngle = self.spinGui.rockAngleThumbwheel.get()
        if (self.spinVar.get() in [2, 3] ) and (self.spinAnglePerFrame < 0):
            self.spinAnglePerFrame = - self.spinGui.anglePerFrameThumbwheel.get()
        else:
            self.spinAnglePerFrame = self.spinGui.anglePerFrameThumbwheel.get()
        self.numOfFrameRock = abs(self.rockAngle / self.spinAnglePerFrame)
        self.spinAxis = self.spinGui.vectorGUI.vector
        self.spinRotationMat = self.makeSpinRotationMat(self.spinAxis, self.spinAnglePerFrame)
