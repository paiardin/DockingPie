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
# $Header: /mnt/raid/services/cvs/DejaVu2/Qt/Camera.py,v 1.4.4.1 2017/07/13 22:25:45 annao Exp $
#
# $Id: Camera.py,v 1.4.4.1 2017/07/13 22:25:45 annao Exp $
#
import numpy

from PySide import QtGui, QtCore, QtOpenGL
from math import pow, exp, sqrt

import DejaVu2
#from DejaVu2.Qt.EventHandler import EventManager
from DejaVu2.Camera import StandardCameraBase
from DejaVu2.Transformable import Transformable
from DejaVu2.Insert2d import Insert2d

from mglutil.gui import widgetsOnBackWindowsCanGrabFocus

from opengltk.OpenGL import GL
from opengltk.OpenGL.GL import *
from opengltk.extent.utillib import glTrackball
from opengltk.OpenGL.GLU import gluUnProject

class Camera(QtOpenGL.QGLWidget, StandardCameraBase):
    keyPressSignal = QtCore.Signal(str)
    def __init__(self, master, screenName, viewer, num, check=1,
                 cnf={}, **kw):
        
        #print 'CAMERA MASTER', id(master)
        self.name = 'Camera'+str(num)
        self.num = num
        self.uniqID = self.name+viewer.uniqID

        fmt = QtOpenGL.QGLFormat()
        fmt.setDoubleBuffer(kw.get('double', True))
        fmt.setOverlay(kw.get('overlay', True))
        fmt.setDepth(kw.get('depth', True))
        fmt.setStencil(kw.get('stencil', True))
        fmt.setAccum(kw.get('accum', True))
        stereo = kw.get('stereo', False)
        if stereo=='native': stereo = True
        else: stereo = False
        fmt.setStereo(stereo)
        fmt.setSampleBuffers(kw.get('sampleBuffers', True))
        #print 'MASTER', master
        sharedContext = kw.get('sharecontext', None)
        if sharedContext is None and len(viewer.cameras):
             # share context with the default camera.
             c0 = viewer.cameras[0]
             sharecontextWidget = viewer.cameras[0]
             QtOpenGL.QGLWidget.__init__(self, fmt, #c0.context(),
                                         parent=master,
                                         shareWidget=c0)
             #import pdb
             #pdb.set_trace()
        else:
             QtOpenGL.QGLWidget.__init__(self, fmt, parent=master)
        self.setAutoBufferSwap(0)
        print "CONTEXTS", self.name, self.getContext()
        #self.gridLayout = QtGui.QGridLayout(master)
        #self.gridLayout.addWidget(self, 0, 0, 2, 1)
        #self.setLayout(self.gridLayout)
        
        # make sure we get mouse move event even no button is pressed
        self.setMouseTracking(True)
        #self._motionAxis = None
        
        self._width = self.width()
        self._height = self.height()

        from opengltk.exception import GLerror
        try:
            self.makeCurrent()
            GL.glAccum(GL.GL_LOAD, 1.0)
            viewer.accumBuffersError = False
        except GLerror:
            viewer.accumBuffersError = True

        if num==0:
            viewer.redrawTimer = QtCore.QTimer()
            self.connect(viewer.redrawTimer, QtCore.SIGNAL("timeout()"),
                         viewer.ReallyRedraw)
            
        StandardCameraBase.__init__(self, master, screenName, viewer, num, **kw)

##         if hasattr(viewer, 'accumBuffersError') and viewer.accumBuffersError:
##             # if we tried to create a context with accumulation buffers before
##             # and it failed we set accum to False
##             kw['accum'] = 0
##         else:
##             # if the user did not specify accumulation we turn them on by defaut
##             if not kw.has_key('accum'):
##                 kw['accum'] = 1

##         if not kw.has_key('ident'):
##             kw['ident'] = self.uniqID
## #            kw['ident'] = 'camera%d' % num

##         if not kw.has_key('sharelist') and len(viewer.cameras):
##             # share list with the default camera.
##             cam = viewer.cameras[0]
##             kw['sharelist'] = cam.uniqID

##         if not kw.has_key('sharecontext') and len(viewer.cameras):
##             # share context with the default camera.
##             cam = viewer.cameras[0]
##             kw['sharecontext'] = cam.uniqID
##             if not hasattr(cam, 'shareCTXWith'): cam.shareCTXWith = []
##             cam.shareCTXWith.append(self)

##         #self.frame.master.protocol("WM_DELETE_WINDOW", self.hide)

##         cfg = 0

##         if 'width' in cnf.keys():
##             self.width = cnf['width']
##             cfg = 1
##         else:
##             cnf['width'] = self.width = 406

##         if 'height' in cnf.keys():
##             self.height = cnf['height']
##             cfg = 1
##         else:
##             cnf['height'] = self.height = 406

##         self.rootx = 320
##         if 'rootx' in cnf.keys():
##             self.rootx = cnf['rootx']
##             del cnf['rootx']
##             cfg = 1
##         self.rooty = 180
##         if 'rooty' in cnf.keys():
##             self.rooty = cnf['rooty']
##             del cnf['rooty']
##             cfg = 1
##         if cfg: self.Geometry()

##         if 'side' in cnf.keys():
##             side = cnf['side']
##             del cnf['side']
##         else: side = 'top'

##         # after this line self.master will be set to self frame
##         from Tkinter import TclError
##         from opengltk.exception import GLerror
##         try:
##             Tkinter.Widget.__init__(self, self.frame, 'togl', cnf, kw)
##             try:
##                 GL.glAccum(GL.GL_LOAD, 1.0)
##                 viewer.accumBuffersError = False
##             except GLerror:
##                 viewer.accumBuffersError = True
##                 #viewer.accumBuffersError = False
##         except TclError, e:
##             print "Warning: disabling accumulation buffers",e
##             viewer.accumBuffersError = True
##             kw.pop('accum')
##             Tkinter.Widget.__init__(self, self.frame, 'togl', cnf, kw)

        self.preventIntelBug_BlackTriangles()
        
        # x,y coordinates used to compute delta during mouse move
        self._oldx = self._oldy = 0
        # x,y coordinate when button is first pressed
        self._orginalx = self._orginaly = 0
        self._isPressed = False
        self._hasMoved = False # will be set to true if cursor moves after
                              # button press

        # dictionaries of functions to be called upon mouse button events
        # we create a dictionary for each button event
        d1 = self._mousePressActions = {}
        d2 = self._mouseMoveActions = {}
        d3 = self._mouseReleaseNoMotionActions = {}
        d4 = self._mouseReleaseWithMotionActions = {}
        
        # - the first key in these dicts is the button number QtCore.Qt.[
        #   LeftButton, MidButton, RightButton, and combinations].
        # - the next key is the modifiers

        for d in [d1,d2,d3,d4]:
            for button in [
                QtCore.Qt.LeftButton,
                QtCore.Qt.MidButton,
                QtCore.Qt.RightButton,
                QtCore.Qt.LeftButton | QtCore.Qt.MidButton,
                QtCore.Qt.LeftButton |  QtCore.Qt.RightButton,
                QtCore.Qt.MidButton | QtCore.Qt.RightButton,
                QtCore.Qt.LeftButton | QtCore.Qt.MidButton | QtCore.Qt.RightButton ]:
                d1 = {}
                d[int(button)] = d1
                for mod in [
                    QtCore.Qt.NoModifier, 
                    QtCore.Qt.ShiftModifier,
                    QtCore.Qt.ControlModifier,
                    QtCore.Qt.AltModifier,
                    QtCore.Qt.MetaModifier,
                    QtCore.Qt.KeypadModifier,
                    QtCore.Qt.ShiftModifier | QtCore.Qt.ControlModifier,
                    QtCore.Qt.ShiftModifier | QtCore.Qt.AltModifier,
                    QtCore.Qt.ShiftModifier | QtCore.Qt.MetaModifier,
                    QtCore.Qt.ShiftModifier | QtCore.Qt.KeypadModifier,
                    QtCore.Qt.ControlModifier | QtCore.Qt.AltModifier,
                    QtCore.Qt.ControlModifier | QtCore.Qt.MetaModifier,
                    QtCore.Qt.ControlModifier | QtCore.Qt.KeypadModifier,
                    QtCore.Qt.AltModifier | QtCore.Qt.MetaModifier,
                    QtCore.Qt.AltModifier | QtCore.Qt.KeypadModifier,
                    ]:
                    d1[int(mod)] = None

        self._mouseWheelActions = {}
        self._mouseMoveNoButtonActions = {}
        for mod in [
            QtCore.Qt.NoModifier, 
            QtCore.Qt.ShiftModifier,
            QtCore.Qt.ControlModifier,
            QtCore.Qt.AltModifier,
            QtCore.Qt.MetaModifier,
            QtCore.Qt.KeypadModifier,
            QtCore.Qt.ShiftModifier | QtCore.Qt.ControlModifier,
            QtCore.Qt.ShiftModifier | QtCore.Qt.AltModifier,
            QtCore.Qt.ShiftModifier | QtCore.Qt.MetaModifier,
            QtCore.Qt.ShiftModifier | QtCore.Qt.KeypadModifier,
            QtCore.Qt.ControlModifier | QtCore.Qt.AltModifier,
            QtCore.Qt.ControlModifier | QtCore.Qt.MetaModifier,
            QtCore.Qt.ControlModifier | QtCore.Qt.KeypadModifier,
            QtCore.Qt.AltModifier | QtCore.Qt.MetaModifier,
            QtCore.Qt.AltModifier | QtCore.Qt.KeypadModifier,
            ]:
            self._mouseWheelActions[int(mod)] = None
            self._mouseMoveNoButtonActions[int(mod)] = None
            
        # create initial bindings
        self._mouseReleaseNoMotionActions[int(QtCore.Qt.LeftButton)][int(QtCore.Qt.NoModifier)] = self.pick
        self._mouseMoveActions[int(QtCore.Qt.LeftButton)][int(QtCore.Qt.NoModifier)] = self.rotateScene
        self._mouseMoveActions[int(QtCore.Qt.MidButton)][int(QtCore.Qt.NoModifier)] = self.translateSceneXY
        self._mouseMoveActions[int(QtCore.Qt.RightButton)][int(QtCore.Qt.NoModifier)] = self.handleZoomCamera

        self._mouseMoveActions[int(QtCore.Qt.LeftButton)][int(QtCore.Qt.ShiftModifier)] = self.addToSelection
        self._mouseReleaseWithMotionActions[int(QtCore.Qt.LeftButton)][int(QtCore.Qt.ShiftModifier)] = self.endSelectionRectangle

        self._mouseWheelActions[int(QtCore.Qt.NoModifier)] = self.handleScaleSlab
        self._mouseWheelActions[int(QtCore.Qt.ShiftModifier)] = self.translateSceneZ
        self._mouseWheelActions[int(QtCore.Qt.ControlModifier)] = self.handleTranslateSlab

        self._mouseMoveNoButtonActions[int(QtCore.Qt.NoModifier)] = self.showTools
        #self._mouseMoveActions[int(QtCore.Qt.RightButton)][int(QtCore.Qt.NoModifier)] = self.translateSceneZ
        #self._mouseReleaseWithMotionActions[QtCore.Qt.NoModifier][QtCore.Qt.LeftButton] = self.
        #print self._mouseReleaseNoMotionActions

        # create a trackball rotation calculator
        self.tb = glTrackball(0.8, 2.0, 97)


    ##
    ## programmatic methods for changing the camera
    ##
    def zoomCamera(self, delta):
        """Zoom camera by a factor of delta that can be positive or negative
        by changing the FOV
        """
        #import pdb; pdb.set_trace()
        lFovyFactor = self.fovy / self.fovyNeutral
        if lFovyFactor > 1:
            lFovyFactor = pow(lFovyFactor, 1.5)
        lDistanceFactor = (self.far - self.near)  * .001
        lDistanceFactor = exp(lDistanceFactor) - 1
        lFactor = lDistanceFactor * lFovyFactor * .5

        if delta>0:
            delta = -37.0*lFactor
        elif delta < 0:
            delta = 37.0*lFactor
        else:
            return

        x,y,z = self.direction
        d = self.direction*delta/sqrt(x*x+y*x+z*z)
        self._setLookFrom(self.lookFrom + d)
        self.lookAt = self.lookAt + d
        self.fog.Set( start=self.fog.start-delta, end=self.fog.end-delta)
        self.Set( far=self.far-delta)
        ## print 'FAC', lFactor
        ## if (lFactor<0.005):
        ##     self.viewer.rootObject.Set(lineWidth=5, redo=True)
        ##     for g in self.viewer.rootObject.AllObjects():
        ##         g.Set(lineWidth=5, inheritLineWidth=0)
        ##     print 'LW 5'
        ## elif (lFactor<0.01):
        ##     self.viewer.rootObject.Set(lineWidth=4, redo=True)
        ##     for g in self.viewer.rootObject.AllObjects():
        ##         g.Set(lineWidth=4, inheritLineWidth=0)
        ##     print 'LW 4'
        ## elif (lFactor<0.03):
        ##     self.viewer.rootObject.Set(lineWidth=3, redo=True)
        ##     for g in self.viewer.rootObject.AllObjects():
        ##         g.Set(lineWidth=3, inheritLineWidth=0)
        ##     print 'LW 3'
        ## elif (lFactor<0.06):
        ##     self.viewer.rootObject.Set(lineWidth=2, redo=True)
        ##     for g in self.viewer.rootObject.AllObjects():
        ##         g.Set(lineWidth=2, inheritLineWidth=0)
        ##     print 'LW 2'
        ## else:
        ##     self.viewer.rootObject.Set(lineWidth=1, redo=True)
        ##     for g in self.viewer.rootObject.AllObjects():
        ##         g.Set(lineWidth=1, inheritLineWidth=0)
        ##     print 'LW 1'

        #print 'DELTA', delta, self.far
        ## fovy = self.fovy - self.fovy * -0.0018 * delta
        ## print delta, self.fovy, fovy, self.fovy-fovy
        ## if (fovy < 80.) and (fovy > 0.):
        ##     if (fovy<5.0):
        ##         self.viewer.rootObject.Set(lineWidth=5, redo=True)
        ##         for g in self.viewer.rootObject.AllObjects():
        ##             g.Set(lineWidth=5, inheritLineWidth=0)
        ##         print 'LW 5',
        ##     elif (fovy<10.0):
        ##         self.viewer.rootObject.Set(lineWidth=4, redo=True)
        ##         for g in self.viewer.rootObject.AllObjects():
        ##             g.Set(lineWidth=4, inheritLineWidth=0)
        ##         print 'LW 4',
        ##     elif (fovy<20.0):
        ##         self.viewer.rootObject.Set(lineWidth=3, redo=True)
        ##         for g in self.viewer.rootObject.AllObjects():
        ##             g.Set(lineWidth=3, inheritLineWidth=0)
        ##         print 'LW 3',
        ##     elif (fovy<30.0):
        ##         self.viewer.rootObject.Set(lineWidth=2, redo=True)
        ##         for g in self.viewer.rootObject.AllObjects():
        ##             g.Set(lineWidth=2, inheritLineWidth=0)
        ##         print 'LW 2',
        ##     elif (fovy<35.0):
        ##         self.viewer.rootObject.Set(lineWidth=1, redo=True)
        ##         for g in self.viewer.rootObject.AllObjects():
        ##             g.Set(lineWidth=1, inheritLineWidth=0)
        ##         print 'LW 1',
        ##     else:
        ##         self.viewer.rootObject.Set(lineWidth=1, redo=True)
        ##         for g in self.viewer.rootObject.AllObjects():
        ##             g.Set(lineWidth=1, inheritLineWidth=0)
        ##         print 'LW 1',
        ##     self.fovy = fovy
            
        
    def scaleSlab(self, delta):
        """change near and far clipping planes to make slab thicker or thinner
        """
        new_near = max(0.1, self.near - 0.5*delta)
        new_far = self.far + 0.5*delta
        if new_near < new_far:
            self.Set(near=new_near, far=new_far)
            #print 'scale slab', delta, new_near, new_far
        
    def translateSlab(self, delta):
        new_near = max(0.1, self.near + 0.5*delta)
        new_far = self.far + 0.5*delta
        self.Set(near=new_near, far=new_far)
        #print 'translate slab', delta, new_near, new_far
        
    ##
    ## GUI callbacks for changing the camera
    ##
    def getDelta(self, dx, dy, e, normalize=True):
        if isinstance(e, QtGui.QMouseEvent):
            #if self._motionAxis==0: delta = dx
            #else: delta = dy
            if abs(dx) > abs(dy): delta = dx
            else: delta = dy
        elif isinstance(e, QtGui.QWheelEvent):
            delta = e.delta() / 120
        if normalize:
            if delta > 0: delta = 1
            else: delta = -1
        #print dx, dy, delta, isinstance(e, QtGui.QMouseEvent)
        return delta
    
    def handleZoomCamera(self, x, y, dx, dy, e ):
        self.zoomCamera( self.getDelta(dx,dy,e,normalize=False))
        self.viewer.OneRedraw()

    def handleScaleSlab(self, x, y, dx, dy, e):
        self.scaleSlab( self.getDelta(dx,dy,e) )
        self.viewer.OneRedraw()

    def handleTranslateSlab(self, x, y, dx, dy, e):
        self.translateSlab( self.getDelta(dx,dy,e) )
        self.viewer.OneRedraw()

    def showTools(self, x, y, e):
        #print 'showTools', x, y
        pass
    ##
    ## selection drag GUI methods
    ##
    def initSelectionRectangle(self, x, y):
        if __debug__:
            if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        if isinstance(self.viewer.currentObject, Insert2d):
            return
        if self.stereoMode != 'MONO':
            return
        self.selectDragRect = 1
        self.afid=None
        self.fill=0
        self.Activate()
        glPushMatrix()
        self.BuildTransformation()  # camera transformation
        glDrawBuffer(GL_FRONT)
#        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
#        glDisable(GL_CULL_FACE)
        glLineWidth( 1.0 )
        glDisable( GL_LIGHTING )
        glColor3f( float(self.selectionColor[0]), float(self.selectionColor[1]),
                   float(self.selectionColor[2]) )
        glDisable(GL_DEPTH_TEST)
        glEnable(GL_COLOR_LOGIC_OP)
        glLogicOp(GL_XOR)
        x2 = self._x1 = x
        y2 = self._y1 = self._height - y

        self.unProj_model = glGetDoublev(GL_MODELVIEW_MATRIX)
        self.unProj_proj = glGetDoublev(GL_PROJECTION_MATRIX)
        self.unProj_view = glGetIntegerv(GL_VIEWPORT)

        from opengltk.extent import _glulib as glulib
        self._P1 = gluUnProject( (self._x1, self._y1, 0.5),
                                 self.unProj_model, self.unProj_proj,
                                 self.unProj_view )
        self._P2 = gluUnProject( (self._x1, y2, 0.5),
                                 self.unProj_model, self.unProj_proj,
                                 self.unProj_view )
        self._P3 = gluUnProject( (x2, y2, 0.5),
                                 self.unProj_model, self.unProj_proj,
                                 self.unProj_view )
        self._P4 = gluUnProject( (x2, self._y1, 0.5),
                                 self.unProj_model, self.unProj_proj,
                                 self.unProj_view )

        # draw first rectangle
        #self.history = [ ('draw', self._P3) ]
        self.drawRect( self._P1, self._P2, self._P3, self._P4 )
        glPopMatrix()


    def addToSelection(self, x, y, dx, dy, e):
        if not hasattr(self, '_firstCallDone'):
            self.initSelectionRectangle(x, y)
            self._firstCallDone = True
        else:
            if not self.selectDragRect: return
            self.Activate()
            glPushMatrix()
            self.BuildTransformation()  # camera transformation
            #print 'AAAAAAAAAAA2 drawSelectionRectangle'

            # draw over previous rectangle
            #self.history.append( ('hide', self._P3) )
            glDisable(GL_DEPTH_TEST)
            if self.fill:
                self.drawRect( self._P1, self._P2, self._P3, self._P4, 1 )
                self.fill = 0
            else:
                self.drawRect( self._P1, self._P2, self._P3, self._P4, 0 )

            # draw new rectangle
            x2 = x
            y2 = self._height - y

            self.unProj_model = glGetDoublev(GL_MODELVIEW_MATRIX)

            self._P2 = gluUnProject( (self._x1, y2, 0.5),
                                     self.unProj_model, self.unProj_proj,
                                     self.unProj_view)
            self._P3 = gluUnProject( (x2, y2, 0.5),
                                     self.unProj_model, self.unProj_proj,
                                     self.unProj_view)
            self._P4 = gluUnProject( (x2, self._y1, 0.5),
                                     self.unProj_model, self.unProj_proj,
                                     self.unProj_view)
            #self.history.append( ('draw', self._P3) )
            if self._P1[0] < self._P3[0]:
                self.setCursor(QtCore.Qt.DragCopyCursor)
                ## cursor_px = QtGui.QPixmap('/mgl/ms1/people/sanner/python/lazy/PmvApp/GUI/Icons/selectionAdd.png')
                ## cursor_px.setMask(cursor_px.mask())
                ## cursor = QtGui.QCursor(cursor_px)
                ## self.setCursor(cursor)
                
            else:
                self.setCursor(QtCore.Qt.DragLinkCursor)
                ## cursor_px = QtGui.QPixmap('/mgl/ms1/people/sanner/python/lazy/PmvApp/GUI/Icons/selectionRemove.png')
                ## cursor_px.setMask(cursor_px.mask())
                ## cursor = QtGui.QCursor(cursor_px)
                ## self.setCursor(cursor)
                
            self.drawRect( self._P1, self._P2, self._P3, self._P4, self.fill )
            glPopMatrix()
            glFlush()

    def endSelectionRectangle(self, x, y, dx, dy, e):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        #print "Camera.endSelectionRectangle"

        self.setCursor(QtCore.Qt.ArrowCursor)

        del self._firstCallDone
        if not self.selectDragRect: return
        self.selectDragRect = 0
        self.Activate()
        glPushMatrix()
        self.BuildTransformation()  # camera transformation

        # this line is required for last rectangle to disappear !
        # not sure why
        # oct 2001: apparently not necessary
        #glEnable(GL_COLOR_LOGIC_OP)

        #self.history.append( ('hide', self._P3) )
        self.drawRect( self._P1, self._P2, self._P3, self._P4 )
            
        #self.drawRect( self._P1, self._P2, self._P3, self._P4 )
        glPopMatrix()
        
        glEnable(GL_DEPTH_TEST)
        glDisable(GL_COLOR_LOGIC_OP)
        #glEnable(GL_LIGHTING)
        if self.viewer is not None:
            self.viewer.enableOpenglLighting()
        glDrawBuffer(GL_BACK)
        pick = self.DoPick(x, y, self._x1, self._height-self._y1, event=e)
        del self._x1
        del self._y1
        del self._P1
        del self._P2
        del self._P3
        del self._P4
        if self.viewer and len(pick.hits):
            self.viewer.processPicking(pick)


    def rotateScene(self, x, y, dx, dy, e):
        self._rotateObject(self.viewer.rootObject, x, y, dx, dy, e)
        
    def rotateScene(self, x, y, dx, dy, e):
        self._rotateObject(self.viewer.currentObject, x, y, dx, dy, e)

    def _rotateObject(self, obj, x, y, dx, dy, e):
        # Rotate like a trackball
	self.tb.update(self._oldx, self._oldy, e.x(), e.y(), self._width, self._height, 1)
        self.viewer.rootObject.ConcatRotation(self.tb.mat)
        self.viewer.Redraw()


    def translateSceneXY(self, x, y, dx, dy, e):
        self._translateObjectXY(self.viewer.rootObject, x, y, dx, dy, e)

    def translateObjectXY(self, x, y, dx, dy, e):
        self._translateObjectXY(self.viewer.currentObject, x, y, dx, dy, e)

    def _translateObjectXY(self, obj, x, y, dx, dy, e):
        lFovyFactor = self.fovy / self.fovyNeutral
        if lFovyFactor > 1:
            lFovyFactor = pow(lFovyFactor, 1.5)
        lDistanceFactor = (self.lookFrom[2] - 0.0)  * .001

        lDistanceFactor = exp(lDistanceFactor) - 1
        lFactor = lDistanceFactor * lFovyFactor        
        t =  (dx * lFactor, dy*lFactor, 0.0 )
        obj.ConcatTranslation( t )
        self.viewer.Redraw()

    ## def scaleScene(self, x, y, dx, dy, e):
    ##     self._scaleObject(self.viewer.rootObject, x, y, dx, dy, e)

    ## def scaleObject(self, x, y, dx, dy, e):
    ##     self._scaleObject(self.viewer.currentObject, x, y, dx, dy, e)

    ## def _scaleObject(self, obj, x, y, dx, dy, e):
    ##     if abs(dx) > abs(dy): delta = dx
    ##     else: delta = dy
    ##     oldPivot = self.currentObject.SetPivot( self.lastPickedPixe )
    ##     obj = self.currentObject
    ##     self.ScaleCurrentObject(self, None, None, None, delta)


    ## def translateCameraZ(self, x, y, dx, dy, e):
    ##     # move camera to cerate zoom in perspective
    ##     # clipping planes and depthcueing are conserved
    ##     # BUG: - when near gets below 0 it is clamped at 0 but does
    ##     #      not recover when we move back
    ##     #      - the motion should be done along the vector from
    ##     #      lookFrom to lookAt rather thatn on the Z axis
    ##     if abs(dx) > abs(dy): delta = dx
    ##     else: delta = dy
    ##     t =  (0, 0, -delta )
    ##     #self.ConcatTranslation( t )
    ##     self.lookAt[2] -= delta
    ##     self.lookFrom[2] -= delta
    ##     #print "ASDSADASA", self.lookFrom, self.lookAt, self.fog.start, self.fog.end
    ##     self.fog.Set(start=self.fog.start-delta, end=self.fog.end-delta)
    ##     self.Set(near=self.near-delta, far=self.far-delta)
    ##     self.viewer.Redraw()

    def translateSceneZ(self, x, y, dx, dy, e):
        self._translateObjectZ(self.viewer.rootObject, x, y, dx, dy, e)

    def translateObjectZ(self, x, y, dx, dy, e):
        self._translateObjectZ(self.viewer.currentObject, x, y, dx, dy, e)

    def _translateObjectZ(self, obj, x, y, dx, dy, e):
        if isinstance(obj, Transformable):
            lFovyFactor = self.fovy / self.fovyNeutral
            if lFovyFactor > 1:
                lFovyFactor = pow(lFovyFactor, 1.5)
            lDistanceFactor = (self.far - self.near)  * .001
            lDistanceFactor = exp(lDistanceFactor) - 1
            lFactor = lDistanceFactor * lFovyFactor * .5
            if abs(dx) > abs(dy): delta = dx
            else: delta = dy
            t =  (0.0, 0.0, -delta*lFactor)
            obj.ConcatTranslation( t )
            self.viewer.Redraw()


    def pick(self, x, y, dx, dy, e):
        pick = self.DoPick(x, y)
        self.viewer.processPicking(pick)
        
    def enterEvent(self, e):
        #print "enter"
        self.setFocus()

    def leaveEvent(self, e):
        #print "leave"
        self.clearFocus()

    def Activate(self):
        self.makeCurrent()

    def getContext(self):
        return self.context()
        
    def SwapBuffers(self):
        self.swapBuffers()

    def wheelEvent(self, e):
        #numDegrees = event.delta() / 8
        #numSteps = event.delta() / 120
        func = self._mouseWheelActions[int(e.modifiers())]
        if func:
            func(None, None, e.delta(), e.delta(), e)

    def mousePressEvent(self, e):
        #print "mouse press"
        #print 'func modifier', int(e.buttons()), int(e.modifiers())
        # save the mouse button as it is lost on the release event
        self._isPressed = int(e.buttons())
        func = self._mousePressActions[self._isPressed][int(e.modifiers())]
        if func:
            func(x, y, 0, 0, e)
        self._hasMoved = False
        self._orginalx = self._oldx = e.x()
        self._orginaly = self._oldy = e.y()
        for func in self.onButtonDownCBlist:
            func(e)
        
    def mouseMoveEvent(self, e):
        x = e.x()
        y = e.y()
        if int(e.buttons()) == QtCore.Qt.NoButton :
            func = self._mouseMoveNoButtonActions[int(e.modifiers())]
            if func:
                func(x, y, e)
        else:
            delta_x = x - self._oldx
            delta_y = self._oldy - y
            self._hasMoved = True
            #if self._motionAxis is None:
            #    if abs(delta_x) > abs(delta_y):
            #        self._motionAxis = 0 # X axis
            #    else:
            #        self._motionAxis = 1
            #print 'MOVE', self._isPressed, int(e.modifiers())
            func = self._mouseMoveActions[self._isPressed][int(e.modifiers())]
            if func:
                func(x, y, delta_x, delta_y, e)
            self.update()
            self._oldx = x
            self._oldy = y

    def mouseReleaseEvent(self, e):
        #self._motionAxis = None
        x = e.x()
        y = e.y()
        delta_x = x - self._oldx
        delta_y = self._oldy - y
        if self._hasMoved:
            func = self._mouseReleaseWithMotionActions[self._isPressed][int(e.modifiers())]
            #print "after moving", func, self._isPressed, int(e.modifiers())
            if func:
                func(x, y, delta_x, delta_y, e)
        else:
            #print 'func modifier', self._isPressed, int(e.modifiers())
            func = self._mouseReleaseNoMotionActions[self._isPressed][int(e.modifiers())]
            if func:
                func(x, y, delta_x, delta_y, e)

        for func in self.onButtonUpCBlist:
            func(e)

    def keyPressEvent(self, e):
        #print "key press", e.key(), e.isAutoRepeat(), e.text()
        shift = control = alt = meta = keypad = False
        if e.modifiers() == QtCore.Qt.NoModifier:
            noModifier = True
        else:
            noModifier = False
            if e.modifiers() & QtCore.Qt.ShiftModifier: shift = True
            if e.modifiers() & QtCore.Qt.ControlModifier: control = True
            if e.modifiers() & QtCore.Qt.AltModifier: alt = True
            if e.modifiers() & QtCore.Qt.MetaModifier: meta= True
            if e.modifiers() & QtCore.Qt.KeypadModifier: keypad = True
            

        #print 'MODIOFIERS', noModifier, shift, control, alt, meta, keypad  
        if noModifier:
            if e.text()=='a': self.viewer.AutoDepthCue()
            elif e.text()=='A': self.viewer.AutoDepthCue()
            elif e.text()=='c': self.viewer.Center_cb()
            elif e.text()=='C': self.viewer.Center_cb()
            elif e.text()=='d': self.viewer.ToggleDepth()
            elif e.text()=='D': self.viewer.ToggleDepth()
            elif e.text()=='n': self.viewer.Normalize_cb()
            elif e.text()=='N': self.viewer.Normalize_cb()
            elif e.text()=='o': self.viewer.SSAO_cb_arg()
            elif e.text()=='O': self.viewer.SSAO_cb_arg()
            elif e.text()=='r': self.viewer.Reset_cb()
            elif e.text()=='R': self.viewer.Reset_cb()
            elif e.text()=='q': self.saveLargeImage()
            elif e.text()=='f': self.keyPressSignal.emit('f')
            elif e.text()=='F': self.keyPressSignal.emit('f')
        self.viewer.Redraw()


    ## def keyReleaseEvent(self, e):
    ##     print "key release", e.key(), e.isAutoRepeat(), e.text()


    def focusInEvent(self, e):
        #print 'focus In'
        pass

    def focusOutEvent(self, e):
        #print 'focus Out'
        pass
        
    def storeGeometry(self):
        #print 'ASASASASA'
        g = self.geometry()
        self.rootx = g.left()
        self.rooty = g.top()
        #self._width = self.width()
        #self._height = self.height()
        
    def initializeGL(self):
        glClearColor(0.0, 0.0, 0.0, 1.0)
        glClearDepth(1.0)

    def paintGL(self):
        self.exposeEvent = True
        self.viewer.needsRedraw = True
        self.viewer.ReallyRedraw()
        
    def resizeGL(self, widthInPixels, heightInPixels):
        #print 'RESIZE', widthInPixels, heightInPixels, self.exposeEvent, self.viewer.autoRedraw
        # if viewer is in autoRedraw mode the next redraw will handle it
        #print self._width, widthInPixels, self._height, heightInPixels, self.width(), self.height()
        self._width = widthInPixels
        self._height = heightInPixels
        if not self.exposeEvent and self.viewer.autoRedraw:
            self.viewer.Redraw()
            self.exposeEvent = True

        for o in self.viewer.rootObject.AllObjects():
            if o.needsRedoDpyListOnResize or o.scissor:
                self.viewer.objectsNeedingRedo[o] = None
        
        #self.camera.setViewportDimensions(widthInPixels, heightInPixels)
        #self.Activate()
        #glViewport(0, 0, widthInPixels, heightInPixels)

    ## def getGeometry(self):
    ##     raise

    ## def Geometry(self):
    ##     raise
    
    ## def lift(self):
    ##     raise

    ## def Enter_cb(self, event=None):
    ##     raise

    ## def activeStereoSupport(self):
    ##     raise
    def saveLargeImage(self):
        # does not work yet. raises an eception in amera.py", line 3913, in RedrawAASwitch
        # glReadBuffer(GL_BACK)

        pixmap = self.renderPixmap(w=2000, h=2000)#, useContext=false]]])
        pixmap.save('large.png', 'PNG') 
        #pixmap = QPixmap()
        #bytes = QByteArray()
        #buffer(bytes)
        #buffer.open(QIODevice.WriteOnly)
        #pixmap.save(buffer, "PNG") # writes pixmap into bytes in PNG format
