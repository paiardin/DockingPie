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
# $Header: /mnt/raid/services/cvs/DejaVu2/Qt/Viewer.py,v 1.1.1.1.4.1 2017/07/13 22:25:45 annao Exp $
#
# $Id: Viewer.py,v 1.1.1.1.4.1 2017/07/13 22:25:45 annao Exp $
#
import threading

from PySide import QtGui, QtCore

import DejaVu2
from DejaVu2.Qt.Camera import Camera
from DejaVu2.Viewer import ViewerBase

import numpy, math

#class Viewer(QtGui.QDockWidget, ViewerBase):
#class Viewer(QtGui.QMainWindow, ViewerBase):
#class Viewer(QtGui.QWidget, ViewerBase):
#class Viewer(QtGui.QWidget, ViewerBase):
class Viewer(ViewerBase):

    def __init__(self, nogui=0,
                 guiMaster=None,
                 showViewerGUI=True,
                 autoRedraw=True,
                 cnf={}, **kw):
	"""Viewer's constructor
"""
        if guiMaster is not None:
            self.ownMaster = True
        #print 'VIEWER MASTER', id(master)
        #QtGui.QDockWidget.__init__(self, master)
        #QtGui.QMainWindow.__init__(self, master)
        #self.setWindowTitle('DejaVu2_Viewer_Qt')        
        #QtGui.QWidget.__init__(self, master)

        #self.fileMenu = self.menuBar().addMenu(self.tr("&File"))
        #self.fileMenu.hide()

        #self.toolBar = self.addToolBar(self.tr("Focus"))

        #self.redrawTimer = QtCore.QTimer()
        #self.connect(self.redrawTimer, QtCore.SIGNAL("timeout()"),
        #             self.ReallyRedraw)

        # string var used to decide what the trackball is moving
        self.Xform = 'Object'

        self.contourTk = False

        self.spinVar = DejaVu2.defaultSpinningMode
        self.spinCallBack = None

        # Decides if the call to enable GL_LIGHTNING will be considered 
        self.OverAllLightingIsOn = 1
        
        # not sure about this but if it is not there I have 100x3 black box in upper left corner
        #mainLayout = QtGui.QGridLayout()
        #self.setLayout(mainLayout)
        ViewerBase.__init__(self, nogui=0, #screenName=None,
                            guiMaster=None, #classCamera=None, 
                            autoRedraw=True,
                            #verbose=True,
                            cnf=cnf, **kw)

        #objects = QtGui.QDockWidget('ObjectTree', self)
        #from tree import ObjectTree
        #self.objTree = ObjectTree(self, parent=objects)
        #objects.setWidget(self.objTree)
        
        #from dashboard import Dashboard
        #self.objTree = Dashboard(parent=objects)
        #objects.setWidget(self.objTree)
        #mainLayout.addWidget(objects)
        
        #self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, objects)

        #mainLayout.addWidget(self.cameras[0].master)

        # create the material editor
        self.materialEditor = None
    ##     self.materialEditor = MaterialEditor(None, self)
    ##     self.materialEditor.dismissTk.grid_forget()
    ##     self.materialEditor.dismiss()

    ##     # register Object that have their own OpenGL context but need to have
    ##     # the same lighting model and lights
    ##     for l in self.lights:
    ##         l.applyTo = [self.materialEditor]

    ##     self.lightModel.applyTo = [self.materialEditor]
        
    ##     # create the ViewerGui
    ##     if showViewerGUI:
    ##         self.GUI = ViewerGUI(self, self.maxLights, self.maxClipP,
    ##                              nogui=nogui, master=guiMaster)

    ##         #self.GUI.CameraBackgroundColor = self.CurrentCameraBackgroundColor
    ##         #self.GUI.LightColor = self.CurrentLightColor
    ##         #self.GUI.ClipColor = self.CurrentClipColor
    ##         #self.GUI.ObjectFrontColor = self.CurrentObjectFrontColor
    ##         #self.GUI.ObjectBackColor = self.CurrentObjectBackColor
    ##         #self.GUI.LightModelColor = self.LMColor

    ##         self.GUI.addObject(self.rootObject, None)

    ##         self.GUI.bindResetButton( self.Reset_cb)
    ##         self.GUI.bindNormalizeButton( self.Normalize_cb)
    ##         self.GUI.bindCenterButton( self.Center_cb)
    ##         self.GUI.bindDeleteButton( self.Delete_cb)
    ## ##         self.GUI.Exit = self.__del__

            ## if nogui and isinstance(self.GUI.root, Tkinter.Toplevel):
            ##     self.GUI.withdraw()

        #self.GUI.addObject(self.pickVerticesSpheres, None)
        #self.GUI.showPickedVertex.set(self.showPickedVertex)

    def FocusOnBox(self, mini, maxi):
        """Moves camera to provided as [mini, maxi] bounding box"""

        g = numpy.sum( (mini, maxi),0 ) * .5 # center of BoundingBox        
        self.rootObject.Set(scale=(1.,1.,1.))
        lBox = maxi - mini
        self.lHalfObject = max(lBox)/2.
        if self.lHalfObject == 0.:
            self.lHalfObject = 1.
        cam = self.currentCamera
        self.oldLookFrom = cam.lookFrom[2]

        self.diff_fovy =  cam.fovyNeutral - cam.fovy
        cam.lookAt = numpy.array((0.,0.,0.))
        Rmini, Rmaxi = self.rootObject.ComputeBB(self.currentCamera)
        Rg = numpy.sum( (Rmini, Rmaxi), 0) * .5

        self.diffVect = -g-self.rootObject.translation
        if not (g-Rg).any():
            self.diffVect = -Rg

        cam.Set(fov=cam.fovy + self.diff_fovy)
        self.rootObject.ConcatTranslation( self.diffVect[:3])
        lNewDist = self.lHalfObject / math.tan(cam.fovy/2*math.pi/180.0)
        newDist = cam.nearDefault+lNewDist+self.lHalfObject
        dist = self.oldLookFrom + newDist
        cam.lookFrom = numpy.array( ( 0., 0., dist ) )            
        cam.direction = cam.lookAt - cam.lookFrom 

        currentObject = self.currentObject
        self.currentObject = self.rootObject
        self.CenterCurrentObject()
        self.currentObject = currentObject 
        self.currentCamera.AutoDepthCue(object=self.rootObject)
        self.OneRedraw()
        g = ( (0.5*(mini[0]+maxi[0])), (0.5*(mini[1]+maxi[1])),
              (0.5*(mini[2]+maxi[2])) )
        self.rootObject.SetPivot(g)
        self.OneRedraw()
        
    def postNextRedraw(self):
        if self.autoRedraw:
            if self.redrawTimer.isActive():
                self.redrawTimer.stop()
            self.redrawTimer.start(100)


    def startAutoRedraw(self):
        self.autoRedraw = True
        self.redrawTimer.start(100)


    def stopAutoRedraw(self):
        if self.redrawTimer.isActive():
            self.redrawTimer.stop()
        self.autoRedraw = False   


    def checkIfRedrawIsNeeded(self):
        if self.suspendRedraw:
            self.redrawTimer.start(1000)
            #print 'NO REDRAW SUSPEND'
            return False

        if not self.needsRedraw:
            self.redrawTimer.start(100)
            #print 'NO REDRAW no needs'
            return False
                    
        if threading.currentThread().getName()!='MainThread':
            print 'NO REDRAW not main thread'
            #self.redrawTimer.start(100)
            return False

        if self.autoRedraw and not self.redrawTimer.isActive():
            #print 'NO REDRAW no active timer'            
            return False
        return True
    

    def closeEvent(self, event):
        reply = QtGui.QMessageBox.question(
            self, "Confirmation",
            "Are you sure you want to quit?",
            QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.No)
        if reply == QtGui.QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()


    def AddCamera(self, master=None, screenName=None, classCamera=None,
                  stereo='none', num=None, verbose=True, cnf={}, **kw):
	"""Add one more camera to this viewer"""

        if num is None:
            num = len(self.cameras)

        if classCamera is None:
            classCamera = Camera
            name = 'camera '+str(num)
        else:
            name = classCamera.__name__+str(num)
            
        cameraOwnsMaster = False
	if not master:
            cameraOwnsMaster = True
        
        kw['stereo'] = stereo
            
        c = classCamera(master, screenName, self, num, check=1, cnf=cnf, **kw)
        c.show()
        c.fog.Set(enabled=True)
        #if hasattr(c.frame.master,"protocol"):
        #    c.frame.master.protocol("WM_DELETE_WINDOW",self.closeEvent)
        ## c.eventManager.AddCallback('<KeyPress>', self.modifierDown)
        ## c.eventManager.AddCallback('<KeyRelease>', self.modifierUp)
        ## c.eventManager.AddCallback('R', self.Reset_cb_arg)
        ## c.eventManager.AddCallback('r', self.Reset_cb_arg)
        ## c.eventManager.AddCallback('A', self.AutoDepthCue)
        ## c.eventManager.AddCallback('a', self.AutoDepthCue)
        ## c.eventManager.AddCallback('N', self.Normalize_cb_arg)
        ## c.eventManager.AddCallback('n', self.Normalize_cb_arg)
        ## c.eventManager.AddCallback('C', self.Center_cb_arg)
        ## c.eventManager.AddCallback('c', self.Center_cb_arg)
        ## c.eventManager.AddCallback('D', self.Depth_cb_arg)
        ## c.eventManager.AddCallback('d', self.Depth_cb_arg)
        ## c.eventManager.AddCallback('T', self.toggleTransformRootOnly)
        ## c.eventManager.AddCallback('t', self.toggleTransformRootOnly)
        ## c.eventManager.AddCallback('L', self.toggleOpenglLighting)
        ## c.eventManager.AddCallback('l', self.toggleOpenglLighting)
        ## c.eventManager.AddCallback('O', self.SSAO_cb_arg)
        ## c.eventManager.AddCallback('o', self.SSAO_cb_arg)    
        c.ownMaster = cameraOwnsMaster
        
        ## if self.GUI is not None:
        ##     self.GUI.bindModifersToTransformInfo(master)

	self.cameras.append( c )
        #c.frame.config( background = "#900000" )

	# make the trackball transform the current object
        #if self.rootObject:
        #    self.BindTrackballToObject(self.rootObject)

        c.firstRedraw = True
        c.Activate()

	if len(self.cameras)==1:
            from opengltk.OpenGL import GL
            from DejaVu2.Clip import ClippingPlane
            from DejaVu2.Light import Light
            
	    self.currentCamera = c
            c.hasBeenCurrent = 1

            self.SetCurrentObject(self.rootObject)

            #self.GLextensions = string.split(GL.glGetString(GL.GL_EXTENSIONS))
            #self.GLversion = GL.glGetString(GL.GL_VERSION)
            #self.GLvendor = GL.glGetString(GL.GL_VENDOR)
            if verbose:
                print 'OpenGL-based graphics'
                print ' GL_VENDOR', GL.glGetString(GL.GL_VENDOR)
                print ' GL_RENDERER', GL.glGetString(GL.GL_RENDERER)
                print ' GL_VERSION', GL.glGetString(GL.GL_VERSION)
        
            #self.hasOffsetExt = "GL_EXT_polygon_offset" in self.GLextensions

            # Find out how many clipping planes OpenGL defines
            maxClipP = 0
            for i in range(6):
                try:
                    cpn = getattr(GL, "GL_CLIP_PLANE%d"%i)
                    ClippingPlane.clipPlaneNames.append(cpn)
                    maxClipP +=1
                except AttributeError:
                    break
            self.maxClipP = maxClipP

            # Find out how many clipping planes OpenGL "thinks" it supports
            #maxClipP = min(int(GL.glGetDoublev(GL.GL_MAX_CLIP_PLANES)[0]), 6)

            # Find out how many light sources
            maxLights = int(GL.glGetDoublev(GL.GL_MAX_LIGHTS)[0])
            if maxLights > 8:
                print 'WARNING: Reducing number of light sources from %d to 8' % \
                    maxLights
                maxLights = 8
            self.maxLights = maxLights

            # create the light sources
            for i in range(self.maxLights):
                l = Light(i, self)
                self.lights.append( l )
                l.viewer = self

            self.InitLighting(c)

            self.rootObject.clipSide = [1]*self.maxClipP

            for i in range(self.maxClipP):
                #from DejaVu2.csgClip import CsgClippingPlane
                #cp = CsgClippingPlane(self.rootObject, i, self)
                cp = ClippingPlane(self.rootObject, i, self)
                self.clipP.append(cp)

            if self.maxClipP:
                self.currentClip = self.clipP[0]
            else:
                self.currentClip = None

            # see if it PolygonOffset actually implemented
            if self.hasOffsetExt:
                try:
                    self.polyOffset(1.0, 1.0)
                except ValueError:
                    self.viewer.hasOffsetExt = False
                
##         if cnf.has_key("addScenarioButton"):
##             self.addScenarioButton = cnf["addScenarioButton"]
##         else:
##             self.addScenarioButton = True

            # this line causes networks wih viewers not to appear completely
            # (missing connections, which appear only after a window takes focus
            # MS march 26, 03
            #self.master.wait_visibility()
                    
            #print "DejaVu2.enableStereo", DejaVu2.enableStereo
            if DejaVu2.enableStereo is True:
                self.activeStereoSupport = self.currentCamera.activeStereoSupport()
            else:
                self.activeStereoSupport = False

            #print "self.activeStereoSupport", self.activeStereoSupport

            #make sure that at the end of the init the current context is the one of the current camera
            self.currentCamera.Redraw()
            #self.currentCamera.Activate()

            # create fake geometry to find out if VBO are supported
            from DejaVu2.IndexedPolygons import IndexedPolygons
            fake = IndexedPolygons('test', vertexArrayFlag=True,
                                   vertices=((0,0,0),))

            if verbose:
                print 'Enable VBO:', DejaVu2.enableVBO

            self.suspendRedraw = False # set to True to Prevent redrawing

            if self.autoRedraw:
                self.pendingAutoRedrawID = self.redrawTimer.start(10)

	return c


    ## def AddObject(self, obj, parent=None, redraw=True, redo=True):
    ##     ViewerBase.AddObject(self, obj, parent=parent, redraw=redraw, redo=redo)
    ##     self.objTree.addObject(obj.parent, obj, obj.name, str(obj))

##     def _DeleteCamera(self, camera):
##         """Remove the given camera in the right order
## """
##         #print 'Viewer._DeleteCamera ', camera
##         # Check if this camera shareCTX with other camera.
##         if hasattr(camera, 'shareCTXWith'): 
##             while len(camera.shareCTXWith):
##                 cam = camera.shareCTXWith[0]
##                 self._DeleteCamera(cam)
##         camera.destroy()
##         if camera.ownMaster:
##             camera.frame.master.destroy()
##         else:
##             camera.frame.destroy()
##         self.cameras.remove(camera)
##         for c in self.cameras:
##             if hasattr(c, 'shareCTXWith'):
##                 c.shareCTXWith.remove(camera)

        
##     def DeleteCamera(self, camera):
##         """
##         Remove the given camera from the viewer and takes care
##         of the dpyList if camera is cameras[0]
##         """
## #        #delete NPR rendering toplevel
## #        if camera.imCanvastop:
## #            camera.imCanvastop.destroy()
## #            camera.imCanvas = None
## #            camera.imCanvas1 = None
## #            camera.imCanvas2 = None
## #            camera.imCanvas3 = None
            
##         # Remove the camera from the list of cameras associated to
##         # the viewer.
##         camIndex = self.cameras.index(camera)

##         # the current openGL context has been destroyed so
##         # the dpyList need to be destroyed only if the CTX is not
##         # shared by any other camera.
##         if camIndex == 0:
##             self.objectsNeedingRedo = {}
##             for g in self.rootObject.AllObjects():
##                 g.deleteOpenglList()
##                 self.objectsNeedingRedo[g] = None

##         self._DeleteCamera(camera)
##         # If this camera is the current camera
##         if self.currentCamera == camera:
##             if len(self.cameras) == 0:
##                 # There is no more cameras then set currentCamera to None
##                 self.currentCamera = None
##             else:
##                 # Set the current Camera to be the first camera of the list.
##                 self.currentCamera = self.cameras[0]
