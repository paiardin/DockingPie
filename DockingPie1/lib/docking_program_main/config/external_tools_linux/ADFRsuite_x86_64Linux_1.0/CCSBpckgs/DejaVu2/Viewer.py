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
# Date: 2000 Authors: Michel Sanner
#
#    vareille@scripps.edu
#    sanner@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Michel Sanner and TSRI 2000
#
# revision: Guillaume Vareille
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/DejaVu2/Viewer.py,v 1.2.4.2 2017/11/08 23:27:21 annao Exp $
#
# $Id: Viewer.py,v 1.2.4.2 2017/11/08 23:27:21 annao Exp $
#

import warnings, thread, threading, os, numpy, math, string, types, sys, re
from time import time, sleep

from PIL import Image

from opengltk.OpenGL import GL

from mglutil.events import Event, EventHandler
from mglutil.util.callback import CallbackFunction
from mglutil.math.rotax import rotVectToVect

import DejaVu2
from DejaVu2 import viewerConst
from DejaVu2.Geom import Geom
from DejaVu2.Camera import Camera
from DejaVu2.Clip import ClippingPlane
from DejaVu2.Light import Light, LightModel
from DejaVu2.colorTool import OneColor
from DejaVu2.tileRenderer import TRctx
from DejaVu2.Transformable import Transformable 
from DejaVu2.Insert2d import Insert2d 
from DejaVu2.Common2d3dObject import Common2d3dObject
from DejaVu2.cursors import cursorsDict

# count clipping planes

if hasattr( DejaVu2, 'enableStereo') is False:
    DejaVu2.enableStereo = False


class SetCurrentObjectEvent(Event):

    def __init__(self, object):
        Event.__init__(self)
        self.object = object

class SetCurrentCameraEvent(SetCurrentObjectEvent):
    pass

class SetCurrentLightEvent(SetCurrentObjectEvent):
    pass

class SetCurrentClipEvent(SetCurrentObjectEvent):
    pass


class AddObjectEvent(Event):
    def __init__(self, object):
        Event.__init__(self)
        self.object = object


class RemoveObjectEvent(AddObjectEvent):
    pass


class  ReparentObjectEvent(Event):
    
    def __init__(self, object, oldparent):
        Event.__init__(self)
        self.object = object
        self.oldparent = oldparent



class ViewerBase(object):
    """A viewer holds a bunch of object that inherit from Geom and displays
       them in cameras"""


    def getViewerStateDefinitionCode(self, viewerName, indent='',
                                     withMode=True, rootxy=True):
        """Returns a list of strings containing Python code to restore the
state of the viewer. The returned code will restore the state of all non-
transient objects, i.e. light models, light sources, clipping planes,
cameras, fog and root object.
indent is the level of indentation that is desired.
viewerName is the name of the Viewer object in the source code generated.
"""
        lines = []
        
        if withMode is True:
            lines.append(indent+"if mode=='viewer' or mode=='both':\n")
            indent += '    '

        lines.append(indent+"##\n")
        lines.append(indent+"    ## Saving State for Viewer\n")
        lines.append(indent+"%s.TransformRootOnly(%d)\n"%(viewerName, self.redirectTransformToRoot))
        lines.append(indent+"##\n\n")
        lines.append(indent+"## Light Model\n")
        lm = self.lightModel
        if lm._modified:
            lines.append(indent+"state = %s\n"%str(lm.getState()))
            lines.append(indent+"%s.lightModel.Set(**state)\n\n"%(
                viewerName))
        lines.append(indent+"## End Light Model\n\n")

        lines.append(indent+"## Light sources\n")
        for i,l in enumerate(self.lights):
            if not l._modified: continue
            lines.append(indent+"## Light Number %d\n"%i)
            lines.append(indent+"state = %s\n"%str(self.lights[i].getState()))
            lines.append(indent+"%s.lights[%d].Set(**state)\n\n"%(
                viewerName, i))
        lines.append(indent+"## End Light sources %d\n\n"%i)
        
        lines.append(indent+"## Cameras\n")
        for i,c in enumerate(self.cameras):
            if not c._modified: continue
            lines.append(indent+"## Camera Number %d\n"%i)
            lCamState = self.cameras[i].getState()
            if rootxy is False:
                lCamState.pop('rootx')
                lCamState.pop('rooty')
            lines.append(indent+"state = %s\n"%str(lCamState))
            lines.append(indent+"%s.cameras[%d].Set(**state)\n\n"%(
                viewerName, i))
            lines.append(indent+"state = %s\n"%str(
                self.cameras[i].fog.getState()))
            lines.append(indent+"%s.cameras[%d].fog.Set(**state)\n\n"%(
                viewerName, i))
        lines.append(indent+"## End Cameras\n\n")

        lines.append(indent+"## Clipping planes\n")
        for i,c in enumerate(self.clipP):
            if not c._modified: continue
            lines.append(indent+"state = %s\n"%str(c.getState()))
            lines.append(indent+"%s.clipP[%d].Set(**state)\n\n"%(
                viewerName, i))
        lines.append(indent+"## End Clipping planes\n\n")

        lines.append(indent+"## Root object\n")
        if self.rootObject._modified:
            lState = self.rootObject.getState()
            lState.pop('protected')
            lines.append(indent+"state = %s\n"%str(lState))
            lines.append(indent+"%s.rootObject.Set(**state)\n\n"%(
                viewerName))
        lines.append(indent+"## End Root Object\n\n")

        # get material definition of root object
        lines.extend(self.rootObject.getGeomMaterialCode(
            '%s.rootObject'%viewerName, indent))

        # rootObject's clipping planes
        lines.extend(self.rootObject.getGeomClipPlanesCode(
            '%s.rootObject'%viewerName, indent))
        
        return lines
    

    def getObjectsStateDefinitionCode(self, viewerName, indent='', withMode=True, includeBinding=True):
        """Returns a list of strings containing Python code to restore the
state of the geometry objects in the viewer.
indent is the level of indentation that is desired.
viewerName is the name of the Viewer object in the source code generated.
"""
        lines = []
        
        if withMode is True:
            lines.append(indent+"if mode=='objects' or mode=='both':\n")
            indent += '    '
            lines.append(indent+"pass #needed in case there no modif\n")

        lines.append(indent+"##\n")
        lines.append(indent+"## Saving State for objects in Viewer\n")
        lines.append(indent+"##\n\n")
        for obj in self.rootObject.AllObjects():
            if obj == self.rootObject:
                continue
            lines.append(indent+"## Object %s\n"%obj.fullName)
            if obj._modified:
                lines.append(indent+"state = %s\n"%str(obj.getState()))
                lines.append(indent+"obj = %s.FindObjectByName('%s')\n"%(
                    viewerName, obj.fullName))
                lines.append(indent+"if obj:\n")
                lines.append(indent+"    obj.Set(**state)\n\n")
                
                if isinstance(obj, Geom):
                    # get material definition of obj
                    lines.extend(obj.getGeomMaterialCode('obj', indent, includeBinding=includeBinding))

                    # clipping planes for obj
                    lines.extend(obj.getGeomClipPlanesCode('obj', indent))
        lines.append(indent+"## End Object %s\n\n"%obj.fullName)

        return lines


    def InitLighting(self, camera):
	"""Set up some basic lighting (single infinite light source)
"""
        #print "InitLighting"
        camera.Activate()
        self.lightModel = LightModel(self, twoSide=True)

        # Graham light model
        # 'key light': front, top left (80%)
        self.lights[0].Set( enabled=viewerConst.YES,
                            direction=(-0.7, 0.8, 1.4, 1.),
                            diffuse = [0.75, 0.75, 0.75, 1.],
                            specular=(.8, .8, .65, 1.),
                            ambient=(.12, .12, .12, 1.),
                            tagModified=False)

        # 'fill light': front, top right, warmer (15-30%)
        self.lights[1].Set( enabled=viewerConst.YES,
                            direction=(9., .5, 3., 1.),
                            diffuse=(.35,.25,.175,1.),
                            specular=(.0,.0,.0,1.),
                            ambient=(.12, .12, .12, 1.),
                            tagModified=False )

        # 'reflective or back': back, bottom center, cooler (10-30%)
        self.lights[2].Set( enabled=viewerConst.YES,
                            direction=(0., -.45, -.35, 1.),
                            diffuse=(.4,.5,.75,1.),
                            specular=(.0,.0,.0,1.),
                            ambient=(.12, .12, .12, 1.),
                            tagModified=False )

        self.currentLight = self.lights[0]
        self.currentLight.hasBeenCurrent = 1
        #for l in self.lights:
            # activate is needed else lights other than 0 seem to tbe lost
        #    self.cameras[0].Activate()
        #    l.apply()


    def findGeomsByName(self, pat, parent=None):
        """find geoms whose name matches the string pat,
a list of matches is returned, option parent: if parent,
search only children of parent.
see also: FindObjectByName"""

        # here 'root|mol|MSMS-MOL' is interepreter as either of the strings
        # (in FindObjectByName | is interpreted as a parent child separator)
        if parent is None:
            allObjs = self.rootObject.AllObjects()
        else:
            allObjs = parent.AllObjects()
        match = re.compile(pat)
        result = {}
        for object in allObjs:
            if match.search(object.name):
                result[object] = 1
            for c in object.children:
                if match.search(c.name):
                    result[c] = 1
        return result.keys()


    def __init__(self, #master=None,
                 nogui=0, #screenName=None,
                 guiMaster=None, #classCamera=None,
                 autoRedraw=True,
                 #verbose=True,
                 eventHandler=None,
                 cnf={}, **kw):
	"""Viewer's constructor
"""
        if eventHandler:
            assert isinstance(eventHandler, EventHandler)
            self.eventHandler = eventHandler
        else:
            self.eventHandler = EventHandler()

        #print "Viewer.__init__" 
        self.isInitialized=0
        self.suspendRedraw = True # set to True to Prevent redrawing

        self.ownMaster = False  # is set to True if this Viewer creates its master window
        self.redoDspLst = 0       # set to 1 to for Display List re-calculation
        
        if sys.platform == 'osf1V4':
            # for some reason using nested display lists on DecAlpha OpenGL
            # freezes the viewer
            self.useMasterDpyList = 0
        else:
            self.useMasterDpyList = 1
        self.oldUseMasterDpyList = self.useMasterDpyList
        self.singleDpyList = False # when false, any object can have a dpyList
                                   # which is called from the master
                                   # when True, only the master list is built
        self.uniqID = '_'+str(time())

        # timing
        self.numberOfRedraws = 0
        self.timeButonDown = 0
        self.timeButonUp = 0
        self.lastFrameTime = 0
        self.avgFrameRate = 0.  # average frame rate between button down and
                                # button up. Compute in tackball
        self.autoRedraw = autoRedraw

        # identifier returned by last called to after_idle in autoRedraw mode
        # it is set to None of no autoredraw has been posted
        self.pendingAutoRedrawID = None
        
        self.rootObject = None
        self.currentObject = self.rootObject
	self.lastObjectNumber = 1
	self.lastPick = None
	self.lastPickedCamera = None
        self._ButtonPress = False # set to True in Trackball.MouseButtonUp
                                  # and False in Trackball.MouseButtonDown
	self.currentLight = None
	self.lights = []
	self.GUI = None
        self.redirectTransformToRoot = 0

        self.Reset_cb = self.ResetCurrentObject
        self.Normalize_cb = self.NormalizeCurrentObject
        self.Center_cb = self.CenterCurrentObject
        self.Delete_cb = self.DeleteCurrentObject

        self.redrawLock = thread.allocate_lock()
        self.objectsNeedingRedo = {}
        
        self.materialEditor = None

	self.cameras = [ ]

        self.maxClipP = 0
	self.clipP = []
        self.activeClippingPlanes = 0 # number of active clipping planes
        self.cappedGeoms = [] # will hold a list of (geom, clipp) pairs

        self.lights = []
        self.OverAllLightingIsOn = True
        self.currentLight = None

        self.activeScissor = 0 # number of active scissors

        # lists of object not using a dispalylist (managed in ReallyRedraw)
        self.noDpyListOpaque = []
        self.noDpyListTransp = []

	# create the root objects
        root = Geom('root', protected=True, shape=(0,0),
                    inheritPointWidth=0,
                    inheritLighting=0,
                    inheritLineWidth=0,
                    inheritStippleLines=0,
                    inheritStipplePolygons=0,
                    inheritFrontPolyMode=0,
                    inheritBackPolyMode=0,
                    inheritShading=0,
                    inheritSharpColorBoundaries=0,
                    inheritCulling=0,
                    inheritXform=0)
        self.rootObject = root
        self.rootObject.viewer = self
	root.materials[GL.GL_FRONT].InitDefault()
	root.materials[GL.GL_BACK].InitDefault()
	root.inheritMaterial = viewerConst.NO
	root.inheritFrontPolyMode = viewerConst.NO
	root.inheritBackPolyMode = viewerConst.NO
        root.inheritShading = viewerConst.NO
        root.inheritLighting = viewerConst.NO
        root.inheritSharpColorBoundaries = viewerConst.NO
	root.inheritCulling = viewerConst.NO
        self.rootObject = root
        self.currentObject = root
        root.hasBeenCurrent = 1
        self.showPickedVertex = 1
        self.rootObject._modified = False

        # find out if PolygonOffset is supported
        self.hasOffsetExt = True
        self.polyOffset = None   # OpenGL function to call for PolygonOffset
        try:
            self.polyOffset = GL.glPolygonOffset
        except AttributeError:
            try:
                self.polyOffset = GL.glPolygonOffsetEXT
            except AttributeError:
                self.hasOffsetExt = False

        # flags to track modifier positions
        self.kbdModifier = { 'Shift_L':0,  'Shift_R':0,
                             'Control_L':0, 'Control_R':0,
                             'Alt_L':0, 'Alt_R':0,
                             'Meta_L':0, 'Meta_R':0 }

        
        # list of functions to be called when picking occurs
        self.pickingFunctions = []
        
        self.AddPickingCallback(self.unsolicitedPick)
        self.pickLevel = 'vertices' #can also be 'parts'

        self.isInitialized=1
        self.needsRedraw = True
        self.hasRedrawn = False

        # tiled rendering
        self.tileRender = False
        self.tileRenderCtx = None
        self.tileRenderFinalWidth = 0
        self.tileRenderFinalHeight = 0
        self.tileRenderBorder = 0
        self.tileRenderOutput = None
        self.tileRenderBuffer = None
        self.tileCount = None

	#c = self.AddCamera(master=master, #screenName=screenName,
        #                                  #classCamera=classCamera,
        #                   stereo='none', num=0,
        #                   cnf=cnf, **kw)


    def update(self):
        self.master.update()

        
    def saveViewerState(self, filename):
        """Save the source code to restore the state of all the Viewer, i.e.
Cameras, Fog, LightModel, Lights, Clipping Planes and rootObject
"""
        code = self.getViewerStateDefinitionCode('vi')
        f = open(filename, 'w')
        dum = map(lambda x, f=f: f.write("%s"%x), code)
        f.close()


    def saveObjectsStates(self, filename):
        """Save the source code to restore the state of all objects in the
viewer except for rootObject.
"""
        code = self.getObjectsStateDefinitionCode('vi')
        f = open(filename, 'w')
        dum = map(lambda x, f=f: f.write("%s"%x), code)
        f.close()

            
    def saveViewerAndObjectsStates(self, filename):
        """Save the source code to restore the state of all objects in the
Viewer's"""
        f = open(filename, 'w')
        code = self.getViewerStateDefinitionCode('vi')
        dum = map(lambda x, f=f: f.write("%s"%x), code)
        code = self.getObjectsStateDefinitionCode('vi')
        dum = map(lambda x, f=f: f.write("%s"%x), code)
        f.close()


    def startTileRendering(self, border=0, outputFile='test.jpg',
                           width=None, height=None, backBuffer=True,
                           checkeredBackground=False):
        self.stopAutoRedraw()

        lDepthArray = self.currentCamera.GrabZBufferAsArray()
        lDepthArray = lDepthArray.ravel()
        zmin = min(lDepthArray)
        zmax = max(lDepthArray)

        self.tileRender = True
        self.tileRenderCtx = TRctx(zmin, zmax)
        self.tileRenderCtx.checkeredBackground = checkeredBackground
        # set to true to see tiles rendered in foreground
        self.tileRenderCtx.backBuffer = backBuffer
        # set to true to create colors tiles background from black to white
        self.tileRenderOutput = outputFile
        if width is not None:
            self.tileRenderFinalWidth = width
        #    self.tileRenderFinalHeight=None
        if height is not None:
            self.tileRenderFinalHeight = height
        #    self.tileRenderFinalWidth = None
        if height is None and width is None:
            self.tileRender = False

            
    def stopTileRendering(self):
        self.tileRender = False
        self.startAutoRedraw()
        

    def unsolicitedPick(self, pick):
        """treat and unsollicited picking event"""
##         if self.isShift():
##             self.centerCurrentObject(pick)
##         elif self.isControl():
##             self.selectGeom(pick)
##         else:
##             for object, items in pick.hits.items():
##                 print "Object", object.name, pick.type, items
        for object, items in pick.hits.items():
            print "unsolicitedPick Object", object.name, pick.type, items
                
            
    def isShift(self):
        return self.kbdModifier['Shift_L'] or self.kbdModifier['Shift_R']


    def isControl(self):
        return self.kbdModifier['Control_L'] or self.kbdModifier['Control_R']


    def isAlt(self):
        return self.kbdModifier['Alt_L'] or self.kbdModifier['Alt_R']


    def getModifier(self):
        if self.kbdModifier['Shift_L']: return 'Shift_L'
        elif self.kbdModifier['Control_L']: return 'Control_L'
        elif self.kbdModifier['Alt_L']: return 'Alt_L'
        elif self.kbdModifier['Shift_R']: return 'Shift_R'
        elif self.kbdModifier['Control_R']: return 'Control_R'
        elif self.kbdModifier['Alt_R']: return 'Alt_R'
        else: return 'None'

        
    def modifierDown(self, event):
        """track changes in SHIFT, CONTROL, ALT kyes positions"""
        #print "event.keysym", event.keysym

        if event.keysym in ['Shift_L', 'Shift_R', 'Control_L', 'Control_R',
                            'Alt_L', 'Alt_R']:

            self.kbdModifier[event.keysym] = 1

            # grab all event to make sure get the key release event even
            # if the mouse is outside the application
            #self.currentCamera.master.grab_set_global()
            #self.master.grab_set_global()
        

    def modifierUp(self, event):
        """track changes in SHIFT, CONTROL, ALT keys positions"""

        #print 'modifier up', event.keysym
        if event.keysym in ['Shift_L', 'Shift_R', 'Control_L', 'Control_R',
                            'Alt_L', 'Alt_R']:
            self.currentCamera.configure(cursor=cursorsDict['default'])
            self.kbdModifier[event.keysym] = 0
            # release the grab 
            self.currentCamera.master.grab_release()
            #self.master.grab_release()

        
##     def __del__(self):
## 	"""Viewer's destructor"""

## 	self.GUI.root.destroy()
## 	for c in self.cameras:
## 	    c.frame.master.destroy()
    def Exit(self):
        #if sys.platform != "win32":
            # FIX this: the following line causes
            # "PyEval_Restore_Thread: Null tstate" error on Windows
            # when Python interpreter quits.
            # Python2.6.5 , numpy 1.4.1, Togl 1.7
            #if self.GUI:
                #self.GUI.occlusionCamera.destroy()
    
        self.materialEditor._destroy()
        if self.GUI:
            # When "occlusionCamera" class was commented out
            # "PyEval_Restore_Thread: Null tstate" error on Windows returned 
            # ( happens when Python interpreter quits).
            # The following line seems to fix this problem.
            self.GUI.GraphToolpanel.destroy()
            #
            self.GUI.root.destroy()
            
        if self.autoRedraw:
            self.stopAutoRedraw()
            sleep(0.1)
            self.master.update()
        if len(self.cameras) > 1:
            lReversedListOfCameras = self.cameras[:]
            lReversedListOfCameras.reverse()
        else:
            lReversedListOfCameras = self.cameras[:]
        for c in lReversedListOfCameras:
            self._DeleteCamera(c)
        if self.ownMaster:
            self.master.destroy()

    def processPicking(self, pick):
        """call all functions that need to be called after something has been
        picked"""
            
        for func in self.pickingFunctions:
            func(pick)
            

    def AddPickingCallback(self, func):
	"""Add function func to the list of picking callback fucntions"""
        if func not in self.pickingFunctions:
            self.pickingFunctions.append(func)


    def RemovePickingCallback(self, func, modifier=None):
	"""Delete function func from the list of picking callbacks"""
        if type(func)==types.StringType:
            findFunc = self.currentCamera.eventManager.FindFunctionByName
	    func = findFunc(func, self.pickingFunctions)
        if not func:
            return

	self.pickingFunctions.remove(func)


    def SetPickingCallback(self, func):
	"""Replace the functions to be called after a picking event.
        Return the list of functions that got replaced"""

        save = self.pickingFunctions
        # used to restore save callbacks sent back earlier
        if type(func)==types.ListType:
            self.pickingFunctions = func
        else:
            self.pickingFunctions = [func]
        return save


##     def selectGeom(self, pick):
## 	"""Picking call back"""

##         l = len(pick.hits.keys())
##         if l == 0:
##             if pick.camera.currentTransfMode == 'Object':
##                 self.BindTrackballToObject(self.rootObject)
##             else:
##                 self.BindTrackballToObject(self.currentObject)
##             return
            
##         if l > 1:
##             print "WARNING selectPick does not work for multiple objects"""
##             return

##         object = pick.hits.keys()[0]
## 	if isinstance(object, Geom):
##             self.BindTrackballToObject(object)
## 	elif isinstance(object, ClippingPlane):
## 	    self.BindTrackballToClip( object )

## 	elif isinstance(object, Light):
## 	    self.BindTrackballToLight( object )


##     def centerCurrentObject(self, pick):
## 	"""Center the world on the picked vertex"""

##         if pick.type=='parts':
##             self.vf.message("WARNING: centering for picking done on parts is not implmented.")
##             return
        
##         g = numpy.zeros( (3), 'f' )
##         object = self.rootObject
##         inv = object.GetMatrixInverse()
##         for obj, vertInd in pick.hits.items():
##             if len(vertInd)==0:
##                 print 'WARNING: object',o.name,' is not vertex pickable'
##             else:
##                 m = numpy.dot( inv, obj.GetMatrix() )
##                 vert = obj.vertexSet.vertices * m
##                 vertsel = numpy.take( vert, vertInd )
##                 g = g + numpy.sum(vertsel)/len(vertsel)
##         object.SetPivot( g )


    def toggleTransformRootOnly(self, event=None):
        self.TransformRootOnly(not self.redirectTransformToRoot)

        
    def TransformRootOnly(self, value):
        """if value is true, all 3D transformations are redirected to the
        root object, and all objects added subsequently are have their
        transformation redirected automatically"""
        root = self.rootObject
        objs = root.AllObjects()
        objs.remove(root)
        objsTransformableOnly = []
        for o in objs:
            if isinstance(o, Transformable) is True:
                objsTransformableOnly.append(o)
        if value:
            for o in objsTransformableOnly:
                o.RedirectTransformTo(root)
            self.redirectTransformToRoot = 1
        else:
            for o in objsTransformableOnly:
                o.RedirectTransformTo(None)
            self.redirectTransformToRoot = 0
            # since the root transormation might have changed we need to
            # recompute R and Ri for the current Object
            if self.currentObject is not self.rootObject:
                if isinstance(self.currentObject, Transformable):
                    self.currentObject.FrameTransform(self.currentCamera)

        if self.GUI:
            self.GUI.moveRootOnly.set(value)


    def ResetCurrentObject(self):
	"""Reset the current's object transformation
"""
	#self.trackball.autospin = 0
        if self.redirectTransformToRoot:
            self.rootObject.ResetTransformation()
        else:
            self.currentObject.ResetTransformation()
	self.Redraw()


    def Reset_cb_arg(self, event):
	"""Function to be bound to the GUI's reset widget"""
	self.Reset_cb()


    def AutoDepthCue(self, nearOffset=0.0, farOffset=0.0):
        """
"""
        self.currentCamera.AutoDepthCue()
        self.Redraw()


    def NormalizeCurrentObject(self):
        """Normalize the current object or RootObject
"""
        ##
        ## MS Jan 31 2013 I was unable to get the normalize center to work properly on
        ## object when the current object was not the root. The workaround is to
        ## set the root to be the current onject and hide all but current objects to perform
        ## normalize then restore current object
        lObject = self.rootObject
        if self.redirectTransformToRoot or self.currentObject==self.rootObject:
            hidden = None
        else:
            hidden = []
            for g in self.rootObject.VisibleObjects():
                #print g.name , g.visible
                if g is not self.currentObject and g.visible and not self.AreYouMyDescendant(self.currentObject, g):
                    g.visible = False
                    hidden.append(g)
            current = self.currentObject
            self.SetCurrentObject(self.rootObject)
            #print "CURRENT ROOT", self.currentObject
            #print self.rootObject.VisibleObjects()
                
        if isinstance(lObject, Transformable):
            self.rootObject.Set(scale=(1.,1.,1.))
            mini, maxi = lObject.ComputeBB(self.currentCamera)
            g = numpy.add.reduce( (mini, maxi) ) * .5 # center of BoundingBox
            if lObject == self.rootObject:
                self.rootObject.ConcatTranslation( -g[:3] ) #* lObject.Si )
            else:
                self.rootObject.SetTranslation( -g[:3] ) #* lObject.Si )
            lBox = maxi - mini
            lHalfObject = max(lBox)/2.
            if lHalfObject == 0.:
                lHalfObject = 1.
            c = self.currentCamera
            c.Set(fov=c.fovyNeutral)
            lNewDist = lHalfObject / math.tan(c.fovy/2*math.pi/180.0)
            c.lookFrom = numpy.array( ( 0., 0., c.nearDefault+lNewDist+2*lHalfObject ) )
            c.lookAt = numpy.array((0.,0.,0.))
            c.direction = c.lookAt - c.lookFrom
            c.AutoDepthCue(object=lObject)
        else: # Insert2d
            self.lObject.ResetSize()

        if hidden:
            self.SetCurrentObject(current)
            for g in hidden:
                g.visible = True
        self.Redraw()


    def Normalize_cb_arg(self, event):
        """Function to be bound to the GUI's normalize widget
"""
        self.Normalize_cb()


    def CenterCurrentObject(self):
        """Center the current object on the center of gravity of its subtree
in "redirect transform to root" mode, root is recentered on the selected geometry
"""
        #print "self.currentObject", self.currentObject
        if isinstance(self.currentObject, Transformable):
            #mini, maxi = self.currentObject.BoundingBoxSubtree()
            mini, maxi = self.currentObject.ComputeBB(self.currentCamera)
            #g = numpy.add.reduce( (mini, maxi) ) / 2.0 # center of BoundingBox
            g = mini + ( (maxi-mini) *0.5 ) # center of BoundingBox
            g = numpy.concatenate( (g, [1.0]) )
            
            # FIXME won't work with instance matrices
            m = self.currentObject.GetMatrixInverse(root=self.currentObject)
            g = numpy.dot(m, g)[:3]
            if self.redirectTransformToRoot and self.currentObject != self.rootObject:
                self.rootObject.SetPivot( self.currentObject.translation + g )
            else:
                self.currentObject.SetPivot( g )
        else: # Insert2d
            self.currentObject.ResetPosition()
        #self.Redraw()


    def get3DPointFromPick(self,event):
        """Compute 3D point from pick event.

        point, background <- get3DPointFromPick(event)

        Returns a 3D point and a boolean which is True if picking occured on
        background
        """
        if isinstance(self.currentObject, Transformable):
            winx = event.x
            winy = self.currentCamera._height - event.y
            g, background = self.currentCamera.xyzFromPixel(winx, winy)
            g = numpy.concatenate( (g, [1.0]) )
            # FIXME won't work with instance matrices
            if self.redirectTransformToRoot or \
                   self.currentObject==self.rootObject:
                m = self.rootObject.GetMatrixInverse(root=self.rootObject)
                g = numpy.dot(m, g)[:3]
            else:
                m = self.currentObject.GetMatrixInverse(
                    root=self.rootObject)
                g = numpy.dot(m, g)[:3]
            return g, background
                

    def pivotOnPixel(self, event):
        #print "pivotOnPixel"
        g, background = self.get3DPointFromPick(event)
        #print 'pivotOnPixel setting pivot', background, g
        if not background: # not background pick
            if self.redirectTransformToRoot or \
                   self.currentObject == self.rootObject:
                oldPivot = self.rootObject.SetPivot( g )
            else:
                oldPivot = self.currentObject.SetPivot( g )
            return oldPivot
        return None

    def Center_cb_arg(self, event):
	"""Function to be bound to the GUI's center widget"""
	self.Center_cb()

    def _FocusCameraBB(self):
        cam = self.currentCamera
        dtime = time() - self.focusTimeStart
        steps = self.steps
        i = (dtime/self.focusDuration)*steps
        cam.Set(fov=cam.fovy + self.diff_fovy/float(steps))
        self.rootObject.ConcatTranslation( self.diffVect[:3]/float(steps))
        lNewDist = self.lHalfObject / math.tan(cam.fovy/2*math.pi/180.0)
        newDist = cam.nearDefault+lNewDist+self.lHalfObject
        dist = self.oldLookFrom*(steps-i-1.)/steps + newDist*(i+1.)/steps
        cam.lookFrom = numpy.array( ( 0., 0., dist ) )            
        cam.direction = cam.lookAt - cam.lookFrom 
        self.OneRedraw()
        if i<steps:
            self.master.after(5, self._FocusCameraBB)
        elif self.after_cb:
            self.after_cb()
            del self.steps
            del self.diff_fovy
            del self.lHalfObject
            del self.diffVect
        
    def FocusCameraBBT(self, mini, maxi, focusDuration=.5, after_cb=None):
        """Moves camera to provided as [mini, maxi] bounding box"""
        g = numpy.sum( (mini, maxi),0 ) * .5 # center of BoundingBox        
        self.rootObject.Set(scale=(1.,1.,1.))
        lBox = maxi - mini
        self.lHalfObject = max(lBox)/2.
        if self.lHalfObject == 0.:
            self.lHalfObject = 1.
        cam = self.currentCamera
        self.oldLookFrom = cam.lookFrom[2]
        if self.lastFrameTime == 0:
            self.OneRedraw()
        self.steps = int(2/self.lastFrameTime)

        self.diff_fovy =  cam.fovyNeutral - cam.fovy
        cam.lookAt = numpy.array((0.,0.,0.))
        Rmini, Rmaxi = self.rootObject.ComputeBB(self.currentCamera)
        Rg = numpy.sum( (Rmini, Rmaxi), 0) * .5

        self.diffVect = -g-self.rootObject.translation
        if not (g-Rg).any():
            self.diffVect = -Rg
        self.focusTimeStart = time()
        self.focusDuration = focusDuration
        self.after_cb=after_cb
        self._FocusCameraBB()
        

    def FocusOnBoxT(self, obj, mini, maxi):
        """Focuses currentCamera on a box defined by the mini and maxi points
        by transforming object"""
        cb = CallbackFunction(self.afterFocus, obj, mini, maxi)
        self.FocusCameraBBT(mini, maxi, after_cb=cb)


    def afterFocus(self, obj, mini, maxi):
        currentObject = self.currentObject
        self.currentObject = obj
        self.CenterCurrentObject()
        self.currentObject = currentObject 
        self.currentCamera.AutoDepthCue(object=obj)
        self.OneRedraw()
        #self.update()
        g = ( (0.5*(mini[0]+maxi[0])), (0.5*(mini[1]+maxi[1])),
              (0.5*(mini[2]+maxi[2])) )
        print 'after focus set pivot ', g
        obj.SetPivot(g)


    def FocusCameraBB(self, obj, mini, maxi):
        """Moves camera to provided as [mini, maxi] bounding box"""

        #print 'mini, maxi', mini, maxi
        g = numpy.sum( (mini, maxi), 0 ) * .5 # center of BoundingBox

        # reset root scaling
        self.rootObject.Set(scale=(1.,1.,1.))

        # compute half of largest length
        lBox = maxi - mini
        lHalfObject = max(lBox)/2.
        if lHalfObject == 0.:
            lHalfObject = 1.

        #print 'lHalfObject', lHalfObject
        cam = self.currentCamera

        # save looFrom
        oldLookFrom = cam.lookFrom[2]
        #print 'oldLookFrom', oldLookFrom
        
        # find out how many steps we want to use
        if self.lastFrameTime == 0:
            self.OneRedraw()      
        steps = int(2/self.lastFrameTime)
        steps = 1
        
        # compute delta fovy to go to neutral
        diff_fovy =  cam.fovyNeutral - cam.fovy
        #print 'diff_fovy', diff_fovy
        # overwerite camera's look at
        cam.lookAt = numpy.array((0.,0.,0.))

        # compute bounding box of visible scene
        Rmini, Rmaxi = self.rootObject.ComputeBB(self.currentCamera)

        #compute
        Rg = numpy.sum( (Rmini, Rmaxi), 0) *.5
        #print 'Rg', Rg

        # compute difference between center of original BB and desired BB
        diffVect = -g-self.rootObject.translation
        #print 'diffVect', diffVect
        #print 'g-Rg', g-Rg
        if not (g-Rg).any():
            diffVect = -Rg
            #print 'diffVect 1', diffVect

        delta = diffVect[:3]/float(steps)
        toRad = math.pi/180.
        onear = cam.near
        ofar = cam.far
        #print 'onear, ofar', onear, ofar
        for i in range(steps):
            #print i, cam.near, cam.far, delta
            cam.Set(fov=cam.fovy + diff_fovy/float(steps))
            self.rootObject.ConcatTranslation( delta )
            lNewDist = lHalfObject / math.tan( cam.fovy*.5*toRad )
            newDist = cam.nearDefault+lNewDist+lHalfObject
            dist = oldLookFrom*(steps-i-1.)/steps + newDist*(i+1.)/steps
            #print dist, newDist*(i+1.)/steps, oldLookFrom*(steps-i-1.)/steps
            cam.lookFrom = numpy.array( ( 0., 0., dist ) )
            cam.direction = cam.lookAt - cam.lookFrom 
            dc = dist-oldLookFrom
            cam.AutoDepthCue(object=obj)
            #print 'dc', dc
            #cam.Set(near=onear+dc+delta[2], far=ofar+dc+delta[2])
            self.OneRedraw()
            #self.update()
                    
    
    def FocusOnBox(self, obj, mini, maxi, vector=None):
        """Focuses currentCamera on a box defined by the mini and maxi points
        by transforming object"""
        root = self.rootObject
        # save root transform
        cam = self.currentCamera
        #print 'near, far BEFORE', cam.near, cam.far, cam.fog.start, cam.fog.end
        orig = [cam.fovy, cam.lookFrom.copy(), cam.near, cam.far,
                cam.fog.start, cam.fog.end,
                root.rotation[:], root.translation[:],
                root.scale[:], root.pivot[:]]
        
        # R N on obj
        self.suspendRedraw = True
        self.rootObject.ResetTransformation()
        self.FocusCameraBB(obj, mini, maxi)

        # does not work. When the rotation is added the object ends up
        # not centerd
        # rotate so that vector is aligned with -Z (i.e. pointing towards
        # camera
        if vector:
            g = ( (0.5*(mini[0]+maxi[0])), (0.5*(mini[1]+maxi[1])),
                  (0.5*(mini[2]+maxi[2])) )
            obj.SetPivot(g)
            rot = numpy.array(rotVectToVect( vector, (0,0,-1) ), 'f')
            root.ConcatRotation(rot.flatten())

        #bbcenter = (maxi-mini)*0.5
        #sceneCenter = root.ComputeBB()
        #v = 
        # get MAA
        from DejaVu2.scenarioInterface.animations import OrientationMAA
        maa = OrientationMAA(root, 'focus', None, 0, keyPos=[0,0,0,20])

        # restore rot transformation
        cam.Set(fov=orig[0], lookFrom=orig[1], near=orig[2], far=orig[3])
        cam.fog.Set(start=orig[4], end=orig[5])
        root.rotation = orig[6]
        root.translation = orig[7]
        root.scale = orig[8]
        root.pivot = orig[9]
        self.suspendRedraw = False
        self.OneRedraw()
        # mow play maa
        try:
            import Scenario2
            maa.runIn(1., self.master.after)
        except ImportError:
            maa.run()
            
##         from DejaVu2.scenarioInterface import getActor
##         from Scenario2.keyframes import KFAutoCurrentValue, KF
##         actor = getActor(root, 'transformation')
##         val = [matToQuaternion(rotation), translation, scale, pivot]
##         kf0 = KFAutoCurrentValue( 0 , actor)
##         kf1 = KF( 30, val)
##         i = Interval( kf0, kf1, generator=actor.behaviorList)
##         actor.addIntervals( [i] )

        
##         currentObject = self.currentObject
##         self.currentObject = obj
##         self.CenterCurrentObject()
##         self.currentObject = currentObject 
##         #print 'BEFORE', self.currentCamera.near, self.currentCamera.far
##         #self.currentCamera.AutoDepthCue(object=obj)
##         #print 'AFTER', self.currentCamera.near, self.currentCamera.far
##         self.OneRedraw()
##         #self.update()
        g = ( (0.5*(mini[0]+maxi[0])), (0.5*(mini[1]+maxi[1])),
              (0.5*(mini[2]+maxi[2])) )
        obj.SetPivot(g)

    def DeleteCurrentObject(self):
        """Delete the current object
"""
        #print "DeleteCurrentObject", self

        if self.currentObject.protected is True:
            return

        self.RemoveObject(self.currentObject)
        self.SetCurrentObject(self.rootObject)
        self.Redraw()


    def ToggleDepth(self):
        """Toggles depthcueing
"""
        #print "ToggleDepth"
        self.currentCamera.ToggleDepth()
        self.Redraw()


    def Depth_cb_arg(self, event):
	"""Toggle depthcueing flag
"""
        #print "Depth_cb_arg"
        self.ToggleDepth()


    def SSAO_cb_arg(self, event=None):
        camera = self.currentCamera
        if camera.ssao is False: # we are switching to SSAO
            camera.notSsaoNear = camera.near # save camera.near
            camera.Set(near=2.0)
        else:
            camera.Set(near=camera.notSsaoNear)
            
        self.currentCamera.ssao = not self.currentCamera.ssao
        if self.GUI is not None:
            self.GUI.ssaoTk.set(self.currentCamera.ssao)
        self.Redraw()

    def CurrentCameraBackgroundColor(self, col):
	"""Method to be bound to a GUI widget"""

	if col:
#	    col = OneColor(col)
	    self.currentCamera.Set(color = col)
	    self.currentCamera.fog.Set( color= col )
	    self.currentCamera.Redraw()


    def FindObjectByName(self, fullName):
        """find a geometry using a name
returns the corresponding unique object
see also: findGeomsByName 
"""
        names = string.split(fullName,'|')
        obj = self.rootObject
        for n in names[1:]:
            no = None # find the child with name n
            for o in obj.children:
                if o.name==n:
                    no = o
                    break
            if no is None: # failed
               # print "WARNING: object %s not found in %s\n"%(n,fullName)
                return None
            obj = no # iterate with child named n
        return obj


    def ensureUniqueName(self, aName, parent=None, local=True, exclude=None):
        """return a localy unique name if "local" is True
return a globaly unique name if "local" is False
default is local
exclude can handle one geometry wich will be excluded from the name comparison
"""
        if parent is None or local is False:
            parent = self.rootObject

        if local is True:
            lChildrenList = parent.children
        else:
            lChildrenList = parent.AllObjects()
            
        for lChild in lChildrenList :
            if aName == lChild.name :
                if exclude is None or exclude != lChild:
                    break
        else: #if we didn't break
           return aName 

        lNamebase = aName
        for lChild in lChildrenList :
            if aName == lChild.name :
                lNamebase = aName   
                if lNamebase == None :
                    namebase = ""
                # to create a nice name base, we get rid of the last digits                
                while lNamebase and lNamebase[-1].isdigit() :
                    lNamebase = lNamebase[:-1]
                break
        #print "lNamebase" , lNamebase

        # maybe lNamebase is good enough 
        for lChild in lChildrenList :                                            
            if lNamebase == lChild.name :
                if exclude is None or exclude != lChild:
                    break # get out of loop for lChild
        else: # if we didn't break on the loop for lChild
            return lNamebase 
                 
        # let's complete the name now
        for i in range(len(lChildrenList)) :
            lName = lNamebase + str(i+1)
            for lChild in lChildrenList :                                            
                if lName == lChild.name :
                    if exclude is None or exclude != lChild:
                        break # get out of loop for lChild
            else: # if we didn't break on the loop for lChild
                break # get out of loop for i
            
        #print "lName" , lName
        return lName


    def AddObjectAndParents(self, obj, parent=None, 
                            redraw=True, redo=True, local=True):
        """Same as AddObject but add parents that are not yet in the viewer
return the highest parent added.
local specifies if the names will be locally or globaly unique.
"""
        assert obj.viewer is None
               
        # find last parents before a geom belonging to this viewer
        lListOfAdditions = []
        p = obj
        while p:
            lListOfAdditions.append(p)
            if p.parent is None:
                break
            if p.parent.viewer==self:
                break
            elif p.viewer is not None:
                raise RuntimeError('geom %s already belong to a viewer'%p.name)
            p = p.parent

        # add the object and its parents to the viewer
        for geom in lListOfAdditions[::-1]: # going backwards
            if geom.viewer is None:
                if geom.replace is False:
                    geom.name = \
                        self.ensureUniqueName(geom.name, geom.parent, local)
                else:
                    if geom.parent is None:
                        lparent = self.rootObject
                    else:
                        lparent = geom.parent
                    for c in lparent.children:
                        if geom.name == c.name:
                            if hasattr(c,'node'):                               
                                c.node().nodeOnlyGeomRemoval(c)
                            break;

                self.AddObject(geom, parent=geom.parent)
            elif geom.viewer != self:
                raise RuntimeError('geom %s already belong to a viewer'%geom.name)
        
        return p



    def AddObject(self, obj, parent=None, 
                  redraw=True, redo=True):
        """Add an object to the viewer
- parent: by default a geoemtry is added under the root object or under the
          specified parent geometry
- replace: when True, if a sibling has the same name the sibling is replaced
           when false, if a sibling has the same name the new geom is renamed
           Defaults to True
"""
		# obj refers to the geom

        #assert obj.viewer == None , "ERROR: object already has a viewer"

        if parent is None:
            parent = self.rootObject

        assert parent.viewer == self , \
                        "ERROR: the parent must be in the same viewer"

#        if replace is None:
#            replace = obj.replace
#        if replace is None:
#            replace = True
                        
		#we verify that the full path name is unique
        #if ( replace is None ) or ( replace == False ) :
        #    for child in parent.children :
        #        if child.name == obj.name : 
        #            assert False

        self.redrawLock.acquire()
        try:
            self._AddObject(obj, parent, #replace, 
                            redo=redo)
            obj.OnAddGeomToViewer()
            
        finally:
            self.redrawLock.release()
            if redraw:
                self.Redraw()
  

    def _AddObject(self, obj, parent, #replace, 
                   redo=False):
	"""Add an object to the viewer
- parent:  geometry under which obj will be added
- replace: when True, if a sibling has the same name the sibling is replaced
           when false, if a sibling has the same name the new geom is renamed
"""

        #assert obj.viewer == None , "ERROR: object already has a viewer"
        assert parent.viewer == self , \
                        "ERROR: the parent must be in the same viewer"
        
        #we verify that the full path name is unique
        #if ( replace is None ) or ( replace == False ) :
        #    for child in parent :
        #        if child.name == name :
        #            assert False , "ERROR: one of the child has already the same name"

        ### WARNING this is not thread save when redraw redo are True
        ### and obj.dpyList is None

        if not ( isinstance(obj, Geom) or isinstance(obj, Insert2d) ):
            raise AttributeError('Objects has to inherit from Geom')
        if not ( isinstance(parent, Geom) or isinstance(parent, Insert2d) ):
            raise AttributeError('Parent Objects has to inherit from Geom')

        assert obj.replace in [False, True, 'force']
#        self.redrawLock.acquire()

        ## check the name of obj
        ## if it is unique among the children of parent then we just add it
        ## if there already is a sibling with the same name we test 'replace'
        ##   if replace is true:
        ##       obj will replace the node with the same name
        ##   else:
        ##       obj will be renamed with a unique name by adding a _##### 
        ##       string at the end of its name where ##### is the number of
        ##       geoms added to the viewer
        ##
        if not obj.name: # create unique object name
            obj.name = str(self.lastObjectNumber)
        else:
            # look for a child of parent with the same name
            for child in parent.children:
                if obj.name==child.name:
                    if obj.replace is not False: # replace child by obj
                        # copy all potential children to the new node
                        obj.children = child.children
                        child.children = []
                        for c in obj.children:
                            c.parent = obj
                        self.redrawLock.release()
                        # remove the node from the DejaVu2 tree of geometry
                        # and update tree widget in TreeViewer
                        try:
                            if obj.replace == 'force':
                                 child.protected = False
                            self.RemoveObject(child, redo=redo)
                        finally:
                            self.redrawLock.acquire()

                    else: # make a unique name for the obj
                        if obj.name[-6:]=='_':
                            try: # if it had a number already replace it
                                int(obj.name[-5:])
                                obj.name = obj.name[:-6]+'_%05d'% \
                                           self.lastObjectNumber
                            except: # else append the number
                                obj.name = obj.name+ '_%05d'% \
                                           self.lastObjectNumber
                        else:
                            obj.name = obj.name+ '_%05d'% self.lastObjectNumber
                        
                    break

        # create  the object's full name
        obj.fullName = parent.fullName+'|'+obj.name

        # increment object counter (even in case of replacement, i.e.
        # high water mark)
        self.lastObjectNumber = self.lastObjectNumber + 1

        # add the object as a child of parent
        parent.children.append(obj)
        obj.parent = parent
        
        if self.GUI and obj.listed:
            self.GUI.addObject(obj,parent)
        self.eventHandler.dispatchEvent(AddObjectEvent(obj))

#        if obj.immediateRendering:
#            obj.updateParentsForImmediateRendering(1)

        obj.clipSide = [1]*self.maxClipP

        obj.setViewer(self, not obj.dpyList and redo)

        if isinstance(obj, Transformable):
            # DST & MS modification to enable adding geoms that have children 
            #obj.viewer = self
            if self.redirectTransformToRoot:
                root = obj.parent
                while root.parent:
                    root = root.parent
                obj.RedirectTransformTo(root)
            if obj.transparent:
                obj._setTransparent(1)
        # DST & MS modification to enable adding geoms that have children 
##          if not obj.dpyList and redo:
##              self.objectsNeedingRedo.append(obj)
##              #obj.RedoDisplayList()
#        self.redrawLock.release()


    def _RemoveObjectRec(self, obj ):
        """Remove an object from the viewer
"""
        #print "_RemoveObjectRec", obj

        if obj.protected:
            raise RuntimeError('attempt to remove a protected object from the DejaVu2 viewer: %s'%obj)
            return
        if obj == self.currentObject:
            if obj == self.rootObject:
                p = obj
            else:
                p = obj.parent
            self.SetCurrentObject(p)
            # the selected node will be displayed in TreeView during the 
            # SetCurrentObject process
        children = obj.children[:]
        for o in children:
            self._RemoveObjectRec(o)

        if obj.transparent and isinstance(obj, Transformable):
            for c in self.cameras:
                c.delButtonUpCB(obj.sortPoly_cb, silent=True)
                
        if len(obj.children)==0 and obj.parent is not None:
            obj.parent.children.remove(obj)
            if obj.dpyList: # destroy display list
                self.currentCamera.Activate()
                obj.deleteOpenglList()
                if isinstance(obj, Transformable):
                    obj.deleteOpenglTemplateList()
            obj.viewer = None
            if self.GUI:
                self.GUI.deleteObject(obj)
            self.eventHandler.dispatchEvent(RemoveObjectEvent(obj))
            if obj in self.noDpyListOpaque:
                self.noDpyListOpaque.remove(obj)
            elif obj in self.noDpyListTransp:
                self.noDpyListTransp.remove(obj)
            if obj in self.objectsNeedingRedo.keys():
                del self.objectsNeedingRedo[obj]

            for func in self.pickingFunctions[:]:
                if func.im_self == obj:
                    self.pickingFunctions.remove(func)
                    break
        else:
            assert obj == self.rootObject
            obj.ResetTransformation()


    def RemoveObject(self, obj, redo=False):
        """Remove an object from the viewer
"""
        #print "RemoveObject", obj
        self.redrawLock.acquire()
        try:
            self._RemoveObjectRec(obj)
            #obj.parent = None
        finally:
            self.redrawLock.release()
        if redo:
            self.deleteOpenglList()
            self.Redraw()


    def AreYouMyDescendant(self, parent, child):
        #print "Is", child.name , "a descendant of" , parent.name , "?"

        assert parent != child , "ERROR: parent and child are the same"   
        if child == self.rootObject :
            return False       
        for lChild in parent.children :
            if lChild == child :
                #print "yes"
                return True
            elif self.AreYouMyDescendant(lChild, child) is True :
                return True 
        return False


    def ReparentObject(self, obj, newparent, objectRetainsCurrentPosition=False):
	"""Reparent an object
An AssertionError will be raised if ...
"""
        #print "ReparentObject", obj, newparent

        assert newparent.viewer is None or newparent.viewer==self, \
               "Parent is in another viewer"

        assert self.AreYouMyDescendant(obj , newparent) is False , \
            "ERROR: a parent can't be a descendant of it's parent"
        
        assert obj != newparent , "ERROR: trying to reparent with itself"
        assert obj != self.rootObject , "ERROR: trying to reparent the root"
                            
        if obj.parent==newparent:
            return

        # make sure the name of obj is not used by another child of newparent
        for child in newparent.children:
            if obj.name==child.name:
                obj.name = obj.name+'_%05d'%self.lastObjectNumber
                self.lastObjectNumber += 1

        if obj == self.currentObject:
            self.SetCurrentObject(self.rootObject)
            wasCurrent = True
        else:
            wasCurrent = False

        if objectRetainsCurrentPosition in (1, True):
            m = obj.GetMatrix(root=self.rootObject)
            mInv = newparent.GetMatrixInverse(root=self.rootObject)
            lMatrice = numpy.dot(mInv, m)
            obj.SetTransformation(lMatrice, transpose=True)

        oldparent = obj.parent
        # remove this object from the list of children in parent
        obj.parent.children.remove(obj)
        # add it to the lits of children of the new parent
        newparent.children.append(obj)

        # update the parent attribute
        obj.parent = newparent
        ##for anObject in obj.AllObjects():
        ##    if self.GUI:
        ##        self.GUI.occlusionCamera.chooser.remove(anObject.fullName)
        ##        anObject.fullName = anObject.BuildFullName()
        ##        val = (anObject.fullName,'no comment available', anObject)
        ##        self.GUI.occlusionCamera.chooser.add(val)
        
        # update the Tree widget
        if self.GUI:
            self.GUI.deleteObject(obj)
            self.GUI.addObject(obj, newparent)

        # if the object cas the current one set it back
        if wasCurrent:
            self.SetCurrentObject(obj)
        self.eventHandler.dispatchEvent(ReparentObjectEvent(obj, oldparent))
        
        if objectRetainsCurrentPosition in (0, False) \
          and obj.immediateRendering is False:
            if obj not in self.objectsNeedingRedo.keys():
                self.objectsNeedingRedo[obj] = None
            self.deleteOpenglList() #to allow instance matrices of new parent to be shawn
            self.Redraw()
        
## 	if newparent in self.rootObject.AllObjects():
##             if obj.viewer :
##                 self.stopAutoRedraw()
##                 self.RemoveObject(obj)
##                 self.AddObject(obj, parent=newparent)
##                 self.startAutoRedraw()
##             else:
##                 obj.parent = parent

##             #obj.parent.children.remove(obj)
## 	    #self.GUI.deleteObject(obj)
## 	    #newparent.children.append(obj)
## 	    #obj.parent = newparent
##             #obj.fullName = obj.parent.name+'.'+obj.name
## 	    #self.GUI.addObject(obj, newparent)
## 	else:
## 	    raise ValueError( 'ERROR: parent object not found for %s' % \
## 			      obj.name )


    def beforeRedraw(self):
        pass
        
    def afterRedraw(self):
        pass

    def beforeRebuildDpyLists(self):
        pass

    def afterRebuildDpyLists(self):
        pass


    def doRedraw(self):
        self.needsRedraw = 0
        t1 = time()
        self.numberOfRedraws = self.numberOfRedraws + 1
        for c in self.cameras:
            # added by AG 01/12/2006 to suspend a camera to be redraw
            # needed for ARViewer
            if c.suspendRedraw: continue
            if not c.visible: continue
            c.Activate()

            ## added to port to Mountain Lion, else light are not turned on !
            if c.firstRedraw:
                self.lightModel.apply()
                [l.apply() for l in self.lights]
                c.firstRedraw = False
                
            if c.exposeEvent:
                c.storeGeometry()
                c.exposeEvent = False

            if self.tileRender:
                ## SSAO does not work with tile rendering so we turn it off
                ## do the tile rendering, then render the AO mask in a single
                ## tile, scale and blur it and composite it with the large
                ## image
                ssaoWasOn = False # use to restore SSAO mode
                if c.ssao:
                    c.ssao = False # turn SSAO off
                    ssaoWasOn = True # remember it was on

                # do regular tile rendering
                c.lift()
                tr = self.tileRenderCtx
                print 'tile size', c._width, c._height, self.tileRenderBorder
                tr.tileSize(c._width, c._height, self.tileRenderBorder)
                w = self.tileRenderFinalWidth 
                h = self.tileRenderFinalHeight
                if w is not None and h is None:
                    h = int(w*c._height/float(c._width))
                elif w is None and h is not None:
                    w = int(h*c._width/float(c._height))
                self.tileRenderFinalHeight = h

                if c.stereoMode.startswith('SIDE_BY_SIDE'):
                    self.tileRenderFinalWidth = w/2
                    numOfPasses = 2
                    self.tileRenderBuffer = [None, None]
                    im = [None, None]
                else:
                    self.tileRenderFinalWidth = w
                    numOfPasses = 1
                    self.tileRenderBuffer = [None]
                    im = [None]

                for i in range(numOfPasses):
                    if c.stereoMode.endswith('_CROSS'):
                        if i == 0:
                            c.imageRendered = 'RIGHT_EYE'
                        else:
                            c.imageRendered = 'LEFT_EYE'
                    elif c.stereoMode.endswith('_STRAIGHT'):
                        if i == 0:
                            c.imageRendered = 'LEFT_EYE'
                        else:
                            c.imageRendered = 'RIGHT_EYE'

                    print 'New tileRenderBuffer'
                    self.tileRenderBuffer[i] = numpy.zeros( self.tileRenderFinalWidth * h * 3,
                                                           numpy.uint8)
                    tr.oneTileRenderBuffer = numpy.zeros(
                        c._width * c._height * 3, numpy.uint8)

                    print 'image size', self.tileRenderFinalWidth, self.tileRenderFinalHeight
                    tr.imageSize(self.tileRenderFinalWidth,
                                 self.tileRenderFinalHeight)

                    tr.imageBuffer(GL.GL_RGB, GL.GL_UNSIGNED_BYTE,
                                   self.tileRenderBuffer[i])
                    #tr.rowOrder(113) # 112

                    c.Activate()

                    if c.projectionType==c.PERSPECTIVE:
                        self.tileRenderCtx.perspective(c.fovy,
                                                   float(c._width)/float(c._height),
                                                       c.near, c.far)
                    else:
                        self.tileRenderCtx.ortho(c.left, c.right,
                                                 c.bottom, c.top, 
                                                 c.near, c.far)
                    moreTiles = 1
                    self.tileCount = 0
                    GL.glPixelStorei(GL.GL_PACK_ALIGNMENT, 1)
                    if tr.backBuffer:
                        GL.glReadBuffer(GL.GL_BACK)
                        c.swap = False
                    else:
                        GL.glReadBuffer(GL.GL_FRONT)
                        c.swap = True
                    # protect projection matrix
                    GL.glMatrixMode(GL.GL_PROJECTION)
                    GL.glPushMatrix()
                    GL.glMatrixMode(GL.GL_MODELVIEW)
                    tr.setup()
                    print tr.Columns,tr.Rows
                    if tr.checkeredBackground:
                        col =  0.0
                        colIncr = 1.0/(tr.Columns*tr.Rows)
                        originalBackground = c.backgroundColor
                    while moreTiles:
                        #print tileCount
                        tr.beginTile()
                        if tr.checkeredBackground:
                            #print 'color', col
                            c.Set(color = (col, col, col))
                            col += colIncr
                        c.InitGL()
                        if c.renderMode == GL.GL_RENDER:
                            c.Redraw()
                        self.tileCount+=1
                        #sleep(0.2)
                        moreTiles = tr.endTile()
                        print 'rendered tile %d'%self.tileCount
                        #raw_input("hit enter")
                    c.swap = True
                    if tr.checkeredBackground:
                        c.Set(color = originalBackground)

                    print "%d tiles drawn"%self.tileCount
                    # restore projection matrix
                    GL.glMatrixMode(GL.GL_PROJECTION)
                    GL.glPopMatrix()
                    GL.glMatrixMode(GL.GL_MODELVIEW)

                    im[i] = Image.fromstring('RGB', ( self.tileRenderFinalWidth,
                                                  self.tileRenderFinalHeight),
                                          self.tileRenderBuffer[i]) 
                    if os.name != 'nt': #sys.platform!='win32':
                        im[i] = im[i].transpose(Image.FLIP_TOP_BOTTOM)
                try:
                  if c.stereoMode.startswith('SIDE_BY_SIDE'):
                    imageFinal = Image.new('RGB', (w, self.tileRenderFinalHeight) )
                    imageFinal.paste( im[0], (0 , 0) )
                    imageFinal.paste( im[1], (self.tileRenderFinalWidth , 0) )
                    tileRenderedImage = imageFinal
                  else:
                    tileRenderedImage = im[0]

                  # if SSAO was on we will now do an additional redraw with
                  # SSAO on only_ssao mode to create the mask
                  if ssaoWasOn:
                      self.tileRender = False
                      c.ssao = True
                      old = c.SSAO_OPTIONS['only_ssao'][0]
                      c.SSAO_OPTIONS['only_ssao'][0] = True
                      c.Redraw()
                      # get the mask
                      aoim = c.GrabFrontBuffer(lock=False)
                      # scale it to cover match the large image
                      scale = self.tileRenderFinalWidth/float(c._width)
                      w, h = aoim.size
                      aoims = aoim.resize( (int(w*scale), int(h*scale)) )
                      # blur the mask
                      from PIL import ImageFilter, ImageChops
                      aoimsb = aoims.filter(ImageFilter.BLUR)
                      # multiply the large iamge and the mask
                      result = ImageChops.multiply(tileRenderedImage, aoimsb)
                      result.save(self.tileRenderOutput)
                      # restore SSAO mode
                      c.SSAO_OPTIONS['only_ssao'][0] = old
                  else:
                      tileRenderedImage.save(self.tileRenderOutput)
                except KeyError, e:
                  print "saving failed: unknown extension", e

                self.stopTileRendering()
                print 'done rendering tiles'

            # redraw full scenes
            c.InitGL()
            if c.renderMode == GL.GL_RENDER:
                c.Redraw()

        #self.needsRedraw = 0
        t2 = time()
        self.lastFrameTime = t2-t1


    def ReallyRedraw(self):
	"""Redraw the scene"""

        #print 'ReallyRedraw', self.pendingAutoRedrawID
        if self.checkIfRedrawIsNeeded() is False:
            #print 'NOOOOO'
            return True

        if len(self.cameras)==0:
            return True
        #import pdb
        #pdb.set_trace()

        #print 'ReallyRedraw', id(self)#self.pendingAutoRedrawID
        self.pendingAutoRedrawID = None

        #print 'Acquire =======================', id(self.redrawLock)
        self.redrawLock.acquire()
        #print 'Gotten'
        self.suspendRedraw = True

        try:
            # activate OpenGL context
            #c = self.cameras[0]
            c = self.currentCamera
            c.Activate()
            stillNeedRedraw = {}
            if len(self.objectsNeedingRedo):
                #free redraw lock in case beforeRebuildDpyLists tries to acquire it
                self.redrawLock.release()
                self.beforeRebuildDpyLists()
                self.redrawLock.acquire()

                for o in self.objectsNeedingRedo.keys():
                    if hasattr(o, 'justGotRedone'):
                        continue
                    if not o.visible:
                        stillNeedRedraw[o] = self.objectsNeedingRedo[o]
                        continue
                    o.justGotRedone = True
                    if o.immediateRendering \
                       or (self.tileRender and o.needsRedoDpyListOnResize):

                        # in immediate mode we can destroy any existing
                        # display list for this object
                        if o.dpyList is not None:
                            o.deleteOpenglList()
                        if o.transparent: #o.isTransparent():
                            if o in self.noDpyListOpaque:
                                self.noDpyListOpaque.remove(o)
                            if o not in self.noDpyListTransp:
                                self.noDpyListTransp.append(o)
                        else:
                            if o not in self.noDpyListOpaque:
                                self.noDpyListOpaque.append(o)
                            if o in self.noDpyListTransp:
                                self.noDpyListTransp.remove(o)
                    else:
                        self.currentCamera.Activate()

                        if o.viewer == None:
                            import pdb;pdb.set_trace()
                        
                        o.RedoDisplayList()

                for o in self.objectsNeedingRedo.keys():
                    # test because if object is more than once in list
                    # we only have to delete attribute once
                    if hasattr(o, 'justGotRedone'):
                        del o.justGotRedone

                # free redraw lock in case beforeRebuildDpyLists tried to acquire it
                self.redrawLock.release()
                self.afterRebuildDpyLists()
                self.redrawLock.acquire()

            self.objectsNeedingRedo = stillNeedRedraw
            
            # free redraw lock in case beforeRedraw tried to acquire it
            self.redrawLock.release()
            self.beforeRedraw()
            self.redrawLock.acquire()

            self.doRedraw()

            # free redraw lock in case afterRedraw tried to acquire it
            self.redrawLock.release()
            self.afterRedraw()
            self.redrawLock.acquire()

        finally:
            #print '================= Release'
            self.redrawLock.acquire(0)
            self.redrawLock.release()
        self.suspendRedraw = False

        self.postNextRedraw()
        #if self.autoRedraw:
        #    self.master.after(10, self.postNextRedrawCheck)
        return True


    def OneRedraw(self):

        if self.suspendRedraw is True:
            warnings.warn("invalid call to OneRedraw() when redraw is suspended")
            return

        # calls itself until a redraw was actually done
        if threading.currentThread().getName()!='MainThread':
            return

        restoreAutoRedraw = False
        if self.autoRedraw:
            #print 'stopping from OneRedraw'
            self.stopAutoRedraw()
            restoreAutoRedraw = True
            #raise RuntimeError('OneRedraw should NOT be called while Viewer is in AutoRedraw mode')
            
        self.needsRedraw = 1
        while True:
            done = self.ReallyRedraw()
            if done:
                break

        if restoreAutoRedraw:
            self.startAutoRedraw()
            
        
    def Redraw(self):
	"""Redraw the scene"""
        #print 'Redraw'
	#for c in self.cameras: c.Redraw()
        self.redrawLock.acquire()
        self.needsRedraw = True
        self.hasRedrawn = False
        self.redrawLock.release()
        

    def SetCurrentObject(self, obj):
        """Select an object to be the current object"""

        obj.hasBeenCurrent = 1
        
        if isinstance(obj, Transformable):
            # build R, Ri and Si
            obj.FrameTransform(self.currentCamera)

            if obj.scissor:
                self.currentScissor = obj

        if self.GUI:
            self.GUI.SetCurrentObject(obj)

        # generate an event
        self.eventHandler.dispatchEvent(SetCurrentObjectEvent(obj))

        # shouldn't the screen reflect the new selection ?
        #self.Redraw()


    def SetCurrentTexture(self, obj):
	"""Select a texture to be the current one"""

	# build R, Ri and Si
	obj.FrameTransform(self.currentCamera)
	self.currentTexture = obj


    def SetCurrentClip(self, obj):
	"""Select an Clipping plane to be the current one"""

        obj.hasBeenCurrent = 1
	# build R, Ri and Si
	obj.object.FrameTransform(self.currentCamera)
	obj.FrameTransform(self.currentCamera)
	oldcp = self.currentClip
	self.currentClip = obj
        if self.GUI:
            self.GUI.SetCurrentClip(oldcp)
        # generate an event
        self.eventHandler.dispatchEvent(SetCurrentClipEvent(obj))


    def SetCurrentCamera(self, obj):
	"""Select a Camera to be the current one"""

        obj.hasBeenCurrent = 1
	if obj == self.currentCamera: return
	self.currentCamera = obj

        if self.GUI:
            self.GUI.SetCurrentCamera()
        # generate an event
        self.eventHandler.dispatchEvent(SetCurrentCameraEvent(obj))


    def SetCurrentLight(self, obj):
	"""Select a Light to be the current one"""

        obj.hasBeenCurrent = 1
        obj.FrameTransform(self.currentCamera)
	if obj == self.currentLight: return
	self.currentLight = obj

	# update GUI here
        if self.GUI:
            self.GUI.SetCurrentLight(obj)
        # generate an event
        self.eventHandler.dispatchEvent(SetCurrentLightEvent(obj))


    def RotateCurrentObject(self, event, matrix, transXY, transZ):
        """Apply the rotation matrix
"""
        if isinstance(self.currentObject, Transformable):
            self.currentObject.ConcatRotation(matrix)
            self.Redraw()


    def TranslateCurrentObjectXY(self, event, matrix, transXY, transZ):
        """Add the translation trans to the current object
"""
#        c = self.currentCamera
#        t =  (transXY[0]*0.05, transXY[1]*0.05, 0.0 )
#        self.currentObject.ConcatTranslation( t )
#        self.Redraw()

        c = self.currentCamera
        lFovyFactor = c.fovy / c.fovyNeutral
        if lFovyFactor > 1:
            lFovyFactor = math.pow(lFovyFactor, 1.5)

        #lDistanceFactor = (c.far - c.near)  * .001
        # FIXME this will not work if we rotate the camera
        # using c.lookFrom[2] - c.lookAt[2] so taht translation is the same
        # regardless of distance of camera to scene
        lDistanceFactor = (c.lookFrom[2] - 0.0)  * .001

        lDistanceFactor = math.exp(lDistanceFactor) - 1
        lFactor = lDistanceFactor * lFovyFactor
        t =  (transXY[0] * lFactor, transXY[1]*lFactor, 0.0 )
        if c.currentTransfMode == 'Clip':
            self.currentClip.ConcatTranslation( t )
        else:
            self.currentObject.ConcatTranslation( t )
        self.Redraw()


    def screenTranslateCurrentObjectXY(self, event, matrix, transXY, transZ):
        """Add the translation trans to the current object
"""
        c = self.currentCamera
        lFovyFactor = c.fovy / c.fovyNeutral
        if lFovyFactor > 1:
            lFovyFactor = math.pow(lFovyFactor, 1.5)
        lDistanceFactor = (c.far - c.near)  * .001
        lDistanceFactor = math.exp(lDistanceFactor) - 1
        lFactor = lDistanceFactor * lFovyFactor
        t =  (transXY[0] * lFactor, transXY[1]*lFactor, 0.0 )
        if c.currentTransfMode == 'Clip':
            self.currentClip.screenConcatTranslation( t )
        else:
            # not implemented yet
            self.currentObject.ConcatTranslation( t )
        self.Redraw()


    def TranslateCurrentObjectZ(self, event, matrix, transXY, transZ):
        """Add the translation trans to the current object
"""
        if isinstance(self.currentObject, Transformable):
            c = self.currentCamera
            lFovyFactor = c.fovy / c.fovyNeutral
            if lFovyFactor > 1:
                lFovyFactor = math.pow(lFovyFactor, 1.5)
            lDistanceFactor = (c.far - c.near)  * .001
            lDistanceFactor = math.exp(lDistanceFactor) - 1
            lFactor = lDistanceFactor * lFovyFactor * .5
            t =  (0.0, 0.0, transZ*lFactor)
            if c.currentTransfMode == 'Clip':
                self.currentClip.ConcatTranslation( (transZ*lFactor, 0., 0.) )
            else:
                self.currentObject.ConcatTranslation( t )
            self.Redraw()


    def screenTranslateCurrentObjectZ(self, event, matrix, transXY, transZ):
        """Add the translation trans to the current object
"""
        if isinstance(self.currentObject, Transformable):
            c = self.currentCamera
            lFovyFactor = c.fovy / c.fovyNeutral
            if lFovyFactor > 1:
                lFovyFactor = math.pow(lFovyFactor, 1.5)
            lDistanceFactor = (c.far - c.near)  * .001
            lDistanceFactor = math.exp(lDistanceFactor) - 1
            lFactor = lDistanceFactor * lFovyFactor * .5
            t =  (0.0, 0.0, transZ*lFactor)
            if c.currentTransfMode == 'Clip':
                self.currentClip.screenConcatTranslation( t )
            else:
                # not implemented yet
                self.currentObject.ConcatTranslation( t )
            self.Redraw()


    def ScaleCurrentObject(self, event, matrix, transXY, transZ):
	"""Scale the current object"""
        # set pivot pivot
        if self.redirectTransformToRoot or \
               self.currentObject == self.rootObject:
            oldPivot = self.rootObject.SetPivot( self.lastPickedPixel )
        else:
            oldPivot = self.currentObject.SetPivot( self.lastPickedPixe )

        # scale
        scaleFactor = max(0.001, 1.0 + transZ*0.01)
        if scaleFactor < 0.001: scaleFactor = 0.001
	self.currentObject.ConcatScale(scaleFactor)

        # restore pivot
        if self.redirectTransformToRoot or \
               self.currentObject == self.rootObject:
            self.rootObject.SetPivot( oldPivot )
        else:
            self.currentObject.SetPivot( oldPivot )
	self.Redraw()


    def CurrentObjectFrontColor(self, prop, value):
	"""Method to be bound to a GUI widget"""

	if value:
	    self.currentObject.Set( materials = (value,), propName=prop )


    def CurrentObjectBackColor(self, prop, value):
	"""Method to be bound to a GUI widget"""

	if value:
	    self.currentObject.Set( materials = (value,), propName=prop,
                                    polyFace=GL.GL_BACK)


    def CurrentObjectOpacity(self, val):
	"""Method to be bound to a GUI widget"""

	if val:
	    self.currentObject.Set( opacity = val )


    def BindTrackballToClip(self, c, allCameras=None):
	"""Bind trackball to the current clipping plane"""

        self.useMasterDpyList = self.oldUseMasterDpyList

	if c is None: return
	if not isinstance(c, ClippingPlane):
	    raise AttributeError ('first parameter has to be an instance\
 of ClippingPlane')

	self.SetCurrentClip(c)

	if allCameras:
	    for c in self.cameras:
                c.bindAllActions('Clip')
##  		c.trackball.B2motion = self.RotateCurrentClipPlane
##  		c.trackball.B3motion = self.TranslateCurrentClipPlane
##  		c.trackball.ShiftB3motion = c.trackball.NoFunc
##  		c.trackball.ShiftB2motion = self.ScaleCurrentClipPlane
		self.Reset_cb = self.ResetCurrentClipPlane
		c.currentTransfMode = 'Clip'
	else:
	    c = self.currentCamera
            c.bindAllActions('Clip')
##  	    c.trackball.B2motion = self.RotateCurrentClipPlane
##  	    c.trackball.B3motion = self.TranslateCurrentClipPlane
##  	    c.trackball.ShiftB3motion = c.trackball.NoFunc
##  	    c.trackball.ShiftB2motion = self.ScaleCurrentClipPlane
	    self.Reset_cb = self.ResetCurrentClipPlane
	    c.currentTransfMode = 'Clip'

        self.Xform.set(c.currentTransfMode)
        if self.GUI:
            self.GUI.BindTrackballToClip(c, allCameras)


    def ResetCurrentClipPlane(self):
	"""Reset the current clipping plane's transformation"""
	self.currentClip.ResetTransformation()
	self.Redraw()


    def RotateCurrentClipPlane(self, event, matrix, transXY, transZ):
	"""Apply rotation matrix to the current clipping plane"""
	c = self.currentClip
	c.ConcatRotation(matrix)
	self.Redraw()


#    def TranslateCurrentClipPlane(self, event, matrix, transXY, transZ):
#	"""Apply translation trans to the current clipping plane"""
#
##        clip = self.currentClip
##        t =  (transXY[0]*0.05, transXY[1]*0.05, 0.0 )
##        clip.ConcatTranslation( t )
##        self.Redraw()
#
#        c = self.currentCamera
#        lFovyFactor = c.fovy / c.fovyNeutral
#        if lFovyFactor > 1:
#            lFovyFactor = math.pow(lFovyFactor, 1.5)
#        lDistanceFactor = (c.far - c.near)  * .001
#        lDistanceFactor = math.exp(lDistanceFactor) - 1
#        lFactor = lDistanceFactor * lFovyFactor
#        t =  (transXY[0] * lFactor, transXY[1]*lFactor, 0.0 )
#        self.currentClip.ConcatTranslation( t )
#        self.Redraw()


    def ScaleCurrentClipPlane(self, event, matrix, transXY, transZ):
	"""Scale clipping plane representation"""
        scaleFactor = max(0.001, 1.0 + transZ*0.01)
        if scaleFactor < 0.001: scaleFactor = 0.001
	self.currentClip.ConcatScale(scaleFactor)
	self.Redraw()


    def CurrentClipColor(self, col):
	"""Method to be bound to a GUI widget"""

	if col:
	    self.currentClip.Set( color = col )
	    self.currentCamera.Redraw()


    def BindTrackballToCamera(self, c, allCameras=None):
	"""Bind the trackball to the current camera"""

        self.useMasterDpyList = self.oldUseMasterDpyList

	if not isinstance(c, Camera):
	    raise AttributeError('first parameter has to be an instance of \
Camera')
        self.SetCurrentCamera(c)
	if allCameras:
	    for c in self.cameras:
                c.bindAllActions('Camera')
##  		c.trackball.B2motion = self.RotateCurrentCamera
##  		c.trackball.B3motion = self.TranslateCurrentCamera
##  		c.trackball.ShiftB3motion = c.trackball.NoFunc
##                  c.trackball.ShiftB2motion = self.ScaleCurrentCamera
		self.Reset_cb = self.ResetCurrentCamera
		c.currentTransfMode = 'Camera'
	else:
            c.bindAllActions('Camera')
##  	    c.trackball.B2motion = self.RotateCurrentCamera
##  	    c.trackball.B3motion = self.TranslateCurrentCamera
##  	    c.trackball.ShiftB3motion = c.trackball.NoFunc
##  	    c.trackball.ShiftB2motion = self.ScaleCurrentCamera
	    self.Reset_cb = self.ResetCurrentCamera
	    c.currentTransfMode = 'Camera'

        self.Xform.set(c.currentTransfMode)
        if self.GUI:
            self.GUI.BindTrackballToCamera(c, allCameras)


##      def RotateCameraLookAt(self, event, matrix):
##          self.currentCamera.ConcatLookAtRot(matrix)
##  	for l in self.lights:
##  	    if self.currentCamera.renderMode == GL.GL_SELECT:
##                  continue
##              if l.positional==viewerConst.YES:
##                  l.posFlag = 1
##                  l.spotFlag = 1
##              else:
##                  l.dirFlag = 1
##          self.Redraw()

##      # not implemented now; will need to fix light if it is used (see fn above)
##      def TranslateCameraLookAt(self, event, trans):
##          self.currentCamera.ConcatLookAtTrans(trans)
##          self.Redraw()
        

    def RotateCurrentCamera(self, event, matrix, transXY, transZ):
	"""Apply the rotation matrix to camera.  If lights, the flag is set
        to redraw the lights in a new position."""
	#self.currentCamera.ConcatLookFromRot(matrix)
        #rot = numpy.transpose(numpy.reshape(matrix, (4,4)))
        #rot = numpy.reshape(rot, (16,))
        self.currentCamera.ConcatRotation(matrix)
	for l in self.lights:
	    if self.currentCamera.renderMode == GL.GL_SELECT:
                continue
            if l.positional==viewerConst.YES:
                l.posFlag = 1
                l.spotFlag = 1
            else:
                l.dirFlag = 1
            #l.ConcatRotationDir(matrix)
            #l.ConcatRotationDir(rot)
##              for type in ('dir', 'pos', 'spot'):
##                  setattr(l, '%sFlag'%type, 1)
	self.Redraw()


    def TranslateCurrentCamera(self, event, matrix, transXY, transZ):
        """Add the translation trans to the current object
"""
        #t =  (- transXY[0]*0.05, - transXY[1]*0.05, -transZ )
        t =  (0, 0, -transZ )
        self.currentCamera.ConcatTranslation( t )
        self.Redraw()


    def TranslateXYCurrentCamera(self, event, matrix, transXY, transZ):
        """Add the translation trans to the current object
"""
        t =  (transXY[0]*0.05, transXY[1]*0.05, 0. )
        self.currentCamera.ConcatTranslation( t )
        self.Redraw()


    def translateCurrentCameraMouseWheel(self, event):
        """translate the current camera via the mouse wheel
"""
        #print "translateCurrentCameraMouseWheel", event.num
        if event.num == 4:
            t =  (0, 0, -10 )
        else:
            t =  (0, 0, 10 )
        self.currentCamera.ConcatTranslation( t )
        self.Redraw()

        
    def ScaleCurrentCamera(self, event, matrix, transXY, transZ):
	"""Scale the current object
"""
        #print "ScaleCurrentCamera", event, matrix, transXY, transZ
        # set picked pixel as pivot
        c = self.currentCamera
        value = c.fovy * .01 * transZ
        fovy = c.fovy - value
        print c.fovy, transZ, value, fovy
        if (fovy < 180.) and (fovy > 0.):
            c.Set(fov=fovy)
        self.Redraw()


    def scaleCurrentCameraMouseWheel(self, event):
        """Scale the current object
"""
        #print "scaleCurrentCameraUp", event.num
        if os.name == 'nt': #sys.platform == 'win32':
            if event.delta > 0:
                lEventNum = 4
            else:
                lEventNum = 5
        else:
            lEventNum = event.num
        c = self.currentCamera
        value = c.fovy * .1
        if lEventNum == 4:
            fovy = c.fovy - value
            c.Set(fov=fovy)
        else:
            fovy = c.fovy + value
            if fovy < 180.:
                c.Set(fov=fovy)
        self.Redraw()


    def ResetCurrentCamera(self):
	"""Reset the current camera's transformation"""
	self.currentCamera.ResetTransformation()
	for l in self.lights:
	    if self.currentCamera.renderMode == GL.GL_SELECT:
                continue
            if l.positional==viewerConst.YES:
                l.posFlag = 1
                l.spotFlag = 1
            else:
                l.dirFlag = 1
 	self.Redraw()


    def BindTrackballToLight(self, light, allCameras=None):
	"""Bind the trackball to the current ligth source"""
	
        self.useMasterDpyList = self.oldUseMasterDpyList

	if not isinstance(light, Light):
	    raise AttributeError('first parameter has to be an instance of \
Light')
        self.SetCurrentLight(light)

	if light.positional==viewerConst.NO:  # directional light
	    if allCameras:
		for c in self.cameras:
                    c.bindAllActions('Light')
##  		    c.trackball.B2motion = self.RotateCurrentDLight
##  		    c.trackball.B3motion = c.trackball.NoFunc
##  		    c.trackball.ShiftB3motion = c.trackball.NoFunc
##  		    c.trackball.ShiftB2motion = c.trackball.NoFunc
		    self.Reset_cb = self.ResetCurrentDLight
		    c.currentTransfMode = 'Light'
	    else:
		c = self.currentCamera
                c.bindAllActions('Light')
##  		c.trackball.B2motion = self.RotateCurrentDLight
##  		c.trackball.B3motion = c.trackball.NoFunc
##  		c.trackball.ShiftB3motion = c.trackball.NoFunc
##  		c.trackball.ShiftB2motion = c.trackball.NoFunc
		self.Reset_cb = self.ResetCurrentDLight
		c.currentTransfMode = 'Light'

	else:  # positional light
	    light.FrameTransform(self.currentCamera)

        self.Xform.set(c.currentTransfMode)
        if self.GUI:
            self.GUI.BindTrackballToLight(light, allCameras)


    def RotateCurrentDLight(self, event, matrix, transXY, transZ):
	self.currentLight.ConcatRotationDir(matrix)
	self.Redraw()


    def ResetCurrentDLight(self):
	self.currentLight.Set(direction = (1,1,1,1))
	self.Redraw()


    def CurrentLightColor(self, prop, value):
	"""Method to be bound to a GUI widget"""

	if value:
	    if prop=='ambient': self.currentLight.Set( ambient = value )
	    elif prop=='diffuse': self.currentLight.Set( diffuse = value )
	    elif prop=='specular': self.currentLight.Set( specular = value )
	    self.currentCamera.Redraw()


    def LMColor(self, col):
	"""Method to be bound to a GUI widget"""

	if col:
	    self.lightModel.Set( ambient=col, tagModified=False )
	    self.currentCamera.Redraw()


    def BindTrackballToTexture(self, o, allCameras=None):
        """Bind trackball to the texture of the current object
"""

        # we cannot use masterDpyList else Xform matrix for texture does not
        # get updated in DrawOneObject (since it is embbeded in masterDpyList
        self.oldUseMasterDpyList = self.useMasterDpyList
        self.useMasterDpyList = 0

        if o is None:
            self.Xform.set(self.currentCamera.currentTransfMode)
            return
        if not isinstance(o, Geom):
            raise AttributeError ('first parameter has to be an instance of Geom')
        if not o.texture:
            self.Xform.set(self.currentCamera.currentTransfMode)
            return
        self.SetCurrentTexture(o.texture)

        if allCameras:
            for c in self.cameras:
                c.bindAllActions('Texture')
                #c.trackball.B2motion = self.RotateCurrentTexture
                #c.trackball.B3motion = self.TranslateCurrentTextureXY
                #c.trackball.ShiftB3motion = self.TranslateCurrentTextureZ
                #c.trackball.ShiftB2motion = self.ScaleCurrentTexture
                self.Reset_cb = self.ResetCurrentTexture
                c.currentTransfMode = 'Texture'
        else:
            c = self.currentCamera
            c.bindAllActions('Texture')
            #c.trackball.B2motion = self.RotateCurrentTexture
            #c.trackball.B3motion = self.TranslateCurrentTextureXY
            #c.trackball.ShiftB3motion = self.TranslateCurrentTextureZ
            #c.trackball.ShiftB2motion = self.ScaleCurrentTexture
            self.Reset_cb = self.ResetCurrentTexture
            c.currentTransfMode = 'Texture'

        self.Xform.set(c.currentTransfMode)
        if self.GUI:
            self.GUI.BindTrackballToTexture(o, allCameras)


    def ResetCurrentTexture(self):
	"""Reset the current tesxture's transformation"""
	self.currentTexture.ResetTransformation()
	self.Redraw()


    def ScaleCurrentTexture(self, event, matrix, transXY, transZ):
	"""Scale the current object"""
        scaleFactor = max(0.001, 1.0 + transZ*0.01)
        if scaleFactor < 0.001: scaleFactor = 0.001
	self.currentTexture.ConcatScale(1.0/scaleFactor)
	self.Redraw()


    def RotateCurrentTexture(self, event, matrix, transXY, transZ):
	"""Apply the rotation matrix"""
	self.currentTexture.ConcatRotation(matrix)
	self.Redraw()


    def TranslateCurrentTextureXY(self, event, matrix, transXY, transZ):
	"""Add the translation trans to the current Texture"""
        t =  (-transXY[0]*0.05, -transXY[1]*0.05, 0.0 )
	self.currentTexture.ConcatTranslation(t)
	self.Redraw()


    def TranslateCurrentTextureZ(self, event, matrix, transXY, transZ):
	"""Add the translation trans to the current Texture"""
        t =  (0.0, 0.0, -transZ*0.05)
	self.currentTexture.ConcatTranslation(t)
	self.Redraw()


    def BindTrackballToScissor(self, o, allCameras=None):
        """Bind trackball to the scissor of the current object
"""
        if o is None:
            self.Xform.set(self.currentCamera.currentTransfMode)
            return
        if not isinstance(o, Geom):
            raise AttributeError('first parameter has to be an instance of Geom')
        if not o.scissor:
            self.Xform.set(self.currentCamera.currentTransfMode)
            return

        self.useMasterDpyList = self.oldUseMasterDpyList
        self.currentScissor = o
	
        if allCameras:
            for c in self.cameras:
                c.bindAllActions('Scissor')
                #c.trackball.B2motion = c.trackball.NoFunc
                #c.trackball.B3motion = self.TranslateCurrentScissor
                #c.trackball.ShiftB3motion = self.AspectRatioScissor
                #c.trackball.ShiftB2motion = self.ScaleCurrentScissor
                self.Reset_cb = self.ResetCurrentScissor
                c.currentTransfMode = 'Scissor'
        else:
            c = self.currentCamera
            c.bindAllActions('Scissor')
            #c.trackball.B2motion = c.trackball.NoFunc
            #c.trackball.B3motion = self.TranslateCurrentScissor
            #c.trackball.ShiftB3motion = self.AspectRatioScissor
            #c.trackball.ShiftB2motion = self.ScaleCurrentScissor
            self.Reset_cb = self.ResetCurrentScissor
            c.currentTransfMode = 'Scissor'

        self.Xform.set(c.currentTransfMode)
        if self.GUI:
            self.GUI.BindTrackballToScissor(o, allCameras)


    def ResetCurrentScissor(self):
	"""Reset the current tesxture's transformation"""
	self.currentScissor.Set( scissorX=0, scissorY=0,
                                 scissorW=200, scissorH=200 )
	self.Redraw()


    def AspectRatioScissor(self, event, matrix, transXY, transZ):
        """modify scissor aspect ratio"""
        scale = max(0.001, 1.0 + transZ*0.01)
        sc = self.currentScissor
        sc.Set(scissorAspectRatio = sc.scissorAspectRatio*scale)
	self.Redraw()

        
    def ScaleCurrentScissor(self, event, matrix, transXY, transZ):
	"""Scale the current object"""
        scaleFactor = max(0.001, 1.0 + transZ*0.01)
        if scaleFactor < 0.001: scaleFactor = 0.001
        sc = self.currentScissor
	sc.Set( scissorW=int(sc.scissorW*scaleFactor),
                scissorH=int(sc.scissorH*scaleFactor) )
	self.Redraw()


    def TranslateCurrentScissor(self, event, matrix, transXY, transZ):
	"""Add the translation trans to the current Scissor"""
        sc = self.currentScissor
	sc.Set( scissorX=int(sc.scissorX+transXY[0]),
                scissorY=int(sc.scissorY+transXY[1]) )
	self.Redraw()


    def deleteOpenglList(self):
        #import traceback;traceback.print_stack() 
        #print "Viewer.deleteOpenglList"
        self.redoDspLst = 0

        for c in self.cameras:
            if c is not None and c.dpyList is not None:
                c.Activate()
                #currentcontext = self.currentCamera.getContext()
                currentcontext = c.getContext()
                if currentcontext != c.dpyList[1]:
                    warnings.warn("""deleteOpenglList failed because the current context is the wrong one""")
                    #print "currentcontext != self.dpyList[1]", currentcontext, self.dpyList[1]
                else:
                    #print '-%d'%c.dpyList[0], currentcontext, "glDeleteLists Viewer"
                    GL.glDeleteLists(c.dpyList[0], 1)
                    c.dpyList = None


    def deleteOpenglListAndCallRedraw(self):
        self.deleteOpenglList()
        self.Redraw()


#    def deleteAllTheDisplayListAndCallRedraw(self):
#        # except the templates
#        for o in self.rootObject.AllObjects():
#            if o.dpyList is not None:
#                o.RedoDisplayList()
#        self.deleteOpenglList()
#        self.Redraw()


    def removeAllTheDisplayListsExceptTemplatesAndVBO(self):
        self.objectsNeedingRedo = {}
        for g in self.rootObject.AllObjects():
            if not (hasattr(g, 'vertexArrayFlag') and \
               g.vertexArrayFlag is True and \
               hasattr(g, 'vertexSet') is True and \
               DejaVu2.enableVBO is True):
                g.deleteOpenglList()
                self.objectsNeedingRedo[g] = None        
        self.deleteOpenglList()


    def deleteOpenglListAndCallRedrawAndCallDisableGlLighting(self):
        #print "deleteOpenglListAndCallRedrawAndCallDisableGlLighting"
        GL.glDisable(GL.GL_LIGHTING)
        self.deleteOpenglListAndCallRedraw()


    def transformedCoordinatesWithInstances(self, hits):
        """ hist is pick.hits = {geom: [(vertexInd, intance),...]}
This function will use the instance information to return a list of transformed
coordinates
"""
        vt = []
        for geom, values in hits.items():
            coords = geom.vertexSet.vertices.array
            for vert, instance in values:
                M = geom.GetMatrix(geom.LastParentBeforeRoot(), instance[1:])
                pt = coords[vert]
                ptx = M[0][0]*pt[0]+M[0][1]*pt[1]+M[0][2]*pt[2]+M[0][3]
                pty = M[1][0]*pt[0]+M[1][1]*pt[1]+M[1][2]*pt[2]+M[1][3]
                ptz = M[2][0]*pt[0]+M[2][1]*pt[1]+M[2][2]*pt[2]+M[2][3]
                vt.append( (ptx, pty, ptz) )
        return vt


    def enableOpenglLighting(self):
        #print "enableOpenglLighting", self.OverAllLightingIsOn.get()
        if self.OverAllLightingIsOn == 1:
            #print "Opengl lighting enabled"
            GL.glEnable(GL.GL_LIGHTING)


    def toggleOpenglLighting(self, event=None):
        #print "toggleOpenglLighting"
        self.OverAllLightingIsOn = not self.OverAllLightingIsOn
        self.deleteOpenglListAndCallRedrawAndCallDisableGlLighting()


    ## NEW METHODS
        
