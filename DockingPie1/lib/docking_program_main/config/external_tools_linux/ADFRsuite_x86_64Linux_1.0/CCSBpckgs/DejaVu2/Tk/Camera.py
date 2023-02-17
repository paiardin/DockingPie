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
# $Header: /mnt/raid/services/cvs/DejaVu2/Tk/Camera.py,v 1.1.1.1.4.1 2017/07/13 22:22:07 annao Exp $
#
# $Id: Camera.py,v 1.1.1.1.4.1 2017/07/13 22:22:07 annao Exp $
#
import Tkinter, Pmw, os, types

from opengltk.OpenGL import GL

import DejaVu2
from DejaVu2.Tk import loadTogl
from DejaVu2.Tk.EventHandler import EventManager
from DejaVu2.Tk.Trackball import Trackball
from DejaVu2.Camera import StandardCameraBase
from DejaVu2.viewerFns import checkKeywords
from DejaVu2.cursors import cursorsDict

from mglutil.gui import widgetsOnBackWindowsCanGrabFocus


## class DynamicComboBox(Pmw.ComboBox):

##     def __init__(self, viewer, parent = None, **kw):
##         if __debug__:
##          if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
##         self.viewer = viewer
##         Pmw.ComboBox.__init__(self, parent, **kw)
        
##     def _postList(self, event=None):
##         if __debug__:
##          if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
##         geomNames = []
##         for g in self.viewer.rootObject.AllObjects():
##             if isinstance(g, IndexedPolygons):
##                 geomNames.append(g.fullName)
        
##         self.component('scrolledlist').setlist(geomNames)
##         Pmw.ComboBox._postList(self, event)



class Camera(Tkinter.Widget, Tkinter.Misc, StandardCameraBase):

    def __init__(self, master, screenName, viewer, num, check=1,
                 cnf={}, **kw):
        
# Once togl will allow to create widgets with same context
# but for now we ignore it ==> no more multiple cameras
#                kw['samecontext'] = 1
#
# Togl Tk widget options as of Version 1.5
# -height [DEFAULT_WIDTH=400], -width [DEFAULT_HEIGHT=400],
# -rgba [true]
# -redsize [1], -greensize [1], -bluesize [1]
# -double [false]
# -depth [false]
# -depthsize [1]
# -accum [false]
# -accumredsize [1], -accumgreensize [1], -accumbluesize [1], -accumalphasize [1]
# -alpha [false], -alphasize [1]
# -stencil [false], -stencilsize [1]
# -auxbuffers [0]
# -privatecmap [false]
# -overlay [false]
# -stereo [false]
# -time [DEFAULT_TIME = "1"]
# -sharelist [NULL]
# -sharecontext [NULL]
# -ident [DEFAULT_IDENT = ""]

        self.name = 'Camera'+str(num)
        self.num = num
        self.uniqID = self.name+viewer.uniqID

        if __debug__:
            if check:
                apply( checkKeywords, (self.name,self.initKeywords), kw)

        if not kw.has_key('double'): kw['double'] = 1
#        if not kw.has_key('overlay'): kw['overlay'] = 1
        if not kw.has_key('depth'): kw['depth'] = 1

        if not kw.has_key('stencil'): kw['stencil'] = 1

        if hasattr(viewer, 'accumBuffersError') and viewer.accumBuffersError:
            # if we tried to create a context with accumulation buffers before
            # and it failed we set accum to False
            kw['accum'] = 0
        else:
            # if the user did not specify accumulation we turn them on by defaut
            if not kw.has_key('accum'):
                kw['accum'] = 1

        if not kw.has_key('stereo'): kw['stereo'] = 'none'

        if not kw.has_key('ident'):
            kw['ident'] = self.uniqID
#            kw['ident'] = 'camera%d' % num

        if not kw.has_key('sharelist') and len(viewer.cameras):
            # share list with the default camera.
            cam = viewer.cameras[0]
            kw['sharelist'] = cam.uniqID

        if not kw.has_key('sharecontext') and len(viewer.cameras):
            # share context with the default camera.
            cam = viewer.cameras[0]
            kw['sharecontext'] = cam.uniqID
            if not hasattr(cam, 'shareCTXWith'): cam.shareCTXWith = []
            cam.shareCTXWith.append(self)

##         if master is None:
##             from os import path
##             from opengltk.OpenGL import Tk
##             toglInstallDir = path.dirname(path.abspath(Tk.__file__))
##             tclIncludePath = master.tk.globalgetvar('auto_path')
##             master.tk.globalsetvar('auto_path', toglInstallDir + ' ' +
##                                    tclIncludePath)
##             master.tk.call('package', 'require', 'Togl')

        self.frameBorderWidth = 3
        self.frame = Tkinter.Frame(master, bd=self.frameBorderWidth)

        #self.frame.master.protocol("WM_DELETE_WINDOW", self.hide)

        cfg = 0

        if 'width' in cnf.keys():
            self._width = cnf['width']
            cfg = 1
        else:
            cnf['width'] = self._width = 406

        if 'height' in cnf.keys():
            self._height = cnf['height']
            cfg = 1
        else:
            cnf['height'] = self._height = 406

        self.rootx = 320
        if 'rootx' in cnf.keys():
            self.rootx = cnf['rootx']
            del cnf['rootx']
            cfg = 1
        self.rooty = 180
        if 'rooty' in cnf.keys():
            self.rooty = cnf['rooty']
            del cnf['rooty']
            cfg = 1
        if cfg: self.Geometry()

        if 'side' in cnf.keys():
            side = cnf['side']
            del cnf['side']
        else: side = 'top'

        self.frame.pack(fill=Tkinter.BOTH, expand=1, side=side)
        self.defFrameBack = self.frame.config()['background'][3]

        toglVersion = loadTogl(self.frame)
		#print "Togl Version:", toglVersion

        # after this line self.master will be set to self frame
        from Tkinter import TclError
        from opengltk.exception import GLerror
        try:
            Tkinter.Widget.__init__(self, self.frame, 'togl', cnf, kw)
            try:
                GL.glAccum(GL.GL_LOAD, 1.0)
                viewer.accumBuffersError = False
            except GLerror:
                viewer.accumBuffersError = True
                #viewer.accumBuffersError = False
        except TclError, e:
            print "Warning: disabling accumulation buffers",e
            viewer.accumBuffersError = True
            kw.pop('accum')
            Tkinter.Widget.__init__(self, self.frame, 'togl', cnf, kw)

        self.preventIntelBug_BlackTriangles()
        
        # create a TK-event manager for this camera
        self.eventManager = EventManager(self)
        self.eventManager.AddCallback('<Map>', self.Map)
        self.eventManager.AddCallback('<Expose>', self.Expose)
        self.eventManager.AddCallback('<Configure>', self.Expose)
        self.eventManager.AddCallback('<Enter>', self.Enter_cb)

        StandardCameraBase.__init__(self, master, screenName, viewer, num, **kw)
        self.pack(side='left', expand=1, fill='both')

        if os.name == 'nt': #sys.platform == 'win32':
            self.eventManager.AddCallback(
                "<MouseWheel>", viewer.scaleCurrentCameraMouseWheel)
                #"<MouseWheel>", viewer.translateCurrentCameraMouseWheel)
            #self.bind("<MouseWheel>", viewer.scaleCurrentCameraMouseWheel)
        else:
            self.eventManager.AddCallback(
                "<Button-4>", viewer.scaleCurrentCameraMouseWheel)
                #"<Button-4>", viewer.translateCurrentCameraMouseWheel)
            self.eventManager.AddCallback(
                "<Button-5>", viewer.scaleCurrentCameraMouseWheel)
                #"<Button-5>", viewer.translateCurrentCameraMouseWheel)
            #self.bind("<Button-4>", viewer.scaleCurrentCameraMouseWheel)
            #self.bind("<Button-5>", viewer.scaleCurrentCameraMouseWheel)

        self.AddTrackball()
        for binding in ['Object', 'Insert2d', 'Camera', 'Clip', 'Light',
                        'Texture', 'Scissor']:
            self.actions[binding]['None'] = self.trackball.NoFunc


    def Map(self, *dummy):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Cause the opengl widget to redraw itself."""
        self.Expose()


##     def Configure(self, *dummy):
##         """Cause the opengl widget to redraw itself."""
##         self.Activate()
##         self._width = self.winfo_width()
##         self._height = self.winfo_height()
##         self.rootx = self.winfo_rootx()
##         self.rooty = self.winfo_rooty()
##         # FIXME .. this is not symetrical ... should we just recompute left, right,
##         # top and botom here ???
##         if self.projectionType==self.ORTHOGRAPHIC:
##             self.PerspectiveToOrthogonal()
##         else:
##             self.SetupProjectionMatrix()

##         self.Expose()


    def Expose(self, event=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """set the camera's exposeEvent so that at the next redraw
the camera width and height are updated
"""
        if not self.viewer.isInitialized:
            self.after(100, self.Expose)
        else:
            # if viewer is in autoRedraw mode the next redraw will handle it
            if not self.exposeEvent and self.viewer.autoRedraw:
                self.viewer.Redraw()
                self.exposeEvent = True

            for o in self.viewer.rootObject.AllObjects():
                if o.needsRedoDpyListOnResize or o.scissor:
                    self.viewer.objectsNeedingRedo[o] = None


    def Set(self, check=1, redo=1, **kw):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Set various camera parameters
"""
        StandardCameraBase.Set(self, check=check, redo=redo, **kw)

        w = kw.get( 'width')
        if not w is None:
            assert type(w)==types.IntType
            if w > 0:
                self._width = w
            else:
                raise AttributeError('width has to be > 0')

        h = kw.get( 'height')
        if not h is None:
            assert type(h)==types.IntType
            if h > 0:
                # we disable redraw because autoRedraw could cause ReallyRedraw
                # to set camera._width and camera._height BEFORE te window gets
                # resized by Tk
                self._height = h
            else:
                raise AttributeError('height has to be > 0')

        px = kw.get( 'rootx')
        if not px is None:
            px = max(px, 0)
            master = self.frame.master
            while hasattr(master, 'wm_maxsize') is False:
                master = master.master
            xmax, ymax = master.wm_maxsize()
            px = min(px, xmax -100)
            self.rootx = px

        py = kw.get( 'rooty')
        if not py is None:
            py = max(py, 0)
            master = self.frame.master
            while hasattr(master, 'wm_maxsize') is False:
                master = master.master
            xmax, ymax = master.wm_maxsize()
            py = min(py, ymax - 100)
            self.rooty = py

        # if width, height of position of camera changed apply changes
        if w or h or px or py:
            self.Geometry()
            self.viewer.update()


    def AddTrackball(self, size=0.8, rscale=2.0, tscale=0.05,
                     sscale=0.01, renorm=97):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Add a Trackball to this camera with default bindings
"""
        #print "AddTrackball"

        self.trackball = Trackball(self, size, rscale, tscale, sscale, renorm )

    
    def setCursor(self, buttonNum):
        modDict = self.viewer.kbdModifier
        xform = self.viewer.Xform.get()
        if modDict['Shift_L']:
            action = self.mouseButtonActions[xform][buttonNum]['Shift']
        elif modDict['Control_L']:
            action = self.mouseButtonActions[xform][buttonNum]['Control']
        elif modDict['Alt_L']:
            action = self.mouseButtonActions[xform][buttonNum]['Alt']
        else:
            action = self.mouseButtonActions[xform][buttonNum]['None']
            
        #print 'AAAAAAAAAAAA', action#, modDict
        if cursorsDict.has_key(action):
            self.configure(cursor=cursorsDict[action])


    def bindAllActions(self, actionDict):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        #print "bindAllActions", actionDict
        for b in (1,2,3):
            d = self.mouseButtonActions[actionDict][b]
            for mod in self.mouseButtonModifiers:
                self.bindActionToMouseButton(d[mod], b, mod, actionDict)


    # FIXME picking and other mouse actions should be unified
    # Press events register the motion and release callbacks
    # this will enable to process picking and others in the same way
    def bindActionToMouseButton(self, action, buttonNum, modifier='None',
                                actionDict='Object'):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """registers callbacks to trigger picking when the specified mouse
button is pressed. buttonNum can be 1,2,3 or 0 to remove call backs
"""
        #print "bindActionToMouseButton", buttonNum, action, modifier, actionDict

        ehm = self.eventManager
        func = self.actions[actionDict][action]

        # find button and modifier to which actions is bound currently
        oldnum, oldmod = self.findButton(action, actionDict)

        if action in ['picking', 'addToSelection', 'removeFromSelection',
                      'pivotOnPixel']:
            if oldnum: # picking action was bound to a mouse button
                # remove the picking callbacks from previous picking button
                if modifier=='None': modifier1=''
                else: modifier1=modifier+'-'
                ev = '<'+modifier1+'ButtonPress-'+str(oldnum)+'>'
                ehm.RemoveCallback(ev, func)
                    
            # remove all motion call back on buttonNum for all modifiers
            if modifier=='None': modifier1=''
            else: modifier1=modifier+'-'
            setattr(self.trackball, modifier1+'B'+str(buttonNum)+'motion',
                    self.trackball.NoFunc)

            # now bind picking to buttonNum for all modifiers
            self.mouseButtonActions[actionDict][buttonNum][modifier] = action
            if modifier=='None': modifier1=''
            else: modifier1=modifier+'-'
            ev = '<'+modifier1+'ButtonPress-'+str(buttonNum)+'>'
            ehm.AddCallback(ev, func)

        else: # not picking
            if oldnum: # action was bound to a mouse button
                # remove the picking callbacks from previous picking button
                if oldmod=='None': mod1=''
                else: mod1=oldmod+'-'
                ev = '<'+mod1+'ButtonPress-'+str(oldnum)+'>'
                oldaction = self.mouseButtonActions[actionDict][buttonNum][oldmod]
                if oldaction=='picking' or oldaction=='pivotOnPixel':
                    ehm.RemoveCallback(ev, self.actions[actionDict][oldaction])

                # remove motion callback on oldbuttonNum for oldmodifier
                if oldmod=='None': mod=''
                else: mod=oldmod
                setattr(self.trackball, mod+'B'+str(oldnum)+'motion',
                        self.trackball.NoFunc)
                self.mouseButtonActions[actionDict][oldnum][oldmod] = 'None'

            # remove picking callback on buttonNum for modifier
            if modifier=='None': mod=''
            else: mod=modifier+'-'
            oldaction = self.mouseButtonActions[actionDict][buttonNum][modifier]
            if oldaction=='picking' or oldaction=='pivotOnPixel':
                ev = '<'+mod+'ButtonPress-'+str(buttonNum)+'>'
                ehm.RemoveCallback(ev, self.actions[actionDict][oldaction])

            
            # now bind picking to buttonNum for all modifiers
            if modifier=='None': mod=''
            else: mod=modifier
            setattr(self.trackball, mod+'B'+str(buttonNum)+'motion', func)
            self.mouseButtonActions[actionDict][buttonNum][modifier] = action


    def __del__(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Destroy the camera
"""
        #print "StandardCamera.__del__", self
        self.frame.master.destroy()


    def storeGeometry(self):
        self.rootx = self.winfo_rootx()
        self.rooty = self.winfo_rooty()
        self._width = self.winfo_width()
        self._height = self.winfo_height()
        
        
    def getGeometry(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """return the posx, posy, width, height of the window containing the camera"""
        geom = self.winfo_geometry()
        size, x, y = geom.split('+')
        w, h = size.split('x')
        return int(x), int(y), int(w), int(h)


    def Geometry(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """
Resize the Tk widget holding the camera to:
     self._widthxself._height+self.rootx+self.rooty
This only applies when the camera is in its own top level.
"""

        # get a handle to the master of the frame containing the Togl widget
        window = self.frame.master

        # if window is a toplevel window
        
        # replaced this test by test for geometry method
        #if isinstance(, Tkinter.Tk) or \
        #   isinstance(self.frame.master, Tkinter.Toplevel):
        if hasattr(window, 'geometry') and callable(window.geometry):
            # we have to set the geoemtry of the window to be the requested
            # size plus 2 times the border width of the frame containing the
            # camera
            off = 2*self.frameBorderWidth
            geom = '%dx%d+%d+%d' % (self._width+off, self._height+off,
                                    self.rootx, self.rooty)
            window.geometry(geom)


    def SelectCamera(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """First pick in a non current camera selects it"""

        curCam = self.viewer.currentCamera
        if curCam == self: return 0
        curCam.frame.config( background = self.defFrameBack )
        self.frame.config( background = "#900000" )
        if self.viewer:
            self.viewer.currentCamera = self
            self.viewer.BindTrackballToObject(self.viewer.currentObject)

        self.SetupProjectionMatrix()
        return 1


    def lift(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """brings the window containing the camera in front of others"""
        window = self.frame.master
        if isinstance(window, Tkinter.Tk) or \
           isinstance(window, Tkinter.Toplevel):
            self.frame.master.lift()
        else:
            m = self.master
            while m.master:
                m = m.master
            m.lift()

    
    def Enter_cb(self, event=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Call back function trigger when the mouse enters the camera
"""
        if widgetsOnBackWindowsCanGrabFocus is False:
            lActiveWindow = self.focus_get()
            if    lActiveWindow is not None \
              and ( lActiveWindow.winfo_toplevel() != self.winfo_toplevel() ):
                return

        self.focus_set()
        self.tk.call(self._w, 'makecurrent')
        if not self.viewer:
            return
        self.SelectCamera()
        if not hasattr(self.viewer.GUI, 'Xform'):
            return
        self.viewer.Xform.set(self.currentTransfMode)
        if self.currentTransfMode=='Object':
            self.viewer.GUI.enableNormalizeButton(Tkinter.NORMAL)
            self.viewer.GUI.enableCenterButton(Tkinter.NORMAL)
        else:
            self.viewer.GUI.enableNormalizeButton(Tkinter.DISABLED)
            self.viewer.GUI.enableCenterButton(Tkinter.DISABLED)


    def SwapBuffers(self):
        if __debug__:
            if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.tk.call(self._w, 'swapbuffers')

    def Activate(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Make this Opengl widget the current destination for drawing."""
        self.tk.call(self._w, 'makecurrent')

    def getContext(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Make this Opengl widget the current destination for drawing."""
        return self.tk.call(self._w, 'contexttag')


    def activeStereoSupport(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        try:
            root = Tkinter.Tk()
            root.withdraw()
            loadTogl(root)
            Tkinter.Widget(root, 'togl', (), {'stereo':'native'} )
            #currentcontext = self.tk.call(self._w, 'contexttag')
            #print "StandardCamera.activeStereoSupport currentcontext", currentcontext

            if GL.glGetIntegerv(GL.GL_STEREO)[0] > 0:
                lReturn = True
            else:
                #print 'Stereo buffering is not available'
                lReturn = False

            Tkinter.Widget(root, 'togl', (), {'stereo':'none'} )
            return lReturn

        except Exception, e:
            if str(e)=="Togl: couldn't get visual":
                #print 'Stereo buffering is not available:', e
                lReturn = False
            elif str(e)=="Togl: couldn't choose stereo pixel format":
                #print 'Stereo buffering is not available:', e
                lReturn = False
            else:
                print 'Togl error:',e
                lReturn = False
            Tkinter.Widget(root, 'togl', (), {'stereo':'none'} )
            return lReturn
