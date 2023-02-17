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
# $Header: /mnt/raid/services/cvs/DejaVu2/Tk/Viewer.py,v 1.1.1.1.4.1 2017/07/13 22:22:07 annao Exp $
#
# $Id: Viewer.py,v 1.1.1.1.4.1 2017/07/13 22:22:07 annao Exp $
#
import Tkinter, threading

import DejaVu2
from DejaVu2.Geom import Geom
from DejaVu2.MaterialEditor import MaterialEditor
from DejaVu2.Viewer import ViewerBase, SetCurrentObjectEvent
from DejaVu2.Common2d3dObject import Common2d3dObject

from DejaVu2.Tk.ViewerGUI import ViewerGUI
from DejaVu2.Tk.Camera import Camera
from DejaVu2.Tk import loadTogl


class Viewer(Tkinter.Widget, Tkinter.Misc, ViewerBase):

    def __init__(self, master=None, nogui=0, screenName=None,
                 guiMaster=None, classCamera=None, showViewerGUI=True,
                 autoRedraw=True, verbose=True, cnf={}, **kw):
	"""Viewer's constructor
"""
        Tkinter.Widget.__init__(self, None, None)
        name = 'DejaVu2_Viewer_Tk'
        if guiMaster is not None and master is None:
            master = guiMaster
            master = master.winfo_toplevel()
            #master = Tkinter.Toplevel(master)
            self.ownMaster = True

        withdrawDefaultRoot = Tkinter._default_root is None

        if not master:
            # DISABLE DISPLAY handling as it would require creating a new
            # Tk, and this creates many problems, including that Pmw uses
            # Tkinter._default_root in many places
            #if screenName is None:
            #    screenName = os.environ.get("DISPLAY", None)
            #if screenName is not None:
            #    master = Tk(screenName=screenName)
            #    self.ownMaster = True
            #    #print 'new Tk', repr(master)
            #else:
            #    master = Toplevel()
            #    self.ownMaster = True
                #print 'new Tk', repr(master)
                
            master = Tkinter.Toplevel()
            self.ownMaster = True

            if guiMaster is None and self.ownMaster:
                guiMaster = Tkinter.Toplevel(master)
            master.title(name)
        else:
            assert isinstance(master, Tkinter.BaseWidget) or isinstance(master, Tkinter.Tk)

        if withdrawDefaultRoot:
            Tkinter._default_root.withdraw()

        # string var used to decide what the trackball is moving
        self.Xform = Tkinter.StringVar(value="Object")

        self.contourTk = Tkinter.BooleanVar()

        self.spinVar = Tkinter.IntVar()
        self.spinVar.set(DejaVu2.defaultSpinningMode)
        self.spinCallBack = None

        # Decides if the call to enable GL_LIGHTNING will be considered 
        self.OverAllLightingIsOn = Tkinter.IntVar()
        self.OverAllLightingIsOn.set(1)
        
        #print 'master', repr(master)
        #print 'guiMaster', repr(guiMaster)
        self.master = master

        loadTogl(master)

        ViewerBase.__init__(self, nogui=0, screenName=None,
                            guiMaster=None, classCamera=None, 
                            autoRedraw=True, verbose=True, cnf={}, **kw)

        self.BindTrackballToObject(self.rootObject)

        # create the material editor
        self.materialEditor = MaterialEditor(None, self)
        self.materialEditor.dismissTk.grid_forget()
        self.materialEditor.dismiss()

        # register Object that have their own OpenGL context but need to have
        # the same lighting model and lights
        for l in self.lights:
            l.applyTo = [self.materialEditor]

        self.lightModel.applyTo = [self.materialEditor]
        
        # create the ViewerGui
        if showViewerGUI:
            self.GUI = ViewerGUI(self, self.maxLights, self.maxClipP,
                                 nogui=nogui, master=guiMaster)

            #self.GUI.CameraBackgroundColor = self.CurrentCameraBackgroundColor
            #self.GUI.LightColor = self.CurrentLightColor
            #self.GUI.ClipColor = self.CurrentClipColor
            #self.GUI.ObjectFrontColor = self.CurrentObjectFrontColor
            #self.GUI.ObjectBackColor = self.CurrentObjectBackColor
            #self.GUI.LightModelColor = self.LMColor

            self.GUI.addObject(self.rootObject, None)

            self.GUI.bindResetButton( self.Reset_cb)
            self.GUI.bindNormalizeButton( self.Normalize_cb)
            self.GUI.bindCenterButton( self.Center_cb)
            self.GUI.bindDeleteButton( self.Delete_cb)
    ##         self.GUI.Exit = self.__del__

            if nogui and isinstance(self.GUI.root, Tkinter.Toplevel):
                self.GUI.withdraw()

        #self.GUI.addObject(self.pickVerticesSpheres, None)
        #self.GUI.showPickedVertex.set(self.showPickedVertex)

        if self.autoRedraw:
            self.pendingAutoRedrawID = self.master.after_idle(self.ReallyRedraw)


    def dialog(self):
#        t="Do you Wish to Quit?"
#        from SimpleDialog import SimpleDialog
#        d = SimpleDialog(self.currentCamera, text=t,
#                                     buttons=["Quit","Cancel"],
#                                     default=0, title="Quit?")
#        
#        ok=d.go()
        ok = tkMessageBox.askokcancel("Quit?","Do you Wish to Quit?")
        if ok:
            self.quit_cb()
        else:
            return


    def BindTrackballToObject(self, obj, allCameras=None):
	"""Bind trackball to the current object"""

        self.useMasterDpyList = self.oldUseMasterDpyList

	if not isinstance(obj, Common2d3dObject):
	    raise AttributeError('first parameter has to be an instance of \
Comon2d3dObject')

	self.SetCurrentObject(obj)

	if allCameras:
	    for c in self.cameras:
                if isinstance(obj, Geom):
                    c.bindAllActions('Object')
                elif isinstance(obj, Insert2d):
                    c.bindAllActions('Insert2d')
##  		c.trackball.B2motion = self.RotateCurrentObject
##  		c.trackball.B3motion = self.TranslateCurrentObjectXY
##  		c.trackball.ShiftB3motion = self.TranslateCurrentObjectZ
##  		c.trackball.ShiftB2motion = self.ScaleCurrentObject
		self.Reset_cb = self.ResetCurrentObject
		c.currentTransfMode = 'Object'
	else:
	    c = self.currentCamera
            if isinstance(obj, Geom):
                c.bindAllActions('Object')
            elif isinstance(obj, Insert2d):
                c.bindAllActions('Insert2d')
##  	    c.trackball.B2motion = self.RotateCurrentObject
##  	    c.trackball.B3motion = self.TranslateCurrentObjectXY
##  	    c.trackball.ShiftB3motion = self.TranslateCurrentObjectZ
##  	    c.trackball.ShiftB2motion = self.ScaleCurrentObject
	    self.Reset_cb = self.ResetCurrentObject
	    c.currentTransfMode = 'Object'
            
        if self.GUI:
            self.GUI.BindTrackballToObject(obj, allCameras)

        # generate an event
        self.eventHandler.dispatchEvent(SetCurrentObjectEvent(obj))


    def quit_cb(self):
        self.currentCamera.frame.master.destroy()


    def AddCamera(self, master=None, screenName=None, classCamera=None,
                  stereo='none', num=None, cnf={}, **kw):
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
##             if screenName is None:
##                 if os.environ.has_key("DISPLAY"):
##                     screenName = os.environ["DISPLAY"]
##                 if screenName is not None:
##                     master = Tk(screenName=screenName)
##                 else:
##                     master = Tk()
            master = Tkinter.Toplevel(self.master)
            master.title(name)
            cameraOwnsMaster = True
            
        #else:
        #    assert isinstance(master, Widget)

        # simulate the setting of TCLLIPATH
        # we shoudl really check is the directory where we know togl is
        # isn't already in that path
        # path = master.tk.globalgetvar('auto_path')
        # master.tk.globalsetvar('auto_path', os.path.join(sys.exec_prefix,
        #                                                 'lib')+ ' '+ path )
        # master.tk.call('package', 'require', 'Togl')


        #else:
        #    assert isinstance(master, Widget)

##         # simulate the setting of TCLLIPATH
##         # we should really check is the directory where we know togl is
##         # isn't already in that path
##         path = master.tk.globalgetvar('auto_path')
## #        master.tk.globalsetvar('auto_path', os.path.join(sys.exec_prefix,
## #                                                         'lib')+ ' '+ path )
##         import opengltk
##         toglpath = os.path.join(opengltk.__path__[0], 'OpenGL/Tk/ ')
##         master.tk.globalsetvar('auto_path', (toglpath,) + path )
##         #master.tk.globalsetvar('auto_path', toglpath)
##         #print master.tk.globalgetvar('auto_path')
##         master.tk.call('package', 'require', 'Togl')

        #kw['takefocus'] = 1
        kw['stereo'] = stereo
        c = classCamera(master, screenName, self, num, check=1, cnf=cnf, **kw)
        if hasattr(c.frame.master,"protocol"):
            c.frame.master.protocol("WM_DELETE_WINDOW",self.dialog)
        c.eventManager.AddCallback('<KeyPress>', self.modifierDown)
        c.eventManager.AddCallback('<KeyRelease>', self.modifierUp)
        c.eventManager.AddCallback('R', self.Reset_cb_arg)
        c.eventManager.AddCallback('r', self.Reset_cb_arg)
        c.eventManager.AddCallback('A', self.AutoDepthCue)
        c.eventManager.AddCallback('a', self.AutoDepthCue)
        c.eventManager.AddCallback('N', self.Normalize_cb_arg)
        c.eventManager.AddCallback('n', self.Normalize_cb_arg)
        c.eventManager.AddCallback('C', self.Center_cb_arg)
        c.eventManager.AddCallback('c', self.Center_cb_arg)
        c.eventManager.AddCallback('D', self.Depth_cb_arg)
        c.eventManager.AddCallback('d', self.Depth_cb_arg)
        c.eventManager.AddCallback('T', self.toggleTransformRootOnly)
        c.eventManager.AddCallback('t', self.toggleTransformRootOnly)
        c.eventManager.AddCallback('L', self.toggleOpenglLighting)
        c.eventManager.AddCallback('l', self.toggleOpenglLighting)
        c.eventManager.AddCallback('O', self.SSAO_cb_arg)
        c.eventManager.AddCallback('o', self.SSAO_cb_arg)    
        c.ownMaster = cameraOwnsMaster
        
        if self.GUI is not None:
            self.GUI.bindModifersToTransformInfo(master)

	self.cameras.append( c )
	if len(self.cameras) == 1:
	    self.currentCamera = c
            c.hasBeenCurrent = 1
	    #self.trackball = c.trackball
	    c.frame.config( background = "#900000" )

	# make the trackball transform the current object
        if self.rootObject:
            self.BindTrackballToObject(self.rootObject)

        c.firstRedraw = True
        c.Activate()

	return c

    
    def _DeleteCamera(self, camera):
        """Remove the given camera in the right order
"""
        #print 'Viewer._DeleteCamera ', camera
        # Check if this camera shareCTX with other camera.
        if hasattr(camera, 'shareCTXWith'): 
            while len(camera.shareCTXWith):
                cam = camera.shareCTXWith[0]
                self._DeleteCamera(cam)
        camera.destroy()
        if camera.ownMaster:
            camera.frame.master.destroy()
        else:
            camera.frame.destroy()
        self.cameras.remove(camera)
        for c in self.cameras:
            if hasattr(c, 'shareCTXWith'):
                c.shareCTXWith.remove(camera)

        
    def DeleteCamera(self, camera):
        """
        Remove the given camera from the viewer and takes care
        of the dpyList if camera is cameras[0]
        """
#        #delete NPR rendering toplevel
#        if camera.imCanvastop:
#            camera.imCanvastop.destroy()
#            camera.imCanvas = None
#            camera.imCanvas1 = None
#            camera.imCanvas2 = None
#            camera.imCanvas3 = None
            
        # Remove the camera from the list of cameras associated to
        # the viewer.
        camIndex = self.cameras.index(camera)

        # the current openGL context has been destroyed so
        # the dpyList need to be destroyed only if the CTX is not
        # shared by any other camera.
        if camIndex == 0:
            self.objectsNeedingRedo = {}
            for g in self.rootObject.AllObjects():
                g.deleteOpenglList()
                self.objectsNeedingRedo[g] = None

        self._DeleteCamera(camera)
        # If this camera is the current camera
        if self.currentCamera == camera:
            if len(self.cameras) == 0:
                # There is no more cameras then set currentCamera to None
                self.currentCamera = None
            else:
                # Set the current Camera to be the first camera of the list.
                self.currentCamera = self.cameras[0]
    
            
    def postNextRedraw(self):
        if self.autoRedraw:
            self.pendingAutoRedrawID = self.master.after(10, self.ReallyRedraw)

    def startAutoRedraw(self):
        self.autoRedraw = True
        self.pendingAutoRedrawID = self.master.after(10, self.ReallyRedraw)


    def stopAutoRedraw(self):
        if self.pendingAutoRedrawID:
            #print 'cancelling autoRedraw_after', self.pendingAutoRedrawID
            # cancelling does not seem to work (i.e. the callback gets called
            # even if we cancelit), so we force the id to None
            # and check for ID None is ReallyRedraw
            # we still have to cancel else its lows down and freezes eventually
            self.master.after_cancel(self.pendingAutoRedrawID)
            self.pendingAutoRedrawID = None
        self.autoRedraw = False
        #  This was used before I added pendingAutoRedrawID
        #from time import sleep
        #sleep(0.2)
        #self.master.update()
        

    def checkIfRedrawIsNeeded(self):
        if self.pendingAutoRedrawID:
            self.master.after_cancel(self.pendingAutoRedrawID)

        if self.suspendRedraw:
            self.master.after(100, self.ReallyRedraw)
            return False

        if not self.needsRedraw:
            self.master.after(10, self.ReallyRedraw)
            return False
                    
        if threading.currentThread().getName()!='MainThread':
            self.master.after(10, self.ReallyRedraw)
            return False

        if self.autoRedraw and self.pendingAutoRedrawID is None:
            return False


