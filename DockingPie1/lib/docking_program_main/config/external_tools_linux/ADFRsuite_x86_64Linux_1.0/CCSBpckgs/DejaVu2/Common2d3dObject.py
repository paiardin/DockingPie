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
# Date: October 2006 Authors: Guillaume Vareille, Michel Sanner
#
#    vareille@scripps.edu
#    sanner@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Guillaume Vareille, Michel Sanner and TSRI
#
# Revision:
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/DejaVu2/Common2d3dObject.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#
# $Id: Common2d3dObject.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#
# modified by Stefano Forli 2009

import types
import numpy
import warnings

from opengltk.OpenGL import GL

from DejaVu2.viewerFns import checkKeywords

class Common2d3dObject:
    """Base Class inherited by Insert2d and Geom
"""

    keywords = [
        'immediateRendering', # set to 1 to avoid using dpyList
        'listed',    # 0/1 when set geometry appears in object list
        'name',
        'needsRedoDpyListOnResize',
        'pickable',
        'protected', # 0/1 when set geometry cannot be deleted
        'replace',
        'tagModified', # use False to avoid toggling _modified
        'visible',
        'animatable'
    ]

    _redoFlags = {
        'redrawFlag' : 1,
        'updateOwnGuiFlag' : 2,
        'redoViewerDisplayListFlag' : 4,
        'redoDisplayListFlag' : 8,
        'redoTemplateFlag' : 16,
        'redoDisplayChildrenListFlag' : 32
    }

    def __init__(self, name='Common2d3dObject', check=1, **kw):
        """Constructor
"""
        #self.newList = GL.glGenLists(1)

        self.name = name
        self.viewer = None
        self.immediateRendering = False
        self.needsRedoDpyListOnResize = False
        self.protected = False      # False == not protected against deletion
        self.replace = True     # this is used to store in a geometry the 
          # desired behavior when adding the geom to the Viewer.
          # in AddObject, if replace is not specified, geom.replace will be
          # used
        self.isExpandedInObjectList = 1 # set to 1 when children of that geom
                                # are visible in ViewerGUI object list
        
        self.parent = None      # parent object
        self.fullName = name    # object's name including parent's name
        self.children = []      # list of children
        self.listed = True         # When true the name of this object appears
                                # in the GUI's object

        self.pickable = True       # when set the object can be picked

        self.pickNum = 0
        self.hiddenInCamera = {} # key:Camera object, value dummy

        self.dpyList = None     # display list for drawing
        self.pickDpyList = None # display list for picking
        
        # FIXME: quick hack. Flag set by SetCurrentXxxx. When set object's
        # transformation and material are saved in log file
        self.hasBeenCurrent = False


        self.depthMask = 1 # 0: zbuffer is readOnly for transparent objects
                           # 1: zbuffer is read write for transparent objects
        self.srcBlendFunc = GL.GL_SRC_ALPHA
        #self.dstBlendFunc = GL.GL_ONE #GL.GL_ONE_MINUS_SRC_COLOR
        # ONE_MINUS_SRC_ALPHA works for colorMapEditor
        self.dstBlendFunc = GL.GL_ONE_MINUS_SRC_ALPHA 

        self.transparent = False #viewerConst.NO
        
        self.visible = True

        m = numpy.reshape( numpy.identity(4), (1,16) )
        self.instanceMatricesFortran = m.astype('f')

        self.ownGui = None

        self.animatable = True

        apply( self.Set, (check, 0), kw)

        self._modified = False  # used to decide if we should save state


    def Set(self, check=1, redo=1, updateOwnGui=True, **kw):
        """set data for this object
check=1 : verify that all the keywords present can be handle by this func 
redo=1 : append self to viewer.objectsNeedingRedo
updateOwnGui=True : allow to update owngui at the end this func
"""
        if __debug__ and check:
            apply( checkKeywords, (self.name, self.keywords), kw)

        if self.viewer:
            self.viewer.currentCamera.Activate()

        tagModified = kw.pop( 'tagModified', True )
        assert tagModified in (True, False)
        self._modified = tagModified

        redoFlags = 0
        # this allows to grab directly the correct function rather than going through the full Set function
        for k, v in kw.items():
            if k == 'transparent':
                continue
            lFunctionName = '_set' + k[0].upper() + k[1:]
            if hasattr(self, lFunctionName):
                lFunction = getattr(self, lFunctionName)
                redoFlags |= lFunction(v)
                kw.pop(k)

        if len(kw) == 0:
            return self.redoNow(redo, updateOwnGui, redoFlags)

        val = kw.pop( 'name', None)
        if not val is None:
            assert type(val) == types.StringType
            if self.parent:
                # make sure name is unique amongst siblings
                sameName = None
                for o in self.parent.children:
                    if o.name==val and o!=self:
                        raise RuntimeError("name %s already used by child of %s geometry"%(val, self.parent.name))
            oldFullName = self.fullName
            self.name = val
            self.fullName = self.BuildFullName()
            if self.viewer and self.viewer.GUI:
                # rename object in tree viewer
                self.viewer.GUI.renameObject(self, val)
                #self.viewer.GUI.occlusionCamera.chooser.remove(oldFullName)
                #self.viewer.GUI.occlusionCamera.chooser.add( (self.fullName, 'no comment available', self) )

        val = kw.pop( 'protected', None)
        if val is not None:
            if val in [True, False, 0, 1]:
                allObjects = self.AllObjects()
                for obj in allObjects:
                    obj.protected = val
                if self.viewer and self.viewer.GUI:
                    if self.viewer.currentObject==self:
                        self.viewer.SetCurrentObject(self)

        val = kw.pop( 'listed', None)
        if val is not None:
            if val in [True, False]:
                self.listed = val

        val = kw.pop( 'visible', None)
        if  not val is None:
            if val in (True, False):
                if val != self.visible and self.viewer:
                    redoFlags |= self._redoFlags['redoViewerDisplayListFlag']
                self.visible = val
            else: raise ValueError("visible can only be True or False")

        val = kw.pop( 'pickable', None)
        if val is not None:
            if val in (True, False):
                self.pickable = val
            else: raise ValueError("pickable can only be True or False")
            
        val = kw.pop( 'immediateRendering', None)
        if val is not None:
            assert val in [False,True], "only False or True are possible"
            self.immediateRendering = val
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        val = kw.pop( 'replace', None)
        if val is not None:
            self.replace = val

        val = kw.pop( 'animatable', None)
        if val is not None:
            if val in (True, False):
                self.animatable = val
	    else: raise ValueError("visible can only be True or False")

        return self.redoNow(redo, updateOwnGui, redoFlags)


    def redoNow(self, redo, updateOwnGui, redoFlags, setCurrent=True):
        if redo and redoFlags and self.viewer:
            if (redoFlags & self._redoFlags['redoTemplateFlag']) != 0:
                self.redoTemplate()
                redoFlags |= self._redoFlags['redrawFlag']
            if (redoFlags & self._redoFlags['redoDisplayListFlag']) != 0:
                if self not in self.viewer.objectsNeedingRedo.keys():
                    self.viewer.objectsNeedingRedo[self] = None
                redoFlags |= self._redoFlags['redrawFlag']
            if (redoFlags & self._redoFlags['redoDisplayChildrenListFlag']) != 0:
                lObjectsNeedingRedo = self.viewer.objectsNeedingRedo.keys()
                for child in self.AllObjects():
                    if child not in lObjectsNeedingRedo:
                        self.viewer.objectsNeedingRedo[child] = None
                redoFlags |= self._redoFlags['redrawFlag']
            if (redoFlags & self._redoFlags['redoViewerDisplayListFlag']) != 0:
                self.viewer.deleteOpenglList()
                redoFlags |= self._redoFlags['redrawFlag']
            if updateOwnGui is True and self.ownGui is not None \
              and (redoFlags & self._redoFlags['updateOwnGuiFlag']) != 0:
                self.updateOwnGui()
            if redoFlags & self._redoFlags['redrawFlag']:
                self.viewer.Redraw()

            if setCurrent and self.viewer and self.viewer.GUI:
                if self.viewer.currentObject==self:
                    # reflect changes in DejaVu2 GUI
                    self.viewer.GUI.SetCurrentObject(self)

        return redoFlags


    def _Hide(self):
        """Transient signal handler"""
        self.visible = False
        self.viewer.deleteOpenglList()
        self.viewer.Redraw()


    def _Remove(self, a, b):
        """Transient signal handler"""
        
        self.protected=False
        self.viewer.RemoveObject(self)


    def getState(self, full = 1):
        """return a dictionary describing this object's state
This dictionary can be passed to the Set method to restore the object's state
"""
        state = {
            'immediateRendering':self.immediateRendering,
            'listed':self.listed,
            'name':self.name,
            'needsRedoDpyListOnResize': self.needsRedoDpyListOnResize,
            'pickable':self.pickable,
            'protected':self.protected,
            'replace':self.replace,
            #'tagModified': ,
            'visible':self.visible,
            'animatable': self.animatable
            }
        return state


    def getVisible(self):
        #print "getVisible", self.name
        obj = self
        while obj.visible:
            if obj.parent:
                obj = obj.parent
            else:
                break
        return obj.visible


    def deleteOpenglList(self):
        #print "Common2d3dObject.deleteOpenglList", self

        self.redoDspLst = 0

        currentcontext = self.viewer.currentCamera.getContext()
        if self.dpyList is not None:
            if currentcontext != self.dpyList[1]:
                warnings.warn("""deleteOpenglList failed because the current context is the wrong one""")
                #print "currentcontext != self.dpyList[1]", currentcontext, self.dpyList[1]
                #import traceback;traceback.print_stack()
            else:
                #print '-%d'%self.dpyList[0], currentcontext,"glDeleteLists Common2d3dObject"
                GL.glDeleteLists(self.dpyList[0], 1)
                self.dpyList = None

        if self.pickDpyList is not None:
            if currentcontext != self.pickDpyList[1]:
                warnings.warn("""deleteOpenglList failed because the current context is the wrong one""")
                #print "currentcontext != self.pickDpyList[1]", currentcontext, self.pickDpyList[1]
            else:
                #print '-%d'%self.pickDpyList[0], currentcontext, "glDeleteLists Common2d3dObject2", 
                GL.glDeleteLists(self.pickDpyList[0], 1)
                self.pickDpyList = None


    def setViewer(self, viewer, buildDisplayList):
        """set viewer in all children of geom object
"""
        #print "Common2d3dObject.setViewer", self, buildDisplayList
        #import traceback;traceback.print_stack()
        assert viewer is not None
        self.viewer = viewer
        
        # first, we make sure the template exists for this object if it uses template
        if hasattr(self, 'makeTemplate') is True:
            if self.templateDSPL is None:
                viewer.currentCamera.Activate()
                self.makeTemplate()                
        
        if buildDisplayList:
            viewer.objectsNeedingRedo[self] = None

        for o in self.children:
            # set viewer in all children of geom object
            o.setViewer(viewer, buildDisplayList)


    def getDepthMask(self):
        return self.depthMask


    def BuildFullName(self):
        """Returns the full name of current geom """
        if self.parent:
            return self.parent.fullName +'|'+self.name
        else:
            return self.name


    def LastParentBeforeRoot(self):
        """return the last parent of a node before the root object is reached
        """
        obj = self
        
        #this function can be call even if the geom doesn't belong to a viewer
        if self.viewer is None: 
            while obj.parent:
                obj = obj.parent
        else:
            while obj.parent and obj.parent!=self.viewer.rootObject:
                obj = obj.parent
        return obj


    def _visibleObjs(self, geom, result):
        if len(geom.vertexSet):
            result.append(geom)
        for child in geom.children:
            if child.visible:
                result = self._visibleObjs(child, result)
        return result

    
    def VisibleObjects(self):
        """Return a list of visible objects that contain vertices in the object
        sub tree. This function is recursive and respects visibility inheritence"""
        if not self.visible: return []
        result = []
        result = self._visibleObjs(self, result)
        return result
    

    def AllVisibleObjects(self):
        """Return a list of all visible children of this object in the object tree"""
        if not self.visible: return
        obj = [self] + filter(lambda x: x.visible, self.children)
        if len(obj) == 1: return obj
        oi = 1
        while oi != len(obj):
            obj = obj + obj[oi].children
            oi = oi + 1
        return obj
    
        
    def AllObjects(self):
        """Return a list of all children of this object in the object tree"""

        obj = [self] + self.children
        if len(obj) == 1: return obj
        oi = 1
        while oi != len(obj):
            obj = obj + obj[oi].children
            oi = oi + 1
        return obj


    def RedoDisplayList(self):
        #print "Common2d3dObject.RedoDisplayList", self

        if hasattr(self, 'dpyList') and self.dpyList is not None:
            lNewList = self.dpyList[0]
        else:
            lNewList = GL.glGenLists(1)
            self.viewer.deleteOpenglList()

        #lNewList = self.newList
        lContext = self.viewer.currentCamera.getContext()
        #print "lNewList Common2d3dObject.RedoDisplayList", lNewList, lContext, self.name
        self.dpyList = ( lNewList,
                         self.viewer.currentCamera.getContext()
                       )
        GL.glNewList(self.dpyList[0], GL.GL_COMPILE)
        #print '+%d'%self.dpyList[0], lContext, "glNewList Common2d3dObject"
        if hasattr(self, 'Draw'):
            status = self.Draw()
        else:
            status = 0
        #print '*%d'%GL.glGetIntegerv(GL.GL_LIST_INDEX), "glEndList Common2d3dObject"
        GL.glEndList()
        if status == 0:
            self.deleteOpenglList()


    def DisplayFunction(self):
        """ display function. may be re-implemented by subclass
"""
        #print "Common2d3dObject.DisplayFunction", self
        if self.dpyList is not None:
            currentcontext = self.viewer.currentCamera.getContext()
            if currentcontext != self.dpyList[1]:
                warnings.warn("""DisplayFunction failed because the current context is the wrong one""")
                #print "currentcontext != self.dpyList[1]", currentcontext, self.dpyList[1]
            else:
                #print '#%d'%self.dpyList[0], currentcontext, "glCallList Common2d3d"
                GL.glCallList(self.dpyList[0])


    def showOwnGui(self):
        #print "showOwnGui", self
        if self.ownGui is None:
            self.createOwnGui()
        if self.ownGui.winfo_ismapped() == 0:
            self.ownGui.deiconify()
        self.ownGui.lift()


    def hideOwnGui(self):
        #print "hideOwnGui", self
        if self.ownGui is not None \
          and self.ownGui.winfo_ismapped() == 1:
            self.ownGui.withdraw()
