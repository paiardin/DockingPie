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
# Date: 2000 Author: Michel F. SANNER
#
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
# $Header: /mnt/raid/services/cvs/DejaVu2/Camera.py,v 1.5.4.2 2017/11/08 23:27:21 annao Exp $
#
# $Id: Camera.py,v 1.5.4.2 2017/11/08 23:27:21 annao Exp $
#

"""Camera Module:

This Module implements the Camera class and the Fog class.
"""

## NOTE about share context and sharelist keyword arguments to Camera
##   By default ane new Camera will share its context with camaera 0
## which creates the OpenGL context
## Passing an argument sharecontext=None will make the new camera not share
## the context but still the display lists will be shared.
## Passing sharelist=None in addition to sharecontext=None will create a
## camera with a completely separate OpenGL context.
## All display list are created in the context of the Camera 0 which is the
## one activated in ReallyRedraw when display list are created.

import os, sys, warnings

from PIL import Image
from PIL import ImageFilter
from PIL import ImageChops

from opengltk.OpenGL.GLU import gluPerspective, gluPickMatrix, gluUnProject, gluErrorString, gluLookAt
from opengltk.extent import _gllib
from opengltk.OpenGL.GL import *
from opengltk.extent.utillib import glCleanRotMat
from opengltk.OpenGL import GL

import DejaVu2
from DejaVu2.Insert2d import Insert2d
from DejaVu2.Spheres import Spheres
from DejaVu2.Ellipsoids import Ellipsoids
from DejaVu2.Cylinders import Cylinders
from DejaVu2 import bitPatterns
from DejaVu2.cursors import cursorsDict
from DejaVu2.Texture import Texture

if hasattr(DejaVu2, 'enableVertexArray') is False:
    DejaVu2.enableVertexArray = False

if hasattr( DejaVu2, 'allowedAntiAliasInMotion') is False:
    DejaVu2.allowedAntiAliasInMotion = 0

if hasattr( DejaVu2, 'enableSelectionContour') is False:
    DejaVu2.enableSelectionContour = False

if hasattr( DejaVu2, 'selectionContourSize') is False:
    DejaVu2.selectionContourSize = 0

if hasattr( DejaVu2, 'selectionContourColor') is False:
    DejaVu2.selectionContourColor = (1., 0., 1., .7)

if hasattr( DejaVu2, 'selectionPatternSize') is False:
    DejaVu2.selectionPatternSize = 6

if hasattr( DejaVu2, 'enableSSAO') is False:
    DejaVu2.enableSSAO = True

sndDeriv = [ -0.125, -0.125, -0.125,
             -0.125,    1.0, -0.125,
             -0.125, -0.125, -0.125]
    
fstDeriveV1 = [-0.125,  -0.25, -0.125,
               0.0  ,    0.0,  0.0,
               0.125,   0.25,  0.125]

fstDeriveV2 = [ 0.125,   0.25,  0.125,
                0.0  ,    0.0,  0.0,
                -0.125,  -0.25, -0.125]

fstDeriveH1 = [-0.125,    0.0, 0.125,
               -0.25 ,    0.0, 0.25,
               -0.125,    0.0, 0.125]

fstDeriveH2 = [ 0.125,    0.0, -0.125,
                0.25 ,    0.0, -0.25,
                0.125,    0.0, -0.125]
    

from time import time
from opengltk.extent.utillib import namedPointsWithNames

import os, math, types, weakref
import numpy

import viewerConst, jitter
from Geom import Geom
from Transformable import Transformable
import colorTool
import viewerFns
#from Trackball import Trackball
from IndexedPolygons import IndexedPolygons
from viewerFns import checkKeywords

import array
matf = array.array('f', [0]*16)


try:
    from opengltk.extent.glextlib import *
except ImportError:
    try:
        from opengltk.extent._glextlib import *
    except Exception, e:
        print 'ERROR: failed to import opengltk.extent._glextlib'
        print e
        #print "no glActiveTextue"
        DejaVu2.enableSSAO = False

class Fog:

    keywords = [
        'tagModified',
        'enabled',
        'start',
        'end',
        'density',
        'mode',
        'color'
        ]
    

    def Reset(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.color = (0.0, 0.0, 0.0, 1.0)
        self.enabled = False
        self.start = 25
        self.end = 40
        self.mode = GL_LINEAR
        self.density = 0.1
        self._modified = False


    def __init__(self, camera):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.Reset()
        self.camera = weakref.ref(camera) # used to activate right context
                             # the alternative would be to have the Set of
                             # fog values be done trough Camera.Set
        self.name = 'Fog'+camera.name

        
    def getState(self):
        """return a dictionary describing this object's state
This dictionary can be passed to the Set method to restore the object's state
"""
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        state = {'enabled':self.enabled,
                 'start':self.start,
                 'end':self.end,
                 'density':self.density,
                 'color':self.color
                }
        mode='GL_LINEAR'
        if self.mode==GL_EXP: mode='GL_EXP'
        elif self.mode==GL_EXP2: mode='GL_EXP2'
        state['mode'] = mode
        return state

    
    def __repr__(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        return '<Fog enabled=%d from %5.2f to %5.2f mode=%d density=%f \
 color=%s>' % \
               (self.enabled, self.start, self.end, self.mode, self.density,
                repr(self.color) )


    def Set(self, check=1, **kw):
        """Set various fog parameters"""
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        if __debug__:
            if check:
                apply( checkKeywords, ("Fog",self.keywords), kw)
        
        val = kw.get( 'tagModified', True )
        assert val in [True, False]
        self._modified = val
        self.camera().Activate()

        val = kw.get( 'enabled')
        if val is not None:
            if val in [False, 0]:
                glDisable(GL_FOG)
            elif val in [True, 1]:
                glFogi(GL_FOG_MODE, self.mode)
                if self.mode == GL_LINEAR:
                    glFogf(GL_FOG_START, float(self.start) )
                    glFogf(GL_FOG_END, float(self.end) )
                else:
                    glFogf(GL_FOG_DENSITY, float(self.density))
                glFogfv(GL_FOG_COLOR, self.color)
                glEnable(GL_FOG)
            else:
                raise ValueError('Bad argument, Only True ot False are possible %s'%val)
            self.enabled = val
            
        val = kw.get( 'start')
        if not val is None:
            if kw.has_key('end'): end = kw.get('end')
            else: end = self.end
            if val < end:
                glFogf(GL_FOG_START, float(val) )
                self.start = val
            else:
                raise AttributeError('start %f has to be smaller than end=%f'%(
                                     val, self.end))

        val = kw.get( 'end')
        if not val is None:
            if val > self.start:
                glFogf(GL_FOG_END, float(val) )
                self.end = val
            else:
                raise AttributeError('end %d has to be larger than start %d'%(
                                     val, self.start))

        val = kw.get( 'density')
        if not val is None:
            if val <= 1.0 and val >= 0.0:
                glFogf(GL_FOG_DENSITY, float(val))
                self.density = val
            else:
                raise AttributeError('density has to be <=1.0 and >= 0.0')

        val = kw.get( 'mode')
        if not val is None:
            if val=='GL_LINEAR': val=GL_LINEAR
            elif val=='GL_EXP': val=GL_EXP
            elif val=='GL_EXP2': val=GL_EXP2
            if val in (GL_LINEAR, GL_EXP, GL_EXP2):
                glFogi(GL_FOG_MODE, int(val))
                self.mode = val
            else:
                raise AttributeError('mode has to be GL_LINEAR,GL_EXP or\
GL_EXP2')

        val = kw.get( 'color')
        if not val is None:
            self.color = colorTool.OneColor( val )
            glFogfv(GL_FOG_COLOR, self.color)



class PickObject:
    """Class to represent the result of picking or drag selection
the keys of hits dictionary are geometries,
the values are lists of 2-tuples (vertexIndexe, instance), where vertexInd
is the index of the a vertex of a face of the geometry and instance is a list
of integer providing the instance matrix index for the geometry and all its
parents.
"""
    
    def __init__(self, mode, camera, type='vertices'):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        assert mode in ['pick', 'drag select']
        assert type in ['vertices', 'parts']
        self.type = type
        self.mode = mode
        self.hits = {}
        self.camera = weakref.ref(camera)
        self.p1 = None # intersection of pick ray with front clip plane
        self.p2 = None # intersection of pick ray with back clip plane
        self.event = None # event that triggered picking
        self.box = (0,0,0,0) # screen coordinates of selection box
        self.operation = 'add' # can be 'remove'
        
    def add(self, object=None, vertex=None, instance=0, operation='add'):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        if self.hits.has_key(object):
            self.hits[object].append( (vertex, instance) )
        else:
            self.hits[object] = [ (vertex, instance) ]
        self.operation = operation
        
    def printInfo(self):
        print 'Pick mode:%s type:%s hits:%d'%(self.mode, self.type, len(self.hits))
        for key, values in self.hits.items():
            print key.name,
            for v in values:
                print 'vertex:%d instance:%s'%(v[0], str(v[1])),
            print
            

from math import sqrt, fabs
import numpy
from IndexedPolygons import IndexedPolygons
from Spheres import Spheres

class StandardCameraBase(Transformable):
    """Class for Opengl 3D drawing window"""

    initKeywords = [
        'height',
        'width',
        'rgba'
        'redsize',
        'greensize',
        'bluesize',
        'double',
        'depth',
        'depthsize',
        'accum',
        'accumredsize',
        'accumgreensize',
        'accumbluesize',
        'accumalphasize',
        'alpha',
        'alphasize',
        'stencil',
        'stencilsize',
        'auxbuffers',
        'privatecmap',
        'overlay',
        'stereo',
        'time',
        'sharelist',
        'sharecontext',
        'ident',
        'rootx',
        'rooty',
        'side',
        'stereoflag',
        'ssao',
        'SSAO_OPTIONS'
        ]
    
    setKeywords = [
        'tagModified',
        'height',
        'width',
        'fov',
        'near',
        'far',
        'color',
        'antialiased',
        'contours',
        'd1ramp',
        'd1scale',
        'd1off',
        'd1cutL',
        'd1cutH',
        'd2scale',
        'd2off',
        'd2cutL',
        'd2cutH',
        'boundingbox',
        'rotation',
        'translation',
        'scale',
        'pivot',
        'direction',
        'lookAt',
        'lookFrom',
        'projectionType',
        'rootx',
        'rooty',
        'stereoMode',
        'sideBySideRotAngle',
        'sideBySideTranslation',
        'TV3DTranslation',
        'TV3DRotAngle',
        'TV3DScaling',
        'suspendRedraw',
        'drawThumbnail',
        'ssao',
        'SSAO_OPTIONS'        
        ]

    PERSPECTIVE = 0
    ORTHOGRAPHIC = 1

    stereoModesList = ['MONO',
                      'SIDE_BY_SIDE_CROSS',
                      'SIDE_BY_SIDE_STRAIGHT',
                      '3DTV',
                      'STEREO_BUFFERS',
                      'COLOR_SEPARATION_RED_BLUE',
                      'COLOR_SEPARATION_BLUE_RED',
                      'COLOR_SEPARATION_RED_GREEN',
                      'COLOR_SEPARATION_GREEN_RED',
                      'COLOR_SEPARATION_RED_GREENBLUE',
                      'COLOR_SEPARATION_GREENBLUE_RED',
                      'COLOR_SEPARATION_REDGREEN_BLUE',
                      'COLOR_SEPARATION_BLUE_REDGREEN'
                     ]

    def getState(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """return a dictionary describing this object's state
This dictionary can be passed to the Set method to restore the object's state
"""
        states = {
            'height':self._height,
            'width':self._width,
            'rootx':self.rootx,
            'rooty':self.rooty,
            'fov':self.fovy,
            'near':self.near,
            'far':self.far,
            'color':self.backgroundColor,
            'antialiased':self.antiAliased,
            'boundingbox':self.drawBB,
            
            'rotation':list(self.rotation),
            'translation':list(self.translation),
            'scale':list(self.scale),
            'pivot':list(self.pivot),
            'direction':list(self.direction),
            'lookAt':list(self.lookAt),
            'lookFrom':list(self.lookFrom),

            'projectionType':self.projectionType,
            'stereoMode':self.stereoMode,
            'sideBySideRotAngle':self.sideBySideRotAngle,
            'sideBySideTranslation':self.sideBySideTranslation,
            'TV3DTranslation': self.TV3DTranslation,
            'TV3DRotAngle': self.TV3DRotAngle,
            'TV3DScaling': self.TV3DScaling,
            
            'suspendRedraw':self.suspendRedraw,
            'drawThumbnail':self.drawThumbnailFlag,

            'contours': self.contours,
            'd1ramp':list(self.d1ramp),
            'd1scale': self.d1scale,
            'd1off': self.d1off,
            'd1cutL': self.d1cutL,
            'd1cutH': self.d1cutH,

            'd2scale': self.d2scale,
            'd2off': self.d2off,
            'd2cutL': self.d2cutL,
            'd2cutH': self.d2cutH,
            'ssao':self.ssao,
            }
        if self.ssao:
            d = {}
#            for k in self.SSAO_OPTIONS:
#                d[k] = self.SSAO_OPTIONS[k][0]
            states['SSAO_OPTIONS'] = self.SSAO_OPTIONS
        return states


    def AutoDepthCue(self, nearOffset=0.0, farOffset=0.0, object=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """
        AutoDepthCue(nearOffset=0.0, farOffset=0.0)
        set fog start and end automatically using the bounding box of the
        specified object.
        if delta is the depth of the bounding box,
        start will be set to near+(nearOffset*delta)
        end will be set to farn+(farOffset*delta)
        """
        #print "StandardCamera.AutoDepthCue"
        #print 'BEFORE', self.near, self.far, self.fog.start, self.fog.end, self.nearDefault, self.farDefault
        self.fog.Set(enabled=True)
        if object is None:
            object = self.viewer.rootObject
        bb = object.ComputeBB(camera=self)
        lf = self.lookFrom
        la = self.lookAt
        v = lf-la
        from math import sqrt
        fl = frustrumlength = sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

        # compute unit vector in viewing direction
        viewDir1 = v/fl

        # project center of bb to view direction
        # and measure length from lookfrom to projected point
        bbCenter = (bb[1]+bb[0])*.5
        u = bbCenter-lf
        u1 = u/sqrt(numpy.sum(u*u))
        cosAlpha = numpy.dot(viewDir1, u1)
        v = u*cosAlpha
        l = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) # length from lookFrom to middle of BB
        
        # subtract from this distance
        halfDiag = bbCenter-bb[1]
        rad = sqrt(numpy.sum(halfDiag*halfDiag))
        d = rad#/sqrt(3)
        
        v = (l-d)*viewDir1
        start = sqrt(numpy.sum(v*v))
        v = (l+d)*viewDir1
        end =  sqrt(numpy.sum(v*v))
        #print 'AAAA1', l, d, start, end
        #print 'AAAA2', bb
        #print 'fog', self.fog.start, self.fog.end
        #print 'clip', self.near, self.far
        #far = -min(bb[0])+frustrumlength
        #near= -max(bb[1])+frustrumlength
        #delta = far-near
        #start = near+delta*nearOffset
        #end = far+delta*farOffset

        #if start < end:
        #    self.fog.Set(start=start, end=end)

        # we make the fog span the Z range of object's BB
        #fstart = l+bb[0][2]-10
        #import pdb; pdb.set_trace()
        fstart = l#+bb[0][2] # fog starts a middle of BB 
        fend = l+bb[1][2]+10
        if fend-fstart < 1.0:
            fend = fstart + 1.0
        if start < end:
            self.fog.Set(start=fstart, end=fend)

        # update camera near and far
        #self.Set(near=self.nearDefault, far=self.farDefault)
        #self.Set(near=self.fog.start*.9, far=self.fog.end*1.1)
        #self.Set(near=self.fog.start*0.9)
        #self.Set(far=self.fog.end*1.1)
        
        #self.Set(near=self.fog.start*0.9, far=self.fog.end*1.1)
        self.Set(near=0.1, far=self.fog.end*1.1)

        #print 'AAAA1', l, bb, self.fog.start, self.fog.end, self.near, self.far
        #print 'AFTER', self.nearDefault, self.farDefault
        #print 'fog', self.fog.start, self.fog.end
        #print 'clip', self.near, self.far


    def GrabFrontBufferAsArray(self, lock=True, buffer=GL.GL_FRONT):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Grabs the front buffer and returns it as Numeric array of
size camera._width*camera._height*3"""
        from opengltk.extent import _gllib as gllib
        width = self._width
        height = self._height
        nar = numpy.zeros(3*width*height, numpy.uint8)

        # get the redraw lock to prevent viewer from swapping buffers
        if lock:
            self.viewer.redrawLock.acquire()
            self.Activate()
        glPixelStorei(GL.GL_PACK_ALIGNMENT, 1)

        current_buffer = int(GL.glGetIntegerv(GL.GL_DRAW_BUFFER)[0])
        # tell OpenGL we want to read pixels from front buffer
	glReadBuffer(buffer)

        glFinish() #was glFlush()
        gllib.glReadPixels( 0, 0, width, height, GL.GL_RGB,
                            GL.GL_UNSIGNED_BYTE, nar)
        glFinish()
        # restore buffer from on which we operate
        glReadBuffer(current_buffer)
        if lock:
            self.viewer.redrawLock.release()
        return nar
    
        
    def GrabFrontBuffer(self, lock=True, buffer=GL.GL_FRONT):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Grabs the detph buffer and returns it as PIL P image"""
        nar = self.GrabFrontBufferAsArray(lock, buffer)

        from PIL import Image
        image = Image.fromstring('RGB', (self._width, self._height),
                                 nar.tostring())
        #if sys.platform!='win32':
        image = image.transpose(Image.FLIP_TOP_BOTTOM)
        return image

    
    def GrabZBufferAsArray(self, lock=True):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Grabs the detph buffer and returns it as a numpy array"""

        from opengltk.extent import _gllib as gllib
        width = self._width
        height = self._height
        nar = numpy.zeros(width*height, 'f')

        if lock:
            # get the redraw lock to prevent viewer from swapping buffers
            self.viewer.redrawLock.acquire()

        glFinish() #was glFlush()
        gllib.glReadPixels( 0, 0, width, height, GL.GL_DEPTH_COMPONENT,
                            GL.GL_FLOAT, nar)
        glFinish()
        if lock:
            self.viewer.redrawLock.release()

        return nar

        
    def GrabZBuffer(self, lock=True, flipTopBottom=True, zmin=None, zmax=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Grabs the detph buffer and returns it as PIL P image"""
        deptharray = self.GrabZBufferAsArray(lock)

        # map z values to unsigned byte
        if zmin is None:
            zmin = min(deptharray)
        if zmax is None:
            zmax = max(deptharray)
        if (zmax!=zmin):
            zval1 = 255 * ((deptharray-zmin) / (zmax-zmin))
        else:
            zval1 = numpy.ones(self._width*self._height, 'f')*zmax*255

        from PIL import Image

        depthImage = Image.fromstring('L', (self._width, self._height),
                                      zval1.astype('B').tostring())

        #if sys.platform!='win32':
        if flipTopBottom is True:
            depthImage = depthImage.transpose(Image.FLIP_TOP_BOTTOM)

        return depthImage
    

    def SaveImage(self, filename, transparentBackground=False):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """None <- cam.SaveImage(filename, transparentBackground=False)
The file format is defined by the filename extension.
Transparent background is only supported in 'png' format.
"""
        im = self.GrabFrontBuffer()
        if transparentBackground:
            errmsg = 'WARNING: transparent background is only supported with the png file format'
            name, ext = os.path.splitext(filename)
            if ext.lower != '.png':
                print errmsg
                filename += '.png'
            def BinaryImage(x):
                if __debug__:
                 if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
                if x==255:
                    return 0
                else:
                    return 255
            # grab z buffer
            z = self.GrabZBuffer()
            # turn 255 (i.e. bg into 0 opacity and everything else into 255)
            alpha = Image.eval(z, BinaryImage)
            im = Image.merge('RGBA', im.split()+(alpha,))
        extension = os.path.splitext(filename)[1]
        kw = {}
        if not extension:
            filename += '.png'
        elif extension in ['.jpg', '.jpeg', '.JPG', '.JPEG']:
            kw = {'quality':95}

        im.save(filename, **kw)


    def ResetTransformation(self, redo=1):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        Transformable.ResetTransformation(self, redo=redo)

        # Point at which we are looking
        self.lookAt = numpy.array([0.0, 0.0, 0.0])

        # Point at which the camera is
        self.lookFrom = numpy.array([0.0, 0.0, 180.0])

        # Vector from lookFrom to lookAt
        self.direction = self.lookAt - self.lookFrom

        # Field of view in y direction
        self.fovyNeutral = 40.
        self.fovy = self.fovyNeutral

        self.projectionType = self.PERSPECTIVE
        self.left = 0.
        self.right = 0.
        self.top = 0.
        self.bottom = 0.

        # Position of clipping planes.
        self.nearDefault = .1
        self.near = self.nearDefault
        self.near_real = self.near
        self.farDefault = 50.
        self.far = self.farDefault
        self.zclipWidth = self.far - self.near

        self.SetupProjectionMatrix()

    def ResetDepthcueing(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()    
        self.fog.Set(start = 25, end = 40)

    def BuildTransformation(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Creates the camera's transformation"""

        fx, fy, fz = self.lookFrom
        ax, ay, az = self.lookAt
        #ux, uy, uz = self.up
        gluLookAt( float(fx), float(fy), float(fz),
                   float(ax), float(ay), float(az),
                   float(0), float(1), float(0))
        #print 'LOOK', self.lookAt, self.lookFrom
        #lMatrix = numpy.array(glGetDoublev(GL_MODELVIEW_MATRIX)).astype('f')
        #print 'Before', lMatrix
        #lRotation, lTranslation, lScale = self.Decompose4x4(lMatrix, cleanup=False)
        #print lRotation
        #print lTranslation
        #glLoadIdentity()
        #glTranslatef(float(lTranslation[0]), float(lTranslation[1]), float(lTranslation[2]))
        #glMultMatrixf(lRotation)
        #lMatrix = numpy.array(glGetDoublev(GL_PROJECTION_MATRIX)).astype('f')
        #print 'After', lMatrix
        
                   #float(ux), float(uy), float(uz) )
## ##         eye = self.lookFrom
## ##         center = self.lookAt
## ##         gluLookAt( eye[0], eye[1], eye[2], center[0], center[1], center[2],
## ##                    0, 1, 0)
## ##         return
## ##          rot = numpy.reshape(self.rotation, (4,4))
## ##          dir = numpy.dot(self.direction, rot[:3,:3])
## ##          glTranslatef(dir[0], dir[1], dir[2])
##         glTranslatef(float(self.direction[0]),float(self.direction[1]),float(self.direction[2]))
##         glMultMatrixf(self.rotation)
##         glTranslatef(float(-self.lookAt[0]),float(-self.lookAt[1]),float(-self.lookAt[2]))
## ##          glTranslatef(self.pivot[0],self.pivot[1],self.pivot[2])
## ##          glMultMatrixf(self.rotation)
## ##          glTranslatef(-self.pivot[0],-self.pivot[1],-self.pivot[2])


    def GetMatrix(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Returns the matrix that is used to transform the whole scene to the
        proper camera view"""

        glPushMatrix()
        glLoadIdentity()
        self.BuildTransformation()
        m = numpy.array(glGetDoublev(GL_MODELVIEW_MATRIX)).astype('f')
        glPopMatrix()
        return numpy.transpose(m)

    def GetMatrixInverse(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Returns the inverse of the matrix used to transform the whole scene
        to the proper camera view"""

        m = self.GetMatrix()
        m = numpy.reshape(numpy.transpose(m), (16,)).astype('f')
        rot, transl, scale = self.Decompose4x4(m)
        sc = numpy.concatenate((numpy.reshape(scale,(3,1)), [[1]]))
        n = numpy.reshape(rot, (4,4))/sc
        tr = numpy.dot(n, (transl[0], transl[1], transl[2],1) )
        n[:3,3] = -tr[:3]
        return n        

##      def ConcatLookAtRot(self, matrix):
##          """Rotates the lookAt point around the lookFrom point."""
##          matrix = numpy.transpose(numpy.reshape( matrix, (4,4) ))

##          rot = numpy.reshape( self.rotation, (4,4))
##          m = numpy.dot(rot, matrix)
##          m = numpy.dot(m, numpy.transpose(rot))
        
##          dir = numpy.dot(self.direction, m[:3,:3])
##          self.lookAt = self.lookFrom + dir
##          self.ConcatRotation(numpy.reshape(matrix, (16,)))
##          self.pivot = self.lookFrom

##      def ConcatLookAtRot(self, matrix):
##          """Rotates the lookAt point around the lookFrom point."""
##          self.SetPivot(self.lookFrom)
##          #print "ConcatLookAtRot", self.pivot
##          self.ConcatRotation(matrix)

##      def ConcatLookFromRot(self, matrix):
##          """Rotates the lookFrom point around the lookAt point."""
##          self.SetPivot(self.lookAt)
##          #print "ConcatLookFromRot", self.pivot
##          self.ConcatRotation(matrix)

        
##      # not implemented
##      def ConcatLookAtTrans(self, trans):
##          pass

    def ConcatRotation(self, matrix):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Rotates the camera around the lookAt point"""

        self._modified = True
        x,y,z = self.lookFrom

        import numpy
        matrix = numpy.reshape( matrix, (4,4) )
        newFrom = numpy.dot((x,y,z,1), matrix )
        self.Set(lookFrom=newFrom[:3])
        
##         glPushMatrix()
##         glLoadIdentity()
##         matrix = numpy.reshape( matrix, (4,4) )

## ##          rot = numpy.reshape( self.rotation, (4,4))
## ##          m = numpy.dot(rot, matrix)
## ##          m = numpy.dot(m, numpy.transpose(rot))
## ##          self.direction = numpy.dot(self.direction, m[:3,:3])
## ##          self.lookFrom = self.lookAt - self.direction
        
##         matrix = numpy.reshape( numpy.transpose( matrix ), (16,) )
##         glMultMatrixf(matrix)
##         glMultMatrixf(self.rotation)
##         self.rotation = numpy.array(glGetDoublev(GL_MODELVIEW_MATRIX)).astype('f')
##         self.rotation = glCleanRotMat(self.rotation).astype('f')
##         self.rotation.shape = (16,)
##         glPopMatrix()
##         #self.pivot = self.lookAt

        self.viewer.deleteOpenglList() # needed to redraw the clippingplanes little frame

    def ConcatTranslation(self, trans, redo=1):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Translate the camera from lookFrom to looAt"""

        self._modified = True
        newFrom = self.lookFrom+trans
        newAt = self.lookAt+trans
        # FIXME to update near far and fog we should compute deltaZ using view direction rather than assuming Z 
        self.Set(lookFrom = newFrom, lookAt=newAt,
                 near=self.near_real+trans[2], far= max(self.nearDefault+.1, self.far+trans[2]))
        self.fog.Set(start=self.fog.start+trans[2], end=self.fog.end+trans[2])

        #print 'GOGO123'
        #self.AutoDepthCue()
        
##         trans = numpy.array(trans)
##         sign = numpy.add.reduce(trans)
##         if sign > 0.0:
##             n = 1 + (math.sqrt(numpy.add.reduce(trans*trans)) * 0.01 )
##             newdir = self.direction*n
##             diff = self.direction-newdir
##             diff = math.sqrt(numpy.add.reduce(diff*diff))
##         else:
##             n = 1 - (math.sqrt(numpy.add.reduce(trans*trans)) * 0.01 )
##             newdir = self.direction*n
##             diff = self.direction-newdir        
##             diff = -math.sqrt(numpy.add.reduce(diff*diff))
##         self.direction = newdir
## #        print self.lookFrom, self.near, self.far, self.fog.start, self.fog.end

##         # update near and far
##         near = self.near_real + diff
##         far = self.far + diff
##         if near < far:
##             self.Set(near=near, far=far, redo=redo)

##         # update fog start and end
##         #self.fog.Set(start=self.fog.start, end=self.fog.end + diff)
##         if self.fog.start < self.fog.end + diff:
##             self.fog.Set(end=self.fog.end + diff)

##         # update viewerGUI
##         if self == self.viewer.currentCamera:
##             self.viewer.GUI.NearFarFog.Set(self.near, self.far,
##                                            self.fog.start, self.fog.end)

##         self.lookFrom = self.lookAt - self.direction
        #self.viewer.Redraw()

        self.viewer.deleteOpenglList() # needed to redraw the clippingplanes little frame

    def ConcatScale(self, scale, redo=1):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Open and close camera's FOV
"""
        #print "Camera.ConcatScale", scale
        self._modified = True
        if scale > 1.0 or self.scale[0] > 0.001:
            #fovy = abs(math.atan( math.tan(self.fovy*math.pi/180.) * scale ) * 180.0/math.pi)
            fovy = self.fovy * scale
            if fovy < 180.:
                self.Set(fov=fovy, redo=redo)

    def __repr__(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        return self.name

    def preventIntelBug_BlackTriangles(self):
        if   (DejaVu2.preventIntelBug_BlackTriangles is None) \
          or (DejaVu2.preventIntelBug_WhiteTriangles is None):
            isIntelOpenGL = GL.glGetString(GL.GL_VENDOR).find('Intel') >= 0
            if isIntelOpenGL is True:
                isIntelGmaRenderer = GL.glGetString(GL.GL_RENDERER).find("GMA") >= 0
                lPreventIntelBugs = isIntelOpenGL and not isIntelGmaRenderer
            else:
                lPreventIntelBugs = False
            if DejaVu2.preventIntelBug_BlackTriangles is None:
                DejaVu2.preventIntelBug_BlackTriangles = lPreventIntelBugs
            if DejaVu2.preventIntelBug_WhiteTriangles is None:
                DejaVu2.preventIntelBug_WhiteTriangles = lPreventIntelBugs

    def ToggleDepth(self):
        """Toggles depthcueing"""
        self.fog.Set(enabled=not self.fog.enabled)
        self.viewer.Redraw()

    def __init__(self, master, screenName, viewer, num, **kw):

        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        #print "StandardCamera.__init__"

        self.master = master
        self.notSsaoNear = False

        # used for unprojection un gluUnProject
        self.unProj_model = None
        self.unProj_proj = None
        self.unProj_view = None
        
        self.swap = True # set to false to prevent camera from swapping after
                         # redraw (used by tile renderer)
        self.suspendRedraw = False # set to True to prevent camera to be redraw
                                   # Use for ARViewer ( Added by AG 01/12/2006)
        self.visible = True
        self.lastBackgroundColorInPhotoMode = (0.,0.,0.,1.)
        
        self.initialized = 0
        Transformable.__init__(self, viewer)

        # FIXME: quick hack. Flag set by SetCurrentXxxx. When set object's
        # transformation and material are saved in log file
        self.hasBeenCurrent = 0
        self._modified = False

        #self.posLog = []

        self.ResetTransformation(redo=0)
        self.currentTransfMode = 'Object' # or 'Clip', 'Camera'

        self.renderMode = GL_RENDER
        self.pickNum = 0
        self.objPick = []

        self.drawBB = 0

        self.drawMode = None  # bit 1: for Opaque objects with no dpyList
                              # bit 2: for objects using dpy list
                              # bit 3: for Transp Object with dpyList
                              
        self.drawTransparentObjects = 0
        self.hasTransparentObjects = 0
        self.selectDragRect = 0
        self.fillSelectionBox = 0

        # variable used for stereographic display
        self.sideBySideRotAngle = 3.
        self.sideBySideTranslation = 0.
        self.imageRendered = 'MONO' # can be LEFT_EYE or RIGHT_EYE
        self.stereoMode = 'MONO' # or 'SIDE_BY_SIDE'

        self.TV3DTranslation = 0.0
        self.TV3DRotAngle = 1.5
        self.TV3DScaling = (0.5, 0.5, 0.25)
        
        self.backgroundColor = (.0,.0,.0, 1.0)

        self.selectionColor = (1.0, 1.0, 0.0, 1.0)
        self.fillDelay = 200 # delay in miliseconds for filling selection box

        ## if os.name == 'nt':
        ##     if kw['stereo'] == 'none':
        ##         kw['stereo'] = 0
        ##     else:
        ##         kw['stereo'] = 1

        #currentcontext = self.getContext()
        #print "StandardCamera.__init__ currentcontext", currentcontext

        if DejaVu2.defaultAntiAlias is None:
            if os.name == 'nt':
                DejaVu2.defaultAntiAlias = 0
            else:
                if sys.platform=='darwin':
                    # disable accumulation buffers vecause they do not resize at the widget resizes :(
                    DejaVu2.defaultAntiAlias = 0
                else:
                    try:
                        from opengltk.extent import _glextlib # tells us the graphic card is not too bad
                        DejaVu2.defaultAntiAlias = 4
                    except:
                        DejaVu2.defaultAntiAlias = 0

        #self.antiAliased = DejaVu2.defaultAntiAlias
        self.antiAliased = 0
        self._wasAntiAliased = 0
        self.drawThumbnailFlag = False
        if self.antiAliased == 0:
            self.accumWeigth = 1.
            self.jitter = None
        else:
            self.accumWeigth = 1./self.antiAliased
            self.jitter = eval('jitter._jitter'+str(self.antiAliased))

        self.newList = GL.glGenLists(1)
        self.dpyList = None

        self.visible = 1
        self.ownMaster = False # set tot true is the master of self.Frame
                               # has to be destroyed when Camera is deleted
        self.exposeEvent = False # set to true on expose events and reset in
                                 # Viewer.ReallyRedraw

        self.firstRedraw = True

        # QT
        ## # create a TK-event manager for this camera
        ## self.eventManager = EventManager(self)

        self.onButtonUpCBlist = []    # list of functions to be called when
        self.onButtonDownCBlist = []  # mouse button is pressed or depressed
                                      # they take 1 argument of type Tk event

        # register funtion to swi
        self.addButtonDownCB(self.suspendAA)
        self.addButtonUpCB(self.restoreAA)

        #self.addButtonDownCB(self.suspendNPR)
        #self.addButtonUpCB(self.restoreNPR)

        self.addButtonDownCB(self.reduceSSAOQuality)
        self.addButtonUpCB(self.restoreSSAOQuality)

        self.addButtonUpCB(self.capClippedGeoms)

        # these are used in bindPickingToMouseButton to bind picking to a
        # given mouse button
        self.mouseButtonModifiers = ['None', 'Shift', 'Control', 'Alt',
                                     'Meta']
        # this keys if self.mouseButtonActions are the objects to which the
        # trackball can be attached i.e. 'Object', 'Insert2d', 'Camera' etc.
        # for each such key there is a dict with buttonnum as key (1,2,3)
        # for each such key there is a dict with keys modifiers and values
        # a string describing an action
        #
        # e.g mouseButtonActions['Object'][3]['Shift'] = 'Ztranslation'
        self.mouseButtonActions = {}
        
        for bindings in ['Object', 'Insert2d', 'Camera', 'Clip', 'Light',
                         'Texture', 'Scissor']:
            self.mouseButtonActions[bindings] = { 1:{}, 2:{}, 3:{} }
            bd = self.mouseButtonActions[bindings]
            for b in (1,2,3):
                d = bd[b]
                for mod in self.mouseButtonModifiers:
                    d[mod] = 'None'

        # initialize actions for object
        bd = self.mouseButtonActions['Object']
        for mod in self.mouseButtonModifiers:
            bd[2][mod] = 'picking'

        bd[1]['None'] = 'rotation'
        bd[1]['Shift'] = 'addToSelection'
        bd[1]['Control'] = 'removeFromSelection'
        bd[2]['None'] = 'None'
        #bd[2]['Alt'] = 'camZtranslation'
        bd[2]['Control'] = 'zoom'
        bd[2]['Shift'] = 'scale'
        bd[3]['None'] = 'XYtranslation'
        bd[3]['Control'] = 'pivotOnPixel'
        bd[3]['Shift'] = 'Ztranslation'

        # initialize actions for Insert2d
        bd = self.mouseButtonActions['Insert2d']
        for mod in self.mouseButtonModifiers:
            bd[1][mod] = 'picking'
        #bd[2]['Alt'] = 'camZtranslation'
        bd[2]['Shift'] = 'zoom'

        # initialize actions for Clip
        bd = self.mouseButtonActions['Clip']
        bd[1]['None'] = 'rotation'
        #bd[2]['None'] = 'rotation'
        #bd[2]['Alt'] = 'camZtranslation'
        bd[2]['Control'] = 'scale'
        bd[2]['Shift'] = 'zoom'
        bd[3]['None'] = 'screenXYtranslation'
        bd[3]['Control'] = 'screenZtranslation'
        bd[3]['Shift'] = 'Ztranslation'

        # initialize actions for Light
        bd = self.mouseButtonActions['Light']
        bd[1]['None'] = 'rotation'
        #bd[2]['None'] = 'rotation'
        #bd[2]['Alt'] = 'camZtranslation'
        bd[2]['Shift'] = 'zoom'

        # initialize actions for Camera
        bd = self.mouseButtonActions['Camera']
        bd[1]['None'] = 'camRotation'
        #bd[2]['None'] = 'camRotation'
        #bd[2]['Alt'] = 'camZtranslation'
        #bd[2]['Control'] = 'zoom'
        bd[3]['Shift'] = 'zoom '
        bd[3]['None'] = 'camXYtranslation'
        bd[3]['Control'] = 'camZtranslation'

        # initialize actions for Texture
        bd = self.mouseButtonActions['Texture']
        bd[1]['None'] = 'rotation'
        #bd[2]['None'] = 'rotation'
        #bd[2]['Alt'] = 'camZtranslation'
        bd[2]['Control'] = 'scale'
        bd[2]['Shift'] = 'zoom'
        bd[3]['None'] = 'XYtranslation'
        bd[3]['Shift'] = 'Ztranslation'

        # initialize actions for Scissor
        bd = self.mouseButtonActions['Scissor']
        #bd[2]['Alt'] = 'camZtranslation'
        bd[2]['Control'] = 'scale'
        bd[2]['Shift'] = 'zoom'
        bd[3]['None'] = 'translation'
        bd[3]['Shift'] = 'ratio'

        # define actionName and callback fucntions equivalence
        self.actions = {
            'Object': {
                'picking':self.initSelectionRectangle,
                'addToSelection':self.initSelectionRectangle,
                'removeFromSelection':self.initSelectionRectangle,
                'rotation':viewer.RotateCurrentObject,
                'scale':viewer.ScaleCurrentObject,
                'XYtranslation':viewer.TranslateCurrentObjectXY,
                'Ztranslation':viewer.TranslateCurrentObjectZ,
                'zoom':viewer.ScaleCurrentCamera,
                'camZtranslation':viewer.TranslateCurrentCamera,                
                'pivotOnPixel':viewer.pivotOnPixel,                
                },
            'Insert2d': {
                'picking':self.SetInsert2dPicking,
                'zoom':viewer.ScaleCurrentCamera,
                'camZtranslation':viewer.TranslateCurrentCamera,                
                },
            'Clip': {
                'picking':None,
                'rotation':viewer.RotateCurrentClipPlane,
                'scale':viewer.ScaleCurrentClipPlane,
                'screenXYtranslation':viewer.screenTranslateCurrentObjectXY,
                'Ztranslation':viewer.TranslateCurrentObjectZ,
                'zoom':viewer.ScaleCurrentCamera,
                'camZtranslation':viewer.TranslateCurrentCamera,
                'screenZtranslation':viewer.screenTranslateCurrentObjectZ,
                },
            'Light': {
                'picking':None,
                'rotation':viewer.RotateCurrentDLight,
                'zoom':viewer.ScaleCurrentCamera,
                'camZtranslation':viewer.TranslateCurrentCamera,                
                },
            'Camera': {
                'picking':None,
                'camRotation':viewer.RotateCurrentCamera,
                'zoom':viewer.ScaleCurrentCamera,
                'zoom ':viewer.ScaleCurrentCamera, #it doesn't work if we use the same name
                'camZtranslation':viewer.TranslateCurrentCamera,
                'camXYtranslation':viewer.TranslateXYCurrentCamera, #it doesn't work if we use the same name
                },
            'Texture': {
                'picking':None,
                'rotation':viewer.RotateCurrentTexture,
                'scale':viewer.ScaleCurrentTexture,
                'XYtranslation':viewer.TranslateCurrentTextureXY,
                'Ztranslation':viewer.TranslateCurrentTextureZ,
                'zoom':viewer.ScaleCurrentCamera,
                'camZtranslation':viewer.TranslateCurrentCamera,                
                },
            'Scissor': {
                'picking':None,
                'ratio':viewer.AspectRatioScissor,
                'scale':viewer.ScaleCurrentScissor,
                'translation':viewer.TranslateCurrentScissor,
                'zoom':viewer.ScaleCurrentCamera,
                'camZtranslation':viewer.TranslateCurrentCamera,                
               },
            }

# light model is define in Viewer for all cameras
#        self.lightModel = None
        self.fog = Fog(self)
        self.fog.Set(color=self.backgroundColor)

        # attributes used to draw black outlines
#        self.imCanvastop = None   # top level window for silhouette rendering
        self.contouredImage = None # will hole the final PIL image
        self.outlineim = None # will hole the final contour
        self.outline = None  # accumulation buffer used for AA contour
        self.contours = False # set to True to enable contouring
        self.d1scale = 0.013
        self.d1off = 4
        self.d1cutL = 0
        self.d1cutH = 60

        self.d1ramp = numpy.arange(0,256,1,'f')
        
        self.d2scale = 0.0  # turn off second derivative
        self.d2off = 1
        self.d2cutL = 150
        self.d2cutH = 255

        self._suspendNPR = False # used during motion
         
        # we save it so it can be reapplied if the projection matrix is recreated
        self.pickMatrix = None

        # prepare contour highlight
        self.contourTextureName = int(glGenTextures(1)[0])
        glPrioritizeTextures(numpy.array([self.contourTextureName]), numpy.array([1.])) # supposedly make this texture fast
        _gllib.glBindTexture(GL_TEXTURE_2D, int(self.contourTextureName) )
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)

        self.stencilTextureName = int(glGenTextures(1)[0])
        # glPrioritizeTextures supposedly make this texture fast
        glPrioritizeTextures(numpy.array([self.stencilTextureName]),
                             numpy.array([1.]))
        _gllib.glBindTexture(GL_TEXTURE_2D, int(self.stencilTextureName) )
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)

        # prepare patterned highlight
        self.textureName = int(glGenTextures(1)[0])
        lTexture = numpy.array(
        (
          (.0,.0,.0,.8),(.0,.0,.0,.8),(.0,.0,.0,.8),(.0,.0,.0,.8),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),
          (.0,.0,.0,.8),(1.,1.,1.,.6),(1.,1.,1.,.6),(.0,.0,.0,.8),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),
          (.0,.0,.0,.8),(1.,1.,1.,.6),(1.,1.,1.,.6),(.0,.0,.0,.8),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),
          (.0,.0,.0,.8),(.0,.0,.0,.8),(.0,.0,.0,.8),(.0,.0,.0,.8),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),
          (.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),
          (.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),
          (.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),
          (.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),(.0,.0,.0,.0),
        ),'f')

        _gllib.glBindTexture(GL_TEXTURE_2D, self.textureName )
        glPrioritizeTextures(numpy.array([self.textureName]), numpy.array([1.])) # supposedly make this texture fast
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
        _gllib.glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 8, 8, 0, GL_RGBA,
                            GL.GL_FLOAT, lTexture)
        #glEnable(GL_TEXTURE_GEN_S)
        #glEnable(GL_TEXTURE_GEN_T)
        #glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_EYE_LINEAR )
        #glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_EYE_LINEAR )
        #glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_BLEND)

        # to protect our texture, we bind the default texture while we don't use it


                             
        self.dsT = 0
        # to protect our texture, we bind the default texture while we don't use it
        self.ssao_method = 0
        self.copy_depth_ssao = True
        if 'ssao' not in kw: 
            kw['ssao'] = False
            self.ssao = False
            self.use_mask_depth = 0
        else :
            self.ssao = kw['ssao']
            
        if 'SSAO_OPTIONS' in kw :
            self.SSAO_OPTIONS = kw['SSAO_OPTIONS']
        else :
            self.use_mask_depth = 0
            self.setDefaultSSAO_OPTIONS()            
        
        self.setShaderSSAO()  
        _gllib.glBindTexture(GL_TEXTURE_2D, 0 ) 
        self.setShaderSelectionContour()

        self.cursorStack = []
        self.currentCursor = 'default'
        

    def addButtonDownCB(self, func):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        assert callable(func)
        if func not in self.onButtonDownCBlist:
            self.onButtonDownCBlist.append(func)


    def delButtonDownCB(self, func, silent=False):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        if func in self.onButtonDownCBlist:
            self.onButtonDownCBlist.remove(func)
        else:
            if not silent:
                print 'WARNING: delButtonDownCB: function not found', func


    def addButtonUpCB(self, func):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        assert callable(func)
        if func not in self.onButtonUpCBlist:
            self.onButtonUpCBlist.append(func)


    def delButtonUpCB(self, func, silent=False):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        if func in self.onButtonUpCBlist:
            self.onButtonUpCBlist.remove(func)
        else:
            if not silent:
                print 'WARNING: delButtonUpCB: function not found', func

    def reduceSSAOQuality(self, event):
        if self.ssao:
            self._oldSSAOSamples = self.SSAO_OPTIONS['samples'][0]
            if len(self.SSAO_OPTIONS['samples']) == 5 :
                self.SSAO_OPTIONS['samples'][-1].set(2)
            else:
                self.SSAO_OPTIONS['samples'][0] = 2
            
    def restoreSSAOQuality(self, event):
        if self.ssao:
            if len(self.SSAO_OPTIONS['samples']) == 5 :
                self.SSAO_OPTIONS['samples'][-1].set(self._oldSSAOSamples)
            else:
                self.SSAO_OPTIONS['samples'][0] = self._oldSSAOSamples
            del self._oldSSAOSamples

            
    def suspendAA(self, event):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Function used to turn off anti-aliasing during motion
"""
        if self.antiAliased == False:
            self.antiAliased = 0
        assert isinstance(self.antiAliased, types.IntType), self.antiAliased
        if self._wasAntiAliased == 0 : # else it is already suspended
            if self.antiAliased <= DejaVu2.allowedAntiAliasInMotion:
                # save the state
                self._wasAntiAliased = self.antiAliased
                # we don't suspend anti alias or anti alias is already suspended
            else:
                # save the state
                self._wasAntiAliased = self.antiAliased
                # we suspend anti alias
                self.antiAliased = 0


    def restoreAA(self, event):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Function used to restore anti-aliasing after motion
"""
        #print "restoreAA"
        self.antiAliased = self._wasAntiAliased
        self._wasAntiAliased = 0


    def capClippedGeoms(self, event):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Function used to cap clipped geoms when button is released"""

        if self.currentTransfMode!='Clip':
            return
        for geom, cp, capg in self.viewer.cappedGeoms:
            if cp != self.viewer.currentClip:
                continue
            vc, fc = cp.getCapMesh(geom)
            capg.Set(vertices=vc, faces=fc)


    def suspendNPR(self, event):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Function used to turn off NPR during motion
"""
        if self.viewer.contourTk.get() is True:
            self._suspendNPR = True
            self._suspendNPR_rootMaterialsPropDiffuse = self.viewer.rootObject.materials[GL.GL_FRONT].prop[1]
            self._suspendNPR_rootMaterialsbindDiffuse = self.viewer.rootObject.materials[GL.GL_FRONT].binding[1]
            #print "self._suspendNPR_rootMaterialsPropDiffuse", self._suspendNPR_rootMaterialsPropDiffuse
            #print "self._suspendNPR_rootMaterialsbindDiffuse", self._suspendNPR_rootMaterialsbindDiffuse
            if    self._suspendNPR_rootMaterialsbindDiffuse == 1 \
              and len(self._suspendNPR_rootMaterialsPropDiffuse) == 1 \
              and len(self._suspendNPR_rootMaterialsPropDiffuse[0]) >= 3 \
              and self._suspendNPR_rootMaterialsPropDiffuse[0][0] == 1 \
              and self._suspendNPR_rootMaterialsPropDiffuse[0][1] == 1 \
              and self._suspendNPR_rootMaterialsPropDiffuse[0][2] == 1 :
                # root color is set to grey
                self.viewer.rootObject.materials[GL.GL_FRONT].prop[1] = numpy.array( ((.5, .5, .5, 1.0 ), ), 'f' )
                self.viewer.rootObject.materials[GL.GL_FRONT].binding[1] = viewerConst.OVERALL
            # lighting is turned on
            self._suspendNPR_OverAllLightingIsOn = self.viewer.OverAllLightingIsOn.get()
            self.viewer.OverAllLightingIsOn = True
            self.viewer.deleteOpenglListAndCallRedrawAndCallDisableGlLighting()


    def restoreNPR(self, event):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Function used to restore NPR after motion
"""
        if self.viewer.contourTk.get() is True:
            # root color is set to what it was
            self.viewer.rootObject.materials[GL.GL_FRONT].prop[1] = self._suspendNPR_rootMaterialsPropDiffuse
            self.viewer.rootObject.materials[GL.GL_FRONT].binding[1] = self._suspendNPR_rootMaterialsbindDiffuse
            # lighting is set to what it was
            self.viewer.OverAllLightingIsOn = self._suspendNPR_OverAllLightingIsOn
            self.viewer.deleteOpenglListAndCallRedrawAndCallDisableGlLighting()
            self._suspendNPR = False


    def pushCursor(self, cursorName):
        self.configure(cursor=cursorsDict[cursorName])
        self.cursorStack.append(self.currentCursor)
        self.currentCursor = cursorName


    def popCursor(self):
        if len(self.cursorStack):
            self.currentCursor = cursorName = self.cursorStack.pop(-1)
            self.configure(cursor=cursorsDict[cursorName])

    
    def findButton(self, action, actionDict):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        for b in (1,2,3):
            d = self.mouseButtonActions[actionDict][b]
            for mod in self.mouseButtonModifiers:
                if d[mod]==action: return b, mod
        return None, None

    
    def _setFov(self, fov):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        if fov > 0.0 and fov < 180.0: self.fovy = fov
        else: raise AttributeError('fov has to be < 180.0 and > 0.0 was %f'%fov)


    def _setLookFrom(self, val):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        mat = numpy.array(val, 'f')
        self.lookFrom = mat
        self.direction = self.lookAt - self.lookFrom
        

    def Set(self, check=1, redo=1, **kw):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Set various camera parameters
"""
        #print "Camera.Set", redo
        if __debug__:
            if check:
                apply( checkKeywords, (self.name,self.setKeywords), kw)
        
        val = kw.get( 'tagModified', True )
        assert val in [True, False]
        self._modified = val

        # we disable redraw because autoRedraw could cause ReallyRedraw
        # to set camera._width and camera._height BEFORE te window gets
        # resized by Tk
        if self.viewer.autoRedraw:
            restoreAutoRedraw = True
            self.viewer.stopAutoRedraw()
        else:
            restoreAutoRedraw = False

        if restoreAutoRedraw:
            self.viewer.startAutoRedraw()

        fov = kw.get( 'fov')
        if not fov is None:
            if fov > 0.0 and fov < 180.0: self.fovy = fov
            else: raise AttributeError('fov has to be < 180.0 and > 0.0 was %f'%fov)

        near = kw.get( 'near', self.near)
        far = kw.get( 'far', self.far)
        self.zclipWidth = far - near

        near = kw.get( 'near')
        if not near is None:
            #print 'NEAR:', near, far,
            if near >= far:
                raise AttributeError('near sould be smaller than far')

            if near > far-self.zclipWidth:
                near = max(far-self.zclipWidth, self.nearDefault)
            # truncate near at 0.1
            self.near = max(near, self.nearDefault)
            self.near_real = near
            #print self.near, self.near_real, self.far, self.zclipWidth
            
        if not far is None:
            #print 'FAR:', self.near, far,
            if far <= self.near:
                self.far = self.near + self.zclipWidth
                #if self.near == self.nearDefault:
                #    self.far = self.near * 1.5
                #else:
                #    raise AttributeError('far sould be larger than near')
            else:
                self.far = far

            #print self.near, self.near_real, self.far, self.zclipWidth


        if fov or near or far:
            self.SetupProjectionMatrix()

        val = kw.get( 'color')
        if not val is None:
            color = colorTool.OneColor( val )
            if color:
                self.backgroundColor = color

        val = kw.get( 'antialiased')
        if not val is None:
            if val in jitter.jitterList:
                self.antiAliased = val
                if val!=0:
                    self.accumWeigth = 1.0/val
                    self.jitter = eval('jitter._jitter'+str(val))
            else: raise ValueError('antiAliased can only by one of', \
                                   jitter.jitterList)

        # NPR parameters

        #first derivative
        val = kw.get( 'd1ramp')
        if not val is None:
                self.d1ramp=val

        val = kw.get( 'd1scale')
        if not val is None:
            assert isinstance(val, float)
            assert val>=0.0
            self.d1scale = val

        val = kw.get( 'd1off')
        if not val is None:
            assert isinstance(val, int)
            self.d1off = val

        val = kw.get( 'd1cutL')
        if not val is None:
            assert isinstance(val, int)
            assert val>=0
            self.d1cutL = val

        val = kw.get( 'd1cutH')
        if not val is None:
            assert isinstance(val, int)
            assert val>=0
            self.d1cutH = val

        #second derivative
        val = kw.get( 'd2scale')
        if not val is None:
            assert isinstance(val, float)
            assert val>=0.0
            self.d2scale = val

        val = kw.get( 'd2off')
        if not val is None:
            assert isinstance(val, int)
            self.d2off = val

        val = kw.get( 'd2cutL')
        if not val is None:
            assert isinstance(val, int)
            assert val>=0
            self.d2cutL = val

        val = kw.get( 'd2cutH')
        if not val is None:
            assert isinstance(val, int)
            assert val>=0
            self.d2cutH = val

        val = kw.get( 'contours')
        if not val is None:
            assert val in [True, False]
            if self.contours != val:
                self.contours = val
                if val is True:
                    self.lastBackgroundColorInPhotoMode = self.backgroundColor
                    self.backgroundColor = (1.,1.,1.,1.)
                    if self.viewer.OverAllLightingIsOn == 1:
                        self.viewer.OverAllLightingIsOn = 0
                else:
                    self.backgroundColor = self.lastBackgroundColorInPhotoMode
                    if self.viewer.OverAllLightingIsOn == 0:
                        self.viewer.OverAllLightingIsOn = 1
                self.fog.Set(color=self.backgroundColor)
                self.viewer.deleteOpenglListAndCallRedrawAndCallDisableGlLighting()
                                                        
#            if val:
#                self.imCanvastop = Tkinter.Toplevel()
#                self.imCanvas = Tkinter.Canvas(self.imCanvastop)
#                self.imCanvas1 = Tkinter.Canvas(self.imCanvastop)
#                self.imCanvas2 = Tkinter.Canvas(self.imCanvastop)
#                self.imCanvas3 = Tkinter.Canvas(self.imCanvastop)
#                self.canvasImage = None
#                self.canvasImage1 = None
#                self.canvasImage2 = None
#                self.canvasImage3 = None
#            elif self.imCanvastop:
#                self.imCanvastop.destroy()
#                self.imCanvas = None
#                self.imCanvas1 = None
#                self.imCanvas2 = None
#                self.imCanvas3 = None

#                if hasattr(self, 'ivi'):
#                    self.ivi.Exit()
#                from opengltk.OpenGL import GL
#                if not hasattr(GL, 'GL_CONVOLUTION_2D'):
#                    print 'WARNING: camera.Set: GL_CONVOLUTION_2D nor supported'
#                    self.contours = False
#                else:
#                    from DejaVu2.imageViewer from PIL import ImageViewer
#                    self.ivi = ImageViewer(name='zbuffer')
                
        val = kw.get( 'boundingbox')
        if not val is None:
            if val in viewerConst.BB_MODES:
                self.drawBB = val
            else: raise ValueError('boundingbox can only by one of NO, \
ONLY, WITHOBJECT')

        val = kw.get( 'rotation')
        if not val is None:
            self.rotation = numpy.identity(4, 'f').ravel()
            mat = numpy.reshape(numpy.array(val), (4,4)).astype('f')
            self.ConcatRotation(mat)

        val = kw.get( 'translation')
        if not val is None:
            self.translation = numpy.zeros( (3,), 'f')
            mat = numpy.reshape(numpy.array(val), (3,)).astype('f')
            self.ConcatTranslation(mat, redo=redo)

        val = kw.get( 'scale')
        if not val is None:
            self.SetScale( val, redo=redo )

        val = kw.get( 'pivot')
        if not val is None:
            self.SetPivot( val )

        val = kw.get( 'direction')
        if not val is None:
            mat = numpy.reshape(numpy.array(val), (3,)).astype('f')
            self.direction = mat

        valLookFrom = kw.get('lookFrom')
        if valLookFrom is not None:
            mat = numpy.array(valLookFrom, 'f')
            assert mat.shape==(3,)
            self.lookFrom = mat
            self.direction = self.lookAt - self.lookFrom

        valLookAt = kw.get('lookAt')
        if valLookAt is not None:
            mat = numpy.array(valLookAt, 'f')
            assert mat.shape==(3,)
            self.lookAt = mat
            self.direction = self.lookAt - self.lookFrom

        valUp = kw.get('up')
        if valUp is not None:
            mat = numpy.array(valUp, 'f')
            assert mat.shape==(3,)
            self.up = mat

##         # compute .translation and .rotation used to build camera transformation
##         # FIXME .. could recompute only if one changed
##         if valLookAt is not None \
##           or valLookFrom is not None \
##           or valUp is not None:
##             glPushMatrix()
##             glLoadIdentity()
##             glMatrixMode(GL_MODELVIEW)
##             fx, fy, fz = self.lookFrom
##             ax, ay, az = self.lookAt
##             ux, uy, uz = self.up
##             gluLookAt( float(fx), float(fy), float(fx),
##                        float(ax), float(ay), float(az),
##                        float(ux), float(uy), float(uz) )
##             lMatrix = numpy.array(glGetDoublev(GL_MODELVIEW_MATRIX)).astype('f')
##             lRotation, lTranslation, lScale = self.Decompose4x4(lMatrix, cleanup=False)
##             glPopMatrix()
##             if numpy.any(lMatrix):
##                 self.rotation = lRotation
##                 self.translation = lTranslation
##             self.lookAt = numpy.reshape(numpy.array(valLookAt), (3,)).astype('f')

##         val = kw.get( 'lookAt')
##         if not val is None:
##             mat = numpy.reshape(numpy.array(val), (3,)).astype('f')
##             self.lookAt = mat
##             self.direction = self.lookAt - self.lookFrom

##         val = kw.get( 'lookFrom')
##         if not val is None:
##             mat = numpy.reshape(numpy.array(val), (3,)).astype('f')
##             self.lookFrom = mat
##             self.direction = self.lookAt - self.lookFrom

        val = kw.get( 'projectionType')
        if val==self.PERSPECTIVE:
            if self.projectionType==self.ORTHOGRAPHIC:
                self.projectionType = val
                self.OrthogonalToPerspective()
        elif val==self.ORTHOGRAPHIC:
            if self.projectionType==self.PERSPECTIVE:
                self.projectionType = val
                self.PerspectiveToOrthogonal()

        val = kw.get( 'sideBySideRotAngle')
        if not val is None:
            assert type(val)==types.FloatType
            self.sideBySideRotAngle = val
            if self.viewer:
                self.viewer.Redraw()

        val = kw.get( 'sideBySideTranslation')
        if not val is None:
            assert type(val)==types.FloatType
            self.sideBySideTranslation = val
            if self.viewer:
                self.viewer.Redraw()

        val = kw.get( 'TV3DTranslation')
        if not val is None:
            assert type(val)==types.FloatType
            self.TV3DTranslation = val
            if self.viewer:
                self.viewer.Redraw()

        val = kw.get( 'TV3DRotAngle')
        if not val is None:
            assert type(val)==types.FloatType
            self.TV3DRotAngle = val
            if self.viewer:
                self.viewer.Redraw()

        val = kw.get( 'TV3DScaling')
        if not val is None:
            assert len(val)==3
            assert type(val[0])==types.FloatType
            assert type(val[1])==types.FloatType
            assert type(val[2])==types.FloatType
            self.TV3DScaling = val
            if self.viewer:
                self.viewer.Redraw()
        
        val = kw.get( 'ssao')
        if not val is None:
            self.ssao = kw["ssao"]
        val = kw.get( 'SSAO_OPTIONS')   
        if not val is None:
            self.SSAO_OPTIONS = kw['SSAO_OPTIONS']

        val = kw.get( 'stereoMode')
        if val is not None:
            assert val in self.stereoModesList

            glDrawBuffer(GL_BACK)

            if self.viewer is None:
                warnings.warn("""Stereo buffers are not present
or not enabled on this system.

enableStereo must be set to True in:
~/.mgltools/(ver_number)/DejaVu2/_dejavurc
"""
)
                val = 'MONO'
            elif self.viewer.activeStereoSupport is False:
                if val == 'STEREO_BUFFERS':
                    warnings.warn("""Stereo buffers are not present
or not enabled on this system.

enableStereo must be set to True in:
~/.mgltools/(ver_number)/DejaVu2/_dejavurc
"""
)
                    val = 'MONO'

            self.stereoMode = val
            if self.viewer:
                self.viewer.Redraw()

        val = kw.get( 'suspendRedraw')
        if val is not None:
            assert val in [True, False]
            self.suspendRedraw = val

        val = kw.get( 'drawThumbnail')
        if val is not None:
            assert val in [True, False]
            self.drawThumbnailFlag = val
        

    def OrthogonalToPerspective(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Compute left, right, top, bottom from field of view"""
        d = self.near + (self.far - self.near)*0.5
        self.fovy = (math.atan(self.top/d) * 360.0) / math.pi
        self.SetupProjectionMatrix()
        

    def PerspectiveToOrthogonal(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Compute left, right, top, bottom from field of view"""

        aspect = self._width / float(self._height)

        fov2 = (self.fovy*math.pi) / 360.0  # fov/2 in radian
        d = self.near + (self.far - self.near)*0.5

        self.top = d*math.tan(fov2)
        self.bottom = -self.top
        self.right = aspect*self.top
        self.left = -self.right
        self.SetupProjectionMatrix()


    def SetupProjectionMatrix(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Setup the projection matrix"""

        if not self.initialized:
            return

        self.Activate()
        #print 'SETUP PROJ', self._width, self._height
        glViewport(0, 0, self._width, self._height)

        glMatrixMode(GL_TEXTURE)
        glLoadIdentity()

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity()

##         if self.viewer.tileRender:
##             print 'FUGU'
##             if self.projectionType==self.PERSPECTIVE:
##                 self.viewer.tileRenderCtx.perspective(self.fovy,
##                            float(self._width)/float(self._height),
##                            self.near, self.far)
##             else:
##                 self.viewer.tileRenderCtx.ortho(self.left, self.right,
##                                                 self.bottom, self.top, 
##                                                 self.near, self.far)
##             print 'near', self.viewer.tileRenderCtx.Near
            
        if self.projectionType==self.PERSPECTIVE:
            
            # protect from bug in mac intel when zomming on molecule 1BX4
            if sys.platform == 'darwin':
                if self.fovy < .003:
                    self.fovy = .003

            gluPerspective(float(self.fovy),
                           float(self._width)/float(self._height),
                           float(self.near), float(self.far))
        else:
            glOrtho(float(self.left), float(self.right),
                    float(self.bottom), float(self.top), 
                    float(self.near), float(self.far))
            
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()


##          from opengltk.OpenGL import GLU

##          GLU.gluLookAt(0., 0., 30., 0.,0.,0., 0.,1.,0.)
##          print 'Mod Mat:', glGetDoublev(GL_MODELVIEW_MATRIX)
##          glLoadIdentity();
    #  The first 6 arguments are identical to the glFrustum() call.
    # 
    #  pixdx and pixdy are anti-alias jitter in pixels. 
    #  Set both equal to 0.0 for no anti-alias jitter.
    #  eyedx and eyedy are depth-of field jitter in pixels. 
    #  Set both equal to 0.0 for no depth of field effects.
    #
    #  focus is distance from eye to plane in focus. 
    #  focus must be greater than, but not equal to 0.0.
    #
    #  Note that AccFrustum() calls glTranslatef().  You will 
    #  probably want to insure that your ModelView matrix has been 
    #  initialized to identity before calling accFrustum().

    def AccFrustum(self, left, right, bottom, top, near, far,
                   pixdx, pixdy, eyedx, eyedy, focus):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        viewport = numpy.array(glGetDoublev (GL_VIEWPORT))

        xwsize = right - left
        ywsize = top - bottom

        dx = -(pixdx*xwsize/viewport[2] + eyedx*near/focus)
        dy = -(pixdy*ywsize/viewport[3] + eyedy*near/focus)

        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glFrustum (float(left + dx), float(right + dx),
                   float(bottom + dy), float(top + dy),
                   float(near), float(far))
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glTranslatef (float(-eyedx), float(-eyedy), 0.0)


    #  The first 4 arguments are identical to the gluPerspective() call.
    #  pixdx and pixdy are anti-alias jitter in pixels. 
    #  Set both equal to 0.0 for no anti-alias jitter.
    #  eyedx and eyedy are depth-of field jitter in pixels. 
    #  Set both equal to 0.0 for no depth of field effects.
    #
    #  focus is distance from eye to plane in focus. 
    #  focus must be greater than, but not equal to 0.0.
    #  Note that AccPerspective() calls AccFrustum().

    def AccPerspective(self, pixdx, pixdy, eyedx, eyedy, focus,
                       left=None, right=None, bottom=None, top=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Build perspective matrix for jitter"""

        from math import pi, cos, sin

        fov2 = self.fovy*pi / 360.0;
        if top is None:
            top = self.near / (cos(fov2) / sin(fov2))
        if bottom is None:
            bottom = -top
        if right is None:
            right = top * float(self._width)/float(self._height)
        if left is None:
            left = -right;
        
        self.AccFrustum (left, right, bottom, top, self.near, self.far,
                         pixdx, pixdy, eyedx, eyedy, focus)


    def InitGL(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Initialize some GL features"""

        self.Activate()
        glEnable(GL_CULL_FACE)
        glEnable(GL_NORMALIZE)  # required if glScale is used,
                                    # else the lighting doesn't work anymore
        glDepthFunc(GL_LESS)
        glEnable(GL_DEPTH_TEST)

        # blend function used for line anti-aliasing
        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA)
        self.initialized = 1
        

## moved the width and height into viewer.ReallyRedraw
##     def ReallyExpose(self, *dummy):
##         """Redraw the widget.
##         Make it active, update tk events, call redraw procedure and
##         swap the buffers.  Note: swapbuffers is clever enough to
##         only swap double buffered visuals."""
##         self.Activate()
##         self._width = self.winfo_width()
##         self._height = self.winfo_height()
##         #self.SetupProjectionMatrix()
##         #self.InitGL()
##         self.Redraw()

    

    def drawOneObjectThumbnail(self, obj):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.ActivateClipPlanes( obj.clipP, obj.clipSide )
        if obj.scissor:
            glEnable(GL_SCISSOR_TEST)
            glScissor(obj.scissorX, obj.scissorY, obj.scissorW, obj.scissorH)

        inst = 0
        v = obj.vertexSet.vertices.array
        v = numpy.reshape(v, (-1,3)).astype('f')

        for m in obj.instanceMatricesFortran:
            inst = inst + 1
            glPushMatrix()
            glMultMatrixf(m)
            if isinstance(obj, Transformable):            
                #v = obj.vertexSet.vertices.array
                v, indices = obj.getVisibleVertices(picking=True)
                if len(v.shape) >2 or v.dtype.char!='f':
                    v = numpy.reshape(v, (-1,3)).astype('f')
                   
                if len(v) >0:
                    #t1 = time()
                    namedPointsWithNames(v, indices)
                    #glVertexPointer(2, GL_FLOAT, 0, v)
                    #glEnableClientState(GL_VERTEX_ARRAY)
                    #glDrawArrays(GL_POINTS, 0, len(v) )
                    #glDisableClientState(GL_VERTEX_ARRAY)

            else: #Insert2d
                obj.pickDraw()

            glPopMatrix()

        for c in obj.clipP: # disable object's clip planes
            c._Disable()
        if obj.scissor:
            glDisable(GL_SCISSOR_TEST)


    def DrawObjThumbnail(self, obj):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """draws an object for picking purposes. If type is 'vertices' a
display list for vertices identification is used. If type is parts
a display list identifying geometric primitives such as triangles or
lines is used
"""
        if isinstance(obj, Transformable):
            glPushMatrix()
            obj.MakeMat()
            self.ActivateClipPlanes( obj.clipPI, obj.clipSide )
            inst = 0
            for m in obj.instanceMatricesFortran:
                glPushMatrix()
                glMultMatrixf(m)
                for child in obj.children:
                    if child.visible:
                        self.DrawObjThumbnail(child)
                glPopMatrix()
                inst = inst + 1
            if obj.visible:
                self.drawOneObjectThumbnail(obj)
            for c in obj.clipPI: # disable object's clip planes that are
                c._Disable()     # inherited by children
            glPopMatrix()     # Restore the matrix


    def RedrawThumbnail(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        glClear(GL_DEPTH_BUFFER_BIT)
        glPushMatrix()
        glPointSize(1.0)
        self.BuildTransformation()  # camera transformation
        obj = self.viewer.rootObject
        obj.MakeMat() 
        if len(obj.clipPI):
            self.ActivateClipPlanes( obj.clipPI, obj.clipSide )
        inst = 0
        for m in obj.instanceMatricesFortran:
            glPushMatrix()
            glMultMatrixf(m)
            for child in obj.children:
                if child.visible:
                    self.DrawObjThumbnail(child)
            glPopMatrix()
            inst = inst + 1
        for c in obj.clipPI: c._Disable()
        glPopMatrix()


    def drawThumbnail(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.Activate()

        glViewport(0, 0, self._width/10, self._height/10)

        glMatrixMode (GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity ()

        # create projection matrix and modeling
#        self.SetupProjectionMatrix()
        if self.projectionType==self.PERSPECTIVE:
            gluPerspective(float(self.fovy),
                           float(self._width)/float(self._height),
                           float(self.near),float(self.far))
        else:
            glOrtho(float(self.left),float(self.right),
                    float(self.bottom), float(self.top), 
                    float(self.near), float(self.far))
            
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        self.RedrawThumbnail()
        glFlush ()

        glMatrixMode (GL_PROJECTION)
        glPopMatrix ()
        glMatrixMode(GL_MODELVIEW)

        
    def DoPick(self, x, y, x1=None, y1=None, type=None, event=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Makes a redraw in GL_SELECT mode to pick objects
"""
        if self.stereoMode != 'MONO':
            return
        if type is None:
            type = self.viewer.pickLevel
        if x1 is None: x1 = x
        if y1 is None: y1 = y
        #print 'pick at', x, y, x1, y1

        self.Activate()
        vp = [0,0,self._width, self._height]

        selectBuf = glSelectBuffer( 1000000 ) 
        # 500000 was insuficient on bigger molecule as 2plv.pdb with MSMS

        y1 = vp[3] - y1
        y = vp[3] - y

        dx = x - x1
        dy = y - y1
        #x = x
        pickWinSize = 10
        
        if math.fabs(dx) < pickWinSize: dx=pickWinSize
        else: x = x1 + (dx/2)

        if math.fabs(dy) < pickWinSize: dy=pickWinSize
        else: y = y1 + (dy/2)

        if dx==pickWinSize and dy==pickWinSize:
            mode = 'pick'
        else: 
            mode = 'drag select'

        abdx = int(math.fabs(dx))
        abdy = int(math.fabs(dy))
        #print 'pick region ', x-math.fabs(dx), y-math.fabs(dy), x+math.fabs(dx), y+math.fabs(dy), mode

        # debug
        glRenderMode (GL_SELECT)
        self.renderMode = GL_SELECT

        glInitNames()
        glPushName(0)

        glMatrixMode (GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity ()

        # create abs(dx)*abs(dy) pixel picking region near cursor location
        gluPickMatrix( x, y, abdx, abdy, vp)

        # we save it so it can be reapplied if the projection matrix is recreated
        self.pickMatrix = glGetFloatv(GL_PROJECTION_MATRIX)

        # create projection matrix and modeling
#        self.SetupProjectionMatrix()
        if self.projectionType==self.PERSPECTIVE:
            gluPerspective(self.fovy,float(self._width)/float(self._height),
                           self.near, self.far)
        else:
            glOrtho(float(self.left),float(self.right),
                    float(self.bottom),float(self.top), 
                    float(self.near), float(self.far))
            
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()

        self.RedrawPick(type);
        glFlush ()

        # debug
        ###self.Activate()
        ###return PickObject('pick')

        self.renderMode = GL_RENDER;
        # get the number of hits

        selectHits = glRenderMode (GL_RENDER)
        #print "hits", selectHits

        #removed because it generates bug in displaylist
        #if selectHits == 0:
        #    # if we pick the background, root becomes selected
        #    self.viewer.SetCurrentObject(self.viewer.rootObject)
        
        # restore projection matrix
        glMatrixMode (GL_PROJECTION)
        glPopMatrix ()

        # unproject points x and y
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()  # restore matrix for unproject
        self.BuildTransformation()
        self.viewer.rootObject.MakeMat()

        self.unProj_model = glGetDoublev(GL_MODELVIEW_MATRIX)
        self.unProj_proj = glGetDoublev(GL_PROJECTION_MATRIX)
        self.unProj_view = glGetIntegerv(GL_VIEWPORT)
        
        p1 = gluUnProject( (x, y, 0.), self.unProj_model, self.unProj_proj,
                           self.unProj_view)
        p2 = gluUnProject( (x, y, 1.), self.unProj_model, self.unProj_proj,
                           self.unProj_view)
        glLoadIdentity()

        if selectHits:
            if mode == 'pick':
                pick = self.handlePick(selectHits, selectBuf, type)
            else:
                if dx > 0: op = 'add'
                else: op = 'remove'
                pick = self.handleDragSelect(selectHits, selectBuf, type, op)
        else:
            pick = PickObject('pick', self)

##          if selectHits and self.viewer.showPickedVertex:
##              self.DrawPickingSphere(pick)
            
        pick.p1 = p1
        pick.p2 = p2
        pick.box = (x, y, x1, y1)
        pick.event = event
        
        self.viewer.lastPick = pick

        #DEBUG used to debug picking
#        from IndexedPolylines import IndexedPolylines
#        l = IndexedPolylines('line', vertices = (self._p1,self._p2), faces=((0,1),) )
#        self.viewer.AddObject(l)
#        return o, parts, self.viewer.lastPickedVertex, self._p1, self._p2 

        return pick

#    ## DEPRECATED, is not using instance for computing coordinates (MS 12/04)
#    def DrawPickingSphere(self, pick):
#        """display a transient sphere at picked vertex"""
#        from warnings import warn
#        warnings.warn('DrawPickingSphere is deprecated',
#                      DeprecationWarning, stacklevel=2)
#        
#        obj = pick.hits.keys()[0]
#        varray = obj.vertexSet.vertices.array
#        coords = numpy.take( varray, pick.hits[obj], axis=0 )
#        p = self.viewer.pickVerticesSpheres
#        mat = obj.GetMatrix( obj.LastParentBeforeRoot() )
#        mat = numpy.transpose(mat)
#        p.Matrix = numpy.reshape(mat, (16, ))
#        p.SetMatrix(mat)
#        self.viewer.pickVerticesSpheres.Set(vertices=coords,
#                                     transient=self.viewer.pickReminiscence,
#                                     redo=redo)
#        self.Redraw()


    def xyzFromPixel( self, winx, winy):
        """compute x, and y values in scene from x and y corrdinates in event
        find z from Zbuffer.

        pt3d, background = xyzFromPixel( self, winx, winy)

        background is True if pick occured on background"""

        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        # unproject pixel
        self.Activate()
        glPushMatrix()
        self.BuildTransformation()  # camera transformation
        nar = numpy.zeros(1, 'f')
        glFinish()
        from opengltk.extent import _gllib as gllib
        gllib.glReadPixels( winx, winy, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT,
                            nar)
        winz = float(nar[0])
        background = winz == 1.0

        #print 'xyzFromPixel: picks',winz, background

        #print "winx, winy, winz", winx, winy, winz
        l3dPoint = gluUnProject( (winx, winy, winz), 
                                 glGetDoublev(GL_MODELVIEW_MATRIX),
                                 glGetDoublev(GL_PROJECTION_MATRIX),
                                 glGetIntegerv(GL_VIEWPORT)
                               )
        #print "l3dPoint", l3dPoint
        glPopMatrix()
        return l3dPoint, background


    def handleDragSelect(self, selectHits, selectBuf, type, operation):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        pick = PickObject('drag select', self, type)
        p = 0
        for i in range(selectHits):
            nb_names = selectBuf[p]
            z1 = float(selectBuf[p+1]) #/ 0x7fffffff
            end = int(p+3+nb_names)
            instance = list(selectBuf[p+3:end-2])
            obj = self.objPick[int(selectBuf[end-2])]
            vertex = selectBuf[end-1]
            #print 'vertex', vertex, obj.name, instance
            pick.add(object=obj, vertex=vertex, instance=instance,
                     operation=operation)
            p = end
        return pick


    def handlePick(self, selectHits, selectBuf, type):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """We are looking for the vertex closest to the viewer"""
        # now handle selection
        mini = 9999999999.9
        p = 0
        vertex = None
        obj = None
        # loop over hits. selectBuf contains for each hit:
        #   - number of names for that hit (num)
        #   - z1
        #   - z2
        #   - instNumRoot, instNumRootParent_1, instNumRootParent_2, ... ,
        #      geomIndex, vertexNum
        # For each parent we have an instance number
        # the geomIndex is the index of the geometry in self.objPick
        # vertexNum is the index of the vertex in the geometry's vertexSet
        #
        pick = PickObject('pick', self, type)
        
        for i in range(selectHits):
            nb_names = selectBuf[p]
            z1 = float(selectBuf[p+1]) #/ 0x7fffffff
            # compute the end of the pick info for this hit
            end = int(p+3+nb_names)
            #print 'vertex', selectBuf[end-1], self.objPick[selectBuf[end-2]], list(selectBuf[p+3:end-2]), z1
            if z1 < mini:
                #end = p+3+nb_names
                mini = z1
                instance = list(selectBuf[p+3:end-2])
                obj = self.objPick[int(selectBuf[end-2])]
                vertex = selectBuf[end-1]
            # advance p to begin of next pick hit data
            p = end
        pick.add(object=obj, vertex=vertex, instance=instance )
        return pick
    

    def drawOneObjectPick(self, obj, type):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        # skip object if it is in list of objects not shown in this Camera
        if obj.hiddenInCamera.has_key(self):
            return

        lIsInstanceTransformable = isinstance(obj, Transformable)
        if lIsInstanceTransformable:
            self.ActivateClipPlanes( obj.clipP, obj.clipSide )
            if obj.scissor:
                glEnable(GL_SCISSOR_TEST)
                glScissor(obj.scissorX, obj.scissorY,
                         obj.scissorW, obj.scissorH)

        obj.pickNum = self.pickNum
        self.objPick.append(obj)

        inst = 0
        
        #if lIsInstanceTransformable and type=='vertices':
        #    v, indices = obj.getVisibleVertices()
        #    v = numpy.reshape(v, (-1,3)).astype('f')

        for m in obj.instanceMatricesFortran:
            glPushName(inst) # instance number
            inst = inst + 1
            glPushName(self.pickNum) # index into self.objPick
            glPushMatrix()
            glMultMatrixf(m)
            if lIsInstanceTransformable:            
                if type=='vertices':
                    v, indices = obj.getVisibleVertices(picking=True)
                    if len(v.shape) >2 or v.dtype.char!='f':
                        v = numpy.reshape(v, (-1,3)).astype('f')
                   
                    if len(v) >0:
                        #t1 = time()
                        namedPointsWithNames(len(v), v.astype('float32'),
                                             indices.astype('int32'))
                        #print 'draw points %d poitns for %s:'%(len(v), obj.name)#, time()-t1
    ##                  i = 0
    ##                  for p in v:
    ##                      glPushName(i)
    ##                      glBegin(GL_POINTS)
    ##                      glVertex3f(p[0], p[1], p[2])
    ##                      glEnd()
    ##                      glPopName()
    ##                      i = i + 1
                    
    # can't be used since we cannot push and pop names
    ##                  glVertexPointer(2, GL_FLOAT, 0, v)
    ##                  glEnableClientState(GL_VERTEX_ARRAY)
    ##                  glDrawArrays(GL_POINTS, 0, len(v) )
    ##                  glDisableClientState(GL_VERTEX_ARRAY)
                     
                elif type=='parts':
                    if obj.pickDpyList:
                        #print "pick displayFunction"
                        currentcontext = self.viewer.currentCamera.getContext()
                        if currentcontext != obj.pickDpyList[1]:
                            warnings.warn("""DisplayFunction failed because the current context is the wrong one""")
                            #print "currentcontext != obj.pickDpyList[1]", currentcontext, obj.pickDpyList[1]
                        else:
                            print '#%d'%obj.pickDpyList[0], currentcontext, "glCallList Camera"
                            glCallList(obj.pickDpyList[0])
                    else:
                        #print "displayFunction"
                        obj.DisplayFunction()
                else:
                    print 'Error: bad type for PickRedraw: ',type
            else: #Insert2d
                obj.pickDraw()

            glPopMatrix()
            glPopName() # instance number

            glPopName() # index into self.objPick
        self.pickNum = self.pickNum + 1

        if lIsInstanceTransformable:
            for c in obj.clipP: # disable object's clip planes
                c._Disable()
            if obj.scissor:
                glDisable(GL_SCISSOR_TEST)


    def DrawObjPick(self, obj, type):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """draws an object for picking purposes. If type is 'vertices' a
        display list for vertices identification is used. If type is parts
        a display list identifying geometric primitives such as triangles or
        lines is used"""

        lIsInstanceTransformable = isinstance(obj, Transformable)
        if lIsInstanceTransformable:
            glPushMatrix()
            obj.MakeMat()
            self.ActivateClipPlanes( obj.clipPI, obj.clipSide )

        inst = 0
        for m in obj.instanceMatricesFortran:
            glPushName(inst) # instance number
            glPushMatrix()
            glMultMatrixf(m)
            for child in obj.children:
                if child.visible:
                    self.DrawObjPick(child, type)
            glPopMatrix()
            glPopName() # instance number
            inst = inst + 1
        
        if obj.pickable:
            self.drawOneObjectPick(obj, type)

        if lIsInstanceTransformable:
            for c in obj.clipPI: # disable object's clip planes that are
                c._Disable()     # inherited by children

            glPopMatrix()     # Restore the matrix

               
    def RedrawPick(self, type):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.pickNum = 0
        self.objPick = [ ]
        glClear(GL_DEPTH_BUFFER_BIT)
        glPushMatrix()
        glPointSize(1.0)
        self.BuildTransformation()  # camera transformation
        obj = self.viewer.rootObject
        obj.MakeMat() 
        if len(obj.clipPI):
            self.ActivateClipPlanes( obj.clipPI, obj.clipSide )
        inst = 0
        for m in obj.instanceMatricesFortran:
            glLoadName(inst)
            glPushMatrix()
            glMultMatrixf(m)
            for child in obj.children:
                if child.visible:
                    self.DrawObjPick(child, type)
            glPopMatrix()
            inst = inst + 1
        for c in obj.clipPI: c._Disable()
        glPopMatrix()

        
    def drawRect(self, P1, P2, P3, P4, fill=0):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        if fill: prim=GL_POLYGON
        else: prim=GL_LINE_STRIP
        glBegin(prim)
        glVertex3f( float(P1[0]), float(P1[1]), float(P1[2]) )
        glVertex3f( float(P2[0]), float(P2[1]), float(P2[2]) )
        glVertex3f( float(P3[0]), float(P3[1]), float(P3[2]) )
        glVertex3f( float(P4[0]), float(P4[1]), float(P4[2] ))
        glVertex3f( float(P1[0]), float(P1[1]), float(P1[2] ))
        glEnd()


    def Insert2dPickingCallBack(self, event):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        # this shoudn't be necessary but somewhere/sometimes the 
        # current object is set without a call to SetCurrentObject
        # So self.viewer.currentCamera.bindAllActions('Insert2d') is no set
        if isinstance(self.viewer.currentObject, Insert2d) is False:
            assert isinstance(self.viewer.currentObject, Transformable), self.viewer.currentObject
            self.viewer.SetCurrentObject(self.viewer.currentObject)
            return

        self.viewer.currentObject.respondToMouseMove(event)


    def SetInsert2dPicking(self, event):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        num = str(self.findButton('picking', 'Insert2d')[0])            
        ehm = self.eventManager
        for mod in self.mouseButtonModifiers:
            if mod == 'None': 
                mod = ''
            else: 
                mod = mod + '-'
            ev = '<' + mod + 'B' + num + '-Motion>'
            #print "ev", ev
            ehm.AddCallback(ev, self.Insert2dPickingCallBack )       
        pick = self.DoPick(event.x, event.y, event=event)
        #if pick and len(pick.hits):
        self.viewer.processPicking(pick)


    def initSelectionRectangle(self, event):
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
        x2 = self._x1 = event.x
        y2 = self._y1 = self._height - event.y

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

        # set "<ButtonMotion-1>" to call draw selection rectangle 
        # add un-drawing of last selection rectangle and call functions
        ehm = self.eventManager
        num = str(event.num)#str(self.findButton('picking', 'Object')[0])
        for mod in self.mouseButtonModifiers:
            if mod=='None': mod=''
            else: mod=mod+'-'
            ev = '<'+mod+'B'+num+'-Motion>'
            ehm.AddCallback(ev, self.drawSelectionRectangle )
            ehm.AddCallback("<"+mod+"ButtonRelease-"+num+">",
                            self.endSelectionRectangle )


    def drawSelectionRectangle(self, event):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
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
        x2 = event.x
        y2 = self._height - event.y

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
        self.drawRect( self._P1, self._P2, self._P3, self._P4, self.fill )
        glPopMatrix()
        glFlush()
        if self.fillSelectionBox:
            if self.afid:
                self.after_cancel(self.afid)
            self.afid = self.after(self.fillDelay, self.DrawFilledSelectionBox)
        

    def DrawFilledSelectionBox(self, event=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.fill = 1
        self.Activate()
        glPushMatrix()
        self.BuildTransformation()  # camera transformation

        # draw over previous rectangle
        self.drawRect( self._P1, self._P2, self._P3, self._P4, 0 )

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
        glDisable(GL_CULL_FACE)
        glColor3f( float(self.selectionColor[0]), float(self.selectionColor[1]),
                      float(self.selectionColor[2] ))
        self.drawRect( self._P1, self._P2, self._P3, self._P4, 1 )
        glFlush()
        glPopMatrix()


    def endSelectionRectangle(self, event):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        #print "Camera.endSelectionRectangle"

        if not self.selectDragRect: return

        if self.afid:
            self.after_cancel(self.afid)

        # remove "<Any-ButtonMotion-1>" to call draw selection rectangle 
        # remove un-drawing of last selection rectangle and call functions
        ehm = self.eventManager
        num = str(event.num)#str(self.findButton('picking', 'Object')[0])
        for mod in self.mouseButtonModifiers:
            if mod=='None': mod=''
            else: mod=mod+'-'
            ev = '<'+mod+'B'+num+'-Motion>'
            ehm.RemoveCallback(ev, self.drawSelectionRectangle )
            ehm.RemoveCallback("<"+mod+"ButtonRelease-"+num+">",
                               self.endSelectionRectangle )

        self.selectDragRect = 0
        self.Activate()
        glPushMatrix()
        self.BuildTransformation()  # camera transformation

        # this line is required for last rectangle to disappear !
        # not sure why
        # oct 2001: apparently not necessary
        #glEnable(GL_COLOR_LOGIC_OP)

        #self.history.append( ('hide', self._P3) )
        self.drawRect(self._P1, self._P2, self._P3, self._P4)

        #self.drawRect(self._P1, self._P2, self._P3, self._P4)
        glPopMatrix()
        
        glEnable(GL_DEPTH_TEST)
        glDisable(GL_COLOR_LOGIC_OP)
        #glEnable(GL_LIGHTING)
        if self.viewer is not None:
            self.viewer.enableOpenglLighting()
        glDrawBuffer(GL_BACK)
        pick = self.DoPick(event.x, event.y, self._x1, self._height-self._y1,
                           event=event)
        del self._x1
        del self._y1
        del self._P1
        del self._P2
        del self._P3
        del self._P4
        if self.viewer and len(pick.hits):
            self.viewer.processPicking(pick)

        #print "len(pick.hits)", len(pick.hits)


    def DoubleSelectPick_cb(self, event):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Handle a double pick event in this camera"""

        self.DoPick(event.x, event.y, event=event)
        self.viewer.BindTrackballToObject(self.viewer.rootObject)


    ## DEPRECATED, is not using instance for computing coordinates
    def CenterPick_cb(self, event):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Handle a pick event. Center the cur. object on the picked vertex"""

        from warnings import warn
        warnings.warn('CenterPick_cbg is deprecated',
                      DeprecationWarning, stacklevel=2)

        pick = self.DoPick(event.x, event.y, event=event)
        if len(pick.hits)==0: return
        
        g = numpy.zeros( (3), 'f' )
        object = self.viewer.currentObject
        inv = object.GetMatrixInverse()
        for obj, vertInd in pick.hits.items():
            if len(vertInd)==0:
                print 'WARNING: object',obj.name,' is not vertex pickable'
            else:
                m = numpy.dot(inv, obj.GetMatrix())
                vert = obj.vertexSet.vertices * m
                vertsel = numpy.take(vert, vertInd, axis=0)
                g = g + numpy.sum(vertsel)/len(vertsel)

        object.SetPivot(g)

##                varray = obj.vertexSet.vertices.array
##                vert = numpy.take(varray, vertInd, axis=0)
##                vert = numpy.concatenate( vert, numpy.ones( (1,len(vert))), 1)
##                vert = numpy.dot(m, vert)
##                obj.SetPivot( (numpy.sum(vert)/len(vert))[:3] )
##                 obj.SetPivot(self.viewer.lastPickedVertex[1])


    def TransfCamera_cb(self, event):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.viewer.BindTrackballToCamera(self.viewer.currentCamera)


    def ActivateClipPlanes(self, cplist, cpside):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """activate a list of clipping planes"""

        if len(cplist)==0: return
        glPushMatrix()
        glLoadIdentity()
        glMultMatrixf(self.rootObjectTransformation)

        glPushAttrib(GL_CURRENT_BIT | GL_LIGHTING_BIT |
                        GL_LINE_BIT)

        glDisable (GL_LIGHTING)
        for c in cplist:         # display them
            if c.visible: 
                # if we push a name here we won't know how many values
                # to skip in handlePick after z1 and z2
                #if self.renderMode == GL_SELECT:
                #    c.pickNum = self.pickNum
                #    self.objPick.append(c)
                #    self.pickNum = self.pickNum + 1
                #    glLoadName(c.pickNum)
                c.DisplayFunction()

        glPopAttrib()

        for c in cplist:         # enable them 
            c._Enable(cpside[c.num])

        glPopMatrix();


    def DrawAllDirectionalLights(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Draw all directional lights"""

        glDisable (GL_LIGHTING)
        for l in self.viewer.lights:
            if l.positional is True or l.visible in (False,0):
                continue
            if self.renderMode == GL_SELECT:
                l.pickNum = self.pickNum
                self.objPick.append(l)
                self.pickNum = self.pickNum + 1
                glLoadName(l.pickNum)
            l.DrawDirectionalLight(self)




    def drawWithCap(self, obj):
        MY_CLIP_PLANE = obj.cap
        # compute clip plane number
        cpn = MY_CLIP_PLANE - GL.GL_CLIP_PLANE0
        #GL.glPushAttrib(GL.GL_ENABLE_BIT | GL.GL_STENCIL_BUFFER_BIT |
        #                GL.GL_POLYGON_BIT | GL_CURRENT_BIT)
        
        ##
        ## SINGLE PASS
        ##
##         GL.glDepthMask(0)
##         GL.glColorMask(0,0,0,0)
##         GL.glDisable(GL_CULL_FACE)
##         GL.glEnable(GL.GL_STENCIL_TEST)
##         GL.glEnable(GL.GL_STENCIL_TEST_TWO_SIDE_EXT)

##         GL.glActiveStencilFaceEXT(GL.GL_BACK)
##         GL.glStencilOp(GL.GL_KEEP,            # stencil test fail
##                     GL.GL_KEEP,            # depth test fail
##                     GL.GL_DECR_WRAP_EXT)  # depth test pass
##         GL.glStencilMask(~0)
##         GL.glStencilFunc(GL_ALWAYS, 0, ~0)

##         GL.glActiveStencilFaceEXT(GL.GL_FRONT)
##         GL.glStencilOp(GL.GL.GL_KEEP,            # stencil test fail
##                     GL.GL.GL_KEEP,            # depth test fail
##                     GL.GL.GL_INCR_WRAP_EXT);  # depth test pass
##         GL.glStencilMask(~0)
##         GL.glStencilFunc(GL.GL_ALWAYS, 0, ~0)
##         obj.Draw()

        ## 2 passes version
        GL.glEnable(GL.GL_CULL_FACE)
        GL.glEnable(GL.GL_STENCIL_TEST)
        glStencilMask(1)
        GL.glEnable(GL.GL_DEPTH_TEST)
        #GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glClear(GL.GL_STENCIL_BUFFER_BIT)
        GL.glColorMask(GL.GL_FALSE, GL.GL_FALSE, GL.GL_FALSE, GL.GL_FALSE)

        # first pass: increment stencil buffer value on back faces
        GL.glStencilFunc(GL.GL_ALWAYS, 0, 0)
        GL.glStencilOp(GL.GL_KEEP, GL.GL_KEEP, GL.GL_INCR)
        GL.glPolygonMode(GL.GL_BACK, GL.GL_FILL)
        GL.glCullFace(GL.GL_FRONT) # render back faces only
        obj.Draw()

        # second pass: decrement stencil buffer value on front faces
        GL.glStencilOp(GL.GL_KEEP, GL.GL_KEEP, GL.GL_DECR)
        GL.glPolygonMode(GL.GL_FRONT, GL.GL_FILL)
        GL.glCullFace(GL.GL_BACK) # render front faces only
        obj.Draw()

        # FIXME this can create problems whent there are multiple geoms
        GL.glClear(GL_DEPTH_BUFFER_BIT)

        # drawing clip planes masked by stencil buffer content
        #_gllib.glBindTexture(GL_TEXTURE_2D, 0)      
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glColorMask(GL.GL_TRUE, GL.GL_TRUE, GL.GL_TRUE, GL.GL_TRUE)
        GL.glDisable(MY_CLIP_PLANE)
        GL.glDisable(GL.GL_CULL_FACE)
        #GL.glStencilFunc(GL.GL_ALWAYS, 0, 1)
        #GL.glStencilFunc(GL.GL_LESS, 0, 1)
        #GL.glStencilFunc(GL.GL_NOTEQUAL, 0, 1)
        GL.glStencilFunc(GL.GL_NOTEQUAL, 0, 1)
        GL.glMaterialfv( GL.GL_FRONT, GL.GL_DIFFUSE, obj.capFrontColor)
        GL.glMaterialfv( GL.GL_BACK, GL.GL_DIFFUSE, obj.capBackColor)
        GL.glBegin(GL.GL_QUADS) #rendering the plane quad. Note, it should be 
                                #big enough to cover all clip edge area.
        GL.glNormal3fv( (1, 0, 0))

        if obj.clipSide[cpn]==-1:
            pts = [(0,  9999, -9999), (0,  9999,  9999),
                   (0, -9999,  9999), (0, -9999, -9999) ]
        else:
            pts = [ (0, -9999, -9999), (0, -9999,  9999),
                    (0,  9999,  9999), (0,  9999, -9999) ]
            
        cp = self.viewer.clipP[cpn]
        xpts = self.transformPoints(cp.translation, cp.rotation, pts)
        #GL.glDrawBuffer(GL.GL_FRONT)
        for p in xpts:
            GL.glVertex3fv( p )
        GL.glEnd()

        #****** Rendering object ********* 

        #GL.glPopAttrib()
        GL.glDisable(GL.GL_STENCIL_TEST)
        GL.glClear(GL_STENCIL_BUFFER_BIT)

        GL.glEnable(MY_CLIP_PLANE) # enabling clip plane again

        obj.SetupGL()
        obj.Draw()
        GL.glStencilFunc(GL.GL_ALWAYS, 0, 0)
        glStencilMask(0)


    def transformPoints(self, trans, rot, points):
        tx,ty,tz = trans
        pos = []
        for xs,ys,zs in points:
            x = rot[0]*xs + rot[4]*ys + rot[8]*zs + tx
            y = rot[1]*xs + rot[5]*ys + rot[9]*zs + ty
            z = rot[2]*xs + rot[6]*ys + rot[10]*zs + tz
            pos.append( [x,y,z] )
        return pos


    def DrawOneObject(self, obj):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """set up clipping planes, activate scissors and call display function
for all instances of obj
"""
        #print "Camera.DrawOneObject", obj
        
        lIsInstancetransformable = isinstance(obj, Transformable)
        
        if lIsInstancetransformable:

                self.ActivateClipPlanes( obj.clipP, obj.clipSide )
        
                if obj.scissor:
                    glDisable(GL_LIGHTING)
                    lViewport = glGetIntegerv(GL_VIEWPORT)
                    GL.glMatrixMode(GL.GL_PROJECTION)
                    GL.glPushMatrix()
                    GL.glLoadIdentity()
                    GL.glOrtho(float(lViewport[0]),
                               float(lViewport[2]),
                               float(lViewport[1]),
                               float(lViewport[3]),
                               -1, 1)
                    GL.glMatrixMode(GL.GL_MODELVIEW)
                    GL.glPushMatrix()
                    GL.glLoadIdentity()
                    glEnable(GL_COLOR_LOGIC_OP)
                    glLogicOp(GL_XOR)            
                    glBegin(GL_LINE_STRIP)
                    glVertex2f( float(obj.scissorX), float(obj.scissorY ))
                    glVertex2f( float(obj.scissorX+obj.scissorW), float(obj.scissorY ))
                    glVertex2f( float(obj.scissorX+obj.scissorW), float(obj.scissorY+obj.scissorH))
                    glVertex2f( float(obj.scissorX), float(obj.scissorY+obj.scissorH ))
                    glVertex2f( float(obj.scissorX), float(obj.scissorY ))
                    glEnd()
                    glDisable(GL_COLOR_LOGIC_OP)
                    GL.glMatrixMode(GL.GL_PROJECTION)
                    GL.glPopMatrix()
                    GL.glMatrixMode(GL.GL_MODELVIEW)
                    GL.glPopMatrix()
                    glEnable(GL_SCISSOR_TEST)
                    glScissor(obj.scissorX, obj.scissorY, obj.scissorW, obj.scissorH)
                    #glEnable(GL_LIGHTING)
                    if self.viewer is not None:
                        self.viewer.enableOpenglLighting()
                    
                if obj.drawBB: obj.DrawBoundingBox()

        # hack to fix this inheritence problem without slowing down
        # the "many object case (see Expensive inheritence pb above)

        # MS oct 11 2001. removed this code since obj.shadModel is set to
        # the proper value when it is inherited (see RenderMode()) and
        # shading is called by SetupGL for each object
        #if not obj.inheritShading and \
        #   obj.shading != GL_NONE:
        #    glShadeModel(obj.shading)

        if (lIsInstancetransformable and obj.drawBB != viewerConst.ONLY) \
           or isinstance(obj, Insert2d):
            if lIsInstancetransformable and obj.texture:
                glActiveTexture(GL_TEXTURE0)
                _gllib.glBindTexture(GL_TEXTURE_2D, 0)
                obj.texture.Setup()

            if lIsInstancetransformable and obj.invertNormals is True \
              and isinstance(obj, Ellipsoids) is False \
              and isinstance(obj, Cylinders) is False \
              and (isinstance(obj, Spheres) is True and DejaVu2.enableVertexArray is True):
                lCwIsOn = True
                GL.glFrontFace(GL.GL_CW)
                #print "GL.glGetIntegerv(GL.GL_FRONT_FACE)",GL.glGetIntegerv(GL.GL_FRONT_FACE)
            else:
                lCwIsOn = False

            mi = 0
            for m in obj.instanceMatricesFortran:
                if lIsInstancetransformable and not obj.inheritMaterial:
                    obj.InitMaterial(mi)
                    obj.InitColor(mi)
                    mi = mi + 1
                glPushMatrix()
                glMultMatrixf(m)
                if obj.immediateRendering or \
                   (self.viewer.tileRender and obj.needsRedoDpyListOnResize):
                    if hasattr(obj, "cap") and obj.cap:
                        self.drawWithCap(obj)
                    else:
                        obj.Draw()
                    #obj.Draw()
                elif self.viewer.singleDpyList:
                    obj.Draw()
                else:
                    obj.DisplayFunction()
                glPopMatrix()

            if lCwIsOn is True:
                GL.glFrontFace(GL.GL_CCW)
                
            #print obj.name, self.glError()
            if lIsInstancetransformable and obj.texture:
                glActiveTexture(GL_TEXTURE0);
                _gllib.glBindTexture(GL_TEXTURE_2D, 0)
                glDisable(obj.texture.dim)

        if lIsInstancetransformable:
            for c in obj.clipP: # disable object's clip planes
                c._Disable()
        
            if obj.scissor:
                glDisable(GL_SCISSOR_TEST)


    def Draw(self, obj, pushAttrib=True):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Draw obj and subtree below it
"""
        #print "Camera.Draw", obj

        assert obj.viewer is not None

        # skip object if it is in list of objects not shown in this Camera
        if obj.hiddenInCamera.has_key(self): return
        
        # if object is not visible, all subtree is hidden --> we can return
        if obj.getVisible() in [False, 0]: #not obj.visible:
            return

        if self.drawMode & 5: # immediateRendring draw mode
            if obj.immediateRendering:
                pass
            elif (self.viewer.tileRender and obj.needsRedoDpyListOnResize):
                pass
            else:
                return

        glPushMatrix()                        # Protect our matrix

        lIsInstanceTransformable = isinstance(obj, Transformable)
        if lIsInstanceTransformable: 
            if not obj.inheritXform:
                glPushMatrix()
                glLoadIdentity()
                glMultMatrixf(self.beforeRootTransformation)
            
            obj.MakeMat()

# VERY expensive and I am not sure what I need it for
        if pushAttrib:
            glPushAttrib(GL_CURRENT_BIT | GL_LIGHTING_BIT
                         | GL_POLYGON_BIT)
            glDisable(GL_LIGHTING)

# required for the transparent box of AutoGrid not to affect the molecule's
# color binding. The box is FLAT shaded and this shade models gets propagated
# to the molecule preventing color interpolation over bonds.
# This shows that thi models of inheritence has serious limitations :(
#        gl.glPushAttrib(GL_LIGHTING_BIT)

        if lIsInstanceTransformable:
            obj.SetupGL() # setup GL properties the can be inherited

            # enable and display object's clipping planes that also
            # clip its children
            self.ActivateClipPlanes( obj.clipPI, obj.clipSide )

            if obj.scissor:
                glEnable(GL_SCISSOR_TEST)
                glScissor(obj.scissorX, obj.scissorY, obj.scissorW, obj.scissorH)

        if self.drawMode == 2:
            # recusive call for children
            # for o in obj.children: self.Draw(o)
            mi = 0
            for m in obj.instanceMatricesFortran:
                if lIsInstanceTransformable and not obj.inheritMaterial:
                    obj.InitMaterial(mi)
                    obj.InitColor(mi)
                    mi = mi + 1
                glPushMatrix()
                glMultMatrixf(m)
                map ( self.Draw, obj.children)
                glPopMatrix()

        # Should be done in a function to setup GL properties that cannot
        # be inherited
        # should actually only be done once since all transparent objects
        # have to be drawn after all opaque objects for this to work
        draw = 1
        transp = obj.transparent #obj.isTransparent()
        if obj.immediateRendering or \
           (self.viewer.tileRender and obj.needsRedoDpyListOnResize):
            if transp:
                if not self.drawMode & 4: # only render if immediate & transp
                    draw = 0
            else:
                if not self.drawMode & 1: # only render if immediate & opaque
                    draw = 0
##                 if draw==1:
##                     print 'drawing immediate opaque object'

        if draw and obj.visible:
            if transp:
                if self.drawTransparentObjects:
                    glEnable(GL_BLEND)
                    if not obj.getDepthMask():
                        glDepthMask(GL_FALSE)
                    glBlendFunc(obj.srcBlendFunc, obj.dstBlendFunc)

                    self.DrawOneObject(obj)
                
                    glDisable(GL_BLEND)
                    glDepthMask(GL_TRUE)
                else:
                    self.hasTransparentObjects = 1
            else: # was: elif not self.drawTransparentObjects:
                self.DrawOneObject(obj)

# VERY Expensif
        if pushAttrib:
            glPopAttrib()

        if lIsInstanceTransformable:
            for c in obj.clipPI: # disable object's clip planes that are
                c._Disable()   # inherited by children

            if obj.scissor:
                glDisable(GL_SCISSOR_TEST)

            if not obj.inheritXform:
                glPopMatrix()

        glPopMatrix()                        # Restore the matrix


    def RedrawObjectHierarchy(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        
        glPushMatrix()
        
        # setup GL state for root object
        obj = self.viewer.rootObject
        
        if self.drawBB: obj.DrawTreeBoundingBox()
        if self.drawBB != viewerConst.ONLY:

            # save the transformation before the root object
            # used by some object who do not want to inherit transform
            m = numpy.array(glGetDoublev(GL_MODELVIEW_MATRIX)).astype('f')
            self.beforeRootTransformation = numpy.reshape( m, (16,) )

            obj.MakeMat()               # root object transformation
            # mod_opengltk
            #m = glGetDoublev(GL_MODELVIEW_MATRIX).astype('f')
            m = numpy.array(glGetDoublev(GL_MODELVIEW_MATRIX)).astype('f')

            self.rootObjectTransformation = numpy.reshape( m, (16,) )
# no reason to push and pop attribute of root object 
#            glPushAttrib(GL_CURRENT_BIT | GL_LIGHTING_BIT |
#                         GL_POLYGON_BIT)

            #obj.InitColor()       # init GL for color with ROOT color
            #obj.InitMaterial()    # init GL for material with ROOT material
            obj.SetupGL()         # setup GL with ROOT properties

            if len(obj.clipPI):
                self.ActivateClipPlanes( obj.clipPI, obj.clipSide )

            # draw all object that do not want a display List
            #print "self.viewer.noDpyListOpaque", self.viewer.noDpyListOpaque
            if len(self.viewer.noDpyListOpaque):
                self.drawMode = 1 # no dispaly list Opaque
                for m in obj.instanceMatricesFortran:
                    glPushMatrix()
                    glMultMatrixf(m)
                    #map( self.Draw, obj.children)
                    map( self.Draw, self.viewer.noDpyListOpaque )
                    glPopMatrix()

            # for "Continuity" (ucsd - nbcr) as we don't have vertexArrays yet
            if hasattr(self, 'vertexArrayCallback'):
                self.vertexArrayCallback()

            # when redraw is triggered by pick event self.viewer has been set
            # to None
            if self.dpyList and self.viewer.useMasterDpyList:
                #print 'calling master list'
                currentcontext = self.getContext()
                if currentcontext != self.dpyList[1]:
                    warnings.warn("""RedrawObjectHierarchy failed because the current context is the wrong one""")
                    #print "currentcontext != self.viewer.dpyList[1]", currentcontext, self.viewer.dpyList[1]
                else:
                    #print '#%d'%self.dpyList[0], currentcontext, "glCallList Camera2"
                    glCallList(self.dpyList[0])

            else:
                #print 'rebuilding master list'
                if self.viewer.useMasterDpyList:
                    lNewList = glGenLists(1)
                    #lNewList = self.newList
                    lContext = self.getContext()
                    #print "lNewList StandardCamera.RedrawObjectHierarchy", lNewList, lContext, self.name
                    self.dpyList = ( 
                          lNewList,
                          lContext
                         )

                    glNewList(self.dpyList[0], GL_COMPILE)
                    #print '+%d'%self.dpyList[0], lContext, "glNewList Camera"

                # call for each subtree with root in obj.children
                # for o in obj.children: self.Draw(o)
                self.drawMode = 2
                self.drawTransparentObjects = 0
                self.hasTransparentObjects = 0
                for m in obj.instanceMatricesFortran:
                    glPushMatrix()
                    glMultMatrixf(m)
                    map( self.Draw, obj.children)
                    glPopMatrix()
                    
                if self.hasTransparentObjects:
                    self.drawTransparentObjects = 1
                    for m in obj.instanceMatricesFortran:
                        glPushMatrix()
                        glMultMatrixf(m)
                        map( self.Draw, obj.children)
                        glPopMatrix()

                if self.viewer.useMasterDpyList:
                    #print '*%d'%GL.glGetIntegerv(GL.GL_LIST_INDEX), "glEndList Camera"
                    glEndList()

                    #print "self.viewer.dpyList", self.viewer.dpyList
                    if self.dpyList:
                        currentcontext = self.getContext()
                        if currentcontext != self.dpyList[1]:
                            warnings.warn("""RedrawObjectHierarchy failed because the current context is the wrong one""")
                            #print "currentcontext != self.viewer.dpyList[1]", currentcontext, self.viewer.dpyList[1]
                        else:
                            #print '#%d'%self.dpyList[0], currentcontext, "glCallList Camera3"
                            glCallList(self.dpyList[0])

            # draw all object that do not want a display List
            #print "self.viewer.noDpyListTransp", self.viewer.noDpyListTransp
            if len(self.viewer.noDpyListTransp):
                self.drawTransparentObjects = 1
                self.drawMode = 4
                for m in obj.instanceMatricesFortran:
                    glPushMatrix()
                    glMultMatrixf(m)
                    #map( self.Draw, obj.children)
                    map( self.Draw, self.viewer.noDpyListTransp )
                    glPopMatrix()
                self.drawTransparentObjects = 0

##            glPopAttrib()
            for c in obj.clipPI:
                c._Disable()

        glPopMatrix()


    def _RedrawCamera(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Actual drawing of lights, clipping planes and objects
"""
        if self.stereoMode.startswith('SIDE_BY_SIDE'):
            if self.viewer.tileRender:
                glPushMatrix()
                # setup GL state for camera
                self.BuildTransformation()  # camera transformation
                self.SetupLights() 
                self.DrawAllDirectionalLights()
                if self.imageRendered == 'RIGHT_EYE':
                    glTranslatef( float(-self.sideBySideTranslation), 0, 0 )
                    glRotatef( float(-self.sideBySideRotAngle), 0., 1., 0.)
                elif self.imageRendered == 'LEFT_EYE':
                    glTranslatef( float(self.sideBySideTranslation), 0, 0 )
                    glRotatef( float(self.sideBySideRotAngle), 0., 1., 0.)
                else:
                    assert False, 'self.imageRendered is not set correctly'
                glViewport(0, 0, int(self._width), int(self._height))
                glScalef( 1., 0.5, 1.)
                self.RedrawObjectHierarchy()
                glPopMatrix()

            elif self.stereoMode.endswith('_CROSS'):
                halfWidth = self._width/2
                
                # setup GL state for camera
                glPushMatrix()
                self.BuildTransformation()  # camera transformation
                self.SetupLights()  # set GL lights; moved from Light object to here
                self.DrawAllDirectionalLights()
                
                # render image for right eye
                self.imageRendered = 'RIGHT_EYE'
                glPushMatrix()
                glTranslatef( float(-self.sideBySideTranslation), 0, 0 )
                glRotatef( float(-self.sideBySideRotAngle), 0., 1., 0.)
                glViewport(0, 0, int(halfWidth), int(self._height) )
                glScalef( .5, 0.5, 1.)
                self.RedrawObjectHierarchy()
                glPopMatrix()

                # render image for left eye
                self.imageRendered = 'LEFT_EYE'
                glPushMatrix()
                glTranslatef( float(self.sideBySideTranslation), 0, 0 )
                glRotatef( float(self.sideBySideRotAngle), 0., 1., 0.)
                glScalef( .5, 0.5, 1.)
                glViewport(int(halfWidth), 0, int(halfWidth), int(self._height)) 
                self.RedrawObjectHierarchy()
                glPopMatrix()

                glPopMatrix()
                glViewport(0, 0, int(self._width), int(self._height))

            elif self.stereoMode.endswith('_STRAIGHT'):
                halfWidth = self._width/2

                # setup GL state for camera
                glPushMatrix()
                self.BuildTransformation()  # camera transformation
                self.SetupLights()  # set GL lights; moved from Light object to here
                self.DrawAllDirectionalLights()

                # render image for left eye
                self.imageRendered = 'LEFT_EYE'
                glPushMatrix()
                glTranslatef( float(self.sideBySideTranslation), 0, 0 )
                glRotatef( float(self.sideBySideRotAngle), 0., 1., 0.)
                glScalef( .5, 0.5, .25)
                glViewport(0, 0, halfWidth, self._height) 
                self.RedrawObjectHierarchy()
                glPopMatrix()

                # render image for right eye
                self.imageRendered = 'RIGHT_EYE'
                glPushMatrix()
                glTranslatef( float(-self.sideBySideTranslation), 0, 0 )
                glRotatef( float(-self.sideBySideRotAngle), 0., 1., 0.)
                glViewport(int(halfWidth), 0, int(halfWidth), int(self._height))
                glScalef( .5, 0.5, .25)
                self.RedrawObjectHierarchy()
                glPopMatrix()

                glPopMatrix()
                glViewport(0, 0, self._width, self._height)

        elif self.stereoMode=='3DTV':
            halfWidth = self._width/2

            # setup GL state for camera
            glPushMatrix()
            self.BuildTransformation()  # camera transformation
            self.SetupLights()  # set GL lights; moved from Light object to here
            self.DrawAllDirectionalLights()

            # render image for left eye
            self.imageRendered = 'LEFT_EYE'
            glPushMatrix()
            glTranslatef( float(self.TV3DTranslation), 0, 0 )
            glRotatef( float(self.TV3DRotAngle), 0., 1., 0.)
            sx, sy, sz = self.TV3DScaling
            glScalef( sx, sy, sz)
            glViewport(0, 0, halfWidth, self._height) 
            self.RedrawObjectHierarchy()
            glPopMatrix()

            # render image for right eye
            self.imageRendered = 'RIGHT_EYE'
            glPushMatrix()
            glTranslatef( float(-self.TV3DTranslation), 0, 0 )
            glRotatef( float(-self.TV3DRotAngle), 0., 1., 0.)
            glViewport(int(halfWidth), 0, int(halfWidth), int(self._height))
            glScalef( sx, sy, sz)
            self.RedrawObjectHierarchy()
            glPopMatrix()

            glPopMatrix()
            glViewport(0, 0, self._width, self._height)

        elif self.stereoMode.startswith('COLOR_SEPARATION'):
            if self.stereoMode.endswith('_RED_BLUE'):
                left = (GL_TRUE, GL_FALSE, GL_FALSE)
                right = (GL_FALSE, GL_FALSE, GL_TRUE)
            elif self.stereoMode.endswith('_BLUE_RED'):
                left = (GL_FALSE, GL_FALSE, GL_TRUE)
                right = (GL_TRUE, GL_FALSE, GL_FALSE)
            elif self.stereoMode.endswith('_RED_GREEN'):
                left = (GL_TRUE, GL_FALSE, GL_FALSE)
                right = (GL_FALSE, GL_TRUE, GL_FALSE)
            elif self.stereoMode.endswith('_GREEN_RED'):
                left = (GL_FALSE, GL_TRUE, GL_FALSE)
                right = (GL_TRUE, GL_FALSE, GL_FALSE)
            elif self.stereoMode.endswith('_RED_GREENBLUE'):
                left = (GL_TRUE, GL_FALSE, GL_FALSE)
                right = (GL_FALSE, GL_TRUE, GL_TRUE)
            elif self.stereoMode.endswith('_GREENBLUE_RED'):
                left = (GL_FALSE, GL_TRUE, GL_TRUE)
                right = (GL_TRUE, GL_FALSE, GL_FALSE)
            elif self.stereoMode.endswith('_REDGREEN_BLUE'):
                left = (GL_TRUE, GL_TRUE, GL_FALSE)
                right = (GL_FALSE, GL_FALSE, GL_TRUE)
            elif self.stereoMode.endswith('_BLUE_REDGREEN'):
                left = (GL_FALSE, GL_FALSE, GL_TRUE)
                right = (GL_TRUE, GL_TRUE, GL_FALSE)

            glPushMatrix()
            # setup GL state for camera
            self.BuildTransformation()  # camera transformation
            self.SetupLights()  # set GL lights; moved from Light object to here
            self.DrawAllDirectionalLights()

            # render image for left eye
            self.imageRendered = 'MONO' #'LEFT_EYE'
            glPushMatrix()
            glTranslatef( float(self.sideBySideTranslation), 0, 0 )
            glRotatef( float(self.sideBySideRotAngle), 0., 1., 0.)
            glViewport(0, 0, int(self._width), int(self._height)) 
            glColorMask(int(left[0]), int(left[1]), int(left[2]), GL_FALSE)
            self.RedrawObjectHierarchy()
            glPopMatrix()
            self.drawHighlight()
            glClear(GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT)

            # render image for right eye
            self.imageRendered = 'RIGHT_EYE'
            glPushMatrix()
            glTranslatef( float(-self.sideBySideTranslation), 0, 0 )
            glRotatef( float(-self.sideBySideRotAngle), 0., 1., 0.)
            glViewport(0, 0, int(self._width), int(self._height))
            glColorMask(int(right[0]),int(right[1]),int(right[2]), GL_FALSE)
            self.RedrawObjectHierarchy()
            glPopMatrix()
            self.drawHighlight()

            glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_FALSE)
            glPopMatrix()

        elif self.stereoMode == 'STEREO_BUFFERS':
            glPushMatrix()
            # setup GL state for camera
            self.BuildTransformation()  # camera transformation
            self.SetupLights()  # set GL lights; moved from Light object to here
            self.DrawAllDirectionalLights()

            # render image for left eye
            self.imageRendered = 'LEFT_EYE'
            glPushMatrix()
            glTranslatef( float(self.sideBySideTranslation), 0, 0 )
            glRotatef( float(self.sideBySideRotAngle), 0., 1., 0.)
            glViewport(0, 0, int(self._width), int(self._height)) 
            self.RedrawObjectHierarchy()
            glPopMatrix()
            self.drawHighlight()

            # render image for right eye
            self.imageRendered = 'RIGHT_EYE'
            glDrawBuffer(GL_BACK_RIGHT)
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT |
                    GL_STENCIL_BUFFER_BIT)
            glPushMatrix()
            glTranslatef( float(-self.sideBySideTranslation), 0, 0 )
            glRotatef( float(-self.sideBySideRotAngle), 0., 1., 0.)
            glViewport(0, 0, int(self._width), int(self._height)) 
            self.RedrawObjectHierarchy()
            glPopMatrix()
            self.drawHighlight()

            glPopMatrix()
            glDrawBuffer(GL_BACK_LEFT)

        else: # default to "Mono mode"
            self.imageRendered = 'MONO'
            glPushMatrix()
            # setup GL state for camera
            self.BuildTransformation()  # camera transformation
            self.SetupLights() 
            self.DrawAllDirectionalLights()
            glViewport(0, 0, int(self._width), int(self._height))
            self.RedrawObjectHierarchy()
            glPopMatrix()


#    def secondDerivative(self, im, width, height):
#        return im.filter(sndDerivK)

    def firstDerivative(self, im, width, height):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        fstDeriveV1K = ImageFilter.Kernel( (3,3), fstDeriveV1, self.d1scale)
        c = im.filter(fstDeriveV1K)

        fstDeriveV2K = ImageFilter.Kernel( (3,3), fstDeriveV2, self.d1scale)
        d = im.filter(fstDeriveV2K)

        fstDeriveH1K = ImageFilter.Kernel( (3,3), fstDeriveH1, self.d1scale)
        e = im.filter(fstDeriveH1K)

        fstDeriveH2K = ImageFilter.Kernel( (3,3), fstDeriveH2, self.d1scale)
        f = im.filter(fstDeriveH2K)

        result1 = ImageChops.add(c, d, 2.)
        result2 = ImageChops.add(e, f, 2.)
        result = ImageChops.add(result1, result2, 2.)
        return result

        ## alternative adding numeric arrays seems slower
##         t1 = time()
##         c = im.filter(fstDeriveV1K)
##         c = numpy.fromstring(c.tostring(), numpy.uint8)

##         d = im.filter(fstDeriveV2K)
##         d = numpy.fromstring(d.tostring(), numpy.uint8)

##         e = im.filter(fstDeriveH1K)
##         e = numpy.fromstring(e.tostring(), numpy.uint8)

##         f = im.filter(fstDeriveH2K)
##         f = numpy.fromstring(f.tostring(), numpy.uint8)
        
##         result = numpy.fabs(c) + numpy.fabs(d) +\
##                  numpy.fabs(e) + numpy.fabs(f)
##         result.shape = im.size
##         result = numpy.array( result )*0.25
##         return result
    
    
        
        
    
    def drawNPR(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        t1 = time()

        if self.viewer.tileRender is False:
            zbuf = self.GrabZBuffer(lock=False, flipTopBottom=False)
            #imageZbuf = Image.fromstring('L', (self._width, self._height), zbuf.tostring() )
            #imageZbuf.save("zbuf.png")
        else:
            zbuf = self.GrabZBuffer(lock=False, flipTopBottom=False,
                                    zmin=self.viewer.tileRenderCtx.zmin,
                                    zmax=self.viewer.tileRenderCtx.zmax
                                    )
            imageFinal = Image.new('L', (self._width+2, self._height+2) )

            #padding
            # the 4 corners
            imageFinal.putpixel((0, 0), zbuf.getpixel((0, 0)) )
            imageFinal.putpixel((0, self._height+1), zbuf.getpixel((0, self._height-1)) )
            imageFinal.putpixel((self._width+1, 0), zbuf.getpixel((self._width-1, 0)) )
            imageFinal.putpixel((self._width+1, self._height+1), zbuf.getpixel((self._width-1, self._height-1)) )

            # the top and bottom line
            for i in range(self._width):
                imageFinal.putpixel((i+1, 0), zbuf.getpixel((i, 0)) )
                imageFinal.putpixel((i+1, self._height+1), zbuf.getpixel((i, self._height-1)) )

            # the left and right columns
            for j in range(self._height):
                imageFinal.putpixel((0, j+1), zbuf.getpixel((0, j)) )
                imageFinal.putpixel((self._width+1, j+1), zbuf.getpixel((self._width-1, j)) )

            # the main picture
            imageFinal.paste( zbuf, (1 , 1) )
            zbuf = imageFinal

        d1strong = d2strong = None

        useFirstDerivative = self.d1scale!=0.0
        useSecondDerivative = self.d2scale!=0.0

        if useFirstDerivative:

            if self.viewer.tileRender:
                d1 = self.firstDerivative(zbuf, self._width+2, self._height+2)
                d1 = d1.crop((1,1,self._width+1,self._width+1))
            else:
                d1 = self.firstDerivative(zbuf, self._width, self._height)

            #d1.save('first.jpg')
            #d1strong = numpy.fromstring(d1.tostring(),numpy.uint8)
            #print 'first', min(d1.ravel()), max(d1.ravel()), numpy.sum(d1.ravel())
            #print self.d1cut, self.d1off
            d1 = numpy.fromstring(d1.tostring(),numpy.uint8)
            #mini = min(d1)
            #maxi = max(d1)
            #print 'd1',mini, maxi
#            d1strong = numpy.clip(d1, self.d1cutL, self.d1cutH)
#            d1strong = d1strong.astype('f')*self.d1off

            #print 'd1',min(d1strong), max(d1strong)
            #print 'first', d1strong.shape, min(d1strong.ravel()), max(d1strong.ravel()), numpy.sum(d1strong.ravel())

            # LOOKUP ramp
            #print "self.d1ramp", len(self.d1ramp), self.d1ramp

            #d1strong = numpy.choose(d1.astype('B'), self.d1ramp )
            d1strong = numpy.zeros(len(d1))
            for index, elt in enumerate(d1.astype('B')):
                 d1strong[index] = self.d1ramp[elt]

            #print 'firstramp', d1strong.shape, min(d1strong.ravel()), max(d1strong.ravel()), numpy.sum(d1strong.ravel())

        if useSecondDerivative:
            #d2 = self.secondDerivative(zbuf, self._width, self._height)
            sndDerivK = ImageFilter.Kernel( (3,3), sndDeriv, self.d2scale)
            d2 = zbuf.filter(sndDerivK)

            if self.viewer.tileRender:
                d2 = d2.crop((1,1,self._width+1,self._width+1))            

            #d2.save('second.jpg')
            #print 'second1', min(d2.ravel()), max(d2.ravel()), numpy.sum(d2.ravel())
            #print self.d2cut, self.d2off
            d2 = numpy.fromstring(d2.tostring(),numpy.uint8)
            #mini = min(d2)
            #maxi = max(d2)
            #print 'd2',mini, maxi
            d2strong = numpy.clip(d2, self.d2cutL, self.d2cutH)
            d2strong = d2strong.astype('f')*self.d2off
            #d2strong = numpy.where(numpy.greater(d2,self.d2cutL),
            #                         d2, 0)
            #d2strong = numpy.where(numpy.less(d2strong,self.d2cutH),
            #                         d2, 0)
            #d2strong += self.d2off
            #print 'd2',min(d2strong), max(d2strong)
            #print 'second2', d2strong.shape, min(d2strong.ravel()), max(d2strong.ravel()), numpy.sum(d2strong.ravel())

        if useFirstDerivative and useSecondDerivative:
            self.outline = numpy.maximum(d1strong, d2strong)
            #self.outline = (d1strong + d2strong)/2
            #self.outlineim = ImageChops.add(d1strong, d2strong)
        elif useFirstDerivative:
            self.outline = d1strong
        elif useSecondDerivative:
            self.outline = d2strong
        else:
            self.outline = None

        ## working OpenGL version
##         zbuf = self.GrabZBuffer(lock=False)
##         self.ivi.setImage(zbuf)

##         deriv = self.ivi.firstDerivative()
##         d1strong = numpy.where(numpy.greater(deriv,self.d1cut),
##                                  self.d1scale*(deriv+self.d1off), 0)

##         deriv2 = self.ivi.secondDerivative()
##         d2strong = numpy.where(numpy.greater(deriv2,self.d2cut),
##                                  self.d2scale*(deriv2+self.d2off), 0)
##         self.outline = numpy.maximum(d1strong, d2strong)
##         self.tk.call(self._w, 'makecurrent')

        #print 'time NPR rendering', time()-t1
        return


    def displayNPR(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        t1 = time()
        if self.outline is None:
            return

        lViewport = glGetIntegerv(GL_VIEWPORT)
        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()
        glOrtho(float(lViewport[0]),
                float(lViewport[0]+lViewport[2]),
                float(lViewport[1]),
                float(lViewport[1]+lViewport[3]),
                -1, 1)
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glLoadIdentity()
        glDisable(GL_DEPTH_TEST)
        glEnable(GL_BLEND)
        glBlendFunc(GL_DST_COLOR, GL_ZERO);  
        glRasterPos2f(0.0, 0.0)

        self.outline = numpy.fabs(self.outline-255)
        self.outline = self.outline.astype('B')
#        glActiveTexture(GL_TEXTURE0);
        GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1)
        _gllib.glDrawPixels(self._width, self._height, 
                            GL_LUMINANCE, GL_UNSIGNED_BYTE, 
                            self.outline )
        
        glMatrixMode(GL_PROJECTION)
        glPopMatrix()
        glMatrixMode(GL_MODELVIEW)
        glPopMatrix()

        #print 'time NPR display', time()-t1

        
    def RedrawAASwitch(self, *dummy):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Redraw all the objects in the scene"""

        glStencilMask(1)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT |
                GL_STENCIL_BUFFER_BIT)
        glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE)
        glStencilFunc ( GL_ALWAYS, 0, 1 ) ;
        glEnable ( GL_STENCIL_TEST ) ;

        if not self.viewer.accumBuffersError and \
               self.antiAliased and self.renderMode==GL_RENDER:
            dist = math.sqrt(numpy.add.reduce(self.direction*self.direction))
            # jitter loop
            if self.antiAliased>0:
                sca = 1./self.antiAliased
                accumContour = None

            for i in range(self.antiAliased):
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

                if self.projectionType == self.PERSPECTIVE:
                    if self.viewer.tileRender:
                        tr = self.viewer.tileRenderCtx
                        left = tr.tileleft
                        right = tr.tileright
                        bottom = tr.tilebottom
                        top = tr.tiletop
                    else:
                        left=right=top=bottom=None

                    self.AccPerspective (self.jitter[i][0], self.jitter[i][1],
                                      0.0, 0.0, dist, left, right, bottom, top)
                    self._RedrawCamera()
                else:
                    W = self.right-self.left # width in world coordinates
                    H = self.top-self.bottom # height in world coordinates
                    glPushMatrix()
                    glTranslatef (self.jitter[i][0]*W/self._width,
                                  self.jitter[i][1]*H/self._height,
                                  0.0)
                    self._RedrawCamera()
                    glPopMatrix()

                if i == 0 and self.ssao :
                    self.copyDepth()

                # we accumulate the back buffer
                glReadBuffer(GL_BACK)
                if i==0:
                    glAccum(GL_LOAD, self.accumWeigth)
                else:
                    glAccum(GL_ACCUM, self.accumWeigth)
                    
                if self.contours and self._suspendNPR is False:
                    self.drawNPR()
                    if self.outline is not None:
                        if accumContour is None:
                            accumContour = self.outline.copy()
                            accumContour *= sca
                        else:
                            accumContour += sca*self.outline
            print 'GL ACCUM'
            glAccum (GL_RETURN, 1.0)

            if self.contours and self._suspendNPR is False:
                self.outline = accumContour

        else: # plain drawing
            if self.ssao : self.copy_depth_ssao = True
            self._RedrawCamera()
    
        if self.contours and self._suspendNPR is False:
            self.drawNPR()


        if self.contours and self._suspendNPR is False:
            glViewport(0, 0, self._width, self._height)
            self.displayNPR()

        if    self.stereoMode.startswith('COLOR_SEPARATION') is False \
          and self.stereoMode != 'STEREO_BUFFERS':     
            self.drawHighlight()

        
        if self.ssao :
            self.drawSSAO()    

        glDisable(GL_STENCIL_TEST)            
        if self.drawThumbnailFlag:
            self.drawThumbnail()

        if self.swap:
            self.SwapBuffers()


    def setShaderSelectionContour(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        #import pdb;pdb.set_trace()
        #if DejaVu2.enableSelectionContour is True:
        #    if GL.glGetString(GL.GL_VENDOR).find('Intel') >= 0:
        #        DejaVu2.enableSelectionContour = False
        #        print "intel (even gma) gpu drivers don't handle properly the stencil with FBO"

        if DejaVu2.enableSelectionContour is True:
            try:
                ## do not try to import this before creating the opengl context, it would fail.
                from opengltk.extent import _glextlib
                DejaVu2.enableSelectionContour = True
            except ImportError:
                DejaVu2.enableSelectionContour = False
                print "could not import _glextlib"

        #import pdb;pdb.set_trace()
        if DejaVu2.enableSelectionContour is True:
            extensionsList = glGetString(GL_EXTENSIONS)
            if extensionsList.find('GL_EXT_packed_depth_stencil') < 0:
                DejaVu2.enableSelectionContour = False
                print "opengl extension GL_EXT_packed_depth_stencil is not present"
        # starting with swig 3.#.# GL_FRAGMENT_SHADER and other GL_###_###
        # constants are imported from opengltk/extent/glextlib.py (from opengltk.extent.glextlib import *)
        if DejaVu2.enableSelectionContour is True:
            f = _glextlib.glCreateShader(GL_FRAGMENT_SHADER)
    
            # This shader performs a 9-tap Laplacian edge detection filter. 
            # (converted from the separate "edges.cg" file to embedded GLSL string)    
            self.fragmentShaderCode = """
uniform float contourSize;
uniform vec4 contourColor;
uniform sampler2D texUnit;
void main(void)
{
    float lOffset = .001 * contourSize ; // (1./512.);
    vec2 texCoord = gl_TexCoord[0].xy;
    if (     ( texCoord [ 0 ] > lOffset ) // to suppress artifacts on frame buffer bounds
          && ( texCoord [ 1 ] > lOffset ) 
          && ( texCoord [ 0 ] < 1. - lOffset )
          && ( texCoord [ 1 ] < 1. - lOffset ) )
    {
        float c  = texture2D(texUnit, texCoord)[0];
        float bl = texture2D(texUnit, texCoord + vec2(-lOffset, -lOffset))[0];
        float l  = texture2D(texUnit, texCoord + vec2(-lOffset,     0.0))[0];
        float tl = texture2D(texUnit, texCoord + vec2(-lOffset,  lOffset))[0];
        float t  = texture2D(texUnit, texCoord + vec2(    0.0,  lOffset))[0];
        float tr = texture2D(texUnit, texCoord + vec2( lOffset,  lOffset))[0];
        float r  = texture2D(texUnit, texCoord + vec2( lOffset,     0.0))[0];
        float br = texture2D(texUnit, texCoord + vec2( lOffset,  lOffset))[0];
        float b  = texture2D(texUnit, texCoord + vec2(    0.0, -lOffset))[0];
        if ( 8. * (c + -.125 * (bl + l + tl + t + tr + r + br + b)) != 0. )
        {
           gl_FragColor = contourColor; //vec4(1., 0., 1., .7) ;
        }
        else
        {
            //due to bizarre ATI behavior
            gl_FragColor = vec4(1., 0., 0., 0.) ;
        }
    }
}
"""

            _glextlib.glShaderSource(f, 1, self.fragmentShaderCode, 0x7FFFFFFF)
            _glextlib.glCompileShader(f)
            lStatus = 0x7FFFFFFF
            lStatus = _glextlib.glGetShaderiv(f, GL_COMPILE_STATUS, lStatus )
            if lStatus == 0:
                print "compile status", lStatus
                charsWritten  = 0
                shaderInfoLog = '\0' * 2048
                charsWritten, infoLog = _glextlib.glGetShaderInfoLog(f, len(shaderInfoLog), 
                                                                     charsWritten, shaderInfoLog)
                print "shaderInfoLog", shaderInfoLog
                DejaVu2.enableSelectionContour = False
                print "selection contour shader didn't compile"
            else:
                self.shaderProgram = _glextlib.glCreateProgram()
                _glextlib.glAttachShader(self.shaderProgram,f)
                _glextlib.glLinkProgram(self.shaderProgram)
                lStatus = 0x7FFFFFFF
                lStatus = _glextlib.glGetProgramiv(self.shaderProgram,
                                                   GL_LINK_STATUS,
                                                   lStatus )
                if lStatus == 0:
                    print "link status", lStatus
                    DejaVu2.enableSelectionContour = False
                    print "selection contour shader didn't link"
                else:
                    _glextlib.glValidateProgram(self.shaderProgram)
                    lStatus = 0x7FFFFFFF
                    lStatus = _glextlib.glGetProgramiv(self.shaderProgram,
                                                       GL_VALIDATE_STATUS,
                                                       lStatus )
                    if lStatus == 0:
                        print "validate status", lStatus
                        DejaVu2.enableSelectionContour = False
                        print "selection contour shader not validated"
                    else:
                        # Get location of the sampler uniform
                        self.texUnitLoc = int(_glextlib.glGetUniformLocation(self.shaderProgram, "texUnit"))
                        #print "selection contour shader successfully compiled and linked"

                        # Get location of the contourSize uniform
                        self.contourSizeLoc = int(_glextlib.glGetUniformLocation(self.shaderProgram, "contourSize"))

                        # Get location of the contourSize uniform
                        self.contourColorLoc = int(_glextlib.glGetUniformLocation(self.shaderProgram, "contourColor"))


                        ## create a framebuffer to receive the texture
                        self.fb = 0
                        self.fb = _glextlib.glGenFramebuffersEXT(1, self.fb)
                        #print "self.fb", self.fb
                        lCheckMessage = _glextlib.glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT)
                        if lCheckMessage != GL_FRAMEBUFFER_COMPLETE_EXT: # 0x8CD5
                            print 'glCheckFramebufferStatusEXT %x'%lCheckMessage
                            DejaVu2.enableSelectionContour = False


    def drawHighlight(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
#        lStencil = numpy.zeros(self._width*self._height, dtype='uint32')
#        GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1)
#        _gllib.glReadPixels( 0, 0, self._width, self._height,
#                             GL.GL_STENCIL_INDEX, GL.GL_UNSIGNED_INT,
#                             lStencil)
#        print "lStencil", lStencil[0] # background is 254 # fill is 255
#        lStencilImage = Image.fromstring('L', (self._width, self._height), 
#                                       lStencil.astype('uint8').tostring() )
#        lStencilImage.save("lStencil.png")
#        zbuf = self.GrabZBuffer(lock=False, flipTopBottom=False)
#        imageZbuf = Image.fromstring('L', (self._width, self._height), zbuf.tostring() )
#        imageZbuf.save("lZbuff.png")

        # to draw only the highlight zone from the stencil buffer
        GL.glEnable(GL.GL_STENCIL_TEST)
        glStencilMask(0)
        glStencilFunc(GL_EQUAL, 1, 1)
        glPushMatrix()
        glLoadIdentity()
        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()
        lViewport = glGetIntegerv(GL_VIEWPORT)
        glOrtho(float(lViewport[0]),
                float(lViewport[0]+lViewport[2]),
                float(lViewport[1]),
                float(lViewport[1]+lViewport[3]),
                -1, 1)
        glDisable(GL_DEPTH_TEST)
        glDisable(GL_LIGHTING)
#        glActiveTexture(GL_TEXTURE0);
        ## here we apply the patterned highlight
        if DejaVu2.selectionPatternSize >= 3:
            #glEnable(GL_COLOR_LOGIC_OP)
            #glLogicOp(GL_AND)
            glEnable(GL_TEXTURE_2D)
            _gllib.glBindTexture(GL_TEXTURE_2D, self.textureName)
            glEnable(GL_BLEND)
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
            glBegin(GL_POLYGON)
            glTexCoord2f(0., 0.)
            glVertex2i(0, 0)
            lTexWidth = float(self._width)/DejaVu2.selectionPatternSize # determine how big the pattern will be 
            lTexHeight = float(self._height)/DejaVu2.selectionPatternSize
            glTexCoord2f(lTexWidth, 0.)
            glVertex2i(self._width, 0)
            glTexCoord2f(lTexWidth, lTexHeight)
            glVertex2i(self._width, self._height)
            glTexCoord2f(0., lTexHeight)
            glVertex2i(0, self._height)
            glEnd()
            glDisable(GL_BLEND)
            glDisable(GL_TEXTURE_2D)
            #glDisable(GL_COLOR_LOGIC_OP)

        if    DejaVu2.enableSelectionContour is True \
          and DejaVu2.selectionContourSize != 0:
            from opengltk.extent import _glextlib
            #import pdb;pdb.set_trace()
            ## copying the current stencil to a texture
            _gllib.glBindTexture(GL_TEXTURE_2D, self.stencilTextureName )
            glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_STENCIL_EXT, 0, 0,
                             self._width, self._height, 0)
            ## switch writting to the FBO (Frame Buffer Object)
            _glextlib.glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, self.fb )
            #lCheckMessage = _glextlib.glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT)
            #print 'glCheckFramebufferStatusEXT 0 %x'%lCheckMessage
            ## attaching the current stencil to the FBO (Frame Buffer Object)
            _glextlib.glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,
                                         GL_DEPTH_ATTACHMENT_EXT,
                                         GL_TEXTURE_2D, self.stencilTextureName,0)
            #lCheckMessage = _glextlib.glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT)
            #print 'glCheckFramebufferStatusEXT 1 %x'%lCheckMessage
            _glextlib.glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,
                                         GL_STENCIL_ATTACHMENT_EXT,
                                         GL_TEXTURE_2D, self.stencilTextureName, 0)
            #lCheckMessage = _glextlib.glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT)
            #print 'glCheckFramebufferStatusEXT 2 %x'%lCheckMessage
            ## attaching the texture to be written as FBO
            _gllib.glBindTexture(GL_TEXTURE_2D, self.contourTextureName )
            _gllib.glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, self._width, self._height, 
                                   0, GL_RGBA, GL.GL_FLOAT, 0 )
            _glextlib.glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,
                                                GL_COLOR_ATTACHMENT0_EXT,
                                                GL_TEXTURE_2D, self.contourTextureName, 0)
            #lCheckMessage = _glextlib.glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT)
            #print 'glCheckFramebufferStatusEXT 3 %x'%lCheckMessage

            ## verifying that everything gets attached to the FBO
            lCheckMessage = _glextlib.glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT)
            if lCheckMessage != GL_FRAMEBUFFER_COMPLETE_EXT: # 0x8CD5
                print 'glCheckFramebufferStatusEXT %x'%lCheckMessage
                DejaVu2.enableSelectionContour = False
                _glextlib.glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0 )
                print 'opengl frame buffer object not available, selection contour is now disabled.'
                #print 'you may need to set enableSelectionContour to False in the following file:'
                #print '~/.mgltools/(version numer)/DejaVu2/_dejavurc'
            else:
                ## writing the stencil to the texture / FBO 
                glClear(GL_COLOR_BUFFER_BIT)
                glBegin(GL_POLYGON)
                glVertex2i(0, 0)
                glVertex2i(self._width, 0)
                glVertex2i(self._width, self._height)
                glVertex2i(0, self._height)
                glEnd()

                ## switch writing to the regular frame buffer (normaly the back buffer)
                _glextlib.glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0 )
                
                ## here we obtain the contour of the stencil copied in the texture and draw it
                from opengltk.extent import _glextlib
                _glextlib.glUseProgram( self.shaderProgram )
                _glextlib.glUniform1i( self.texUnitLoc, 0)
                _glextlib.glUniform1f( self.contourSizeLoc, float(DejaVu2.selectionContourSize))
                _glextlib.glUniform4f( self.contourColorLoc, 
                                       DejaVu2.selectionContourColor[0],
                                       DejaVu2.selectionContourColor[1],
                                       DejaVu2.selectionContourColor[2],
                                       DejaVu2.selectionContourColor[3]
                                     )                
                glEnable(GL_BLEND)
                glStencilFunc(GL_NOTEQUAL, 1, 1) #so the contour is drawn only outside of the highlighted area
                glBegin(GL_POLYGON)
                glTexCoord2i(0, 0)
                glVertex2i(0, 0)
                glTexCoord2i(1, 0)
                glVertex2i(self._width, 0)
                glTexCoord2i(1, 1)
                glVertex2i(self._width, self._height)
                glTexCoord2i(0, 1)
                glVertex2i(0, self._height)
                glEnd()
                glDisable(GL_BLEND)
            _glextlib.glUseProgram( 0 )

        # to protect our texture, we bind the default texture while we don't use it
        _gllib.glBindTexture(GL_TEXTURE_2D, 0)      
        glEnable(GL_DEPTH_TEST)
        glPopMatrix()
        glMatrixMode(GL_MODELVIEW)
        glPopMatrix()
        glStencilMask(1)
        glStencilFunc ( GL_ALWAYS, 0, 1 )

#===============================================================================
#     SSAO Shader
#===============================================================================
    def setTextureSSAO(self):
        self.depthtextureName = int(glGenTextures(1)[0])
        glPrioritizeTextures(numpy.array([self.depthtextureName]),
                             numpy.array([1.]))
        _gllib.glBindTexture (GL_TEXTURE_2D, self.depthtextureName);
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri (GL_TEXTURE_2D, GL_DEPTH_TEXTURE_MODE, GL_LUMINANCE);
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_NONE )
        _gllib.glTexImage2D (GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT32, self._width,
                      self._height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, None);
        #SSAO depthexture
        self.illumtextureName = int(glGenTextures(1)[0])
        glPrioritizeTextures(numpy.array([self.illumtextureName]),
                             numpy.array([1.]))
        _gllib.glBindTexture (GL_TEXTURE_2D, self.illumtextureName);
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri (GL_TEXTURE_2D, GL_DEPTH_TEXTURE_MODE, GL_LUMINANCE);
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_NONE )
        _gllib.glTexImage2D (GL_TEXTURE_2D, 0, GL_LUMINANCE, self._width,
                      self._height, 0, GL_LUMINANCE, GL_UNSIGNED_INT, None);

        #SSAO depthexture
        self.rendertextureName = int(glGenTextures(1)[0])
        glPrioritizeTextures(numpy.array([self.rendertextureName]),
                             numpy.array([1.]))
        _gllib.glBindTexture (GL_TEXTURE_2D, self.rendertextureName);
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_NONE )
        _gllib.glTexImage2D (GL_TEXTURE_2D, 0, GL_RGB, self._width,
                      self._height, 0, GL_RGB, GL_UNSIGNED_INT, None);
        #depth mask texture
        self.depthmasktextureName = int(glGenTextures(1)[0])
        glPrioritizeTextures(numpy.array([self.depthmasktextureName]),
                             numpy.array([1.]))
        _gllib.glBindTexture (GL_TEXTURE_2D, self.depthmasktextureName);
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
        glTexParameteri (GL_TEXTURE_2D, GL_DEPTH_TEXTURE_MODE, GL_LUMINANCE);
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_NONE )
        _gllib.glTexImage2D (GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT32, self._width,
                      self._height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, None);
        #SSAO randomtexture
#
#        self.randomTexture = Texture()
#        from PIL import Image
#        img = Image.open(DejaVu2.__path__[0]+os.sep+"textures"+os.sep+"random_normals.png")
#        self.randomTexture.Set(enable=1, image=img)
        self.randomtextureName = int(glGenTextures(1)[0])
        glPrioritizeTextures(numpy.array([self.randomtextureName]),
                             numpy.array([1.]))
#        _gllib.glBindTexture (GL_TEXTURE_2D, self.randomtextureName);       
#        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
#        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
#        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
#        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
#        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_NONE )
#        _gllib.glTexImage2D (GL_TEXTURE_2D, 0, self.randomTexture.format, 
#                             self.randomTexture._width, self.randomTexture._height, 
#                             0, self.randomTexture.format, GL_UNSIGNED_INT, 
#                             self.randomTexture.image);       
        _gllib.glBindTexture (GL_TEXTURE_2D,0)
        
    def setShaderSSAO(self):
        if __debug__:
            if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()        
        if not DejaVu2.enableSSAO :
            return
        try:
            ## do not try to import this before creating the opengl context, it would fail.
            from opengltk.extent import _glextlib
            DejaVu2.enableSSAO = True
        except ImportError:
            DejaVu2.enableSSAO = False
            print "could not import _glextlib"  
        #import pdb;pdb.set_trace()
        if DejaVu2.enableSSAO is True:
            extensionsList = glGetString(GL_EXTENSIONS)
            if extensionsList.find('GL_EXT_packed_depth_stencil') < 0:
                DejaVu2.enableSSAO = False
                print "opengl extension GL_EXT_packed_depth_stencil is not present"

        if DejaVu2.enableSSAO is True:
#            self.setDefaultSSAO_OPTIONS()
            self.setTextureSSAO()
            f = _glextlib.glCreateShader(GL_FRAGMENT_SHADER)
            sfile = open(DejaVu2.__path__[0]+os.sep+"shaders"+os.sep+"fragSSAO","r")
            lines = sfile.readlines()
            sfile.close()
            self.fragmentSSAOShaderCode=""
            for l in lines :
                self.fragmentSSAOShaderCode+=l
            v = _glextlib.glCreateShader(GL_VERTEX_SHADER)
            sfile = open(DejaVu2.__path__[0]+os.sep+"shaders"+os.sep+"vertSSAO","r")
            lines = sfile.readlines()
            sfile.close()
            self.vertexSSAOShaderCode=""
            for l in lines :
                self.vertexSSAOShaderCode+=l
            lStatus = 0x7FFFFFFF
            _glextlib.glShaderSource(v, 1, self.vertexSSAOShaderCode, lStatus)
            _glextlib.glCompileShader(v)
            _glextlib.glShaderSource(f, 1, self.fragmentSSAOShaderCode, lStatus)
            _glextlib.glCompileShader(f)

            lStatus1 = _glextlib.glGetShaderiv(f, GL_COMPILE_STATUS, lStatus )
            #print "SSAO compile status frag: %s"%lStatus,
            lStatus2 = _glextlib.glGetShaderiv(v, GL_COMPILE_STATUS, lStatus )
            #print " vertex: %s"%lStatus
            if lStatus1 == 0 or lStatus2 == 0:
#                print "compile status", lStatus
                charsWritten  = 0
                shaderInfoLog = '\0' * 2048
                charsWritten, infoLog = _glextlib.glGetShaderInfoLog(f, len(shaderInfoLog), 
                                                                     charsWritten, shaderInfoLog)
                print "shaderInfoLog", shaderInfoLog
                DejaVu2.enableSSAO = False
                print "SSAO shader didn't compile"
            else:
                self.shaderSSAOProgram = _glextlib.glCreateProgram()
                _glextlib.glAttachShader(self.shaderSSAOProgram,v)
                _glextlib.glAttachShader(self.shaderSSAOProgram,f)                
                _glextlib.glLinkProgram(self.shaderSSAOProgram)
                lStatus = 0x7FFFFFFF
                lStatus = _glextlib.glGetProgramiv(self.shaderSSAOProgram,
                                                   GL_LINK_STATUS,
                                                   lStatus )
                if lStatus == 0:
#                    print "link status", lStatus
                    log = ""
                    charsWritten  = 0
                    progInfoLog = '\0' * 2048
                    charsWritten, infoLog = _glextlib.glGetProgramInfoLog(self.shaderSSAOProgram, len(progInfoLog), 
                                                                     charsWritten, progInfoLog)
                    DejaVu2.enableSSAO = False
                    print "SSAO shader didn't link"
                    print "log ",progInfoLog
    
                else:
                    self.SSAO_LOCATIONS = {
                'RandomTexture': _glextlib.glGetUniformLocation( self.shaderSSAOProgram, 'RandomTexture' ),
                'DepthTexture': _glextlib.glGetUniformLocation( self.shaderSSAOProgram, 'DepthTexture' ),
                'RenderedTexture': _glextlib.glGetUniformLocation( self.shaderSSAOProgram, 'RenderedTexture' ),
                'LuminanceTexture': _glextlib.glGetUniformLocation( self.shaderSSAOProgram, 'LuminanceTexture' ),                
                'DepthMaskTexture': _glextlib.glGetUniformLocation( self.shaderSSAOProgram, 'DepthMaskTexture' ),
                'RenderedTextureWidth': _glextlib.glGetUniformLocation( self.shaderSSAOProgram, 'RenderedTextureWidth' ),
                'RenderedTextureHeight': _glextlib.glGetUniformLocation( self.shaderSSAOProgram, 'RenderedTextureHeight' ),                
                'realfar': _glextlib.glGetUniformLocation( self.shaderSSAOProgram, 'realfar' ),
                #'realnear': _glextlib.glGetUniformLocation( self.shaderSSAOProgram, 'realnear' ),
#                'fogS': _glextlib.glGetUniformLocation( self.shaderSSAOProgram, 'fogS' ),
#                'fogE': _glextlib.glGetUniformLocation( self.shaderSSAOProgram, 'fogE' ),
                    }
                        
                    for k in self.SSAO_OPTIONS :
                            self.SSAO_LOCATIONS[k] = _glextlib.glGetUniformLocation( self.shaderSSAOProgram, k )
                    _glextlib.glValidateProgram(self.shaderSSAOProgram)
                    lStatus = 0x7FFFFFFF
                    lStatus = _glextlib.glGetProgramiv(self.shaderSSAOProgram,
                                                       GL_VALIDATE_STATUS,
                                                       lStatus )
                    if lStatus == 0:
                        print "SSAO shader did not validate, status:", lStatus
#                        DejaVu2.enableSSAO = False
                        print "enableSSAO not validated"
#                    else:
                        # Get location 


    def setDefaultSSAO_OPTIONS(self):
        from opengltk.extent import _glextlib
	#can be Intel, NVIDIA
	
        self.SSAO_OPTIONS={'far': [300.0,0.,1000.,"float"],
                'near': [self.near,0.,100.,"float"],
                'method':[self.ssao_method,0,1,"int"],
                'do_noise':[0,0,1,"int"],
                'fog':[0,0,1,"int"],
                'use_fog':[1,0,1,"int"],
                'use_mask_depth':[self.use_mask_depth,0,1,"int"],
		'only_depth':[0,0,1,"int"],
		'mix_depth':[0,0,1,"int"],
                'only_ssao':[0,0,1,"int"],
                'show_mask_depth':[0,0,1,"int"],
                'scale':[1.0,1.0,100.0,"float"],
                'samples': [6,1,8,"int"],
                'rings': [6,1,8,"int"],
                'aoCap': [1.2,0.0,10.0,"float"],
                'aoMultiplier': [100.0,1.,500.,"float"],
                'depthTolerance': [0.0,0.0,1.0,"float"],
                'aorange':[60.0,1.0,500.0,"float"],
                'negative':[0,0,1,"int"],
		'correction':[6.0,0.0,1000.0,"float"],
                #'force_real_far':[0,0,1,"int"],
                }
        vendor = GL.glGetString(GL.GL_VENDOR)
        #if vendor.find("ATI"):
        self.SSAO_OPTIONS["near"][0] = 2.0
        self.Set(near = 2.0)
        self.SSAO_OPTIONS_ORDER = ['use_fog','only_ssao','only_depth','mix_depth',
				'negative','do_noise','near','scale','rings','samples','aoCap','aoMultiplier',
				'aorange','depthTolerance','use_mask_depth']#'far','correction',


    def copyDepth(self):
        _gllib.glBindTexture (GL_TEXTURE_2D, self.depthtextureName);
        glCopyTexImage2D(GL_TEXTURE_2D, 0,GL_DEPTH_COMPONENT, 0, 0, self._width,
                         self._height, 0);        
        self.copy_depth_ssao = False


    def copyBuffer(self, depthmask=False):
        if depthmask :
            _gllib.glBindTexture (GL_TEXTURE_2D, self.depthmasktextureName);
            glCopyTexImage2D(GL_TEXTURE_2D, 0,GL_DEPTH_COMPONENT, 0, 0, self._width,
                      self._height, 0);
            return 
        if self.copy_depth_ssao : 
            _gllib.glBindTexture (GL_TEXTURE_2D, self.depthtextureName);
            glCopyTexImage2D(GL_TEXTURE_2D, 0,GL_DEPTH_COMPONENT, 0, 0, self._width,
                         self._height, 0);
        
        _gllib.glBindTexture (GL_TEXTURE_2D, self.illumtextureName);
        glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, 0, 0, self._width,
                      self._height, 0);
        
        _gllib.glBindTexture (GL_TEXTURE_2D, self.rendertextureName);
        glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 0, 0, self._width,
                      self._height, 0);        
        
        _gllib.glBindTexture (GL_TEXTURE_2D, 0);
    
    def drawSSAO(self,combine=0):
        if __debug__:
            if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        if DejaVu2.enableSSAO is True and self.ssao:
            #if not self.AR.use_mask: 
            self.copyBuffer()
            from opengltk.extent import _glextlib
            
            import numpy
            _glextlib.glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0 )

            _glextlib.glUseProgram( self.shaderSSAOProgram )

            glActiveTexture(GL_TEXTURE3);
            _gllib.glBindTexture (GL_TEXTURE_2D, self.depthtextureName);
            _glextlib.glUniform1i(self.SSAO_LOCATIONS['DepthTexture'],3)
            
            glActiveTexture(GL_TEXTURE4);
            _gllib.glBindTexture (GL_TEXTURE_2D, self.rendertextureName);
            _glextlib.glUniform1i(self.SSAO_LOCATIONS['RenderedTexture'],4)
            
            glActiveTexture(GL_TEXTURE5);
            _gllib.glBindTexture (GL_TEXTURE_2D, self.illumtextureName);        
            _glextlib.glUniform1i(self.SSAO_LOCATIONS['LuminanceTexture'],5)
            
            glActiveTexture(GL_TEXTURE6);
            _gllib.glBindTexture (GL_TEXTURE_2D, self.randomtextureName);        
            _glextlib.glUniform1i(self.SSAO_LOCATIONS['RandomTexture'],6)            
            
            glActiveTexture(GL_TEXTURE7);
            _gllib.glBindTexture (GL_TEXTURE_2D, self.depthmasktextureName);        
            _glextlib.glUniform1i(self.SSAO_LOCATIONS['DepthMaskTexture'],7)

            _glextlib.glUniform1f(self.SSAO_LOCATIONS['RenderedTextureWidth'],self._width)
            _glextlib.glUniform1f(self.SSAO_LOCATIONS['RenderedTextureHeight'],self._height)
            
            _glextlib.glUniform1f(self.SSAO_LOCATIONS['realfar'],self.far)
            #_glextlib.glUniform1f(self.SSAO_LOCATIONS['realnear'],self.near)

            for k in self.SSAO_OPTIONS:
                if len(self.SSAO_OPTIONS[k]) == 5 :
                    self.SSAO_OPTIONS[k][0] = self.SSAO_OPTIONS[k][-1].get()
                val = self.SSAO_OPTIONS[k][0]
                if k == "correction" :
                    #val = self.SSAO_OPTIONS[k][0]
                    self.SSAO_OPTIONS[k][0] = val = self.far / 10.0
                    if len(self.SSAO_OPTIONS[k]) == 5 :
                        val = self.SSAO_OPTIONS[k][-1].get()
                if k == "far" :
                    #val = self.SSAO_OPTIONS[k][0]
                    val = self.far + 230.
                    if len(self.SSAO_OPTIONS[k]) == 5 :
                        val = self.SSAO_OPTIONS[k][-1].get()
                if k == "fog" :
                    self.SSAO_OPTIONS[k][0] = val = int(self.fog.enabled)
                    if len(self.SSAO_OPTIONS[k]) == 5 :
                        self.SSAO_OPTIONS[k][-1].set(val)
                if self.SSAO_OPTIONS[k][3] == "float":
                    _glextlib.glUniform1f(self.SSAO_LOCATIONS[k],val)
                else : 
                    _glextlib.glUniform1i(self.SSAO_LOCATIONS[k],val)
                
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT |
                GL_STENCIL_BUFFER_BIT) 
            self.drawTexturePolygon()
            self.endSSAO()

    def endSSAO(self):
        if DejaVu2.enableSSAO is True :
            from opengltk.extent import _glextlib
            _glextlib.glUseProgram(0)
            glActiveTexture(GL_TEXTURE0);            
            _gllib.glBindTexture (GL_TEXTURE_2D, 0); 
#            _gllib.glBindTexture (GL_TEXTURE_2D,GL_TEXTURE0)

            
            
    def drawTexturePolygon(self):
        glPushMatrix()
        glLoadIdentity()
        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()
        lViewport = glGetIntegerv(GL_VIEWPORT)
        glOrtho(float(lViewport[0]),
                float(lViewport[0]+lViewport[2]),
                float(lViewport[1]),
                float(lViewport[1]+lViewport[3]),
                -1, 1)        
        glEnable(GL_BLEND)

#        glEnable(GL_TEXTURE_2D)
#        if self.dsT == 0 :
#            _gllib.glBindTexture (GL_TEXTURE_2D, self.depthtextureName);
#        elif self.dsT == 1 :
#            _gllib.glBindTexture (GL_TEXTURE_2D, self.illumtextureName);
#        else :
#            _gllib.glBindTexture (GL_TEXTURE_2D, self.rendertextureName);
        glBegin(GL_POLYGON)
        glTexCoord2i(0, 0)
        glVertex2i(0, 0)
        glTexCoord2i(1, 0)
        glVertex2i(self._width, 0)
        glTexCoord2i(1, 1)
        glVertex2i(self._width, self._height)
        glTexCoord2i(0, 1)
        glVertex2i(0, self._height)
        glEnd()
        glDisable(GL_BLEND)  
        _gllib.glBindTexture (GL_TEXTURE_2D, 0);      
        glPopMatrix()
        glMatrixMode(GL_MODELVIEW)
        glPopMatrix()

##     def DrawImage(self, imarray, mode='RGB', format=GL.GL_UNSIGNED_BYTE,
##                   filter=None, swap=False):
##         """Draw an array of pixel in the camera
## """
##         if imarray is None:
##             return

##         GL.glViewport(0, 0, self._width, self._height)
##         GL.glMatrixMode(GL.GL_PROJECTION)
##         GL.glPushMatrix()
##         GL.glLoadIdentity()
##         GL.glOrtho(0, self._width, 0, self._height, -1.0, 1.0)
##         GL.glMatrixMode(GL.GL_MODELVIEW)
##         GL.glPushMatrix()
##         GL.glLoadIdentity()

##         GL.glRasterPos2i( 0, 0)
##         GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1)

##         if not swap:
##             GL.glDrawBuffer(GL.GL_BACK)

##         if filter:
##             GL.glEnable(GL.GL_CONVOLUTION_2D)
##             GL.glConvolutionFilter2D(GL.GL_CONVOLUTION_2D, GL.GL_LUMINANCE,
##                                      3, 3, GL.GL_LUMINANCE, GL.GL_FLOAT,
##                                      filter)

##         if mode=='RGB':
##             _gllib.glDrawPixels( self._width, self._height,
##                                  GL.GL_RGB, format, imarray)
##         elif mode in ['L','P']:
##             _gllib.glDrawPixels( self._width, self._height,
##                                  GL.GL_LUMINANCE, format, imarray)

##         GL.glMatrixMode(GL.GL_PROJECTION)
##         GL.glPopMatrix()
##         GL.glMatrixMode(GL.GL_MODELVIEW)
##         GL.glPopMatrix()
        
##         GL.glDisable(GL.GL_CONVOLUTION_2D)

##         GL.glFlush()

##         if swap:
##             self.tk.call(self._w, 'swapbuffers')


##     def secondDerivative(self, imarray, mode='RGB', format=GL.GL_UNSIGNED_BYTE,
##                          swap=False):
##         sndDeriv = numpy.array([ [-0.125, -0.125, -0.125,],
##                                    [-0.125,    1.0, -0.125,],
##                                    [-0.125, -0.125, -0.125,] ], 'f')
##         self.DrawImage(imarray, mode, format, sndDeriv, swap)
##         if not swap:
##             buffer=GL.GL_BACK
##         deriv = self.GrabFrontBufferAsArray(lock=False, buffer=buffer)
##         return numpy.fabs(deriv).astype('B')
    

##     def firstDerivative(self, imarray, mode='RGB', format=GL.GL_UNSIGNED_BYTE,
##                         swap=False):
##         if not swap:
##             buffer=GL.GL_BACK
##         else:
##             buffer=GL.GL_FRONT
    
##         fstDeriveV = numpy.array([ [-0.125,  -0.25, -0.125],
##                                      [ 0.0  ,    0.0,  0.0  ],
##                                      [0.125,  0.25, 0.125] ], 'f')

##         self.DrawImage(imarray, mode, format, fstDeriveV, swap)
##         derivV = self.GrabFrontBufferAsArray(lock=False, buffer=buffer)
##         if derivV is None:
##             return None
    
##         fstDeriveV = numpy.array([ [ 0.125,   0.25,  0.125],
##                                      [ 0.0  ,    0.0,  0.0  ],
##                                      [-0.125,  -0.25, -0.125] ], 'f')

##         self.DrawImage(imarray, mode, format, fstDeriveV, swap)
##         derivV += self.GrabFrontBufferAsArray(lock=False, buffer=buffer)
##         derivV = numpy.fabs(derivV*0.5)

##         fstDeriveH = numpy.array([ [-0.125,    0.0, 0.125],
##                                      [-0.25 ,    0.0, 0.25  ],
##                                      [-0.125,    0.0, 0.125] ], 'f')
##         self.DrawImage(imarray, mode, format, fstDeriveH, swap)
##         derivH = self.GrabFrontBufferAsArray(lock=False, buffer=buffer)
        
##         fstDeriveH = numpy.array([ [ 0.125,    0.0, -0.125],
##                                      [ 0.25 ,    0.0, -0.25  ],
##                                      [ 0.125,    0.0, -0.125] ], 'f')
##         self.DrawImage(imarray, mode, format, fstDeriveH, swap)
##         derivH += self.GrabFrontBufferAsArray(lock=False, buffer=buffer)
##         derivH = numpy.fabs(derivH*0.5)

##         return derivH+derivV.astype('B')

    
        
    def Redraw(self, *dummy):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        if self.selectDragRect:
            return
        if not self.initialized:
            return
        if not self.visible:
            return

        # MS commented out when adding QT
        #if not self.master.winfo_ismapped():
        #    return

        # was causing very slow update of light source
        # This function shoudl force redrawing and hence not
        # wait for other things to happen
        #self.update_idletasks()

        # At this point We expect no interruption after this
        # and this OpenGL context should remain the current one until
        # we swap buffers
        self.Activate()
        #glPushAttrib(GL_ALL_ATTRIB_BITS)

        if not self.viewer.tileRender:
            # only setup camera projection if we are not rendering tiles
            # when tiles are rendered the projection is set by the tile
            # renderer tr.beginTile() method

            ## removed by MS. Seems unecessary
            #if self.projectionType==self.ORTHOGRAPHIC:
            #    self.PerspectiveToOrthogonal()
            #else:
            #    self.SetupProjectionMatrix()
            self.SetupProjectionMatrix()

        if self.renderMode == GL_RENDER:
            # This activate has to be done here, else we get a
            # GLXBadContextState on the alpha. If we are in 
            # GL_SELECT mode we did already an activate in DoPick

            # here we need to restore the camrea's GLstate
            if self.viewer and len(self.viewer.cameras) > 1:
                #self.SetupProjectionMatrix()
                self.fog.Set(enabled=self.fog.enabled)

        r,g,b,a = self.backgroundColor
        glClearColor(r, g, b, a )

        self.RedrawAASwitch()

        while 1:
            errCode = glGetError()
            if not errCode: break
            errString = gluErrorString(errCode)
            print 'GL ERROR: %d %s' % (errCode, errString)

       
    def glError(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        while 1:
            errCode = glGetError()
            if not errCode: return
            errString = gluErrorString(errCode)
            print 'GL ERROR: %d %s' % (errCode, errString)
            

    def SetupLights(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """ loops over the lights and sets the direction and/or position of a
        light only if it is enabled and the flag is set to positive. Also, if
        there is more than one camera, it always sets the lights that are on"""
        for l in self.viewer.lights:
            if l.enabled is True:
                if len(self.viewer.cameras) > 1:
                    if l.positional is True:
                        glLightfv(l.num, GL_POSITION, l.position)
                        glLightfv(l.num, GL_SPOT_DIRECTION,
                                         l.spotDirection)
                    else:  # directional
                        glLightfv(l.num, GL_POSITION, l.direction)
                else:
                    if l.positional is True:
                        if l.posFlag:
                            glLightfv(l.num, GL_POSITION, l.position)
                            l.posFlag = 0
                        if l.spotFlag:
                            glLightfv(l.num, GL_SPOT_DIRECTION,
                                         l.spotDirection)
                            l.spotFlag = 0
                    else:  #directional
                        if l.dirFlag:
                            #rot = numpy.reshape(self.rotation, (4,4))
                            #dir = numpy.dot(l.direction, rot).astype('f')
                            #glLightfv(l.num, GL_POSITION, dir)
                            glLightfv(l.num, GL_POSITION, l.direction)
                            l.dirFlag = 0
                            #self.posLog.append(l.direction)
                            #print 'setting light to ',l.direction
##                              l.posFlag = 0
##                              l.spotFlag = 0


try:
    import pymedia.video.vcodec as vcodec
    
    class RecordableCamera(StandardCameraBase):
        """Subclass Camera to define a camera supporting saving mpg video
    """

        def __init__(self, master, screenName, viewer, num, check=1,
                     cnf={}, **kw):
            if __debug__:
             if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

            StandardCameraBase.__init__( self, *(master, screenName, viewer,
                                                 num, check), **kw )
            # this attribute can be 'recording', 'autoPaused', 'paused' or 'stopped'
            self.videoRecordingStatus = 'stopped'
            self.videoOutputFile = None
            self.pendingAutoPauseID = None
            self.pendingRecordFrameID = None
            self.encoder = None
            self.autoPauseDelay = 1 # auto pause after 1 second
            self.pauseLength = 15 # if recording pauses automatically
                                  # add self.pauseLength frames will be added
                                  # when recording resumes
            self.videoRecorder = None
            self.nbRecordedframes = 0
            self.vidimage = None

        def setVideoOutputFile(self, filename='out.mpg'):
            if __debug__:
             if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
            # open the file
            self.videoOutputFile = open( filename, 'wb' )
            assert self.videoOutputFile, 'ERROR, failed to open file: '+filename


        def setVideoParameters(self, width=None, height=None, pauseLength=0,
                               autoPauseDelay=1, outCodec = 'mpeg1video'):
            if __debug__:
             if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
            # make sure image is dimensions are even
            if width is None:
                width = self._width
            w = self.videoFrameWidth = width - width%2

            if height is None:
                height = self._height
            h = self.videoFrameHeight = height - height%2
            if h != self._height or w != self._width:
                self.Set(width=w, height=h)

            ## FIXME here we should lock the size of the Camera
            ##

            params= { 'type': 0, 'gop_size': 12, 'frame_rate_base': 125,
                     'max_b_frames': 0, 'width': w, 'height': h,
                     'frame_rate': 3000,'deinterlace': 0, #'bitrate': 9800000,
                     'id': vcodec.getCodecID(outCodec)
                     }
##             params= { 'type': 0, 'gop_size': 12, 'frame_rate_base': 125,
##                      'max_b_frames': 0, 'width': w, 'height': h,
##                      'frame_rate': 2997,'deinterlace': 0, 'bitrate': 2700000,
##                      'id': vcodec.getCodecID( 'mpeg1video' )
##                      }
            if outCodec== 'mpeg1video':
                params['bitrate']= 2700000
            else:
                params['bitrate']= 9800000
            self.encoder = vcodec.Encoder( params )
            self.pauseLength = pauseLength
            self.autoPauseDelay = autoPauseDelay


        def start(self):
            if __debug__:
             if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
            # toggle to recording mode
            if self.videoOutputFile is None or self.encoder is None:
                return
            self.videoRecordingStatus = 'recording'
            self.nbRecordedframes = 0 
            self.viewer.master.after(1, self.recordFrame)


        def recordFrame(self, force=False, imageFile=None):
            if __debug__:
             if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

            if not self.viewer.hasRedrawn and not force:
                return

            root = self.viewer.master

            if self.videoRecordingStatus=='paused':
                self.pendingRecordFrameID = root.after(1, self.recordFrame)
                return
        
##             if self.videoRecordingStatus=='autopaused':
##                 # add frames for pause
##                 #print 'adding %d pause frames --------------', self.pauseLength
##                 imstr = self.lastFrame.tostring()
##                 for i in range(self.pauseLength):
##                     bmpFrame = vcodec.VFrame(
##                         vcodec.formats.PIX_FMT_RGB24, self.lastFrame.size,
##                         (imstr, None, None))
##                     yuvFrame = bmpFrame.convert(vcodec.formats.PIX_FMT_YUV420P)
##                     self.videoOutputFile.write(
##                         self.encoder.encode(yuvFrame).data)
##                 self.videoRecordingStatus = 'recording'

            if self.videoRecordingStatus=='recording' or force:
                if self.pendingAutoPauseID:
                    print 'auto pause id', self.pendingAutoPauseID
                    root.after_cancel(self.pendingAutoPauseID)
                    self.pendingAutoPauseID = None
#                if DejaVu2.enableSSAO is True and self.ssao :#and self.vidimage is not None:
#                    image = self.GrabFrontBuffer(lock=False)
#                else :
                image = self.GrabFrontBuffer(lock=False)
                    # FIXME this resizing can be avoided if camera size is locked
                image = image.resize((self.videoFrameWidth,
                                      self.videoFrameHeight))
                self.lastFrame = image
                bmpFrame = vcodec.VFrame(
                    vcodec.formats.PIX_FMT_RGB24, image.size,
                    (self.lastFrame.tostring(), None, None))
                #print "bmp rate_base", bmpFrame.rate_base, "rate:", bmpFrame.rate, "bitrate", bmpFrame.bitrate, "size:", bmpFrame.size 
                yuvFrame = bmpFrame.convert(vcodec.formats.PIX_FMT_YUV420P)
                #print "yuv rate_base", yuvFrame.rate_base, "rate:", yuvFrame.rate, "bitrate", yuvFrame.bitrate, "size:", yuvFrame.size 
                self.videoOutputFile.write(self.encoder.encode(yuvFrame).data)
                #self.pendingAutoPauseID = root.after(self.autoPauseDelay*1000,
                #                                     self.autoPause)
                if imageFile != None:
                    image.save(imageFile)
                self.viewer.hasRedrawn = False

                self.pendingRecordFrameID = root.after(1, self.recordFrame)
                self.nbRecordedframes += 1

        def Redraw(self):
            if __debug__:
             if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
            StandardCamera.Redraw(self)
            self.viewer.hasRedrawn = True
            self.recordFrame()
            
            
        def autoPause(self):
            if __debug__:
             if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
            #print 'autoPause =========================='
            self.videoRecordingStatus = 'autoPaused'
            if self.pendingRecordFrameID:
                root = self.viewer.master
                root.after_cancel(self.pendingRecordFrameID)    
            self.pendingAutoPauseID = None


        def pause(self):
            if __debug__:
             if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
            #print 'pause =========================='
            self.videoRecordingStatus = 'paused'
            if self.pendingRecordFrameID:
                root = self.viewer.master
                root.after_cancel(self.pendingRecordFrameID)    
            self.pendingAutoPauseID = None


        def stop(self):
            if __debug__:
             if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
            self.videoRecordingStatus = 'stopped'
            self.videoOutputFile.close()
            root = self.viewer.master
            if self.pendingAutoPauseID:
                root.after_cancel(self.pendingAutoPauseID)
            if self.pendingRecordFrameID:
                root.after_cancel(self.pendingRecordFrameID)

    Camera = RecordableCamera

except ImportError:
    Camera = StandardCameraBase
