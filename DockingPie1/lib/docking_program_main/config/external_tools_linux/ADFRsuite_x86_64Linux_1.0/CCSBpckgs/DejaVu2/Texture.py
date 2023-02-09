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
#############################################################################

#
# $Header: /mnt/raid/services/cvs/DejaVu2/Texture.py,v 1.1.1.1.4.2 2017/11/08 23:27:21 annao Exp $
#
# $Id: Texture.py,v 1.1.1.1.4.2 2017/11/08 23:27:21 annao Exp $
#

import sys, os
from PIL import Image
import numpy
import warnings

from opengltk.OpenGL import GL
import viewerConst
import math
from colorTool import OneColor
from Transformable import Transformable
from Geom import Geom
from viewerFns import checkKeywords, getkw

class Texture(Transformable):
    """Base class for textures"""

    keywords = [
        'enable',
        'wrapS',
        'wrapT',
        'magFilter',
        'minFilter',
        'genModS',
        'genModT',
        'genPlaneS',
        'genPlaneT',
        'planeS',
        'planeT',
        'level',
        'auto',
        'envMode',
        'envColor',
        'border',
        'format',
        'image',
        
        ]
    
    def __init__(self, check=1, viewer=None, **kw):
        
        if __debug__:
            if check:
                apply( checkKeywords, ('texture',self.keywords), kw)

        Transformable.__init__(self, viewer)

        self.objects = [] # will list the geometries to which this texture
                          # is applied
        self.name = 'NoName'
        self.enabled = viewerConst.NO
        self.dim = GL.GL_TEXTURE_2D
        self.wrap = [GL.GL_REPEAT, ]*3
        self.magFilter = GL.GL_NEAREST
        self.minFilter = GL.GL_NEAREST
        self.border = 0
        self.genMod = [ GL.GL_OBJECT_LINEAR, GL.GL_OBJECT_LINEAR]
        self.genPlane = [ GL.GL_OBJECT_PLANE, GL.GL_OBJECT_PLANE]
        self.plane = [ [1.0, 0.0, 0.0, 0.0],   # x=0
                       [0.0, 1.0, 0.0, 0.0] ]  # y=0
        self.image = None
        self.width = 0
        self.height = 0
        self.level = 0
        self.format = GL.GL_RGB
        self.auto = 1  # automatic computation of texture indices
        self.envMode = GL.GL_MODULATE
        self.envColor = OneColor( (1.0, 1.0, 0.0, 0.0) )
        
        # as we resize the texture to be of the power of two,
        # this allow to recompute texture coordinate after
        self.resizeRatio = (1, 1)

        apply( self.Set, (0,), kw)


    def Set(self, check=1, **kw):
	"""Set various parameters for texture objects:
- image=im  set an image as the texture to be used. im can be either a numpy
array of bytes with shape (n x m x d) where d can be either 3 or 4.
Im can also be a PIL Image object. If thes image is not of type 'RGB' or
'RGBA' it will be converted to 'RGB' if possible, else ... error !

"The set command will also pad the image to make sure that width and height
are powers of 2." (has been tested for Image but not for numpy array yet - guillaume)
"""
        if __debug__:
            if check:
                apply( checkKeywords, (self.name,self.keywords), kw)

	val = getkw(kw, 'enable')
	if not val is None:
	    assert val in (viewerConst.YES, viewerConst.NO), "enable can only be YES or NO"
	    self.enabled = val

	val = getkw(kw, 'wrapS')
	if not val is None:
	    assert val in (GL.GL_REPEAT, GL.GL_CLAMP)
	    self.wrap[0] = val

	val = getkw(kw, 'wrapT')
	if not val is None:
	    assert val in (GL.GL_REPEAT, GL.GL_CLAMP)
	    self.wrap[1] = val

	val = getkw(kw, 'magFilter')
	if not val is None:
	    assert val in (GL.GL_NEAREST, GL.GL_LINEAR)
	    self.magFilter = val

	val = getkw(kw, 'minFilter')
	if not val is None:
	    assert val in (GL.GL_NEAREST, GL.GL_LINEAR)
	    self.minFilter = val

	val = getkw(kw, 'genModS')
	if not val is None:
	    assert val in (GL.GL_OBJECT_LINEAR, GL.GL_EYE_LINEAR, GL.GL_SPHERE_MAP)
	    self.genMod[0] = val

	val = getkw(kw, 'genModT')
	if not val is None:
	    assert val in (GL.GL_OBJECT_LINEAR, GL.GL_EYE_LINEAR, GL.GL_SPHERE_MAP)
	    self.genMod[1] = val

	val = getkw(kw, 'genPlaneS')
	if not val is None:
	    assert val in (GL.GL_OBJECT_PLANE, GL.GL_EYE_PLANE)
	    self.genPlane[0] = val

	val = getkw(kw, 'genPlaneT')
	if not val is None:
	    assert val in (GL.GL_OBJECT_PLANE, GL.GL_EYE_PLANE)
	    self.genPlane[1] = val

	val = getkw(kw, 'planeS')
	if not val is None:
	    assert len(val)==4 and type(val[0])==type(0.0), "Plane has to be 4Float vector"
	    self.plane[0] = val

	val = getkw(kw, 'planeT')
	if not val is None:
	    assert len(val)==4 and type(val[0])==type(0.0), "Plane has to be 4Float vector"
	    self.plane[1] = val

	val = getkw(kw, 'level')
	if not val is None:
	    assert type(val)==type(0) and val >=0
	    self.level = val

	val = getkw(kw, 'auto')
	if not val is None:
	    assert val in (viewerConst.YES,viewerConst.NO), "auto can only be YES or NO"
	    self.auto = val

	val = getkw(kw, 'envMode')
	if not val is None:
	    assert val in (GL.GL_MODULATE, GL.GL_REPLACE, GL.GL_DECAL, GL.GL_BLEND ), "envMode can only be GL_MODULATE, GL_DECAL, or GL_BLEND"
	    self.envMode = val
	val = getkw(kw, 'envColor')
	if val:
	    col = OneColor(val)
	    if col: self.envColor = col

	b = getkw(kw, 'border')
	f = getkw(kw, 'format')
	im = getkw(kw, 'image')
	if im is not None:

            if isinstance(im, Image.Image):
                lImInstaceImage = True
            else:
                lImInstaceImage = False

            if lImInstaceImage is True:
                width = im.size[0]
                height = im.size[1]
            else:
                height = im.shape[0]
                width = im.shape[1]

            # find smallest power of 2 larger than image size
            dim1 = 1
            dim2 = 1
            while dim1 < width:
                dim1 = dim1 << 1
            while dim2 < height:
                dim2 = dim2 << 1
            self.resizeRatio = ( width / float(dim1),
                                 height / float(dim2) )

            if os.name != 'nt': #sys.platform != 'win32':
                lMaxTextureSize = GL.glGetInteger(GL.GL_MAX_TEXTURE_SIZE)
                #print "width", width, height, dim1, dim2, lMaxTextureSize
                if (dim1 > lMaxTextureSize) or (dim2 > lMaxTextureSize):
                    warnings.warn('texture map too big for this implementation of opengl %d'%lMaxTextureSize)

            if lImInstaceImage is True:
                if im.mode !='RGB' and im.mode !='RGBA':
                    im = im.convert('RGB')
                im = im.transpose(Image.FLIP_TOP_BOTTOM)
                imstr = im.tostring()
                imarr = numpy.fromstring( imstr, 'B')
                if im.mode=='RGB':
                    imarr.shape = (height, width, 3)
                elif im.mode=='RGBA':
                    imarr.shape = (height, width, 4)
                im = imarr
                
            if (dim1 != width) or (dim2 != height):
                if len(im.shape) == 3:
                    newArray = numpy.zeros( (dim2, dim1, len(im[0][0]) ) )
                else:
                    newArray = numpy.zeros( (dim2, dim1 ) )
                for i in range(height):
                    for j in range(width):    
                        newArray[i][j] = im[i][j]
                im = newArray.astype('B')

	    if b:
		assert type(b)==type(0)
	    else: b = 0
	    self.border = b

	    if f:
		assert f in (GL.GL_RGB, GL.GL_RGBA), "format can only be GL_RGB or GL_RGBA"

	    assert type(im).__name__ == 'ndarray'
	    assert im.dtype.char == 'B'
	    if im.shape[-1] == 3:
		if f and f != GL.GL_RGB: raise ValueError("bad image format")
	        self.format = GL.GL_RGB
	    elif im.shape[-1] == 4:
		if f and f != GL.GL_RGBA: raise ValueError("bad image format")
	        self.format = GL.GL_RGBA

                for o in self.objects:
                    o.transparent=1
                    o.inheritMaterial=0
                    if self.viewer:
                        self.viewer.objectsNeedingRedo[o] = None
                    
	    l = len(im.shape)
	    if l==2:
		w=im.shape[0] - 2*b
		q, r  = divmod(math.log(w)/math.log(2), 1)
		if r != 0.0:
		    raise ValueError("Image width must be 2**m +2b")
		self.dim = GL.GL_TEXTURE_1D
		self.image = im
		self.width = im.shape[0]

	    elif l==3:
		w=im.shape[0] - 2*b
		q, r  = divmod(math.log(w)/math.log(2), 1)
		if r != 0.0:
		    raise ValueError("Image width must be 2**m +2b")
		h=im.shape[1] -2*b
		q, r  = divmod(math.log(h)/math.log(2), 1)
		if r != 0.0:
		    raise ValueError("Image height must be 2**m +2b")
		self.dim = GL.GL_TEXTURE_2D
		self.image = im
		self.width = im.shape[1]
		self.height = im.shape[0]

	    else:
		raise ValueError("Bad shape for image")

        if self.viewer:
            self.viewer.deleteOpenglList()
            self.viewer.Redraw()


    def getTextureCoordinatesAfterResizeRatio(self, textureCoordinates):
        """ to be called with textureCoordinates=geom.vertexSet.texCoords
"""
        if self.resizeRatio == (1., 1.):
            return textureCoordinates.copy()
        lenTextureCoordinates = len(textureCoordinates)
        assert lenTextureCoordinates > 0
        assert len(textureCoordinates[0]) == 2
        newTextureCoordinates = []
        for i in range(lenTextureCoordinates):
            newTextureCoordinates.append (
                           ( self.resizeRatio[0] * textureCoordinates[i][0],
                             self.resizeRatio[1] * textureCoordinates[i][1] ) )

        return newTextureCoordinates


    def MakeMatrix(self):
        """Build texture transformation matrix"""

        t = self
        GL.glMatrixMode(GL.GL_TEXTURE)
        GL.glLoadIdentity()
        GL.glTranslatef(float(t.translation[0]),float(t.translation[1]),float(t.translation[2]))
        GL.glMultMatrixf(t.rotation)
        GL.glScalef(float(t.scale[0]),float(t.scale[1]),float(t.scale[2]))
        GL.glMatrixMode(GL.GL_MODELVIEW)


    def dump(self):
	print  '\nenabled: ', self.enabled, \
	       '\ndim: ', self.dim, \
	       '\nwrap: ', self.wrap, \
	       '\nmagFilter: ', self.magFilter, \
	       '\nminFilter: ', self.minFilter, \
	       '\nborder: ', self.border, \
	       '\ngenMod: ', self.genMod, \
	       '\ngenPlane: ', self.genPlane, \
	       '\nplane: ', self.plane, \
	       '\nimage: ', self.image.shape, \
	       '\nwidth: ', self.width, \
	       '\nheight: ', self.height, \
	       '\nlevel: ', self.level, \
	       '\nformat: ', self.format, \
	       '\nauto: ', self.auto, \
	       '\nenvMode: ', self.envMode, \
	       '\nenvColor: ', self.envColor


    def Setup(self):
	"""Setup texture mapping"""
	
	t = self
	if not t.enabled:
	    GL.glDisable(t.dim)
	    return

	t.MakeMatrix()

	GL.glTexEnvf(GL.GL_TEXTURE_ENV,
		     GL.GL_TEXTURE_ENV_MODE, t.envMode)
	if t.envMode==GL.GL_BLEND:
	    GL.glTexEnvf(GL.GL_TEXTURE_ENV,
			 GL.GL_TEXTURE_ENV_COLOR, t.envColor)

	GL.glTexParameterf(t.dim, GL.GL_TEXTURE_WRAP_S, t.wrap[0])
	GL.glTexParameterf(t.dim, GL.GL_TEXTURE_MAG_FILTER, 
			   t.magFilter)
	GL.glTexParameterf(t.dim, GL.GL_TEXTURE_MIN_FILTER, 
			   t.minFilter)
	if t.auto:
 	    GL.glTexGeni(GL.GL_S, GL.GL_TEXTURE_GEN_MODE,
 			 t.genMod[0])
 	    GL.glTexGenfv(GL.GL_S, t.genPlane[0], t.plane[0])
	    GL.glEnable(GL.GL_TEXTURE_GEN_S)
        else:
	    GL.glDisable(GL.GL_TEXTURE_GEN_S)

	if t.dim==GL.GL_TEXTURE_1D:
            from opengltk.extent import _gllib
	    _gllib.glTexImage1D(t.dim, t.level, t.format,
			    t.width, 0, t.format,
			    GL.GL_UNSIGNED_BYTE, t.image)

	elif t.dim==GL.GL_TEXTURE_2D:
	    GL.glTexParameterf(t.dim, GL.GL_TEXTURE_WRAP_T, t.wrap[1])
            from opengltk.extent import _gllib
            # directly call the C function to not have to transform the
            # t.image into a bufarray
            
            #print "t.width, t.height", t.width, t.height
            
	    _gllib.glTexImage2D(t.dim, t.level, t.format,
			    t.width, t.height, t.border!=0, t.format,
			    GL.GL_UNSIGNED_BYTE, t.image)
##             import bufarray
##             bimage = bufarray.Bufarray(im, bufarray.ubyteCtype)
## 	    GL.glTexImage2D(t.dim, t.level,GL.GL_UNSIGNED_BYTE,
##                             t.width, t.height, t.border!=0, t.format,
##                             bimage)
            
## 	    GL.glTexImage2D(t.dim, t.level, t.format,
## 			    t.width, t.height, t.border!=0, t.format,
## 			    GL.GL_UNSIGNED_BYTE, t.image)
            
	    if t.auto:
 		GL.glTexGeni(GL.GL_T, GL.GL_TEXTURE_GEN_MODE,
 			     t.genMod[1])
 		GL.glTexGenfv(GL.GL_T, t.genPlane[1], t.plane[1])
		GL.glEnable(GL.GL_TEXTURE_GEN_T)
            else:
                GL.glDisable(GL.GL_TEXTURE_GEN_T)

 	GL.glEnable(t.dim)


if __name__ == '__main__':

    t = Texture()
    from PIL import Image
    im = Image.open('lena.jpg')
    t.Set(enable=1, image=im)

    # if g is a DejaVu2 geometry
    #g.Set(texture=t)
