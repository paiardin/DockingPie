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
## Copyright (c) MGL TSRI 2016
##
################################################################################

#
# copyright_notice
#

"""gl wrappers
"""

# the function that are checked in the C wrapper are commented out

import numpy
from gl_auto import *
from gl_deprec import *

from opengltk.extent import _gllib
from opengltk import util
from opengltk.util import glGetXXDim ,gltypmap, revtypmap, \
     GLboolean,\
     GLbyte,\
     GLubyte,\
     GLshort,\
     GLushort,\
     GLint,\
     GLuint,\
     GLfloat,\
     GLdouble

## if we do this here by hand we get better doc string than in gl_auto
##
##  def glFogfv( pname, params):
##      """glFogfv( pname, params)
##      pname           params
##      GL_FOG_MODE     GL_LINEAR, GL_EXP, GL_EXP2, GL_FOG_FUNC_SGIS
##      GL_FOG_DENSITY  float
##      GL_FOG_START    float
##      GL_FOG_END      float
##      GL_FOG_INDEX    float
##      GL_FOG_COLOR    4-floats
##      GL_FOG_OFFSET_VALUE_SGIX 4-floats
##      """
##      params = readArray( params, GLfloat)
##      if glGetXXDim[ pname] != len( params):
##          raise TypeError( glGetXXDim[ pname], len( params), 'wrong len( param)')
##      _gllib.glFogfv( pname, params)


##  def glFogfi( pname, params):
##      """glFogfi( pname, params)
##      pname           params
##      GL_FOG_MODE     GL_LINEAR, GL_EXP, GL_EXP2, GL_FOG_FUNC_SGIS
##      GL_FOG_DENSITY  int
##      GL_FOG_START    int
##      GL_FOG_END      int
##      GL_FOG_INDEX    int
##      GL_FOG_COLOR    4-ints
##      GL_FOG_OFFSET_VALUE_SGIX 4-ints
##      """
##      params = readArray( params, GLint)
##      if glGetXXDim[ pname] != len( params):
##          raise TypeError( glGetXXDim[ pname], len( params), 'wrong len( param)')
##      _gllib.glFogfi( pname, params)


##  def glMaterialfv(face, pname, params ):
##      """glMaterialfv(face, pname, params )
##      face: GL_FRONT, GL_BACK, GL_FRONT_AND_BACK
##      pname                     params - sequence
##      GL_SHININESS              float
##      GL_AMBIENT                4-float
##      GL_DIFFUSE                4-float
##      GL_SPECULAR               4-float
##      GL_EMISSION               4-float
##      GL_AMBIENT_AND_DIFFUSE    4-float
##      GL_COLOR_INDEXES          3-float
##      """
##      params = readArray( params, GLfloat)
##      if glGetXXDim[ pname] != len( params):
##          raise TypeError( glGetXXDim[ pname], len( params), 'wrong len( param)')
##      _gllib.glMaterialfv(face, pname, params)

    
##  def glMaterialfi(face, pname, params ):
##      """glMaterialfi(face, pname, params )
##      face: GL_FRONT, GL_BACK, GL_FRONT_AND_BACK
##      pname                     params - sequence
##      GL_SHININESS              int
##      GL_AMBIENT                4-int
##      GL_DIFFUSE                4-int
##      GL_SPECULAR               4-int
##      GL_EMISSION               4-int
##      GL_AMBIENT_AND_DIFFUSE    4-int
##      GL_COLOR_INDEXES          3-int
##      """
##      params = readArray( params, GLint)
##      if glGetXXDim[ pname] != len( params):
##          raise TypeError( glGetXXDim[ pname], len( params), 'wrong len( param)')
##      _gllib.glMaterialfv(face, pname, params)


##################################################################




def glAreTexturesResident( textures):
    """wrapper over _gllib.glAreTexturesResident
    Return None if all textures are resident, or a numpy Array.
    """

    ncell = len( textures)
    
    residences = numpy.zeros( ncell, GLboolean)
    
    return not(_gllib.glAreTexturesResident( ncell, textures, residences) ) \
           and residences or None


def glCallLists( lists):
    """
    lists - Num. array
    """
    _gllib.glCallLists( len( lists), gltypmap[ lists.dtype.char], lists)
    
 
def glClipPlane( plane, equation):
    if 4 != len( equation):
        raise TypeError( len( equation), 'len( equation) must be 4')
    _gllib.glClipPlane( plane, equation)


def glColorPointer( size, pointer, stride=0, gltype=None):
    """
    size - 3 or 4
    pointer - Num. array or Python list
    stride - int
    gltype - GL<type>
    """
    if not gltype:
        if  isinstance(pointer, numpy.ndarray):
            gltype = gltypmap[pointer.dtype.char]
        else:
            raise ValueError( (gltype,), "GL<type> is not known")
    else:
        if  not isinstance(pointer, numpy.ndarray):
            pointer = numpy.array(pointer, revtypmap[gltype])
    _gllib.glColorPointer( size, gltype ,stride, pointer)



def glDrawElements( mode, indices):
    """
    mode - GLenum
    indices - Num. array
    """
    _gllib.glDrawElements( mode, len( indices),
                          gltypmap[ indices.dtype.char], indices)




def glDrawPixels( width, format, pixels):
    """
    width - GLsizei
    format - GLenum
    pixels - numpy array
    """
    npix = len( pixels)
    if npix % width:
        raise TypeError( npix, width, npix % width,
                         'len( pixels) must be a multiple of "width"')
    _gllib.glDrawPixels( width, npix/width, format,
                        gltypmap[ pixels.dtype.char], pixels)



def glEdgeFlagv( flag):
    """
    flag - numpy array
    """
    _gllib.glEdgeFlagv( flag)



def glEdgeFlagPointer( pointer, stride=0):
    """
    pointer - numpy array
    stride - int
    """
    _gllib.glEdgeFlagPointer( stride, pointer)
    

def glFeedbackBuffer( size, type):
    """wrapper over glFeedbackBuffer
    Return a size-GLfloat buffer
    """
    result = numpy.zeros( size, GLfloat)
    _gllib.glFeedbackBuffer( size, type, result)
    return result


def glFogiv( pname, params):
    """
    pname - GLenum
    params - sequence
    """
    if glGetXXDim[ pname] != len( params):
        raise TypeError( glGetXXDim[ pname], len( params), 'wrong len( param)')
    _gllib.glFogiv( pname, params)

def glGetClipPlane( plane):
    """
    plane - GLenum
    return - numpy.array( 4, GLdouble)
"""
    result = numpy.zeros( 4, GLdouble)
    _gllib.glGetClipPlane( plane, result)
    return result



def glGetPolygonStipple( mask):
    """To Wrap: what is the expact output format ??"""
    result = numpy.zeros( 128,  GLubyte)
    _gllib.glGetPolygonStipple( result)
    return result



def glIndexPointer( pointer, stride=0, ctype=None):
    """
    pointer - numpy array or python list
    stride - int
    ctype - GL<type>
    """
    if ctype is None:
        if  isinstance(pointer, numpy.ndarray):
            _gllib.glIndexPointer( gltypmap[pointer.dtype.char],
                                   stride, pointer)
        else:
            raise ValueError( (ctype,), "GL<type> is not known")
    else:
        if  not isinstance(pointer, numpy.ndarray):
            pointer = numpy.array(pointer, revtypmap[ctype])
        _gllib.glColorPointer(ctype ,stride, pointer)



def glInterleavedArrays( format, pointer, stride=0):
    """
    format - Glenum
    pointer - numpy array 
    stride - int
    """
    _gllib.glInterleavedArrays( format, stride, pointer)
    

def glLightfv( light, pname, params):
    """
    light - Glenum
    pname - GLenum
    params - numpy array/ Python list
    """
    if len( params) != glGetXXDim[ pname]:
        raise TypeError( glGetXXDim[ pname], len( params),
                         'params has the wrong size')
    _gllib.glLightfv( light, pname, params)
    

def glLightiv( light, pname, params):
    """
    light - GLenum
    pname - GLenum
    params - numpy array/ Python list
    """
    if len( params) != glGetXXDim[ pname]:
        raise TypeError( glGetXXDim[ pname], len( params),
                         'params has the wrong size')
    _gllib.glLightiv( light, pname, params)


def glLoadMatrixd( m):
    """
    m - sequence( GLdouble, 16)
    """
    if 16 != len( m):
        raise TypeError( len( m), 'len(m) must be 16')
    _gllib.glLoadMatrixd( m)
    

def glLoadMatrixf( m):
    """
    m - sequence( GLfloat, 16)
    """
    if 16 != len( m):
        raise TypeError( len( m), 'len( m, GLfloat) must be 16')
    _gllib.glLoadMatrixf( m)


def __mapstride2order( stride, lpoints):
    if stride:
        if lpoints % stride:
            raise TypeError( (stride, lpoints), '(stride, point) uncompatible')
        else:
            return lpoints / stride
    else:
        return lpoints / datadim
    

def glMap1d( target, u1, u2, points, stride=0):
    """
    target - GLenum
    u1, u2 - float
    points - sequence( GLdouble)
    stride - int
Note: C 'order' is deducted from (target, len( points), stride)
    """
    datadim = util.mapDim[ target]
    if stride and stride < datadim:
        raise ValueError( (target, datadim, stride),
                          '(target, dim( target), stride) uncompatible')
    _gllib.glMap1d( target, u1, u2, stride,
                   __mapstride2order( stride, len( points)), points)
    

def glMap1f( target, u1, u2, points, stride=0):
    """
    target - GLenum
    u1, u2 - float
    points - sequence( GLfloat)
    stride - int
Note: C 'order' is deducted from (target, len( points), stride)
    """
    datadim = util.mapDim[ target]
    if stride and stride < datadim:
        raise ValueError( (target, datadim, stride),
                          '(target, dim( target), stride) uncompatible')
    _gllib.glMap1f( target, u1, u2, stride,
                   __mapstride2order( stride, len( points)), points)
    

def glMap2d( target, u1, u2, v1, v2, points, ustride=0, vstride=0):
    """
    target - GLenum
    u1, u2, v1, v2 - float
    points - sequence( GLdouble)
    ustride, vstride - int
Note: C '[uv]order' are deducted from (target, len( points), [uv]stride)
    """
    datadim = util.mapDim[ target]
    for stride in (ustride, vstride):
        if stride and stride < datadim:
            raise ValueError( (target, datadim, stride),
                              '(target, dim( target), stride) uncompatible')
    lpoints = len( points)
    _gllib.glMap2d( target,
                   u1, u2, ustride, __mapstride2order( ustride, lpoints),
                   v1, v2, vstride, __mapstride2order( vstride, lpoints),
                   points)
    

def glMap2f( target, u1, u2, v1, v2, points, ustride=0, vstride=0):
    """
    target - GLenum
    u1, u2, v1, v2 - float
    points - sequence( GLdouble)
    ustride, vstride - int
Note: C '[uv]order' are deducted from (target, len( points), [uv]stride)
    """
    datadim = util.mapDim[ target]
    for stride in (ustride, vstride):
        if stride and stride < datadim:
            raise ValueError( (target, datadim, stride),
                              '(target, dim( target), stride) uncompatible')
    lpoints = len( points)
    _gllib.glMap2f( target,
                   u1, u2, ustride, __mapstride2order( ustride, lpoints),
                   v1, v2, vstride, __mapstride2order( vstride, lpoints),
                   points)
    


def glMultMatrixd( m):
    """
    m - sequence( GLdouble, 16)
    """
    if 16 != len( m):
        raise TypeError( len( m), 'len( m, GLdouble) must be 16')
    _gllib.glMultMatrixd( m)
    

def glMultMatrixf( m):
    """
    m - sequence( GLfloat, 16)
    """

    if 16 != len( m):
        raise TypeError( len( m), 'len( m, GLfloat) must be 16')
    _gllib.glMultMatrixf( m)


    
def glNormalPointer( pointer, stride=0, ctype=None):
    """
    pointer - numpy array/ list
    stride - int
    """
    if ctype is None:
        if  isinstance(pointer, numpy.ndarray):
            _gllib.glNormalPointer( gltypmap[pointer.dtype.char],
                                   stride, pointer)
        else:
            raise ValueError( (ctype,), "GL<type> is not known")
    else:
        if  not isinstance(pointer, numpy.ndarray):
            pointer = numpy.array(pointer, revtypmap[ctype])
        _gllib.glNormalPointer(ctype ,stride, pointer)

def glPolygonStipple( mask):
    """mask - sequence (128, GLubyte) """
    
    if len(mask) != 128:
        raise TypeError( len( mask), 128,
                         'mask must have the size 128 GLubytes')
    
    if isinstance(mask, numpy.ndarray):
        if mask.dtype.char != 'B': #numpy.uint8
            mask = mask.astype('B')
    else:
        mask = numpy.array(mask, 'B')
        
    return _gllib.glPolygonStipple( mask)


def glPrioritizeTextures( textures, priorities):
    """
    textures - seq( GLuint)
    priorities - seq( GLclampf)
    """
    n = len( textures)
    if n != len( priorities):
        raise TypeError( len( textures), len( priorities),
                         'texture and priorities must have the same size')
    _gllib.glPrioritizeTextures( n, textures, priorities)


def glReadPixels( x, y, width, format, pixels):
    """
    x, y - GLint
    width - GLsizei
    format - GLenum
    pixels - numpy array
    """
    npix = len( pixels)
    if npix % width:
        raise TypeError( npix, width, npix % width,
                         'len( pixels) must be a multiple of "width"')
    _gllib.glReadPixels( x, y, width, npix/width, format,
                        gltypmap[ pixels.dtype.char], pixels)


def glSelectBuffer( size):
    """
    size - int
    """
    result = numpy.zeros( size, GLuint)
    _gllib.glSelectBuffer( size, result)
    return result


def glTexCoordPointer( size, pointer, stride=0, ctype=None):
    """
    size - 1, 2, 3 or 4
    pointer - numpy array or Python list
    stride - int
    ctype - GL<type>
    """
    if ctype is None:
        if  isinstance(pointer, numpy.ndarray):
            _gllib.glTexCoordPointer( size, gltypmap[pointer.dtype.char],
                                   stride, pointer)
        else:
            raise ValueError( (ctype,), "GL<type> is not known")
    else:
        if  not isinstance(pointer, numpy.ndarray):
            pointer = numpy.array(pointer, revtypmap[ctype])
            _gllib.glTexCoordPointer(size, ctype ,stride, pointer)
    

def glTexImage1D( target, level, internalformat, border, format, pixels):
    """
    pixels - numpy array
    """
    width = len( pixels)
    _gllib.glTexImage1D( target, level, internalformat, width, border,
                        format, gltypmap[ pixels.dtype.char], pixels)


def glTexImage2D( target, level, internalformat, width, height, border, format,
                  pixels):
    """
    pixels - numpy array
    """
    pixlen = len( pixels)
    if pixlen % width:
        raise TypeError( (width, len( pixels)),
                         'len( pixels) must be a multiple of width')
    length = pixlen / width
    _gllib.glTexImage2D( target, level, internalformat, width,  height, border,
                        format, gltypmap[ pixels.dtype.char], pixels)


def glTexSubImage1D( target, level, xoffset, width, format, pixels):
    """
    pixels - numpy array
    """
    width = len( pixels)
    _gllib.glTexSubImage1D( target, level, xoffset, width,
                           format, gltypmap[ pixels.ctype], pixels)


def glTexSubImage2D( target, level, xoffset, yoffset, width, height,format,
                     pixels):
    """
    pixels - numpy array
    """
    pixlen = len( pixels)
    if pixlen % width:
        raise TypeError( (width, len( pixels)),
                         'len( pixels) must be a multiple of width')
    length = pixlen / width
    _gllib.glTexSubImage2D( target, level, xoffset, yoffset, width,height,
                           format, gltypmap[ pixels.dtype.char], pixels)


def glVertexPointer( size, pointer, stride=0, ctype=None):
    """
    size - 2, 3 or 4
    pointer - numpy array or Python list
    stride - int
    ctype - Ctype
    """
    if ctype is None:
        if  isinstance(pointer, numpy.ndarray):
            _gllib.glVertexPointer( size, gltypmap[pointer.dtype.char],
                                   stride, pointer)
        else:
            raise ValueError( (ctype,), "GL<type> is not known")
    else:
        if  not isinstance(pointer, numpy.ndarray):
            pointer = numpy.array(pointer, revtypmap[ctype])
            _gllib.glVertexPointer(size, ctype ,stride, pointer)
            

__all__ = [
    'glAreTexturesResident',
    'glCallLists',
    'glClipPlane',
    'glColorPointer',
    'glDeleteTextures',
    'glDrawElements',
    'glDrawPixels',
    'glEdgeFlagv',
    'glEdgeFlagPointer',
    'glFeedbackBuffer',
#    'glFogfv',
#    'glFogiv',
    'glGetClipPlane',
    'glGetPolygonStipple',
#    'glLightfv',
#    'glLightiv',
    'glLoadMatrixd',
    'glLoadMatrixf',
    'glMultMatrixd',
    'glMultMatrixf',
    'glNormalPointer',
    'glPolygonStipple',
    'glPrioritizeTextures',
    'glReadPixels',
    'glSelectBuffer',
    'glTexCoordPointer',
    'glTexImage1D',
    'glTexImage2D',
    'glTexSubImage1D',
    'glTexSubImage2D',
    'glVertexPointer',
    ]

import gl_auto
__all__ += gl_auto.__all__
import gl_deprec
__all__ += gl_deprec.__all__
