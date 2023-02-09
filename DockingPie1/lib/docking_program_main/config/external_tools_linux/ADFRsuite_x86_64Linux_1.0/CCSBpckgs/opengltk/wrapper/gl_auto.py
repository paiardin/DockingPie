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

from opengltk.extent import _gllib
from opengltk.util import glGetXXDim, GLbyte, GLdouble, GLfloat, GLint, GLdouble, GLfloat, GLushort, GLshort, GLboolean, GLuint, GLint, GLubyte
import numpy

def glGetBooleanv( pname):
    """

return numpy array ( GLboolean)
"""
    result = numpy.zeros( glGetXXDim[ pname], GLboolean)
    _gllib.glGetBooleanv( pname, result)
    return result

def glGetBoolean( pname):
    """
    
return numpy array( GLboolean) or GLboolean if singleton
"""
    result = glGetBooleanv( pname)
    if 1 == len( result):
        return result[ 0]
    else:
        return result

def glGetIntegerv( pname):
    """

return numpy array( GLint)
"""
    result = numpy.zeros( glGetXXDim[ pname], GLint)
    _gllib.glGetIntegerv( pname, result)
    return result

def glGetInteger( pname):
    """
    
return numpy array( GLint) or GLint if singleton
"""
    result = glGetIntegerv( pname)
    if 1 == len( result):
        return result[ 0]
    else:
        return result

def glGetFloatv( pname):
    """

return numpy array( GLfloat)
"""
    result = numpy.zeros( glGetXXDim[ pname], GLfloat)
    _gllib.glGetFloatv( pname, result)
    return result

def glGetFloat( pname):
    """
    
return numpy array( GLfloat) or GLfloat if singleton
"""
    result = glGetFloatv( pname)
    if 1 == len( result):
        return result[ 0]
    else:
        return result

def glGetDoublev( pname):
    """

return numpy array( GLdouble)
"""
    result = numpy.zeros( glGetXXDim[ pname], GLdouble)
    _gllib.glGetDoublev( pname, result)
    return result

def glGetDouble( pname):
    """
    
return numpy array ( GLdouble) or GLdouble if singleton
"""
    result = glGetDoublev( pname)
    if 1 == len( result):
        return result[ 0]
    else:
        return result

def glGetLightfv(light, pname):
    """light - GLenum

return numpy array( GLfloat)
"""
    result = numpy.zeros( glGetXXDim[ pname], GLfloat)
    _gllib.glGetLightfv(light, pname, result)
    return result

def glGetLightf(light, pname):
    """
    light - GLenum
return numpy array( GLfloat) or GLfloat if singleton
"""
    result = glGetLightfv(light, pname)
    if 1 == len( result):
        return result[ 0]
    else:
        return result

def glGetLightiv(light, pname):
    """light - GLenum

return numpy array ( GLint)
"""
    result = numpy.zeros( glGetXXDim[ pname], GLint)
    _gllib.glGetLightiv(light, pname, result)
    return result

def glGetLighti(light, pname):
    """
    light - GLenum
return numpy array ( GLint) or GLint if singleton
"""
    result = glGetLightiv(light, pname)
    if 1 == len( result):
        return result[ 0]
    else:
        return result

def glGetMaterialfv(face, pname):
    """face - GLenum

return numpy array( GLfloat)
"""
    result = numpy.zeros( glGetXXDim[ pname], GLfloat)
    _gllib.glGetMaterialfv(face, pname, result)
    return result

def glGetMaterialf(face, pname):
    """
    face - GLenum
return numpy array( GLfloat) or GLfloat if singleton
"""
    result = glGetMaterialfv(face, pname)
    if 1 == len( result):
        return result[ 0]
    else:
        return result

def glGetMaterialiv(face, pname):
    """face - GLenum

return numpy array( GLint)
"""
    result = numpy.zeros( glGetXXDim[ pname], GLint)
    _gllib.glGetMaterialiv(face, pname, result)
    return result

def glGetMateriali(face, pname):
    """
    face - GLenum
return numpy array( GLint) or GLint if singleton
"""
    result = glGetMaterialiv(face, pname)
    if 1 == len( result):
        return result[ 0]
    else:
        return result

def glGetTexEnvfv(target, pname):
    """target - GLenum

return numpy array( GLfloat)
"""
    result = numpy.zeros( glGetXXDim[ pname], GLfloat)
    _gllib.glGetTexEnvfv(target, pname, result)
    return result

def glGetTexEnvf(target, pname):
    """
    target - GLenum
return numpy array( GLfloat) or GLfloat if singleton
"""
    result = glGetTexEnvfv(target, pname)
    if 1 == len( result):
        return result[ 0]
    else:
        return result

def glGetTexEnviv(target, pname):
    """target - GLenum

return numpy array( GLint)
"""
    result = numpy.zeros( glGetXXDim[ pname], GLint)
    _gllib.glGetTexEnviv(target, pname, result)
    return result

def glGetTexEnvi(target, pname):
    """
    target - GLenum
return numpy array( GLint) or GLint if singleton
"""
    result = glGetTexEnviv(target, pname)
    if 1 == len( result):
        return result[ 0]
    else:
        return result

def glGetTexGendv(coord, pname):
    """coord - GLenum

return numpy array( GLdouble)
"""
    result = numpy.zeros( glGetXXDim[ pname], GLdouble)
    _gllib.glGetTexGendv(coord, pname, result)
    return result

def glGetTexGend(coord, pname):
    """
    coord - GLenum
return numpy array( GLdouble) or GLdouble if singleton
"""
    result = glGetTexGendv(coord, pname)
    if 1 == len( result):
        return result[ 0]
    else:
        return result

def glGetTexGenfv(coord, pname):
    """coord - GLenum

return numpy array( GLfloat)
"""
    result = numpy.zeros( glGetXXDim[ pname], GLfloat)
    _gllib.glGetTexGenfv(coord, pname, result)
    return result

def glGetTexGenf(coord, pname):
    """
    coord - GLenum
return numpy array( GLfloat) or GLfloat if singleton
"""
    result = glGetTexGenfv(coord, pname)
    if 1 == len( result):
        return result[ 0]
    else:
        return result

def glGetTexGeniv(coord, pname):
    """coord - GLenum

return numpy array( GLint)
"""
    result = numpy.zeros( glGetXXDim[ pname], GLint)
    _gllib.glGetTexGeniv(coord, pname, result)
    return result

def glGetTexGeni(coord, pname):
    """
    coord - GLenum
return numpy array( GLint) or GLint if singleton
"""
    result = glGetTexGeniv(coord, pname)
    if 1 == len( result):
        return result[ 0]
    else:
        return result

def glGetTexLevelParameterfv(target, level, pname):
    """target - GLenum
    level - GLint

return numpy array( GLfloat)
"""
    result = numpy.zeros( glGetXXDim[ pname], GLfloat)
    _gllib.glGetTexLevelParameterfv(target, level, pname, result)
    return result

def glGetTexLevelParameterf(target, level, pname):
    """
    target - GLenum
    level - GLint
return numpy array( GLfloat) or GLfloat if singleton
"""
    result = glGetTexLevelParameterfv(target, level, pname)
    if 1 == len( result):
        return result[ 0]
    else:
        return result

def glGetTexLevelParameteriv(target, level, pname):
    """target - GLenum
    level - GLint

return numpy array( GLint)
"""
    result = numpy.zeros( glGetXXDim[ pname], GLint)
    _gllib.glGetTexLevelParameteriv(target, level, pname, result)
    return result

def glGetTexLevelParameteri(target, level, pname):
    """
    target - GLenum
    level - GLint
return numpy array( GLint) or GLint if singleton
"""
    result = glGetTexLevelParameteriv(target, level, pname)
    if 1 == len( result):
        return result[ 0]
    else:
        return result

def glGetTexParameterfv(target, pname):
    """target - GLenum

return numpy array( GLfloat)
"""
    result = numpy.zeros( glGetXXDim[ pname], GLfloat)
    _gllib.glGetTexParameterfv(target, pname, result)
    return result

def glGetTexParameterf(target, pname):
    """
    target - GLenum
return numpy array( GLfloat) or GLfloat if singleton
"""
    result = glGetTexParameterfv(target, pname)
    if 1 == len( result):
        return result[ 0]
    else:
        return result

def glGetTexParameteriv(target, pname):
    """target - GLenum

return numpy array( GLint)
"""
    result = numpy.zeros( glGetXXDim[ pname], GLint)
    _gllib.glGetTexParameteriv(target, pname, result)
    return result

def glGetTexParameteri(target, pname):
    """
    target - GLenum
return numpy array( GLint) or GLint if singleton
"""
    result = glGetTexParameteriv(target, pname)
    if 1 == len( result):
        return result[ 0]
    else:
        return result

def glColor3bv( v):
    """
    v - seq( GLbyte, 3)
    """
    if 3 != len( v):
        raise TypeError( len( v), "3-array expected")
    _gllib.glColor3bv( v)

def glColor3dv( v):
    """
    v - seq( GLdouble, 3)
    """
    if 3 != len( v):
        raise TypeError( len( v), "3-array expected")
    _gllib.glColor3dv( v)

def glColor3fv( v):
    """
    v - seq( GLfloat, 3)
    """
    if 3 != len( v):
        raise TypeError( len( v), "3-array expected")
    _gllib.glColor3fv( v)

def glColor3iv( v):
    """
    v - seq( GLint, 3)
    """
    if 3 != len( v):
        raise TypeError( len( v), "3-array expected")
    _gllib.glColor3iv( v)

def glColor3sv( v):
    """
    v - seq( GLshort, 3)
    """
    if 3 != len( v):
        raise TypeError( len( v), "3-array expected")
    _gllib.glColor3sv( v)

def glColor3ubv( v):
    """
    v - seq( GLubyte, 3)
    """
    if 3 != len( v):
        raise TypeError( len( v), "3-array expected")
    _gllib.glColor3ubv( v)

def glColor3uiv( v):
    """
    v - seq( GLuint, 3)
    """
    if 3 != len( v):
        raise TypeError( len( v), "3-array expected")
    _gllib.glColor3uiv( v)

def glColor3usv( v):
    """
    v - seq( GLushort, 3)
    """
    if 3 != len( v):
        raise TypeError( len( v), "3-array expected")
    _gllib.glColor3usv( v)

def glColor4bv( v):
    """
    v - seq( GLbyte, 4)
    """
    if 4 != len( v):
        raise TypeError( len( v), "4-array expected")
    _gllib.glColor4bv( v)

def glColor4dv( v):
    """
    v - seq( GLdouble, 4)
    """
    if 4 != len( v):
        raise TypeError( len( v), "4-array expected")
    _gllib.glColor4dv( v)

def glColor4fv( v):
    """
    v - seq( GLfloat, 4)
    """
    if 4 != len( v):
        raise TypeError( len( v), "4-array expected")
    _gllib.glColor4fv( v)

def glColor4iv( v):
    """
    v - seq( GLint, 4)
    """
    if 4 != len( v):
        raise TypeError( len( v), "4-array expected")
    _gllib.glColor4iv( v)

def glColor4sv( v):
    """
    v - seq( GLshort, 4)
    """
    if 4 != len( v):
        raise TypeError( len( v), "4-array expected")
    _gllib.glColor4sv( v)

def glColor4ubv( v):
    """
    v - seq( GLubyte, 4)
    """
    if 4 != len( v):
        raise TypeError( len( v), "4-array expected")
    _gllib.glColor4ubv( v)

def glColor4uiv( v):
    """
    v - seq( GLuint, 4)
    """
    if 4 != len( v):
        raise TypeError( len( v), "4-array expected")
    _gllib.glColor4uiv( v)

def glColor4usv( v):
    """
    v - seq( GLushort, 4)
    """
    if 4 != len( v):
        raise TypeError( len( v), "4-array expected")
    _gllib.glColor4usv( v)

def glEvalCoord1dv( v):
    """
    v - seq( GLdouble, 1)
    """
    if 1 != len( v):
        raise TypeError( len( v), "1-array expected")
    _gllib.glEvalCoord1dv( v)

def glEvalCoord1fv( v):
    """
    v - seq( GLfloat, 1)
    """
    if 1 != len( v):
        raise TypeError( len( v), "1-array expected")
    _gllib.glEvalCoord1fv( v)

def glEvalCoord2dv( v):
    """
    v - seq( GLdouble, 2)
    """
    if 2 != len( v):
        raise TypeError( len( v), "2-array expected")
    _gllib.glEvalCoord2dv( v)

def glEvalCoord2fv( v):
    """
    v - seq( GLfloat, 2)
    """
    if 2 != len( v):
        raise TypeError( len( v), "2-array expected")
    _gllib.glEvalCoord2fv( v)

def glRasterPos2dv( v):
    """
    v - seq( GLdouble, 2)
    """
    if 2 != len( v):
        raise TypeError( len( v), "2-array expected")
    _gllib.glRasterPos2dv( v)

def glRasterPos2fv( v):
    """
    v - seq( GLfloat, 2)
    """
    if 2 != len( v):
        raise TypeError( len( v), "2-array expected")
    _gllib.glRasterPos2fv( v)

def glRasterPos2iv( v):
    """
    v - seq( GLint, 2)
    """
    if 2 != len( v):
        raise TypeError( len( v), "2-array expected")
    _gllib.glRasterPos2iv( v)

def glRasterPos2sv( v):
    """
    v - seq( GLshort, 2)
    """
    if 2 != len( v):
        raise TypeError( len( v), "2-array expected")
    _gllib.glRasterPos2sv( v)

def glRasterPos3dv( v):
    """
    v - seq( GLdouble, 3)
    """
    if 3 != len( v):
        raise TypeError( len( v), "3-array expected")
    _gllib.glRasterPos3dv( v)

def glRasterPos3fv( v):
    """
    v - seq( GLfloat, 3)
    """
    if 3 != len( v):
        raise TypeError( len( v), "3-array expected")
    _gllib.glRasterPos3fv( v)

def glRasterPos3iv( v):
    """
    v - seq( GLint, 3)
    """
    if 3 != len( v):
        raise TypeError( len( v), "3-array expected")
    _gllib.glRasterPos3iv( v)

def glRasterPos3sv( v):
    """
    v - seq( GLshort, 3)
    """
    if 3 != len( v):
        raise TypeError( len( v), "3-array expected")
    _gllib.glRasterPos3sv( v)

def glRasterPos4dv( v):
    """
    v - seq( GLdouble, 4)
    """
    if 4 != len( v):
        raise TypeError( len( v), "4-array expected")
    _gllib.glRasterPos4dv( v)

def glRasterPos4fv( v):
    """
    v - seq( GLfloat, 4)
    """
    if 4 != len( v):
        raise TypeError( len( v), "4-array expected")
    _gllib.glRasterPos4fv( v)

def glRasterPos4iv( v):
    """
    v - seq( GLint, 4)
    """
    if 4 != len( v):
        raise TypeError( len( v), "4-array expected")
    _gllib.glRasterPos4iv( v)

def glRasterPos4sv( v):
    """
    v - seq( GLshort, 4)
    """
    if 4 != len( v):
        raise TypeError( len( v), "4-array expected")
    _gllib.glRasterPos4sv( v)

def glTexCoord1dv( v):
    """
    v - seq( GLdouble, 1)
    """
    if 1 != len( v):
        raise TypeError( len( v), "1-array expected")
    _gllib.glTexCoord1dv( v)

def glTexCoord1fv( v):
    """
    v - seq( GLfloat, 1)
    """
    if 1 != len( v):
        raise TypeError( len( v), "1-array expected")
    _gllib.glTexCoord1fv( v)

def glTexCoord1iv( v):
    """
    v - seq( GLint, 1)
    """
    if 1 != len( v):
        raise TypeError( len( v), "1-array expected")
    _gllib.glTexCoord1iv( v)

def glTexCoord1sv( v):
    """
    v - seq( GLshort, 1)
    """
    if 1 != len( v):
        raise TypeError( len( v), "1-array expected")
    _gllib.glTexCoord1sv( v)

def glTexCoord2dv( v):
    """
    v - seq( GLdouble, 2)
    """
    if 2 != len( v):
        raise TypeError( len( v), "2-array expected")
    _gllib.glTexCoord2dv( v)

def glTexCoord2fv( v):
    """
    v - seq( GLfloat, 2)
    """
    if 2 != len( v):
        raise TypeError( len( v), "2-array expected")
    _gllib.glTexCoord2fv( v)

def glTexCoord2iv( v):
    """
    v - seq( GLint, 2)
    """
    if 2 != len( v):
        raise TypeError( len( v), "2-array expected")
    _gllib.glTexCoord2iv( v)

def glTexCoord2sv( v):
    """
    v - seq( GLshort, 2)
    """
    if 2 != len( v):
        raise TypeError( len( v), "2-array expected")
    _gllib.glTexCoord2sv( v)

def glTexCoord3dv( v):
    """
    v - seq( GLdouble, 3)
    """
    if 3 != len( v):
        raise TypeError( len( v), "3-array expected")
    _gllib.glTexCoord3dv( v)

def glTexCoord3fv( v):
    """
    v - seq( GLfloat, 3)
    """
    if 3 != len( v):
        raise TypeError( len( v), "3-array expected")
    _gllib.glTexCoord3fv( v)

def glTexCoord3iv( v):
    """
    v - seq( GLint, 3)
    """
    if 3 != len( v):
        raise TypeError( len( v), "3-array expected")
    _gllib.glTexCoord3iv( v)

def glTexCoord3sv( v):
    """
    v - seq( GLshort, 3)
    """
    if 3 != len( v):
        raise TypeError( len( v), "3-array expected")
    _gllib.glTexCoord3sv( v)

def glTexCoord4dv( v):
    """
    v - seq( GLdouble, 4)
    """
    if 4 != len( v):
        raise TypeError( len( v), "4-array expected")
    _gllib.glTexCoord4dv( v)

def glTexCoord4fv( v):
    """
    v - seq( GLfloat, 4)
    """
    if 4 != len( v):
        raise TypeError( len( v), "4-array expected")
    _gllib.glTexCoord4fv( v)

def glTexCoord4iv( v):
    """
    v - seq( GLint, 4)
    """
    if 4 != len( v):
        raise TypeError( len( v), "4-array expected")
    _gllib.glTexCoord4iv( v)

def glTexCoord4sv( v):
    """
    v - seq( GLshort, 4)
    """
    if 4 != len( v):
        raise TypeError( len( v), "4-array expected")
    _gllib.glTexCoord4sv( v)

def glVertex2dv( v):
    """
    v - seq( GLdouble, 2)
    """
    if 2 != len( v):
        raise TypeError( len( v), "2-array expected")
    _gllib.glVertex2dv( v)

def glVertex2fv( v):
    """
    v - seq( GLfloat, 2)
    """
    if 2 != len( v):
        raise TypeError( len( v), "2-array expected")
    _gllib.glVertex2fv( v)

def glVertex2iv( v):
    """
    v - seq( GLint, 2)
    """
    if 2 != len( v):
        raise TypeError( len( v), "2-array expected")
    _gllib.glVertex2iv( v)

def glVertex2sv( v):
    """
    v - seq( GLshort, 2)
    """
    if 2 != len( v):
        raise TypeError( len( v), "2-array expected")
    _gllib.glVertex2sv( v)

def glVertex3dv( v):
    """
    v - seq( GLdouble, 3)
    """
    if 3 != len( v):
        raise TypeError( len( v), "3-array expected")
    _gllib.glVertex3dv( v)

def glVertex3fv( v):
    """
    v - seq( GLfloat, 3)
    """
    if 3 != len( v):
        raise TypeError( len( v), "3-array expected")
    _gllib.glVertex3fv( v)

def glVertex3iv( v):
    """
    v - seq( GLint, 3)
    """
    if 3 != len( v):
        raise TypeError( len( v), "3-array expected")
    _gllib.glVertex3iv( v)

def glVertex3sv( v):
    """
    v - seq( GLshort, 3)
    """
    if 3 != len( v):
        raise TypeError( len( v), "3-array expected")
    _gllib.glVertex3sv( v)

def glVertex4dv( v):
    """
    v - seq( GLdouble, 4)
    """
    if 4 != len( v):
        raise TypeError( len( v), "4-array expected")
    _gllib.glVertex4dv( v)

def glVertex4fv( v):
    """
    v - seq( GLfloat, 4)
    """
    if 4 != len( v):
        raise TypeError( len( v), "4-array expected")
    _gllib.glVertex4fv( v)

def glVertex4iv( v):
    """
    v - seq( GLint, 4)
    """
    if 4 != len( v):
        raise TypeError( len( v), "4-array expected")
    _gllib.glVertex4iv( v)

def glVertex4sv( v):
    """
    v - seq( GLshort, 4)
    """
    if 4 != len( v):
        raise TypeError( len( v), "4-array expected")
    _gllib.glVertex4sv( v)

def glDeleteTextures( seq):
    """
    vseq - sequence( GLuint)
    """
    return _gllib.glDeleteTextures( len( seq), seq)

def glGenTextures( n):
    """
return - sequence( GLuint, n)
    """
    result = numpy.zeros( n, GLuint)
    _gllib.glGenTextures( n, result)
    return result

def glIndexdv( v):
    """
    v - seq( GLdouble, 1)
"""
    if 1 != len( v):
        raise TypeError( len( v), "1-array expected")
    _gllib.glIndexdv( v)

def glIndexfv( v):
    """
    v - seq( GLfloat, 1)
"""
    if 1 != len( v):
        raise TypeError( len( v), "1-array expected")
    _gllib.glIndexfv( v)

def glIndexiv( v):
    """
    v - seq( GLint, 1)
"""
    if 1 != len( v):
        raise TypeError( len( v), "1-array expected")
    _gllib.glIndexiv( v)

def glIndexsv( v):
    """
    v - seq( GLshort, 1)
"""
    if 1 != len( v):
        raise TypeError( len( v), "1-array expected")
    _gllib.glIndexsv( v)

def glIndexubv( v):
    """
    v - seq( GLubyte, 1)
"""
    if 1 != len( v):
        raise TypeError( len( v), "1-array expected")
    _gllib.glIndexubv( v)

def glRectdv( v1, v2):
    """
    v1, v2 - seq( GLdouble, 2)
"""
    if 2 != len( v1):
        raise TypeError( len( v1), "2-array expected for v1")
    if 2 != len( v2):
        raise TypeError( len( v2), "2-array expected for v2")
    _gllib.glRectdv( v1, v2)

def glRectfv( v1, v2):
    """
    v1, v2 - seq( GLfloat, 2)
"""
    if 2 != len( v1):
        raise TypeError( len( v1), "2-array expected for v1")
    if 2 != len( v2):
        raise TypeError( len( v2), "2-array expected for v2")
    _gllib.glRectfv( v1, v2)

def glRectiv( v1, v2):
    """
    v1, v2 - seq( GLint, 2)
"""
    if 2 != len( v1):
        raise TypeError( len( v1), "2-array expected for v1")
    if 2 != len( v2):
        raise TypeError( len( v2), "2-array expected for v2")
    _gllib.glRectiv( v1, v2)

def glRectsv( v1, v2):
    """
    v1, v2 - seq( GLshort, 2)
"""
    if 2 != len( v1):
        raise TypeError( len( v1), "2-array expected for v1")
    if 2 != len( v2):
        raise TypeError( len( v2), "2-array expected for v2")
    _gllib.glRectsv( v1, v2)

def glLightfv(light, pname, parms):
    """light - GLenum
    parms - sequence
"""
    if len( parms) != glGetXXDim[ pname]:
        raise TypeError( len( parms), glGetXXDim[ pname],
                         "wrong size of parms")
    _gllib.glLightfv(light, pname, parms)

def glLightiv(light, pname, parms):
    """light - GLenum
    parms - sequence
"""
    if len( parms) != glGetXXDim[ pname]:
        raise TypeError( len( parms), glGetXXDim[ pname],
                         "wrong size of parms")
    _gllib.glLightiv(light, pname, parms)

def glMaterialfv(face, pname, parms):
    """face - GLenum
    parms - sequence
"""
    if len( parms) != glGetXXDim[ pname]:
        raise TypeError( len( parms), glGetXXDim[ pname],
                         "wrong size of parms")
    _gllib.glMaterialfv(face, pname, parms)

def glMaterialiv(face, pname, parms):
    """face - GLenum
    parms - sequence
"""
    if len( parms) != glGetXXDim[ pname]:
        raise TypeError( len( parms), glGetXXDim[ pname],
                         "wrong size of parms")
    _gllib.glMaterialiv(face, pname, parms)

def glLightModelfv( pname, parms):
    """
    parms - sequence
"""
    if len( parms) != glGetXXDim[ pname]:
        raise TypeError( len( parms), glGetXXDim[ pname],
                         "wrong size of parms")
    _gllib.glLightModelfv( pname, parms)

def glLightModeliv( pname, parms):
    """
    parms - sequence
"""
    if len( parms) != glGetXXDim[ pname]:
        raise TypeError( len( parms), glGetXXDim[ pname],
                         "wrong size of parms")
    _gllib.glLightModeliv( pname, parms)

def glFogfv( pname, parms):
    """
    parms - sequence
"""
    if len( parms) != glGetXXDim[ pname]:
        raise TypeError( len( parms), glGetXXDim[ pname],
                         "wrong size of parms")
    _gllib.glFogfv( pname, parms)

def glFogiv( pname, parms):
    """
    parms - sequence
"""
    if len( parms) != glGetXXDim[ pname]:
        raise TypeError( len( parms), glGetXXDim[ pname],
                         "wrong size of parms")
    _gllib.glFogiv( pname, parms)

def glTexEnvfv(target, pname, parms):
    """target - GLenum
    parms - sequence
"""
    if len( parms) != glGetXXDim[ pname]:
        raise TypeError( len( parms), glGetXXDim[ pname],
                         "wrong size of parms")
    _gllib.glTexEnvfv(target, pname, parms)

def glTexEnviv(target, pname, parms):
    """target - GLenum
    parms - sequence
"""
    if len( parms) != glGetXXDim[ pname]:
        raise TypeError( len( parms), glGetXXDim[ pname],
                         "wrong size of parms")
    _gllib.glTexEnviv(target, pname, parms)

def glTexGenfv(coord, pname, parms):
    """coord - GLenum
    parms - sequence
"""
    if len( parms) != glGetXXDim[ pname]:
        raise TypeError( len( parms), glGetXXDim[ pname],
                         "wrong size of parms")
    _gllib.glTexGenfv(coord, pname, parms)

def glTexGeniv(coord, pname, parms):
    """coord - GLenum
    parms - sequence
"""
    if len( parms) != glGetXXDim[ pname]:
        raise TypeError( len( parms), glGetXXDim[ pname],
                         "wrong size of parms")
    _gllib.glTexGeniv(coord, pname, parms)

def glTexParameterfv(target, pname, parms):
    """target - GLenum
    parms - sequence
"""
    if len( parms) != glGetXXDim[ pname]:
        raise TypeError( len( parms), glGetXXDim[ pname],
                         "wrong size of parms")
    _gllib.glTexParameterfv(target, pname, parms)

def glTexParameteriv(target, pname, parms):
    """target - GLenum
    parms - sequence
"""
    if len( parms) != glGetXXDim[ pname]:
        raise TypeError( len( parms), glGetXXDim[ pname],
                         "wrong size of parms")
    _gllib.glTexParameteriv(target, pname, parms)

def glPixelMapfv(map, values):
    """
    map - GLenum
    values - seq( GLfloat)
    """
    _gllib.glPixelMapfv(map, values)

def glPixelMapuiv(map, values):
    """
    map - GLenum
    values - seq( GLuint)
    """
    _gllib.glPixelMapuiv(map, values)

def glPixelMapusv(map, values):
    """
    map - GLenum
    values - seq( GLushort)
    """
    _gllib.glPixelMapusv(map, values)

__all__ = ["glGetBoolean",
    "glGetBooleanv",
    "glGetInteger",
    "glGetIntegerv",
    "glGetFloat",
    "glGetFloatv",
    "glGetDouble",
    "glGetDoublev",
    "glGetLightf",
    "glGetLightfv",
    "glGetLighti",
    "glGetLightiv",
    "glGetMaterialf",
    "glGetMaterialfv",
    "glGetMateriali",
    "glGetMaterialiv",
    "glGetTexEnvf",
    "glGetTexEnvfv",
    "glGetTexEnvi",
    "glGetTexEnviv",
    "glGetTexGend",
    "glGetTexGendv",
    "glGetTexGenf",
    "glGetTexGenfv",
    "glGetTexGeni",
    "glGetTexGeniv",
    "glGetTexLevelParameterf",
    "glGetTexLevelParameterfv",
    "glGetTexLevelParameteri",
    "glGetTexLevelParameteriv",
    "glGetTexParameterf",
    "glGetTexParameterfv",
    "glGetTexParameteri",
    "glGetTexParameteriv",
    "glColor3bv",
    "glColor3dv",
    "glColor3fv",
    "glColor3iv",
    "glColor3sv",
    "glColor3ubv",
    "glColor3uiv",
    "glColor3usv",
    "glColor4bv",
    "glColor4dv",
    "glColor4fv",
    "glColor4iv",
    "glColor4sv",
    "glColor4ubv",
    "glColor4uiv",
    "glColor4usv",
    "glEvalCoord1dv",
    "glEvalCoord1fv",
    "glEvalCoord2dv",
    "glEvalCoord2fv",
    "glRasterPos2dv",
    "glRasterPos2fv",
    "glRasterPos2iv",
    "glRasterPos2sv",
    "glRasterPos3dv",
    "glRasterPos3fv",
    "glRasterPos3iv",
    "glRasterPos3sv",
    "glRasterPos4dv",
    "glRasterPos4fv",
    "glRasterPos4iv",
    "glRasterPos4sv",
    "glTexCoord1dv",
    "glTexCoord1fv",
    "glTexCoord1iv",
    "glTexCoord1sv",
    "glTexCoord2dv",
    "glTexCoord2fv",
    "glTexCoord2iv",
    "glTexCoord2sv",
    "glTexCoord3dv",
    "glTexCoord3fv",
    "glTexCoord3iv",
    "glTexCoord3sv",
    "glTexCoord4dv",
    "glTexCoord4fv",
    "glTexCoord4iv",
    "glTexCoord4sv",
    "glVertex2dv",
    "glVertex2fv",
    "glVertex2iv",
    "glVertex2sv",
    "glVertex3dv",
    "glVertex3fv",
    "glVertex3iv",
    "glVertex3sv",
    "glVertex4dv",
    "glVertex4fv",
    "glVertex4iv",
    "glVertex4sv",
    "glDeleteTextures",
    "glGenTextures",
    "glIndexdv",
    "glIndexfv",
    "glIndexiv",
    "glIndexsv",
    "glIndexubv",
    "glRectdv",
    "glRectfv",
    "glRectiv",
    "glRectsv",
    "glLightfv",
    "glLightiv",
    "glMaterialfv",
    "glMaterialiv",
    "glLightModelfv",
    "glLightModeliv",
    "glFogfv",
    "glFogiv",
    "glTexEnvfv",
    "glTexEnviv",
    "glTexGenfv",
    "glTexGeniv",
    "glTexParameterfv",
    "glTexParameteriv",
    "glPixelMapfv",
    "glPixelMapuiv",
    "glPixelMapusv"]
