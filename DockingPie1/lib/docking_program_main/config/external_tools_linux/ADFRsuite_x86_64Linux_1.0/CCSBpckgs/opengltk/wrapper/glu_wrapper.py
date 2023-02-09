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

"""glu wrappers
"""
import numpy
from opengltk import util
from opengltk.extent import _gllib, _glulib
from opengltk.util import gltypmap, revtypmap ,\
     GLdouble, GLfloat, GLint, readModProjView
from opengltk.ccallback import swigptrcallback
import gl_wrapper


def gluBuild1DMipMaps( target, component, format, data):
    """
    target - GLenum
    component - GLint
    format - GLenum
    data - numpy array
    """
    _glulib.gluBuild1DMipMaps( target, component, len( data), format,
                              gltypmap[ data.dtype.char], data)


def gluBuild2DMipMaps( target, component, width, format, data):
    """
    target - GLenum
    component - GLint
    width - GLsizei
    format - GLenum
    data - numpy array
    """
    lendat = len( data)
    if lendat % width:
        raise TypeError( (width, len( data)),
                         'len( data) nota multiple of "width"')
    _glulib.gluBuild2DMipMaps( target, component, width, lendat/width, format,
                              gltypmap[ data.dtype.char], data)


def gluGetNurbsProperty( nurb, property):
    """
    nurb - GLUnurbs*
    property - GLenum
return - GLfloat
"""
    data = numpy.zeros( 1, 'f')
    _glulib.gluGetNurbsProperty( nurb, property, data)
    return data[ 0]


def gluGetTessProperty( tess, property):
    """
    tess - GLUtesselator*
    property - GLenum
return - GLdouble
"""
    data = numpy.zeros( 1, 'd')
    _glulib.gluGetTessProperty( nurb, property, data)
    return data[ 0]


def gluLoadSamplingMatrices( nurb, model, perspective, view):
    """
    nurb - GLUnurb
    model - sequence of 16 floats
    perspective - sequence of 16 floats
    view - sequence of 4 floats
    """
    model, perspective, view = readModProjView( model, perspective, view)
    _glulib.gluLoadSamplingMatrices( nurb, model, perspective, view)
    

def gluNurbsCallback( element, which, callback):
    return _glulib.gluNurbsCallback( element, which,
                                    swigptrcallback( util,
                                                     ('void', ['GLenum']),
                                                     callback))

def gluNurbsCurve( nurb, type, knots, order, control, stride):
    """
    nurb - GLUnurbs
    type - GLenum
    knots - seq( GLfloat)
    order - GLint
    control - seq( GLfloat)
    stride - GLint
    """

    lknot = len( knots)
    assert stride
    if lknot % stride:
        raise TypeError( (len( knots), stride),
                         "(len( knots), stride) not compatible")
    nknot = lknot / stride
    if len( control) != nknot - order:
        raise TypeError( (len( knots), order, stride, len( control)),
                         '(len( knots), order, stride, len( control))'
                         ' not compatible')
    
    _glulib.gluNurbsCurve( nurb, nknot, knots, stride,
                          control, order, type)


def gluNurbsSurface( nurb, type, sKnots, tKnots, sOrder, tOrder, control,
                     sStride, tStride):
    """
    nurb - GLUnurbs
    type - GLenum
    sKnots, tKnots - seq( GLfloat)
    sOrder, tOrder - GLint
    control - seq( GLfloat)
    sStride, tStride - GLint
    """
    assert sStride
    assert tStride

    lsknot = len( sKnots)
    if lsknot % sStride:
        raise TypeError( (len( sKnots), sStride),
                         "(len( sKnots), sStride) not compatible")
    nsknot = lsknot / sStride

    
    ltknot = len( tKnots)
    if ltknot % sStride:
        raise TypeError( (len( tKnots), sStride),
                         "(len( tKnots), sStride) not compatible")
    ntknot = ltknot / sStride

    
    if len( control) != (nsknot - sOrder)*(ntknot - tOrder):
        raise TypeError( (len( control), (nsknot, sOrder), (ntknot, tOrder)),
                         '(len( control), (nsknot, sOrder), (ntknot, tOrder))'
                         ' not compatible')
    
    _glulib.gluNurbsSurface( nurb, nsknot, sKnots, ntknot, tKnots,
                            sStride, tStride, control, sOrder, tOrder, type)


def gluPickMatrix( x, y, delX, delY, view):
    """
    x, y, delX, delY - GLdouble
    """
    if len( view) != 4:
        raise TypeError( view, len( view), '4 expected')
    if not isinstance(view , numpy.ndarray):
        view = numpy.array(view, GLint)
    _glulib.gluPickMatrix( x, y, delX, delY, view)


def gluProject( obj, model, proj, view):
    """
    obj - seq( 3, GLdouble)
    model - seq( gl.MODELVIEW_MATRIX, GLfloat)
    perspective - seq( gl.PROJECTION_MATRIX, GLfloat)
    view - seq( gl.VIEWPORT, GLfloat)
    return - numpy array (shape: (3,), type: GLdouble)
    """
    win = numpy.zeros( 3, GLdouble)
    gluProjectv( obj, model, proj, view, win[0:1], win[ 1:2], win[ 2:3])
    return win

def gluProjectv( obj, model, proj, view, win):
    """
    obj - seq( 3, GLdouble)
    model - seq( gl.MODELVIEW_MATRIX, GLfloat)
    perspective - seq( gl.PROJECTION_MATRIX, GLfloat)
    view - seq( gl.VIEWPORT, GLfloat)
    win - numpy array (shape: (3,), type: GLdouble): gets the new coordinates
    """

    if 3 != len( obj):
        raise TypeError( len( obj), "obj not 3-array")
    model, proj, view = readModProjView( model, proj, view)
    _glulib.gluProject( obj[ 0], obj[ 1], obj[ 2], model, proj, view,
                     win[0:1], win[ 1:2], win[ 2:3])


def gluPwlCurve( nurb, type, data, stride):
    """
    nurb - GLUnurbs*
    type - GLenum
    data - seq( GLfloat)
    stride - GLint
    """
    assert stride
    ldata = len( data)
    if ldata % stride:
        raise TypeError( (len( data), stride),
                         "len( data) not a multiple of stride")
    _glulib.gluPwlCurve( nurb, len( data) / stride, data, stride, type)


def gluQuadricCallback( element, which, callback):
    return _glulib.gluQuadricCallback( element, which,
                                      swigptrcallback( util,
                                                       ('void', ['GLenum']),
                                                       callback))


def gluScaleImage( format, wIn, dataIn, wOut, hOut, typeOut):
    """
    format - GLenum
    wIn - GLsizei
    dataIn - numpy array / list
    wOut, hOut - GLsizei
    typeOut - GLenum
    return - numpy array(shape: wOut*hOut, type; typeOut)
    """
    lin = len( dataIn)
    if lin % wIn:
        raise TypeError( (len( dataIn), wIn),
                         "len( dataIn) not multiple of wIn")

    gltout = typeOut
    typeOut = revtypmap[ typeout]
    result = numpy.zeros( wOut * hOut, typeOut)
    _glulib.gluScaleImage( format,
                          wIn, lin / wIn, gltypmap[ dataIn.dtype.char], dataIn,
                          wOut, hOut, gltout, result)
    return result

        
def gluTessCallback( element, which, callback):
    return _glulib.gluTessCallback( element, which,
                                   swigptrcallback( util,
                                                    ('void', ['GLenum']),
                                                    callback))


def gluTessNormaldv( tess, normale):
    """
    tess - GLUtesselator*
    normale - seq( 3, GLdouble)
    """
    _glulib.gluTessNormal( tess, normale[ 0], normale[ 1], normale[ 2])


def gluTessNormalfv( tess, normale):
    """
    tess - GLUtesselator*
    normale - seq( 3, GLfloat)
    """
    _glulib.gluTessNormal( tess, normale[ 0], normale[ 1], normale[ 2])


def gluTessVertex( tess, location, data):
    """
    tess - GLUtesselator*
    location - seq( 3, GLdouble)
    data - GLvoid*
    """
    if 3 != len( location):
        raise TypeError( len( location), "location not a 3-array")
    _glulib.gluTessVertex( tess, location, data)


def gluUnProject( obj, model, proj, view):
    """
    obj - seq( 3, GLdouble)
    model - seq( gl.MODELVIEW_MATRIX, GLfloat)
    perspective - seq( gl.PROJECTION_MATRIX, GLfloat)
    view - seq( gl.VIEWPORT, GLfloat)
return - numpy.array ( 3, GLdouble)
    """
    win = numpy.zeros( 3, GLdouble)
    gluUnProjectv( obj, model, proj, view, win )
    return win[:3]


def gluUnProjectv( obj, model, proj, view, win):
    """
    obj - seq( 3, GLdouble)
    model - seq( gl.MODELVIEW_MATRIX, GLfloat)
    perspective - seq( gl.PROJECTION_MATRIX, GLfloat)
    view - seq( gl.VIEWPORT, GLfloat)
    win - numpy array( 3, GLdouble): gets the new coordinates
    """

    if 3 != len( obj):
        raise TypeError( len( obj), "obj not 3-array")
    model, proj, view = readModProjView( model, proj, view)
    _glulib.gluUnProject( obj[ 0], obj[ 1], obj[ 2], model, proj, view,
                         win[0:1], win[ 1:2], win[ 2:9])



__all__ = [
    'gluBuild1DMipMaps',
    'gluBuild2DMipMaps',
    'gluGetNurbsProperty',
    'gluGetTessProperty',
    'gluLoadSamplingMatrices',
    'gluNurbsCallback',
    'gluNurbsCurve',
    'gluNurbsSurface',
    'gluPickMatrix',
    'gluProject',
    'gluProjectv',
    'gluPwlCurve',
    'gluQuadricCallback',
    'gluScaleImage',
    'gluTessCallback',
    'gluTessNormaldv',
    'gluTessNormalfv',
    'gluTessVertex',
    'gluUnProject',
    'gluUnProjectv',
    ]
