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

"""glx wrappers
"""

import numpy
from opengltk.extent import _glxlib


def glXChooseVisual( dpy, screen, attribpairs):
    """
    dpy - Display*
    screen - int
    attribpairs - seq( (attribute, value))
    """
    from types import IntType
    larray = numpy.zeros( 2 * len( attribpairs) + 1, numpy.int32)
    idx = 0
    for assg in attribpairs:
        if isinstance( assg, IntType):
            larray[ idx] = assg
        else:
            larray[ idx] = assg[ 0]
            idx += 1
            larray[ idx] = assg[ 1]
        idx += 1

    return _glxlib.glXChooseVisual( dpy, screen, larray)


def glXGetConfig( dpy, vis, attrib):
    """
    dpy - Display*
    vis - XVisualInfo*
    attrib - int
return - int
    """
    value = numpy.zeros( 1, numpy.int32)
    res = _glxlib.glXGetConfig( dpy, vis, attrib, value)
    if res:
        from exception import Glxerror
        raise Glxerror( res)
    return value[ 0]


def glXQueryExtension( dpy):
    """
    dpy - Display*
return - bool, int, int: support, errorBase, eventBase
"""
    errorBase = numpy.zeros( 1, numpy.int32)
    eventBase = numpy.zeros( 1, numpy.int32)
    res = _glxlib.glXQueryExtension( dpy, errorBase, eventBase)
    return res, errorBase[ 0], eventBase[ 0]


def glXQueryVersion( dpy):
    """
    dpy - Display*
return - bool, int, int: support, major, minor
"""
    major = numpy.zeros( 1, numpy.int32)
    minor = numpy.zeros( 1, numpy.int32)
    res = _glxlib.glXQueryVersion( dpy, major, minor)
    return res, major[ 0], minor[ 0]


def XOpenDisplay( name=None):
    if name is None:
        from os import environ
        name = environ[ 'DISPLAY']
    return _glxlib.XOpenDisplay( name)

__all__ = [
    'glXChooseVisual',
    'glXGetConfig',
    'glXQueryExtension',
    'glXQueryVersion',
    'XOpenDisplay',
    ]

