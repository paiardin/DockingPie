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



__all__ = ["glGetMapdv",
    "glGetMapfv",
    "glGetMapiv",
    "glGetPixelMapfv",
    "glGetPixelMapuiv",
    "glGetPixelMapusv",
    "glGetPointerv",
    "glGetTexImage"]

from warnings import warn

from opengltk.extent import _gllib


class ToWrap( RuntimeWarning):
    pass


def glGetMapdv( target, query, v):
    """To Wrap: """
    #warn( desc, ToWrap)
    return _gllib.glGetMapdv( target, query, v)

def glGetMapfv( target, query, v):
    """To Wrap: """
    #warn( desc, ToWrap)
    return _gllib.glGetMapfv( target, query, v)

def glGetMapiv( target, query, v):
    """To Wrap: """
    #warn( desc, ToWrap)
    return _gllib.glGetMapiv( target, query, v)

def glGetPixelMapfv( map, values):
    """To Wrap: use glGet to get value sz to return"""
    #warn( desc, ToWrap)
    return _gllib.glGetPixelMapfv( map, values)

def glGetPixelMapuiv( map, values):
    """To Wrap: use glGet to get value sz to return"""
    #warn( desc, ToWrap)
    return _gllib.glGetPixelMapuiv( map, values)

def glGetPixelMapusv( map, values):
    """To Wrap: use glGet to get value sz to return"""
    #warn( desc, ToWrap)
    return _gllib.glGetPixelMapusv( map, values)

def glGetPointerv( pname, params):
    """To Wrap: convert cobject to appropriate Numeric array (need size...)"""
    #warn( desc, ToWrap)
    return _gllib.glGetPointerv( pname, params)

def glGetTexImage( target, level, format, type, pixels):
    """To Wrap: """
    #warn( desc, ToWrap)
    return _gllib.glGetTexImage( target, level, format, type, pixels)
