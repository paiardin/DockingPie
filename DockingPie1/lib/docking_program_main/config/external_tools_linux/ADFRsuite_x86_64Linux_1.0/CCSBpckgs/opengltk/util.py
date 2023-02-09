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

'''util module
'''

__all__ = (
    #'glarray',
    'attachCurrentThread',
    'detachCurrentThread',
    )
import numpy
from extent._utillib import *
from extent import utillib, gllib

import threading

attachLock = threading.Lock()
attachedThread = None

def attachCurrentThread( blocking=1):
    global attachedThread
    attachLock.acquire( blocking)
    assert attachedThread is None
    utillib.attachCurrentThread()
    attachedThread = threading.currentThread()

def detachCurrentThread():
    global attachedThread
    assert attachedThread == threading.currentThread()
    utillib.detachCurrentThread()
    attachedThread = None
    attachLock.release()


def sizedict( types):
    lst = []
    for type in types:
        lst.append((numpy.ones(1, type).itemsize, type))
    lst.sort()
    result = {}
    for sz, typ in lst:
        if not result.has_key( sz):
            result[ sz] = typ
    return result

def importstripped( globs, module, preflist):
    import re
    if hasattr( module, '__all__'):
        symbs = module.__all__
    else:
        exppat = re.compile( '__')
        symbs = [x for x in dir( module) if not exppat.match( x)]
    subpat = re.compile( '^(%s)' % reduce( lambda x, y: '%s|%s' % (x, y),
                                           preflist[ 1:], preflist[ 0]))
    for symb in symbs:
        globs[ subpat.sub( '', symb)] = getattr( module, symb)


for mtyp in [
    'GLbitfield',
    'GLboolean',
    'GLubyte',
    'GLuint',
    'GLushort',
    ]:
    globals()[ mtyp] = sizedict( [numpy.uint8,
                                  numpy.uint16,
                                  numpy.uint32,
                                  ])[ getattr( utillib, 'sizeof_%s' % mtyp)]

for mtyp in [
    'GLbyte',
    'GLint',
    'GLshort',
    ]:
    globals()[ mtyp] = sizedict( [numpy.int8,
                                  numpy.int16,
                                  numpy.int32,
                                  ])[ getattr( utillib, 'sizeof_%s' % mtyp)]
for mtyp in [
    'GLclampd',
    'GLclampf',
    'GLdouble',
    'GLfloat',
    ]:
    globals()[ mtyp] = sizedict( [numpy.float32,
                                  numpy.float64,
                                  ])[ getattr( utillib, 'sizeof_%s' % mtyp)]
    for mtyp in ['GLenum', 'GLsizei']:
        globals()[ mtyp] = sizedict([numpy.uint32
                                    ])[getattr(utillib, 'sizeof_%s' % mtyp)] 

gltypmap = {
    GLbyte: gllib.GL_BYTE,
    GLubyte: gllib.GL_UNSIGNED_BYTE,
    GLshort: gllib.GL_SHORT,
    GLushort: gllib.GL_UNSIGNED_SHORT,
    GLint: gllib.GL_INT,
    GLuint: gllib.GL_UNSIGNED_INT,
    GLfloat: gllib.GL_FLOAT,
    GLdouble: gllib.GL_DOUBLE,
    }


revtypmap = {}
for ctype, econst in gltypmap.items():
    revtypmap[ econst] = ctype

# GL_[234]_BYTES not supported ... ?


glXXDim = {}

if __debug__:
    for param in (
        "GL_ACCUM_ALPHA_BITS",
        "GL_ACCUM_BLUE_BITS",
        "GL_ACCUM_GREEN_BITS",
        "GL_ACCUM_RED_BITS",
        "GL_ALPHA_BIAS",
        "GL_ALPHA_BITS",
        "GL_ALPHA_SCALE",
        "GL_ALPHA_TEST",
        "GL_ALPHA_TEST_FUNC",
        "GL_ALPHA_TEST_REF",
        "GL_ATTRIB_STACK_DEPTH",
        "GL_AUTO_NORMAL",
        "GL_AUX_BUFFERS",
        "GL_BLEND",
        "GL_BLEND_DST",
        "GL_BLEND_EQUATION_EXT",
        "GL_BLEND_SRC",
        "GL_BLUE_BIAS",
        "GL_BLUE_BITS",
        "GL_BLUE_SCALE",
        "GL_CLIENT_ATTRIB_STACK_DEPTH",
        "GL_CLIP_PLANE0",
        "GL_CLIP_PLANE1",
        "GL_CLIP_PLANE2",
        "GL_CLIP_PLANE3",
        "GL_CLIP_PLANE4",
        "GL_CLIP_PLANE5",
        "GL_COLOR_ARRAY",
        "GL_COLOR_ARRAY_SIZE",
        "GL_COLOR_ARRAY_STRIDE",
        "GL_COLOR_ARRAY_TYPE",
        "GL_COLOR_LOGIC_OP",
        "GL_COLOR_MATERIAL",
        "GL_COLOR_MATERIAL_FACE",
        "GL_COLOR_MATERIAL_PARAMETER",
        "GL_CONSTANT_ATTENUATION",
        "GL_CULL_FACE",
        "GL_CULL_FACE_MODE",
        "GL_CURRENT_INDEX",
        "GL_CURRENT_RASTER_DISTANCE",
        "GL_CURRENT_RASTER_INDEX",
        "GL_CURRENT_RASTER_POSITION_VALID",
        "GL_DECAL",
        "GL_DEPTH_BIAS",
        "GL_DEPTH_BITS",
        "GL_DEPTH_CLEAR_VALUE",
        "GL_DEPTH_FUNC",
        "GL_DEPTH_SCALE",
        "GL_DEPTH_TEST",
        "GL_DEPTH_WRITEMASK",
        "GL_DITHER",
        "GL_DOUBLEBUFFER",
        "GL_DRAW_BUFFER",
        "GL_EDGE_FLAG",
        "GL_EDGE_FLAG_ARRAY",
        "GL_EDGE_FLAG_ARRAY_STRIDE",
        "GL_FOG",
        "GL_FOG_DENSITY",
        "GL_FOG_END",
        "GL_FOG_HINT",
        "GL_FOG_INDEX",
        "GL_FOG_MODE",
        "GL_FOG_START",
        "GL_FRONT_FACE",
        "GL_GREEN_BIAS",
        "GL_GREEN_BITS",
        "GL_GREEN_SCALE",
        "GL_INDEX_ARRAY",
        "GL_INDEX_ARRAY_STRIDE",
        "GL_INDEX_ARRAY_TYPE",
        "GL_INDEX_BITS",
        "GL_INDEX_CLEAR_VALUE",
        "GL_INDEX_LOGIC_OP",
        "GL_INDEX_MODE",
        "GL_INDEX_OFFSET",
        "GL_INDEX_SHIFT",
        "GL_INDEX_WRITEMASK",
        "GL_LIGHT0",
        "GL_LIGHT1",
        "GL_LIGHT2",
        "GL_LIGHT3",
        "GL_LIGHT4",
        "GL_LIGHT5",
        "GL_LIGHT6",
        "GL_LIGHT7",
        "GL_LIGHTING",
        "GL_LIGHT_MODEL_LOCAL_VIEWER",
        "GL_LIGHT_MODEL_TWO_SIDE",
        "GL_LINEAR_ATTENUATION",
        "GL_LINE_SMOOTH",
        "GL_LINE_SMOOTH_HINT",
        "GL_LINE_STIPPLE",
        "GL_LINE_STIPPLE_PATTERN",
        "GL_LINE_STIPPLE_REPEAT",
        "GL_LINE_WIDTH",
        "GL_LINE_WIDTH_GRANULARITY",
        "GL_LIST_BASE",
        "GL_LIST_INDEX",
        "GL_LIST_MODE",
        "GL_LOGIC_OP_MODE",
        "GL_MAP1_COLOR_4",
        "GL_MAP1_GRID_SEGMENTS",
        "GL_MAP1_INDEX",
        "GL_MAP1_NORMAL",
        "GL_MAP1_TEXTURE_COORD_1",
        "GL_MAP1_TEXTURE_COORD_2",
        "GL_MAP1_TEXTURE_COORD_3",
        "GL_MAP1_TEXTURE_COORD_4",
        "GL_MAP1_VERTEX_3",
        "GL_MAP1_VERTEX_4",
        "GL_MAP2_COLOR_4",
        "GL_MAP2_INDEX",
        "GL_MAP2_NORMAL",
        "GL_MAP2_TEXTURE_COORD_1",
        "GL_MAP2_TEXTURE_COORD_2",
        "GL_MAP2_TEXTURE_COORD_3",
        "GL_MAP2_TEXTURE_COORD_4",
        "GL_MAP2_VERTEX_3",
        "GL_MAP2_VERTEX_4",
        "GL_MAP_STENCIL",
        "GL_MATRIX_MODE",
        "GL_MAX_CLIENT_ATTRIB_STACK_DEPTH",
        "GL_MAX_ATTRIB_STACK_DEPTH",
        "GL_MAX_CLIP_PLANES",
        "GL_MAX_EVAL_ORDER",
        "GL_MAX_LIGHTS",
        "GL_MAX_LIST_NESTING",
        "GL_MAX_MODELVIEW_STACK_DEPTH",
        "GL_MAX_NAME_STACK_DEPTH",
        "GL_MAX_PIXEL_MAP_TABLE",
        "GL_MAX_PROJECTION_STACK_DEPTH",
        "GL_MAX_TEXTURE_SIZE",
        "GL_MAX_TEXTURE_STACK_DEPTH",
        "GL_MODELVIEW_STACK_DEPTH",
        "GL_MODULATE",
        "GL_NAME_STACK_DEPTH",
        "GL_NORMAL_ARRAY",
        "GL_NORMAL_ARRAY_STRIDE",
        "GL_NORMAL_ARRAY_TYPE",
        "GL_NORMALIZE",
        "GL_PACK_ALIGNMENT",
        "GL_PACK_LSB_FIRST",
        "GL_PACK_ROW_LENGTH",
        "GL_PACK_SKIP_PIXELS",
        "GL_PACK_SKIP_ROWS",
        "GL_PACK_SWAP_BYTES",
        "GL_PERSPECTIVE_CORRECTION_HINT",
        "GL_PIXEL_MAP_A_TO_A_SIZE",
        "GL_PIXEL_MAP_B_TO_B_SIZE",
        "GL_PIXEL_MAP_G_TO_G_SIZE",
        "GL_PIXEL_MAP_I_TO_A_SIZE",
        "GL_PIXEL_MAP_I_TO_B_SIZE",
        "GL_PIXEL_MAP_I_TO_G_SIZE",
        "GL_PIXEL_MAP_I_TO_I_SIZE",
        "GL_PIXEL_MAP_I_TO_R_SIZE",
        "GL_PIXEL_MAP_R_TO_R_SIZE",
        "GL_PIXEL_MAP_S_TO_S_SIZE",
        "GL_POINT_SIZE",
        "GL_POINT_SIZE_GRANULARITY",
        "GL_POINT_SMOOTH",
        "GL_POINT_SMOOTH_HINT",
        "GL_POLYGON_OFFSET_FACTOR",
        "GL_POLYGON_OFFSET_UNITS",
        "GL_POLYGON_OFFSET_FILL",
        "GL_POLYGON_OFFSET_LINE",
        "GL_POLYGON_OFFSET_POINT",
        "GL_POLYGON_SMOOTH",
        "GL_POLYGON_SMOOTH_HINT",
        "GL_POLYGON_STIPPLE",
        "GL_PROJECTION_STACK_DEPTH",
        "GL_READ_BUFFER",
        "GL_REPLACE",
        "GL_QUADRATIC_ATTENUATION",
        "GL_RED_BIAS",
        "GL_RED_BITS",
        "GL_RED_SCALE",
        "GL_RENDER_MODE",
        "GL_RGBA_MODE",
        "GL_SCISSOR_TEST",
        "GL_SHADE_MODEL",
        "GL_SHININESS",
        "GL_SPOT_EXPONENT",
        "GL_SPOT_CUTOFF",
        "GL_STENCIL_BITS",
        "GL_STENCIL_CLEAR_VALUE",
        "GL_STENCIL_FAIL",
        "GL_STENCIL_FUNC",
        "GL_STENCIL_PASS_DEPTH_FAIL",
        "GL_STENCIL_PASS_DEPTH_PASS",
        "GL_STENCIL_REF",
        "GL_STENCIL_TEST",
        "GL_STENCIL_VALUE_MASK",
        "GL_STENCIL_WRITEMASK",
        "GL_STEREO",
        "GL_SUBPIXEL_BITS",
        "GL_TEXTURE_1D",
        "GL_TEXTURE_1D_BINDING",
        "GL_TEXTURE_2D",
        "GL_TEXTURE_2D_BINDING",
        "GL_TEXTURE_COORD_ARRAY",
        "GL_TEXTURE_COORD_ARRAY_SIZE",
        "GL_TEXTURE_COORD_ARRAY_STRIDE",
        "GL_TEXTURE_COORD_ARRAY_TYPE",
        "GL_TEXTURE_ENV_MODE",
        "GL_TEXTURE_GEN_MODE",
        "GL_TEXTURE_GEN_Q",
        "GL_TEXTURE_GEN_R",
        "GL_TEXTURE_GEN_S",
        "GL_TEXTURE_GEN_T",
        "GL_TEXTURE_MAG_FILTER",
        "GL_TEXTURE_MIN_FILTER",
        "GL_TEXTURE_PRIORITY",
        "GL_TEXTURE_RESIDENT",
        "GL_TEXTURE_STACK_DEPTH",
        "GL_TEXTURE_WRAP_S",
        "GL_TEXTURE_WRAP_T",
        "GL_UNPACK_ALIGNMENT",
        "GL_UNPACK_LSB_FIRST",
        "GL_UNPACK_ROW_LENGTH",
        "GL_UNPACK_SKIP_PIXELS",
        "GL_UNPACK_SKIP_ROWS",
        "GL_UNPACK_SWAP_BYTES",
        "GL_VERTEX_ARRAY",
        "GL_VERTEX_ARRAY_SIZE",
        "GL_VERTEX_ARRAY_STRIDE",
        "GL_VERTEX_ARRAY_TYPE",
        "GL_ZOOM_X",
        "GL_ZOOM_Y",
        # add glGetTexLevelParameter pnames (all on 1)
        ): # not exhausitve... see Xxdim
        if hasattr(gllib, param):
            glXXDim[ getattr(gllib, param)] = 1

for param in (
    'GL_DEPTH_RANGE',
    'GL_LINE_WIDTH_RANGE',
    'GL_MAP1_GRID_DOMAIN',
    'GL_MAP2_GRID_SEGMENTS',
    'GL_MAX_VIEWPORT_DIMS',
    'GL_POINT_SIZE_RANGE',
    'GL_POLYGON_MODE'
    ):
    if hasattr(gllib, param):
        assert not glXXDim.has_key( param), (param, glXXDim[ param])
        glXXDim[ getattr(gllib, param)] = 2

for param in (
    'GL_COLOR_INDEXES',
    'GL_CURRENT_NORMAL',
    'GL_SPOT_DIRECTION',
    ):
    if hasattr(gllib, param):
        assert not glXXDim.has_key( param), (param, glXXDim[ param])
        glXXDim[ getattr(gllib, param)] = 3

for param in (
    'GL_ACCUM_CLEAR_VALUE',
    'GL_AMBIENT',
    'GL_BLEND_COLOR_EXT',
    'GL_COLOR_CLEAR_VALUE',
    'GL_COLOR_WRITEMASK',
    'GL_CURRENT_COLOR',
    'GL_CURRENT_RASTER_COLOR',
    'GL_CURRENT_RASTER_POSITION',
    'GL_CURRENT_RASTER_TEXTURE_COORDS',
    'GL_CURRENT_TEXTURE_COORDS',
    'GL_DIFFUSE',
    'GL_EMISSION',
    'GL_EYE_PLANE',
    'GL_FOG_COLOR',
    'GL_LIGHT_MODEL_AMBIENT',
    'GL_MAP2_GRID_DOMAIN',
    'GL_OBJECT_PLANE',
    'GL_POSITION',
    'GL_SCISSOR_BOX',
    'GL_SPECULAR',
    'GL_TEXTURE_ENV_COLOR',
    'GL_TEXTURE_BORDER_COLOR',
    'GL_VIEWPORT',
    ):
    if hasattr(gllib, param):
        assert not glXXDim.has_key( param), (param, glXXDim[ param])
        glXXDim[ getattr(gllib, param)] = 4

for param in (
    'GL_MODELVIEW_MATRIX',
    'GL_PROJECTION_MATRIX',
    'GL_TEXTURE_MATRIX',
    ):
    if hasattr(gllib, param):
        assert not glXXDim.has_key( param), (param, glXXDim[ param])
        glXXDim[ getattr(gllib, param)] = 16


class Xxdim:
    def __getitem__( self, key):
        if __debug__:
            try:
                return glXXDim[ key]
            except KeyError:
                from warnings import warn
                warn( 'key %s not in glXXXDim, return default %i' % (key,
                                                                     default),
                      RuntimeWarning)
                return 1
        else:
            return glXXDim.get( key, 1)

glGetXXDim = Xxdim()


mapDim = {
    gllib.GL_MAP1_COLOR_4: 4,
    gllib.GL_MAP1_INDEX: 1,
    gllib.GL_MAP1_NORMAL: 3,
    gllib.GL_MAP1_TEXTURE_COORD_1: 1,
    gllib.GL_MAP1_TEXTURE_COORD_2: 2,
    gllib.GL_MAP1_TEXTURE_COORD_3: 3,
    gllib.GL_MAP1_TEXTURE_COORD_4: 4,
    gllib.GL_MAP1_VERTEX_3: 3,
    gllib.GL_MAP1_VERTEX_4: 4,
    gllib.GL_MAP2_COLOR_4: 4,
    gllib.GL_MAP2_INDEX: 1,
    gllib.GL_MAP2_NORMAL: 3,
    gllib.GL_MAP2_TEXTURE_COORD_1: 1,
    gllib.GL_MAP2_TEXTURE_COORD_2: 2,
    gllib.GL_MAP2_TEXTURE_COORD_3: 3,
    gllib.GL_MAP2_TEXTURE_COORD_4: 4,
    gllib.GL_MAP2_VERTEX_3: 3,
    gllib.GL_MAP2_VERTEX_4: 4,
    }

def readModProjView( model, projection, view):
    if len( model) != glGetXXDim[ gllib.GL_MODELVIEW_MATRIX]:
        raise TypeError( (len( model),
                          glGetXXDim[ gllib.GL_MODELVIEW_MATRIX]),
                         'len( model) wrong')
    if len( projection) != glGetXXDim[ gllib.GL_PROJECTION_MATRIX]:
        raise TypeError( (len( projection),
                          glGetXXDim[ gllib.GL_PROJECTION_MATRIX]),
                         'len( projection) wrong')
    if len( view) != glGetXXDim[ gllib.GL_VIEWPORT]:
        raise TypeError( (len( view), glGetXXDim[ gllib.GL_VIEWPORT]),
                         'len( view) wrong')
    return model, projection, view
