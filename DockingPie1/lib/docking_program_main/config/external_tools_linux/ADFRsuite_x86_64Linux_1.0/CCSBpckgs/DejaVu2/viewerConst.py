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
# $Header: /mnt/raid/services/cvs/DejaVu2/viewerConst.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#
# $Id: viewerConst.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $ #
import numpy
from opengltk.OpenGL import GL

FPRECISION = numpy.float32
IPRECISION = numpy.int32

INHERIT = -1

NO = False
YES = True

# binding modes for normals and colors
OVERALL = 10
PER_VERTEX = 11
PER_PART = 12
PER_INSTANCE = 13

## FIXME these constances are defined here AND in Materials.py !
propConst = ( GL.GL_AMBIENT, GL.GL_DIFFUSE,
	      GL.GL_EMISSION, GL.GL_SPECULAR, GL.GL_SHININESS )
propNum = { GL.GL_AMBIENT:0, GL.GL_DIFFUSE:1,
            GL.GL_EMISSION:2, GL.GL_SPECULAR:3, 
            GL.GL_SHININESS:4
            }
# constants for properties status and computation
NONE = 20
UNKNOWN = 21
TOO_MANY = 22
TOO_FEW = 23
COMPUTED = 24
SET = 25

# constants for properties computation
NO_COMPUTATION = 30
AUTO = 31

#drawBB values
NO = NO
ONLY = 41   # BB only
WITHOBJECT = 42   # BB and object
BB_MODES = (NO, ONLY, WITHOBJECT)

OUTLINED = 50

Front_POLYGON_MODES_keys = ('point', 'line', 'fill', 'outlined', 'inherit')
Front_POLYGON_MODES_values = (GL.GL_POINT, GL.GL_LINE, GL.GL_FILL, OUTLINED, INHERIT)
Back_POLYGON_MODES_keys = ('point', 'line', 'fill', 'outlined', 'inherit', 'as front')
Back_POLYGON_MODES_values = (GL.GL_POINT, GL.GL_LINE, GL.GL_FILL, OUTLINED, INHERIT, GL.GL_FRONT_AND_BACK)
POLYGON_MODES = dict(zip(Back_POLYGON_MODES_keys, Back_POLYGON_MODES_values))
POLYGON_MODES_REV = dict(zip(Back_POLYGON_MODES_values, Back_POLYGON_MODES_keys))
#POLYGON_MODES['as front'] = GL.GL_FRONT_AND_BACK
#POLYGON_MODES[GL.GL_FRONT_AND_BACK] = 'as front'

POINTS_PRIMITIVES = (GL.GL_POINTS, )
LINES_PRIMITIVES = (GL.GL_LINES, GL.GL_LINE_LOOP,
		    GL.GL_LINE_STRIP )
POLYGON_PRIMITIVES = (GL.GL_TRIANGLES, GL.GL_TRIANGLE_STRIP,
		      GL.GL_TRIANGLE_FAN, GL.GL_QUADS,
		      GL.GL_QUAD_STRIP, GL.GL_POLYGON)
PRIMITIVES = POINTS_PRIMITIVES + LINES_PRIMITIVES + POLYGON_PRIMITIVES
MINIMUM_LENGTH = {
 GL.GL_POINTS:1,
 GL.GL_LINES:2, GL.GL_LINE_LOOP:3, GL.GL_LINE_STRIP:3,
 GL.GL_TRIANGLES:3, GL.GL_TRIANGLE_STRIP:4, GL.GL_TRIANGLE_FAN:4,
 GL.GL_QUADS:4, GL.GL_QUAD_STRIP:5, GL.GL_POLYGON:3
}

SHADINGS = {'flat':GL.GL_FLAT, 'smooth':GL.GL_SMOOTH, 'inherit':INHERIT}
SHADINGS_REV = {GL.GL_FLAT:'flat', GL.GL_SMOOTH:'smooth', INHERIT:'inherit'}

CULLINGS_keys = ('none', 'back', 'front', 'front_and_back', 'inherit')
CULLINGS_values = (GL.GL_NONE, GL.GL_BACK, GL.GL_FRONT,
                   GL.GL_FRONT_AND_BACK, INHERIT)
CULLINGS = dict(zip(CULLINGS_keys, CULLINGS_values))
CULLINGS_REV = dict(zip(CULLINGS_values, CULLINGS_keys))


srcBFnames_keys = ( 'GL_ZERO', 'GL_ONE',
               'GL_DST_COLOR', 'GL_ONE_MINUS_DST_COLOR',
               'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA',
               'GL_DST_ALPHA', 'GL_ONE_MINUS_DST_ALPHA',
               'GL_SRC_ALPHA_SATURATE' )
srcBFnames_values = ( GL.GL_ZERO, GL.GL_ONE,
               GL.GL_DST_COLOR, GL.GL_ONE_MINUS_DST_COLOR,
               GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA,
               GL.GL_DST_ALPHA, GL.GL_ONE_MINUS_DST_ALPHA,
               GL.GL_SRC_ALPHA_SATURATE )
srcBFnames = dict(zip(srcBFnames_keys, srcBFnames_values))
srcBFnames_REV = dict(zip(srcBFnames_values, srcBFnames_keys))

