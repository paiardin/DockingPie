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
# $Header: /mnt/raid/services/cvs/DejaVu2/Polylines.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#
# $Id: Polylines.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#

from opengltk.OpenGL import GL
from opengltk.extent import _gllib as gllib

from Geom import Geom
import datamodel, viewerConst, types
from viewerFns import checkKeywords

import DejaVu2

class Polylines(Geom):
    """Class for sets of lines"""

    keywords = Geom.keywords + [
        'type',
        ]

    def __init__(self, name=None, check=1, **kw):

        if not kw.get('shape'):
            kw['shape'] = (0,0,3)    # default shape for line's vertex set

        if kw.has_key('type') is False:
            kw['type'] = GL.GL_LINES
        
        apply( Geom.__init__, (self, name, check), kw)


    def Set(self, check=1, redo=1, updateOwnGui=True, **kw):
        """set data for this object:
check=1 : verify that all the keywords present can be handle by this func 
redo=1 : append self to viewer.objectsNeedingRedo
updateOwnGui=True : allow to update owngui at the end this func
"""
        redoFlags = apply( Geom.Set, (self, check, 0), kw )

        type = kw.pop('type', None)
        if type is not None:
            self._PrimitveType(type)

        return self.redoNow(redo, updateOwnGui, redoFlags)


    def MaterialBindingMode(self, num, face=GL.GL_FRONT, mode=None):
	"""Figure out how materials should be used to color the object
        we overwrite the method from Geom because the number of parts
        corresponds to the nuber of lines which is
        vertexSet.vertices.array.shape[0]"""

        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        if num is None:
            num = self.materials[GL.GL_FRONT].diff
        elif type(num) is types.StringType:
            num = getattr(self.materials[GL.GL_FRONT], num)

        if face is None:
            face = GL.GL_FRONT

        f = face
        if face == GL.GL_FRONT_AND_BACK:
            f = GL.GL_FRONT
            self.frontAndBack = True

        lNumOfMaterials = self.materials[f].prop[num].shape[0]

        if mode is None:
            if lNumOfMaterials == len(self.vertexSet.vertices):
                self.materials[f].binding[num] = viewerConst.PER_VERTEX
            elif lNumOfMaterials == self.vertexSet.vertices.array.shape[0]:
                self.materials[f].binding[num] = viewerConst.PER_PART
            elif lNumOfMaterials == len(self.instanceMatricesFortran):
                self.materials[f].binding[num] = viewerConst.PER_INSTANCE
            elif lNumOfMaterials > 0:
                self.materials[f].binding[num] = viewerConst.OVERALL
            else:
                self.materials[f].binding[num] = viewerConst.INHERIT
        else: # a mode is  requested
            if mode==viewerConst.PER_VERTEX and lNumOfMaterials >= len(self.vertexSet.vertices):
                self.materials[f].binding[num] = viewerConst.PER_VERTEX
            elif lNumOfMaterials == self.vertexSet.vertices.array.shape[0]:
                self.materials[f].binding[num] = viewerConst.PER_PART
            elif mode==viewerConst.OVERALL and lNumOfMaterials >= 1:
                self.materials[f].binding[num] = viewerConst.OVERALL
            elif mode==viewerConst.PER_INSTANCE and lNumOfMaterials >= len(self.instanceMatricesFortran):
                self.materials[f].binding[num] = viewerConst.PER_INSTANCE
            else:
                self.materials[f].binding[num] = viewerConst.INHERIT


    def __repr__(self):
        return '<%s> %s with %d lines' % (self.__class__,
                          self.name, len(self.vertexSet) )


    def _PrimitveType(self, type):
	"""Find out out what type of lines
"""
	assert type in viewerConst.LINES_PRIMITIVES, type
	self.fixedLength = self.vertexSet.vertices.ashape[1]
	if self.fixedLength == 2 and (type is None or type==GL.GL_LINES):
	    self.primitiveType = GL.GL_LINES
	elif type:
	    if self.fixedLength > viewerConst.MINIMUM_LENGTH[type]:
		self.primitiveType = type
	    else:
		self.fixedLength = viewerConst.NO
		self.primitiveType = GL.GL_LINE_STRIP
	else:
            self.primitiveType = GL.GL_LINE_STRIP


    def Draw(self):
        c = self.vertexSet.vertices.array
        if len(c)==0: return
        GL.glDisable(GL.GL_LIGHTING)

        if self.materials[GL.GL_FRONT]:
            mat = self.materials[GL.GL_FRONT]
            binding = self.materials[GL.GL_FRONT].binding[mat.diff]
        else:
            binding = viewerConst.OVERALL        

        if binding == viewerConst.INHERIT:
            for i in xrange(c.shape[0]): #loop over lines
                GL.glPushName(i)
                GL.glBegin(self.primitiveType)
                for v in c[i]:
                    #GL.glVertex3dv(v.tolist())
                    #gllib.glVertex3dv(v)
                    gllib.glVertex3fv(v)
                GL.glEnd()
                GL.glPopName()

        elif binding == viewerConst.OVERALL:
            if self.materials[GL.GL_FRONT]:
                col = mat.prop[mat.diff][0]
            GL.glColor4fv(col)
            for i in xrange(c.shape[0]): #loop over lines
                GL.glPushName(i)
                GL.glBegin(self.primitiveType)
                for v in c[i]:
                    #GL.glVertex3dv(v.tolist())
                    #gllib.glVertex3dv(v)
                    gllib.glVertex3fv(v)
                GL.glEnd()
                GL.glPopName()

        elif binding == viewerConst.PER_VERTEX:
            if self.materials[GL.GL_FRONT]:
                col = mat.prop[mat.diff]
            vi = 0
            for i in xrange(c.shape[0]):
                GL.glPushName(i)
                GL.glBegin(self.primitiveType)
                for v in c[i]:
                    GL.glColor4fv(col[vi])
                    vi = vi + 1
                    #GL.glVertex3dv(v.tolist())
                    #gllib.glVertex3dv(v)
                    gllib.glVertex3fv(v)
                GL.glEnd()
                GL.glPopName()

        elif binding == viewerConst.PER_PART: # i.e. line
            if self.materials[GL.GL_FRONT]:
                col = mat.prop[mat.diff]
            for i in xrange(c.shape[0]):
                GL.glColor4fv(col[i])
                GL.glPushName(i)
                GL.glBegin(self.primitiveType)
                for v in c[i]:
                    #GL.glVertex3dv(v.tolist())
                    #gllib.glVertex3dv(v)
                    gllib.glVertex3fv(v)
                GL.glEnd()
                GL.glPopName()

        #glEnable(GL_LIGHTING)
        if self.viewer is not None:
            self.viewer.enableOpenglLighting()
        return True
