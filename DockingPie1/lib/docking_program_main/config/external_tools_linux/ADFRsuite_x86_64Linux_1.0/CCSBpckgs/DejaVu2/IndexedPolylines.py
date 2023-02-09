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

from inspect import isclass

from DejaVu2.IndexedGeom import IndexedGeom
import datamodel
import viewerConst
from DejaVu2.viewerFns import checkKeywords


class IndexedPolylines(IndexedGeom):
    """Set of lines sharing vertices"""

    def __init__(self, name=None, check=1, **kw):
        #print "IndexedPolylines.__init__", name

        #kw['vertexArrayFlag'] = True
        #kw['immediateRendering'] = True
        apply( IndexedGeom.__init__, (self, name, check), kw)


    def Set(self, check=1, redo=1, updateOwnGui=True, **kw):
        """set data for this object: Set polygon's vertices, faces, normals or materials
check=1 : verify that all the keywords present can be handle by this func 
redo=1 : append self to viewer.objectsNeedingRedo
updateOwnGui=True : allow to update owngui at the end this func
"""
        redoFlags = apply( IndexedGeom.Set, (self, check, 0), kw )

        if    hasattr(self, 'faceSet') \
          and len(self.faceSet.faces.array) > 0 \
          and len(self.faceSet.faces.array[0]) > 2 :
            # register functions to compute normals
            self.VertexNormalFunction(self.ComputeVertexNormals)
            self.vertexSet.normals.ComputeMode( viewerConst.AUTO )
            self.FaceNormalFunction(self.ComputeFaceNormals)
            self.faceSet.normals.ComputeMode( viewerConst.AUTO )
            self.GetNormals()
            from opengltk.OpenGL.GL import GL_LINE_STRIP
            self._PrimitiveType(type=GL_LINE_STRIP)

        return self.redoNow(redo, updateOwnGui, redoFlags)


    def Add(self, check=1, redo=1, **kw):
	"""Add polygon's vertices, faces, normals or materials"""

        if __debug__:
            if check:
                apply( checkKeywords, (self.name,self.keywords), kw)

        apply( IndexedGeom.Add, (self, 0, 0), kw )

        if self.viewer and redo:
            if self.redoDspLst:
                self.viewer.objectsNeedingRedo[self] = None
