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


# $Header: /mnt/raid/services/cvs/DejaVu2/IndexedPolygons.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#
# $Id: IndexedPolygons.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#

import os
import numpy

from opengltk.OpenGL import GL
from opengltk.extent import _gllib

import DejaVu2
from DejaVu2.IndexedPolylines import IndexedPolylines
from DejaVu2.Geom import Geom
from DejaVu2.IndexedGeom import IndexedGeom
import datamodel, viewerConst
from viewerFns import checkKeywords
if hasattr(DejaVu2, 'enableVertexArray') is False:
    DejaVu2.enableVertexArray = False


def IndexedPolygonsFromFile(filename, name=None):
    if __debug__:
     if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
    """IndexedPolygons <- indexedPolygonsFromFile(filename, name=None)
Create an IndexedPolygons object from a file saved using the
geom.saveIndexedPolygons method. 
The extension of filename (if any) will be remove.
filename.vert and filename.face should be present.
The name arguemnt will be the name of the Geom object. If omitted, filename
will be used.
"""
    filename = os.path.splitext(filename)[0]

    from string import split
    f = open(filename+'.indpolvert')
    data = map(split, f.readlines())
    f.close()
    verts = map( lambda x: (float(x[0]), float(x[1]), float(x[2])), data )
    norms = map( lambda x: (float(x[3]), float(x[4]), float(x[5])), data )

    f = open(filename+'.indpolface')
    data = map(split, f.readlines())
    f.close()
    faces = []
    fnorms = []
    for line in data:
        faces.append( map( int, line[:-3]) )
        fnorms.append( map( float, line[-3:]) )
    if name is None:
        name = filename

    faces = numpy.array(faces)
    if min(faces.ravel())==1:
        faces = faces - 1
        
    pol = IndexedPolygons(name, vertices=verts, vnormals=norms, faces=faces)
    # FIXME face normals have to be set after the constructor else they are
    # overwritten with computed ones
    pol.Set(fnormals=fnorms, shading=GL.GL_FLAT)
    return pol


class IndexedPolygons(IndexedPolylines):
    """Set of polygons sharing vertices
"""

    def __init__(self, name=None, check=1, redo=1, replace=True, **kw):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        #print "IndexedPolygons.__init__", name

        kw['replace'] = replace
        #if DejaVu2.enableVertexArray is True:
        #    kw['vertexArrayFlag'] = True
        #    kw['immediateRendering'] = True

        apply( IndexedPolylines.__init__, (self, name, check), kw)


    def Set(self, check=1, redo=1, updateOwnGui=True, setCurrent=True, **kw):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """set data for this object: Set polygon's vertices, faces, normals or materials
check=1 : verify that all the keywords present can be handle by this func 
redo=1 : append self to viewer.objectsNeedingRedo
updateOwnGui=True : allow to update owngui at the end this func
"""
        redoFlags = apply( IndexedPolylines.Set, (self, check, 0), kw )

        if    hasattr(self, 'faceSet') \
          and len(self.faceSet.faces.array) > 0 \
          and len(self.faceSet.faces.array[0]) > 2 :
            self._PrimitiveType()

        return self.redoNow(redo, updateOwnGui, redoFlags, setCurrent)


    def delete(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.vertexSet.normals.Compute = None
        self.faceSet.normals.Compute = None

    
    def sortPoly(self, order=-1):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """None <- sortPoly(order=-1)
Sorts the geometry polygons according to z values of polygon's
geomtric centers. Order=-1 sorts by furthest z first, order=1 sorts
by closest z first"""
        # FIXME will not work with instance matrices
        mat = self.GetMatrix()
        mat = numpy.reshape(mat, (4,4))
        vt = self.vertexSet.vertices*mat
        if vt is None:
            return
        triv = numpy.take(vt, self.faceSet.faces.array, axis=0)
        trig = numpy.sum(triv,1)/3.
        trigz = trig[:,2]  #triangle's center of gravity z value
        
        ind = numpy.argsort(trigz) # sorted indices
        
        if len(self.faceSet.faces.array):
            faces = numpy.take(self.faceSet.faces.array, ind[::order], axis=0)
        
            if self.shading==GL.GL_FLAT: # we also need to re-arrange the
                                       # face normals
                if self.normals is None:
                    normals = None
                else:
                    if len(self.normals)>1:
                        normals = numpy.take(self.normals, ind[::order], axis=0)
                    else:
                        normals = self.normals
            else:
                normals = None

            # if materials are define per face with have to reshuffle materials too
            for i in range(5):
                if self.materials[GL.GL_FRONT].binding[i]==viewerConst.PER_PART:
                    prop = self.materials[GL.GL_FRONT].prop[i]
                    self.materials[GL.GL_FRONT].prop[i]= numpy.take(prop, ind[::order], axis=0)
                if self.materials[GL.GL_BACK].binding[i]==viewerConst.PER_PART:
                    prop = self.materials[GL.GL_BACK].prop[i]
                    self.materials[GL.GL_BACK].prop[i]= numpy.take(prop, ind[::order], axis=0)
            self.Set(faces=faces, fnormals=normals)


    def asIndexedPolygons(self, run=1, removeDupVerts=0, **kw):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Should return an IndexedPolygons object if this object can be
        represented using triangles, else return None. run=0 returns 1
        if this geom can be represented as an IndexedPolygon and None if not
        run=1 returns the IndexedPolygon object."""

        if run==0:
            return 1 # yes, I can be represented as IndexedPolygons
        
        if removeDupVerts==0:
            return self
        else:
            vu, fu = self.removeDuplicatedVertices()
            return IndexedPolygons('copy_'+self.name, vertices=vu, faces=fu,
                                   visible=1, 
                                   invertNormals=self.invertNormals)


    def writeToFile(self, filename):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """swriteToFile(filename)
Creates a .vert and a face file describing this indexed polygons geoemtry.
only vertices, and noramsl are saved the the .vert file (x y z nx ny nz)
and 0-based topoly in the .face file (i, j, k, ... ).
"""
        verts = self.getVertices()
        vnorms = self.getVNormals()
        fnorms = self.getFNormals()
        faces = self.getFaces()

        file = open(filename+'.indpolvert', 'w')
        map( lambda v, n, f=file: \
             f.write("%f %f %f %f %f %f\n"%tuple(tuple(v)+tuple(n))),
             verts, vnorms)
        file.close()

        file = open(filename+'.indpolface', 'w')
        for v, face in zip(fnorms, faces):
            map( lambda ind, f=file: f.write("%d "%ind ), face )
            file.write("%f %f %f\n"%tuple(v))
        file.close()

