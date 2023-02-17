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
# $Header: /mnt/raid/services/cvs/DejaVu2/Points.py,v 1.3.2.1 2017/07/13 22:28:32 annao Exp $
#
# $Id: Points.py,v 1.3.2.1 2017/07/13 22:28:32 annao Exp $
#

from opengltk.OpenGL import GL
from opengltk.extent.utillib import glDrawIndexedGeom

#from Geom import Geom
from .IndexedGeom import IndexedGeom
import datamodel, viewerConst
import numpy
from viewerFns import checkKeywords

class Points(IndexedGeom):
    """Class for sets of spheres"""

    keywords = IndexedGeom.keywords + [
        'centers',
        ]

    def __init__(self, name=None, check=1, **kw):

        v = kw.get( 'centers')
        if v:
            kw['vertices'] = v     # rename centers in vertices for Geom.__init
        elif not kw.get('shape'):
            kw['shape'] = (0,3)    # default shape for sphere set

        if not kw.get('frontPolyMode'):
            kw['frontPolyMode'] = GL.GL_POINT
        if not kw.get('inheritFrontPolyMode'):
            kw['inheritFrontPolyMode'] = viewerConst.NO
        if not kw.get('lighting'):
            kw['lighting'] = viewerConst.NO
        if not kw.get('inheritLighting'):
            kw['inheritLighting'] = viewerConst.NO

        IndexedGeom.__init__(self, name, check, **kw)
        
        assert len(self.vertexSet.vertices.ashape)==2

        self._modified = False


    def Set(self, check=1, redo=1, updateOwnGui=True, **kw):
        """set data for this object
check=1 : verify that all the keywords present can be handle by this func 
redo=1 : append self to viewer.objectsNeedingRedo
updateOwnGui=True : allow to update owngui at the end this func
"""
        # Exceptionnaly this has to be before the call to Geom.Set
        v = kw.get( 'centers')
        if v:
            kw['vertices'] = v     # rename centers in vertices for Geom.__init

        redoFlags = apply( IndexedGeom.Set, (self, check, 0), kw)

        return self.redoNow(redo, updateOwnGui, redoFlags)


    def Add(self, check=1, redo=1, **kw):
	"""Add spheres"""

        if __debug__:
            if check:
                apply( checkKeywords, (self.name,self.keywords), kw)

	v = kw.get( 'centers')
	if v:
	    kw['vertices'] = v     # rename centers in vertices for Geom.__init
            self.redoDspLst = 1

        apply(IndexedGeom.Add, (self, 0, 0), kw)
        
        if self.viewer and redo:
            if self.redoDspLst:
                self.viewer.objectsNeedingRedo[self] = None
#                self.RedoDisplayList()


    def Draw(self):

        centers = self.vertexSet.vertices.array
        if len(centers)==0:
            return

        if len(self.faceSet.faces)==0:
            faces = [[x] for x in xrange(len(centers))]
        else:
            faces = self.faceSet.faces.array

        if self.materials[GL.GL_FRONT] and \
               not self.inheritMaterial:
            fpProp = self.materials[GL.GL_FRONT].prop[:5]
            fpBind = self.materials[GL.GL_FRONT].binding[:5]
        else:
            fpProp = None
            fpBind = None

        if self.materials[GL.GL_BACK] and \
               not self.inheritMaterial:
            bpProp = self.materials[GL.GL_BACK].prop[:5]
            bpBind = self.materials[GL.GL_BACK].binding[:5]
        else:
            bpProp = None
            bpBind = None

        texCoords = None

        if self.normals is None:
            GL.glDisable(GL.GL_LIGHTING)

        status = glDrawIndexedGeom(
            GL.GL_POINTS,
            self.vertexSet.vertices.array,
            faces,
            self.normals,
            texCoords,
            fpProp, bpProp, fpBind, bpBind,
            self.frontAndBack, 1)
        return status


# FIXME picking is probably not working onthis geometry
class CrossSet(Points):

    keywords = Points.keywords + [
        'offset',
        ]
    
    def __init__(self, name=None, check=1, **kw):

        v = kw.get( 'centers')
        if v:
            kw['vertices'] = v     # rename centers in vertices for Geom.__init

        # has to be fone before Points.__init__
        offset = kw.get( 'offset')
        if offset != None: self.offset = offset
        else: self.offset = 0.3

#        if not kw.get('transparent'):
#            kw['transparent'] = viewerConst.YES

        apply( Points.__init__, (self, name, check), kw)

        self._modified = False

        
    def Set(self, check=1, redo=1, updateOwnGui=True, **kw):
        """set data for this object
check=1 : verify that all the keywords present can be handle by this func 
redo=1 : append self to viewer.objectsNeedingRedo
updateOwnGui=True : allow to update owngui at the end this func
"""

        # Exceptionnaly this has to be before the call to Geom.Set
        v = kw.get( 'centers')
        if v:
            kw['vertices'] = v     # rename centers in vertices for Geom.__init
            
        redoFlags = apply( Points.Set, (self, check, 0), kw)

        offset = kw.get( 'offset')
        if offset is not None:
            self.offset = offset
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        return self.redoNow(redo, updateOwnGui, redoFlags)


    def getCoordsAndFaces(self, c):
        c = numpy.array(c)
        x1 = c - (self.offset, 0, 0)
        x2 = c + (self.offset, 0, 0)
        y1 = c - (0, self.offset, 0)
        y2 = c + (0, self.offset, 0)
        z1 = c - (0, 0, self.offset)
        z2 = c + (0, 0, self.offset)
        l = len(c)
        l2 = 2*l
        l3 = 3*l
        l4 = 4*l
        l5 = 5*l
        l6 = 6*l
        f = numpy.zeros((l*3, 2),'int32')
        f[:l,0] = xrange(l)
        f[:l,1] = xrange(l, l2)
        f[l:l2,0] = xrange(l2, l3)
        f[l:l2,1] = xrange(l3, l4)
        f[l2:l3,0] = xrange(l4, l5)
        f[l2:l3,1] = xrange(l5, l6)
        return numpy.concatenate( (x1,x2,y1,y2,z1,z2) ), f


    def getCoords(self, c):
        c = numpy.array(c)
        x1 = c - (self.offset, 0, 0)
        x2 = c + (self.offset, 0, 0)
        y1 = c - (0, self.offset, 0)
        y2 = c + (0, self.offset, 0)
        z1 = c - (0, 0, self.offset)
        z2 = c + (0, 0, self.offset)
        return numpy.concatenate( (x1,x2,y1,y2,z1,z2) )
    

    def getFaces(self, length=None):
        if length is None:
            length = len(self.vertexSet.vertices.array)
        off = []
        for i in range(1,6): off.append(i*length)
        faces = []
        for i in xrange(length):
            faces.append( (i, i+off[0]) )
            faces.append( (i+off[1], i+off[2]) )
            faces.append( (i+off[3], i+off[4]) )
        return faces


    def Draw(self):
        from time import time
        t0 = time()
        #centers = self.getCoords(self.vertexSet.vertices.array)
        #if len(centers)==0:
        #    return
        #    #faces = numpy.zeros((0,3), 'i')
        #else: faces = self.getFaces(len(self.vertexSet.vertices.array))

        centers, faces = self.getCoordsAndFaces(self.vertexSet.vertices.array)
        if len(centers)==0:
            return

        if self.materials[GL.GL_FRONT] and \
           not self.inheritMaterial:
            fpProp = self.materials[GL.GL_FRONT].prop[:5]
            fpBind = self.materials[GL.GL_FRONT].binding[:5]
        else:
            fpProp = None
            fpBind = None

        if self.materials[GL.GL_BACK] and \
           not self.inheritMaterial:
            bpProp = self.materials[GL.GL_BACK].prop[:5]
            bpBind = self.materials[GL.GL_BACK].binding[:5]
        else:
            bpProp = None
            bpBind = None

        texCoords = None

        if self.lighting:
            if (self.invertNormals) and (self.normals is not None):
                norms = - self.normals
            else:
                norms = self.normals
        else:
            norms = None

        #GL.glDisable(GL.GL_DEPTH_TEST)

        if self.disableStencil is True:
            GL.glDisable(GL.GL_STENCIL_TEST)

        status = glDrawIndexedGeom(
            GL.GL_LINES,
            centers.astype('f'),
            faces,
            norms,
            texCoords,
            fpProp, bpProp, fpBind, bpBind,
            self.frontAndBack, 1)        

        if self.disableStencil is True:
            GL.glEnable(GL.GL_STENCIL_TEST)

        #print 'time for faces', time()-t0

        #GL.glEnable(GL.GL_DEPTH_TEST)

        return status
