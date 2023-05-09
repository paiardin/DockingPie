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

########################################################################
#
# Date: 2000 Authors: Michel F. SANNER, Daniel Stoffler, Guillaume Vareille
#
#    vareille@scripps.edu
#    sanner@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Michel F. SANNER, Daniel Stoffler, Guillaume Vareille and TSRI
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/DejaVu2/IndexedGeom.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#
# $Id: IndexedGeom.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#


import numpy, types
import warnings

from opengltk.OpenGL import GL
from opengltk.extent import _gllib
from opengltk.extent.utillib import glDrawIndexedGeom
from geomutils.geomalgorithms import  TriangleNormals

import DejaVu2
from DejaVu2.datamodel import FaceSet
from DejaVu2.viewerFns import checkKeywords
from DejaVu2 import viewerConst
from DejaVu2.Geom import Geom


class IndexedGeom(Geom):
    """Geometry specified by a VertexSet and a FaceSet
"""
    keywords = Geom.keywords + [
        'type',
        'faces',
        'fnormals',
        'freshape',
        ]


    def __init__(self, name=None, check=1, **kw):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        #self.outlList = GL.glGenLists(1)

        self.numOfHighlightedIndex = 0

        if not kw.get('shape'):
            kw['shape'] = (0,3)    # default shape for sphere set

        self.faceSet = FaceSet( shape= (0,0) )
        self.pickingFaces = None
        apply( Geom.__init__, (self, name, check), kw)

        self._modified = False

        
    def getState(self, full=False):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        state = Geom.getState(self, full)
        if full:
            state['faces'] = self.getFaces()
            state['fnormals'] = self.getFNormals()
            
        return state


    def getFaces(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """returns a handle to the faces array"""
        if len(self.faceSet.faces.array)== 0:
            return []
        return self.faceSet.faces.array


    def getFNormals(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """returns a handle to the face normals"""
        if self.faceSet.normals.status == viewerConst.NONE:
            self.faceSet.normals.GetProperty()

        if len(self.faceSet.normals.array) == 0:
            return []
        return self.faceSet.normals.array


    def splitFacesOnHighlightedVertices(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        #print "splitFacesOnHighlightedVertices", self.name
        if len(self.highlight) > 0:
            lFacesWithHighlightedVertices = []
            lFacesWithoutHighlightedVertices = []
            for face in self.faceSet.faces.array:
               for lVertexIndex in face:
                   if self.highlight[lVertexIndex]:
                       lFacesWithHighlightedVertices.append(face)
                       break
               else: # we didn't break
                   lFacesWithoutHighlightedVertices.append(face)
            lFaces = lFacesWithHighlightedVertices + lFacesWithoutHighlightedVertices
            self.numOfHighlightedIndex = len(lFacesWithHighlightedVertices)
            self.faceSet.faces.SetValues(lFaces)
        else:
            self.numOfHighlightedIndex = 0


    def removeFacesWithoutHighlightedVertices(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        if self.numOfHighlightedIndex > 0:
            self.Set(faces=self.faceSet.faces.array[:self.numOfHighlightedIndex])

#    def removeFacesWithoutHighlightedVertices(self):
#        if len(self.highlight) > 0:
#            lFacesWithHighlightedVertices = []
#            for face in self.faceSet.faces.array:
#               for lVertexIndex in face:
#                   if self.highlight[lVertexIndex]:
#                       lFacesWithHighlightedVertices.append(face)
#                       break
#            self.Set(faces=lFacesWithHighlightedVertices)


    def _FixedLengthFaces(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """sets self.fixedLength to the number of vertices perface if all faces
        have the same number of vertices, else self.fixedLength=0.
        Check if there are negative indices finishing faces lists"""

        ind = self.faceSet.faces.array
        min = numpy.minimum.reduce( numpy.minimum.reduce (ind) )
        if min > -1 and ind.shape[1] < 5: self.fixedLength = ind.shape[1]
        else: self.fixedLength = False

    def getVisibleVertices(self, picking=False):
        # return a list of vertices currently used to daw this geoemetry
        # the list is used for picking
        if picking and self.pickingVertices is not None:
            indices = numpy.unique(self.pickingFaces.flatten())
            return self.pickingVertices[indices], indices
        else:
            indices = numpy.unique(self.faceSet.faces.array.flatten())
            return self.vertexSet.vertices.array[indices], indices


    def _PrimitiveType(self, type=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Set the geometric primitives type for indexed geometries
        Type can be: None, GL_LINES, GL_LINE_STRIP, GL_LINE_LOOP
                           GL_TRIANGLES, GL_QUADS, GL_POLYGON, GL_TRIANGLE_FAN
"""
        #print "IndexedGeom._PrimitiveType", self, type
        #
        # - GL_POINTS, GL_TRIANGLE_STRIP and GL_QUAD_STRIP are not
        #   considered because they are not indexed geometries
        #
        # - GL_LINES, GL_TRIANGLES, GL_QUADS, are NOT pickable but fast
        #
        # - GL_LINE_STRIP, GL_LINE_LOOPS, GL_POLYGON, GL_TRIANGLE_FAN
        #   return per primitive picking info
        #
        assert type in viewerConst.PRIMITIVES+(None,)
        if len(self.faceSet)==0: return
        self._FixedLengthFaces()
        old = self.primitiveType

        # no type has been given
        # so use the most efficient primitive
        if type is None:
            if not self.pickableVertices:
                if self.fixedLength==2:
                    self.primitiveType = GL.GL_LINES
                elif self.fixedLength==3:
                    self.primitiveType = GL.GL_TRIANGLES
                elif self.fixedLength==4:
                    self.primitiveType = GL.GL_QUADS
                else:
                    self.primitiveType = GL.GL_POLYGON
                    self.pickableVertices = True # it will pickable
            else:
                if self.fixedLength==2:
                    # MS DEC 02: we make it a line strip so the display list
                    # build by the *DSPL function will let pick parts
                    self.primitiveType = GL.GL_LINE_STRIP
                elif self.fixedLength==3:
                    self.primitiveType = GL.GL_TRIANGLES
                elif self.fixedLength==4:
                    self.primitiveType = GL.GL_QUADS
                else:
                    self.primitiveType = GL.GL_POLYGON
                
        else: # type has been provided
            if type == GL.GL_LINES:
                if self.fixedLength==2: 
                    self.primitiveType = GL.GL_LINES
                    self.pickableVertices = False
                else:
                    raise AttributeError('Bad faces for GL.GL_LINES')

            elif type == GL.GL_TRIANGLES:
                if self.fixedLength==3: 
                    self.primitiveType = GL.GL_TRIANGLES
                    self.pickableVertices = False
                else:
                    raise AttributeError('Bad faces for GL.GL_TRIANGLES')

            elif type == GL.GL_QUADS:
                if self.fixedLength==4:
                    self.primitiveType = GL.GL_QUADS
                    self.pickableVertices = False
                else: raise AttributeError('Bad faces for GL.GL_QUADS')

            elif type == GL.GL_TRIANGLE_FAN:
                self.primitiveType = GL.GL_QUADS
                self.pickableVertices = True

            else:
                self.primitiveType = type
                self.pickableVertices = True

        if old != self.primitiveType:
            self.redoDspLst = 1


    def Add(self, check=1, redo=1, **kw):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """add faces (polygon or lines) to this object
"""
        #print "IndexedGeom.Add"

        if __debug__:
            if check:
                apply( checkKeywords, (self.name,self.keywords), kw)
            
        t = kw.get( 'type')
        f = kw.get( 'faces')
        fn = kw.get( 'fnormals')
        if f is not None and len(f):
            self.redoDspLst = 1
            self.faceSet.faces.AddValues( f )

        if fn is not None and len(fn):
            self.redoDspLst = 1
            self.faceSet.faces.AddValues(fn)

        Geom.Add(self, check=0, redo=0,
                 vertices = kw.get( 'vertices'),
                 vnormals = kw.get( 'vnormals'),
                 materials = kw.get( 'materials'),
                 polyFace = kw.get( 'polyFace'),
                 matBind = kw.get( 'matBind'),
                 propName = kw.get( 'propName') )

        if f is not None:
            pf = kw.get( 'polyFace')
            pn = kw.get( 'propName')
            mbm = kw.get( 'matBind')
            self._PrimitiveType(t)
            self.MaterialBindingMode(pn, face=pf, mode=mbm)
            
        if f is not None or fn is not None:
            if self.shading==GL.GL_FLAT:
                self.GetNormals()

        if self.viewer and redo:
            if self.redoDspLst:
                self.viewer.objectsNeedingRedo[self] = None
#                self.RedoDisplayList()


    def Set(self, check=1, redo=1, updateOwnGui=True, **kw):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """set data for this object: add faces (polygon or lines) to this object
check=1 : verify that all the keywords present can be handle by this func 
redo=1 : append self to viewer.objectsNeedingRedo
updateOwnGui=True : allow to update owngui at the end this func
"""
        #print "IndexedPolygons.Set"

        #import pdb; pdb.set_trace()

        redoFlags = 0

        # Exceptionnaly this has to be before the call to Geom.Set
        # because we want to override the treatment of it by Geom.Set
        invertNormals = kw.get('invertNormals')
        if invertNormals is not None:
            kw.pop('invertNormals')
            if self.invertNormals != invertNormals:
                self.invertNormals = invertNormals
                redoFlags |= self._redoFlags['redoDisplayListFlag']

                if hasattr(self, 'vertexArrayFlag') \
                  and self.vertexArrayFlag is True \
                  and len(self.vertexSet.normals.array) > 0:
                   if self.invertNormals:
                       lnormals = -self.vertexSet.normals.array
                   else:
                       lnormals = self.vertexSet.normals.array
                   if hasattr(DejaVu2, 'enableVBO') and DejaVu2.enableVBO is True:
                       from opengltk.extent import _glextlib 
                       _glextlib.glBindBufferARB(_glextlib.GL_ARRAY_BUFFER_ARB,
                                              int(self.vertexArrayFlagBufferList[1]))
                       _glextlib.glBufferDataARB(_glextlib.GL_ARRAY_BUFFER_ARB,
                             4*len(self.vertexSet.normals.array)*len(self.vertexSet.normals.array[0]),
                                              lnormals,
                                              _glextlib.GL_STATIC_DRAW_ARB)
                       _gllib.glNormalPointer(GL.GL_FLOAT, 0, 0)
                   elif DejaVu2.enableVertexArray is True:
                       _gllib.glNormalPointer(GL.GL_FLOAT, 0, lnormals)

        # Exceptionnaly this has to be before the call to Geom.Set
        # because we want to complete the treatment of it by Geom.Set
        highlight = kw.get('highlight', None)

        redoFlags |= apply( Geom.Set, (self, check, 0), kw)

        if kw.has_key('faces'):
            self.numOfHighlightedIndex = 0
            self.faceSet.faces.SetValues( [] )
            if hasattr(self, 'edges'): del self.edges
            if hasattr(self, 'faceEdges'): del self.faceEdges

        t = kw.get( 'type')
        f = kw.get( 'faces')
        reshape = kw.get( 'freshape')
        fn = kw.get( 'fnormals')

        if not f is None:
            try:
                len(f)
            except TypeError:
                raise TypeError ("faces should be sequences of integers")

            if len(f)==1 and len(f[0])==0:  # handle [[]]
                f = []

            ind = numpy.array(f)
            if len(ind.ravel()) > 0:
                m = numpy.minimum.reduce(ind)
                if m.size > 1:
                    m = min(m)
                m = max(0, m)
                if ( m < 0 ):
                    raise ValueError ("vertex index %d out of range" % m)

                m = numpy.maximum.reduce(ind)
                if m.size > 1:
                    m = max(m)
                if ( m >= len(self.vertexSet) ):
                    raise ValueError ("vertex index %d out of range, max %d" %
                                      (m, len(self.vertexSet)-1) )

            redoFlags |= self._redoFlags['redoDisplayListFlag']
            self.faceSet.faces.SetValues( f, reshape)
            assert len(self.faceSet.faces.ashape)==2

        if highlight is not None: #len(self.highlight) > 0:
            # the rest of the treatment was done in Geom.Set
            self.splitFacesOnHighlightedVertices()

        if not fn is None:
            redoFlags |= self._redoFlags['redoDisplayListFlag']
            self.faceSet.normals.SetValues(fn)

        if not f is None or t or kw.get( 'pickableVertices'):
            pf = kw.get( 'polyFace')
            if type(pf) is types.StringType:
                pf = self._strToOGL[pf.lower()]
            pn = kw.get( 'propName')
            mbm = kw.get( 'matBind')
            self._PrimitiveType(t)
            self.MaterialBindingMode(pn, face=pf, mode=mbm)

        if f is not None or fn is not None:
            if self.shading==GL.GL_FLAT:
                redoFlags |= self._redoFlags['redoDisplayListFlag']
                self.faceSet.normals.PropertyStatus(len(self.faceSet))
                if self.lighting:
                    self.GetNormals()

        if self.faceSet.normals.status < 24:
            if self.lighting:
                self.GetNormals()

        if (f is not None or highlight is not None) and \
           hasattr(DejaVu2, 'enableVBO') and \
           DejaVu2.enableVBO is True and \
           self.vertexArrayFlag is True and \
           hasattr(self, 'faceSet') and \
           len(self.faceSet.faces.array) > 0:

            # faces
            #print "contiguous", self.faceSet.faces.array.flags.contiguous
            #if self.faceSet.faces.array.flags.contiguous is False:
            #    self.faceSet.faces.array = numpy.array(self.faceSet.faces.array,copy=1)
            from opengltk.extent import _glextlib
            _glextlib.glBindBufferARB(_glextlib.GL_ELEMENT_ARRAY_BUFFER,
                                      int(self.vertexArrayFlagBufferList[2]))
            _glextlib.glBufferDataARB(_glextlib.GL_ELEMENT_ARRAY_BUFFER,
                  4*len(self.faceSet.faces.array)*len(self.faceSet.faces.array[0]),
                            self.faceSet.faces.array, _glextlib.GL_STATIC_DRAW_ARB)
            _glextlib.glBindBufferARB(_glextlib.GL_ELEMENT_ARRAY_BUFFER, 0)

        return self.redoNow(redo, updateOwnGui, redoFlags)


    def ComputeVertexNormals(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Compute the vertex normals"""
        v = self.vertexSet.vertices.array
        f = self.faceSet.faces.array
        if len(v) > 2 and len(f) > 1:
            return TriangleNormals( v, f[:,:3], 'PER_VERTEX')
        else: return None


    def ComputeFaceNormals(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Compute the face normals"""

        v = self.vertexSet.vertices.array
        f = self.faceSet.faces.array
        if len(v) > 2 and len(f) > 0:
            return TriangleNormals( v, f[:,:3], 'PER_FACE')
        else: return None


    def VertexNormalFunction(self, func=None, args=()):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Set the function used to compute vertices normals"""
        if func is None: return self.vertexSet.normals.Compute
        assert callable(func)
        self.vertexSet.normals.ComputeFunction( func, args )


    def FaceNormalFunction(self, func=None, args=()):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Set the function used to compute faces normal"""

        if func is None: return self.faceSet.normals.Compute
        assert callable(func)
        self.faceSet.normals.ComputeFunction( func, args )

            
    def DisplayFunction(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """display a set of indexed geometric primitives"""
        
        if self.dpyList:

#            print "DisplayFunction", self.dpyList, self.fullName

            lDrawOutline = (self.getDrawOutlineMode('front'), self.getDrawOutlineMode('back'))
            if (lDrawOutline[0] or lDrawOutline[1]) and self.viewer.hasOffsetExt:

                outl = self.outline

                if   self.GetPolyMode('front') == GL.GL_FILL \
                  or self.GetPolyMode('back') == GL.GL_FILL:

                    mode = GL.GL_POLYGON_OFFSET_FILL

                    GL.glEnable(mode)
                    self.viewer.polyOffset( outl.factor, outl.unit)
                    Geom.DisplayFunction(self)
                    GL.glDisable(mode)

                    GL.glPolygonMode(GL.GL_FRONT_AND_BACK, GL.GL_LINE)
                    if not outl.colorAsMaterial:
                        if outl.lighting:
                            GL.glMaterialfv( GL.GL_FRONT_AND_BACK,
                                             GL.GL_EMISSION,
                                             outl.color )
                        else:
                            GL.glDisable(GL.GL_LIGHTING)
                            GL.glColor4fv (outl.color)

                    GL.glLineWidth(outl.lineWidth)

                    if lDrawOutline[0] is False or lDrawOutline[1] is False:
                        GL.glEnable(GL.GL_CULL_FACE)
                        if lDrawOutline[0]:
                            GL.glCullFace(GL.GL_BACK)
                        elif lDrawOutline[1]:
                            GL.glCullFace(GL.GL_FRONT)
                    else:
                        GL.glDisable(GL.GL_CULL_FACE)

                    if outl.dpyList:
                        currentcontext = self.viewer.currentCamera.getContext()
                        if currentcontext != outl.dpyList[1]:
                            warnings.warn("""DisplayFunction failed because the current context is the wrong one""")
                            #print "currentcontext != outl.dpyList[1]", currentcontext, outl.dpyList[1]
                        else:
                            #print '#%d'%outl.dpyList[0], currentcontext, "glCallList IndexedGeom"
                            GL.glCallList(outl.dpyList[0])

                    GL.glEnable(GL.GL_CULL_FACE)
                    GL.glEnable(GL.GL_LIGHTING)

                else:
                    Geom.DisplayFunction(self)
            else:
                Geom.DisplayFunction(self)


    def Draw(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """ draw geom
"""
        #print "IndexedGeom.Draw", self.name
        if self.vertexArrayFlag is True \
          and DejaVu2.enableVertexArray is True:
            if not (hasattr(self.vertexSet, "texCoords") \
              and self.vertexSet.texCoords.status >= viewerConst.COMPUTED):
                return self.drawVertexArray()

        if len(self.faceSet) and len(self.vertexSet):
            if self.materials[GL.GL_FRONT] and \
                   not self.inheritMaterial:
                mat = self.materials[GL.GL_FRONT]
                fpProp = []
                fpBind = []
                for propInd in range(4):
                    b, p = mat.GetProperty(propInd)
                    fpProp.append(p)
                    fpBind.append(b)
                fpProp.append(mat.prop[4])
                fpBind.append(mat.binding[4])
            else:
                fpProp = None
                fpBind = None

            if self.materials[GL.GL_BACK] and \
               not self.inheritMaterial:
                mat = self.materials[GL.GL_BACK]
                bpProp = []
                bpBind = []
                for propInd in range(4):
                    b, p = mat.GetProperty(propInd)
                    bpProp.append(p)
                    bpBind.append(b)
                bpProp.append(mat.prop[4])
                bpBind.append(mat.binding[4])

            else:
                bpProp = None
                bpBind = None

            texCoords = None
            if hasattr(self.vertexSet, "texCoords"):
                if self.vertexSet.texCoords.status >= viewerConst.COMPUTED:
                    texCoords = self.vertexSet.texCoords.array

            if self.lighting and self.normals is not None:
                if self.invertNormals:
                    norms = - self.normals
                else:
                    norms = self.normals
            else:
                norms = None

            from DejaVu2 import preventIntelBug_BlackTriangles
            if preventIntelBug_BlackTriangles:
                preventIntelBug = 1
            else:
                preventIntelBug = 0

            lsharpColorBoundaries = self.getSharpColorBoundaries()

            if self.disableStencil is True:
                GL.glDisable(GL.GL_STENCIL_TEST)

            status = glDrawIndexedGeom(
                self.primitiveType,
                self.vertexSet.vertices.array,
                self.faceSet.faces.array,
                norms,
                texCoords,
                fpProp, bpProp, fpBind, bpBind,
                self.frontAndBack, 1,
                lsharpColorBoundaries,
                preventIntelBug,
                highlight=self.highlight,
                )

            if self.disableStencil is True:
                GL.glEnable(GL.GL_STENCIL_TEST)

            return status


    def drawVertexArray(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """ drawVertexArray
"""
        #print "drawVertexArray", self.name
        if hasattr(self, 'faceSet') and len(self.faceSet.faces.array) > 0:
            # vertices
            GL.glEnableClientState(GL.GL_VERTEX_ARRAY)
            
            # normals
            if len(self.vertexSet.normals.array) > 0:
                GL.glEnableClientState(GL.GL_NORMAL_ARRAY)
        
            # colors
            if hasattr(self, 'colorPointerIsOn') and self.colorPointerIsOn is True:
                GL.glEnableClientState(GL.GL_COLOR_ARRAY)
                from DejaVu2 import preventIntelBug_WhiteTriangles
                if preventIntelBug_WhiteTriangles:
                    GL.glColorMaterial(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE)
                else:
                    GL.glColorMaterial(GL.GL_FRONT_AND_BACK, GL.GL_AMBIENT_AND_DIFFUSE)
                GL.glEnable( GL.GL_COLOR_MATERIAL )
                
            if DejaVu2.enableVBO is True:
                from opengltk.extent import _glextlib
                # vertices
                _glextlib.glBindBufferARB(_glextlib.GL_ARRAY_BUFFER_ARB,
                                 int(self.vertexArrayFlagBufferList[0]))
                _gllib.glVertexPointer(len(self.vertexSet.vertices.array[0]),
                                              GL.GL_FLOAT, 0, 0)
                # normals
                if len(self.vertexSet.normals.array) > 0:
                    _glextlib.glBindBufferARB(_glextlib.GL_ARRAY_BUFFER_ARB,
                                              int(self.vertexArrayFlagBufferList[1]))
                    _gllib.glNormalPointer(GL.GL_FLOAT, 0, 0)
                # colors
                if hasattr(self, 'colorPointerIsOn') and self.colorPointerIsOn is True:
                    _glextlib.glBindBufferARB(_glextlib.GL_ARRAY_BUFFER_ARB,
                                              int(self.vertexArrayFlagBufferList[3]))
                    _gllib.glColorPointer(4, GL.GL_FLOAT, 0, 0)
            else:
                # vertices
                _gllib.glVertexPointer(len(self.vertexSet.vertices.array[0]),
                                   GL.GL_FLOAT, 0, self.vertexSet.vertices.array)
                # normals
                if len(self.vertexSet.normals.array) > 0:
                    _gllib.glNormalPointer(GL.GL_FLOAT, 0, self.vertexSet.normals.array)
                # colors
                if hasattr(self, 'colorPointerIsOn') and self.colorPointerIsOn is True:
                    _gllib.glColorPointer(4, GL.GL_FLOAT, 0, self.colorArray)

            # Draw faces
            if self.primitiveType == GL.GL_LINE_STRIP:
                lPrimitiveType = GL.GL_LINES
            elif self.primitiveType == GL.GL_TRIANGLES:
                #print "triangles' length:", len(self.faceSet.faces.array[0])
                lPrimitiveType = GL.GL_TRIANGLES
            elif self.primitiveType == GL.GL_QUADS:
                #print "quads' length:", len(self.faceSet.faces.array[0])
                lPrimitiveType = GL.GL_QUADS
            else:
                #print "what's that ?" , self.primitiveType
                lPrimitiveType = self.primitiveType

            lNumOfNonHighlightedIndices = len(self.faceSet.faces.array) - self.numOfHighlightedIndex
            if DejaVu2.enableVBO is True:
                _glextlib.glBindBufferARB(_glextlib.GL_ELEMENT_ARRAY_BUFFER, 
                                      int(self.vertexArrayFlagBufferList[2])) #this protect from unexplained segfault
                if self.disableStencil is False and self.numOfHighlightedIndex > 0:
                    # highlighted
                    GL.glStencilFunc(GL.GL_ALWAYS, 1, 1)
                    _gllib.glDrawElements(lPrimitiveType,
                                  self.numOfHighlightedIndex*len(self.faceSet.faces.array[0]),
                                  GL.GL_UNSIGNED_INT,
                                  0 
                                 )
                    GL.glStencilFunc(GL.GL_ALWAYS, 0, 1)

                    # non highlighted
                    _gllib.glDrawElements(lPrimitiveType,
                                  lNumOfNonHighlightedIndices * len(self.faceSet.faces.array[0]),
                                  GL.GL_UNSIGNED_INT,
                                  self.numOfHighlightedIndex * len(self.faceSet.faces.array[0]) * 4
                                 )
                else:
                    _gllib.glDrawElements(lPrimitiveType,
                                  len(self.faceSet.faces.array)*len(self.faceSet.faces.array[0]),
                                  GL.GL_UNSIGNED_INT,
                                  0 
                                 )
                _glextlib.glBindBufferARB(_glextlib.GL_ELEMENT_ARRAY_BUFFER, 0 ) #this protect from unexplained segfault
            else:
                if self.disableStencil is False and self.numOfHighlightedIndex > 0:
                    # highlighted
                    GL.glStencilFunc(GL.GL_ALWAYS, 1, 1)
                    _gllib.glDrawElements(lPrimitiveType,
                                  self.numOfHighlightedIndex*len(self.faceSet.faces.array[0]),
                                  GL.GL_UNSIGNED_INT,
                                  self.faceSet.faces.array 
                                 )
                    GL.glStencilFunc(GL.GL_ALWAYS, 0, 1)

                    # non highlighted
                    _gllib.glDrawElements(lPrimitiveType,
                                  lNumOfNonHighlightedIndices * len(self.faceSet.faces.array[0]),
                                  GL.GL_UNSIGNED_INT,
                                  self.faceSet.faces.array[self.numOfHighlightedIndex: ]
                                 )
                else:
                    _gllib.glDrawElements(lPrimitiveType,
                                  len(self.faceSet.faces.array)*len(self.faceSet.faces.array[0]),
                                  GL.GL_UNSIGNED_INT,
                                  self.faceSet.faces.array
                                 )
            if hasattr(self, 'colorPointerIsOn') and self.colorPointerIsOn is True:
                GL.glDisable( GL.GL_COLOR_MATERIAL )
                GL.glDisableClientState(GL.GL_COLOR_ARRAY)
            if len(self.vertexSet.normals.array) > 0:
                GL.glDisableClientState(GL.GL_NORMAL_ARRAY)
            GL.glDisableClientState(GL.GL_VERTEX_ARRAY)

        return 1


    def RedoDisplayList(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
            #print "IndexedGeom.RedoDisplayList", self.name
##          if __debug__:
##              print 'IndexedGeom RedoDisplayList for', self.fullName

        Geom.RedoDisplayList(self)

        if len(self.faceSet) and len(self.vertexSet) \
          and (   self.primitiveType == GL.GL_TRIANGLES \
               or self.primitiveType == GL.GL_QUADS \
               or self.primitiveType == GL.GL_POLYGON):

            # we always build this, that way we don't have to built on demand
            outl = self.outline
            if outl.colorAsMaterial:
                if self.materials[GL.GL_FRONT] and \
                       not self.inheritMaterial:
                    mat = self.materials[GL.GL_FRONT]
                    fpProp = []
                    fpBind = []
                    for propInd in range(4):
                        b, p = mat.GetProperty(propInd)
                        fpProp.append(p)
                        fpBind.append(b)
                    fpProp.append(mat.prop[4])
                    fpBind.append(mat.binding[4])
                else:
                    fpProp = None
                    fpBind = None

                if self.materials[GL.GL_BACK] and \
                   not self.inheritMaterial:
                    mat = self.materials[GL.GL_BACK]
                    bpProp = []
                    bpBind = []
                    for propInd in range(4):
                        b, p = mat.GetProperty(propInd)
                        bpProp.append(p)
                        bpBind.append(b)
                    bpProp.append(mat.prop[4])
                    bpBind.append(mat.binding[4])

                else:
                    bpProp = None
                    bpBind = None
            else:
                fpProp = bpProp = fpBind = bpBind = None

            texCoords = None

            if outl.lighting:
                norms = self.normals
            else:
                norms = None

            # WARNING: if texture, fpProp, bpProp, fpBind, bpBind,
            # are not passed (either None or arrays) we get a segmentation
            # fault if the surface has many colors (i.e. color MSMS by atom
            # type and dispaly outline seg faults)

            # calling with too many arguments segaults too
            # just add  None, None, None, None, after the line with colors

            if hasattr(outl, 'dpyList') and outl.dpyList is not None:
                lNewList = outl.dpyList[0]
            else:
                lNewList = GL.glGenLists(1)
                self.viewer.deleteOpenglList()

            #print "lNewList IndexedGeom.RedoDisplayList", lNewList, self.name
            lCurrentContext = self.viewer.currentCamera.getContext()
            outl.dpyList = ( lNewList, lCurrentContext)
                             
            GL.glNewList(outl.dpyList[0], GL.GL_COMPILE)
            #print '+%d'%outl.dpyList[0], lCurrentContext, "glNewList IndexedGeom"
            status=glDrawIndexedGeom(
                GL.GL_TRIANGLES,
                self.vertexSet.vertices.array,
                self.faceSet.faces.array,
                norms,
                texCoords,
                fpProp, bpProp, fpBind, bpBind,
                self.frontAndBack,
                1)  # 1 means use diffuse component if no lighting
            #print '*%d'%GL.glGetIntegerv(GL.GL_LIST_INDEX), "glEndList IndexedGeom"
            GL.glEndList()
            if not status:
                #print '-%d'%outl.dpyList[0], "glDeleteLists IndexedGeom"
                GL.glDeleteLists(outl.dpyList[0], 1)
                outl.dpyList = None


    def _setTransparent(self, val):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        #print "IndexedGeom._setTransparent", self.name, val
        t1 = self.materials[GL.GL_FRONT].fixOpacity()
        t2 = self.materials[GL.GL_BACK].fixOpacity()
        lDetectedTransparency =  t1 or t2

        if val == 'implicit':
            val = lDetectedTransparency \
                   or ( (self.getDrawOutlineMode('front') is True or \
                         self.getDrawOutlineMode('back') is True)
                        and self.outline.color[3]<1.0)
        if val is True:
            val = 1
        elif val is False:
            val = 0
        assert val in [0,1], "only 0 or 1 are possible"
        self.transparent = val

        if self.viewer:
            if val in (1, True):
                for c in self.viewer.cameras:
                    c.addButtonUpCB(self.sortPoly_cb)
            else:
                for c in self.viewer.cameras:
                    if self.sortPoly_cb in c.onButtonUpCBlist:
                        c.delButtonUpCB(self.sortPoly_cb)

        return self._redoFlags['redoDisplayListFlag'] | self._redoFlags['redoViewerDisplayListFlag'] # if container we force rebuilding main dpyList


    def buildEdgeList(self):
        """
        builds and returns a list of edges at pairs of vertex indices
        also return a dict with key face and value list of edge indices
        """
        edges = []
        edgeKey = {}
        en = 0
        faces = self.getFaces()
        for f in faces:
            s1,s2,s3 = f
            if s1<s2:
                edges.append( (s1,s2) )
                edgeKey['%d,%d'%(s1,s2)] = en
                en += 1
            if s2<s3:
                edges.append( (s2,s3) )
                edgeKey['%d,%d'%(s2,s3)] = en
                en += 1
            if s3<s1:
                edges.append( (s3,s1) )
                edgeKey['%d,%d'%(s3,s1)] = en
                en += 1

        faceEdges = map( lambda x: [], faces)
        for fn, f in enumerate(faces):
            s1,s2,s3 = f
            key = '%d,%d'%(s1,s2)
            if edgeKey.has_key(key):
                faceEdges[fn].append(edgeKey[key])
            else:
                key = '%d,%d'%(s2,s1)
                faceEdges[fn].append(edgeKey[key])

            key = '%d,%d'%(s2,s3)
            if edgeKey.has_key(key):
                faceEdges[fn].append(edgeKey[key])
            else:
                key = '%d,%d'%(s3,s2)
                faceEdges[fn].append(edgeKey[key])

            key = '%d,%d'%(s3,s1)
            if edgeKey.has_key(key):
                faceEdges[fn].append(edgeKey[key])
            else:
                key = '%d,%d'%(s1,s3)
                faceEdges[fn].append(edgeKey[key])

        return edges, faceEdges


    def getEdges(self):
        if not hasattr(self, 'edges') or not hasattr(self, 'faceEdges'):
           edges, faceEdges = self.buildEdgeList()
           self.edges = edges
           self.faceEdges = faceEdges

        return self.edges, self.faceEdges

    
    def removeDuplicatedVertices(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """find duplicated vertices and remove them, re-index face list"""
        # hash vertices
        d = {}
        for vert in self.vertexSet.vertices.array:
            d['%f,%f,%f'%tuple(vert)] = []

        # build list of unique vertices and lookup table
        lookup = {}
        nv = []
        nn = []
        i = 0
        for k in d.keys():
            nv.append(eval(k))
            lookup[k] = i
            i = i + 1

        # new facelist
        v = self.vertexSet.vertices.array
        nflist = []
        for face in self.faceSet.faces.array:
            nf = []
            for vind in face:
                nf.append(lookup['%f,%f,%f'%tuple(v[vind])])
            nflist.append(nf)

        return nv, nflist


#    def modifiedFacesAndVerticesForSharpColorBoundaries(self):
##        self.vertexSet.vertices.array,
##        self.faceSet.faces.array,
##        self.normals,
##        fpProp, bpProp, fpBind, bpBind,
##        self.frontAndBack,
#
#        if not \
#           (    ( self.faceSet.faces.array.shape[1] == 3 ) \
#            and ( fpBind ) \
#            and (   ( fpBind [ 0 ] == PER_VERTEX ) \
#                 or ( fpBind [ 1 ] == PER_VERTEX ) \
#                 or ( fpBind [ 2 ] == PER_VERTEX ) \
#                 or ( fpBind [ 3 ] == PER_VERTEX ) \
#                 or ( fpBind [ 4 ] == PER_VERTEX ) \
#                 or (    ( not self.frontAndBack )
#                     and ( backMatBind ) 
#                     and (   ( backMatBind [ 0 ] == PER_VERTEX )
#                          or ( backMatBind [ 1 ] == PER_VERTEX )
#                          or ( backMatBind [ 2 ] == PER_VERTEX )
#                          or ( backMatBind [ 3 ] == PER_VERTEX )
#                          or ( backMatBind [ 4 ] == PER_VERTEX ) 
#                         )
#                    )
#                )
#           ):
#            return
#
#        if self.normals:
#  {
#                if (lennorm[0] == lencoord) 
#                        normBinding = PER_VERTEX;
#    else if (lennorm[0] == lenind[0]) 
#                        normBinding = PER_PART;
#    else if (lennorm[0] == 1) 
#            normBinding = OVERALL;
#    else normBinding = NONE;
#  }
#  else normBinding = NONE;
#
#

