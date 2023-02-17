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
# Date: 2000 Authors: Michel Sanner, Kevin Chan
#
#    vareille@scripps.edu
#    sanner@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Michel Sanner and TSRI 2000
#
# revision: Guillaume Vareille
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/DejaVu2/triangle_strip.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#
# $Id: triangle_strip.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#


from opengltk.OpenGL import GL
from opengltk.extent import _gllib as gllib

from DejaVu2.Geom import Geom
from DejaVu2.IndexedPolygons import IndexedPolygons
from DejaVu2.viewerFns import checkKeywords
import DejaVu2.datamodel, DejaVu2.viewerConst
import numpy
from DejaVu2.colorTool import glMaterialWithCheck, resetMaterialMemory

class Triangle_strip(IndexedPolygons):
    """ Class to draw a triangle strip or multiple triangle strip 
    geometries. """

    keywords = IndexedPolygons.keywords + [
        'stripBegin', 'stripEnd',
        'fnormals',
        ]


    def __init__(self, name=None, check=1, **kw):
        """ Constructor:
Takes an array of vertices and splits it into separate triangle strips
at specific vertices specified by stripBegin, which is and array of
indices. Generates faceSet. In this class, a face corresponds to a
triangle strip.  Calls IndexedPolygons constructor.  Calls
MatBindingMode() to set the correct mode specific to a triangle strip,
overwriting the binding mode set by MaterialBindingMode in the Geom
class. Initializes frontPolyMode to FILL.  Initializes currentCol and
currentMat."""

        self.stripBegin = None
        self.stripEnd = None
        self.normals = None

        apply( IndexedPolygons.__init__, (self, name, 0), kw )

        for face in self.materials.keys():
            for j in range(5):
                self.MatBindingMode(j, face)

        self.frontPolyMode = GL.GL_FILL
        #self.inheritFrontPolyMode = 0

        self.currentCol = []
        for i in range(4):
            self.currentCol.append(-1.0)
        self.currentMat = []
        for i in range(2):
            mat = []
            for j in range(5):
                rgba = []
                for k in range(4):
                    rgba.append(-1.0)
                mat.append(rgba)
            self.currentMat.append(mat)


    def Set(self, check=1, redo=1, updateOwnGui=True, **kw):
        """set data for this object:
check=1 : verify that all the keywords present can be handle by this func 
redo=1 : append self to viewer.objectsNeedingRedo
updateOwnGui=True : allow to update owngui at the end this func
"""
        sb = kw.get('stripBegin')
        se = kw.get('stripEnd')
        v = kw.get('vertices')

        if not sb is None or not v is None:
            # Modify self.stripBegin only if another set of vertices is
            # given or antother set of stripBegin.
            self.stripBegin = sb
            if sb is None:
                self.stripBegin = [0, len(v)]

            if (se and len(self.stripBegin) != len(se)) or se is None:
                self.stripEnd = map(lambda x,sB = self.stripBegin  :
                                    sB[x+1], range(len(self.stripBegin[:-1])))
                self.stripBegin = self.stripBegin[:-1]
            else:
                self.stripEnd = se
                
        faces = []
        # only self.stripBegin no sef.stripEnd or
        # len(self.stripBeg)!= len(sel.stripBeg)

        if self.stripBegin and self.stripEnd:
                for i in range(len(self.stripBegin)):
                    strip_ind = ()
                    begStrip = self.stripBegin[i]
                    endStrip = self.stripEnd[i]
                    for j in range(begStrip, endStrip):
                        strip_ind = strip_ind + (j,)
                    faces.extend([strip_ind])
##              else:
##                  for i in range(0, len(self.stripBegin)-1):
##                      strip_ind = ()
##                      for j in range(self.stripBegin[i], self.stripBegin[i+1]):
##                          strip_ind = strip_ind + (j,)
##                      print 'strip_ind', strip_ind
##                      faces.extend([strip_ind])

            
##          if self.stripBegin:
##              for i in range(0, len(self.stripBegin)-1):
##                  strip_ind = ()
##                  for j in range(self.stripBegin[i], self.stripBegin[i+1]):
##                      strip_ind = strip_ind + (j,)
##                  print 'strip_ind', strip_ind
##                  faces.extend([strip_ind])
                
                
        kw['faces'] = faces
	self.fnormals = kw.get( 'fnormals')

        redoFlags = apply( IndexedPolygons.Set, (self, check, 0), kw )

        return self.redoNow(redo, updateOwnGui, redoFlags)
	
    
    def MatBindingMode(self, propNum, face=GL.GL_FRONT):
        """ For the given face and material property corresponding to propNum,
   	the method sets the binding mode.  If there is one value, the mode is
        OVERALL; if the number of values equals the number of vertices, strips,
	or faces(triangles), the mode is PER_VERTEX, PER_PART, or PER_FACE,
        respectively. Else, mode is none. """

        #NONE, OVERALL, PER_VERTEX, PER_PART, PER_FACE = -1, 10, 11, 12, 13
        OVERALL, PER_VERTEX, PER_PART, PER_FACE = 10, 11, 12, 13
        num = propNum
	f = face

	nn = self.materials[f].prop[num].shape[0]
	self.inheritMaterial = 0

        if nn == 1:
            self.materials[f].binding[num] = OVERALL
        elif nn == len(self.vertexSet.vertices):
            self.materials[f].binding[num] = PER_VERTEX
        elif hasattr(self, 'faceSet') and nn == len(self.faceSet):
            self.materials[f].binding[num] = PER_PART
        elif hasattr(self, 'IndexedFaceSet') and nn==len(self.IndexedFaceSet):
            self.materials[f].binding[num] = PER_FACE
            self.shading = GL.GL_FLAT
            self.GetNormals()
        else:
            self.materials[f].binding[num] = -1
            self.inheritMaterial = 1

    def buildIndexedFaceSet(self):
        """Build the set of indices describing the strips"""
        
        f = []
        for stripNum in range(len(self.stripBegin)):
            start = self.stripBegin[stripNum]
            end = self.stripEnd[stripNum]
            for i in range(start, end-3,2):
                f.extend([(i, i+1, i+2), (i+2, i+1, i+3)])
            if (end-start) % 2== 1:
                f.extend([(end-3, end-2, end-1)])
        
##          if not self.stripEnd is None:
##              for stripNum in range(len(self.stripBegin)):
##                  print 'stripNum', stripNum
##                  start = self.stripBegin[stripNum]
##                  end = self.stripEnd[stripNum]
##                  for i in range(start, end-3,2):
##                      f.extend([(i, i+1, i+2), (i+2, i+1, i+3)])
##                  if (end-start) % 2== 1:
##                      f.extend([(end-3, end-2, end-1)])
##          else:
##              for stripNum in range(1,len(self.stripBegin)):
##                  print 'stripNum'
##                  for i in range(self.stripBegin[stripNum-1],
##                                 self.stripBegin[stripNum]-3, 2):
##                      f.extend([(i, i+1, i+2), (i+2, i+1, i+3)])
##                  if (self.stripBegin[stripNum]-self.stripBegin[stripNum-1]) % 2== 1:
##                      f.extend([(self.stripBegin[stripNum]-3,
##                                 self.stripBegin[stripNum]-2,
##                                 self.stripBegin[stripNum]-1)])
        self.IndexedFaceSet = DejaVu2.datamodel.FaceSet( f, (0,0) )

        
    def GetNormals(self):
        """ Gets the proper normals for smooth or flat shading.  Calls
        buildIndexedFaceSet(). Sets face normals or computes them if not given.
	Sets object normals to the vertex normals for smooth shading or to the
	face normals for flat shading. If shading none, normals are none."""
        if not self.stripBegin: return
        if self.shading==GL.GL_NONE:
            self.normals = None
        else:
            if hasattr(self, 'faceSet'):
                self.StripFaceSet = self.faceSet
                self.buildIndexedFaceSet()
                self.faceSet = self.IndexedFaceSet

		if self.fnormals is not None:
		    self.faceSet.normals.SetValues(self.fnormals)
                
		else:
                    self.FaceNormalFunction( self.ComputeFaceNormals )
                    self.faceSet.normals.ComputeMode( DejaVu2.viewerConst.AUTO )
                
                if self.shading==GL.GL_FLAT:
                    if hasattr(self, 'faceSet'):
                        self.normals = self.faceSet.normals.GetProperty()
                    else: self.normals = None
                elif self.shading==GL.GL_SMOOTH:
                    self.normals = self.vertexSet.normals.GetProperty()
                self.faceSet = self.StripFaceSet
            else: self.normals = None
        

    def isNewColor(self, c=None):
        """ Compares new color c to the current color.  If the same, method
        returns 0.  If the new color is different, the current color gets the
        values of the new color, and the method returns 1. """

        if c == None or not c.all():
            for i in self.currentCol:
                i = -1.0 # set an impossible color
            return 0
        elif abs(c[0]-self.currentCol[0]) < 0.0001 and \
             abs(c[1]-self.currentCol[1]) < 0.0001 and \
             abs(c[2]-self.currentCol[2]) < 0.0001 and \
             abs(c[3]-self.currentCol[3]) < 0.0001:
                return 0
        else:
            self.currentCol[0] = c[0]
            self.currentCol[1] = c[1]
            self.currentCol[2] = c[2]
            self.currentCol[3] = c[3]
            return 1

##      def isNewMaterial(self, face, prop, c):
##          """ For the given face (face) and property number (prop), the method
##          compares the new material value c to the current material. If
##          the same, method returns 0.  If different, the current material gets
##          the new material value, and the method returns 1. """

##          f = not(face==GL.GL_FRONT)
##          if not c:
##              for i in range(2):
##                  for j in range(5):
##                      for k in self.currentMat[i][j]:
##                          k = -1.0
##              return 0

##          elif abs(c[0]-self.currentMat[f][prop][0]) < 0.0001 and \
##               abs(c[1]-self.currentMat[f][prop][1]) < 0.0001 and \
##               abs(c[2]-self.currentMat[f][prop][2]) < 0.0001 and \
##               abs(c[3]-self.currentMat[f][prop][3]) < 0.0001:
##              return 0

##          else:
##              self.currentMat[f][prop][0] = c[0]
##              self.currentMat[f][prop][1] = c[1]
##              self.currentMat[f][prop][2] = c[2]
##              self.currentMat[f][prop][3] = c[3]
##              return 1


    def DisplayFunction(self):
        """ Either executes the present display list or creates a display
        list to display the triangle strip. """
        
        # if in select mode, switches to IndexedPolygons
        if self.viewer.currentCamera.renderMode == GL.GL_SELECT:
            temp = self.primitiveType
            if self.frontPolyMode == GL.GL_FILL:
                self.primitiveType = GL.GL_POLYGON
            elif self.frontPolyMode == GL.GL_LINE:
                self.primitiveType = GL.GL_LINE_LOOP
            elif self.frontPolyMode == GL.GL_POINT:
                self.primitiveType = GL.GL_POINTS
            self.faceSet = self.IndexedFaceSet
            IndexedPolygons.DisplayFunction(self)
            self.faceSet = self.StripFaceSet
            self.primitiveType = temp
            return
        
        if self.dpyList:
            Geom.DisplayFunction(self)


    def Draw(self):

        OVERALL, PER_VERTEX, PER_PART, PER_FACE = 10, 11, 12, 13
        #NONE, OVERALL, PER_VERTEX, PER_PART, PER_FACE = -1, 10, 11, 12, 13
        
        propConst = DejaVu2.viewerConst.propConst
        noCol = 1
        vert = self.vertexSet.vertices.array
        if len(vert)==0: return
        if not self.stripBegin or len(self.stripBegin)==0: return


        if self.materials[GL.GL_FRONT] and not self.inheritMaterial:
            frontMat = fpProp = []
            frontMatBind = fpBind = []
            mat = self.materials[GL.GL_FRONT]
            for propInd in range(4):
                b, p = mat.GetProperty(propInd)
                fpProp.append(p)
                fpBind.append(b)
            fpProp.append(mat.prop[4])
            fpBind.append(mat.binding[4])

        else:
            frontMat = None
            frontMatBind = None

        if self.materials[GL.GL_BACK] and not self.inheritMaterial:
            backMat = bpProp = []
            backMatBind = bpBind = []
            mat = self.materials[GL.GL_BACK]
            for propInd in range(4):
                b, p = mat.GetProperty(propInd)
                bpProp.append(p)
                bpBind.append(b)
            bpProp.append(mat.prop[4])
            bpBind.append(mat.binding[4])

        else:
            backMat = None
            backMatBind = None

        if not self.frontAndBack is None:
            face = GL.GL_FRONT
        else:
            face = GL.GL_FRONT_AND_BACK

        n = self.normals
        if self.normals is None:    # overall color for no normals or lighting
            if frontMat:
                if frontMatBind[noCol] == OVERALL:
                    GL.glColor4fv( frontMat[noCol][0] )
        else:
            if len(self.normals)==1:             # overall normal
                GL.glNormal3dv(n[0])
            if frontMat:
                for j in range(5):               # overall materials
                    if frontMatBind[j] == OVERALL:
                        glMaterialWithCheck( face, propConst[j],
                                             frontMat[j][0] )
            if backMat and not self.frontAndBack:
                for j in range(5):
                    if backMatBind[j] == OVERALL:
                        glMaterialWithCheck( GL.GL_BACK, propConst[j],
                                             backMat[j][0] )

        self.isNewColor()
        #self.isNewMaterial(0,0,0)

        # loop over each strip
        sB = self.stripBegin
        sE = self.stripEnd
        self.sE = sE
        self.sB= sB
        for stripNum in range(len(sB)):
            #print 'stripNum',stripNum
            GL.glPushName(stripNum)
            GL.glBegin(GL.GL_TRIANGLE_STRIP)

            # per part material properties
            if frontMat:
                if frontMatBind[noCol] == PER_PART:
                    if self.isNewColor(c=frontMat[noCol][stripNum]):
                        GL.glColor4fv(frontMat[noCol][stripNum])

            if n is not None:
                if frontMat:
                    for j in range(5):
                        if frontMatBind[j]==PER_PART:
                            glMaterialWithCheck( face,
                                                 propConst[j],
                                                 frontMat[j][stripNum] )

                if backMat and not self.frontAndBack:
                    for j in range(5):
                        if backMatBind[j] ==  PER_PART:
                            glMaterialWithCheck( GL.GL_BACK,
                                                 propConst[j],
                                                 backMat[j][stripNum] )

            #   loop over each vertex in a strip
            i = 0
            for ind in range(sB[stripNum],sE[stripNum]):
                if ind >= len(vert):
                    print 'ERROR',ind
                v = vert[ind]
                # normals for flat or smooth shading
                if n is not None:
                    if self.shading==GL.GL_FLAT:
                        if i > 1:
                            GL.glNormal3dv(n[sB[stripNum]+i-(2*(stripNum+1))])
                    elif self.shading==GL.GL_SMOOTH:
                        GL.glNormal3fv(n[ind])
                    else:
                        pass

                # per face (per triangle) material properties
                if n is None:
                    if frontMat:
                        if frontMatBind[noCol] == PER_FACE:
                            if i > 1:
                                if self.isNewColor(c=frontMat[noCol][sB[stripNum]+i-(2*(stripNum+1))]):
                                    GL.glColor4fv(frontMat[noCol][sB[stripNum]+i-(2*(stripNum+1))])

                else:
                    if frontMat:
                        for k in range(5):
                            if frontMatBind[k] == PER_FACE:
                                if i > 1:
                                    glMaterialWithCheck( face,
                                                         propConst[k],
                               frontMat[k][sB[stripNum]+i-(2*(stripNum+1))] )

                    if backMat and not self.frontAndBack:
                        for k in range(5):
                            if backMatBind[k] == PER_FACE:
                                if i > 1 and i%2==0:
                                    glMaterialWithCheck( GL.GL_BACK,
                                                         propConst[k],
                                backMat[k][sB[stripNum]+i-(2*(stripNum+1))] )


                #  per vertex material properties
                if n is None:
                    if frontMat:
                        if frontMatBind[noCol] == PER_VERTEX:
                            if self.isNewColor(c=frontMat[noCol][sB[stripNum]+i]):
                                GL.glColor4fv(frontMat[noCol][sB[stripNum]+i])

                else:
                    if frontMat:
                        for k in range(5):
                            if frontMatBind[k] == PER_VERTEX:
                                glMaterialWithCheck( face,
                                                     propConst[k],
                                            frontMat[k][sB[stripNum]+i] )
                    if backMat and not self.frontAndBack:
                        for k in range(5):
                            if backMatBind[k] == PER_VERTEX:
                                glMaterialWithCheck( GL.GL_BACK,
                                                     propConst[k],
                                             backMat[k][sB[stripNum]+i] )

                # draw vertex
                #GL.glVertex3dv(v)
                gllib.glVertex3fv(v)
                i = i + 1

            GL.glEnd()
            GL.glPopName()
        return 1
    


if __name__=='__main__':
    import pdb, numpy
    from DejaVu2 import Viewer
    vi = Viewer()
    vert = [(0, 1, 0), (1, 0, 1), (1, 1, 0), (2, 0, -1), (2, 1, 0), (3, 0, 1),
            (3, 1, 0), (4, 0, -1), (4, 1, 0), (5, 0, 1), (5, 3, 0), (6, -2, 1)]
    v1 = numpy.array(vert)
    v2 = v1 + 3.0
    v3 = numpy.concatenate( (v1,v2) )
    colors = [(0, 0, 1), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 0), (0, 1, 0),
              (0, 0, 1), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 0), (0, 1, 0)]
    
    v4 = v2 + 3.0
    v5 = numpy.concatenate( (v3,v4) )
    strip3 = Triangle_strip(vertices = v5, stripBegin=[0,12],
                            materials = [(0, 0, 1), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 0), (0, 1, 0),
              (0, 0, 1), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 0), (0, 1, 0)])
    #strip.Set(materials = [(1, 0, 0)])
    vi.AddObject(strip3)

    strip2 = Triangle_strip(vertices = v3, stripBegin=[0,12,24], materials = colors+colors[0:8])
    #strip2.Set(materials = [(0, 1, 0), (0, 0, 1)])
    vi.AddObject(strip2)
    
    strip = Triangle_strip(vertices = vert, materials = colors[2:12])
    #strip.Set(materials = [(1, 0, 0)])
    vi.AddObject(strip)
