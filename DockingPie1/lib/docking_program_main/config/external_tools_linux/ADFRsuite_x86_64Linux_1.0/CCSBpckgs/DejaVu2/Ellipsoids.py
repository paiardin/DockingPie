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
# Authors: Michel F. SANNER, Daniel Stoffler
#
#    sanner@scripps.edu
#    stoffler@scripps.edu
#    
# Copyright: M. Sanner, Daniel Stoffler TSRI 2000
#
#############################################################################


#
# $Header: /mnt/raid/services/cvs/DejaVu2/Ellipsoids.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#
# $Id: Ellipsoids.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#

import warnings
import numpy, math

from opengltk.OpenGL import GL, GLU
from opengltk.extent.utillib import glDrawSphereSet, extractedGlutSolidSphere

from DejaVu2.Geom import Geom
import DejaVu2.datamodel as datamodel
import DejaVu2.viewerConst as viewerConst
from DejaVu2.viewerFns import checkKeywords
from DejaVu2.colorTool import glMaterialWithCheck, resetMaterialMemory

class Ellipsoids(Geom):
    """Class for sets of spheres"""

    if glDrawSphereSet:
        fastSpheres = 1
    else:
        fastSpheres = 0

    keywords = Geom.keywords + [
        'centers',
        'quality',
        'scaling',
        'orientation',
        'slices',
        'stacks'
        ]

    def __init__(self, name=None, check=1, **kw):

        v = kw.get('centers')
        if v:
            kw['vertices'] = v     # rename centers in vertices for Geom.__init
        elif not kw.get('shape'):
            kw['shape'] = (0,3)    # default shape for sphere set

        self.culling = GL.GL_BACK
        self.inheritCulling = 0
        
        self.frontPolyMode = GL.GL_FILL
        self.inheritFrontPolyMode = viewerConst.NO
        self.lighting = viewerConst.YES
        
        self.slices = 5
        self.stacks = 5
        self.templateDSPL = None # (displayList, openglContext)

        self.radius = 1.0

        apply( Geom.__init__, (self, name, check), kw )
        assert len(self.vertexSet.vertices.ashape)==2
        
        self._modified = False


    def Set(self, check=1, redo=1, updateOwnGui=True, **kw):
        """set data for this object
check=1 : verify that all the keywords present can be handle by this func 
redo=1 : append self to viewer.objectsNeedingRedo
updateOwnGui=True : allow to update owngui at the end this func
"""
        redoFlags = 0

        # Exceptionnaly this has to be before the call to Geom.Set
        v = kw.get( 'centers')
        if v:
            kw['vertices'] = v     # rename centers in vertices for Geom.__init

        # Exceptionnaly this has to be before the call to Geom.Set
        # because we want to override the treatment of it by Geom.Set
        invertNormals = kw.get('invertNormals')
        if invertNormals is not None:
            kw.pop('invertNormals')
            if self.invertNormals != invertNormals:
                self.invertNormals = invertNormals
                self.chooseTemplate()
                redoFlags |= self._redoFlags['redoDisplayListFlag']

        redoFlags |= apply( Geom.Set, (self, check, 0), kw)

        if hasattr(self.vertexSet, 'scaling') is False:
            self.vertexSet.scaling = datamodel.ScalarProperties('scaling',
                           shape=(0,3), datatype=viewerConst.FPRECISION)

        if hasattr(self.vertexSet, 'orientation') is False:
            self.vertexSet.orientation = datamodel.ScalarProperties('orientation',
                           shape=(0, 4,4), datatype=viewerConst.FPRECISION)

        scal = kw.get( 'scaling')
        if scal != None:
            self.scaling = numpy.array(scal).astype('f')
            self.scaling.shape = (-1, 3)
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        orient = kw.get( 'orientation')
        if orient != None:
            self.orientation = numpy.array(orient).astype('f')
            self.orientation.shape = (-1, 4, 4)
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        qual = kw.get( 'quality')
        if qual != None:
            if qual > 2:
                self.slices = qual
                self.stacks = qual
                if self.templateDSPL is not None:
                    redoFlags |= self._redoFlags['redoTemplateFlag']
                    redoFlags |= self._redoFlags['redoDisplayListFlag']

        slices = kw.get( 'slices')
        if slices != None:
            if slices > 2:
                self.slices = slices
                if self.templateDSPL is not None:
                    redoFlags |= self._redoFlags['redoTemplateFlag']
                    redoFlags |= self._redoFlags['redoDisplayListFlag']

        stacks = kw.get( 'stacks')
        if stacks != None:
            if stacks > 2:
                self.stacks = stacks
                if self.templateDSPL is not None:
                    redoFlags |= self._redoFlags['redoTemplateFlag']
                    redoFlags |= self._redoFlags['redoDisplayListFlag']

        return self.redoNow(redo, updateOwnGui, redoFlags)


    def Add(self, check=1, redo=1, **kw):
        """Add spheres"""

        if __debug__:
            if check:
                apply( checkKeywords, (self.name,self.keywords), kw)

        v = kw.get( 'centers')
        if v:
            kw['vertices'] = v     # rename centers in vertices for Geom.__init
            apply( Geom.Add, (self,0,0), kw)

        if self.viewer and redo:
            if self.redoDspLst:
                self.viewer.objectsNeedingRedo[self] = None
#                self.RedoDisplayList()


    def deleteTemplate(self):
        #print "Ellipsoids.deleteTemplate"
        # it is asumed the right OpenGL context is active
        assert self.templateDSPL is not None
        currentcontext = self.viewer.currentCamera.getContext()
        if currentcontext != self.templateDSPL[1]:
            import traceback;traceback.print_stack()
            warnings.warn('deleteTemplate failed because the current context is the wrong one')
            print "currentcontext != self.templateDSPL[1]", currentcontext, self.templateDSPL[1]
        else:
            #print '-%d'%self.templateDSPL[0], currentcontext, "glDeleteLists Ellipsoids"
            #print '-%d'%(self.templateDSPL[0]+1), currentcontext, "glDeleteLists Ellipsoids"
            #print '-%d'%(self.templateDSPL[0]+2), currentcontext, "glDeleteLists Ellipsoids"
            GL.glDeleteLists(self.templateDSPL[0], 3)
            self.templateDSPL = None


    def makeTemplate(self):
        #print "Ellipsoids.makeTemplate"
        # it is asumed the right OpenGL context is active
        assert self.templateDSPL is None
        lFirstList = GL.glGenLists(3)
        #print "lFirstList Ellipsoids.makeTemplate", lFirstList, self.name
        lCurrentContext = self.viewer.currentCamera.getContext()
        self.templateDSPL = ( lFirstList, lCurrentContext )

        GL.glNewList(lFirstList+1, GL.GL_COMPILE)
        #print '+%d'%(lFirstList+1), lCurrentContext, "glNewList Ellipsoids1"
        extractedGlutSolidSphere(1, self.slices, self.stacks, 0)
        #print '*%d'%GL.glGetIntegerv(GL.GL_LIST_INDEX), "glEndList Ellipsoids1"
        GL.glEndList()
        
        GL.glNewList(lFirstList+2, GL.GL_COMPILE)
        #print '+%d'%(lFirstList+2), lCurrentContext, "glNewList Ellipsoids2"
        extractedGlutSolidSphere(1, self.slices, self.stacks, 1)
        #print '*%d'%GL.glGetIntegerv(GL.GL_LIST_INDEX), "glEndList Ellipsoids2"
        GL.glEndList()

        self.chooseTemplate()
        

    def redoTemplate(self):
        self.deleteTemplate()
        self.makeTemplate()


    def chooseTemplate(self):
        GL.glNewList(self.templateDSPL[0], GL.GL_COMPILE)
        #print '+%d'%self.templateDSPL[0], "glNewList Ellipsoids0"
        if self.invertNormals:
            #print "GLU_INSIDE reversed normals"
            #print '#%d'%(self.templateDSPL[0]+2), "glCallList Ellipsoids2"
            GL.glCallList(self.templateDSPL[0]+2)
        else:
            #print "GLU_OUTSIDE regular normals"
            #print '#%d'%(self.templateDSPL[0]+1), "glCallList Ellipsoids1"
            GL.glCallList(self.templateDSPL[0]+1)
        #print '*%d'%GL.glGetIntegerv(GL.GL_LIST_INDEX), "glEndList Ellipsoids0"
        GL.glEndList()


    def Draw(self):
        """ Draw function of the geom
return status 0 or 1
If you want fast rendering, you need to set self.templateDSPL
using MakeTemplate.
"""
        #print "Ellipsoids.Draw"
        
        assert self.templateDSPL is not None

        currentcontext = self.viewer.currentCamera.getContext()
        if currentcontext != self.templateDSPL[1]:
            warnings.warn("""draw failed because the current context is the wrong one""")
            #print "currentcontext != self.templateDSPL[1]", currentcontext, self.templateDSPL[1]
            return 0

        centers = self.vertexSet.vertices.array
        if len(centers) == 0: return
        scaling = self.scaling
        orientation = self.orientation
        vertices = self.vertexSet.vertices.array
        
        if len(vertices) != len(scaling) or len(vertices) != len(orientation):
            return

        if self.inheritMaterial:
            fp = None
            bp = None
        else:
            fp = self.materials[GL.GL_FRONT]
            if not self.frontAndBack:
                bp = self.materials[GL.GL_BACK]
                face = GL.GL_FRONT
            else:
                bp = None
                face = GL.GL_FRONT_AND_BACK

        for i in xrange(len(vertices)):
            GL.glPushName(i)

            if fp:
                for m in (0,1,2,3,4):
                    if fp.binding[m] != viewerConst.OVERALL:
                        glMaterialWithCheck( face,
                                             viewerConst.propConst[m],
                                             fp.prop[m][i] )
            if bp:
                for m in (0,1,2,3,4):
                    if bp.binding[m] != viewerConst.OVERALL:
                        glMaterialWithCheck( face,
                                             viewerConst.propConst[m],
                                             bp.prop[m][i] )

            GL.glPushMatrix()
            GL.glTranslatef(float(vertices[i][0]),
                            float(vertices[i][1]),
                            float(vertices[i][2]))
            GL.glMultMatrixf( orientation[i].ravel() )
            GL.glScalef(float(scaling[i][0]),
                        float(scaling[i][1]),
                        float(scaling[i][2]))
            #print '#%d'%self.templateDSPL[0], "glCallList Ellipsoids0"
            GL.glCallList(self.templateDSPL[0])
            GL.glPopMatrix()

            GL.glPopName()

        return True


    def asIndexedPolygons(self, run=1, quality=None, **kw):
        """ run=0 returns 1 if this geom can be represented as an
        IndexedPolygon and None if not. run=1 returns the IndexedPolygon
        object."""

        if run==0:
            return 1 # yes, I can be represented as IndexedPolygons
        
        if quality is None:
            quality = 2
        
        # get centers
        centers = self.vertexSet.vertices.array

        # get radii
        if self.oneRadius == viewerConst.NO:
            radii = self.vertexSet.radii.array
        else:
            radii = numpy.ones( centers.shape[0] ) * self.radius
        
        # create template sphere
        S = TriangulateIcosByEdgeCenterPoint(quality=quality)
        tmpltVertices = S.getVertices(quality=quality)
        tmpltFaces = S.getFaces(quality=quality)
        tmpltNormals = S.getVNormals(quality=quality)

        # these lists will store the data for the new spheres
        vertices = []
        faces = []
        normals = []

        # loop over spheres
        for i in range(len(centers)):
            vert = numpy.array(tmpltVertices[:])*radii[i] + centers[i]
            vertices.extend(list(vert))
            fac = numpy.array(tmpltFaces[:]) + i*len(tmpltVertices)
            faces.extend(list(fac))
            norm = numpy.array(tmpltNormals[:])
            normals.extend(list(norm))

        from DejaVu2.IndexedPolygons import IndexedPolygons
        sphGeom = IndexedPolygons("sph", vertices=numpy.array(vertices),
                               faces=faces, vnormals=numpy.array(normals),
                               visible=1, invertNormals=self.invertNormals)

        # copy Spheres materials into sphGeom
        matF = self.materials[GL.GL_FRONT]
        matB = self.materials[GL.GL_BACK]
        sphGeom.materials[GL.GL_FRONT].binding = matF.binding[:]
        sphGeom.materials[GL.GL_FRONT].prop = matF.prop[:]
        sphGeom.materials[GL.GL_BACK].binding = matB.binding[:]
        sphGeom.materials[GL.GL_BACK].prop = matB.prop[:]

        if sphGeom.materials[GL.GL_FRONT].binding[1] == viewerConst.PER_VERTEX:
            newprop = []
            index = 0
            cnt = 0
            for i in range(len(vertices)):
                newprop.append(sphGeom.materials[GL.GL_FRONT].prop[1][index])
                cnt = cnt + 1
                if cnt == len(tmpltVertices):
                    index = index + 1
                    cnt = 0
            
            sphGeom.materials[GL.GL_FRONT].prop[1] = newprop         
        return sphGeom


class TriangulateIcos:
    """Base class to compute vertices, faces and normals of a sphere based
    on icosahedral subdivision. Subclassed will implement different
    subdivision methods.
    A quality can be passed to the constructur which will trigger the
    precomputation of spheres of quality 0 to quality.
    To access the data, use getVertices(quality=val), getFaces(quality=val),
    getVNormals(quality=val) where val is the quality level """
    

    def __init__(self, quality=None):
        if quality is None:
            quality = 5 # set default to 5
        self.quality = quality

        self.vertices=[] # stores vertices
                         # face lists are created dynamically later on
                         # normals == vertices
        
        X = 0.525731112119133606 # X coord
        Z = 0.850650808352039932 # Y coord

        # build initial icosahedron (lowest quality)
        self.vertices = [
            [-X, 0., Z], [X, 0., Z], [-X, 0., -Z], [X, 0., -Z],
            [0., Z, X], [0., Z, -X], [0., -Z, X], [0., -Z, -X],
            [Z, X, 0.], [-Z, X, 0.], [Z, -X, 0.], [-Z, -X, 0.]
            ]

        self.facesQ0 = [
            [11,6,0], [9,11,0], [0,6,1], [6,10,1], [9,5,2], [7,11,2], [5,3,2],
            [8,10,3], [5,8,3], [0,1,4], [9,4,5], [4,8,5], [7,10,6], [2,3,7],
            [4,1,8], [0,4,9], [8,1,10], [7,3,10], [7,6,11], [9,2,11]
            ]


    def subsample(self, vertices, faces, quality):
        lenF = len(faces)
        lenV = len(vertices)

        for i in xrange(lenF):
            v0 = vertices[faces[i][0]]
            v1 = vertices[faces[i][1]]
            v2 = vertices[faces[i][2]]
            self.subdivideVert(v0, v1, v2)

        self.subdivideFaces(faces, lenV, quality)
    

    def normalize(self, v):
        d = math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
        if d == 0.0:
            print 'Zero length vector!'
            return
        return [v[0] / d, v[1] / d, v[2] / d]


    def getVertices(self, quality=0):
        """ has to be implemented by subclass """
        pass

    def getVNormals(self, quality=0):
        """ has to be implemented by subclass """
        pass
    
    def getFaces(self, quality=0):
        return getattr(self, 'facesQ%d'%quality)


class TriangulateIcosByEdgeCenterPoint(TriangulateIcos):
    
    def __init__(self, quality=None):
        TriangulateIcos.__init__(self, quality)
    
        if self.quality > 0:
            for qu in range(1, self.quality+1):
                self.subsample(self.vertices,
                               getattr(self, 'facesQ%d'%(qu-1,) ),
                               qu)


    def subdivideVert(self, v0, v1, v2):
        # called by subsample
        v01 = []
        v12 = []
        v20 = []
        
        for i in range(3):
            v01.append(v0[i] + v1[i])
            v12.append(v1[i] + v2[i])
            v20.append(v2[i] + v0[i])

        v01=self.normalize(v01)
        v12=self.normalize(v12)
        v20=self.normalize(v20)

        self.vertices.append(v01)
        self.vertices.append(v12)
        self.vertices.append(v20)
        

    def subdivideFaces(self, faces, lenV, quality):
        # called by subsample
        newFaces = []
        
        for i in xrange(len(faces)):
            j = i
            j = j * 3
            f0 = faces[i][0]
            f1 = faces[i][1]
            f2 = faces[i][2]
            f01 = j+lenV
            f12 = j+lenV+1
            f20 = j+lenV+2
            newFaces.append([f0, f01, f20])
            newFaces.append([f01, f12, f20])
            newFaces.append([f01, f1, f12])
            newFaces.append([f20, f12, f2])

        # dynamically create a self.facesQ<quality>
        setattr(self, 'facesQ%d'%quality, newFaces[:])


    def getVertices(self, quality=0):
        # the vertex list is very big, since vertices are added to this
        # list after every subsampling. Thus, only what is needed is returned
        v = 12 # vertices of icosahedron
        f = 20 # faces of icosahedron
        
        for i in range(quality):
            v = v+f*3
            f = f*4
        return self.vertices[:v]


    def getVNormals(self, quality=0):
        # normals == vertices
        self.normals = self.getVertices(quality=quality)[:]
        return self.normals


class TriangulateIcosByFaceCenterPoint(TriangulateIcos):
    """ This class subdivides each face in 3 new faces by putting a center
    in the middle of a face triangle."""

    def __init__(self, quality=None):
        TriangulateIcos.__init__(self, quality)
        
        if self.quality > 0:
            for qu in range(1, self.quality+1):
                self.subsample(self.vertices,
                               getattr(self, 'facesQ%d'%(qu-1,) ),
                               qu)


    def subdivideVert(self, v0, v1, v2):
        # called by subsample
        v012 = []
                
        for i in range(3):
            v012.append(v0[i] + v1[i] + v2[i])

        v012 = self.normalize(v012)
        self.vertices.append(v012)


    def subdivideFaces(self, faces, lenV, quality):
        # called by subsample
        newFaces = []
        
        for i in xrange(len(faces)):
            f0 = faces[i][0]
            f1 = faces[i][1]
            f2 = faces[i][2]
            f012 = i+lenV
            newFaces.append([f0, f1, f012])
            newFaces.append([f1, f2, f012])
            newFaces.append([f2, f0, f012])

        # dynamically create a self.facesQ<quality>
        setattr(self, 'facesQ%d'%quality, newFaces[:])


    def getVertices(self, quality=0):
        # the vertex list is very big, since vertices are added to this
        # list after every subsampling. Thus, only what is needed is returned
        v = 12 # vertices of icosahedron
        f = 20 # faces of icosahedron
        
        for i in range(quality):
            v = v+f
            f = f*3
        return self.vertices[:v]


    def getVNormals(self, quality=0):
        # normals == vertices
        return self.getVertices(quality=quality)



## rotcrn = numpyarray( [
##     [0.62831, -0.67899, -0.37974, 0],
##     [0.36979,  0.69011, -0.62210, 0],
##     [0.68446,  0.25045,  0.68468, 0],
##     [0, 0, 0, 1]], 'f')
## transcrn = [9.26883, 9.78728, 6.96709]
## scalecrn = [8.0957, 11.7227, 16.2550]

## rotcv = numpyidentity(4, 'f')
## transcv = [0.,0.,0.]
## scalecv = [7.8895, 5.5336, 3.3147]

## from DejaVu2.Ellipsoids import Ellipsoids
## el1 = Ellipsoids('ellipsoids1', centers=[transcv,transcrn],
##                 scaling = [scalecv, scalecrn],
##                 orientation = [rotcv, rotcrn],
##                 materials=[[1,0,0],[0,1,0]], inheritMaterial=0,
##                 quality=30)
## mv.GUI.VIEWER.AddObject(el1)

## USING A SPHERE FOR 1 ELLIPSOID
#from DejaVu2.Spheres import Spheres
#s = Spheres('ellipsoid1', vertices=[[0.,0.,0.]], radii=[1.], quality=30)
#mv.GUI.VIEWER.AddObject(s)
#s.setMatrixComponents( rotcrn, transcrn, scalecrn)
