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
# $Header: /mnt/raid/services/cvs/DejaVu2/Spheres.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#
# $Id: Spheres.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#

import warnings
import numpy, math

from opengltk.OpenGL import GL, GLU
from opengltk.extent.utillib import glDrawSphereSet, extractedGlutSolidSphere

import DejaVu2
from DejaVu2.IndexedGeom import IndexedGeom
from DejaVu2 import datamodel, viewerConst
from DejaVu2.viewerFns import checkKeywords
from DejaVu2.colorTool import glMaterialWithCheck, resetMaterialMemory
from DejaVu2.IndexedPolygons import IndexedPolygons

if hasattr(DejaVu2, 'enableVertexArrayNonVBO') is False:
    DejaVu2.enableVertexArrayNonVBO = False

try:
    from UTpackages.UTimposter import utimposterrend
    UTImposterRendererFound = True
except ImportError:
    #warnings.warn('UTpackages.UTimposter not found')
    UTImposterRendererFound = False

    
class GLUSpheres(IndexedGeom):
    """Class for sets of spheres"""

    if glDrawSphereSet:
        fastSpheres = 1
    else:
        fastSpheres = 0

    keywords = IndexedGeom.keywords + [
        'centers',
        'radii',
        'quality',
        'slices', #deprecated
        'stacks'  #deprecated
        ]

    def getState(self, full=False):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        state = IndexedGeom.getState(self, full)
        state['quality'] = self.quality

        if full:
            rad = self.vertexSet.radii.array
            if len(rad):
                state['radii'] =  rad
            
        return state


    def __init__(self, name=None, check=1, **kw):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        #print "Spheres.__init__"
        v = kw.get('centers')
        if v is not None:
            kw['vertices'] = v     # rename centers in vertices for Geom.__init
        elif not kw.get('shape'):
            kw['shape'] = (0,3)    # default shape for sphere set

        self.templateDSPL = None # (displayList, openglContext)
        #self.firstList = GL.glGenLists(3)

        self.culling = GL.GL_BACK
        self.inheritCulling = 0

        self.frontPolyMode = GL.GL_FILL
        self.inheritFrontPolyMode = viewerConst.NO

        self.oneRadius = viewerConst.YES
        self.radius = 1.0

        self.quality = None

        #self.immediateRendering = True

        apply( IndexedGeom.__init__, (self, name, check), kw )
        assert len(self.vertexSet.vertices.ashape)==2
        
        self._modified = False


    def Set(self, check=1, redo=1, updateOwnGui=True, **kw):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """set data for this object
check=1 : verify that all the keywords present can be handle by this func 
redo=1 : append self to viewer.objectsNeedingRedo
updateOwnGui=True : allow to update owngui at the end this func
"""
        redoFlags = 0
        
        # Exceptionnaly this has to be before the call to Geom.Set
        v = kw.pop( 'centers', None)
        if v is not None:
            kw['vertices'] = v # rename centers in vertices for Geom.__init

        # Exceptionnaly this has to be before the call to Geom.Set
        # because we want to override the treatment of it by Geom.Set
        invertNormals = kw.pop('invertNormals', None)
        if invertNormals is not None:
            if self.invertNormals != invertNormals:
                self.invertNormals = invertNormals
                self.chooseTemplate()
                redoFlags |= self._redoFlags['redoDisplayListFlag']

        redoFlags |= apply( IndexedGeom.Set, (self, check, 0), kw)
        if len(kw) == 0:
            return self.redoNow(redo, updateOwnGui, redoFlags)

        rad = kw.pop('radii', None)
        if rad is not None:
            if type(rad).__name__ in ('float','int'):
                self.oneRadius = viewerConst.YES
                self.radius = rad
                self.vertexSet.radii = datamodel.ScalarProperties('radii',
                           shape=(0,), datatype=viewerConst.FPRECISION)
            else: # type(rad).__name__ in ('list', 'tuple', 'ndarray'):
                if len(rad)==1:
                    self.oneRadius = viewerConst.YES
                    self.radius = rad[0]
                self.vertexSet.radii = datamodel.ScalarProperties('radii', rad,
                              datatype=viewerConst.FPRECISION)
        elif hasattr(self.vertexSet, 'radii') is False:
            self.vertexSet.radii = datamodel.ScalarProperties(
                                        'radii',
                                        shape=(0,), 
                                        datatype=viewerConst.FPRECISION)

        if rad is not None or v is not None:
            redoFlags |= self._redoFlags['redoDisplayListFlag']
            self.vertexSet.radii.PropertyStatus(len(self.vertexSet))
            if self.vertexSet.radii.status < viewerConst.COMPUTED:
                self.oneRadius = viewerConst.YES
            else:
                self.oneRadius = viewerConst.NO

        kw.pop( 'stacks', None) # deprecated
        kw.pop( 'slices', None) # deprecated
        quality = kw.pop('quality', None)
        if quality is not None or self.quality not in [1,2,3,4,5]:
            lOldQuality = self.quality
            if quality in [1,2,3,4,5]:
                self.quality = quality
            else:
                if len(self.vertexSet.vertices.array) < 500:
                    self.quality = 5
                elif len(self.vertexSet.vertices.array) < 5000:
                    self.quality = 4
                elif len(self.vertexSet.vertices.array) < 50000:
                    self.quality = 3
                elif len(self.vertexSet.vertices.array) < 100000:
                    self.quality = 2
                else:
                    self.quality = 1
            if lOldQuality != self.quality:
                if self.templateDSPL is not None:
                    redoFlags |= self._redoFlags['redoTemplateFlag']
                    redoFlags |= self._redoFlags['redoDisplayListFlag']

        return self.redoNow(redo, updateOwnGui, redoFlags)



    def deleteTemplate(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        #print "Spheres.deleteTemplate", self.templateDSPL
        # it is asumed the right OpenGL context is active
        if GL.glGetIntegerv(GL.GL_LIST_INDEX) == [0]:
            assert self.templateDSPL is not None
            currentcontext = self.viewer.currentCamera.getContext()
            if currentcontext != self.templateDSPL[1]:
                import traceback;traceback.print_stack()
                warnings.warn('deleteTemplate failed because the current context is the wrong one')
                print "currentcontext != self.templateDSPL[1]", currentcontext, self.templateDSPL[1]
                pass
            else:
                #print '-%d'%self.templateDSPL[0], currentcontext, "glDeleteLists Spheres0"
                #print '-%d'%(self.templateDSPL[0]+1), currentcontext, "glDeleteLists Spheres1"
                #print '-%d'%(self.templateDSPL[0]+2), currentcontext, "glDeleteLists Spheres2"
                GL.glDeleteLists(self.templateDSPL[0], 3)
                self.templateDSPL = None


    def makeTemplate(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        #print "Spheres.makeTemplate", self.quality
        # it is asumed the right OpenGL context is active
        # make sure we are not already in a newlist
        if GL.glGetIntegerv(GL.GL_LIST_INDEX) == [0]:
            assert self.templateDSPL is None
            lFirstList = GL.glGenLists(3)
            #lFirstList = self.firstList
            #print "Spheres.makeTemplate", lFirstList
            #print "lFirstList Spheres.makeTemplate", lFirstList, self.name
            lCurrentContext = self.viewer.currentCamera.getContext()
            self.templateDSPL = ( lFirstList, lCurrentContext )

            if (hasattr(DejaVu2, 'enableVBO') and DejaVu2.enableVBO) \
              or DejaVu2.enableVertexArrayNonVBO is True :
                lSphereGeom = self.asIndexedPolygons(run=1, quality=self.quality-1,
                                                     centers=((0,0,0),), radii=((1.,),) )
                GL.glNewList(lFirstList+1, GL.GL_COMPILE)
                lSphereGeom.Draw()
                GL.glEndList()
                lSphereGeom.Set(invertNormals=True)
                GL.glNewList(lFirstList+2, GL.GL_COMPILE)
                lSphereGeom.Draw()
                GL.glEndList()
            else:
                quality = self.quality * 5
                GL.glNewList(lFirstList+1, GL.GL_COMPILE)
                #print '+%d'%(lFirstList+1), lCurrentContext, "glNewList Spheres1"
                extractedGlutSolidSphere(1, quality, quality, 0)
                #print '*%d'%GL.glGetIntegerv(GL.GL_LIST_INDEX), "glEndList Spheres1"
                GL.glEndList()
                GL.glNewList(lFirstList+2, GL.GL_COMPILE)
                #print '+%d'%(lFirstList+2), lCurrentContext, "glNewList Spheres2"
                extractedGlutSolidSphere(1, quality, quality, 1)
                #print '*%d'%GL.glGetIntegerv(GL.GL_LIST_INDEX), "glEndList Spheres2"
                GL.glEndList()

            self.chooseTemplate()


    def redoTemplate(self):
        if __debug__:
            if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        
        if self.viewer:
            lSuspendRedraw = self.viewer.suspendRedraw
            self.viewer.suspendRedraw = True
        self.deleteTemplate()
        self.makeTemplate()
        if self.viewer:
            self.viewer.suspendRedraw = lSuspendRedraw


    def chooseTemplate(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        # make sure we are not already in a newlist
        if GL.glGetIntegerv(GL.GL_LIST_INDEX) == [0]:
            GL.glNewList(self.templateDSPL[0], GL.GL_COMPILE)
            #print '+%d'%self.templateDSPL[0], "glNewList Spheres0"
            if self.invertNormals:
                #print "GLU_INSIDE reversed normals"
                #print '#%d'%(self.templateDSPL[0]+2), "glCallList Spheres2"
                GL.glCallList(self.templateDSPL[0]+2)
            else:
                #print "GLU_OUTSIDE regular normals"
                #print '#%d'%(self.templateDSPL[0]+1), "glCallList Spheres1"
                GL.glCallList(self.templateDSPL[0]+1)
            #print '*%d'%GL.glGetIntegerv(GL.GL_LIST_INDEX), "glEndList Spheres0"
            GL.glEndList()
        

    def Add(self, check=1, redo=1, **kw):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Add spheres"""

        if __debug__:
            if check:
                apply( checkKeywords, (self.name,self.keywords), kw)

	v = kw.get( 'centers')
	if v:
	    kw['vertices'] = v     # rename centers in vertices for Geom.__init
        apply( IndexedGeom.Add, (self,0,0), kw)

	rad = kw.get( 'radii')
	if rad:
	    if type(rad).__name__ == 'float':
		self.oneRadius = viewerConst.YES
		self.radius = rad
	    else:
		self.vertexSet.radii.AddValues( rad )

	if rad or v:
            self.redoDspLst=1
	    self.vertexSet.radii.PropertyStatus(len(self.vertexSet))
	    if self.vertexSet.radii.status < viewerConst.COMPUTED:
		self.oneRadius = viewerConst.YES
	    else:
		self.oneRadius = viewerConst.NO

        if self.viewer and redo:
            if self.redoDspLst and self not in self.viewer.objectsNeedingRedo:
                self.viewer.objectsNeedingRedo[self] = None
#                self.RedoDisplayList()


    def Draw(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Draw function of the geom
return status 0 or 1
If you want fast rendering, you need to set self.templateDSPL
using MakeTemplate.
"""
        #print "Spheres.Draw", self.name

        assert self.templateDSPL is not None

        currentcontext = self.viewer.currentCamera.getContext()
        if currentcontext != self.templateDSPL[1]:
            # happens al lthe time on windows
            #import traceback;traceback.print_stack()
            #warnings.warn("""draw failed because the current context is the wrong one""")
            ##print "currentcontext != self.templateDSPL[1]", currentcontext, self.templateDSPL[1]
            return 0
            
        centers = self.vertexSet.vertices.array
        if len(centers) == 0: 
            return 0

        # handle overall binding of material
        if self.inheritMaterial:
            fp = None
            fpProp = None
            bp = None
            bpProp = None
        else:
            mat = self.materials[GL.GL_FRONT]
            fpProp = []
            bpProp = None
            for propInd in range(4):
                b, p = mat.GetProperty(propInd)
                fpProp.append(p)
            fpProp.append(mat.prop[4])

            fp = self.materials[GL.GL_FRONT]
            #colorFront = numpy.array(self.materials[GL.GL_FRONT].prop[1], copy=1)
            colorFront = numpy.array(fpProp[1], copy=1)
            if self.frontAndBack:
                bp = None
                face = GL.GL_FRONT_AND_BACK
            else:
                bp = self.materials[GL.GL_BACK]
                face = GL.GL_FRONT
                matb = self.materials[GL.GL_BACK]
                bpProp = []
                for propInd in range(4):
                    b, p = matb.GetProperty(propInd)
                    bpProp.append(p)
                bpProp.append(matb.prop[4])

        if fp:
            for m in (0,1,2,3,4):
                if fp.binding[m] == viewerConst.OVERALL:
                    glMaterialWithCheck( face,
                                         viewerConst.propConst[m],
                                         fpProp[m][0])
            if fp.binding[1] == viewerConst.OVERALL:
                GL.glColor4fv(colorFront[0])

            if fp:
                for m in (0,1,2,3,4):
                    if fp.binding[m] != viewerConst.OVERALL:
                        glMaterialWithCheck( face,
                                             viewerConst.propConst[m],
                                             fpProp[m][0])

                if fp.binding[1] != viewerConst.OVERALL:
                    GL.glColor4fv(colorFront[0])
            if bp:
                for m in (0,1,2,3,4):
                    if bp.binding[m] != viewerConst.OVERALL:
                        glMaterialWithCheck( GL.GL_BACK,
                                             viewerConst.propConst[m],
                                             bp.prop[m][0])

        #print self.name
        #if fp: print fp.prop[1], fp.binding
        #else: print
        
        if self.fastSpheres:
            #print "self.fastSpheres", self.fastSpheres
            if self.oneRadius == viewerConst.NO:
                radii = self.vertexSet.radii.array
                #FIXME: quick fix because can be called from base class Set
                # method after centers have been set BUT before radii have been
                # set
                if len(self.vertexSet.vertices) != len(radii):
                    return 0
            else:
                radii = numpy.ones( centers.shape[0] ) * self.radius
            radii.shape = (-1,1)
            coords = numpy.concatenate ( (centers, radii), 1 )

##             if not self.inheritMaterial:
##                 mat = self.materials[GL.GL_FRONT]
##                 fpProp = []
##                 for propInd in range(4):
##                     b, p = mat.GetProperty(propInd)
##                     fpProp.append(p)
##                 fpProp.append(mat.prop[4])
##                 #fpProp = self.materials[GL.GL_FRONT].prop[:5]
##             else:
##                 fpProp = None

            #print 'FUGU OVERWRITE COLOR', fpProp
            #import numpy
            #GL.glMaterialfv(GL.GL_FRONT, GL.GL_AMBIENT, numpy.array((.6,.6,.6,1), 'f'))
            #GL.glMaterialfv(GL.GL_FRONT, GL.GL_DIFFUSE, numpy.array((1.,1.,1.,1), 'f'))
            #GL.glMaterialfv(GL.GL_FRONT, GL.GL_SPECULAR, numpy.array((.4,.4,.4,1), 'f'))
            #GL.glMaterialfv(GL.GL_FRONT, GL.GL_EMISSION, numpy.array((0,0,0,1), 'f'))
            #GL.glMaterialf(GL.GL_FRONT, GL.GL_SHININESS, 1.)
            if len(self.faceSet.faces)==0:
                faces = range(len(coords))
            else:
                faces = self.faceSet.faces.array
            status = glDrawSphereSet( 
                self.templateDSPL[0],
                coords.astype('f'),
                faces,
                fpProp, bpProp,
                frontAndBack=int(self.frontAndBack),
                noLightCol=1,
                highlight=self.highlight,
                )
            #print "Spheres, status: ", status
            return status
        else:
            resetMaterialMemory()
            #print "SLOW Spheres"
            if self.oneRadius == viewerConst.NO:
                radii = self.vertexSet.radii.array
            else:
                radii = numpy.ones( centers.shape[0] ) * self.radius

            if len(self.vertexSet.vertices) != len(radii):
                return 0
            
            for i in xrange(centers.shape[0]):
                GL.glPushName(i)
                GL.glPushMatrix()
                GL.glTranslatef(float(centers[i][0]),
                                float(centers[i][1]),
                                float(centers[i][2]))
                if not self.oneRadius:
                    GL.glScalef(float(radii[i]),float(radii[i]),float(radii[i]))
                else:
                    GL.glScalef(float(self.radius), float(self.radius), float(self.radius))
                #print '#%d'%self.templateDSPL[0], "glCallList Spheres0"
                if fp:
                    for m in (0,1,2,3,4):
                        if fp.binding[m] != viewerConst.OVERALL:
                            glMaterialWithCheck( face,
                                                 viewerConst.propConst[m],
                                                 fp.prop[m][0], geom=self)
                GL.glCallList(self.templateDSPL[0])
                GL.glPopMatrix()
                GL.glPopName()
            return 1


    def asIndexedPolygons(self, 
                          run=1,
                          quality=None,
                          centers=None,
                          radii=None,
                          **kw):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """implement the sphere as an icosaedre.
run=0 returns 1 if this geom can be represented as an
IndexedPolygon and None if not. run=1 returns the IndexedPolygon object.
"""

        #print "Spheres.asIndexedPolygons", quality
        if run==0:
            return 1 # yes, I can be represented as IndexedPolygons
        
        if quality in [1,2,3,4,5]:
            quality = quality
        elif quality < 0 and self.quality in [2,3,4,5]:
            quality = self.quality - 2
        else:
            quality = self.quality - 1

        # get centers
        if centers is None:
            centers = self.vertexSet.vertices.array

        # get radii
        if radii is None:
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

        #print "Spheres.asIndexedPolygons out", quality
        return sphGeom



if UTImposterRendererFound:
    class UTSpheres(GLUSpheres):

        keywords = list(GLUSpheres.keywords)
        for kw in ['quality', 'stacks', 'slices']:
            keywords.remove(kw)


        def makeTemplate(self):
            if __debug__:
             if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
            pass


        def __init__(self, name=None, check=1, **kw):
            if __debug__:
             if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
            apply( GLUSpheres.__init__, (self, name, check), kw)


        def __init__(self, name=None, check=1, **kw):
            if __debug__:
             if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
            self.imposterRenderer = utimposterrend.ImposterRenderer()

            apply(GLUSpheres.__init__, (self, name, 0), kw )
            
            self.immediateRendering = 1


        def Set(self, check=1, redo=1, updateOwnGui=True, **kw):
            if __debug__:
             if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
            """set data for this object
check=1 : verify that all the keywords present can be handle by this func 
redo=1 : append self to viewer.objectsNeedingRedo
updateOwnGui=True : allow to update owngui at the end this func
"""
            v = kw.get( 'centers')
            mat = kw.has_key('materials')
            rad = kw.get( 'radii')

            redoFlags = apply( GLUSpheres.Set, (self, 0, 0), kw)

            if mat or rad or v:
                self.redoDspLst=1
                coords = self.vertexSet.vertices.array
                rad = self.vertexSet.radii.array
                mat = self.materials[GL.GL_FRONT].prop[1]
                self.initializeImposterRenderer(coords, rad, mat)

            return self.redoNow(redo, updateOwnGui, redoFlags)


        def initializeImposterRenderer(self, coords, rad, col):
            if __debug__:
             if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
            status = self.imposterRenderer.initRenderer()
            #print "initializeImposterRenderer status:", status
            if not status:
                print "Could not initialize the imposter renderer\n"
                return False
            self.imposterRenderer.initSubRenderers(0)
            self.imposterRenderer.clear()
            #brp = utimposterrend.BallRendererPtr(self.imposterRenderer.m_BallRenderer)

            #print "coords:", coords
            #print "radius:", rad
            #print "colors:", col
            if type(coords) == numpy.ndarray:
                coords = coords.tolist()
            if type(rad)== numpy.ndarray: 
                rad = rad.tolist()
            if type(col) == numpy.ndarray:
                col = col.tolist()
            #print "adding spheres...", len(rad)
            brp = self.imposterRenderer.m_BallRenderer
            if len(col) == len(coords):
                for i in range(len(coords)):
                    #print "adding ball", i
                    brp.addBall(coords[i][0], coords[i][1], coords[i][2], rad[i], col[i][0], col[i][1], col[i][2])
            elif len(rad) == len(coords):
                for i in range(len(coords)):
                    #print "adding ball", i
                    brp.addBall(coords[i][0], coords[i][1], coords[i][2], rad[i], col[0][0], col[0][1], col[0][2])
            else:
                for i in range(len(coords)):
                    #print "adding ball", i
                    brp.addBall(coords[i][0], coords[i][1], coords[i][2], 1., col[0][0], col[0][1], col[0][2])
            #print "done"    
                                                             
            return True


        def Draw(self):
            if __debug__:
             if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
            self.imposterRenderer.renderBuffer(True, False, None, None, 0, False, 1.0)

Spheres = GLUSpheres  # UTSpheres GLUSpheres


class TriangulateIcos:
    """Base class to compute vertices, faces and normals of a sphere based
    on icosahedral subdivision. Subclassed will implement different
    subdivision methods.
    A quality can be passed to the constructur which will trigger the
    precomputation of spheres of quality 0 to quality.
    To access the data, use getVertices(quality=val), getFaces(quality=val),
    getVNormals(quality=val) where val is the quality level """
    

    def __init__(self, quality=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
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
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        lenF = len(faces)
        lenV = len(vertices)

        for i in xrange(lenF):
            v0 = vertices[faces[i][0]]
            v1 = vertices[faces[i][1]]
            v2 = vertices[faces[i][2]]
            self.subdivideVert(v0, v1, v2)

        self.subdivideFaces(faces, lenV, quality)
    

    def normalize(self, v):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        d = math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
        if d == 0.0:
            print 'Zero length vector!'
            return
        return [v[0] / d, v[1] / d, v[2] / d]


    def getVertices(self, quality=0):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """ has to be implemented by subclass """
        pass

    def getVNormals(self, quality=0):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """ has to be implemented by subclass """
        pass
    
    def getFaces(self, quality=0):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        return getattr(self, 'facesQ%d'%quality)


class TriangulateIcosByEdgeCenterPoint(TriangulateIcos):
    
    def __init__(self, quality=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        #print "TriangulateIcosByEdgeCenterPoint.__init__"
        TriangulateIcos.__init__(self, quality)
    
        if self.quality > 0:
            for qu in range(1, self.quality+1):
                self.subsample(self.vertices,
                               getattr(self, 'facesQ%d'%(qu-1,) ),
                               qu)


    def subdivideVert(self, v0, v1, v2):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        #print "TriangulateIcosByEdgeCenterPoint.subdivideVert"
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
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        #print "TriangulateIcosByEdgeCenterPoint.subdivideFaces"
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
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        #print "TriangulateIcosByEdgeCenterPoint.getVertices"
        # the vertex list is very big, since vertices are added to this
        # list after every subsampling. Thus, only what is needed is returned
        v = 12 # vertices of icosahedron
        f = 20 # faces of icosahedron
        
        for i in range(quality):
            v = v+f*3
            f = f*4
        return self.vertices[:v]


    def getVNormals(self, quality=0):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        #print "TriangulateIcosByEdgeCenterPoint.getVNormals"
        # normals == vertices
        self.normals = self.getVertices(quality=quality)[:]
        return self.normals


class TriangulateIcosByFaceCenterPoint(TriangulateIcos):
    """ This class subdivides each face in 3 new faces by putting a center
    in the middle of a face triangle."""

    def __init__(self, quality=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        TriangulateIcos.__init__(self, quality)
        
        if self.quality > 0:
            for qu in range(1, self.quality+1):
                self.subsample(self.vertices,
                               getattr(self, 'facesQ%d'%(qu-1,) ),
                               qu)


    def subdivideVert(self, v0, v1, v2):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        # called by subsample
        v012 = []
                
        for i in range(3):
            v012.append(v0[i] + v1[i] + v2[i])

        v012 = self.normalize(v012)
        self.vertices.append(v012)


    def subdivideFaces(self, faces, lenV, quality):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
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
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        # the vertex list is very big, since vertices are added to this
        # list after every subsampling. Thus, only what is needed is returned
        v = 12 # vertices of icosahedron
        f = 20 # faces of icosahedron
        
        for i in range(quality):
            v = v+f
            f = f*3
        return self.vertices[:v]


    def getVNormals(self, quality=0):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        # normals == vertices
        return self.getVertices(quality=quality)

