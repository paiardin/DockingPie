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

#pmv -i -d ms 1crn.pdb
#execfile('testcsg3.py')

#pmv -i -m ms 1crn.pdb cv.pdb blurcv1crn_net.py
#pmv -i -m cpk 1crn.pdb cv.pdb
#pmv -i -m ms 1crn.pdb cv.pdb 
#execfile('testcsg1.py')
#execfile('testcsg2.py')

## TODO
##  try  using IndexedGeomDSPL to build self.dpyListCSG
##  OpenCSGGeom.clearPrimitives leaks memory
##  make PythonPrimitive::render() handle the Python global lock
##  add support for instanceMatrices

from OpenCSG import opencsglib as OpenCSG
import numpy
## from DejaVu2.Spheres import Spheres

class DejaVu2Primitive(OpenCSG.PythonPrimitive):

    def __init__(self, geom):
        apply( OpenCSG.PythonPrimitive.__init__, (self, self.render, OpenCSG.Intersection, 0))
        # does not work for some reason
        #OpenCSG.PythonPrimitive(self.render, OpenCSG.Intersection, 1)
        self.geom = geom
        self.dpyListCSG = None

            
    def redoDisplayListCSG(self):
        if self.dpyListCSG is not None:
            GL.glDeleteLists(1, self.dpyListCSG)
        
        g = self.geom
        self.dpyListCSG = GL.glGenLists(1)
        GL.glNewList(self.dpyListCSG, GL.GL_COMPILE)
##         if isinstance(g, Spheres):
##             g.DisplayFunction()
##         else:
        self.drawpolygons()
        GL.glEndList()


    def drawpolygons(self):
        g = self.geom
        vertices = g.getVertices()
        faces = g.getFaces()
        normals = g.getFNormals()
        GL.glDisable(GL.GL_CULL_FACE)
        for i,f in enumerate(faces):
            GL.glBegin(GL.GL_POLYGON)
            GL.glNormal3fv(normals[i])
            for vi in f:
                GL.glVertex3fv(vertices[vi])
            GL.glEnd()
            i+=1

            
    def render(self, mode='render'):
        # call with mode='csg' to render simple shape to setup Zbuffer for CSG
        # call with mode='render' to render by calling geom's draw function
        if self.geom:
            #import traceback
            #print traceback.print_stack()
            #print self.geom
            #print "========================================================="
            root = self.geom.viewer.rootObject

            instance = [0]
            p = self.geom.parent
            while p:
                instance.append(0)
                p = p.parent

            #mat = numpy.array(GL.glGetDoublev(GL.GL_MODELVIEW_MATRIX)).astype('f')
            #print 'mat OK', mat
            GL.glPushMatrix()
            GL.glLoadIdentity()
            self.geom.viewer.currentCamera.BuildTransformation()
            self.geom.BuildMat(self.geom, root, True, instance)

            #mat = numpy.array(GL.glGetDoublev(GL.GL_MODELVIEW_MATRIX)).astype('f')
            #print 'mat PB', mat
            #print 'render ', mode, self.geom
            if mode=='csg':
                if self.dpyListCSG is None:
                    self.redoDisplayListCSG()
                GL.glCallList(self.dpyListCSG)
            elif mode=='render':
                obj = self.geom
                if not obj.inheritMaterial:
                    obj.InitMaterial(0)
                    obj.InitColor(0)
                obj.DisplayFunction()
                
            GL.glPopMatrix()



from DejaVu2.Geom import Geom
from opengltk.OpenGL import GL
from DejaVu2.viewerFns import checkKeywords

class OpenCSGGeom(Geom):

    keywords = Geom.keywords + [
        'primitives',
        'algo',
        'depthalgo',
        ]
    def __init__(self, name=None, check=1, **kw):

        if __debug__:
            if check:
                apply( checkKeywords, (name,self.keywords), kw)
        apply( Geom.__init__, (self, name, 0), kw )
        self.primitives = OpenCSG.PrimitiveVector() # C++ primitives
        self.pyprimitives = [] # python subclasses used to call python implementation or render
        self.algo = OpenCSG.Goldfeather
        self.depthalgo = OpenCSG.DepthComplexitySampling

        self.Set(culling='none', algo=OpenCSG.Goldfeather,
                 depthalgo=OpenCSG.DepthComplexitySampling)
        
        
    def clearPrimitives(self):
##         for p in self.pyprimitives:
##             if p.dpyListCSG:
##                 print 'AAAAAAAAAAAAAAA', p.geom
##                 print 'AAAAAAAAAAAAAAA', p.dpyListCSG, p.geom.dpyList
##                 GL.glDeleteLists(1, p.dpyListCSG )
##                 p.dpyListCSG = None
        self.primitives.clear()
        self.pyprimitives = []


    def setPrimitives(self, *args):
        self.clearPrimitives()

        for g in args:
            assert isinstance(g, Geom)
            prim = DejaVu2Primitive(g)
            #self.primitives.append(prim)
            OpenCSG.PrimitiveVector_add(self.primitives, prim)
            self.pyprimitives.append(prim)


    def Set(self, check=1, redo=1, **kw):
	"""Set primitives"""

        if __debug__:
            if check:
                apply( checkKeywords, (self.name,self.keywords), kw)

        apply( Geom.Set, (self, 0, 0), kw)

        p = kw.get( 'primitives')
	if p:
            assert isinstance(p, OpenCSG.PythonPrimitiveVector)
            self.primitives = p

        a = kw.get( 'algo')
	if a:
            if a =='automatic':
                a = OpenCSG.Automatic
            elif a== 'goldfeather':
                a = OpenCSG.Goldfeather
            elif a == 'scs':
                a = OpenCSG.SCS
            assert a in (OpenCSG.Automatic, OpenCSG.Goldfeather, OpenCSG.SCS)
            self.algo = a

        d = kw.get( 'depthalgo')
	if d:
            if d =='DepthComplexitySampling':
                d = OpenCSG.DepthComplexitySampling
            elif d== 'NoDepthComplexitySampling':
                d = OpenCSG.NoDepthComplexitySampling
            elif d == 'OcclusionQuery':
                d = OpenCSG.OcclusionQuery
            assert d in (OpenCSG.DepthComplexitySampling,
                         OpenCSG.NoDepthComplexitySampling,
                         OpenCSG.OcclusionQuery)
            self.depthalgo = d


    def Draw(self):
        GL.glEnable(GL.GL_DEPTH_TEST);
        GL.glClear( GL.GL_STENCIL_BUFFER_BIT)
        GL.glDisable(GL.GL_FOG)
        #GL.glPolygonMode(GL.GL_FRONT_AND_BACK, GL.GL_FILL)
        
        OpenCSG.render(self.primitives, self.algo, self.depthalgo)

        GL.glDepthFunc(GL.GL_EQUAL)

        # FIXME should only enable fog if it is on in camera
        GL.glEnable(GL.GL_FOG)
        self.SetupGL()
        for p in self.pyprimitives:
            p.render()

        GL.glDepthFunc(GL.GL_LESS);


o=0
a = 20
vertices = [ (o, o, o), (o, a, o), (a, a, o), (a, o, o),
             (o, o, a), (o, a, a), (a, a, a), (a, o, a) ]
facesInward = [ (3,2,1,0), (7,6,2,3), (5,6,7,4), (0,1,5,4), (2,6,5,1), (4,7,3,0)]
facesOutward = [ (0,1,2,3), (3,2,6,7), (4,7,6,5), (4,5,1,0), (1,5,6,2), (0,3,7,4)]
faces = facesOutward
normals = [(0,0,1), (-1,0,0), (0,0,-1), (1,0,0), (0,-1,0), (0,1,0)]
from DejaVu2.IndexedPolygons import IndexedPolygons
clipBox = IndexedPolygons('clipBox', vertices=vertices, faces=faces, fnormals=normals,
                          frontPolyMode='line', shading='flat')
self.readMolecule('/home/annao/python/dev23/1crn.pdb', ask=0, parser=None, log=0)
#self.readMolecule('/home/annao/python/dev23/cv.pdb', ask=0, parser=None, log=0)
self.browseCommands('msmsCommands', commands=None, log=0, package='Pmv')
#self.computeMSMS("1crn;cv", 'MSMS-MOL', perMol=1, density=4.6, log=0, pRadius=1.5)
self.computeMSMS("1crn", 'MSMS-MOL', perMol=1, density=4.6, log=0, pRadius=1.5)
srf1 = self.Mols[0].geomContainer.geoms['MSMS-MOL']
#srf2 = self.Mols[1].geomContainer.geoms['MSMS-MOL']

cpk1 = self.Mols[0].geomContainer.geoms['cpk']
#cpk2 = self.Mols[1].geomContainer.geoms['cpk']

self.colorByResidueType("1crn;cv", ['MSMS-MOL'], log=0)
#self.displayMSMS("1crn;cv", negate=True, only=False, surfName=['MSMS-MOL'], log=0, nbVert=1)
#self.displayCPK("1crn;cv", cpkRad=0.0, scaleFactor=1.0, only=False, negate=False, quality=17)
#self.showMolecules(['cv'], negate=True, log=0)
#self.colorByAtomType("1crn;cv", ['cpk'], log=0)

srf1.Set(frontPolyMode='line', visible=False)
#srf2.Set(frontPolyMode='line')

err = OpenCSG.glewInit()
print "glewInit status: ", err

g = OpenCSGGeom('inter')

g.setPrimitives(srf1, clipBox)
#g.setPrimitives(srf1, srf2)
#g.setPrimitives(cpk1, clipBox)
g.Set(immediateRendering=True, lighting=True)

self.GUI.VIEWER.AddObject(clipBox)
self.GUI.VIEWER.AddObject(g)
#g.primitives[1].setOperation(OpenCSG.Subtraction)
