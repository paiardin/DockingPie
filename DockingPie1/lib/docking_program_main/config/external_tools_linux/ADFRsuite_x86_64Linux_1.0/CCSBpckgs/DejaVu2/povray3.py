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
# $Header: /mnt/raid/services/cvs/DejaVu2/povray3.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#
# $Id: povray3.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#

"""
Povray module: driver to generate PovRay scenes form a DejaVu2 scene

The driver currently handles: camera parameters (background color, perpective parameters, position), light sources, and Spheres, IndexedPolygons and IndexedPolylines geometries.
The projection is still somewhat approximative ! still a 10 translation that is not explained :(.
"""

from opengltk.OpenGL import GL
#from opengltk.extent.utillib import glTriangleNormals
from geomutils.geomalgorithms import  TriangleNormals
from DejaVu2 import viewerConst
import numpy

class PovRay:
    """Driver for Povray v3.x
    """
    
    def __init__(self, includes = []):
        """Create a PovRay file generator class for a given camera"""
        self.camera = None
        self.entries = ["""
//
// Povray scene file written by Pmv version 0.1
//
// POV-Ray can be retrieved at: http://www.povray.org
// N.Guex, 1995-1999
//
"""]
        for i in includes:
            self.entries.append('#include "%s"\n'%i)


    def write(self, filename):
        """void <- write(filename)  dumps povray scene description"""
        self.filename = filename
        f = open(filename, 'w')
        for entry in self.entries:
            f.write(entry)
        f.close()


    def clear(self):
        """void <- clear()  empties povray scene description"""
        self.entries = []


    def coord3(self, x, y=None, z=None):
        """string <- coord3(x, y=None, z=None) returns a povray vector string
        this function inverts the z coordinates.
        It can be called with a sequence of 3 or with 3 values.
        """
        if y is None:
            return ' <%f, %f, %f> '% (x[0],x[1],-x[2])
        else:
            return ' <%f, %f, %f> '% (x, y, -z)


    def vec3(self, x, y=None, z=None):
        """string <- vec3(x, y=None, z=None) returns a povray vector string
        this function does not invert the z coordinates.
        It can be called with a sequence of 3 or with 3 values.
        """
        if y is None:
            return ' <%.2f, %.2f, %.2f> '% (x[0],x[1],x[2])
        else:
            return ' <%.2f, %.2f, %.2f> '% (x, y, z)


    def vec4(self, x, y=None, z=None, t=None):
        """string <- vec4(x, y=None, z=None, t=None) returns a povray 4-vector
        this function does not invert the z coordinates.
        It can be called with a sequence of 3 and a y-value or with 4 values.
        """
        if z is None:
            return ' <%.2f, %.2f, %.2f, %.2f> '% (x[0],x[1],x[2],y)
        else:
            return ' <%.2f, %.2f, %.2f, %.2f> '% (x, y, z, t)


    def matrix(self, mat):
        """string <- matrix(mat) returns a 4x4 matrix as a povray matrix"""
        
        str = 'matrix < %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f >\n' % (
            mat[0][0], mat[0][1], mat[0][2],
            mat[1][0], mat[1][1], mat[1][2],
            mat[2][0], mat[2][1], mat[2][2],
            mat[0][3], mat[1][3], mat[2][3] )
        return str


    def addCamera(self, camera, scaleLight=2.0):
        """void <- addCamera(camera)  adds the given camera to the povray scene
        handles: camera background color
                 projection matrix (only perspective currently
                 light sources
        """
                 
        # doesn't handle camera transformation yet
        str = 'background { color rgb ' + self.vec3(camera.backgroundColor)
        str = str + '}\n\n'
        self.entries.append( str )

        str = 'camera { perspective'
        # I can't remember why I have to move 10 more than direction vector
        # but it seems familiar
        str = str + '\n  location' + self.coord3(camera.lookFrom[0],
                                                 camera.lookFrom[1],
                                                 camera.lookFrom[2]+10)
        str = str + '\n  look_at' + self.coord3(camera.lookAt)
        str = str + '\n  angle %f'% ( camera.fovy+10.0, )
        str = str + '\n  up <0,1,0> // required for 1/1 aspect ratio p 277'
        str = str + '\n  right <%d/%d,0,0>'%(camera.width,camera.height)
        str = str + '\n}\n\n'
        self.entries.append( str )
        for l in camera.viewer.lights:
            if l.enabled: self.addLight(l, scaleLight)


    def addLight(self, light, scale=2.0):
        """void <- addLight(light) add a light source to the povray scene
        """
        # doesn't handle light transformation yet
        d = light.direction
        str = 'light_source {' + self.coord3(d[0], d[1], d[2])
#        str = str + 'color rgb' + self.vec3(light.diffuse) + '}\n\n'
        str = str + 'color rgb ' + self.vec3(scale*light.diffuse[0],
                                             scale*light.diffuse[1],
                                             scale*light.diffuse[2])
        str += 'parallel '
        str = str + '}\n\n'
        self.entries.append( str )


    def addTexture(self, texture, col):
        """void <-  addTexture(texture, col) Add texture to an object.
        texture is a dictionnary of texture properties,
        col is used as pigment is pimgment is not in textures
        """
        str = '  texture {\n'
        if not 'pigment' in texture.keys():
            str = str + '      pigment { color rgb %s }\n' %self.vec3( col )
        elif texture['pigment']=='':
            str = str + '      pigment { color rgb %s }\n' %self.vec3( col )
        else:
            str = str + '      pigment { %s }\n' % texture['pigment']
        for k,v in texture.items():
            if k=='pigment': continue
            if v is None:
                str = str + '      %s\n' % (k,)
            else:
                str = str + '      %s { %s }\n' % (k,v)
        str = str +  '  }\n'
        self.entries.append( str )


    def endShape(self):
        self.entries.append( '}\n' );

        
    def addGeoms(self, geometries,
                 texture = {'finish':'specular 1 roughness 0.001 ambient 0.3'},
                 bondRad = 0.15):
        """void <- addGeoms(geometries, texture) adds a list of geometries to
        the povray scene
        only visible geometries are added
        """
        for g in geometries:
            if g.visible and len(g.vertexSet):
                self.entries.append('// geometry %s\n//\n'%g.name)
                self.addGeom(g, texture, bondRad)


    def addGeom(self, geom,
                texture = {'finish':'specular 1 roughness 0.001 ambient 0.3'},
                bondRad = 0.15):
        """void <- addGeom(geometry, texture) adds a geometries to the povray
        scene
        Spheres, IndexedPolygons and IndexedPolylines are handled
        """
        from DejaVu2.Spheres import Spheres
        from DejaVu2.IndexedPolygons import IndexedPolygons
        from DejaVu2.IndexedPolylines import IndexedPolylines
        from DejaVu2.Cylinders import Cylinders
        
        if isinstance(geom, Spheres):
            self.entries.append('// Object %s\n//\n' % geom.name)
            self.addSpheres(geom, texture)
        elif isinstance(geom, IndexedPolygons) and len(geom.faceSet):
            self.entries.append('// Object %s\n//\n' % geom.name)
            self.addIndexedPolgygons(geom, texture)
        elif isinstance(geom, IndexedPolylines) and len(geom.faceSet):
            self.entries.append('// Object %s\n//\n' % geom.name)
            self.addIndexedPolylines(geom, texture, bondRad)
        elif isinstance(geom, Cylinders) and len(geom.faceSet):
            self.entries.append('// Object %s\n//\n' % geom.name)
            self.addCylinders(geom, texture)
        else:
            print 'WARNING: %s the geometry is not yet supported'%geom.__class__

    def addIndexedPolylines(self, geom, texture, bondRad):
        """void <- addIndexedPolylines(geom, texture)
        """
        mat = geom.GetMatrix()
        v = geom.vertexSet.vertices*mat
        c = geom.materials[GL.GL_FRONT].prop[1]
        lines = geom.faceSet.faces.array

        for i in xrange(len(v)):
            if len(c)==len(v): col = c[i]
            else: col = c[0]
            self.entries.append('sphere{%s,%f\n' % (self.coord3(v[i]),
                                                    bondRad) )
            self.addTexture( texture, col )
            self.endShape()

        for j in xrange(len(lines)):
            l = lines[j]
            if len(c) == len(v):
                col1 = c[l[0]]
                col2 = c[l[1]]
            else: col1 = col2 = c[0]
            
            if numpy.sum(col1-col2) < 0.0001:
                p2 = v[l[1]]
                oneCyl = 1
            else:
                p2 =  (v[l[1]]+v[l[0]])/2.
                oneCyl = 0

            self.entries.append('cylinder{%s,%s, %f open\n' % \
                                (self.coord3(v[l[0]]),
                                 self.coord3(p2),
                                 bondRad) )

            self.addTexture( texture, col1 )
            self.endShape()

            if not oneCyl:
                self.entries.append('cylinder{%s, %s, %f open\n' % \
                                    (self.coord3(p2),
                                     self.coord3(v[l[1]]),
                                     bondRad))
                self.addTexture( texture, col2 )
                self.endShape()


    def addIndexedPolgygons(self, geom, texture):
        """void <- addIndexedPolgygons(geom, texture)
        """
        mat = geom.GetMatrix()
        vt = geom.vertexSet.vertices*mat

        # FIXME need to add instance
        # FIXME need to handle flat shading
        # add vertices
        lines = ["mesh2 {\n\tvertex_vectors {\n\t\t%d,\n"%len(vt)]
        for v in vt:
            lines.append("\t\t%s,\n"%self.coord3(v))
        lines.append("\t}\n")

        # add normals
        normals = geom.getVNormals()
        lines.append("\tnormal_vectors {\n\t\t%d,\n"%len(normals))
        for n in normals:
            lines.append("\t\t%s,\n"%(self.coord3(n)))
        lines.append("\t}\n")

        # handle colors
        colors = geom.materials[GL.GL_FRONT].prop[1]
            
        colBinding = 'Overall'
        faces = geom.getFaces()
        if len(colors)==len(vt):
            colBinding = 'Per Vertex'
        elif len(colors)==len(faces):
            colBinding = 'Per Face'

        print len(colors), len(faces), len(vt), colBinding
        if colBinding!='Overall':
            lines.append("\ttexture_list {\n\t\t%d,\n"%len(colors))
            for c in colors:
                lines.append(
                    "\t\ttexture { pigment { color rgb<%6.3f, %6.3f, %6.3f> }\n"% tuple(c[:3]))
                lines.append("\t\t\tfinish { specular 1 roughness 0.001 ambient 0.3 } }\n")
            lines.append("\t}\n")
       
        # add faces
        lines.append("\tface_indices {\n\t\t%d,\n"%len(faces))
        faceNumberLine = len(lines)-1
        nbfaces = 0
        if colBinding=='Overall':
            for t in faces:
                for ti in range(len(t)-2):
                    lines.append("\t\t<%d,%d,%d>\n"%(t[0], t[ti+1], t[ti+2]))
                    nbfaces += 1
                    
        elif colBinding=='Per Face':
            for fn,t in enumerate(faces):
                for ti in range(len(t)-2):
                    lines.append("\t\t<%d,%d,%d>,%d\n"%(
                        t[0], t[ti+1], t[ti+2], fn))
                    nbfaces += 1

        elif colBinding=='Per Vertex':
            for t in faces:
                for ti in range(len(t)-2):
                    lines.append("\t\t<%d,%d,%d>,%d,%d,%d\n"%(
                        t[0], t[ti+1], t[ti+2], t[0], t[ti+1], t[ti+2]))
                    nbfaces += 1
        lines.append("\t}\n")

        lines[faceNumberLine] = "\tface_indices {\n\t\t%d,\n"%nbfaces
        
        if colBinding=='Overall':
            lines.append("\tpigment { color rgb<%6.3f, %6.3f, %6.3f> }\n"% tuple(colors[0][:3]))
            lines.append("\t\tfinish { specular 1 roughness 0.001 ambient 0.3 }\n")

        lines.append("}\n")

        self.entries.extend(lines)
        
##         mat = geom.GetMatrix()
##         v = geom.vertexSet.vertices*mat
##         #v = geom.vertexSet.vertices.array
##         c = geom.materials[GL.GL_FRONT].prop[1]
##         n = geom.normals
##         tri = geom.faceSet.faces.array

##         fn = TriangleNormals( v, tri, 'PER_FACE')
        
##         colBinding = 0 #Overall
##         if not color:
##             if len(c)==len(v): colBinding = 1 # Per Vertex
##             elif len(c)==len(tri): colBinding = 2 # Per Face
##             else:
##                 col = c[0] # Overall
##         else:
##             col = color

##         vnorm = n
##         for j in xrange(len(tri)):
##             t = tri[j]
##             if len(tri)==len(n):
##                 # per face normal
##                 pass

##             # tri to fix the normals of singular vertices
##             # if dot(vertex normal, face normal) < 0 negate vertex normal
## ##              vnorm = [0,0,0]
## ##              for i in (0,1,2):
## ##                  norm = n[t[0]]
## ##                  dot = norm[0]*fn[j][0] + norm[1]*fn[j][1] + norm[2]*fn[j][2]
## ##                  if dot < 0.0:
## ##                      vnorm[i] = 0.0 - norm
## ##                  else:
## ##                      vnorm[i] = norm
                    
##             if not color:
##                 if colBinding==1:
##                     col = (c[t[0]] + c[t[1]] + c[t[2]])*0.33333333
##                 elif colBinding==2:
##                     col = c[j]

            
##             self.entries.append('smooth_triangle{\n')
##             self.entries.append('  %s,%s\n'% (self.coord3(v[t[0]]),
##                                               self.coord3(vnorm[t[0]])) )
##             self.entries.append('  %s,%s\n'% (self.coord3(v[t[1]]),
##                                               self.coord3(vnorm[t[1]])) )
##             self.entries.append('  %s,%s\n'% (self.coord3(v[t[2]]),
##                                               self.coord3(vnorm[t[2]])) )
##             self.addTexture( texture, col )
##             self.endShape()

##             if len(t)==4:
##                 self.entries.append('smooth_triangle{\n')
##                 self.entries.append('  %s,%s\n'% (self.coord3(v[t[0]]),
##                                                   self.coord3(vnorm[t[0]])) )
##                 self.entries.append('  %s,%s\n'% (self.coord3(v[t[2]]),
##                                                   self.coord3(vnorm[t[2]])) )
##                 self.entries.append('  %s,%s\n'% (self.coord3(v[t[3]]),
##                                                   self.coord3(vnorm[t[3]])) )
##                 self.addTexture( texture, col )
##                 self.endShape()


    def addSpheres(self, geom, texture):
        """void <- addSpheres(geom, texture)
        """
        mat = geom.GetMatrix()
        v = geom.vertexSet.vertices*mat
        #v = geom.vertexSet.vertices.array
        c = geom.materials[1028].prop[1]
        if geom.oneRadius == viewerConst.NO:
            radii = geom.vertexSet.radii.array

        if geom.inheritMaterial:
            fp = None # FIXME
            bp = None
        else:
            fp = geom.materials[GL.GL_FRONT]
            if not geom.frontAndBack:
                bp = geom.materials[GL.GL_BACK]
                face = GL.GL_FRONT
            else:
                bp = None
                face = GL.GL_FRONT_AND_BACK
        
        self.entries.append('#declare unitSphere = sphere { <0.00, 0.00, 0.00> , 1 }\n\n')
        for i in range(len(geom.vertexSet)):
            self.entries.append('object { unitSphere\n')
            if not geom.oneRadius:
                self.entries.append('  scale %f\n'%(radii[i][0]/2.) )
            else:
                self.entries.append('  scale %f\n'%geom.radius )
            self.entries.append('  translate %s\n'%self.coord3(v[i]) )
            if fp:
                col = fp.prop[1][i][:3]
            else:
                col = (1.,1.,1.)
            self.addTexture( texture, col )
            self.endShape()


    def addCylinders(self, geom, texture):
        """void <- addSpheres(geom, texture)
        """
        mat = geom.GetMatrix()
        v = geom.vertexSet.vertices*mat
        #v = geom.vertexSet.vertices.array
        c = geom.materials[1028].prop[1]
        radii = geom.vertexSet.radii.array

        if geom.inheritMaterial:
            fp = None # FIXME
            bp = None
        else:
            fp = geom.materials[GL.GL_FRONT]
            if not geom.frontAndBack:
                bp = geom.materials[GL.GL_BACK]
                face = GL.GL_FRONT
            else:
                bp = None
                face = GL.GL_FRONT_AND_BACK

        f = geom.getFaces()
        for i in range(len(geom.faceSet)):
            for j in range(len(f[i])-1):
                self.entries.append('cone { \n')
                if not geom.oneRadius:
                    r1 = radii[f[i][j]]/2.
                    r2 = radii[f[i][j+1]]/2.
                else:
                    r1 = r2 = radii[0]/2.

                self.entries.append('\t%s, %f, %s, %f\n'%(
                     self.coord3(v[f[i][j]]), r1,
                     self.coord3(v[f[i][j+1]]), r2))

                if fp:
                    col = c[f[i][j]][:3]
                else:
                    col = (1.,1.,1.)
                self.entries.append("\tpigment { color rgb<%6.3f, %6.3f, %6.3f> }\n"% tuple(col))
                self.entries.append("\tfinish { specular 1 roughness 0.001 ambient 0.3 }\n")
                self.endShape()


if __name__ == '__main__':
    from DejaVu2 import Viewer
    vi = Viewer()
    from DejaVu2.Spheres import Spheres
    s = Spheres('test', centers = ( (0.0, 0., 0.), (3.,0.,0.) ) )
    s.Set(quality = 4)
    vi .AddObject(s)
    from DejaVu2.povray3 import PovRay
    p = PovRay()
    p.addCamera(vi.cameras[0])
    p.addLight(vi.lights[0])
    p.addGeom(s)
    p.write('test.pov')
    
    
