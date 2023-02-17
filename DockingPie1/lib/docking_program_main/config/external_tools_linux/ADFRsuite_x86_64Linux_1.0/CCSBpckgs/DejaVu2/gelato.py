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
# adapted From povray3.py
#############################################################################

#
# $Header: /mnt/raid/services/cvs/DejaVu2/gelato.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#
# $Id: gelato.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#

"""
Povray module: driver to generate PovRay scenes form a DejaVu2 scene

The driver currently handles: camera parameters (background color, perpective parameters, position), light sources, and Spheres, IndexedPolygons and IndexedPolylines geometries.
The projection is still somewhat approximative ! still a 10 translation that is not explained :(.
"""
#from ARViewer import util
from opengltk.OpenGL import GL
#from opengltk.extent.utillib import glTriangleNormals
from geomutils.geomalgorithms import  TriangleNormals
from DejaVu2 import viewerConst
from mglutil.math import rotax
from numpy import matrix
import numpy


PROJ={	1:"orthographic",
	0:"perspective"
	}


class Shader:
    ShaderDic={"plastic" 	: [	"\"float Ka\", 0",
					"\"float Kd\", 1",
					"\"float Ks\", 0.75",
					"\"float roughness\", 0.05"],

		"glass" 	: [	"\"float Ka\", 0.2",
					"\"float Kd\", 0.25",
					"\"float Ks\", 0.4",
					"\"float roughness\", 0.1",
					"\"string envname\",\"reflection\"",	
					"\"float Kr\", 0.5",
					"\"float samples\", 4"],

		"clay" 		: [	"\"float Ka\", 0.9",
					"\"float Kd\", 1",
					"\"float roughness\", 0.2"],

		"shinyplastic"	: [	"\"float Ks\", 0.9",
					"\"float Kd\", 1",
					"\"float eta\", 1.5",
					"\"string envname\", \"reflection\""],

		"plasticss"	: [	"\"float Ks\", 0.9",
					"\"float Kd\", 1",
					"\"float Kss\", 1"],
		
		"metal" 	: [	"\"float Ka\", 1.0",
					"\"float Kd\", 0.5",
					"\"float Ks\", 0.8",
					"\"float roughness\", 0.2",
					"\"string envname\",\"reflection\"",
					"\"float Kr\", 0.5",
					"\"float samples\", 4"],
	
		"ambocclude"	: [	"\"string occlusionname\",\"localocclusion\"",
					"\"float samples\", 256",
					"\"float bias\", 0.01"],

		"greenmarble"	: [	"\"float Ka\", 0.1",
					"\"float Kd\", 0.6",
					"\"float Ks\", 0.4",
					"\"float roughness\", 0.1",
					"\"float veinfreq\", 1",
					"\"float sharpness\", 25",
					"\"float shadingfreq\", 1"],

		"oak"		: [	"\"float Ka\", 1",
					"\"float Kd\", 1",
					"\"float Ks\", 0.25",
					"\"float roughness\", 0.25",
					"\"float divotdepth\", 0.5",
					"\"float ringy\", 1",
					"\"float ringfreq\", 8",
					"\"float ringunevenness\", 0.3",
					"\"float ringnoise\", 0.02",
					"\"float ringnoisefreq\", 1",
					"\"float trunkwobble\", 0.15",
					"\"float trunkwobblefreq\", 0.025",
					"\"float angularwobble\", 0.15",
					"\"float angularwobblefreq\", 0.025",
					"\"float grainy\", 1",
					"\"float grainfreq\", 8",
					"\"float shadingfreq\", 1"],

		"screen"	: [	"\"float Ka\", 1.0",
					"\"float Kd\", 0.75",
					"\"float Ks\", 0.4",
					"\"float roughness\", 0.1",
					"\"float sfreq\", 10",
					"\"float tfreq\", 10",
					"\"float sdensity\", 10",
					"\"float tdensity\", 10"],

		"soapbubble"	: [	"\"float Ka\", 0.2",
					"\"float Kd\", 0.25",
					"\"float Ks\", 0.4",
					"\"float roughness\", 0.1",
					"\"string envname\", \"reflection\"",
					"\"float Kr\", 0.5",
					"\"float samples\", 4",
					"\"float fakerefract\", 1"],

		"shownormals"	: [	"\"float bias\", 1"],

		"showfacing"	: [	"\"float Ka\", 0.25",
					"\"float Kd\", 0.75"]
	   }

    def __init__(self):
	self.name=""
	self.Type="surface"
	self.parameters=[]
	

    def vec3(self, x, y=None, z=None):
        """string <- vec3(x, y=None, z=None) returns a povray vector string
        this function does not invert the z coordinates.
        It can be called with a sequence of 3 or with 3 values.
        """
        if y is None:
            return ' (%.2f, %.2f, %.2f) '% (x[0],x[1],x[2])
        else:
            return ' (%.2f, %.2f, %.2f) '% (x, y, z)
    
    def addParameters(typ,name,value):
	if typ == "float" : self.parameters.append('"float %s", %f' % (typ,name,value))
	elif typ == "string" : self.parameters.append('"string %s", %s' % (typ,name,value))
	elif typ == "color" : self.parameters.append('"color %s", %s' % (typ,name,self.vec3(value)))
	
class Gelato:
    """Driver for Povray v3.x
    """
    
    def __init__(self, shaderD=None,includes = []):
        """Create a PovRay file generator class for a given camera"""
        self.camera = None
	self.PRECISION = 6
	if shaderD== None : self.shaderDic=Shader.ShaderDic
	else : self.shaderDic=shaderD
        self.entries = ["""
# 
# Molecular graphics export from PMV 1.5.4
# Requires NVIDIA Gelato 2.1, PYG format
#
###Some usefull Function######
####nurbscyl comes from the Technical REference guide of gelato
####create a cylinder
def nurbscyl():
    uknot = (0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4)
    vknot = (0, 0, 1, 1)
    Pw = ( 1, 0, 0, 1, 1, 1, 0, 1, 0, 2, 0, 2, -1, 1, 0, 1, -1, 0, 0, 1,
          -1, -1, 0, 1, 0, -2, 0, 2, 1, -1, 0, 1, 1, 0, 0, 1, 1, 0, -3, 1, 1,
          1, -3, 1, 0, 2, -6, 2, -1, 1, -3, 1, -1, 0, -3, 1, -1, -1, -3, 1,
          0, -2, -6, 2, 1, -1, -3, 1, 1, 0, -3, 1 )
    Scale(0.2,0.2,0.3)
    Patch (9, 3, uknot, 0, 4, 2, 2, vknot, 0, 1, "vertex hpoint Pw", Pw)
def vmd_cylinder():
    uknot = (0,0,0,0.25,0.25,0.5,0.5,0.75,0.75,1,1,1,)
    vknot = (0, 0, 1, 1)
    Pw = (0.0429226,0,0,1,0.0303509,0.0303509,0,0.707107,
-1.87621e-09,0.0429226,0,1,-0.0303509,0.0303509,
0,0.707107,-0.0429226,-3.75241e-09,0,1,-0.0303509,-0.0303509,
0,0.707107,5.62862e-09,-0.0429226,0,1,0.0303509,-0.0303509,
0,0.707107,0.0429226,2.38805e-09,0,1,0.0429226,0,0.110961,1,
0.0303509,0.0303509,0.0784616,0.707107,-1.87621e-09,0.0429226,
0.110961,1,-0.0303509,0.0303509,0.0784616,0.707107,-0.0429226,
-3.75241e-09,0.110961,1,-0.0303509,-0.0303509,0.0784616,0.707107,
5.62862e-09,-0.0429226,0.110961,1,0.0303509,-0.0303509,0.0784616,
0.707107,0.0429226,2.38805e-09,0.110961,1,)
    Patch (9, 3, uknot, 0, 1, 2, 2, vknot, 0, 1, "vertex hpoint Pw", Pw)
"""]

    def write(self, filename):
        """void <- write(filename)  dumps gelato scene description"""
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
            return ' (%f, %f, %f) '% (x[0],x[1],(-x[2]))
        else:
            return ' (%f, %f, %f) '% (x, y, (-z))

    def vec2(self, x, y=None):
        """string <- vec3(x, y=None, z=None) returns a povray vector string
        this function does not invert the z coordinates.
        It can be called with a sequence of 3 or with 3 values.
        """
        if y is None:
            return ' (%.2f, %.2f) '% (x[0],x[1])
        else:
            return ' (%.2f, %.2f) '% (x, y)

    def vec3(self, x, y=None, z=None):
        """string <- vec3(x, y=None, z=None) returns a povray vector string
        this function does not invert the z coordinates.
        It can be called with a sequence of 3 or with 3 values.
        """
        if y is None:
            return ' (%.2f, %.2f, %.2f) '% (x[0],x[1],x[2])
        else:
            return ' (%.2f, %.2f, %.2f) '% (x, y, z)


    def vec4(self, x, y=None, z=None, t=None):
        """string <- vec4(x, y=None, z=None, t=None) returns a povray 4-vector
        this function does not invert the z coordinates.
        It can be called with a sequence of 3 and a y-value or with 4 values.
        """
        if z is None:
            return ' (%.2f, %.2f, %.2f, %.2f) '% (x[0],x[1],x[2],y)
        else:
            return ' (%.2f, %.2f, %.2f, %.2f) '% (x, y, z, t)


    def matrix16(self, mat):
        """string <- matrix(mat) returns a 16, matrix"""
        
        str = '( %f, %f, %f,%f, %f, %f, %f,%f, %f, %f, %f,%f, %f, %f, %f ,%f)' % (
            mat[0], mat[1], mat[2],mat[3],
            mat[4], mat[5], mat[6],mat[7],
            mat[8], mat[9], mat[10],mat[11],
            mat[12], mat[13], mat[14]*-1.,mat[15])
        return str

    def set_transform(self, matrix):
                return ('SetTransform %s\n' % self.matrix16(matrix))
		    
    def append_transform(self, matrix):
                return ('AppendTransform %s\n' % self.matrix16(matrix))

    def addRender(self):
	self.entries.append("Render()")

    def addCamera(self, camera, scaleLight=2.0):
        """void <- addCamera(camera)  adds the given camera to the povray scene
        handles: camera background color
                 projection matrix (only perspective currently
                 light sources
        """
                 
        # doesn't handle camera transformation yet
	r,t,s=camera.viewer.rootObject.Decompose4x4(camera.GetMatrix())
	astr=''
	astr = astr + '\nAttribute(\"int[2] resolution\", ('+str(int(camera.width))+','+str(int(camera.height))+'))'
	#astr = astr + '\nAttribute (\"int[2] spatialquality\", (4, 4))'
	proj= PROJ[camera.projectionType]
	astr = astr + '\nAttribute(\"string projection\",  \"'+proj+'\")'
	if (camera.projectionType == 1):#ORTHO
		astr = astr + '\nAttribute ("float[4] screen",'+self.vec4(camera.left,camera.right,camera.bottom,camera.top)+')'
		#astr = astr + ('\nAttribute("float[4] screen", (%f, %f, %f, %f))'%(camera.left,camera.right,camera.bottom,camera.top)) #have to be adjust
		astr = astr + '\nAttribute ("float far", '+str(camera.far)+')'
		astr = astr + '\nAttribute ("float near", '+str(camera.near)+')'
		#astr = astr + '\nTranslate'+self.vec3(camera.translation)
	else :
		#astr = astr + '\nAttribute ( "float shadingquality", 0.25 )'
		#astr = astr + '\nAttribute ( "int[2] spatialquality", (1, 1) )'
		astr = astr + '\nAttribute ("float fov", '+str(camera.fovy)+')'
	#str = str + '\nCamera(\"'+camera.name+'\")'
	astr = astr + '\nTranslate '+self.coord3(t)
        self.entries.append( astr )
	self.entries.append("\nWorld()\n")
	#self.addBackground(camera)
	#self.entries.append("Shader (\"surface\", \"plastic\", \"float Ks\", 0.9, \"float Kd\", 1)\n")
        for l in camera.viewer.lights:
            if l.enabled: self.addLight(l, scaleLight)

    def addBackground(self,camera):
	astr="#####Background#######\n"
	astr = astr + "PushAttributes()\n"
	astr = astr + "Shader(\"surface\", \"constant\")\n"
	astr = astr + ("Attribute(\"color C\", %s)\n"%self.vec3(camera.backgroundColor))
	astr = astr + "Input(\"backplane.pyg\")\n"
	astr = astr + "PopAttributes()\n"
	self.entries.append( astr )

    def addLight1(self):
	self.entries.append("""Light ("amb1", "ambientlight", "float intensity", 0.1)
#Light ("spt1", "spotlight", "float intensity", 10.0,"point from",(0,0,20),"point to",(-0.628775, 10.648891, 26.800826))
Light ("dpt1", "distantlight", "float intensity",0.8,"color lightcolor" , (1,1,1))
#PushTransform ()
#Translate (-0.628775, -2.648891, 20.800826)
#Light ("pt2", "pointlight", "float intensity", 10.0)
#PopTransform ()
#PushTransform ()
#Translate (0, 1, 20)
#Light ("pt1", "pointlight", "float intensity", 1.0)
#PopTransform ()
""")
	

    def addLight(self, light, scale=2.0):
        """void <- addLight(light) add a light source to the povray scene
        """
	gv=[0.,0.,1]
	v=light.direction[0:3]
	mat=numpy.array(rotax.rotVectToVect(gv,v))
	astr="PushTransform ()\n"
	astr= astr + ("AppendTransform %s\n"%self.matrix16(mat.reshape(16,)))
	astr= astr + "Light (\""+light.name.replace(" ","")+"\", \"distantlight\",\"float intensity\",0.5,"
	astr= astr + "\"float __nonspecular\","+str(1-light.specular[0])+","
	astr= astr + ("\"color lightcolor\" , %s)\n" % self.vec3(scale*light.diffuse[0],scale*light.diffuse[1], scale*light.diffuse[2]))

	astr= astr + "PopTransform()\n" 
        self.entries.append( astr )

    def addShader(self,shader):
	if shader == "soapbubble" : astr='Shader("surface","glass"'
	else : astr='Shader("surface","%s"' % shader
	for i in self.shaderDic[shader]:
		astr = astr +","+ i 
	astr = astr + ")\n"
	self.entries.append(astr)
	if shader == 'shinyplastic' or shader == 'metal' or shader == 'soapbubble' or shader == 'glass' :
		self.entries.append("Attribute (\"string geometryset\", \"+reflection\")\n")
	if shader == 'ambocclude' :
		self.entries.append("""Attribute ("string geometryset", "+localocclusion")
Attribute ("float occlusion:maxpixeldist", 20)
Attribute ("float occlusion:maxerror", 0.25)
""")
	
    def addShader_old(self,shader):
	if shader == '' : self.entries.append("")
	if shader == 'clay' : 
		self.entries.append("Shader (\"surface\", \"clay\", \"float Ka\", 0.9, \"float Kd\", 1,\"float roughness\", 0.2)\n")
	if shader == 'plastic' :
		#self.entries.append("Shader (\"surface\", \"plastic\", \"float Ks\", 0.9, \"float Kd\", 1)\n")
		self.entries.append("Shader(\"surface\", \"plastic\", \"float Ka\", 0, \"float Kd\", 1, \"float Ks\", 0.75, \"float roughness\", 0.05)\n")
	elif shader == 'shinyplastic' :
		self.entries.append("Attribute (\"string geometryset\", \"+reflection\")\n")
		self.entries.append("Shader (\"surface\", \"shinyplastic\", \"float Ks\", 0.9, \"float Kd\", 1,\"float eta\", 1.5,\"string envname\", \"reflection\")\n")
	elif shader == 'plasticss':
		self.entries.append("Shader (\"surface\", \"plasticss\", \"float Ks\", 0.9, \"float Kd\", 1,\"float Kss\", 1)\n")
	elif shader == 'metal' :	
		self.entries.append("Attribute (\"string geometryset\", \"+reflection\")\n")
		self.entries.append("Shader (\"surface\", \"metal\", \"float Ka\", 1.0, \"float Kd\", 0.5,\"float Ks\", 0.8,\"float roughness\", 0.2,\"string envname\", \"reflection\",\"float Kr\", 0.5,\"float samples\", 4)\n")
	elif shader == 'occlusion' :
		self.entries.append("""Shader ("surface", "ambocclude", "string occlusionname", "localocclusion", "float samples", 256, "float bias", 0.01)
Attribute ("string geometryset", "+localocclusion")
Attribute ("float occlusion:maxpixeldist", 20)
Attribute ("float occlusion:maxerror", 0.25)
""")
	elif shader == 'screen' :
		self.entries.append("Shader (\"surface\", \"screen\", \"float Ka\", 1.0, \"float Kd\", 0.75,\"float Ks\", 0.4,\"float roughness\", 0.1, \"float sfreq\", 10, \"float tfreq\", 10, \"float sdensity\", 10, \"float tdensity\", 10)\n")
	elif shader == 'glass' :
		self.entries.append("Attribute (\"string geometryset\", \"+reflection\")\n")
		self.entries.append("Shader (\"surface\", \"glass\", \"float Ka\", 0.2, \"float Kd\", 0.25,\"float Ks\", 0.4,\"float roughness\", 0.1,\"string envname\", \"reflection\",\"float Kr\", 0.5,\"float samples\", 4)\n")
	elif shader == 'soapbubble' :
		self.entries.append("Attribute (\"string geometryset\", \"+reflection\")\n")
		self.entries.append("Shader (\"surface\", \"glass\", \"float Ka\", 0.2, \"float Kd\", 0.25,\"float Ks\", 0.4,\"float roughness\", 0.1,\"string envname\", \"reflection\",\"float Kr\", 0.5,\"float samples\", 4,\"float fakerefract\", 1)\n")
	elif shader == 'shownormals' :
		self.entries.append("Shader (\"surface\", \"shownormals\", \"float bias\", 1)\n")
	elif shader == 'showfacing' :
		self.entries.append("Shader (\"surface\", \"showfacing\", \"float Ka\", 0.25, \"float Kd\", 0.75)\n")

    def addColor(self, col):
        """
	Attribute ( "color C", (0.75, 0.75, 0.75) )
        """
	astr=("Attribute ( \"color C\", %s )\n"%self.vec3( col ))
	self.entries.append( astr )
#        str = '  texture {\n'
#        if not 'pigment' in texture.keys():
#            str = str + '      pigment { color rgb %s }\n' %self.vec3( col )
#        elif texture['pigment']=='':
#            str = str + '      pigment { color rgb %s }\n' %self.vec3( col )
#        else:
#            str = str + '      pigment { %s }\n' % texture['pigment']
#        for k,v in texture.items():
#            if k=='pigment': continue
#            if v is None:
#                str = str + '      %s\n' % (k,)
#            else:
#                str = str + '      %s { %s }\n' % (k,v)
#        str = str +  '  }\n'
       


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
                bondRad = 0.15,
		interpolation='linear'):
        """void <- addGeom(geometry, texture) adds a geometries to the povray
        scene
        Spheres, IndexedPolygons and IndexedPolylines are handled
        """
        from DejaVu2.Spheres import Spheres
        from DejaVu2.IndexedPolygons import IndexedPolygons
        from DejaVu2.IndexedPolylines import IndexedPolylines
        from DejaVu2.Cylinders import Cylinders
        print texture
        if isinstance(geom, Spheres):
            self.entries.append('#####Object %s#######\n' % geom.name)
            self.addSpheres(geom, texture)
        elif isinstance(geom, IndexedPolygons) and len(geom.faceSet):
            self.entries.append('#####Object %s#######\n\n' % geom.name)
            self.addIndexedPolgygons(geom, texture,interpolation)
#        elif isinstance(geom, IndexedPolylines) and len(geom.faceSet):
#            self.entries.append('// Object %s\n//\n' % geom.name)
#            self.addIndexedPolylines(geom, texture, bondRad)
        elif isinstance(geom, Cylinders) and len(geom.faceSet):
            self.entries.append('#####Object %s#######\n' % geom.name)
            self.addCylinders(geom, texture)
#        else:
#            print 'WARNING: %s the geometry is not yet supported'%geom.__class__

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

    def join_list(self, l):
                return ', '.join([str(i) for i in l])

    def mesh_geometry(self, name, single_sided, interpolation, nverts,\
                        verts, points, normals = None, vertexcolor = None, holes = None,transform = None):
		#"""
		#interp parameter may be "linear" or "catmull-clark" to indicate Catmull-Clark subdivision
		#nverts[0..n f aces - 1] contains the number of vertices in each face.
		#array verts, whose length is the sum of all the entries in nverts, contains the vertex indices of all faces.
		#"""
		points*=numpy.array((1.,1.,-1.),'f')
		
		astr=""
                if (transform != None ) : 
			astr=astr + self.append_transform(transform)
			#astr=astr + 'Scale (1.,1.,-1)\n'
                if (single_sided):
                        astr= astr +'Attribute ("int twosided", 0)\n'
		

                astr= astr + ('Mesh ("%s", (%s), (%s), "vertex point P", (%s)' %
                        (interpolation, self.join_list(nverts), self.join_list(verts), self.join_list(points.flatten())))

                if (normals != None):
			normals*=numpy.array((1.,1.,-1.),'f')
                        astr = astr +(', "vertex normal N", (%s)' % self.join_list(normals.flatten()))
		#else :
		#	astr = astr +(', "vertex normal N", (')
		#	for i in points :
		#		astr = astr +('surfacenormal %s,'%self.vec3(i))
		#	astr = astr + ')'
                if (vertexcolor != None):
                        astr= astr +(', "vertex color C", (%s)' % self.join_list(vertexcolor))

                if (holes != None):
                        astr= astr +(', "int[%d] holes", (%s)' %
                                (len(holes), self.join_list(holes)))

                astr= astr +')\n'
		return astr
    def sortPoly(self, geom,vt, order=-1):
        """None <- sortPoly(order=-1)
Sorts the geometry polygons according to z values of polygon's
geomtric centers. Order=-1 sorts by furthest z first, order=1 sorts
by closest z first"""
        # FIXME will not work with instance matrices
        #mat = geom.GetMatrix()
        #mat = numpy.reshape(mat, (4,4))
        #vt = geom.vertexSet.vertices*mat
        if vt is None:
            return
        triv = numpy.take(vt, geom.faceSet.faces.array, axis=0)
        trig = numpy.sum(triv,1)/3.
        trigz = trig[:,2]  #triangle's center of gravity z value
        
        ind = numpy.argsort(trigz) # sorted indices
        
        if len(geom.faceSet.faces.array):
            faces = numpy.take(geom.faceSet.faces.array, ind[::order],axis=0)
	    n = geom.getFNormals()
	    n = geom.faceSet.normals * geom.GetMatrix()
	    normals = numpy.take(n, ind[::order], axis=0)
            #
            #if geom.shading==GL.GL_FLAT: # we also need to re-arrange the
            #                          # face normals
            #	if geom.normals is None:
            #        normals = None
           # 	else:
           #         if len(geom.normals)>1:
            #normals = numpy.take(geom.normals, ind[::order], axis=0)
           #         else:
           #             normals = geom.normals
           # else:
           #     normals = None
	    #normals = None
            #geom.Set(faces=faces, fnormals=normals)
	    return faces.copy(),normals.copy()

    def addIndexedPolgygons(self, geom, texture, interpolation='linear'):
        """void <- addIndexedPolgygons(geom, texture)
	self.mesh_geometry(name, transform, single_sided, interpolation, nverts, verts, points, normals)
        """
	smooth=True
	#interpolation = 'catmull-clark'#'linear' - need correct normal to use linear interpolation
        single_sided = False
	name=geom.name
	#mat=matrix(geom.GetMatrix())
	#transform = numpy.array(mat.T).reshape(16,)
        transform = mat =  geom.GetMatrix()
	transform = numpy.array(transform.transpose()).reshape(16,)
	transform = None
        #vertices point
	vt = geom.vertexSet.vertices*mat
	#vt = geom.getVertices()#geom.vertexSet.vertices*mat #> use set transform of gelato : maybe compare performance
	points = vt#*numpy.array((1.,1.,-1.),'f')

	#vertex indice corrrespond to faces in Pmv
	#before getting the face aneed to sort them according the orientation in the viewer
	#f = geom.getFaces().copy()
	#normals = geom.getFNormals()
	faces,normals = self.sortPoly(geom,vt,order=1)
	faces,normals = self.sortPoly(geom,vt,order=1)
	#geom.sortPoly()
	#faces = geom.getFaces().copy()
	#print (f == geom.getFaces()).all()
	verts = faces
	
	#number of vertices for each faces
	nverts = []
	for i in faces :
		nverts.append(len(i))

        #lines = ["mesh2 {\n\tvertex_vectors {\n\t\t%d,\n"%len(vt)]
        #for v in vt:
        #    lines.append("\t\t%s,\n"%self.coord3(v))
        #lines.append("\t}\n")

        # add normals
        normals = geom.getVNormals()
	mat=geom.GetMatrix()[:3,:3]
	normals = numpy.dot(normals,mat.transpose())

        # handle colors
        colors = geom.materials[GL.GL_FRONT].prop[1]
        vtcolor = False    
	vertexcolor = []
        colBinding = 'Overall'
        if len(colors)==len(points):
            colBinding = 'Per Vertex'
	    vtcolor = True
	    for i in colors : vertexcolor.append(i[0:3])
	    vertexcolor = numpy.array(vertexcolor).flatten()
        elif len(colors)==len(faces):
            colBinding = 'Per Face'
	    vtcolor = True
	    for i in colors : vertexcolor.append(i[0:3])
	    vertexcolor = numpy.array(vertexcolor).flatten()
	else : vertexcolor = None
	# color
	self.entries.append('PushAttributes ()\n')
	if texture['pigment'] == "occlusion" or texture['pigment'] == "glass" : vertexcolor = None
	if texture['pigment'] != "" : self.addShader(texture['pigment'])
	if not vtcolor : self.addColor(colors[0][0:3])
	# alpha
        alpha = colors[0][3]
	if (alpha < 1.0):
           alpha = round(alpha, self.PRECISION)
           self.entries.append('Attribute (\"color opacity\", %s)\n' % self.vec3(alpha, alpha, alpha))
	mgeom=self.mesh_geometry(name, single_sided, interpolation, nverts, verts.flatten(), points,normals=normals, vertexcolor =vertexcolor,transform=transform)
	self.entries.append(mgeom)
	self.entries.append('PopAttributes ()\n')


    def addSpheres(self, geom, texture):
        """void <- addSpheres(geom, texture)
r.PushAttributes()
r.Scale(1.7,1.7,1.7)
r.Attribute ( "color C", (0., 0., 0.1) )
r.Translate ( 0, 0, 30.)
r.Sphere(1.0,-1.0,1.0,360.0)
r.PopAttributes()
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
        if texture['pigment'] != '' : self.addShader(texture['pigment'])
        for i in range(len(geom.vertexSet)):
	    self.entries.append('PushAttributes()\n')
            if fp:
                col = fp.prop[1][i][:3]
		alpha = fp.prop[1][i][3]
            else:
                col = (1.,1.,1.)
		alpha = 1.0
            self.addColor( col )
	    # alpha
	    if (alpha < 1.0):
              alpha = round(alpha, self.PRECISION)
              self.entries.append('Attribute (\"color opacity\", %s)\n' % self.vec3(alpha, alpha, alpha))
            self.entries.append('Translate %s\n'%self.coord3(v[i]) )
            if not geom.oneRadius:
                self.entries.append('Scale(%f,%f,%f)\n'%((radii[i][0]),(radii[i][0]),(radii[i][0])) )
            else:
                self.entries.append('Scale(%f,%f,%f)\n'%(geom.radius,geom.radius,geom.radius) )

	    self.entries.append("Sphere (1.0,-1.0,1.0,360.0)\n")
	    self.entries.append("PopAttributes()\n")

            #self.endShape()


    def addOneCylinder(self,v1,v2,color,alpha = 1.0):
	gv=[0.,0.,1]
	v=v2-v1
	mat=numpy.array(rotax.rotVectToVect(gv,v))
	# first translate on v1	
	#then rotate
	astr="PushTransform ()\n"
	astr = astr + "Attribute (\"color C\", %s)\n"% self.vec3(color)
	if (alpha < 1.0):
              alpha = round(alpha, self.PRECISION)	
              self.entries.append('Attribute (\"color opacity\", %s)\n' % self.vec3(alpha, alpha, alpha))
	astr = astr + "Translate %s\n" % self.coord3(v1)
	astr = astr + ("AppendTransform %s\n"%self.matrix16(mat.reshape(16,)))
	astr = astr + "nurbscyl ()\n"
	astr = astr + "PopTransform ()\n"
	return astr
	
	
    def addCylinders(self, geom, texture):
        """void <- addSpheres(geom, texture)
	PushTransform ()
	Attribute ("color C", %s)
	Translate (x, y, z)
	Rotate (a, x, y, z) OR AppendTransform (m16);
	nurbscyl ()
	PopTransform ()
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
	if texture['pigment'] != '' : self.addShader(texture['pigment'])
        f = geom.getFaces()
        for i in range(len(geom.faceSet)):
            for j in range(len(f[i])-1):
                #	self.entries.append('cone { \n')
                if not geom.oneRadius:
                    r1 = radii[f[i][j]]/2.
                    r2 = radii[f[i][j+1]]/2.
                else:
                    r1 = r2 = radii[0]/2.

                #self.entries.append('\t%s, %f, %s, %f\n'%(
                #     self.coord3(v[f[i][j]]), r1,
                #     self.coord3(v[f[i][j+1]]), r2))
                if fp:
                    col = c[f[i][j]][:3]
    		    alpha=c[f[i][j]][3]
                else:
                    col = (1.,1.,1.)
	            alpha = 1.0
		self.entries.append(self.addOneCylinder(v[f[i][j]],v[f[i][j+1]],col))
		self.entries.append(self.addOneCylinder(v[f[i][j+1]],v[f[i][j]],col))
                #self.entries.append("\tpigment { color rgb<%6.3f, %6.3f, %6.3f> }\n"% tuple(col))
                #self.entries.append("\tfinish { specular 1 roughness 0.001 ambient 0.3 }\n")
                #self.endShape()
#########From blendergelato.py#############################




if __name__ == '__main__':
    from DejaVu2 import Viewer
    vi = Viewer()
    from DejaVu2.Spheres import Spheres
    s = Spheres('test', centers = ( (0.0, 0., 0.), (3.,0.,0.) ) )
    s.Set(quality = 15)
    vi .AddObject(s)
    from DejaVu2.povray3 import PovRay
    p = PovRay()
    p.addCamera(vi.cameras[0])
    p.addLight(vi.lights[0])
    p.addGeom(s)
    p.write('test.pov')
    
    
