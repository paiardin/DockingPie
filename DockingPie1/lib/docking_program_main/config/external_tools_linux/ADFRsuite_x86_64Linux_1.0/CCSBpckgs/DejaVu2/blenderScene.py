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

from opengltk.OpenGL import GL
#from opengltk.extent.utillib import glTriangleNormals
from geomutils.geomalgorithms import  TriangleNormals

import numpy, math, string, warnings
from mglutil.math import rotax
from mglutil.util.misc import ensureFontCase
from DejaVu2.Geom import Geom
from DejaVu2.IndexedGeom import IndexedGeom
from DejaVu2.Transformable import Transformable
from DejaVu2.IndexedPolylines import IndexedPolylines
from DejaVu2.IndexedPolygons import IndexedPolygons
from DejaVu2.Spheres import Spheres
from DejaVu2.Cylinders import Cylinders
from DejaVu2 import viewerConst
from DejaVu2.colorTool import glMaterialWithCheck
from DejaVu2.GleObjects import GleExtrude
from DejaVu2.glfLabels import GlfLabels
from DejaVu2.DataOutput import *
from DejaVu2 import colorTool
import Tkinter, Pmw
from mglutil.gui.InputForm.Tk.gui import InputFormDescr, InputForm

class OutputBlender(OutputNode):
    def __init__(self):
        OutputNode.__init__(self)
        self.geomName = None
        self.Transformable = Transformable()
        self.completeFlag  = 1     # is set in getBlender()
        self.saveNormalsFlag   = 0 # is set in getBlender()
        self.colorPerVertexFlag = True # set in getBlender()
        self.usePROTOFlag    = 0   # if set, instanceMatrice geoms are not
                                   # saved as geoms, but re-used with PROTO 
        self.sphereQuality = 2     # default quality for sphere subsampling
        self.cylinderQuality = 10  # default quality for cyl. subsampling
        
        
    def getBlender(self, root, complete=1, normals=0, 
                 colorPerVertex=True,
                 usePROTO=0, sphereQuality=2, cylinderQuality=10):
        """ this method returns a list of strings describing DejaVu2 Geoms in
        blender python API format """

        # this is the method the user should call
        
        self.completeFlag = complete # if 1, header and footer will be added
        self.saveNormalsFlag = normals # if 1, normals are added
        self.colorPerVertexFlag = colorPerVertex # if False, color per face
        self.usePROTOFlag = usePROTO # if set to 1, instance geoms are saved
                                     # with one PROTO. Else: data for these
                                     # geoms will be prepared which blows up
                                     # the file size

        self.sphereQuality = sphereQuality # default is 2
        self.cylinderQuality = cylinderQuality # default is 10

        self.output = []
        
        if self.completeFlag:
            self.output.extend(self.getCopyright())
            self.output.extend(self.getFileHeader1())

            
        self.loopGeoms(root)

        return self.output


    def NodeRepr(self, geom, matrix):
        # called recursively from loopGeom()
        if not self.checkGeom(geom): return

        # call the progress bar update
        self.configureProgressBar(labeltext='parse '+geom.name)
        self.updateProgressBar()
        
        # IndexedGeoms, Spheres and Cylinders have to be treated differently
        if isinstance(geom, IndexedPolygons) or \
           isinstance(geom, IndexedPolylines):
            self.output.extend( self.doit(geom, matrix) )

        elif isinstance(geom, Spheres):
            # convert Spheres geom into IndexedPolygons
            sphGeom = geom.asIndexedPolygons(quality=self.sphereQuality) 
            sphGeom.viewer = geom.viewer
            sphGeom.parent = geom.parent
            sphGeom.name = geom.name
            sphGeom.fullName = geom.fullName
            sphGeom.instanceMatrices = geom.instanceMatrices
            sphGeom.instanceMatricesIndex = geom.instanceMatricesIndex
            sphGeom.VRML2CreatedPROTOForThisGeom = 0
            self.output.extend( self.doit(sphGeom, matrix) )

        elif isinstance(geom, Cylinders):
            # convert Cylinders geom into IndexedPolygons
            cylGeom = geom.asIndexedPolygons(quality=self.cylinderQuality) 
            cylGeom.viewer = geom.viewer
            cylGeom.parent = geom.parent
            cylGeom.name = geom.name
            cylGeom.fullName = geom.fullName
            cylGeom.instanceMatrices = geom.instanceMatrices
            cylGeom.instanceMatricesIndex = geom.instanceMatricesIndex
            cylGeom.VRML2CreatedPROTOForThisGeom = 0
            self.output.extend( self.doit(cylGeom, matrix) )
            
        elif isinstance(geom, GlfLabels):
            # convert Cylinders geom into IndexedPolygons
            glfGeom = geom.asIndexedPolygons() 
            glfGeom.viewer = geom.viewer
            glfGeom.parent = geom.parent
            glfGeom.name = geom.name
            glfGeom.fullName = geom.fullName
            glfGeom.instanceMatrices = geom.instanceMatrices
            glfGeom.instanceMatricesIndex = geom.instanceMatricesIndex
            glfGeom.VRML2CreatedPROTOForThisGeom = 0
            self.output.extend( self.doit(glfGeom, matrix) )


    def doit(self, geom, matrix):
        # gets called from NodeRepr()
        blender = []
        
        # if self.usePROTO is set to 1: don't convert instance matrices
        # geoms into geoms, but use the vrml2 USE 
        if self.usePROTOFlag:
            # create all the necessary data to be put in the PROTO which goes
            # in the header of the vrml2 file
            if geom.VRML2CreatedPROTOForThisGeom == 0:
                name = self.getGeomName(geom)
                blender.append("PROTO "+name+" [ ] {\n")
                identityMatrix = numpy.identity(4).astype('f')
                blender.extend( self.doitReally(geom, identityMatrix) )
                blender.append("}\n")

                # now insert this into the header of self.output
                for i in range(len(vrml2)):
                    self.output.insert(i+1,blender[i])
                geom.VRML2CreatedPROTOForThisGeom = 1 # don't add this geom
                                                      # to the header next time
                # this PROTO flag will be deleted when we leave the recursive
                
            # and add it as USE to the body
            blender = [] # clear the list because this data has been added
            blender.extend( self.doitUsePROTO(geom, matrix) )
            return blender

        else:
            # here we save all the data for all geoms
            return self.doitReally(geom, matrix)


    def doitReally(self, geom, matrix):
        # add header for geom
        blender = []

        blender.extend(self.getGeomHeader(geom))
        mat, colors = self.getMaterial(geom)

        # add texture if applicable:
        #if geom.texture:
        #    blender.extend( self.getTexture(geom) )

        # add coords, faces, etc
        blender.extend( self.getGeomDescr(geom) )

        # add texCoord Coordinates is applicable
        #if geom.texture:
        #    blender.extend( self.getTexCoords(geom) )

        # add texCoordsIndex if applicable
        #if hasattr(geom.vertexSet,'texCoordsIndex'):
        #     blender.extend( self.getTexCoordsIndex(geom))
        
        # add normals if applicable
        #if self.saveNormalsFlag and isinstance(geom, IndexedPolygons):
        #    blender.extend( self.getNormals(geom) )

        # add colors per vertex if applicable
        if colors is not None and len(colors): 
            blender.extend( self.getColors(geom, colors) )

        # add transformations for geom
        #blender.extend( self.getTransforms(matrix) )

        # add closing bracket for Transform{}
        #blender.append("      }\n")
        return blender
    

    def doitUsePROTO(self, geom, matrix):
        # FIXME: this works currently only with geoms that are not grouped
        # i.e. it doesnt work with secondary structures, they will be saved
        # as PROTO too, but also for each instanceMatrix (->redundant)

        vrml2 = []
        
        geom.instanceMatricesIndex = 0
        name = string.split(self.getGeomHeader(geom)[0])[1]
        vrml2.append("    Transform {\n")
        vrml2.append("      children  "+name+" { }\n")

        # add transformations for geom
        vrml2.extend( self.getTransforms(matrix) )

        # add closing bracket for Transform{}
        vrml2.append("      }\n")
        return vrml2
    


    def getFileHeader1(self):
        blender=[]
        blender.append("#!BPY\n")
	blender.append("""import sys, os, os.path, struct, math, string
import Blender
import bpy
from Blender import *
from Blender.Mathutils import *
from Blender import Object
from Blender import Material
from Blender import Window, Scene, Draw
import BPyMesh

def createsNmesh(name,vertices,vnormals,faces,mat=None):
	vlist = []
	me=bpy.data.meshes.new(name)
	me.verts.extend(vertices)	# add vertices to mesh
	me.faces.extend(faces)          # add faces to the mesh (also adds edges)
	#smooth face : the vertex normals are averaged to make this face look smooth
	for face in me.faces:
		face.smooth=1
	me.calcNormals()
	if mat == None :
		mat = Material.New('test')
		mat.R = 0.8
		mat.G = 1.0
		mat.B = 1.0
	me.materials=[mat]
	ob = Blender.Object.New("Mesh")
	ob.link(me)
	return ob,me

sc=Blender.Scene.GetCurrent()

""")
	blender.append("\n")
        return blender


    def getCopyright(self):
        blender=[]
	blender.append("""
#############################################################################
#Blender Scene Export v1.1 by Ludovic Autin
# Date: August 2009   Author: Ludovic Autin
#
# autin@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Ludovic Autin and TSRI
#
#
#############################################################################
""")   
	blender.append("\n")
        return blender


    def getGeomName(self, geom):
        g = geom
        name = "Pmv_"
        while g != geom.viewer.rootObject:
            # g.name can contain whitespaces which we have to get rid of
            gname = string.split(g.name)
            ggname = "" 
            for i in gname:
                ggname = ggname + i
            name = name + string.strip(ggname)+"AT"+\
                   string.strip(str(g.instanceMatricesIndex))+ '_'
            g = g.parent
        name=string.replace(name,"-","_")
	return name


    def getGeomHeader(self, geom):
        # generates geom name
        blender=[]
        g = geom
        name = self.getGeomName(geom)
        blender.append("#Geom "+name+"\nname=\""+name+"\"\n")
	blender.append("mat=None\n")
        return blender


    def getMaterial(self, geom):
        blender=[]
	name = self.getGeomName(geom)
        mat = geom.materials[GL.GL_FRONT].prop[:]
        geom.materials[GL.GL_FRONT].colorIndex = None # will be used later on
        colors = None

        # if only 1 color present, skip this all and use the ambient definition
        # below
        if len(mat[1])> 1:
            colors = mat[1]

            # The ZCorp printer software doesn't support color_per_face,
            # but Pmv does. So, we create a colorIndex list for special cases
            # However, the user can still choose between coloring
            # per face and per vertex

            # FIXME: test for primitive type, i.e. tri_strip or quad_strip
            # currently this works only for tri_strips
            if isinstance(geom, GleExtrude): # special case!
                faces = geom.faceSet.faces.array # triangle_strips
                ifaces = geom.getFaces() # indexed geom

                # if the user forces to save color per vertex:
                if self.colorPerVertexFlag is True:
                    colorIndex = numpy.zeros( (ifaces.shape[0], \
                                                 ifaces.shape[1]) )
                    c = 0
                    cc = 0
                    for face in faces:
                        for j in range(len(face)-2): # -2 because of tri_strip
                            colorIndex[cc] = c
                            cc = cc + 1
                        c = c + 1
                    geom.materials[GL.GL_FRONT].colorIndex = colorIndex
                    
            elif isinstance(geom, IndexedPolygons):
                mat[1]=[mat[1][0]]
                vertices = geom.getVertices()
                faces = geom.getFaces()

                # if current colors are per face:
                if len(colors) != len(vertices) and len(colors) == len(faces):
                    # if the user forces colors per vertices:
                    if self.colorPerVertexFlag is True:
                        colorIndex = numpy.zeros( (faces.shape[0], \
                                                     faces.shape[1]) )
                        c = 0
                        for face in faces:
                            for f in face:
                                colorIndex[c] = c 
                            c = c + 1
                        geom.materials[GL.GL_FRONT].colorIndex = colorIndex

                # if current colors are per vertex
                else:
                    # if the user forces colors per face:
                    if self.colorPerVertexFlag is False:
                        # code from Michel Sanner follows (thanks Michel!):
                        vcol = geom.materials[1028].prop[1]
                        tri = geom.faceSet.faces.array
                        verts= geom.vertexSet.vertices.array
                        colors = []
                        for t in tri:
                            s1,s2,s3 = t
                            col = ( (vcol[s1][0]+vcol[s2][0]+vcol[s3][0])/3.,
                                    (vcol[s1][1]+vcol[s2][1]+vcol[s3][1])/3.,
                                    (vcol[s1][2]+vcol[s2][2]+vcol[s3][2])/3. )
                            colors.append( col)



        ambInt =  '%.5f'%mat[0][0][0]
        difCol =  '%.5f'%mat[1][0][0]+" "+'%.5f'%mat[1][0][1]+" "+\
                  '%.5f'%mat[1][0][2]
        emCol =   '%.5f'%mat[2][0][0]+" "+'%.5f'%mat[2][0][1]+" "+\
                  '%.5f'%mat[2][0][2]
        specCol = '%.5f'%mat[3][0][0]+" "+'%.5f'%mat[3][0][1]+" "+\
                  '%.5f'%mat[3][0][2]
        shin =    '%.5f'%mat[4][0]
        trans =    `1-mat[5][0]`

	blender.append("mat = Material.New('mat"+name+"')")
	blender.append("mat.R = '%.5f'n" %mat[1][0][0])
	blender.append("mat.G = '%.5f'n" %mat[1][0][1])
	blender.append("mat.B = '%.5f'n" %mat[1][0][2])

	#vrml2.append("          material        Material {\n")
        #vrml2.append("            ambientIntensity      "+ambInt+"\n")
        #vrml2.append("            diffuseColor          "+difCol+"\n")
        #vrml2.append("            emissiveColor         "+emCol+"\n")
        #vrml2.append("            specularColor         "+specCol+"\n")
        #vrml2.append("            shininess             "+shin+"\n")
        #vrml2.append("            transparency          "+trans+"\n")
        #vrml2.append("          }\n")
        return blender, colors


    def getGeomDescr(self, geom):
        blender = []
	name = self.getGeomName(geom)
        if isinstance(geom, IndexedPolygons):
	    #blender.append("vertices=\n")
            # add vertices
            blender.extend( self.getCoords(geom) )
            # add face indices
            blender.extend( self.getFaces(geom) )
	    blender.extend( self.getNormals(geom) )
		
            # add color indices if applicable
            #if geom.materials[GL.GL_FRONT].colorIndex is not None and len(geom.materials[GL.GL_FRONT].colorIndex):
            #    vrml2.extend( self.getColorIndex(geom) )
            blender.append(""+name+",mesh"+name+"=createsNmesh(name,vertices,vnormals,faces,mat=mat)\n")
	    blender.append("sc.link("+name+")\n")
        #elif isinstance(geom, IndexedPolylines):
        #    vrml2.append("        geometry        IndexedLineSet {\n")
        #    # add vertices
        #    vrml2.extend( self.getCoords(geom) )
        #    # add face indices
        #    vrml2.extend( self.getFaces(geom) )

        return blender


    ## added by A Gillet 04/13/2006
    def getTexture(self,geom):
        """ return PixelTexture Node
        PixelTexture {
             image 0  0 0  # exposedField SFImage
             repeatS True  # field SFBool
             repeatT True  # field SFBool
             }

        the value of the image field specifies image size and pixel values
        for a texture image
            width (in pixel)
            height (in pixel)
            number of 8-bit bytes for each pixel
               recognize values are:
                   0 disable texturing for shape
                   1 Grayscale
                   2 Grayscale with alpha
                   3 RGB
                   4 RGB with alpha
        (Info taken from Book " VRML 2.0 source book by Andrea L. Ames,
        David R. Nadeau and John L. Moreland ")

        """


        vrml2=[]
        tex = geom.texture
        dims = tex.image.shape
        
        vrml2.append("\n")
        vrml2.append("      texture PixelTexture {\n")
        width  = dims[0]
        if len(dims) == 3:
            height = dims[1]
            num_byte = dims[2]
        elif len(dims)==1:
            height = 1
            num_byte = len(tex.image[0])
        elif len(dims)==2:
            height = 1
            num_byte = dims[1]
            
        vrml2.append("          image "+`width`+" "+`height`+" "+`num_byte`+"\n")

        if len(dims) == 3:
            # we have a 2D texture (image)
            countW =0
            for r in tex.image: # row
                for c in r:     # column
                    istring = "0x"
                    for i in range(3):
                        hexa = "%X"%c[i]
                        if len(hexa)==1: hexa = "0"+hexa
                        istring = istring+ hexa
                    istring = istring + " "
                    vrml2.append(istring)
            vrml2.append("\n")
            vrml2.append("      }\n")

        else:
            # we have a 1-dimensional array
            for line in tex.image:
                istring = "0x"
                for i in range(len(line)):
                    hexa = "%X"%line[i]
                    if len(hexa)==1: hexa = "0"+hexa
                    istring = istring+ hexa
                istring = istring + "\n"
                vrml2.append(istring)        
            vrml2.append("      }\n")
        return vrml2




##  Daniel Stoffler code
##     def getTexture(self, geom):
##         vrml2=[]
##         dims = geom.texture.image.shape
##         vrml2.append("\n")
##         vrml2.append("      texture PixelTexture {\n")
##         # FIXME : what are real dimensions of image?
##         # I never tested this for images larger than one-dimensional array
##         vrml2.append("          image "+`dims[0]`+" "+`1`+" "+\
##                      `len(geom.texture.image[0])`+"\n")
##         for line in geom.texture.image:
##             istring = "0x"
##             for i in range(len(line)):
##                 hexa = "%X"%line[i]
##                 if len(hexa)==1: hexa = "0"+hexa
##                 istring = istring+ hexa
##             istring = istring + "\n"
##             vrml2.append(istring)
##         vrml2.append("      }\n")
##         return vrml2


    def getCoords(self, geom):
        blender=[]
        vertices = geom.getVertices()
        blender.append("vertices=[")
        for vert in vertices:
            vstring = "[%.5f,%.5f,%.5f]," % (vert[0],vert[1],vert[2])
            blender.append(vstring)
	blender[-1][-1].replace(",","")
	blender.append("]\n")
        return blender

    def getFaces(self, geom):
        #print "getFaces"
        blender=[]
        faces = geom.getFaces()
        blender.append("faces=[")
        for face in faces:
            facestring = "[%d,%d,%d]," % (face[0],face[1],face[2])
            blender.append(facestring)

        #for face in faces:
        #    facestring = "                          "
        #    if geom.invertNormals: # reverse faces
        #        facestring = facestring + `face[0]` + ", "
        #        for i in range(len(face)-1,0,-1): 
        #            facestring = facestring + `face[i]` + ", "
        #    else:
        #        for f in face:
        #            facestring = facestring + `f` + ", "#
	#
        #    facestring = facestring + "-1,\n"
        #    blender.append(facestring)
	blender[-1][-1].replace(",","")
	blender.append("]\n")
        return blender

    def getColorIndex(self, geom):
        # only called if len(colors) != len(faces)
        vrml2 = []
        colorIndex = geom.materials[GL.GL_FRONT].colorIndex
        vrml2.append("	  colorIndex	[\n")
        for cI in colorIndex:
            cIstring = "         "
            for c in cI:
                cIstring = cIstring + `c` +", "
            cIstring = cIstring + "-1,\n"
            vrml2.append(cIstring)
        vrml2.append("	                ]\n")
        # clean up the geom object
        del geom.materials[GL.GL_FRONT].colorIndex
        return vrml2


    def getNormals(self, geom):
        if geom.invertNormals:
            fn = -1.0 * geom.normals
	    vn = -1.0 * geom.getVNormals()
        else:
            fn = geom.normals
	    vn = geom.getVNormals()
        blender=[]
        blender.append("vnormals=[")
        for n in vn:
            blender.append("[%d,%d,%d]," % (n[0],n[1],n[2]))
	blender[-1][-1].replace(",","")
	blender.append("]\n")
        return blender

    def getTexCoords(self, geom):
        vrml2=[]
        vrml2.append("\n")
        vrml2.append("      texCoord TextureCoordinate {\n")
        vrml2.append("        point [\n")
        for p in geom.vertexSet.texCoords.array:
            if len(p) == 1: # one dimension array
                vrml2.append("          "+`p[0]`+" 0,\n")
            else:
                vrml2.append("          "+`p[0]`+" "+`p[1]`+",\n")
        vrml2.append("          ]\n")
        vrml2.append("        }\n")
        return vrml2

    def getTexCoordsIndex(self, geom):
        vrml2=[]
        texCoordsIndex = geom.vertexSet.texCoordsIndex.array
        vrml2.append("	  texCoordIndex	[\n")

        for face in texCoordsIndex:
            indexstring = "                          "
            for i in face:
                indexstring = indexstring + `i` + ", "

            indexstring = indexstring + "-1,\n"
            vrml2.append(indexstring)
        vrml2.append("	                ]\n")
        return vrml2

    def blenderColor(self,col):
	    if max(col)<=1.0: col = map( lambda x: x*255, col)
	    return col

    def getColors(self, geom, colors):
	name = self.getGeomName(geom)
        blender=[]
        blender.append("\n")
	

        blender.append("colors=[") #color foreach vertices
        for c in colors:
	    c=self.blenderColor(c)
            cstring = "[%d,%d,%d]," % (c[0],c[1],c[2])
            blender.append(cstring)
	blender[-1][-1].replace(",","")
	blender.append("]\n\n")
	
	blender.append("mesh"+name+".vertexColors = 1  # enable vertex colors \n")
	blender.append("for f in mesh"+name+".faces:\n")
        blender.append("       for i, v in enumerate(f):\n")
        blender.append("               col= f.col[i]\n")
        blender.append("               col.r= colors[v.index][0]\n")
        blender.append("               col.g= colors[v.index][1]\n")
        blender.append("               col.b= colors[v.index][2]\n")
	blender.append("mesh"+name+".materials[0].setMode(\"VColPaint\")\n")
        return blender	


    def getTransforms(self, matrix):
        vrml2=[]
        #mymatrix = matrix.__copy__()
        #mymatrix = numpy.reshape(mymatrix, (16,))
        mymatrix = numpy.reshape(matrix, (16,))
        rot,trans,scale=self.Transformable.Decompose4x4(mymatrix)
        r = rotax.mat_to_quat(rot)
        r[3]=r[3]*math.pi/180.0 #convert to rad
        r[0] = round(r[0],6)
        r[1] = round(r[1],6)
        r[2] = round(r[2],6)
        r[3] = round(r[3],6)
        vrml2.append("      translation "+`trans[0]`+" "+`trans[1]`+" "+\
                           `trans[2]`+"\n")
        vrml2.append("      rotation "+`r[0]`+" "+`r[1]`+" "+\
                           `r[2]`+" "+`r[3]`+"\n")
        vrml2.append("      scale "+`scale[0]`+" "+`scale[1]`+" "+\
                               `scale[2]`+"\n")
        return vrml2

    def write(self,entries, filename):
        """void <- write(filename)  dumps gelato scene description"""
        self.filename = filename
        f = open(filename, 'w')
        for entry in entries:
            f.write(entry)
        f.close()

#from DejaVu2 import blenderScene
#outB=blenderScene.OutputBlender()
#m=self.getMolFromName('1crn')
#g=m.geomContainer.geoms['MSMS-MOL']
#g=m.geomContainer.masterGeom
#test=outB.getBlender(g,complete=1)
#outB.write(test,"test_blender.py")




class BlenderGUI(DatGUI):
    """This is the gui for OutputBlender:
    - Save normals adds normals to all geoms to be saved as vrml2
    - Invert normals inverts the normals for all geoms to be saved as vrml2
    - colorPerVertex: True by default. If set to False, color per face is used
    - Using PROTO will define a prototype geom which can be reused. This
      is usefull to lower the file size when instanceMatrices are applied.
    - Sphere quality is the subsampling of the spheres. Default value is 2
      lowest allowed value is 0. Higher values increase the file size
      significantly.
    - Cylinder quality is the subsampling of cylinders. Default value is 10
      lowest allowed value is 3.
    """

    def __init__(self, master=None, title=None):
        DatGUI.__init__(self, master, title)

        self.saveNormals    = Tkinter.IntVar()
        self.colorPerVertex = Tkinter.IntVar()
        self.colorPerVertex.set(1) # on by default
        self.usePROTO         = Tkinter.IntVar()

        self.createForm()

    def createForm(self):
        row = 0
        self.idf.append({'name':'savNormals',
                         'widgetType':Tkinter.Checkbutton,
                         'wcfg':{'text':'Save Normals            ',
                                 'variable':self.saveNormals,
                                 'command':self.saveNormals_cb},
                         'gridcfg':{'sticky':'w',
                                    'columnspan':2, 'row':row, 'column':0}
                         })
        row+= 1
        self.idf.append({'widgetType':Tkinter.Frame,
                 'wcfg':{'relief':'sunken','borderwidth':2,'height':2},
                 'gridcfg':{'columnspan':2, 'row':row, 'column':0},
                 })

        row += 1
        self.idf.append({'name':'colorPerVertex',
                         'widgetType':Tkinter.Checkbutton,
                         'wcfg':{'text':'color per vertex',
                                 'variable':self.colorPerVertex,
                                 'command':self.colorPerVertex_cb},
                         'gridcfg':{'sticky':'w',
                                    'columnspan':2, 'row':row, 'column':0}
                         })
        
        row+= 1
        self.idf.append({'widgetType':Tkinter.Frame,
                 'wcfg':{'relief':'sunken','borderwidth':2,'height':2},
                 'gridcfg':{'columnspan':2, 'row':row, 'column':0},
                 })

        row+= 1
        self.idf.append({'name':'usePROTO',
                         'widgetType':Tkinter.Checkbutton,
                         'wcfg':{'text':'Use PROTO for instance\n'+\
                                 'matrices to lower file size',
                                 'variable':self.usePROTO,
                                 'command':self.usePROTO_cb},
                         'gridcfg':{'sticky':'w',
                                    'columnspan':2, 'row':row, 'column':0}
                         })
        row+= 1
        self.idf.append({'widgetType':Tkinter.Frame,
                 'wcfg':{'relief':'sunken','borderwidth':2,'height':2},
                 'gridcfg':{'columnspan':2, 'row':row, 'column':0},
                 })

        row+= 1
        self.idf.append({'widgetType':Tkinter.Label,
                 'wcfg':{'text':'Sphere quality'},
                 'gridcfg':{'sticky':'w','columnspan':2, 'row':row,
                            'column':0},
                 })
        
        self.idf.append({'name':'inpSphQual',
                         'widgetType':Tkinter.Entry,
                         'defaultValue':'2',
                         'wcfg':{'font':(ensureFontCase('Courier'),10),
                                 'width':5,'textvariable':self.sphInput,
                                 'command':self.sphereQuality_cb,
                                 },
                         'gridcfg':{'sticky':'e',
                                    'columnspan':2, 'row':row, 'column':1 }
                         })
        row+= 1
        self.idf.append({'widgetType':Tkinter.Label,
                 'wcfg':{'text':'Cylinder quality'},
                 'gridcfg':{'sticky':'w','columnspan':2, 'row':row,
                            'column':0},
                 })

        self.idf.append({'name':'inpCylQual',
                         'widgetType':Tkinter.Entry,
                         'defaultValue':'10',
                         'wcfg':{'font':(ensureFontCase('Courier'),10),
                                 'width':5,'textvariable':self.cylInput,
                                 'command':self.cylinderQuality_cb,
                                 },
                         'gridcfg':{'sticky':'e',
                                    'columnspan':2, 'row':row, 'column':1 }
                         })

        row+= 1
	self.idf.append({'widgetType':Tkinter.Button,
                         'text':'OK',
                         'wcfg':{},
                         'gridcfg':{'sticky':'wens',
                                    'columnspan':1, 'row':row, 'column':0},
                         'command': self.OK_cb})


	self.idf.append({'widgetType':Tkinter.Button,
                         'text':'Cancel',
                         'wcfg':{},
                         'gridcfg':{'sticky':'wens',
                                    'columnspan':1, 'row':row, 'column':1},
                         'command': self.Cancel_cb})
        

    
    def saveNormals_cb(self):
        pass


    def colorPerVertex_cb(self):
        pass


    def usePROTO_cb(self):
        pass


    def getValues(self):
        vals = {}
        vals['saveNormals'] = self.saveNormals.get()
        co = self.colorPerVertex.get()
        if co == 1:
            co = True
        else:
            co = False
        vals['colorPerVertex'] = co
        vals['usePROTO'] = self.usePROTO.get()
        vals['sphereQuality'] = self.sphereQuality
        vals['cylinderQuality'] = self.cylinderQuality
        return vals

