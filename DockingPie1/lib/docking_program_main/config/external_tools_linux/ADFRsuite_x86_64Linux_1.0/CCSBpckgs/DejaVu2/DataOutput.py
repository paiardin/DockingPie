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
# Date: Jun 2002   Author: Daniel Stoffler
#
# stoffler@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Daniel Stoffler and TSRI
#
# Revision: Guillaume Vareille, Alex Gillet
#
#############################################################################


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

import Tkinter, Pmw
from mglutil.gui.InputForm.Tk.gui import InputFormDescr, InputForm


##########################################################################
#
#  functions for generating and parsing SMF format:
#
#  http://www.csit.fsu.edu/~burkardt/data/smf/smf.txt
#  (also DejaVu2/FileFormats/smf.txt)
#
#  use my Michael Garland's Qslim program 
#
##########################################################################


def IndexedPolgonsAsSMFString(geometry):
    """listOfStrings <-- IndexedPolgonsAsSMFString(geometry)
For a given IndexedPolygons geoemtry this function generates the textual SMF
description.  Currently supports vertices (v), faces (f), normals(n) and
colors (c).  Normals and color binding mstring are generated automatically.
Texture indices (r) are not used yet but could be used to store any integer
property.
"""
    assert isinstance(geometry, IndexedPolygons)
    # we can only handle triangles or quads
    assert geometry.faceSet.faces.array.shape[1] in [3,4]

    lines = []
    lines.append('begin\n')
    # vertices
    verts = geometry.vertexSet.vertices.array
    lines.extend( map(lambda x: ('v %f %f %f\n'%tuple(x)), verts) )

    # faces
    tri = geometry.faceSet.faces.array
    lines.extend( map(lambda x:
                      ('f %d %d %d\n'%(x[0]+1,x[1]+1,x[2]+1)), tri) )

    # normals
    normals = geometry.normals
    if len(normals)==len(verts):
        lines.append('bind n vertex\n')
    elif len(normals)==len(tri):
        lines.append('bind n face\n')
    else:
        normals = None
    if normals:
        lines.extend( map(lambda x:
                          ('n %f %f %f\n'%tuple(x)), normals) )

    # colors
    if geometry.materials[1028].binding[1]==11: #per vertex color
        lines.append('bind c vertex\n')
        cols = geometry.materials[1028].prop[1]
        lines.extend( map(lambda x:
                          ('c %f %f %f\n'%tuple(x[:3])), cols) )
    lines.append('end\n')
    return lines

    
def writePolygonsAsSMF(geometry, filename):
    """Write a IndexedPolygons geometry to a file"""

    f = open(filename, 'w')
    map( lambda x, f=f: f.write(x), IndexedPolgonsAsSMFString(geometry) )
    f.close()


def ParseSMFString(stringList):
    """v,f,n,c,r <-- ParseSMFString(stringList)
Parse an ascii SMF file and returnes a list of 3D vertices, triangular faces
(0-based indices, faces with more edges generate warnings), normals, colors
and 2D texture coordinates.
"""
    
    vertices = []
    faces = []
    normals = []
    colors = []
    textures = []
    fi = 0
    for l in stringList:
        w = l.split()
        if w[0]=='v': # vertices
            vertices.append( [float(w[1]),float(w[2]),float(w[3])] )
        elif w[0]=='f': # faces
            if len(w) > 4:
                warnings.warn("face %d has more than 3 edges"%fi);
            faces.append( [int(w[1])-1,int(w[2])-1,int(w[3])-1] )
            fi += 1
        if w[0]=='n': # normal vectors
            normals.append( [float(w[1]),float(w[2]),float(w[3])] )
        if w[0]=='c': # colors
            colors.append( [float(w[1]),float(w[2]),float(w[3])] )
        if w[0]=='r': # 2D texture indices
            textures.append( [float(w[1]),float(w[2])] )
    return vertices, faces, normals, colors, textures


def readSMF(filename):
    """Read an SMF ascii file and return an IndexedPolygons geometry"""

    f = open(filename)
    data = f.readlines()
    f.close()
    v, f, n, c, t = ParseSMFString(data)
    return v, f, n, c, t


class OutputNode:

    """This base class recursively loops over a given DejaVu2 geom tree and
    computes the transformation matrices at each level."""

    def __init__(self):
        self.output = [] # list describing the formated data
        self.lenGeoms = 0 # used for the progress bar to determine how many
                          # visible geoms are in the viewer

    def configureProgressBar(self, **kw):
        # this method is to be implemented by the user from outside
        pass


    def updateProgressBar(self, progress=None):
        # this method is to be implemented by the user from outside
        pass


    def countGeoms(self, obj):
        for child in obj.children:
            if child.visible:
                check = self.checkGeom(child)
                if check:
                    self.lenGeoms = self.lenGeoms + 1
                self.countGeoms(child)
            

    def loopGeoms(self, obj):
        """ this calls the recursive method """

## #commented out 2.7.02 DS, we discard now root transformation
##          # we reset the scale of rootObject to 1,1,1
##          if obj == obj.viewer.rootObject:
##              scaleFactor = obj.scale
##              obj.SetScale([1.,1.,1.])

        # initialize the progress bar to be ratio
        self.lenGeoms = 0 # initialize this to 0
        self.countGeoms(obj) # determine lenght of visible geoms
        self.lenGeoms = self.lenGeoms + 1 # FIXME: dont know why whe
                                              # always have 1 too few
        self.configureProgressBar(init=1, mode='increment',
                                  granularity=1,
                                  progressformat='ratio',
                                  labeltext='parse geoms to vrml2',
                                  max=self.lenGeoms)
        GL.glPushMatrix()
        GL.glLoadIdentity()
        self._loopGeomsRec(obj) # call recursive method
        GL.glPopMatrix()

## #commented out 2.7.02 DS, we discard now root transformation
##          # now we set the scale factor back to what it was
##          if obj == obj.viewer.rootObject:
##              obj.SetScale(scaleFactor)


    def _loopGeomsRec(self, obj):
        """ recursive method """

	GL.glPushMatrix()

        # we discard root object transformation:
        if obj is not obj.viewer.rootObject:
            if hasattr(obj, 'MakeMat'):
                obj.MakeMat()

        obj.VRML2CreatedPROTOForThisGeom = 0 # flag used in vrml2 doit()
            
        for i in range(len(obj.instanceMatricesFortran)):
            GL.glPushMatrix()
            GL.glMultMatrixf(obj.instanceMatricesFortran[i])
            obj.instanceMatricesFortranIndex = i # flag used in stl and vrml2 doit()
            
            matrix = GL.glGetDoublev(GL.GL_MODELVIEW_MATRIX)
            self.NodeRepr(obj, matrix)

            for child in obj.children:
                
                if child.visible:
                    self._loopGeomsRec(child)

            GL.glPopMatrix()

	GL.glPopMatrix()     # Restore the matrix
        del obj.instanceMatricesFortranIndex # clean up object
        del obj.VRML2CreatedPROTOForThisGeom 


    def NodeRepr(self, obj, matrix):
        """ This method, to be implemented by sublcass, should generate
        stl, vrml2 etc descriptions of DejaVu2 geoms"""
        pass


    def checkGeom(self, geom):
        # decides which geoms will be output. return 0 means don't save this
        if geom is None:
            return 0
        elif not isinstance(geom, Geom):
            return 0
        elif not geom.visible:
            return 0
        elif len(geom.vertexSet)==0:
            return 0
        elif isinstance(geom, IndexedGeom) and len(geom.faceSet) == 0:
            return 0
        else: return 1


class OutputSTL(OutputNode):
    """ generates a list of strings describing DejaVu2 geoms in STL
    (stereolithography) format (don't mix this with standard template
    library).
        """
# THIS STUFF HERE WAS USED IN GLEObject
##      def getSTL(self, reverse=0, **kw):
##          """ returns a string describing vertices and normals in the STL format.
##          """
##          from OpenGL import GL

##          # force per face normal calculation
##          if self.shading != GL.GL_FLAT:
##              oldmode = self.shading
##              oldnorm = self.normals
##              self.shading = GL.GL_FLAT
##              self.GetNormals()
##              norms = self.normals
##              self.shading = oldmode
##              self.normals = oldnorm

##          faces = self.getFaces()
##          if len(faces) < 1:
##              raise RuntimeError("No faces found in geometry ",self.name)

##          if faces.shape[1] > 3:
##              raise RuntimeError("More than 3 indices per face in ",self.name)

##          vert = self.getVertices()
##          if len(vert)<3:
##              raise RuntimeError("Less that 3 vertices found in geometry ",
##                                 self.name)

##          fn = norms
##          if len(fn)!=len(faces):
##              raise RuntimeError("Number of face normals does not match \
##              number of faces in geoemtry ",self.name)

##          if reverse:
##              fn = -1.0*fn
##              l=[]
##              for face in faces:
##                  l.append([face[0],face[2],face[1]])
##              faces = numpy.array(l)
            
##          stl = []
##          for i in xrange(len(faces)):
##              stl.append("    facet normal %f %f %f\n"%tuple(fn[i]))
##              stl.append("        outer loop\n")
##              fa = faces[i]
##              stl.append("            vertex %f %f %f\n"%tuple(vert[fa[0]]))
##              stl.append("            vertex %f %f %f\n"%tuple(vert[fa[1]]))
##              stl.append("            vertex %f %f %f\n"%tuple(vert[fa[2]]))
##              stl.append("        endloop\n")
##              stl.append("    endfacet\n")

##          return stl

    def __init__(self):
        OutputNode.__init__(self)
        self.Transformable = Transformable()
        self.sphereQuality = -1     # set in getSTL()
        self.cylinderQuality = -1  # set in getSTL()
        
        
    def getSTL(self, root, filename, sphereQuality=0,
               cylinderQuality=0):

        self.sphereQuality = sphereQuality
        self.cylinderQuality = cylinderQuality
        self.output = []
        self.output.append("solid %s\n"%filename)
        self.loopGeoms(root)
        self.output.append("endsolid %s\n"%filename)
        return self.output


    def NodeRepr(self, geom, matrix):
        # called recursively from loopGeom()
        if not self.checkGeom(geom): return

        # different geoms have to be treated different
        
        # IndexedGeoms are saved
        if isinstance(geom, IndexedPolygons):
            self.doitIndexedGeoms(geom, matrix)

        elif isinstance(geom, Spheres):
            self.doitSpheres(geom, matrix)

        elif isinstance(geom, Cylinders):
            self.doitCylinders(geom, matrix)

        # Lines are obviously not supported (no volume)
        elif isinstance(geom, IndexedPolylines):
            return


    def doitIndexedGeoms(self, geom, matrix):
        # gets called from NodeRepr()

        vert = geom.getVertices()
        faces = geom.getFaces()

        if faces.shape[1]==3:
            fn = TriangleNormals( vert, faces, 'PER_FACE')
        else:
            fn = TriangleNormals( vert, faces[:,:3], 'PER_FACE')
        
        #if geom is quads, or higher, make triangles
        if faces.shape[1] > 3:
            newfaces=[]
            newfn = []
            n=0 # counter for normals
            for face in faces:

                for i in range(len(face)):
                    # stop, if face contains -1 (for example, BSP Tree objects
                    # might 'fill up' faces with -1 to make uniform arrays
                    if face[i] == -1:
                        break
                    newfaces.append( [face[0],face[i-1],face[i]])
                    newfn.append(fn[n])

                # increment normal counter
                n = n + 1

            # finally, set new faces
            faces = newfaces
            fn = numpy.array(newfn).astype('f')

        self.doit(vert, faces, fn, matrix, invertNormals=geom.invertNormals)
        

    def doitSpheres(self, geom, matrix):
        geom = geom.asIndexedPolygons(quality=self.sphereQuality)
        vert = geom.getVertices()
        faces = geom.getFaces()
        fn = TriangleNormals( vert, faces, 'PER_FACE')
        self.doit(vert, faces, fn, matrix, invertNormals=geom.invertNormals)
 

    def doitCylinders(self, geom, matrix):
        geom = geom.asIndexedPolygons(quality=self.cylinderQuality)
        vert = geom.getVertices()
        faces = geom.getFaces()
        fn = TriangleNormals( vert, faces, 'PER_FACE')
        self.doit(vert, faces, fn, matrix, invertNormals=geom.invertNormals)


    def doit(self, vert, faces, fn, matrix, invertNormals=False):
        if invertNormals:
            fn = -1.0*fn
            l=[]
            for face in faces:
                l.append([face[0],face[2],face[1]])
            faces = numpy.array(l)

        coords = numpy.array(vert).astype('f')
        matrix = numpy.array(matrix).astype('f')
        matrix.shape = (4,4) # this changed after python2.3!!
        one = numpy.ones( (coords.shape[0], 1), \
                            coords.dtype.char )
        c = numpy.concatenate( (coords, one), 1 )

        # apply the matrix to the vertices
        newCoords = numpy.dot(c, \
                            numpy.transpose(matrix))[:, :3]

        self.output.extend(self.makestl(newCoords, faces, fn))


    def makestl(self, vert, faces, fn):
        stl = []
        
        for i in xrange(len(faces)):
            stl.append("    facet normal %f %f %f\n"%tuple(fn[i]))
            stl.append("        outer loop\n")
            fa = faces[i]
            stl.append("            vertex %f %f %f\n"%tuple(vert[fa[0]]))
            stl.append("            vertex %f %f %f\n"%tuple(vert[fa[1]]))
            stl.append("            vertex %f %f %f\n"%tuple(vert[fa[2]]))
            stl.append("        endloop\n")
            stl.append("    endfacet\n")

        return stl


class OutputVRML2(OutputNode):
    """ generates a list of strings describing DejaVu2 geoms in VRML 2.0
        format. Usage: OutputVRML2.getVRML2(geom, complete=0/1, normals=0/1)
        geom represents the geoms to be converted
        complete=1 will add vrml2 header and footer
        normals=1 will add normals
        """

    def __init__(self):
        OutputNode.__init__(self)
        self.geomName = None
        self.Transformable = Transformable()
        self.completeFlag  = 1     # is set in getVRML2()
        self.saveNormalsFlag   = 0 # is set in getVRML2()
        self.colorPerVertexFlag = True # set in getVRML2()
        self.usePROTOFlag    = 0   # if set, instanceMatrice geoms are not
                                   # saved as geoms, but re-used with PROTO 
        self.sphereQuality = -1     # default quality for sphere subsampling
        self.cylinderQuality = -1  # default quality for cyl. subsampling
        
        
    def getVRML2(self, root, complete=1, normals=0, 
                 colorPerVertex=True,
                 usePROTO=0, sphereQuality=0, cylinderQuality=0):
        """ this method returns a list of strings describing DejaVu2 Geoms in
        VRML2 format """

        # this is the method the user should call
        
        self.completeFlag = complete # if 1, header and footer will be added
        self.saveNormalsFlag = normals # if 1, normals are added
        self.colorPerVertexFlag = colorPerVertex # if False, color per face
        self.usePROTOFlag = usePROTO # if set to 1, instance geoms are saved
                                     # with one PROTO. Else: data for these
                                     # geoms will be prepared which blows up
                                     # the file size

        self.sphereQuality = sphereQuality # default is 0
        self.cylinderQuality = cylinderQuality # default is 0

        self.output = []
        
        if self.completeFlag:
            self.output.extend(self.getFileHeader1())
            self.output.extend(self.getCopyright())
            self.output.extend(self.getFileHeader2())
            
        self.loopGeoms(root)

        if self.completeFlag:
            self.output.extend(self.getFileTrailer())

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
            sphGeom.instanceMatricesFortran = geom.instanceMatricesFortran
            sphGeom.instanceMatricesFortranIndex = geom.instanceMatricesFortranIndex
            sphGeom.VRML2CreatedPROTOForThisGeom = 0
            self.output.extend( self.doit(sphGeom, matrix) )

        elif isinstance(geom, Cylinders):
            # convert Cylinders geom into IndexedPolygons
            cylGeom = geom.asIndexedPolygons(quality=self.cylinderQuality) 
            cylGeom.viewer = geom.viewer
            cylGeom.parent = geom.parent
            cylGeom.name = geom.name
            cylGeom.fullName = geom.fullName
            cylGeom.instanceMatricesFortran = geom.instanceMatricesFortran
            cylGeom.instanceMatricesFortranIndex = geom.instanceMatricesFortranIndex
            cylGeom.VRML2CreatedPROTOForThisGeom = 0
            self.output.extend( self.doit(cylGeom, matrix) )
            
        elif isinstance(geom, GlfLabels):
            # convert Cylinders geom into IndexedPolygons
            glfGeom = geom.asIndexedPolygons() 
            glfGeom.viewer = geom.viewer
            glfGeom.parent = geom.parent
            glfGeom.name = geom.name
            glfGeom.fullName = geom.fullName
            glfGeom.instanceMatricesFortran = geom.instanceMatricesFortran
            glfGeom.instanceMatricesFortranIndex = geom.instanceMatricesFortranIndex
            glfGeom.VRML2CreatedPROTOForThisGeom = 0
            self.output.extend( self.doit(glfGeom, matrix) )


    def doit(self, geom, matrix):
        # gets called from NodeRepr()
        vrml2 = []
        
        # if self.usePROTO is set to 1: don't convert instance matrices
        # geoms into geoms, but use the vrml2 USE 
        if self.usePROTOFlag:
            # create all the necessary data to be put in the PROTO which goes
            # in the header of the vrml2 file
            if geom.VRML2CreatedPROTOForThisGeom == 0:
                name = self.getGeomName(geom)
                vrml2.append("PROTO "+name+" [ ] {\n")
                identityMatrix = numpy.identity(4).astype('f')
                vrml2.extend( self.doitReally(geom, identityMatrix) )
                vrml2.append("}\n")

                # now insert this into the header of self.output
                for i in range(len(vrml2)):
                    self.output.insert(i+1,vrml2[i])
                geom.VRML2CreatedPROTOForThisGeom = 1 # don't add this geom
                                                      # to the header next time
                # this PROTO flag will be deleted when we leave the recursive
                
            # and add it as USE to the body
            vrml2 = [] # clear the list because this data has been added
            vrml2.extend( self.doitUsePROTO(geom, matrix) )
            return vrml2

        else:
            # here we save all the data for all geoms
            return self.doitReally(geom, matrix)


    def doitReally(self, geom, matrix):
        # add header for geom
        vrml2 = []

        vrml2.extend(self.getGeomHeader(geom))
        vrml2.extend(self.getShape())
        vrml2.extend(self.getAppearance())
        mat, colors = self.getMaterial(geom)
        vrml2.extend(mat)

        # add texture if applicable:
        if geom.texture:
            vrml2.extend( self.getTexture(geom) )

        vrml2.append("        }\n")
        vrml2.append("\n")

        # add coords, faces, etc
        vrml2.extend( self.getGeomDescr(geom) )

        # add texCoord Coordinates is applicable
        if geom.texture:
            vrml2.extend( self.getTexCoords(geom) )

        # add texCoordsIndex if applicable
        if hasattr(geom.vertexSet,'texCoordsIndex'):
             vrml2.extend( self.getTexCoordsIndex(geom))
        
        # add normals if applicable
        if self.saveNormalsFlag and isinstance(geom, IndexedPolygons):
            vrml2.extend( self.getNormals(geom) )

        # add colors per vertex if applicable
        if colors is not None and len(colors): 
            vrml2.extend( self.getColors(geom, colors) )

        # add closing brackets for geom
        vrml2.append("            }\n")
        vrml2.append("         }\n")

        # add transformations for geom
        vrml2.extend( self.getTransforms(matrix) )

        # add closing bracket for Transform{}
        vrml2.append("      }\n")
        return vrml2
    

    def doitUsePROTO(self, geom, matrix):
        # FIXME: this works currently only with geoms that are not grouped
        # i.e. it doesnt work with secondary structures, they will be saved
        # as PROTO too, but also for each instanceMatrix (->redundant)

        vrml2 = []
        
        geom.instanceMatricesFortranIndex = 0
        name = string.split(self.getGeomHeader(geom)[0])[1]
        vrml2.append("    Transform {\n")
        vrml2.append("      children  "+name+" { }\n")

        # add transformations for geom
        vrml2.extend( self.getTransforms(matrix) )

        # add closing bracket for Transform{}
        vrml2.append("      }\n")
        return vrml2
    


    def getFileHeader1(self):
        vrml2=[]
        vrml2.append("#VRML V2.0 utf8 Python Molecular Viewer Geom\n")
        vrml2.append("\n")
        return vrml2


    def getCopyright(self):
        vrml2=[]
        vrml2.append("WorldInfo {\n")
        
        vrml2.append('     title ""\n')
        vrml2.append("     info [\n")
        vrml2.append('         "Copyright (c) 2002 D. Stoffler, M.F. Sanner and A.J. Olson"\n')
        vrml2.append('         "Molecular Graphics Lab"\n')
        vrml2.append('         "The Scripps Research Institute, La Jolla, CA"\n')
        vrml2.append('         "VRML2 file generated with the Python Molecular Viewer:"\n')
        vrml2.append('         "http://www.scripps.edu/~sanner/python/pmv/"\n')
        vrml2.append("     ]\n")
        vrml2.append("}\n")
        return vrml2


    def getFileHeader2(self):
        vrml2=[]
        vrml2.append("Group {\n")
        vrml2.append("  children    [\n")
        return vrml2


    def getFileTrailer(self):
        vrml2=[]
        vrml2.append("  ]\n")
        vrml2.append("}\n")
        return vrml2


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
                   string.strip(str(g.instanceMatricesFortranIndex))+ '_'
            g = g.parent
        return name


    def getGeomHeader(self, geom):
        # generates geom name
        vrml2=[]
        g = geom
        name = self.getGeomName(geom)
        vrml2.append("  DEF "+name+" Transform {\n")
        return vrml2


    def getShape(self):
        vrml2=[]
        vrml2.append("      children    Shape {\n")
        return vrml2


    def getAppearance(self):
        vrml2=[]
        vrml2.append("        appearance        Appearance {\n")
        return vrml2


    def getMaterial(self, geom):
        vrml2=[]

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

        vrml2.append("          material        Material {\n")
        vrml2.append("            ambientIntensity      "+ambInt+"\n")
        vrml2.append("            diffuseColor          "+difCol+"\n")
        vrml2.append("            emissiveColor         "+emCol+"\n")
        vrml2.append("            specularColor         "+specCol+"\n")
        vrml2.append("            shininess             "+shin+"\n")
        vrml2.append("            transparency          "+trans+"\n")
        vrml2.append("          }\n")
        return vrml2, colors


    def getGeomDescr(self, geom):
        vrml2 = []

        if isinstance(geom, IndexedPolygons):
            vrml2.append("        geometry        IndexedFaceSet {\n")
            # add vertices
            vrml2.extend( self.getCoords(geom) )
            # add face indices
            vrml2.extend( self.getFaces(geom) )
            # add color indices if applicable
            if geom.materials[GL.GL_FRONT].colorIndex is not None and len(geom.materials[GL.GL_FRONT].colorIndex):
                vrml2.extend( self.getColorIndex(geom) )
                    
        elif isinstance(geom, IndexedPolylines):
            vrml2.append("        geometry        IndexedLineSet {\n")
            # add vertices
            vrml2.extend( self.getCoords(geom) )
            # add face indices
            vrml2.extend( self.getFaces(geom) )

        return vrml2


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
        vrml2=[]
        vertices = geom.getVertices()
        vrml2.append("          coord               Coordinate {\n")
        vrml2.append("            point        [\n")
        for vert in vertices:
            vstring = "                          "+'%.5f'%vert[0]+" "+\
                      '%.5f'%vert[1]+" "+'%.5f'%vert[2]+",\n"
            vrml2.append(vstring)
        vrml2.append("                         ]\n")
        vrml2.append("	  }\n")
        return vrml2


    def getFaces(self, geom):
        #print "getFaces"
        vrml2=[]
        faces = geom.getFaces()
        vrml2.append("	  coordIndex	[\n")

        for face in faces:
            facestring = "                          "
            if geom.invertNormals: # reverse faces
                facestring = facestring + `face[0]` + ", "
                for i in range(len(face)-1,0,-1): 
                    facestring = facestring + `face[i]` + ", "
            else:
                for f in face:
                    facestring = facestring + `f` + ", "

            facestring = facestring + "-1,\n"
            vrml2.append(facestring)
        vrml2.append("	                ]\n")
        return vrml2


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
        else:
            fn = geom.normals

        vrml2=[]
        vrml2.append("	  normal Normal	{\n")
        vrml2.append("             vector [\n")
        for n in fn:
            vrml2.append("                %.5f"%n[0]+" %.5f"%n[1]+\
                         " %.5f"%n[2]+" \n")
        vrml2.append("	                  ]\n")
        vrml2.append("                  }\n")
        return vrml2


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


    def getColors(self, geom, colors):
        vrml2=[]
        vrml2.append("\n")
        vrml2.append("     colorPerVertex %s\n"%string.upper(
            str(self.colorPerVertexFlag))) # TRUE or FALSE ... capital letters
        vrml2.append("     color Color {\n")
        vrml2.append("            color [\n")
        for c in colors:
            cstring = '                  %.3f'%c[0]+" "+'%.3f'%c[1]+\
                      " "+'%.3f'%c[2]+",\n"
            vrml2.append(cstring)
        vrml2.append("                  ]\n")
        vrml2.append("                }\n")
        return vrml2


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


class DatGUI:
    """ basic gui for DataOutput """
    
    def __init__(self, master=None, title=None):
        self.master = master
        self.root = None
        self.title = title
        if self.title is None:
            self.title = 'Options Panel'

        if self.master is None:
            self.master = Tkinter.Frame()
            self.master.pack()
   
        self.sphInput       = Tkinter.StringVar()
        self.cylInput       = Tkinter.StringVar() 

        self.sphereQuality = -1    # values range from -1 to 5
        self.cylinderQuality = -1 # values range from -1 to 5

        self.readyToRun = 0       # set to 1 in OK_cb, set to 0 in Cancel_cb

        self.idf = InputFormDescr(title=self.title)


    def sphereQuality_cb(self, event=None):
        val = self.sphInput.get()
        if len(val) == 0 or val is None:
            val = self.sphereQuality
        try:
            val = int(val)
            if val < -1:
                val = -1
            elif val > 5:
                val = 5
            self.sphereQuality = val
            self.sphInput.set(str(self.sphereQuality))
        except ValueError:
            self.sphInput.set(str(self.sphereQuality))


    def cylinderQuality_cb(self, event=None):
        val = self.cylInput.get()
        if len(val) == 0 or val is None:
            val = self.cylinderQuality
        try:
            val = int(val)
            if val < -1:
                val = -1
            elif val > 5:
                val = 5
            self.cylinderQuality = val
            self.cylInput.set(str(self.cylinderQuality))
        except ValueError:
            self.cylInput.set(str(self.cylinderQuality))


    def OK_cb(self):
        self.readyToRun = 1
        self.sphereQuality_cb()
        self.cylinderQuality_cb()
        self.master.grab_release()
        self.master.quit()
        self.optionsForm.withdraw()


    def Cancel_cb(self):
        self.readyToRun = 0
        self.master.grab_release()
        self.master.quit()
        self.optionsForm.withdraw()
        

    def displayPanel(self, create):
        self.readyToRun = 0
        if create == 0:
            self.optionsForm.deiconify()
        else:
            self.optionsForm = InputForm(self.master, self.root,
                                         descr=self.idf,okcancel=0)
        # grab the focus, i.e. the program stops until OK or Cancel is pressed
        self.master.grab_set()
        self.master.mainloop()


    def getValues(self):
        vals = {}
        vals['sphereQuality'] = self.sphereQuality
        vals['cylinderQuality'] = self.cylinderQuality
        return vals
        

class STLGUI(DatGUI):
    """ this is the gui for OutputSTL
    - Save normals adds normals to all geoms to be saved as vrml2
    - Invert normals inverts the normals for all geoms to be saved as vrml2
    - Sphere quality is the subsampling of the spheres. Default value is 2
      lowest allowed value is 0. Higher values increase the file size
      significantly.
    - Cylinder quality is the subsampling of cylinders. Default value is 10
      lowest allowed value is 3."""
 
    def __init__(self, master=None, title=None):
        DatGUI.__init__(self, master, title)

        self.createForm()

    def createForm(self):

        self.idf.append({'widgetType':Tkinter.Label,
                 'wcfg':{'text':'Sphere quality'},
                 'gridcfg':{'sticky':'w','columnspan':1, 'row':0, 'column':0},
                 })

        self.idf.append({'name':'inpSphQual',
                         'widgetType':Tkinter.Entry,
                         'defaultValue':self.sphereQuality,
                         'wcfg':{'font':(ensureFontCase('Courier'),10),
                                 'width':5,'textvariable':self.sphInput,
                                 'command':self.sphereQuality_cb,
                                 },
                         'gridcfg':{'sticky':'e',
                                    'columnspan':1, 'row':0, 'column':1 }
                         })

        self.idf.append({'widgetType':Tkinter.Label,
                 'wcfg':{'text':'Cylinder quality'},
                 'gridcfg':{'sticky':'w','columnspan':1, 'row':1, 'column':0},
                 })

        self.idf.append({'name':'inpCylQual',
                         'widgetType':Tkinter.Entry,
                         'defaultValue':self.cylinderQuality,
                         'wcfg':{'font':(ensureFontCase('Courier'),10),
                                 'width':5,'textvariable':self.cylInput,
                                 'command':self.cylinderQuality_cb,
                                 },
                         'gridcfg':{'sticky':'e',
                                    'columnspan':1, 'row':1, 'column':1 }
                         })
        
	self.idf.append({'widgetType':Tkinter.Button,
                         'text':'OK',
                         'wcfg':{},
                         'gridcfg':{'sticky':'wens',
                                    'columnspan':1, 'row':2, 'column':0},
                         'command': self.OK_cb})


	self.idf.append({'widgetType':Tkinter.Button,
                         'text':'Cancel',
                         'wcfg':{},
                         'gridcfg':{'sticky':'wens',
                                    'columnspan':1, 'row':2, 'column':1},
                         'command': self.Cancel_cb})


    def invertNormals_cb(self):
        pass


class VRML2GUI(DatGUI):
    """This is the gui for OutputVRML2:
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


