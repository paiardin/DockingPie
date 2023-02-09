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

#########################################################################
#
# Date: DEC 2003  Author: Daniel Stoffler
#
#       stoffler@scripps.edu
#
# Copyright: Daniel Stoffler and TSRI
#
# revision: Guillaume Vareille
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/DejaVu2/DataInput.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#
# $Id: DataInput.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#

import string
import os
from struct import unpack
import numpy

try:
    from geomutils.geomalgorithms import detectObjFileContent, readObjFileGroup
except:
    pass

from DejaVu2.IndexedPolygons import IndexedPolygons
from DejaVu2.Materials import Materials

def pythonStringFromCString(aCString):
    lPythonString = ''
    for c in aCString:
        if c == '\x00' or c == '':
            break
        lPythonString += c
    assert c == '\x00' or c =='', c
    return lPythonString


def readOBJ(filename):
    if (filename is None) or (filename == ''):
        return None

    lGroupNames = numpy.zeros(128*256, typecode = "c")
    lGroupNames.shape = (128,256)    
    lLibNames = numpy.zeros(128*256, typecode = "c")
    lLibNames.shape = (128,256)    
    lMaterialNames = numpy.zeros(128*256, typecode = "c")
    lMaterialNames.shape = (128,256)    
    out = detectObjFileContent(filename, lGroupNames, lLibNames, lMaterialNames)
    lGeoms = []
    if out[0] is True:
        lMaterialNames = lMaterialNames[:out[3]]
        lNumOfLibs = out[2]

        materialDict = {'default': Materials()}
        for i in range(len(lMaterialNames)):
            materialName = pythonStringFromCString(lMaterialNames[i])
            materialDict[materialName] = Materials()        

        libNames = []
        for i in range(lNumOfLibs):
            libNames.append(pythonStringFromCString(lLibNames[i]))
            materialDict.update(readMTL(libNames[-1]))

        materialList = []
        for i in range(len(lMaterialNames)):
            materialName = pythonStringFromCString(lMaterialNames[i])
            materialList.append(materialDict[materialName])          

        lNumOfGeoms = out[1]
        for i in range(lNumOfGeoms):
            lGroupName = pythonStringFromCString(lGroupNames[i])
            out = readObjFileGroup(filename, lGroupName, None,
                                   None, None, None, None, None )
            if out[0] is True and out[2] > 0 and out[3] > 0:
                vertices = numpy.zeros(out[2]*3, typecode = "f")
                vertices.shape = (out[2],3)
                faces = numpy.zeros(out[3]*3, typecode = "i")
                faces.shape=(out[3],3)
                textureVertices = numpy.zeros(out[4]*2, typecode = "f")
                textureVertices.shape = (out[4],2)
                textureFaces = numpy.zeros(out[5]*3, typecode = "i")
                textureFaces.shape=(out[5],3)
                if ( len(libNames) > 0 ) and ( len(lMaterialNames) > 1 ):
                    triangleMaterialIndices = numpy.zeros(len(faces)*1, typecode = "i")
                else:
                    triangleMaterialIndices = None
                out = readObjFileGroup(filename, lGroupName, lMaterialNames,
                                       vertices, faces, 
                                       textureVertices, textureFaces,
                                       triangleMaterialIndices
                                      )

                if out[0] is True:
                    if len(textureVertices)>0:
                        geom = IndexedPolygons(vertices=vertices,
                                       faces=faces,
                                       textureCoords=textureVertices)
                    else:
                        geom = IndexedPolygons(vertices=vertices,
                                               faces=faces)

                    if triangleMaterialIndices is not None:
                        breaked = False
                        for triangleMaterialIndex in triangleMaterialIndices:
                            if triangleMaterialIndex != triangleMaterialIndices[0]:
                                breaked = True
                                break
                        if breaked is False:
                            if triangleMaterialIndex == 0:
                                # use inherited
                                pass
                            else:
                                geom.Set(materials=[materialList[triangleMaterialIndex].prop[1],],
                                     propName='diffuse',
                                     inheritMaterial=False)
                        else:
                            triangleMaterials = []
                            for triangleMaterialIndex in triangleMaterialIndices:
                                 triangleMaterials.append(materialList[triangleMaterialIndex].prop[1])  
                            geom.Set(materials=triangleMaterials,
                                     propName='diffuse',
                                     inheritMaterial=False)
                               


                    if (lGroupName == 'default'):
                        nameSplit = os.path.splitext(filename)
                        geom.name = os.path.split(nameSplit[0])[-1]
                    else:
                        geom.name = lGroupName
                    #print "geom.name", geom.name
                    lGeoms.append(geom)
        return lGeoms
    else:
        return None


def readMTL(filename):
    if (filename is None) or (filename == ''):
        return None
    file = open(filename, "r")

    materials = {}
    material = None
    for line in file:
        if line.startswith('#'):
            continue
        tokens = line.split()
        if not tokens:
            continue
        if tokens[0] == 'newmtl':
            material = Materials()
            materials[tokens[1]] = material
        elif material is None:
            warnings.warn("missing newmtl statement in mtl file")
        elif tokens[0] == 'Ka':
            material.Set(ambient=(eval(tokens[1]),eval(tokens[2]),eval(tokens[3]),1.))
        elif tokens[0] == 'Kd':
            material.Set(diffuse=(eval(tokens[1]),eval(tokens[2]),eval(tokens[3]),1.))
        elif tokens[0] == 'Ks':
            material.Set(specular=(eval(tokens[1]),eval(tokens[2]),eval(tokens[3]),1.))
        elif tokens[0] == 'Ns':
            material.Set(shininess=eval(tokens[1]))
        elif tokens[0] == 'd' or tokens[0] == 'Tr':
            material.Set(opacity=eval(tokens[1]))
        elif tokens[0] == 'map_Kd':
            # load the texture map
            pass

    return materials


def readAnySTL(filename):
    if (filename is None) or (filename == ''):
        return

    nameSplit = os.path.splitext(filename)
    name = os.path.split(nameSplit[0])[-1]

    P = ReadASCIISTL()
    l = P.read(filename)
    geoms = P.doit(l)
    if len(geoms) > 0:
        geoms[0].name = name
        for i in range(1, len(geoms)):
            geoms[i].name = name + '+' + str(i)
        return geoms

    P = ReadBinarySTL()
    geom = P.doit(filename)
    if geom is not None:
        geom.name = name
        geoms = [geom] 
        return geoms
    else:
        return None



class ReadASCIISTL:
    """This class parses a ASCII STL file and converts it into a DejaVu2
IndexedPolygon geometry.
Usage:

    R = ReadASCIISTL()         # instanciate object
    lines = R.read(filename)   # read ASCII STL file
    geoms = R.doit(lines)      # parse file, build geoms
"""
    

    def __init__(self):
        self.geoms = []  # list storing all Indexed Geoms that will be built
                         # in the doit() method


    def read(self, filename):
        """For your convenience, a read lines method"""
        f = open(filename, "r")
        lines = f.readlines()
        f.close()
        return lines


    def doit(self, lines):
        """Parse ASCII STL files, build DejaVu2 IndexedPolygons and return them
        """
        
        self.geoms = []
        faces = []
        vertices = []
        normals = []
        name = 'IndexedGeom'
        faceIndex = 0 # counter
        

        # parse ASCII STL
        for line in lines:
            spl = string.split(string.lower(line))
            if len(spl)==0: continue
            if spl[0] == 'solid':
                name = spl[1]
                continue

            elif spl[0] == 'facet' and spl[1] == 'normal':
                normals.append( [float(spl[2]), float(spl[3]), float(spl[4])] )
                continue

            elif spl[0] == 'vertex':
                vertices.append( [float(spl[1]), float(spl[2]), float(spl[3])])
                continue

            elif spl[0] == 'endfacet':
                faces.append( [faceIndex, faceIndex+1, faceIndex+2] )
                faceIndex = faceIndex + 3
                continue

            elif spl[0] == 'endsolid':
                geom = IndexedPolygons(name, vertices=vertices, faces=faces,
                                       fnormals = normals)
                self.geoms.append(geom)
                faces = []
                vertices = []
                normals = []
                name = 'IndexedGeom'
                faceIndex = 0 # counter
                continue

            else:
                continue
            
        return self.geoms


class ReadBinarySTL:

    # (c) Daniel Stoffler, Biozentrum, University of Basel, Switzerland,
    # August 2005
    
    """This class parses a Binary STL file and converts it into a DejaVu2
IndexedPolygon geometry.
Usage:

    R = ReadBinarySTL()        # instanciate object
    geoms = R.doit(filename)   # parse file, build geoms
"""
    

    def doit(self, filename):
        """Parse ASCII STL files, build DejaVu2 IndexedPolygons and return them
        """
        # open file to read in binary mode
        f = open(filename, 'rb')
        # read the 84 bytes of the header
        header = f.read(84)

        # get the comment
        comment = unpack("80c", header[0:80])

        # get the total number of facets
        lenFaces = unpack("1i", header[80:84])[0]
                
        data = f.read() # read to end
        f.close()

        # organize the data
        normals = []
        vertices = []
        faces = []

        j = 0 # counter for reading binary data
        f = 0 # counter for generating faces

        for i in range(lenFaces):
            normal = unpack("3f", data[j:j+12])
            normals.append(normal)
            faces.append((f, f+1, f+2))
            verts = unpack("9f", data[j+12:j+48])
            vertices.append(verts[0:3])
            vertices.append(verts[3:6])
            vertices.append(verts[6:9])
            #spacer = unpack("2c", data[j+48:j+50])

            # increment binary counter
            j = j + 50
            # increment faces counter
            f = f + 3

        name = "IndexedGeom"
        geom = IndexedPolygons(name, vertices=vertices, faces=faces,
                               fnormals = normals)
        return geom


    def info(self):
        txt = """Binary STL files consist of a 80 byte header line that can be
interpreted as a comment string. The following 4 bytes interpreted as a long
integer give the total number of facets. What follows is a normal and 3
vertices for each facet, each coordinate represented as a 4 byte floating
point number (12 bytes in all). There is a 2 byte spacer between each facet.
The result is that each facet is represented by 50 bytes, 12 for the normal,
36 for the 3 vertices, and 2 for the spacer. 
""" 

        return txt
        




if __name__ == '__main__':
#    P = ReadASCIISTL()
#    l = P.read('test.stl')
#    geoms = P.doit(l)
#    from DejaVu2 import Viewer
#    vi = Viewer()
#    vi.AddObject(geoms[0])
    #geoms = readOBJ('sample.obj')
    geoms = readAnySTL('test.stl')
    print "len(geoms)", len(geoms)
    from DejaVu2 import Viewer
    vi = Viewer()
    for g in geoms:
        vi.AddObject(g)

    #import pdb;pdb.set_trace()

        
