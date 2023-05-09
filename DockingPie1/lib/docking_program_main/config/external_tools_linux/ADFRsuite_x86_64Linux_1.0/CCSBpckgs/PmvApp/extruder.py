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
# Author: Sophie COON, Kevin CHAN, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2014
#
#############################################################################

# $Header: /mnt/raid/services/cvs/PmvApp/extruder.py,v 1.2.4.1 2017/07/13 20:55:28 annao Exp $
#
# $Id: extruder.py,v 1.2.4.1 2017/07/13 20:55:28 annao Exp $
#

import numpy
import math
from PmvApp import Ribbon
from DejaVu2.Shapes import Rectangle2D

from opengltk.OpenGL import GL
#from opengltk.extent.utillib import glTriangleNormals
from geomutils.geomalgorithms import TriangleNormals
# Gotta move this class in DejaVu2 not dependent on MolKit2 any longer

class Sheet2D:
    """ Class implementing a set of method to compute a path3D, normals,
    binormals and transformation matrices given a 2 sets of control points
    coordinates. (ctrl points and torsion ctrl points)
    """
    def compute(self, coords , isHelix, nbrib = 2,
                nbchords = 10,  width = 1.5, offset = 1.2, off_c = 0.5):
        """ """
        self.nrib = nbrib
        self.coords = coords
        natoms = len(self.coords)
        self.width = width
        self.chords = nbchords
        self.offset = offset
        self.isHelix = isHelix
        # Create a ribbon2D object
        self.smooth = Ribbon.ribbon2D(nbrib, width, nbchords,
                                      offset, natoms, coords,
                                      isHelix, off_c)
        self.oldsmooth = self.smooth
        self.addFirstLastPoints()
        self.verts2D_flat = numpy.array(numpy.reshape(self.smooth,
                                                          (-1,4))[:,:3])
        
        path = (self.smooth[0,:,:3] + self.smooth[1,:,:3])*0.5
        self.path = path.astype('f')
        self.faces2D = self.computeFaces2D()
        self.binormals = self.computeBinormals()
        self.normals = self.computeNormals()
        self.matrixTransfo= self.buildTransformationMatrix(
            self.path, self.normals, self.binormals)

    def interpolateSmoothArray(self, first, last, nPts):
        """Insert self.chord-1 points into self.smooth such that each residue
        is represented by exactly self.nchords points"""

        # create array to hold new points
        beg = numpy.ones( (self.nrib, nPts, 4) ).astype('f')

        # set all points to the coordinates of first point in each rib
        beg = beg * numpy.reshape( self.smooth[:,first], (self.nrib,1,4) )
        # compute vectors from first to second point for each rib
        # and divide by number of intervals
        vec =  (self.smooth[:,last] -  self.smooth[:,first]) / nPts

        # resize vector
        vec1 = numpy.resize( vec[0], (nPts, 4) )
        vec2 = numpy.resize( vec[1], (nPts, 4) )
        vec = numpy.concatenate( (vec1, vec2) )
        vec.shape = (self.nrib, nPts, 4)

        scale = numpy.arange(nPts)
        scale.shape = (1,nPts,1)
        
        # add it 0,1,2,3 ... nchord-1 times to the new points
        return beg + vec*scale
    
    def addFirstLastPoints(self):
        last = self.smooth.shape[1]
        
        ar2 = self.interpolateSmoothArray(last-2, last-1, self.chords)
        self.smooth = numpy.concatenate((self.smooth[:, :-2, :],
                                           ar2, self.smooth[:, -1:, :] ),
                                          1)
        
    def computeFaces2D(self):
        f = []
        n = self.smooth.shape[1]
        faces2D = numpy.array([(x,x+n,x+n+1,x+1) for x in range(n-1)]).astype('i')
        return faces2D


    def computeBinormals(self):
        f = numpy.array(self.faces2D[:,:3])
        binorm = TriangleNormals(self.verts2D_flat,
                                      f, 'PER_VERTEX')
        binorm[self.smooth.shape[1]-1,:] = binorm[self.smooth.shape[1]-2,:]
        binorm1 = binorm[:self.smooth.shape[1]]
        return binorm1

    def computeNormals(self):
        normals = self.smooth[1,:,:3] - self.smooth[0,:,:3]
        normals = normals/self.width
        normals = normals.astype('f')
        return normals

    def buildTransformationMatrix(self, path, normals, binormals):
        matrixTransfo = numpy.zeros( (len(path), 3, 3) ).astype('f')
        for i in range(len(path)):
            matrixTransfo[i][0] = normals[i]
            matrixTransfo[i][1] = binormals[i]
            matrixTransfo[i][2] = path[i]
        return matrixTransfo


def getAverage(elements):
    """Function to get the average of a list of elements"""
    average = numpy.array([0,0,0])
    for el in elements:
        average = average + numpy.array(el)
    average = average/len(elements)
    average = list(average)
    return average



######################################################################
###                                                                ###
###                 CLASSES TO DEFINE EXTRUSION                    ###
###                                                                ###
######################################################################

class ExtrudeObject:
    """
    Base class to take a shape and extrude it along a 3D path.
    """
    
    def __init__(self, path3D, matrix, shape, cap1=0, cap2=0, arrow=0,
                 larrow=3, warrow=2):
        """
        Constructor: Takes as arguments a path3D, matrix, 2D shape, and
        optional cap1, cap2, and arrow.  Calls getextrudeVertices() and
        getFaces() and stores the return values in self.vertices, self.vnormals
        , and self.faces.
        """
        self.path3D = path3D
        self.matrix = matrix
        self.shape = shape
        self.arrow = arrow
        if not isinstance(self.shape, Rectangle2D):
            self.arrow = 0
        self.larrow = larrow
        self.warrow = warrow
        self.cap1 = cap1
        self.cap2 = cap2
        if self.arrow: self.cap2 = 0    # no end cap if arrow
        self.norms = matrix[:,0]
        self.vertices, self.vnormals = self.getextrudeVertices()
        self.faces = self.getFaces()

        # Here need to create a numeric array for each properties:
        # colors
        self.colors = numpy.ones((len(self.faces), 3),'f')
        # opacities
        self.opacities = numpy.ones((len(self.faces),),'f')


    def getextrudeVertices(self):
        """ get the extruded vertices and normals.  If an arrow is specified,
        extra vertices are added.  If caps are specified,
        the end vertices must be duplicated and used as extra vertices."""
        ls = self.shape.lenShape
        length = len(self.path3D)

        f = g = 0
        # f adds another section, a section is ls vertices.
        # g adds the center pt of cap
        if self.cap1: f, g = f+1, g+1   # add a section and a middle point.
        if self.cap2: f, g = f+1, g+1
        if self.arrow:f = f+3           # add 3 sections:
                                        # - 1 junction between the arrow and
                                        # the rest of the geom.
                                        # - 2 for the caps of the arrow.

        # get the right size pts and normals arrays for caps and arrows
        pts = numpy.zeros( ((length+f)*ls+g,3),'f')
        ptsNorm = numpy.zeros( ((length+f)*ls+g,3),'f')

        # extrude 2D shape along 3D path
        for i in range(length):
            pts[i*ls:i*ls+ls] = numpy.dot(self.shape.contpts,
                                            self.matrix[i]).astype('f')
            newNorm = numpy.dot(self.shape.contnorm,
                                             self.matrix[i]).astype('f')
            ptsNorm[i*ls:i*ls+ls] = newNorm

        # add 3*ls vertices for the arrow
        if self.arrow:
            pts = self.addArrow(pts)

        # add duplicate vertices and centers for cap1 and/or cap2 at the end
        # of the array of pts
        if self.cap1 and not self.cap2:
            pts[-ls-1:-1] = pts[0:ls]            
            pts[-1] = self.path3D[0]  # center of cap1
            # normals
            n = TriangleNormals(pts[-4:], [(3, 1, 0)])
            for k in range(ls+1):
                ptsNorm[-ls-1+k] = n[0]
                
        elif self.cap2 and not self.cap1:
            pts[-ls-1:-1] = pts[(length-1)*ls:length*ls]
            pts[-1] = self.path3D[-1]  # center of cap2
            # normals
            n = TriangleNormals(pts[-4:], [(0, 1, 3)])
            for k in range(ls+1):
                ptsNorm[-ls-1+k] = n[0]

        elif self.cap1 and self.cap2:
            pts[-2*ls-2:-ls-2] = pts[0:ls]   # vertices for cap1
            pts[-ls-2] = self.path3D[0] # center of cap1
            # normals
            n = TriangleNormals(pts[-5-ls:-ls-1], [(3, 1, 0)])
            for k in range(ls+1):
                ptsNorm[-2*ls-2+k] = n[0]

            pts[-ls-1:-1] =pts[(length-1)*ls:length*ls]  # vertices for cap2
            pts[-1] = self.path3D[-1]       # center of cap2
            # normals
            n = TriangleNormals(pts[-4:], [(0, 1, 3)])
            for k in range(ls+1):
                ptsNorm[-ls-1+k] = n[0]

        return pts, ptsNorm


    def addArrow(self, pts):
        """ Compute the vertices for an arrowhead.
        self.larrow  specifies how far back the arrow goes along the path3D.
        self.warrow is the width in real units.
        """
        #                                         @     @
        #                                        /|    /|
        #                                       / |   / |   
        #       *-----*-----*-----*-----*      *  |  /  |  
        #      /|    /|    /|    /|    /|      |  | /   |
        #     / |   / |   / |   / |   / |      |  @/    @          
        #    *-----*-----*-----*-----*  |   *  | //    /
        #    |  |  |  |  |  |  |  |  |  |  /|  |//    /
        #    |  *  |  *  |  *  |  *  |  * / |  */    /
        #    | /   | /   | /   | /   | / @  |  @    /
        #    |/    |/    |/    |/    |/  |  |  |   /
	#    *-----*-----*-----*-----*   |  *  |  /
	#                                | /   | /
        #                                |/    |/
        #                            |   @     @
	#   <    BLOC1->  lp-lenarrow>   |
	#                         <lp-lenarrow>
        #                           BLOC2      |
	#                                     <v[lp-lenarrow:] translated >
        #                                             BLOC3
        # lp = len(self.path3D)
        # *: vertices results of the extrusion of the shape2D along the path3D
        # @: those vertices translated .


        lenarrow = self.larrow
        widtharrow = self.warrow

        ls = self.shape.lenShape
        lv = len(pts)
        path = self.path3D
        lp = len(path)
        lnew = lp-lenarrow
        if lnew<0:
            self.arrow = 0
            # the arrow is too long:  no arrow added
            return pts

        # Compute the translation vector that will be used to translate the
        # vertices used to build the arrowhead.

        d = numpy.zeros( lenarrow+1, 'f') # the firs section will be
                                            # duplicate to build the caps.
        for n in xrange(1, lenarrow+1):
            d[n] = d[n-1] + self.dist(path[-n], path[-n-1])
        d = d*widtharrow/d[lenarrow] - abs(self.shape.contpts[0][0])

        # We apply the translation to the arrow vertices 
        # and reorganize them in the array of points so the array contains
        # three blocks of vertices( see the drawing above):
        # A total of 3*ls vertices have been added
        # Build bloc3
        for s in xrange(lenarrow+1):
            #v = [0, 0, 0]
            uvector = self.norms[-s-1]
            v = (uvector*d[s]).astype('f')
            
            for j in range(ls):
                if j < (ls)/2 or j >= ls:
                    pts[(lp-s+2)*ls+j] = pts[(lp-s-1)*ls+j]-v
                elif j >= (ls)/2 and j < ls:
                     pts[(lp-s+2)*ls+j] = pts[(lp-s-1)*ls+j]+v

        # build bloc2
        # duplicate vertices for corner of arrow
        pts[lnew*ls:(lnew+1)*ls] = pts[(lnew-1)*ls:lnew*ls]
        pts[(lnew+1)*ls:(lnew+2)*ls] = pts[(lnew+2)*ls:(lnew+3)*ls]
        return pts
            

    def dist(self, p1, p2):
        """ Calculate distance between two pts."""
        dx = p2[0]-p1[0]
        dy = p2[1]-p1[1]
        dz = p2[2]-p1[2]

        distance = math.sqrt(pow(dx, 2)+pow(dy, 2)+pow(dz, 2))
        return distance

    def getFaces(self):
        """ get the list of faces for the extrusion,
        and add faces for caps."""
        faces = []
        ls = self.shape.lenShape
        lv = len(self.vertices)

        # build faces with 2Dshape that has 2 normals per vertex
        if self.shape.vertDup == 1:
            # You build the faces between two points of the path3D
            # The number of faces  == number of non duplicated points of the
            # shape2D.
            if not self.arrow:
                faces = self.buildFaces(faces=faces,start=0,
                                        end=len(self.path3D)-1,
                                        ls=ls, dup=1)
            else:
                #ARROWS: 3*ls vertices have been added but only 2 faces.
                lp = len(self.path3D)
                # 1- BLOC1: cf drawing above.
                faces = self.buildFaces(faces = faces, start = 0,
                                        end = lp-(self.larrow+1),
                                        ls = ls, dup=1)

                #2- BLOC2 : Only two faces to close the back of the arrow.
                i = lp - self.larrow
                faces.append((i*ls, i*ls+2, i*ls+ls+2, i*ls+ls))
                faces.append((i*ls+ls/2, i*ls+ls/2+2,
                              i*ls+ls/2+ls+2, i*ls+ls+ls/2))

                #3- BLOC3: cf drawing above.
                faces = self.buildFaces(faces=faces,
                                        start=(lp-(self.larrow+1))+3,
                                        end=lp+3-1, ls=ls, dup=1)
                self.fixNormals(faces) # get correct normals for arrow


            # CAPS: get faces for the caps:
            # ls+1 vertices ( a section and a middle face point),
            # have been added to the regular vertices for the
            # caps and normals have been computed for those vertices.
            # The vertices are duplicated except for the middle face point,
            # but the normals are the same That is why below we are only using
            # every other vertices to get the faces.

            if (self.cap1 and not self.cap2) :
                # Only the beginning cap
                # the index of the first caps vertices is:
                # the last vertex  - (ls + 1)
                fc1 = lv - (ls + 1)
                for j in xrange(0,ls-2,2):
                    faces.append( (fc1+j+2, fc1+j, fc1+ls, fc1+ls))
                faces.append( (fc1, fc1+ls-2, fc1+ls, fc1+ls))

            elif (self.cap2 and not self.cap1):
                # Only one cap the end cap:
                # the index of the first caps vertices is:
                # the last vertex  - (ls + 1)
                fc2 = lv - (ls + 1)
                for j in xrange(0,ls-2,2):
                    faces.append((fc2+j, fc2+j+2, fc2+ls, fc2+ls))
                faces.append( (fc2+ls-2, fc2, fc2+ls, fc2+ls))
                
            elif self.cap1 and self.cap2:    # both caps
                # the index of the first cap vertex is :
                # last vertex (lv) - (2*(ls+1)) because 1st caps:
                fc1 = lv - (2*(ls+1))
                for j in xrange(0,ls-2,2):
                    faces.append((fc1+j+2, fc1+j, fc1+ls, fc1+ls))
                faces.append( (fc1, fc1+ls-2, fc1+ls, fc1+ls))

                #cap2:
                # the index of the second cap vertex is :
                # last vertex (lv) - ((ls+1))
                fc2 = lv - (ls+1)
                for j in xrange(0,ls-2,2):
                    faces.append((fc2+j, fc2+j+2,  fc2+ls, fc2+ls))
                faces.append( (fc2+ls-2,fc2, fc2+ls, fc2+ls))
                
        # not duplicate vertices; 1 normal per vertex
        else:
            if not self.arrow:
                faces = self.buildFaces(faces=faces, start=0,
                                        end=len(self.path3D)-1,
                                        ls=ls, dup=0)
            else:
                #ARROWS: 3*ls vertices have been added but only 2 faces.
                lp = len(self.path3D)
                # 1- BLOC1: cf drawing above.
                faces = self.buildFaces(faces=faces, start=0,
                                        end=lp-(self.larrow+1),
                                        ls=ls, dup=0)
                
                #2- BLOC2 : Only two faces to close the back of the arrow.
                #VERIFY!!!!
                i = lp - self.larrow
                faces.append((i*ls, i*ls+1, i*ls+ls+1, i*ls+ls))
                faces.append((i*ls+ls/2, i*ls+ls/2+1,
                              i*ls+ls/2+ls+1, i*ls+ls+ls/2))

                #3- BLOC3: cf drawing above.
                faces = self.buildFaces(faces=faces,
                                        start=(lp-(self.larrow+1))+3,
                                        end=lp+3-1, ls=ls, dup=0)
                #self.fixNormals(faces) # get correct normals for arrow

            # add caps
            if self.cap1 and not self.cap2:
                # first cap only
                fc1 = lv - (ls+1)
                for i in range(ls-1):
                    faces.append((fc1+i+1, fc1+i, fc1+ls, fc1+ls))
                faces.append( (fc1, fc1+ls-1, fc1+ls, fc1+ls))

            elif self.cap2 and not self.cap1:    # second cap only
                fc2 = lv - (ls+1)
                for i in range(ls-1):
                    faces.append((fc2+i+1, fc2+i, fc2+ls, fc2+ls))
                faces.append( (fc2, fc2+ls-1, fc2+ls, fc2+ls))

            elif self.cap1 and self.cap2:    # both caps
                # the index of the first cap vertex is :
                # last vertex (lv) - (2*(ls+1)) because 1st caps:
                fc1 = lv - (2*(ls+1))
                for j in xrange(ls-1):
                    faces.append((fc1+j+1, fc1+j, fc1+ls, fc1+ls))
                faces.append( (fc1, fc1+ls-1, fc1+ls, fc1+ls))

                #cap2:
                # the index of the second cap vertex is :
                # last vertex (lv) - ((ls+1))
                fc2 = lv - (ls+1)
                for j in xrange(ls-1):
                    faces.append((fc2+j, fc2+j+1,  fc2+ls, fc2+ls))
                faces.append( (fc2+ls-1,fc2, fc2+ls, fc2+ls))
        return faces


    def buildFaces(self, faces, start, end, ls, dup):
        """ Method that returns the faces for the points in the path3D
        between start and end. Ls is the length of the shape and dup
        indicates if the vertices are duplicated."""
        if dup:            
            jrange = range(1,ls-2,2)
            #                         b+6
            #       b+ls+6 | /b+ls+7   | /b+7      
            #              |/          |/            
            #              *-----------*
            #   b+ls+5    /|      b+5 /|
            #          | / |       | / |
            #          |/  | b+ls  |/  |  b+0  where b=i*ls
            #          *-----/-----*   | /     and i is the section number.
            #         /|   |/     /|   |/      When you build a facet he 
            # b+ls+4 / |   *-----/-|---*       normals must go out.
            #          |  /|    b+4|  /| 
            #          | / |       | / | 
            #          |/  b+ls+1  |/  b+1
            #          *-----------*       
            #         /|          /|
            # b+ls+3 / |         / |
            #          |       b+3  b+2
            #        b+ls+2 
            
        else:
            jrange = range(ls-1)
            # The vertices are not duplicated.
            #            
            #         
            #       b+ls+3 *-----------*b+3
            #             /|          /|
            #            / |         / |
            #           /  |        /  |        where b=i*ls
            #    b+ls+2*-----------*b+2|        and i is the section number.
            #          |   |       |   |
            #          |   *-------|---*b
            #          |  /b+ls    |  / 
            #          | /         | /  
            #          |/          |/  
            #          *-----------*       
            #         b+ls+1       b+1

            
        for i in xrange(start, end):
            # i represent the section number you are connecting.
            for j in jrange:
                # j represents the vertex number.
                # Build the (ls/2)-1 faces
                faces.append(((i*ls)+j, (i*ls)+(j+1),
                              ((i+1)*ls)+(j+1),((i+1)*ls)+j))
            # Build the last face.
            faces.append((i*ls, (i+1)*ls,
                         (i+1)*ls+(ls-1), (i*ls)+(ls-1)))
        return faces
        
      
    def fixNormals(self, faces):
        """ fixes normals on arrow by using glQuadNormals, which calculates
        normals."""
        ls = self.shape.lenShape
        path = self.path3D
        lp = len(path)
        la = self.larrow
        lnew = lp-la
        if self.shape.vertDup==1: numf = (ls-1)/2
        else: numf = ls-1
        arrfaces = faces[-(la+3)*numf:]
        # glQuadNormals need all the vertices of the secondarystructure but
        # only the faces belonging to the arrow.
        n = self.glQuadNormals(self.vertices, arrfaces)
        # Replacing the normals in the array by the new one.
        self.vnormals[lnew*ls:(lp+3)*ls] = n[lnew*ls:(lp+3)*ls]


    def glQuadNormals(self, vertices, faces):
        """ gets normals for quads by converting a quad face into two triangle
        faces and calling TriangleNormals. """
        faces = numpy.array(faces)
        f = numpy.concatenate((faces[:,-2:], faces[:,:1]), 1)
        F = numpy.concatenate((faces[:,:3], f), 0)
        n = TriangleNormals(vertices, F, 'PER_VERTEX')
        return n


    def setResidueProperties(self, properties, propName, resIndices):
        """ with the given list of pairs of residue indices and corresponding
        colors, the method sets the color per residue in the self.colors array.
        """

        nfaces = self.shape.lenShape
        if self.shape.vertDup==1: nfaces = nfaces/2
        properties = numpy.array(properties).astype('f')
        self.properties = getattr(self, propName)
        for r in range(len(resIndices)):
            st = int(resIndices[r][0]*nfaces)
            en = int(resIndices[r][1]*nfaces)
            for f in range(st, en):
                self.properties[f] = properties[r][:3]


    def setStripProperty(self, properties, propName):
        """ set the color array so that each side strip of the structure has
        one color. """

        nfaces = self.shape.lenShape
        if self.shape.vertDup==1: nfaces = nfaces/2

        assert len(properties)==nfaces
        c = 0
        if self.cap1: c = c+nfaces
        if self.cap2: c = c+nfaces
        for i in range(-c, 0):
            self.colors[i] = (1., 1., 1.)   # color caps white
            self.opacities[i] = 1.   # color caps white
        self.properties = getattr(self, propName)
        if self.arrow:
            for i in range(len(self.properties)-c):
                s = i%nfaces
                if i >= len(self.properties)-(self.larrow+2)*nfaces:
                    self.properties[i] = properties[s-2][:3]
                elif i == len(self.properties)-(self.larrow+2)*nfaces-1:
                    self.properties[i] = properties[s+1][:3]
                else:
                    self.properties[i] = properties[s][:3]
        else:
            for i in range(len(self.properties)-c):
                s = i%nfaces
                self.properties[i] = properties[s][:3]
        

from DejaVu2.Shapes import Circle2D

class ExtrudeCirclesWithRadii(ExtrudeObject):
    """
    Class to extrude a circle of varying radius
    """
    
    def __init__(self, path3D, matrix,
                 radii, quality=12, vertDup=0, firstDup=0,
                 cap1=0, cap2=0, arrow=0, larrow=3, warrow=2):
        """
        Constructor: Takes as arguments a path3D, matrix, 2D shape, and
        optional cap1, cap2, and arrow.  Calls getextrudeVertices() and
        getFaces() and stores the return values in self.vertices, self.vnormals
        , and self.faces.
        """
        assert len(path3D)==len(radii)
        self.radii = radii
        self.quality = quality
        self.vertDup = vertDup
        self.firstDup = firstDup
        shape = Circle2D(1.0, quality, vertDup, firstDup)
        ExtrudeObject.__init__(self, path3D, matrix, shape, cap1=cap1,
                               cap2=cap2, arrow=arrow, larrow=larrow,
                               warrow=warrow)


    def getextrudeVertices(self):
        """ get the extruded vertices and normals.  If an arrow is specified,
        extra vertices are added.  If caps are specified,
        the end vertices must be duplicated and used as extra vertices."""
        ls = self.shape.lenShape
        ctpts = numpy.array(self.shape.contpts)
        spts = ctpts.copy()
        length = len(self.path3D)

        f = g = 0
        # f adds another section, a section is ls vertices.
        # g adds the center pt of cap
        if self.cap1: f, g = f+1, g+1   # add a section and a middle point.
        if self.cap2: f, g = f+1, g+1
        if self.arrow:f = f+3           # add 3 sections:
                                        # - 1 junction between the arrow and
                                        # the rest of the geom.
                                        # - 2 for the caps of the arrow.

        # get the right size pts and normals arrays for caps and arrows
        pts = numpy.zeros( ((length+f)*ls+g,3),'f')
        ptsNorm = numpy.zeros( ((length+f)*ls+g,3),'f')

        # extrude 2D shape along 3D path
        quality = self.quality
        vertDup = self.vertDup
        firstDup = self.firstDup
        for i in range(length):
            radius = self.radii[i]
            #shape = Circle2D(radius, quality, vertDup, firstDup)
            spts[:,:2] = ctpts[:,:2]*radius # scale x, y for Shape2D
            pts[i*ls:i*ls+ls] = numpy.dot(spts, self.matrix[i]).astype('f')
            newNorm = numpy.dot(self.shape.contnorm,
                                  self.matrix[i]).astype('f')
            ptsNorm[i*ls:i*ls+ls] = newNorm

        # add 3*ls vertices for the arrow
        if self.arrow:
            pts = self.addArrow(pts)

        # add duplicate vertices and centers for cap1 and/or cap2 at the end
        # of the array of pts
        if self.cap1 and not self.cap2:
            pts[-ls-1:-1] = pts[0:ls]            
            pts[-1] = self.path3D[0]  # center of cap1
            # normals
            n = TriangleNormals(pts[-4:], [(3, 1, 0)])
            for k in range(ls+1):
                ptsNorm[-ls-1+k] = n[0]
                
        elif self.cap2 and not self.cap1:
            pts[-ls-1:-1] = pts[(length-1)*ls:length*ls]
            pts[-1] = self.path3D[-1]  # center of cap2
            # normals
            n = TriangleNormals(pts[-4:], [(0, 1, 3)])
            for k in range(ls+1):
                ptsNorm[-ls-1+k] = n[0]

        elif self.cap1 and self.cap2:
            pts[-2*ls-2:-ls-2] = pts[0:ls]   # vertices for cap1
            pts[-ls-2] = self.path3D[0] # center of cap1
            # normals
            n = TriangleNormals(pts[-5-ls:-ls-1], [(3, 1, 0)])
            for k in range(ls+1):
                ptsNorm[-2*ls-2+k] = n[0]

            pts[-ls-1:-1] =pts[(length-1)*ls:length*ls]  # vertices for cap2
            pts[-1] = self.path3D[-1]       # center of cap2
            # normals
            n = TriangleNormals(pts[-4:], [(0, 1, 3)])
            for k in range(ls+1):
                ptsNorm[-ls-1+k] = n[0]

        return pts, ptsNorm


class ExtrudeSSElt(ExtrudeObject):
    """ Class to take a shape and extrude it along the 3D path defined by
    a control coordinate and torsion coordinate 
    """

    def __init__(self, ssElt, shape, gapEnd = 0, gapBeg = 0, cap1=0,
                 cap2=0, arrow = 0,larrow = 3, warrow = 2):
        """Constructor's arguments:
        - secondary structure Elt.
        - shape.
        - etc..."""
        self.ssElt = ssElt
        self.gapBeg = gapBeg
        self.gapEnd = gapEnd
        self.lengthPath = len(self.ssElt.sheet2D.path)
        self.residuesInChain = self.ssElt.sheet2D.resInSheet
        self.residuesInChain.nbSS = range(len(self.residuesInChain))
        self.chords = self.ssElt.sheet2D.chords
        # last residue index in chain.
        self.lastResIndex = self.residuesInChain[-1].nbSS
        self.indexStart = self.ssElt.start.nbSS

        # first point and last point index in the path3D of the first residue
        #of the secondarystructure.
        self.fromForStart, self.toForStart = self.getResPts(self.indexStart)

##         self.indexEnd = self.residuesInChain.index(self.ssElt.end)
        self.indexEnd = self.ssElt.end.nbSS
        # first point and last point index in the path3D of the last residue
        #of the secondarystructure.
        self.fromForEnd, self.toForEnd = self.getResPts(self.indexEnd)

        path3D = self.ssElt.sheet2D.path[self.fromForStart:self.toForEnd]
        matrix = self.ssElt.sheet2D.matrixTransfo[self.fromForStart:
                                                  self.toForEnd]
        
        ExtrudeObject.__init__(self, path3D, matrix, shape, cap1=cap1,
                               cap2=cap2, arrow=arrow, larrow=larrow,
                               warrow=warrow)


    def getResIndexFromPts(self, respts):
        """ return the index of the residue to which a point in the path
        belongs."""
        
        if respts < (self.chords/2 + 1 ):
            # first residue
            resIndex = 0
        elif respts > ((self.lengthPath-1)-\
                        (self.chords + self.chords/2)):
            # last residue
            resIndex = self.lastResIndex
        else:
            # all the other nbchords
            resIndex = (respts-(2+ self.chords/2))/self.chords + 1
        

        return resIndex


    def getResPts(self, residueindex):
        """ return the index of the first and the last point in the
        Sheet2D.path for the residue whose index is specified"""
        # a residue is represented in the path3D by chords points.
        # first residue represented by nbchords/2 + 1
        # last residue represented by nbchords+nbchords/2
        # all other by nbchords.
        if residueindex == 0:
            fromPts = 0
            toPts = self.chords/2 + 2

        elif residueindex == self.lastResIndex:
            fromPts = (residueindex-1) * self.chords + self.chords/2+1
            toPts = self.lengthPath-1

        else:
            fromPts = (residueindex-1) * self.chords + self.chords/2+1
            toPts = fromPts + self.chords +1

        toPts = toPts-self.gapEnd
        fromPts = fromPts + self.gapBeg
        return fromPts,toPts


    def getExtrudeResidues(self, resSet, gapBefore=False, gapAfter=False):
        """ Get faces for the specified residues in ResSet """

        if gapBefore:
            resSet = resSet[1:]
            #print ss, 'REMOVING 1st RES'
        if gapAfter:
            resSet = resSet[:-1]
            #print ss, 'REMOVING last RES'
        ls = self.shape.lenShape
        nfaces = ls
        if self.shape.vertDup==1: nfaces = (ls)/2
        residueFaces = []
        resSet.sort()
        
        residueFacesDict = {}
        for r in resSet:
            residueFacesDict[r] = []
            rindex = r.nbSS
            start, end = self.getResPts(rindex)
            # 1st and last points in the path3D corresponding to r

            first = (start-self.fromForStart)*nfaces
            last =  (end-self.fromForStart-1)*nfaces
            # relative to the secondarystructure.

            # faces for arrow (Note: assumes arrow shorter than last residue)
            if r==self.ssElt.end and self.arrow:
                #Only two faces are added when an arrow is added.
                last = last+2
            residueFaces.extend(self.faces[first:last])
            residueFacesDict[r].extend(self.faces[first:last])

            # faces for cap1 or cap2, if residue is beginning or end
            if r==self.ssElt.start:
                if self.cap1 and not self.cap2:
                    residueFaces.extend(self.faces[-nfaces:])
                    residueFacesDict[r].extend(self.faces[-nfaces:])
                    
                elif self.cap1 and self.cap2:
                    residueFaces.extend(self.faces[-2*nfaces:-nfaces])
                    residueFacesDict[r].extend(self.faces[-2*nfaces:-nfaces])
            if r==self.ssElt.end:
                if self.cap2:
                    residueFaces.extend(self.faces[-nfaces:])
                    residueFacesDict[r].extend(self.faces[-nfaces:])
                    
        return residueFaces, residueFacesDict


    def getExtrudeProperties(self, ResSet, propName):
        """ Get the colors for the specified residues in ResSet """

        ls = self.shape.lenShape
        nfaces = ls
        if self.shape.vertDup==1: nfaces = (ls)/2
        # res are the residues used to build the sheet 2D. 
        prop = []
        resProperty = getattr(self, propName).tolist()

        # sort residues in secondarystructure according to nbSS
        ResSet.sort()
        for r in ResSet:
            rindex = r.nbSS
            start, end = self.getResPts(rindex)
            # 1st and last points in the whole path3D corresponding to r

            first = (start-self.fromForStart)*nfaces
            last = (end-self.fromForStart-1)*nfaces
            # 1st and last points index relative to the secondarystructure.

            # colors for arrow (Note: assumes arrow shorter than last residue)
            if r==self.ssElt.end and self.arrow:
                # 3 corresponds to the sections added for the arrow.
                last = last+2
            prop.extend(resProperty[first:last])

            # colors for cap1 or cap2, if residue is beginning or end
            if r==self.ssElt.start:
                if self.cap1 and not self.cap2:
                    prop.extend(resProperty[-nfaces:])
                elif self.cap1 and self.cap2:
                    prop.extend(resProperty[-2*nfaces:-nfaces])
            if r==self.ssElt.end:
                if self.cap2:
                    prop.extend(resProperty[-nfaces:])
                    
        return prop


    def setResProperties(self, propVect, propName, ResSet):
        """ sets each residue in 'ResSet' to the color in 'colors' by calling
        the setResidueColors method of extrusion, which is an instance of the
        Extrude class. """
        
        assert len(propVect)==len(ResSet)
        ls = self.shape.lenShape
        nfaces = ls
        if self.shape.vertDup==1: nfaces = (ls)/2
        firstLastIndices = []
        properties = []
        i=0
        for r in ResSet:
            assert r in self.ssElt.residues
            rindex = r.nbSS
            start, end = self.getResPts(rindex)
            first = start-self.fromForStart
            last  =  end-self.fromForStart-1
            # relative to the secondarystructure.

            # color arrow same color as residue
            if r==self.ssElt.end and self.arrow:
                # When an arrow is added the nbSS of faces is increased by
                # 2. If we add 1 to End this will mean that 4 faces have been
                # added .
                last = last + 0.5  

            firstLastIndices.append((first, last))
            propTmp = propVect[i]
            properties.append(propTmp)
            # color the caps same color as residue
            if r==self.ssElt.start:
                if self.cap1 and not self.cap2:
                    firstLastIndices.append((-1, 0))
                    properties.append(propTmp)
                elif self.cap1 and self.cap2:
                    firstLastIndices.append((-2, -1))
                    properties.append(propTmp)
            if r==self.ssElt.end:
                if self.cap2:
                    firstLastIndices.append((-1, 0))
                    properties.append(propTmp)
            i = i+1

        ExtrudeObject.setResidueProperties(self,properties,
                                           propName,firstLastIndices )


    def getResIndexFromExtrudeVertex(self, vertex):
        """ takes a vertex from the extrusion and returns the index of the
        residue the vertex belongs to """
        assert vertex < len(self.vertices)
        ls = self.shape.lenShape
        # find the section number of the path3D to which the picked
        # vertex belongs to
        pickedEstimate = vertex/ls
        # number of section that have been added to the path3D.
        extraSection = 0
        # len(self.paths3D)-1 because there is a section shared by
        # two residues.
        if pickedEstimate>=len(self.path3D)-1:
            if self.arrow:
                extraSection = 3
                if pickedEstimate <= (len(self.path3D)-1)+extraSection:
                    # What comes after normal points are the arrow sections
                    # then the cap1 there are no cap2.
                    # vertex is in arrow, p is the point index if no
                    # Arrow added
                    pickedSection = len(self.path3D) - 2

                else:
                    pointsLeft = vertex-(len(self.path3D)+extraSection)*ls
                    # remove from the vertex index the whole path3D and the
                    # extraSection added for the arrow to see if the vertex
                    # is in caps. The cap extraSection have been added at the
                    # end.
                    pickedEstimate = pointsLeft/(ls+1)
                    pickedSection = pickedEstimate
            else:
                pointsLeft = vertex-(len(self.path3D)+extraSection)*ls
                pickedEstimate = pointsLeft/(ls+1)
                pickedSection = pickedEstimate

            if pickedEstimate==1:
                pickedSection = len(self.path3D)-2
                # vertex is in cap2
                
            if pickedEstimate<=0:
                if self.cap1:
                    pickedSection = 1                    # vertex is in cap1
                elif self.cap2:
                    pickedSection = len(self.path3D)-2   # vertex is in cap2
                else:
                    pickedSection = pickedEstimate
                    
        else: pickedSection = pickedEstimate
        index = self.getResIndexFromPts(pickedSection+self.fromForStart)
        
        resIndex = index - self.indexStart

        # remove first section because only represented by self.chords/2 + 1
        p = ((self.lengthPath-1)-(self.chords + self.chords/2))
        if pickedSection < self.chords/2+2:
            resIndex = 0

        elif pickedSection >((self.lengthPath-1)-\
                        (self.chords + self.chords/2)):
            resIndex = self.lastResIndex

        if resIndex==-1: resIndex=0
        return resIndex
        
    
from mglutil.math import crossProduct,norm

def getAtsRes(atoms,listesel):
    missingAts = False
    listeCoord = []
    for astr in listesel :
        AT = atoms.objectsFromString(astr)  
        if len(AT) :
            ATcoord = numpy.array(AT[0].coords)
            listeCoord.append(ATcoord)
        else : missingAts = True
    return listeCoord,missingAts

def ExtrudeNA(chain):
    """Computes ribbons for DNA/RNA"""
    coord = []
    coord.append(chain.DNARes[0].atoms[0].coords)
    NA_type = chain.DNARes[0].type.strip()                        
    atoms = chain.DNARes[0].atoms
    missingAts = False
    normal = numpy.array([0.,1.,0.])
    if NA_type in ['A', 'G']:
        listesel = ['N9.*','C8.*','C4.*']
        listeCoord,missingAts = getAtsRes(atoms,listesel)
        if not missingAts :
            N9 =  listeCoord[0]#numpy.array(atoms.objectsFromString('N9.*')[0].coords)
            C8 =  listeCoord[1]#numpy.array(atoms.objectsFromString('C8.*')[0].coords)
            C4 =  listeCoord[2]#numpy.array(atoms.objectsFromString('C4.*')[0].coords)
            N9_C8 = C8-N9
            N9_C4 = C4-N9
            normal = numpy.array(crossProduct(N9_C8, N9_C4, normal=True))
    else:
        listesel = ['N1.*','C2.*','C6.*']
        listeCoord = []
        listeCoord,missingAts = getAtsRes(atoms,listesel)
        if not missingAts :
            N1 =  listeCoord[0]#numpy.array(atoms.objectsFromString('N1.*')[0].coords)
            C2 =  listeCoord[1]#numpy.array(atoms.objectsFromString('C2.*')[0].coords)
            C6 =  listeCoord[2]#numpy.array(atoms.objectsFromString('C6.*')[0].coords)
            N1_C2 = C2-N1
            N1_C6 = C6-N1
            normal = numpy.array(crossProduct(N1_C2, N1_C6, normal=True))
    base_normal = numpy.array(chain.DNARes[0].atoms[0].coords)
    coord.append((base_normal + normal).tolist())

    for res in chain.DNARes[1:]:
        normal = numpy.array([0.,1.,0.])
        if res.atoms.objectsFromString('P.*'):
            P_coord = res.atoms.objectsFromString('P.*')[0].coords
            coord.append(P_coord)
        else: # this in case last residue does not have P
            P_coord = res.atoms[0].coords
        NA_type = res.type.strip()      
        atoms = res.atoms
        if NA_type in ['A', 'G']:
            listesel = ['N9.*','C8.*','C4.*']
            listeCoord,missingAts = getAtsRes(atoms,listesel)
            if not missingAts :
                N9 =  listeCoord[0]#numpy.array(atoms.objectsFromString('N9.*')[0].coords)
                C8 =  listeCoord[1]#numpy.array(atoms.objectsFromString('C8.*')[0].coords)
                C4 =  listeCoord[2]#numpy.array(atoms.objectsFromString('C4.*')[0].coords)
                N9_C8 = C8-N9
                N9_C4 = C4-N9
                normal = numpy.array(crossProduct(N9_C8, N9_C4, normal=True))
        else:
            listesel = ['N1.*','C2.*','C6.*']
            listeCoord = []
            listeCoord,missingAts = getAtsRes(atoms,listesel)
            if not missingAts :
                N1 =  listeCoord[0]#numpy.array(atoms.objectsFromString('N1.*')[0].coords)
                C2 =  listeCoord[1]#numpy.array(atoms.objectsFromString('C2.*')[0].coords)
                C6 =  listeCoord[2]#numpy.array(atoms.objectsFromString('C6.*')[0].coords)
                N1_C2 = C2-N1
                N1_C6 = C6-N1
                normal = numpy.array(crossProduct(N1_C2, N1_C6, normal=True))
        base_normal = numpy.array(P_coord)
        coord.append((base_normal + normal).tolist())
        
    chain.sheet2D['ssSheet2D'] = Sheet2D()
    chain.sheet2D['ssSheet2D'].compute(coord, len(chain.DNARes)*(False,), 
                             width = 2.0,off_c = 0.9,offset=0.0, nbchords=4)
    chain.sheet2D['ssSheet2D'].resInSheet = chain.DNARes

