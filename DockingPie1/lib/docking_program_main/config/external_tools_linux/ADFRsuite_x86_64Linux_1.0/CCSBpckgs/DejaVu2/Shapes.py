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
# Copyright: M. Sanner TSRI 2000
#
#############################################################################
#
# $Header: /mnt/raid/services/cvs/DejaVu2/Shapes.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#
# $Id: Shapes.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#


import numpy, math
from math import sqrt

from opengltk.OpenGL import GL


"""
This module implements a set of Classes to describe 2D shape that can be used
to  perform an extrusion along a 3D path.
"""
   

class Shape2D:
    """
    Shape2D is a base class describing a 2D shape.
    """
    def __init__(self, contpts, contnorm, vertDup=0, firstDup=0):
        """
        contpts   : list of 3D coordinates of the control points describing the
                    shape2D.
        contnorm  : list of 3D vector specifying the normals at each vertex
                    of the shape2D.
        vertDup   : Flag specifying whether or not to duplicate each vertex and
                    normal.
        firstDup  : Flag specifying whether or not to duplicate the first
                    first vertex.
        """
        self.vertDup = vertDup
        self.firstDup = firstDup
        self.contpts = list(contpts)
        self.contnorm = list(contnorm)
        if vertDup:
            firstDupInd = 2
        else:
            firstDupInd = 1

        if firstDup:
            cont = self.contpts
            newcont = cont + cont[:firstDupInd]
            self.contpts = numpy.array(newcont)

            contn = self.contnorm
            newcontn = contn + contn[:firstDupInd]
            self.contnorm = numpy.array(newcontn)
            
        self.lenShape = len(self.contpts)
        
    def capGeom(self, vertices, vertInds, lastInd, flip=False):
        """
        generate vertices, normals and face indice (as quads) for the shape
        when the shape is extruded. This only works for convex shapes.
        If ther are more than 4 vertices in the shape a vertex is added at the
        center of the shape and the quads with repea5ting 4th vertex are used
        the create triangles for the cap.
        """
        l = len(vertices)
        end = lastInd
        if hasattr(vertices, 'tolist'): vertices = vertices.tolist()
        x1,y1,z1 = vertices[0]
        x2,y2,z2 = vertices[1]
        x3,y3,z3 = vertices[2]
        A = [x2-x1, y2-y1, z2-z1]
        B = [x3-x1, y3-y1, z3-z1]
        N = [A[1]*B[2]-B[1]*A[2],
             A[2]*B[0]-B[2]*A[0],
             A[0]*B[1]-B[0]*A[1]]
        if flip:
            n = -1.0/sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2])
        else:
            n = 1.0/sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2])
        normal = [N[0]*n, N[1]*n, N[2]*n]

        if l == 3:
            if flip:
                faces = [[end, end+1, end+2, end]]
            else:
                faces = [[end, end+2, end+1, end]]
            return vertices[:], normal, faces
        elif l==4:
            if flip:
                faces = [[end, end+1, end+2, end+3]]
            else:
                faces = [[end+3, end+2, end+1, end]]
            return vertices[:], normal, faces

        else:
            center = numpy.sum(vertices, 0)/l
            if hasattr(vertices, 'tolist'): vertices = vertices.tolist()
            faces = []
            for i in range(l-1):
                if flip:
                    faces.append( [end, end+1+i, end+2+i, end] )
                else:
                    faces.append( [end, end+2+i, end+1+i, end] )
            if flip:
                faces.append( [end, end+2+i, end+1, end] )
            else:
                faces.append( [end, end+1, end+2+i, end] )
            
            return [center]+vertices, normal, faces

    def returnStringRepr(self):
        """ This method returns the string to be evaluated to create
        the object"""
        cname = self.__class__.__name__
        before = 'from %s import %s'%(self.__module__, cname)
        st = "%s(%s,%s, vertDup=%s, firstDup=%s)"%(cname,
                                                   self.contpts,
                                                   self.contnorm,
                                                   self.vertDup,
                                                   self.firstDup)
        
        return before, st
    
class Triangle2D(Shape2D):
    """ Class derived from Shape2D describing a Triangle."""
    def __init__(self, side = 1.0, vertDup=0, firstDup = 0 ):
        self.side = side
        x = side/2
        y1 = -side/(2*math.sqrt(3))
        y2 = side/math.sqrt(3)
        if not vertDup:
            pts = ( (-x, y1, 1), (0, y2, 1), (x, y1, 1) )
            norms = ((-math.sqrt(3)/2, -.5, 0), (0, 1, 0),
                      (math.sqrt(3)/2, -.5, 0),  )
        else:
            pts   = ( (-x, y1 , 1)    , (-x, y1, 1),
                      ( 0, y2 , 1)    , ( 0, y2, 1),
                      ( x, y1 , 1)    , ( x, y1, 1) )            
            norms = ( (0, -2*x, 0),    ( y1-y2,  x , 0),
                      ( y1-y2,  x , 0),( y2-y1,  x , 0),
                      ( y2-y1,  x , 0),(0, -2*x, 0),)
        
        Shape2D.__init__(self, pts, norms, vertDup=vertDup, firstDup=firstDup)

    def returnStringRepr(self):
        """ This method returns the string to be evaluated to create
        the object"""
        cname = self.__class__.__name__
        before = 'from %s import %s'%(self.__module__, cname)
        st = "%s(side=%s, vertDup=%s, firstDup=%s)"%(cname,
                                                     self.side,
                                                     self.vertDup,
                                                     self.firstDup)
        
        return before, st

class Ellipse2D(Shape2D):
    """ Class derived from Shape2D describing a Ellipse """

        
    def __init__(self, demiGrandAxis ,
                 demiSmallAxis, quality=6, vertDup=0, firstDup=0):
        """demiGrandAxis is 1/2 the width of the ellipse
        demiSmallAxis is 1/2 the height of the ellipse"""
        self.quality = quality
        self.demiGrandAxis = demiGrandAxis
        self.demiSmallAxis = demiSmallAxis
        
        circle = numpy.zeros( (quality,3) ).astype('f')
        
        circleNormal = numpy.zeros( (quality,3) ).astype('f')

        # when points are duplicated:
        # norm0 = (x(i-1)-xi,y(i-1)-y, z(i-1)-z(i)) cross (0,0,1) 
        # norm1 = (0, 0, 1) cross (x(i+1)-xi,y(i+1)-y, z(i+1)-z(i))

        # x = y1*z2 - y2*z2
        # y = z1*x2 - z2*x1
        # z = x1*y2 - x2*y1

        for i in range( quality ):
            circle[i][0] = 2*demiGrandAxis*math.cos( i*2*math.pi/quality)
            circle[i][1] = -2*demiSmallAxis*math.sin( i*2*math.pi/quality)
            circle[i][2] = 1

            circleNormal[i][0] = math.cos( i*2*math.pi/quality )
            circleNormal[i][1] = -math.sin( i*2*math.pi/quality )
            circleNormal[i][2] = 0

        if vertDup:
            pts  = numpy.zeros( (quality*2, 3)) .astype('f')
            norm = numpy.zeros( (quality*2, 3)) .astype('f')
            # index for pts and norm
            ptsInd = 0
            for ind in range(quality):
                if ind == 0:
                    prev = quality-1
                    next = ind+1
                elif ind == quality-1:
                    next = 0
                    prev = ind-1
                else:
                    next = ind + 1
                    prev = ind - 1

                # Compute the Vprev vector and the Vnext vector
                Vprev = circle[prev]-circle[ind]
                Vnext = circle[next]-circle[ind]

                n0 = [ Vprev[1], -Vprev[0], 0 ]
                n1 = [-Vnext[1],  Vnext[0], 0 ]
                norm[ptsInd], norm[ptsInd+1] = n0, n1
                pts[ptsInd], pts[ptsInd+1] = circle[ind], circle[ind]

                ptsInd = ptsInd + 2
        
            circle = pts
            circleNormal = norm

        Shape2D.__init__(self, circle, circleNormal, vertDup=vertDup,
                         firstDup=firstDup)

    def returnStringRepr(self):
        """ This method returns the string to be evaluated to create
        the object"""
        cname = self.__class__.__name__
        before = 'from %s import %s'%(self.__module__, cname)
        st = "%s(%s, %s, quality=%s,  vertDup=%s, firstDup=%s)"%(cname,
                                                                self.demiGrandAxis,
                                                                self.demiSmallAxis,
                                                                self.quality,
                                                                self.vertDup,
                                                                self.firstDup)
        
        return before, st


class Circle2D(Ellipse2D):
    """ Class derived from Ellipse2D describing a Circle."""
    def __init__(self, radius, quality=12, vertDup=0, firstDup=0):
        self.radius = radius
        Ellipse2D.__init__(self, radius, radius, quality=quality,
                           firstDup=firstDup, vertDup=vertDup)

    def returnStringRepr(self):
        """ This method returns the string to be evaluated to create
        the object"""
        cname = self.__class__.__name__
        before = 'from %s import %s'%(self.__module__, cname)
        st = "%s(%s, quality=%s, vertDup=%s, firstDup=%s)"%(cname,
                                                            self.radius,
                                                            self.quality,
                                                            self.vertDup,
                                                            self.firstDup)
        return before, st


class Rectangle2D(Shape2D):
    """ Class derived from Shape2D describing a Rectangle """
    def __init__(self, width, height, vertDup=0, firstDup=0):
        self.width = width
        self.height = height
        if not vertDup:
            pts = ( (-width, -height, 1), (-width, height, 1),
                    (width, height, 1), (width, -height, 1) )
            norms = ( (-1,-1,0), (-1,1,0), (1,1,0), (1,-1,0) )
        else:
            pts = ( (-width, -height, 1), (-width, -height, 1),
                    (-width, height, 1), (-width, height, 1),
                    (width, height, 1), (width, height, 1),
                    (width, -height, 1), (width, -height, 1) )
            norms = ( (0,-1,0), (-1,0,0), (-1,0,0), (0,1,0),
                      (0,1,0), (1,0,0), (1,0,0), (0,-1,0) )
            

        Shape2D.__init__(self, pts, norms, vertDup=vertDup, firstDup=firstDup)

    def returnStringRepr(self):
        """ This method returns the string to be evaluated to create
        the object"""
        before = 'from %s import %s'%(self.__module__, self.__class__.__name__)
        st = "%s(%s, %s, vertDup=%s, firstDup=%s)"%(self.__class__.__name__,
                                                    self.width, self.height,
                                                    self.vertDup,self.firstDup)
        return before, st


class Square2D(Rectangle2D):
    """ Class derived from Shape2D describing a Square """
    def __init__(self, side, vertDup=0, firstDup=0):
        self.side = side
        Rectangle2D.__init__(self, side, side, vertDup=vertDup,
                             firstDup=firstDup)
    def returnStringRepr(self):
        """ This method returns the string to be evaluated to create
        the object"""
        before = 'from %s import %s'%(self.__module__, self.__class__.__name__)
        st = "%s(%s, vertDup=%s, firstDup=%s)"%(self.__class__.__name__,
                                                       self.side,
                                                       self.vertDup, self.firstDup)
        return before, st

