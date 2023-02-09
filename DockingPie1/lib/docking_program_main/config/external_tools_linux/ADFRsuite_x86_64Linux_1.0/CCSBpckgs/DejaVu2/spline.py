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
# $Header: /mnt/raid/services/cvs/DejaVu2/spline.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#
# $Id: spline.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#

"""
reimplementation of lincrv.c in python MS April 2000
original code from Graphics Gem V pp220-222  Ken Shoemake, 1994
"""
import numpy
from copy import deepcopy

# Perform a generic vector unary operation.

def lerp(t, a0, a1, p0, p1):
    t0 = (a1-t)/(a1-a0)
    t1=1-t0
    p = t0*p0 + t1*p1
    return p.astype('f')

"""
    DialASpline(t,a,p,m,n,work,Cn,interp,val) computes a point val at
    parameter t on a spline with knot values a and control points p. The
    curve will have Cn continuity, and if interp is TRUE it will interpolate
    the control points.
    Possibilities include Langrange interpolants, Bezier curves, Catmull-Rom
    interpolating splines, and B-spline curves. Points have m coordinates, and
    n+1 of them are provided. The work array must have room for n+1 points.
"""
def DialASpline( t, a, p, m=3, n=4, work=None, Cn=2, interp=1):

    if Cn>n-1: Cn = n-1         # Anything greater gives one polynomial
    k=0
    while ( t>a[k]): k=k+1      # Find enclosing knot interval
    h=k
    while (t==a[k]): k = k+1    # May want to use fewer legs
    if k>n:
        k = n
        if h>k:
            h = k
    h = 1+Cn - (k-h)
    k=k-1
    lo = k-Cn
    hi = k+1+Cn

    if interp:                # Lagrange interpolation steps
        drop=0
        if lo<0:
            lo = 0
            drop = drop + Cn-k
            if hi-lo<Cn:
                drop = drop + Cn-hi
                hi = Cn
        
        if hi>n:
            hi = n
            drop = drop + k+1+Cn-n
            if hi-lo<Cn:
                drop = drop + lo-(n-Cn)
                lo = n-Cn
        
        for i in range(lo, hi+1):
            work[i] = p[i][::-1]

        for j in range(1, Cn+1):
            for i in range(lo, hi-j+1):
##                  work[i] = lerp(t,a[i],a[i+j],work[i],work[i+1])
##                  print work[i]
                t0 = (a[i+j]-t)/(a[i+j]-a[i])
                t1 = 1-t0
                work[i] = (t0*work[i]+t1*work[i+1]).astype('f')
        h = 1+Cn-drop
    else:                    # Prepare for B-spline steps
        if lo<0:
            h = h+lo
            lo = 0
        
        for i in range(lo, lo+h+1):
            work[i] = p[i][::-1]

        if h<0: h = 0

    for j in range(0, h):
        tmp = 1+Cn-j
        for i in range( h-1, j-1, -1 ):
            #work[lo+i+1] = lerp(t,a[lo+i],a[lo+i+tmp],work[lo+i],work[lo+i+1])
            t0 = (a[lo+i+tmp]-t)/(a[lo+i+tmp]-a[lo+i])
            t1 = 1-t0
            work[lo+i+1] = (t0*work[lo+i]+t1*work[lo+i+1]).astype('f')
    val = deepcopy(work[lo+h][::-1])

    return tuple(val)


class SplineObject:
    """ This class implements a set of methods to compute a spline using
    dialASpline function"""

    def __init__(self, coords, name="", nbchords=4, interp='interpolation',
                 continuity=2, closed=0):
        """
        The constructor of a SplineObjects takes the following arguments:
        coords    : set of coordinates which are the control points of
                    the spline
        name      : Name of the spline (default: '')
        nbchords  : Number of point per control points (nbchords=4) 
        interp    : flag set to 'interpolation' will interpolate the control
                    points is set to 'approximation' it will approximate
                    them.)
        continuity: Define the continuity of the curve (default=2)
        closed    : Flag if set to 1 the computed spline will be closed.
        """
        
        self.name = name
        self.nbchords = nbchords
        self.interp = interp
        self.continuity = continuity
        self.closed = closed
        # If closed==1 some control points needs to be added.
        if closed == 1:
            # the last point:
            llPoint = coords[-2]
            lPoint = coords[-1]
            # 1st Point
            middlePoint = coords[0]
            # 2nd Point
            fPoint = coords[1]
            ffPoint = coords[2]

            oldCoords = coords
            # Now create the new array of control points coords:
            coords.insert(0, lPoint)
            coords.insert(0, llPoint)
            coords.append(middlePoint)


            coords.append(fPoint)
            coords.append(ffPoint)

        self.ctlPts = numpy.array(coords,'f')
        smooth = self.computePath3D()
        if closed:
            index = 2*nbchords-1
            self.smooth = smooth[index:-index]
            
        else:
            self.smooth = self.addMorePoints(smooth)
            
    def computePath3D(self):
        # ctlPts, nbchords, interp, Cn):
        # 1- Compute the path3D.
        # (using the DialASpline function from DejaVu2.spline)
        work = numpy.zeros( (len(self.ctlPts),3), 'f' )
        Kts, interp = self.computeKnots() 
        #print 'Flavor PLY: interp=%d, Cn=%d'%(interp, Cn)
        nbPts = (len(self.ctlPts)-1)*self.nbchords
        tinc = 1.0/nbPts
        smooth = []
        for t in range(nbPts+1):
            #print t
            val = DialASpline(t*tinc, Kts, self.ctlPts, m=3,
                              n=len(self.ctlPts)-1,
                              work=work, Cn=self.continuity, interp=interp)
            smooth.append(val)
        #print 'end'
        return smooth

    def addMorePoints(self, smooth):
        # Add two points one at the beginning and one at the end so the
        # the extruded geometries includes the first and last atom.

        # Add the first point
        pt1 = numpy.array(smooth[0])
        pt2 = numpy.array(smooth[1])
        vec1 = pt1 - pt2
        firstPoint = tuple(pt1 +  vec1)
        smooth.insert(0, firstPoint)

        # Add the last point:
        pt1 = numpy.array(smooth[-1])
        pt2 = numpy.array(smooth[-2])
        vec2 = pt2 - pt1
        lastPoint = tuple(pt1-vec2)
        smooth.append(lastPoint)

        return smooth

    def computeKnots(self):
        #, interp, length):
        length = len(self.ctlPts)
        if self.interp == 'interpolation':
            #print 'in interpolation'
            Kts = []
            interp = 1
            inc = 1.0/(length-1)
            for i in xrange(length):
                Kts.append( i*inc )
            Kts.append(999999999999999.99)

        elif self.interp =='approximation':
            #print 'in  approximation'
            Kts = [0,1]
            interp = 0
            Kts = Kts*((length)/2)
            Kts.sort()
            Kts.append(999999999999999.99)

        #print 'Kts',Kts
        return Kts,interp
    
    



if __name__ == '__main__':
    import numpy
    ctlPts = numpy.array ( (
        ( 16.967 , 12.784 , 4.338),
        ( 13.856 , 11.469 , 6.066),
        ( 13.660 , 10.707 , 9.787),
        ( 10.646 , 8.991 , 11.408),
        ( 9.448 , 9.034 , 15.012),
        ( 8.673 , 5.314 , 15.279),
        ( 8.912 , 2.083 , 13.258),
        ( 5.145 , 2.209 , 12.453),
        ( 5.598 , 5.767 , 11.082),
        ( 8.496 , 4.609 , 8.837),
        ( 6.500 , 1.584 , 7.565),
        ( 3.545 , 3.935 , 6.751),
        ( 5.929 , 6.358 , 5.055),
        ( 7.331 , 3.607 , 2.791),
        ( 3.782 , 2.599 , 1.742),
        ( 2.890 , 6.285 , 1.126),
        ( 5.895 , 6.489 , -1.213),
        ( 4.933 , 3.431 , -3.326),
        ( 2.792 , 5.376 , -5.797),
        ( 5.366 , 8.191 , -6.018),
        ( 3.767 , 10.609 , -3.513),
        ( 6.143 , 13.513 , -2.696),
        ( 8.114 , 13.103 , .500),
        ( 6.614 , 16.317 , 1.913),
        ( 3.074 , 14.894 , 1.756),
        ( 4.180 , 11.549 , 3.187),
        ( 5.879 , 13.502 , 6.026),
        ( 2.691 , 15.221 , 7.194),
        ( .715 , 12.045 , 6.657),
        ( 2.986 , 9.994 , 8.950),
        ( 4.769 , 12.336 , 11.360),
        ( 8.140 , 11.694 , 9.635),
        ( 10.280 , 14.760 , 8.823),
        ( 12.552 , 15.877 , 6.036),
        ( 15.930 , 17.454 , 6.941),
        ( 18.635 , 18.861 , 4.738),
        ( 21.452 , 16.969 , 6.513),
        ( 22.019 , 13.242 , 7.020),
        ( 21.936 , 12.911 , 10.809),
        ( 18.504 , 12.312 , 12.298),
        ( 17.924 , 13.421 , 15.877),
        ( 17.334 , 10.956 , 18.691),
        ( 13.564 , 11.573 , 18.836),
        ( 13.257 , 10.745 , 15.081),
        ( 15.445 , 7.667 , 15.246),
        ( 13.512 , 5.395 , 12.878)
        ), 'f')

    import  pdb
    work = numpy.zeros( (len(ctlPts),3), 'f' )

    lagKts = []
    inc = 1.0/(len(ctlPts)-1)
    for i in xrange(len(ctlPts)):
        lagKts.append( i*inc )
    lagKts.append(999999999999999.99)

    interp = 1
    Cn = 2
    tinc = 1.0/460;
    vals = []
    #print 'Start....'
    print 'Flavor PLY: interp=1, Cn=2'
    print len(ctlPts)
    for t in range(460+1):
        val = DialASpline(t*tinc, lagKts, ctlPts, m=3, n=len(ctlPts)-1,
                          work=work, Cn=Cn, interp=interp)

    ##      val = DialASpline(t*tinc, lagKts, ctlPts, m=3, n=45,
    ##                        work=work, Cn=2, interp=1)
        vals.append( val )
        #print '%9.6f %9.6f %9.6f'%(val[0],val[1],val[2])
    #print 'Done'

    ##  from DejaVu2.Polylines import Polylines
    ##  p = Polylines('pyspline', vertices = (vals,))
    ##  self.GUI.VIEWER.AddObject(p)

    #c test
    ##  f = open("/mgl/home/saxisa/mvprojet/pmv/spline/a")
    ##  data = f.readlines()
    ##  f.close()
    ##  import string
    ##  data = map(lambda x: string.split(x), data)
    ##  data = data[1:]
    ##  data = map(lambda x: [float(x[0]),float(x[1]),float(x[2])], data)
    ##  l = Polylines('cspline', vertices = (data,))
    ##  self.GUI.VIEWER.AddObject(l)
