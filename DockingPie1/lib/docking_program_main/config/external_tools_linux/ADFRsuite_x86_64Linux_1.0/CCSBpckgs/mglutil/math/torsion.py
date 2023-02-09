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

import numpy

#taken from Pmv/measureCommands.py

def torsion( x1, x2, x3, x4):
    """
    Compute the torsion angle between x1, x2, x3, x4.
    All coordinates are cartesian; result is in degrees.
    Raises a ValueError if angle is not defined.
    """
    from math import sqrt, acos
    import numpy
    
    tang=0.0
    x1 = numpy.array(x1, 'f')
    x2 = numpy.array(x2, 'f')
    x3 = numpy.array(x3, 'f')
    x4 = numpy.array(x4, 'f')
    
    assert x1.shape == (3, )
    assert x2.shape == (3, )
    assert x3.shape == (3, )
    assert x4.shape == (3, )

    a = x1-x2
    b = x3-x2
    c = vvmult(a, b)

    a = x2-x3
    b = x4-x3
    d = vvmult(a, b)

    dd=sqrt(numpy.sum(c*c))
    de=sqrt(numpy.sum(d*d))

    if dd<0.001 or de<0.001:
        raise ValueError ( 'Torsion angle undefined, degenerate points')

    vv = numpy.dot(c, d) / (dd*de);
    if vv<1.0: tang=vv
    else: tang= 1.0
    if tang<-1.0: tang=-1.0
    tang = acos(tang)
    tang = tang*57.296

    b = vvmult(c, d)
    if numpy.dot(a, b) > 0.0: tang = -tang
    return tang

def vvmult( a, b):
    """
    Compute a vector product for 3D vectors
    """
    res = numpy.zeros(3, 'f')
    res[0] = a[1]*b[2] - a[2]*b[1]
    res[1] = a[2]*b[0] - a[0]*b[2]
    res[2] = a[0]*b[1] - a[1]*b[0]
    return res
