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

#
#  collection of standard 4x4 rotation matrices
#  about the X, Y and Z axis of 10, 30, 45, 90, 180 degrees
#
rotations = {
}

import math, numpy
from mglutil.math.rotax import rotax

orig = numpy.array( (0,0,0), 'f')
X = numpy.array( (1,0,0), 'f')
Y = numpy.array( (0,1,0), 'f')
Z = numpy.array( (0,0,1), 'f')

for angle in [1, 5, 10, 30, 45, 90, 180]:
    rotations['X'+str(angle)] = rotax( orig, X, angle*math.pi/180.)
    rotations['X-'+str(angle)] = rotax( orig, X, -angle*math.pi/180.)

for angle in [1, 5, 10, 30, 45, 90, 180]:
    rotations['Y'+str(angle)] = rotax( orig, Y, angle*math.pi/180.)
    rotations['Y-'+str(angle)] = rotax( orig, Y, -angle*math.pi/180.)

for angle in [1, 5, 10, 30, 45, 90, 180]:
    rotations['Z'+str(angle)] = rotax( orig, Z, angle*math.pi/180.)
    rotations['Z-'+str(angle)] = rotax( orig, Z, -angle*math.pi/180.)
