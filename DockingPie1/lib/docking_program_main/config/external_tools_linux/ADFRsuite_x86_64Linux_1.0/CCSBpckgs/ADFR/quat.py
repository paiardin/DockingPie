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
# Copyright: M. Sanner and TSRI 2015
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/ADFR/quat.py,v 1.2 2016/12/07 00:38:32 sanner Exp $
#
# $Id: quat.py,v 1.2 2016/12/07 00:38:32 sanner Exp $
#

from math import sqrt, sin, cos, acos, pi
TWOPI = 2*pi
from random import uniform

def randomQuat():
    t1 = uniform(0., TWOPI)
    x0 = uniform(0., 1.)
    r1 = sqrt( 1. - x0 )
    x = sin( t1 ) * r1
    y = cos( t1 ) * r1
    t2 = uniform(0., TWOPI)
    r2 = sqrt( x0 )
    z = sin( t2 ) * r2
    w = cos( t2 ) * r2
    return x, y, z, w

APPROX_ZERO=1.0E-6

def normAxisAngle( x, y, z, a ):
    """
    Normalise the 3D rotation axis or vector nx,ny,nz
    """
    mag3 = sqrt( x*x + y*y + z*z)
    if mag3 > APPROX_ZERO:
        inv_mag3 = 1. / mag3
        nx = inv_mag3 * x
        ny = inv_mag3 * y
        nz = inv_mag3 * z
        return nx, ny,nz, a
    else:
        return 1., 0., 0., 0.
    
def quatToAxisAngle( x, y, z, w ):
    """
    Convert the quaternion components (x,y,z,w) of the quaternion q,
    to the corresponding rotation-about-axis components (nx,ny,nz,ang)
    """
    if w >= 1. or w <= -1.: return 1., 0., 0., 0. # axis angle identity

    angle = 2. * acos(w)
    inv_sin_half_angle =  1. / sin( angle / 2. )
    nx = x * inv_sin_half_angle
    ny = y * inv_sin_half_angle
    nz = z * inv_sin_half_angle
    ## by convention, angles should be in the range -PI to +PI.
    if angle > pi:
        angle = (angle%TWOPI)-TWOPI
    elif angle < -pi:
        angle = (angle%TWOPI)+TWOPI
    return normAxisAngle( nx, ny, nz, angle )

def axisAngleToQuat( x, y, z, a ):
    """
    Normalize the rotation-about-axis vector 
    and convert the rotation-about-axis components (nx,ny,nz,ang)
    to the corresponding quaternion components (x,y,z,w)
    Originally was named convertRotToQuat( )
    """
    nmag = sqrt( x*x + y*y + z*z )
    if nmag <= APPROX_ZERO: return 1., 0., 0., 0. # error

    if a > pi:
        a = (a%TWOPI)-TWOPI
    elif a < -pi:
        a = (a%TWOPI)+TWOPI

    inv_nmag = 1. / nmag
    hqang = 0.5 * a
    s     = sin( hqang )
    nx = s * x * inv_nmag #      /* Normalize axis */
    ny = s * y * inv_nmag #      /* Normalize axis */
    nz = s * z * inv_nmag #      /* Normalize axis */
    nw  = cos( hqang )
    return nx, ny, nz, nw

def randomQuatByAmount(amount):
    """
    returns a quaternion from a random axis and specified angle
    amount is an angle in radians
    """
    q = randomQuat()
    x, y, z, a = quatToAxisAngle( *q ) # will not be 3-element normalized
    return axisAngleToQuat( x, y, z, amount )

def qmult(ql, qr):
    x = (ql[3]*qr[0] + ql[0]*qr[3] + ql[1]*qr[2] - ql[2]*qr[1])
    y = (ql[3]*qr[1] + ql[1]*qr[3] + ql[2]*qr[0] - ql[0]*qr[2])
    z = (ql[3]*qr[2] + ql[2]*qr[3] + ql[0]*qr[1] - ql[1]*qr[0])
    w = (ql[3]*qr[3] - ql[0]*qr[0] - ql[1]*qr[1] - ql[2]*qr[2])
    return x, y, z, w

def conjugate(x, y, z, w):
    """
    compute quaternion conjugate
    """
    return -z, -y, -z, w

def qinv(x, y, z, w):
    """
    compute quaternion inverse
    """
    inv_squared_magnitude = 1. / ( x*x + y*y + z*z + w*w )

    invx = -x * inv_squared_magnitude
    invy = -y * inv_squared_magnitude
    invz = -z * inv_squared_magnitude
    invw = w * inv_squared_magnitude

    return invx,invy,invz,invw

def qnorm(x, y, z, w):
    """
    normalize quaternion
    """
    n1 = 1./sqrt(x*x + y*y + z*z + w*w)
    return x*n1, y*n1, z*n1, w*n1
