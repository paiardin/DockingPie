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
import types

def crossProduct (A, B, normal=True):
    """     Return cross product of two vectors A and B
normal: return normalized vector
"""
    res=[ A[1]*B[2] - A[2]*B[1],
          A[2]*B[0] - A[0]*B[2],
          A[0]*B[1] - A[1]*B[0] ]
    if normal:
        return norm(res)
    else:
        return res

def norm (A):
    """     Return normalized vector A.
"""
    if type(A) == types.ListType:
        A=numpy.array(A,'f')
        res= A/numpy.sqrt(numpy.dot(A,A))
        return res.tolist()    
    elif type(A)==numpy.ndarray:    
        return A/numpy.sqrt(numpy.dot(A,A))    
    else:
        print "Need a list or numpy array"
        return None

def getCenter(coords):
    """ get center of all the coords """
    coords=N.array(coords, 'f')    
    return (N.sum(coords)/len(coords)).tolist()
