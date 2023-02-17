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
# Last modified on Tue Sep  4 16:32:29 PDT 2001 by lindy
#
# $Header: /mnt/raid/services/cvs/python/packages/share1.5/mglutil/math/ncoords.py,v 1.2.18.1 2017/07/26 22:35:42 annao Exp $
#

"""ncoords.py - Numeric coordinates

This class is intented to be the base class of a number
of classes which transform and generally operate on lists
of homogeneous coordinates.
"""

import numpy

class Ncoords:
    def __init__(self, refCoords, tolist=1):
        """refCoords is an nx3 list of n points
        
        resultCoords is set up and maintained as homogeneous coords
        if tolist then return the result coords as a python list
        """
        try:
            self.refCoords = numpy.array(numpy.concatenate(
                (refCoords, numpy.ones( (len(refCoords), 1), 'f')), 1))
        except TypeError:
            raise ValueError, "invalid input array"

        self.resultCoords = self.refCoords
        self.tolist = tolist


    def reset(self):
        self.resultCoords = self.refCoords


    def getResultCoords(self):
        """Return the list of result coordinates

        if tolist is set, return an nx3 Python ListType.
        if tolist is not set, return an nx4 numpy array.
        """
        if self.tolist:
            return self.resultCoords[:,:3].tolist()
        else:
            return self.resultCoords


