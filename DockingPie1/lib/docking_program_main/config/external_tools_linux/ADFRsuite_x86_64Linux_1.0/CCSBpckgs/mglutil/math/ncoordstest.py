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
# Last modified on Tue Sep  4 17:02:59 PDT 2001 by lindy
#
# $Header: /mnt/raid/services/cvs/python/packages/share1.5/mglutil/math/ncoordstest.py,v 1.2.18.1 2017/07/26 22:35:42 annao Exp $
#

"""Unit test for ncoords.py

Requirements for ncoords.py:
    A. __init__:
        1. make Numeric, homogenious coordinates out of refCoords
        2. raise ValueError if refCoords is bad
    B. reset():
        3. nx3 slice of resultCoords must be equal to refCoords
    C. getResultCoords:
        4. return nx3 (not nx4) coordinates
        5. return as Numeric.array or ListType accorinding to self.tolist
"""

from mglutil.math.ncoords import Ncoords
import unittest, math
import numpy
import numpy.random as RandomArray


class NcoordsTest(unittest.TestCase):
    def setUp(self):
        """Called for every test."""

        npts = 500
        dim = 3
        self.max = 9999999.
        self.min = -self.max
        self.random_points = RandomArray.uniform(self.min,
                                                 self.max, (npts,dim)).tolist()


    def tearDown(self):
        pass



class InputOutputValues(NcoordsTest):


    def test_constructor_shape(self):
        """__init__         -- make refCoords and resultCoords homogeneous"""
        n = len(self.random_points)
        ncoords = Ncoords( self.random_points) ### tested call ###
        # confirm shape to be nx4
        self.assertEqual( (n, 4), numpy.shape(ncoords.resultCoords))
        self.assertEqual( (n, 4), numpy.shape(ncoords.refCoords))
        # cofirm that the last column is all ones
        self.assertEqual(numpy.ones(n).tolist(),
                         ncoords.resultCoords[:,3].tolist())
        self.assertEqual(numpy.ones(n).tolist(),
                         ncoords.refCoords[:,3].tolist())


    def test_input_error(self):
        """__init__         -- ValueError on bad input"""
        self.assertRaises(ValueError, Ncoords, range(10))
        self.assertRaises(ValueError, Ncoords, [(1,1,1),(1,1)] )


    def test_reset_values(self):
        """reset            -- points equal input values after reset"""
        nc = Ncoords( self.random_points, tolist=1)
        nc.reset() ### tested call ###
        result = nc.getResultCoords()
        # compare input and output point lists
        self.assertEqual( self.random_points, result)


    def test_getResultCoords_shape(self):
        """getResultCoords  -- if tolist: return nx3 ListType"""
        n = len(self.random_points)
        nc = Ncoords(self.random_points, tolist=0)
        nc.tolist=1
        result = nc.getResultCoords() ### tested call ###
        # confirm shape
        self.assertEqual((n, 3), numpy.shape(result))
        # confirm type
        self.assertEqual(type([]), type(result))


    def test_getResultCoords_type(self):
        """getResultCoords  -- if not tolist: return nx4 numpy.array"""
        n = len(self.random_points)
        nc = Ncoords(self.random_points, tolist=1)
        nc.tolist=0
        result = nc.getResultCoords() ### tested call ###
        # confirm shape
        self.assertEqual((n, 4), numpy.shape(result))
        # confirm type
        self.assertEqual(type(numpy.array([])), type(result))



if __name__ == '__main__':
    unittest.main()   

# for example:
#     py mglutil/math/ncoordstest.py -v
# or, to redirect output to a file:
#     py ncoordstest.py -v > & ! /tmp/nct.out






