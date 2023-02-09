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
# Last modified on Thu Jan 31 10:27:11 PST 2002 by lindy
#
# $Header: /mnt/raid/services/cvs/python/packages/share1.5/mglutil/math/rmsdtest.py,v 1.4.18.1 2017/07/26 22:35:42 annao Exp $
#
"""Unit test for rmsd.py

Requirements for rmsd:
   A. RMSDCalculator.__init__
      0. should ..
   B. RMSDCalculator.setRefCoords
      0. should ..
   C. RMSDCalculator.computeRMSD
      1. should return known result with known input
      2. raise ValueError for input of unlike dimensions
      3. for two random sets of points, rmsd(x,y) == rmsd(y,x)
      4. raise ValueError if the reference coords have not been set
   D.

"""

from mglutil.math import rmsd
import unittest, math
import numpy
import numpy.random as RandomArray

class ComputedValues(unittest.TestCase):
    decimals = 4 # decimal places to round to for float comparison
    point_list_0 = numpy.zeros((5,3))
    point_list_1 = numpy.ones( (5,3))

    knowValues = ( (point_list_0, point_list_0, 0.0),
                   (point_list_1, point_list_1, 0.0),
                   (point_list_0, point_list_1, math.sqrt(3.0)),
                   (point_list_1, point_list_0, math.sqrt(3.0)))

    def test_computeRMSD_KnowValues(self):
        """1. should return known result with known input"""
        for ref, input, known in self.knowValues:
            self.assertEqual(known,
                             rmsd.RMSDCalculator(ref).computeRMSD(input))


    def test_computeRMSD_RandomOffset(self):
        """5. offset point by random value returns offset*sqrt(3)"""
        min = -10000.
        max = 10000.
        num_points = 20
        dimension = 3
        point_list_1 = RandomArray.uniform(min, max, (num_points, dimension))
        delta = point_list_1[0][0]
        point_list_2 = point_list_1 + delta
        answer = rmsd.RMSDCalculator(point_list_1).computeRMSD(point_list_2)
        self.assertEqual(
            round(answer, self.decimals),
            round(abs(delta)*math.sqrt(3.0), self.decimals))


    def test_computeRMSD_Random(self):
        """3. for two random sets of points, rmsd(x,y) == rmsd(y,x)"""
        min = -10000.
        max = 10000.
        num_points = 20
        dimension = 3
        point_list_1 = RandomArray.uniform(min, max, (num_points, dimension))
        point_list_2 = RandomArray.uniform(min, max, (num_points, dimension))
        self.assertEqual(
            rmsd.RMSDCalculator(point_list_1).computeRMSD(point_list_2),
            rmsd.RMSDCalculator(point_list_2).computeRMSD(point_list_1))
            

class InputValues(unittest.TestCase):
    point_list_0 = numpy.zeros((3,3))
    point_list_1 = numpy.ones( (4,3)) # different lengths

    def test_computeRMSD_dimensions(self):
        """2. raise ValueError for input of unlike dimensions"""
        ruler = rmsd.RMSDCalculator(self.point_list_0)
        self.assertRaises(ValueError, ruler.computeRMSD, self.point_list_1)

    def test_computeRMSD_noRefCoords(self):
        """4. raise ValueError if the reference coords have not been set"""
        ruler = rmsd.RMSDCalculator()
        self.assertRaises(ValueError, ruler.computeRMSD, self.point_list_1)

if __name__ == "__main__":
    unittest.main()   

# for example: py mglutil/math/rmsdtest.py -v



