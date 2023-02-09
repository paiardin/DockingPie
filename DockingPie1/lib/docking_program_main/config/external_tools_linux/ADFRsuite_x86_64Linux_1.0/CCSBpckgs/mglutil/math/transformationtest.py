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
# Last modified on Wed Aug 22 14:22:48 PDT 2001 by lindy
#
# $Header: /mnt/raid/services/cvs/python/packages/share1.5/mglutil/math/transformationtest.py,v 1.3.18.1 2017/07/26 22:35:42 annao Exp $
#

"""Unit test for transformation.py

Transformation.py defines three classes:
    Quaternion
    UnitQuaternion(Quaternion)
    Transformation(UnitQuaternion)

"""

from mglutil.math.transformation import UnitQuaternion, Quaternion
from mglutil.math.transformation import Transformation
import unittest

from numpy import random as RA
import math

#
# Unit tests for the Quaternion class
#

class QuaternionTest(unittest.TestCase):
    def setUp(self):
        """The Quaternion class is tested through the UnitQuaternion class"""
        pass



class UnitQuaternionTest(unittest.TestCase):
    def setUp(self):
        self.max = 999999999.
        self.min = -self.max

        

class UnitQuaternionKnownValues(UnitQuaternionTest):
##          theta = 360.0*RA.random()
##          knownValues = (
##              # identity
##              ((0., 0., 0., 0.), ( (1.0, 0.0, 0.0, 0.0),
##                                   (0.0, 1.0, 0.0, 0.0),
##                                   (0.0, 0.0, 1.0, 0.0),
##                                   (0.0, 0.0, 0.0, 1.0))),
##              # rotation about x
##              ((0., 0., 0., 0.), ( (1.0, 0.0, 0.0, 0.0),
##                                   (0.0, 1.0, 0.0, 0.0),
##                                   (0.0, 0.0, 1.0, 0.0),
##                                   (0.0, 0.0, 0.0, 1.0))),
##              # rotation about y
##              ((0., 0., 0., 0.), ( (1.0, 0.0, 0.0, 0.0),
##                                   (0.0, 1.0, 0.0, 0.0),
##                                   (0.0, 0.0, 1.0, 0.0),
##                                   (0.0, 0.0, 0.0, 1.0))),
##              # rotation about z
##              ((0., 0., 0., 0.), ( (1.0, 0.0, 0.0, 0.0),
##                                   (0.0, 1.0, 0.0, 0.0),
##                                   (0.0, 0.0, 1.0, 0.0),
##                                   (0.0, 0.0, 0.0, 1.0))),
##              ((0., 0., 0., 0.), ( (1.0, 0.0, 0.0, 0.0),
##                                   (0.0, 1.0, 0.0, 0.0),
##                                   (0.0, 0.0, 1.0, 0.0),
##                                   (0.0, 0.0, 0.0, 1.0))),
##              ((0., 0., 0., 0.), ( (1.0, 0.0, 0.0, 0.0),
##                                   (0.0, 1.0, 0.0, 0.0),
##                                   (0.0, 0.0, 1.0, 0.0),
##                                   (0.0, 0.0, 0.0, 1.0))),
##              ((0., 0., 0., 0.), ( (1.0, 0.0, 0.0, 0.0),
##                                   (0.0, 1.0, 0.0, 0.0),
##                                   (0.0, 0.0, 1.0, 0.0),
##                                   (0.0, 0.0, 0.0, 1.0))),
##              ((0., 0., 0., 0.), ( (1.0, 0.0, 0.0, 0.0),
##                                   (0.0, 1.0, 0.0, 0.0),
##                                   (0.0, 0.0, 1.0, 0.0),
##                                   (0.0, 0.0, 0.0, 1.0))),
##              ((0., 0., 0., 0.), ( (1.0, 0.0, 0.0, 0.0),
##                                   (0.0, 1.0, 0.0, 0.0),
##                                   (0.0, 0.0, 1.0, 0.0),
##                                   (0.0, 0.0, 0.0, 1.0))))

    def tearDown(self):
        pass


    def testKnown00(self):
        """<describe what's being tested here>"""
        pass



class UnitQuaternionProperties(UnitQuaternionTest):

    def testProperties00(self):
        """The product of the conjugate is the conjucate of the product"""
        q1 = Quaternion(tuple(RA.uniform(self.min,
                                         self.max, (1, 4)).tolist()[0]))
        q2 = Quaternion(tuple(RA.uniform(self.min,
                                         self.max, (1, 4)).tolist()[0]))
##          pc = q1.conjugate()*q2.conjugate()
##          qp = q1*q2
##          cp = qp.conjugate()
##          self.assertEqual( pc, cp)
        # the commented code fails with the same error as this line...
        self.assertEqual( q1.conjugate()*q2.conjugate(), (q2*q1).conjugate())
        

    def testProperties01(self):
        """The magnitudue of the product is the product of the magnitudes"""
        q1 = Quaternion(tuple(RA.uniform(self.min,
                                         self.max, (1, 4)).tolist()[0]))
        q2 = Quaternion(tuple(RA.uniform(self.min,
                                         self.max, (1, 4)).tolist()[0]))
        self.assertEqual((q1*q2).magnitude(), q1.magnitude()*q2.magnitude())


    def testProperties02(self):
        """The conjugate of a unit quaternion is it's inverse"""
        q1 = Quaternion(tuple(RA.uniform(self.min,
                                         self.max, (1, 4)).tolist()[0]))
        self.assertEqual( q1.conjugate(), q1.inverse())

    


class TransformationTest(unittest.TestCase):
    def setUp(self):
        pass


    def tearDown(self):
        pass

    def test(self):
        pass



if __name__ == '__main__':
    unittest.main()   

# for example: py mglutil/math/transformationtest.py -v


