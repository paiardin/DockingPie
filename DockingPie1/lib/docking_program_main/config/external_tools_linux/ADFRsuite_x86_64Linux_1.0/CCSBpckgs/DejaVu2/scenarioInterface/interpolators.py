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

from Scenario2.interpolators import VarVectorInterpolator,CompositeInterpolator, \
     FloatVectorInterpolator, RotationInterpolator, FloatScalarInterpolator

class MaterialInterpolator(CompositeInterpolator):
    nbvar = 4
    def __init__(self, firstVal, lastVal, interpolation='linear', interpolators=None, active = True):
        if not interpolators:
            interpolators = [VarVectorInterpolator, # for ambient RGB
                          VarVectorInterpolator, # for specular RGB
                          VarVectorInterpolator, # for emissive RGB
                          FloatVectorInterpolator, # for shininess
                          ]

        
        CompositeInterpolator.__init__(self, firstVal, lastVal, interpolators=interpolators,
                                       interpolation=interpolation)

class LightColorInterpolator(CompositeInterpolator):
    nbvar = 3
    def __init__(self, firstVal, lastVal, interpolation='linear', interpolators=None, active = True):
        if not interpolators:
            ## interpolators = [VarVectorInterpolator, # for ambient
##                              VarVectorInterpolator, # for diffuse
##                              VarVectorInterpolator, # for specular
##                              ]
            
            interpolators = [FloatVectorInterpolator, # for ambient
                             FloatVectorInterpolator, # for diffuse
                             FloatVectorInterpolator, # for specular
                             ]
        CompositeInterpolator.__init__(self, firstVal, lastVal, interpolators=interpolators,
                                       interpolation=interpolation)


class TransformInterpolator(CompositeInterpolator):
    nbvar = 4
    def __init__(self, firstVal, lastVal, interpolation='linear', interpolators=None, active = True):
        if not interpolators:
            interpolators = [RotationInterpolator, # for rotation
                          FloatVectorInterpolator, # for translation
                          FloatVectorInterpolator, # for scale
                          FloatVectorInterpolator# for pivot
                          ]
        
        CompositeInterpolator.__init__(self, firstVal, lastVal, interpolators=interpolators,
                                       interpolation=interpolation)



class RotationConcatInterpolator(FloatScalarInterpolator):

    def getValue(self, fraction, interval):
        # problem .. how to store last value ==> who calls Interval._getValue2KF
        easeIn = False
        easeOut = False
        if self.behaviorList is not None:
            bl = self.behaviorList()
            #print "Rtation ,easeIn-Out",  bl.easeIn , bl.easeOut
            easeIn = bl.easeIn
            easeOut = bl.easeOut

        if not easeIn and not easeOut:
            return FloatScalarInterpolator.getValue(self, fraction, interval)

        fullValue = interval.data['fullValue']
        if fraction == 0:
            interval._lastSetValue = 0

        _d = self.ease(fraction, interval)*fullValue # whole distance from time 0 to time t
        if interval._lastSetValue:
            dd = _d-interval._lastSetValue               # distance increment
        else:
            dd = _d
        #print fraction, dd , _d, interval.data['fullValue'], interval._lastSetValue
        #sum = sum + dd          # check the distance
        interval._lastSetValue = _d
        return dd


