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
# Copyright: M. Sanner TSRI 2016
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/PmvApp/cubicInterpolate.py,v 1.1.4.1 2017/07/13 20:52:52 annao Exp $
# 
# $Id: cubicInterpolate.py,v 1.1.4.1 2017/07/13 20:52:52 annao Exp $
#
import numpy

def CubicInterpolate(p0, p1, p2, p3, t):
    t2 = t*t
    a0 = [p3[0] - p2[0] - p0[0] + p1[0],
          p3[1] - p2[1] - p0[1] + p1[1],
          p3[2] - p2[2] - p0[2] + p1[2]]
    a1 = [p0[0] - p1[0] - a0[0],
          p0[1] - p1[1] - a0[1],
          p0[2] - p1[2] - a0[2]]
    a2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]] 

    return [ a0[0] * t * t2 + a1[0] * t2 + a2[0] * t + p1[0],
             a0[1] * t * t2 + a1[1] * t2 + a2[1] * t + p1[1],
             a0[2] * t * t2 + a1[2] * t2 + a2[2] * t + p1[2]]

def ResampleControlPoints(controlPoints, nbPointsPerRes, identity):
    # for N control points we generate (N-3)*nbPointsPerRes smoth points
    #
    # for example for 3 CA we have 2 dummy added which corresponds to
    # 2 intervals between CA atoms
    #
    #  X - O - O - 0 - X
    nP = len(controlPoints)
    resampledControlPoints = []
    currentPointId = 1
    currentPosition = controlPoints[currentPointId]
    ##distance = DisplaySettings.Instance.DistanceContraint;
    lerpValues = numpy.arange(0.0, nbPointsPerRes)/nbPointsPerRes
    halfLerp = len(lerpValues)/2
    ## Normalize the distance between control points
    interpIds = []
    while True:
        if currentPointId + 2 >= nP:
            break
        for nl, lerpValue in enumerate(lerpValues):
            candidate = CubicInterpolate(controlPoints[currentPointId-1],
                                         controlPoints[currentPointId],
                                         controlPoints[currentPointId+1],
                                         controlPoints[currentPointId+2],
                                         lerpValue)
            resampledControlPoints.append(candidate)
            if nl<halfLerp:
                interpIds.append(identity[currentPointId-1])
            else:
                interpIds.append(identity[currentPointId])
        currentPointId += 1
    
    return numpy.array(resampledControlPoints), numpy.array(interpIds)

from math import sqrt
def Normalize(v):
    n = 1.0/sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
    return v[0]*n, v[1]*n, v[2]*n

def Cross(A, B):
    return [A[1]*B[2]-B[1]*A[2],
            A[2]*B[0]-B[2]*A[0],
            A[0]*B[1]-B[0]*A[1]]

def GetSmoothNormals(controlPoints):
    controlPoints = numpy.array(controlPoints)
    smoothNormals = numpy.zeros( (len(controlPoints), 3), 'd')
    p0 = controlPoints[0]
    p1 = controlPoints[1]
    p2 = controlPoints[2]
    smoothNormals[0][:] = Normalize(Cross(p0 - p1, p2 - p1))
    for i in range(1, len(controlPoints)-1):
        p0 = controlPoints[i - 1]
        p1 = controlPoints[i]
        p2 = controlPoints[i + 1]
        t = Normalize(p2 - p0);
        b = Normalize(Cross(t, smoothNormals[i-1]))
        n = Normalize(Cross(b, t))
        smoothNormals[i][:] = n[0], n[1], n[2]

    smoothNormals[-1][:] = smoothNormals[-2][:]
    return smoothNormals

def GetFrames(path):
    normals = GetSmoothNormals(path)

    caca1 = numpy.zeros( normals.shape, 'd' )
    caca1 = path[1:]-path[:-1]
    caca1 = caca1/numpy.linalg.norm(caca1,axis=1).reshape((len(caca1),1))
    # double last element
    caca1 = numpy.array( caca1.tolist() + [caca1[-1]] )
    
    binormals = numpy.cross(caca1, normals)
    binormals = binormals/numpy.linalg.norm(binormals,axis=1).reshape((len(binormals),1))
    
    return caca1, normals, binormals
