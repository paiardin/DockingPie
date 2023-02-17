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
# $Header: /mnt/raid/services/cvs/DejaVu2/moveGeom.py,v 1.1.4.1 2017/07/13 22:28:33 annao Exp $
#
# $Id: moveGeom.py,v 1.1.4.1 2017/07/13 22:28:33 annao Exp $
#
from math import sqrt, cos, sin, acos, fabs, pi
import numpy

def matToQuaternion(mat):
    # converts rotation matrix to quaternion
    # http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm
    trace = 1 + mat[0] + mat[5] + mat[10]
    if trace > 0.00000001:
      S = sqrt(trace) * 2
      X = ( mat[9] - mat[6] ) / S
      Y = ( mat[2] - mat[8] ) / S
      Z = ( mat[4] - mat[1] ) / S
      W = 0.25 * S
    else:
        if  mat[0] > mat[5] and mat[0] > mat[10]:      
            S  = sqrt( 1.0 + mat[0] - mat[5] - mat[10] ) * 2
            X = 0.25 * S
            Y = (mat[4] + mat[1] ) / S
            Z = (mat[2] + mat[8] ) / S
            W = (mat[9] - mat[6] ) / S
        elif mat[5] > mat[10] : 
            S  = sqrt( 1.0 + mat[5] - mat[0] - mat[10] ) * 2
            X = (mat[4] + mat[1] ) / S
            Y = 0.25 * S
            Z = (mat[9] + mat[6] ) / S
            W = (mat[2] - mat[8] ) / S
        else:		       
            S  = sqrt( 1.0 + mat[10] - mat[0] - mat[5] ) * 2
            X = (mat[2] + mat[8] ) / S
            Y = (mat[9] + mat[6] ) / S
            Z = 0.25 * S
            W = (mat[4] - mat[1] ) / S
    
      #The quaternion is then defined as:
      #Q = | X Y Z W |
##    #normalize:
##    n = sqrt(X*X + Y*Y + Y*Y + W*W);
##    X /= n
##    Y /= n
##    Z /= n
##    W /= n
    return [X, Y, Z, W]


def quatToMatrix(q):
    # converts quaternion to matrix
    # http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToMatrix/index.htm
    x, y,z, w = q
    sqw = w*w
    sqx = x*x
    sqy = y*y
    sqz = z*z

    #invs (inverse square length) is only required if quaternion is not already normalised
    invs = 1 / (sqx + sqy + sqz + sqw)
    m = numpy.zeros((4,4), 'f')
    m[0][0] = ( sqx - sqy - sqz + sqw)*invs  #since sqw + sqx + sqy + sqz =1/invs*invs
    m[1][1] = (-sqx + sqy - sqz + sqw)*invs 
    m[2][2] = (-sqx - sqy + sqz + sqw)*invs 
    
    tmp1 = x*y
    tmp2 = z*w
    m[1][0] = 2.0 * (tmp1 + tmp2)*invs 
    m[0][1] = 2.0 * (tmp1 - tmp2)*invs 
    
    tmp1 = x*z
    tmp2 = y*w
    m[2][0] = 2.0 * (tmp1 - tmp2)*invs 
    m[0][2] = 2.0 * (tmp1 + tmp2)*invs 
    tmp1 = y*z
    tmp2 = x*w
    m[2][1] = 2.0 * (tmp1 + tmp2)*invs 
    m[1][2] = 2.0 * (tmp1 - tmp2)*invs
    mat = m.flatten()
    mat[-1] = 1.
    return mat

def interpolateQuaternion(firstVal, lastVal, fraction):
    qa = firstVal
    qb = lastVal
    #invertVal = False
    #t = self.ease(fraction, interval)
    cosHalfTheta =  qa[0] * qb[0] + qa[1] * qb[1] + qa[2] * qb[2] + qa[3] * qb[3]
    #print "fraction:", fraction, "cosHalfTheta: ", cosHalfTheta
    #we need cosHalfTheta to be positive - invert the quaternion
    if cosHalfTheta < 0:  # ????
        qb[0] = -qb[0]; qb[1] = -qb[1]; qb[2] = -qb[2]; qb[3] = -qb[3] 
        cosHalfTheata = -cosHalfTheta
        #invertVal = True
    if abs(cosHalfTheta) >= 1.0:
        qm = [qa[0], qa[1], qa[2], qa[3]]
        return qm
    #Calculate temporary values.

    halfTheta = acos(cosHalfTheta)

    sinHalfTheta = sqrt(1.0 - cosHalfTheta*cosHalfTheta)
    # we could rotate around any axis normal to qa or qb
    if fabs(sinHalfTheta) < 0.001: #fabs is floating point absolute
        qm = [(qa[0] * 0.5 + qb[0] * 0.5), (qa[1] * 0.5 + qb[1] * 0.5),
              (qa[2] * 0.5 + qb[2] * 0.5), (qa[3]* 0.5 + qb[3] * 0.5)]
        return qm
    ratioA = sin((1 - fraction) * halfTheta) / sinHalfTheta
    ratioB = sin(fraction * halfTheta) / sinHalfTheta
    #calculate Quaternion.
    qm = [(qa[0] * ratioA + qb[0] * ratioB), (qa[1] * ratioA + qb[1] * ratioB),
          (qa[2] * ratioA + qb[2] * ratioB), (qa[3] * ratioA + qb[3] * ratioB)]
    return qm

def interpolateVector(firstVal, lastVal, fraction):
    if fraction <= 0.0:
        return firstVal
    elif fraction >= 1.0:
        return lastVal
    else:
        value = []
        for v1,v2 in zip(firstVal, lastVal):
            value.append(v1+fraction*(v2-v1))
        return value

def interpolateScalar(firstVal, lastVal, fraction):
    if fraction <= 0.0:
        return firstVal
    elif fraction >= 1.0:
        return lastVal
    else:
        valueRange = lastVal - firstVal
    return firstVal + valueRange*fraction

def comparefloats(a, b, precision = 0.0001 ):
    """Compare two float scalars or arrays"""

    aa = numpy.ravel(numpy.array(a))
    bb = numpy.ravel(numpy.array(b))
    if len(aa) == len(bb) == 0: return True
    if len(aa) != len(bb): return False
    diff = numpy.abs(aa-bb)
    if diff.max() > precision: return False

class MoveGeom:
    """This class is designed to:
    - record current transformation(rotation, translation, scale, pivot) of the specified geometry;
    - record camera attributes (fieldOfView, lookFrom, near, far, fog)
    - move the specified object from one position (defined by transformation and camera attributes)
      to another by interpolating the rotation, translation, etc values; 
    """
    
    def __init__(self, obj, rendering=None):
        """
        - obj is a DejaVu Transformable object that belongs to a DejaVu Viewer
        - rendering is a dictionary storing the state of all geometries;
        """
        self.object = obj
        #self.cameraAttributes, self.transformation = self.getCurrentOrient(orient)
        self.rendering = rendering
        self.viewer = self.object.viewer

    def getTransformation(self, obj):
        rotation = matToQuaternion(obj.rotation)
        return [rotation, obj.translation[:], obj.scale[:], obj.pivot[:] ]

    def getCameraAttr(self):
        cam = self.object.viewer.currentCamera
        fieldOfView = cam.fovy
        lookFrom = cam.lookFrom.copy()
        cameraAttributes = {}
        cameraAttributes['fieldOfView'] = fieldOfView
        cameraAttributes['lookFrom'] = lookFrom
        cameraAttributes['near'] = cam.near
        cameraAttributes['far'] = cam.far
        cameraAttributes['fogStart'] = cam.fog.start
        cameraAttributes['fogEnd'] = cam.fog.end
        return cameraAttributes

    def getCurrentOrient(self):
        """
        """
        cameraAttributes = self.getCameraAttr()
        transformation = {}
        for _obj in self.object.AllObjects():
            #if not obj.animatable: continue
            #oname = obj.fullName
            try:
                transformation[_obj] = self.getTransformation(_obj)
            except:
                pass
        return cameraAttributes, transformation

    def setTransformation(self, obj, value):
        rotation = quatToMatrix(value[0])
        redo = False
        if obj != obj.viewer.rootObject:
            # need to rebuild display list if the object is not root 
             redo = True
        obj.Set(rotation=rotation, translation=value[1], scale=value[2],
                pivot=value[3], redo=redo)

    def setCameraAttr(self, cameraAttr):
        cam = self.viewer.currentCamera
        cam._setFov(cameraAttr['fieldOfView'])
        cam.Set(near=cameraAttr['near'])
        cam.Set(far=cameraAttr['far'])
        cam.Set(lookFrom=cameraAttr['lookFrom'])
        cam.fog.Set(start=cameraAttr['fogStart'], end=cameraAttr['fogEnd'])

    def interTransformation(self, obj, transf1, transf2, fraction):
        if transf1.has_key(obj) and transf2.has_key(obj):
            # rotation
            rot1 = transf1[obj][0]
            rot2 = transf2[obj][0]
            qt = interpolateQuaternion(rot1, rot2, fraction)
            #translation
            tr1 = transf1[obj][1]
            tr2 = transf2[obj][1]
            tran = interpolateVector(tr1, tr2, fraction)
            #scale
            sc1 = transf1[obj][2]
            sc2 = transf2[obj][2]
            scale = interpolateVector(sc1, sc2, fraction)
            #pivot
            pt1 = transf1[obj][3]
            pt2 = transf2[obj][3]
            pivot = interpolateVector(pt1, pt2, fraction)
            return qt, tran, scale, pivot

    def interCameraAttributes(self, ca1, ca2, fraction):
        cameraAttr = {}
        #camera
        cameraAttr['far'] = interpolateScalar(ca1['far'], ca2['far'], fraction)
        cameraAttr['fieldOfView'] = interpolateScalar(ca1['fieldOfView'], ca2['fieldOfView'], fraction)
        cameraAttr['fogEnd'] = interpolateScalar(ca1['fogEnd'], ca2['fogEnd'], fraction)
        cameraAttr['fogStart'] = interpolateScalar(ca1['fogStart'], ca2['fogStart'], fraction)
        cameraAttr['lookFrom'] = interpolateVector(ca1['lookFrom'], ca2['lookFrom'], fraction)
        cameraAttr['near'] = interpolateScalar(ca1['near'], ca2['near'], fraction)
        #print "CAMERA FAR at:", fraction, cameraAttr['far']
        return cameraAttr

    def interpolate(self, interpolationTime, obj=None, transf=None, camattr=None, interAllobj=False):
        """
        Interpolate between two object's positions.
        obj --- a Transformable geometry. If not specified(None) - self.object is used)
        tranf  --- a list of two dictionaries of object(s)transformations;
        camattr--- a list of two dictionaries for camera attributes;
        interAllobj --- if True - interpolate the transf values for the children of the specified object as well.
        """
        #self.viewer.stopAutoRedraw()
        self.viewer.suspendRedraw = True
        if obj is None: obj = self.object
        if interAllobj:
            allobjects = obj.AllObjects()
        else:
            allobjects = [obj]

        if transf is not None:
            transf1 = transf[0]
            transf2 = transf[1]
        else:
            transf1 = transf2 = None
        if camattr is not None:
            ca1 = camattr[0]
            ca2 = camattr[1]
        else:
            ca1 = ca2 = None
            
        def _redraw(transf1, transf2, ca1, ca2, fraction):
            if transf1 is not None:
                for _obj in allobjects:
                    qt, tran, scale, pivot = self.interTransformation(_obj, transf1, transf2, fraction)
                    #print tran, scale, pivot
                # set transformation
                self.setTransformation(_obj, [qt, tran, scale, pivot])
            # camera attributes
            if ca1 is not None:
                cameraAttr = self.interCameraAttributes(ca1, ca2, fraction)
                self.setCameraAttr(cameraAttr)
            self.viewer.OneRedraw()

        from time import time
        timeLeft = interpolationTime
        t0 = time()
        self.viewer.suspendRedraw = False
        _redraw(transf1, transf2, ca1, ca2, 0.0)
        redrawTime = time()-t0
        timeLeft -= redrawTime
        fraction = 0.0
        while timeLeft > 0:
            fraction = (interpolationTime-timeLeft) / interpolationTime
            t0 = time()
            _redraw(transf1, transf2, ca1, ca2, fraction)
            redrawTime = time()-t0
            timeLeft -= redrawTime
        if fraction < 1.0:
            fraction = 1.0
            _redraw(transf1, transf2, ca1, ca2, fraction)
        #self.viewer.startAutoRedraw()
