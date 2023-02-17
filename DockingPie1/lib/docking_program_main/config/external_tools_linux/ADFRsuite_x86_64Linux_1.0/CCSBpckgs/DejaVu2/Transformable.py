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
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/DejaVu2/Transformable.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#
# $Id: Transformable.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#

from opengltk.OpenGL import GL
from opengltk.OpenGL.GLU import gluProject
from opengltk.extent.utillib import glCleanRotMat
from opengltk.extent import _gllib as gllib

import numpy, math

import DejaVu2
from mglutil.events import Event

class translateEvent(Event):
    def __init__(self, arg=None, objects=[]):
        """  """
        self.arg = arg
        self.objects = objects


class Transformable(object):
    """Base Class inherited by objects which can be transformed using the
       mouse"""

    def multMat4pt(self, mat, pt):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        ptx = mat[0][0]*pt[0]+mat[0][1]*pt[1]+mat[0][2]*pt[2]+mat[0][3]
        pty = mat[1][0]*pt[0]+mat[1][1]*pt[1]+mat[1][2]*pt[2]+mat[1][3]
        ptz = mat[2][0]*pt[0]+mat[2][1]*pt[1]+mat[2][2]*pt[2]+mat[2][3]
        return (ptx, pty, ptz)


    def ResetTransformation(self, redo=1):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Reset the tranformations (Rotation, translation, pivot, scale)"""
	self.rotation = numpy.identity(4).astype('f')
	self.rotation.shape = (16, )
	self.translation = numpy.zeros( (3,), 'f')
	self.pivot = numpy.zeros( (3,), 'f')
	self.scale = numpy.ones( (3,), 'f')
        if self.viewer:
            if self.viewer.currentObject != self.viewer.rootObject \
               and redo:
                self.viewer.deleteOpenglList()


    def __init__(self, viewer=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Constructor"""

        self.inheritXform = 1                    # set to 0 not to inherit
	self.redirectXform = None                # object to which transf.
	                                         # should be redirected
	self.propagateRedirection = 1

	self.copyXform = []                      # list of objects to be
                                                 # transf. along with me

                                                 # Object's original transf.
                                                 # in OpenGL form (shape (16,))
        self.Matrix = numpy.identity(4).astype('f')
        self.Matrix.shape = (16, )
        self.MatrixRot = self.Matrix
        self.MatrixRotInv = self.Matrix
        self.MatrixScale = numpy.ones((3,), 'f')
        self.MatrixTransl = numpy.zeros( (3,), 'f')
        self.viewer = viewer
        
	self.ResetTransformation(redo=0)               # init object's transf.

	self.R = numpy.identity(4).astype('f') # Object's frame rotation
	self.R.shape = (16, )

	self.Ri = numpy.identity(4) .astype('f') # Inverse of R
	self.Ri.shape = (16, )

	self.Si = numpy.ones( (3, ) ).astype('f') # Inverse of frame's scale
	self.isScalable = 1

        self.immediateRendering = False
        #self.hasChildWithImmediateRendering = False # set to True if a child is not using dpyList
        self.needsRedoDpyListOnResize = False


    def Decompose4x4(self, matrix, cleanup=True):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """ takes a matrix in shape (16,) in OpenGL form (sequential values go
        down columns) and decomposes it into its rotation (shape (16,)),
        translation (shape (3,)), and scale (shape (3,)) """
        m = matrix
        transl = numpy.array((m[12], m[13], m[14]), 'f')
        scale0 = numpy.sqrt(m[0]*m[0]+m[4]*m[4]+m[8]*m[8])
        scale1 = numpy.sqrt(m[1]*m[1]+m[5]*m[5]+m[9]*m[9])
        scale2 = numpy.sqrt(m[2]*m[2]+m[6]*m[6]+m[10]*m[10])
        scale = numpy.array((scale0,scale1,scale2)).astype('f')
        mat = numpy.reshape(m, (4,4))
        rot = numpy.identity(4).astype('f')
        rot[:3,:3] = mat[:3,:3].astype('f')
        rot[:,0] = (rot[:,0]/scale0).astype('f')
        rot[:,1] = (rot[:,1]/scale1).astype('f')
        rot[:,2] = (rot[:,2]/scale2).astype('f')
        if cleanup:
            rot = glCleanRotMat(rot.ravel())
        rot.shape = (16,)
        #rot1 = rot.astype('f')
        return rot, transl, scale


    def setMatrixComponents(self, rot=None, trans=None, scale=None,
                            redo=1):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Define MatrixRot, MatrixTransl, MatrixScale and MatrixRotInv
from a rotation, translation and scale.
rot should be a 4x4 matrix defining a 3D 3x3 rotation 
trans should be a 3D translation vector
scale should be 3-vector of positive number larger than 0.0
"""
        self._modified = True
        if rot is not None:
            assert rot.shape==(4,4)
            self.MatrixRot = rot.ravel()
            RotInv = numpy.transpose(rot)
            self.MatrixRotInv = numpy.reshape(RotInv, (16,))

        if trans is not None:
            assert len(trans)==3
            self.MatrixTransl = trans

        if scale is not None:
            assert len(scale)==3 and scale[0]>0. and scale[1]>0. and scale[2]>0.
            self.MatrixScale = scale

        if redo:
            self.RedoDisplayList()
            self.viewer.Redraw()


    def SetMatrix(self, matrix):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """ takes a 4x4 matrix.  If shape==(16,), it must be in OpenGL form.
        If shape==(4,4), then it is in standard form, with translation vector
        in right column, etc. calls Decompose4x4 to calculate equivalent
        rotation, translation, and scale matrix/vectors and sets the objects
        attributes, setting the object's original transformation to matrix """
        assert matrix.shape==(16,) or matrix.shape==(4,4)
        self._modified = True
        if matrix.shape==(4,4):
            matrix = numpy.reshape(numpy.transpose(matrix), (16,))
        self.Matrix = matrix
        self.MatrixRot, self.MatrixTransl, self.MatrixScale = self.Decompose4x4(matrix)
        RotInv = numpy.transpose(numpy.reshape(self.MatrixRot, (4,4)))
        self.MatrixRotInv = numpy.reshape(RotInv, (16,))
        # 
        self.MatrixRot = self.MatrixRot.astype('f')
        self.MatrixRotInv = self.MatrixRot.astype('f')
        if self != self.viewer.rootObject:
            self.viewer.deleteOpenglList()
            self.viewer.Redraw()


    def RedirectTransformTo(self, object=None, propagate=1):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Apply transformations to object rather than to myself"""
	
	children = []
	if object:
	    assert isinstance(object, Transformable)

	if object:
	    if hasattr(self, 'AllObjects'):
		children = self.AllObjects()
	    else:
		if object in children:
		    raise AttributeError("%s is a child of %s, therefore it \
already inherits the transformation" % (self.name, object.name) )

	self.propagateRedirection = propagate

	while object and object.redirectXform and object.propagateRedirection:
	    object = object.redirectXform
	self.redirectXform = object
	for o in children:
	    if o.redirectXform: o.redirectXform=object


    def MoveWith(self, object=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Add myself to the list of objects inheriting transf. from object"""
	
        if object:
            assert isinstance(object, Transformable)
        
        if object == self or object is None:
            objs = self.viewer.rootObject.AllObjects()
            for o in objs:
                if self in o.copyXform: o.copyXform.remove(self)
            return

        children = []
        if hasattr(object, 'AllObjects'):
            children = object.AllObjects()

        if self in children:
            raise AttributeError("%s is a child of %s, therefore it already \
inherits the transformation" % (self.name, object.name) )

        if self not in object.copyXform: object.copyXform.append(self)


    def _SetPivot(self, new_pivot):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Set the pivot vector directly"""
        pivot = numpy.array ( new_pivot )
        assert pivot.shape == (3,)
        self.pivot = pivot

        
    def SetPivot(self, new_pivot):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Set the center of rotation.
        the pivot point should be expressed in the coordinate system defined
        by all the transformations above this node including this node's
        transformation"""

        self._modified = True
        MatrixRotInv = numpy.reshape(self.MatrixRotInv, (4,4))
        pivotDiff = self.pivot-new_pivot
        pivotTrans = numpy.dot(MatrixRotInv[:3,:3], pivotDiff)
        scal = self.scale
        scalTrans = -self.MatrixScale*pivotDiff
        newPivotTrans = numpy.dot(MatrixRotInv[:3,:3], scalTrans)
        for j in (0,1,2):
            self.translation[j] = self.translation[j] + \
				  pivotTrans[j]
            for i in (0,1,2):
                self.translation[j] = self.translation[j] + \
                     scal[i] * (newPivotTrans[i]) * self.rotation[i*4+j]
        oldPivot = self.pivot
        self.pivot = numpy.array ( new_pivot )
        return oldPivot
    

# Guillaume wonders if we could write MakeMat this way
# ie get rid of the call with (self.MatrixRotInv) then (self.MatrixRot)
# but I need an exemple where self.MatrixRot is not the identity  
#    def MakeMat(self, scale=True):
#	"""Build the matrix for this object in his parent's frame"""
#        gllib.glTranslatef(float(self.pivot[0]),
#                           float(self.pivot[1]),
#                           float(self.pivot[2]))
#
#        gllib.glTranslatef(float(self.translation[0]),
#                           float(self.translation[1]),
#                           float(self.translation[2]))
#        gllib.glMultMatrixf(self.rotation)
#        if scale:
#            gllib.glScalef(float(self.scale[0]),
#                           float(self.scale[1]),
#                           float(self.scale[2]))
#
#        gllib.glTranslatef(float(self.MatrixTransl[0]),
#                           float(self.MatrixTransl[1]),
#                           float(self.MatrixTransl[2]))
#        gllib.glMultMatrixf(self.MatrixRot)
#        gllib.glScalef(float(self.MatrixScale[0]),
#                       float(self.MatrixScale[1]),
#                       float(self.MatrixScale[2]))
#
#        gllib.glTranslatef(float(-self.pivot[0]),
#                           float(-self.pivot[1]),
#                           float(-self.pivot[2]))


    def MakeMat(self, scale=True):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Build the matrix for this object in his parent's frame"""
        gllib.glTranslatef(float(self.translation[0]),
                           float(self.translation[1]),
                           float(self.translation[2]))
        gllib.glTranslatef(float(self.MatrixTransl[0]),
                           float(self.MatrixTransl[1]),
                           float(self.MatrixTransl[2]))
        #print "self.MatrixRot", self.MatrixRot
        gllib.glMultMatrixf(self.MatrixRot)
        gllib.glTranslatef(float(self.pivot[0]),
                           float(self.pivot[1]),
                           float(self.pivot[2]))
        gllib.glMultMatrixf(self.MatrixRotInv)
        gllib.glMultMatrixf(self.rotation)
        gllib.glMultMatrixf(self.MatrixRot)
        if scale:
            gllib.glScalef(float(self.scale[0]),
                           float(self.scale[1]),
                           float(self.scale[2]))
            gllib.glScalef(float(self.MatrixScale[0]),
                           float(self.MatrixScale[1]),
                           float(self.MatrixScale[2]))
        gllib.glTranslatef(float(-self.pivot[0]),
                           float(-self.pivot[1]),
                           float(-self.pivot[2]))


    def BuildMat(self, obj, root, scale, instance):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Build the matrix by which this object is transformed
instance is a list of integer providing instance indices for all parents
"""

        #print 'build Mat', obj.name, instance
        if obj.parent and obj!=root:
            parentInstance =  instance[:-1]
            self.BuildMat(obj.parent, root, scale, parentInstance)
        obj.MakeMat(scale)
        #print 'build Mat end', obj.name, instance
        #print 'build multiply instance', instance[-1], 'for geom', obj
        instance = obj.instanceMatricesFortran[int(instance[-1])]
        gllib.glMultMatrixf(instance)


## FIXME does not work .. Project return weird stuff
##      def Project(self, point):
##          """Apply all transformations to point and return winx,winy,winz"""
##          camera = self.viewer.currentCamera
##          camera.SetupProjectionMatrix()
##  	GL.glViewport(0, 0, camera.width, camera.height)
##  #        mat = numpy.reshape(numpy.transpose(self.GetMatrix(self)), (16,))
##          mat = numpy.reshape(self.GetMatrix(), (16,))
##          print self.name, mat
##  	GL.glPushMatrix()
##  	GL.glLoadIdentity()
##          GL.glMultMatrix(mat)
##          print 'MOD',GL.glGetDoublev(GL.GL_MODELVIEW_MATRIX)
##          print 'PROJ',GL.glGetDoublev(GL.GL_PROJECTION_MATRIX)
##          projPoint = gluProject(point[0],point[1],point[2])
##    	GL.glPopMatrix()
##          return projPoint


    def GetMatrix(self, root=None, instance=None, scale=True, transpose=True):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Returns the matrix by which this object is transformed
scale = False: returns the rotation and translation. no scaling info included
               Used to save the transformed geom --> coords --> new pdb file
instance is a list of integer instance indices for all parents
"""
        if root is None:
            root = self.viewer.rootObject

        if instance is None:
            instance = [0]
            p = self.parent
            while p:
                instance.append(0)
                p = p.parent

        GL.glPushMatrix()
        GL.glLoadIdentity()
        #print 'GetMatrix', instance
        self.BuildMat(self, root, scale, instance)
        #GL.glMultMatrixf(self.instanceMatricesFortran[instanceList[0]]])
        m = numpy.array(GL.glGetDoublev(GL.GL_MODELVIEW_MATRIX)).astype('f')
        GL.glPopMatrix()
        if numpy.alltrue(m==numpy.zeros(16).astype('f')):
            # this happens when Pmv has no GUI
            m = numpy.identity(4)
        if transpose:
            return numpy.transpose(numpy.reshape(m, (4,4)))
        else:
            return numpy.reshape(m, (4,4))


    def GetMatrixInverse(self, root=None, instance=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Returns the inverse of the matrix used to transform this object"""

        if root is None:
            root = self.viewer.rootObject
        m = self.GetMatrix(root, instance)
        m = numpy.reshape(numpy.transpose(m), (16,)).astype('f')
        rot, transl, scale = self.Decompose4x4(m)
        sc = numpy.concatenate((numpy.reshape(scale,(3,1)), [[1]]))
        n = numpy.reshape(rot, (4,4))/sc
        tr = numpy.dot(n, (transl[0], transl[1], transl[2],1) )
        n[:3,3] = -tr[:3].astype('f')
        return n
##  	m[:3,:3] = numpy.transpose(m[:3,:3])
##  	m[:3,3] = -m[:3,3]
##  	s = numpy.identity(4) * 1.0/root.scale[0]
##  	s[3,3] = 1.0
##  	m = numpy.dot(s, m)


    def FrameTransform(self, camera=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Build the R an RI, the object's frame transformation and inverse"""

	GL.glPushMatrix()
	self.Si = numpy.ones( (3, ) )
	GL.glLoadIdentity()
	if hasattr(self, 'parent'):
            if self.inheritXform:
                parent = self.parent
                while (parent):
                    m = numpy.reshape( parent.rotation, (4,4) )
                    upd = numpy.reshape( numpy.transpose(m), (16, ) )
                    GL.glMultMatrixf(upd)                
                    GL.glMultMatrixf(parent.MatrixRotInv)                
                    self.Si = self.Si / parent.scale
                    self.Si = self.Si / parent.MatrixScale
                    # we have to test here because we need to take into
                    # account the first parent that does not inherit while
                    # building R and Ri
                    if not parent.inheritXform:
                        break
                    parent = parent.parent


	if camera:
	    m = numpy.reshape( camera.rotation, (4,4) )
	    upd = numpy.reshape( numpy.transpose(m), (16, ) )
	    GL.glMultMatrixf(upd)
            self.Si = self.Si / camera.scale

	self.Ri = numpy.array(GL.glGetDoublev(GL.GL_MODELVIEW_MATRIX)).astype('f')
	GL.glPopMatrix()
	self.Ri = glCleanRotMat(self.Ri).astype('f')
        self.Ri.shape = (4,4)
	self.R = numpy.reshape( numpy.transpose(self.Ri), (16, ) ).astype('f')
	self.Ri.shape = (16, )

	if self.redirectXform: self.redirectXform.FrameTransform(camera)
	for o in self.copyXform: o.FrameTransform(camera)


    def SetTransformation(self, matrix, transpose=False, redo=True):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Set the transformation matrix, if transpose is False the rotation
should be in FORTRAN style
"""
        if transpose:
            matrix = numpy.array(numpy.transpose(matrix), 'f')
        if not matrix.flags.contiguous:
            matrix = numpy.array(matrix)
        rot, transl, scale = self.Decompose4x4(matrix.ravel())
        self.SetRotation(rot.astype('f'), redo=False)
        self.SetTranslation(transl, redo=False)
        self.SetScale(scale, redo=False)

        if redo and self != self.viewer.rootObject and \
           not self.immediateRendering:
            self.viewer.deleteOpenglList()
            self.viewer.Redraw()


    def SetRotation(self, matrix, redo=1):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Set the rotation matrix to the object [matrix.shape==(16,)]"""
	assert matrix.shape == (16,)
        self._modified = True
##  	if self.redirectXform: obj = self.redirectXform
##  	else: obj = self
	self.rotation = matrix
##  	for o in self.copyXform: o.SetRotation(camera)
        if redo and self.viewer and self != self.viewer.rootObject and \
           not self.immediateRendering:
            self.viewer.deleteOpenglList()
        #self.viewer.Redraw()


    def ConcatRotationRelative(self, matrix):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Apply the rotation matrix to the object (matrix.shape ==(16,)
        Unlike ConcatRotation you just concatenate the rotation of the object
        without considering Ri and R
        """
        self._modified = True
        obj = self
	GL.glPushMatrix()
	GL.glLoadIdentity()
##         GL.glMultMatrixf(obj.rotation)
        GL.glMultMatrixf(matrix)
        GL.glMultMatrixf(obj.rotation)
        m = numpy.array(GL.glGetDoublev(GL.GL_MODELVIEW_MATRIX)).astype('f')
        obj.rotation = m.astype('f')
        obj.rotation.shape = (16, )
        GL.glPopMatrix()


    def ConcatRotation(self, matrix, redo=1):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Apply the rotation matrix to the object [matrix.shape==(16,)]"""
	
	if self.redirectXform: obj = self.redirectXform
	else: obj = self
        obj._modified = True
	GL.glPushMatrix()
	GL.glLoadIdentity()
	GL.glMultMatrixf(obj.Ri)#.astype('f'))
	GL.glMultMatrixf(matrix)
	GL.glMultMatrixf(obj.R)#.astype('f'))
	GL.glMultMatrixf(obj.rotation)

	m = numpy.array(GL.glGetDoublev(GL.GL_MODELVIEW_MATRIX)).astype('f')
        obj.rotation = glCleanRotMat(m).astype('f')
	obj.rotation.shape = (16, )
	GL.glPopMatrix()
        for o in self.copyXform: o.ConcatRotation(matrix)

        ## This code made rotation very slow because it would rebuild the
        ## master dpyList in cases where it was not needed
##         if redo and not self.immediateRendering:
##             vi = self.viewer
##             print 'faga'
##             vi.deleteOpenglList()
           
        vi = self.viewer
        if vi.activeClippingPlanes > 0 or vi.activeScissor > 0 or \
          (vi.currentObject!=vi.rootObject and not \
           vi.redirectTransformToRoot) and redo and \
           not self.immediateRendering:
            vi.deleteOpenglList()
	

    def SetTranslation(self, trans, redo=1):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Set the translation trans to the object"""
	assert trans.shape == (3,)
        self._modified = True
##  	if self.redirectXform: obj = self.redirectXform
##  	else: obj = self
	self.translation = trans
##  	for o in self.copyXform: o.SetTranslation(camera)
        if self.viewer and self != self.viewer.rootObject and redo and \
           not self.immediateRendering:
            self.viewer.deleteOpenglList()
        #self.viewer.Redraw()
        

    def getCumulatedTranslation(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """calculate the translation cumulated from root
"""
        obj = self
        lCumulatedTranslation = (0.,0.,0.)
        while obj != self.viewer.rootObject:
            lCumulatedTranslation += obj.translation
            obj = obj.parent
        lCumulatedTranslation += obj.translation # to add rootObject translation
        return lCumulatedTranslation
       

    def ConcatTranslation(self, trans, redo=1):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """Apply the translation trans to the object"""
        
        if self.redirectXform: obj = self.redirectXform
        else: obj = self
        obj._modified = True
# MS. Feb 22 '99 obj.R does not contain the scaling factor
#        d = obj.Si * numpy.array( trans )
        #d = numpy.array( trans )
        #d = numpy.concatenate( (d, [1.0]) )  # go to homogenous coords
        d = list(trans)+[1.0]
        rot = numpy.reshape( obj.R, (4,4) )
##         obj.translation = obj.translation + \
##                            numpy.dot( rot, d )[:3]
        obj.translation = obj.translation + self.multMat4pt(rot, d)

        for o in self.copyXform: o.ConcatTranslation(trans)
        vi = self.viewer
        if vi.activeClippingPlanes > 0 or vi.activeScissor > 0 or \
          (vi.currentObject!=vi.rootObject and not \
           vi.redirectTransformToRoot) and redo and \
           not self.immediateRendering:
            vi.deleteOpenglList()
        vi.Redraw()
        event = translateEvent('translate',obj)
        vi.eventHandler.dispatchEvent(event)

    def SetScale(self, scale, redo=1):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Set the scaling factor of the object"""
        #print "SetScale", redo
##  	if self.redirectXform: obj = self.redirectXform
##  	else: obj = self
        self._modified = True
        scale = numpy.array(scale)
	if self.isScalable:
	    assert scale.shape == (3,)
	    if scale[0] > 0.00 and scale[1] > 0.00 and scale[2] > 0.00:
		self.scale = scale
##  	for o in self.copyXform: o.SetScale(camera)
        if self.viewer and (self != self.viewer.rootObject or \
           self.viewer.activeClippingPlanes > 0) and redo and \
           not self.immediateRendering:
            self.viewer.deleteOpenglList()
	#self.viewer.Redraw()
        

    def ConcatScale(self, scale, redo=1):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Multiply the object's scale by scale"""

        if self.redirectXform: obj = self.redirectXform
        else: obj = self
        obj._modified = True
        scale = numpy.array([scale])
        if obj.isScalable:
            if scale[0] > 1.0 or (obj.scale[0] > 0.0 and obj.scale[1] > 0.0
                                      and obj.scale[2] > 0.0):
                obj.scale = obj.scale * scale
        for o in self.copyXform: o.ConcatScale(scale) #was o.ConcatScale(camera)
        vi = self.viewer
        if vi.activeClippingPlanes > 0 or vi.activeScissor > 0 or \
          (vi.currentObject!=vi.rootObject and not \
           vi.redirectTransformToRoot) and redo and \
           not self.immediateRendering:
            vi.deleteOpenglList()


    def transformIsIdentity(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """return true if this object has the identity matrix as
        transformation"""
        mat = self.GetMatrix()
        diff = numpy.sum( (mat-numpy.identity(4)).ravel() )
        if math.fabs(diff) < 0.00001:
            return True
        else:
            return False


    def transformationSourceCode(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """generate the source code for setting the current transformation.
        the returned arryaays can be used to g.SetXXX where XXX is
        Rotation, Translation, Scale or Pivot"""
        mat = self.GetMatrix()
        rot, trans, scale = self.Decompose4x4(numpy.transpose(mat).ravel())
        src = []
        src.append("from numpy import array\n")
        src.append("rot = " + repr(rot) + ".astype('f')\n")
        src.append("trans = " + repr(trans) + "\n")
        src.append("scale = " + repr(scale) + "\n")
        src.append("pivot = " + repr(scale) + "\n")
        return src
