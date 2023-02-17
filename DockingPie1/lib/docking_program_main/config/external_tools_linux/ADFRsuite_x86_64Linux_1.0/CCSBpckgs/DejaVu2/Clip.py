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
# $Header: /mnt/raid/services/cvs/DejaVu2/Clip.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#
# $Id: Clip.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#

from opengltk.OpenGL import GL
from opengltk.extent.utillib import glCleanRotMat
from mglutil.math.rotax import rotVectToVect

import numpy
from viewerFns import getkw
from math import sqrt
from Transformable import Transformable
from colorTool import OneColor, glMaterialWithCheck, resetMaterialMemory
from capGeom import capMesh

class ClippingPlane(Transformable):

    clipPlaneNames = [ ]
    #GL.GL_CLIP_PLANE0, GL.GL_CLIP_PLANE1,
    #GL.GL_CLIP_PLANE2, GL.GL_CLIP_PLANE3,
    #GL.GL_CLIP_PLANE4, GL.GL_CLIP_PLANE5 ]

    clipPlaneColors = [ (.4,1.,1.,1), (1.,.6,1.,1), (1.,1,.5,1),
                        (.4,1,.4,1), (.5,.7,1,1), (1,.5,.8,1) ]

    directions = [(1., 0., 0, 1.), (0., 1., 0, 1.), (0., 0., 1, 1.),
                  (1., 1., 0, 1.), (1., 0., 1, 1.), (0., 1., 1, 1.) ]
                  
    def __init__(self, object, num, viewer):
	"""Create a arbitrary clipping plane.  self.translation represents the
        point in the plane's parent's space about which the plane will
        rotate.  self.eqn is a vector of the  4 coefficients for the equation
        of a plane."""

        self.name = 'ClipPlane'+str(num)
        Transformable.__init__(self, viewer)
	self.num = num
	self.id = self.clipPlaneNames[num]
	self.object = object
        self.Reset()
        self.direction = self.directions[num]
        rot = rotVectToVect(self.direction[:3], (1,.0,0.))
        self.planeRot = numpy.array(rot, 'f').flatten()
        self.ConcatRotation(numpy.identity(4,'f').flatten())
        

    def Reset(self):
        # FIXME since clipping planes are added to objects as they are created
        # this seems superfluous
        # the consequence is that one cannot disable a clipping plane without
        # loosing it! so we should have an add clipping plane button which
        # is separated from enabling the clipping plane
        self.hasBeenCurrent = False

	self.color = self.clipPlaneColors[self.num]
	self.lineWidth = 2
	#self.antiAliased = 0
	self.visible = False
        self.enabled = False
        #self.antialiased = False
	self.eqn = numpy.array( [1.0, 0.0, 0.0, 0.0], numpy.float )
        # this is coefficient vector of the equation of a plane Ax+By+Cz+D = 0
	self.n = 1.0
	self.FrameTransform()
	self.polyMode = GL.GL_LINE_STRIP
        self._modified = False
        

    def ResetTransformation(self, redo=1):
	"""Reset the clipping plane's transformation"""

	Transformable.ResetTransformation(self)
	self.eqn = numpy.array( [1.0, 0.0, 0.0, 0.0] )
	self.n = 1.0
        

    def __repr__(self):
	return '<ClippingPlane %s, eqn=%s>' % (self.name,str(self.eqn) )


    def FrameTransform(self, camera=None):
	"""Build the R an RI, the object's frame transformation and inverse"""

	GL.glPushMatrix()
	self.Si = numpy.ones( (3, ) )
	GL.glLoadIdentity()
	m = numpy.reshape( self.object.rotation, (4,4) )
	upd = numpy.reshape( numpy.transpose(m), (16, ) )
	GL.glMultMatrixf(self.object.Ri)
	GL.glMultMatrixf(upd)
        GL.glMultMatrixf(self.object.MatrixRotInv)
        self.Si = self.Si * self.object.Si / (self.object.scale *
                                              self.object.MatrixScale)

	self.Ri = numpy.array(GL.glGetDoublev(GL.GL_MODELVIEW_MATRIX)).astype('f')
	GL.glPopMatrix()
	#self.Ri = numpy.reshape(glCleanRotMat(self.Ri), (4,4) )
        self.Ri = glCleanRotMat(self.Ri)
	self.R = numpy.reshape( numpy.transpose(self.Ri), (16, ) ).astype('f')
	self.Ri = numpy.reshape(self.Ri, (16, )).astype('f')

	if self.redirectXform: self.redirectXform.FrameTransform(camera)
	for o in self.copyXform: o.FrameTransform(camera)


    def _NormalizeN(self):
	eqn = numpy.zeros( (3,), numpy.float )
	for i in (0,1,2):
	    eqn[i] = self.rotation[i]
	n = numpy.add.reduce(eqn*eqn)
	assert n > 0.0
	if n > 0.00001:
	    self.n = 1.0 / sqrt( n )
	    self.eqn[:3] = eqn
	else: 
	    self.n = 1.0
	    self.eqn[:3] = [1.0, 0.0, 0.0]


    def ConcatRotation(self, matrix):
	"""Overwrite rotation methods"""

        self._modified = True
	Transformable.ConcatRotation(self, matrix)
	self.rotation.shape = (4,4)
        #eqn = numpy.array(self.rotation[0,:3]) # because direction is (1.,0.,0.)
        eqn = numpy.dot(self.direction, self.rotation)
        #print '================================================='
        #print 'eqn', eqn, eqn.shape

        #self.n = 1.0 / sqrt( numpy.add.reduce(eqn*eqn) )
        x,y,z = eqn[:3]
        self.n = 1.0 / sqrt( x*x + y*y +z*z)
        self.eqn[:3] = eqn[:3]
        #print 'self.eqn',self.eqn
        self.eqn[3] = -numpy.dot(self.eqn[:3], self.translation)
        #print self.eqn

        #get the value so that the plane is equivalent to having the proper
        #normal vector and a translation from the origin along the vector
	for o in self.copyXform: o._NormalizeN()
	self.rotation.shape = (16, )
        self.viewer.deleteOpenglList()


    def ConcatTranslation(self, trans):
	"""
        Overwrite translation methods
        """
        #print "ConcatTranslation", trans
        self._modified = True
        self.rotation.shape = (4,4)
        lTransInLocalSystem = numpy.dot(
                                         (trans[0], trans[1], trans[2], 1.),
                                         self.rotation,
                                          )[:3]
        #print "ConcatTranslation", lTransInLocalSystem
        self.rotation.shape = (16, )
        self.translation = self.translation + lTransInLocalSystem
        self.eqn[3] = -numpy.dot(self.eqn[:3], self.translation)

        for o in self.copyXform:
            o.ConcatTranslation(lTransInLocalSystem)
        self.viewer.deleteOpenglList()


    def screenConcatTranslation(self, trans):
        """
        Overwrite translation methods
        """
        #print "localConcatTranslation", trans
        self._modified = True
        self.viewer.rootObject.rotation.shape = (4,4)
        self.viewer.rootObject.rotation.transpose()
        lTransInScreenSystem = numpy.dot(
                                         self.viewer.rootObject.rotation,
                                         (trans[0], trans[1], trans[2], 1.),
                                          )[:3]
        self.viewer.rootObject.rotation.transpose()
        self.viewer.rootObject.rotation.shape = (16, )
        self.translation = self.translation + lTransInScreenSystem
        self.eqn[3] = -numpy.dot(self.eqn[:3], self.translation)

        for o in self.copyXform:
            o.ConcatTranslation(lTransInScreenSystem)
        self.viewer.deleteOpenglList()


    def _Enable(self, side):
	"""Activate the clipping plane"""

        eqnt = self.eqn * side
            #eqnt[3] = 0.0 - eqnt[3]
        GL.glClipPlane(self.clipPlaneNames[self.num], eqnt)
        GL.glEnable(self.clipPlaneNames[self.num])


    def _Disable(self):
	"""Deactivate the clipping plane"""
	GL.glDisable(self.clipPlaneNames[self.num])


    def Set(self, **kw):
	"""Set various clipping plane parameters"""

	self.hasBeenCurrent = True # remember the light has been changed

        tagModified = True
        val = getkw(kw, 'tagModified')
        if val is not None:
            tagModified = val
        assert tagModified in [True, False]
        self._modified = tagModified

	val = getkw(kw, 'enabled')
	if not val is None:
	    if val is True:
                self._Enable(1)
	    elif val is False:
                self._Disable()
	    else: raise AttributeError('enable can only be True or False')
            self.enabled = val

	val = getkw(kw, 'name')
	if not val is None: self.name = val

	val = getkw(kw, 'visible')
	if not val is None:
	    if val in [False,True]:
                self.visible = val
	    else:
                raise AttributeError('visible can only be 0 or 1')

	val = getkw(kw, 'color')
	if not val is None:
	    col = OneColor( val )
	    if col:
		self.color = col

	val = getkw(kw, 'lineWidth')
	if not val is None:
	    try:
	        int(val)
	    except:
	        raise ValueError ('lineWidth must be a positive int')
	    if val>=1:
	        self.lineWidth = int(val)
	    else:
	        raise ValueError ('lineWidth must be a positive int')
            

#	val = getkw(kw, 'antialiased')
#	if not val is None:
#	    if val in (True, False) :
#		self.antialiased = val
#	    else: raise ValueError ('antialiased can only be YES or NO')

	val = getkw(kw, 'rotation')
	if not val is None:
            self.rotation = numpy.identity(4, 'f').ravel()
            mat = numpy.reshape(numpy.array(val), (16,)).astype('f')
            self.ConcatRotation(mat)

	val = getkw(kw, 'translation')
	if not val is None:
            self.translation = numpy.zeros( (3,), 'f')
            mat = numpy.reshape(numpy.array(val), (3,)).astype('f')
            self.ConcatTranslation(mat)

	val = getkw(kw, 'scale')
	if not val is None:
            self.SetScale( val )

	val = getkw(kw, 'pivot')
	if not val is None:
            self.SetPivot( val )

	if len(kw):
	    print 'WARNING1: Keyword(s) %s not used' % kw.keys()

        if self.object.viewer:
            self.object.viewer.objectsNeedingRedo[self.object] = None
            self.object.viewer.Redraw()


    def setColor(self, val):
        self.Set(color=val)


    def getState(self):
        """return a dictionary describing this object's state
This dictionary can be passed to the Set method to restore the object's state
"""
        return {'enabled':self.enabled,
                'name':self.name,
                'visible':self.visible,
                'color':self.color,
                'lineWidth':self.lineWidth,
                #'antialiased':self.antialiased,
                'rotation':list(self.rotation),
                'translation':list(self.translation),
                'scale':list(self.scale),
                'pivot':list(self.pivot)
                }


    def DisplayFunction(self):
        """Draw a square with diagonals to represent the clipping plane
"""
        #print "ClippingPlane.DisplayFunction"
	#trans = self.eqn[3]*(self.eqn[:3]*self.n)
        resetMaterialMemory()
        trans = self.translation
	GL.glPushMatrix()
	#GL.glTranslatef(-trans[0],-trans[1],-trans[2])
        GL.glTranslatef(float(trans[0]),
                        float(trans[1]),
                        float(trans[2]))
	GL.glMultMatrixf(self.rotation)
	GL.glScalef(float(self.scale[0]),
                    float(self.scale[1]),
                    float(self.scale[2]))

	if self.polyMode == GL.GL_QUADS:

	    GL.glPushAttrib(GL.GL_CURRENT_BIT | GL.GL_LIGHTING_BIT |
			    GL.GL_POLYGON_BIT)
            GL.glDisable(GL.GL_LIGHTING)

	    GL.glMaterialWithCheck(GL.GL_FRONT_AND_BACK, GL.GL_AMBIENT,
                                self.color)
            if self.viewer is not None:
                self.viewer.enableOpenglLighting()
	    GL.glPolygonMode(GL.GL_FRONT_AND_BACK, GL.GL_FILL)

            GL.glPushMatrix()
            GL.glMultMatrixf(self.planeRot)

	    GL.glBegin (GL.GL_QUADS)
	    GL.glVertex3f (0.0, -5.0, -5.0)
	    GL.glVertex3f (0.0, -5.0,  5.0)
	    GL.glVertex3f (0.0,  5.0,  5.0)
	    GL.glVertex3f (0.0,  5.0, -5.0)
	    GL.glVertex3f (0.0, -5.0, -5.0)
	    GL.glEnd ()
            GL.glPopMatrix()

	    GL.glPopAttrib()

	else:
# MS disabling GL.GL_BLEND breaks display of transparent surfaces
##  	    if self.antialiased==True:
##  		GL.glEnable(GL.GL_LINE_SMOOTH)
##  		GL.glEnable(GL.GL_BLEND)
##  	    else:
##  		GL.glDisable(GL.GL_LINE_SMOOTH)
##  		GL.glDisable(GL.GL_BLEND)

	    GL.glColor4fv (self.color)
	    GL.glLineWidth(self.lineWidth)

            GL.glPushMatrix()
            GL.glMultMatrixf(self.planeRot)
            
	    # could and should be a display list made once for all planes
	    GL.glBegin (GL.GL_LINE_STRIP)
	    GL.glVertex3f (0.0, -5.0, -5.0)
	    GL.glVertex3f (0.0, -5.0,  5.0)
	    GL.glVertex3f (0.0,  5.0,  5.0)
	    GL.glVertex3f (0.0,  5.0, -5.0)
	    GL.glVertex3f (0.0, -5.0, -5.0)
	    GL.glVertex3f (0.0,  5.0,  5.0)
	    GL.glVertex3f (0.0, -5.0,  5.0)
	    GL.glVertex3f (0.0,  5.0, -5.0)
	    GL.glEnd ()
            GL.glPopMatrix()

	GL.glPopMatrix()


    def getCapMesh(self, geom, delta=0.01):
        x,y,z,s = self.eqn
        planeNormal = (x,y,z)
        delta *= -geom.clipSide[self.num]
        planePoint = (-x*(s+delta), -y*(s+delta), -z*(s+delta))
        edges, faceEdges = geom.buildEdgeList()
        vertsC, facesC = capMesh(
            geom.getVertices(), edges, geom.getFaces(), faceEdges,
            planePoint, planeNormal, 0.01)

        vc = []
        fc = []
        for v, f in zip(vertsC, facesC):
            l = len(vc)
            vc.extend(v)
            fc.extend( [ (t[0]+l, t[1]+l, t[2]+l) for t in f] )
        return vc, fc

    
    def ClipCapMesh(self, obj, onOff):
        if not hasattr(obj, 'buildEdgeList'):
            return 
        vi = obj.viewer
        if onOff:
            x,y,z,s = self.eqn
            vc, fc = self.getCapMesh(obj)
            from DejaVu2.IndexedPolygons import IndexedPolygons
            capg = IndexedPolygons(
                'cap%d'%self.num, vertices=vc, faces=fc, 
                fnormals = ( (-x,-y,-z),),
                inheritFrontPolyMode=False, frontPolyMode='fill',
                inheritBackPolyMode=False, backPolyMode='fill',
                inheritCulling=0, culling='none',
                inheritShading=0, shading='flat'
                )
            vi.AddObject(capg, parent=obj)

            obj.capMesh = 1
            vi.cappedGeoms.append( (obj, self, capg) )
        else:
            obj.capMesh = 0
            for g,c,gc in obj.viewer.cappedGeoms:
                #print 'trying to remove', g, self
                if g==obj and c==self:
                    #print 'removing'
                    vi.cappedGeoms.remove( (g,c,gc) )
                    vi.RemoveObject(gc)
                    break

        if obj != vi.rootObject:
            vi.objectsNeedingRedo[obj] = None
        vi.Redraw()
