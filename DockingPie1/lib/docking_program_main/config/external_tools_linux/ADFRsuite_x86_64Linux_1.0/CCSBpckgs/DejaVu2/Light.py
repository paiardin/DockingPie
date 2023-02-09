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

# $Header: /mnt/raid/services/cvs/DejaVu2/Light.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#
# $Id: Light.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#

from opengltk.OpenGL import GL,GLU
from opengltk.extent.utillib import glCleanRotMat

import types
import numpy
from math import sqrt
from colorTool import OneColor
from viewerFns import GetArray, getkw
from Transformable import Transformable


class LightModel:
    """Class for the OpenGL light model"""

    def Set(self, **kw):
	"""Set various light model parameters"""

        self.viewer.currentCamera.Activate()

        tagModified = True
        val = getkw(kw, 'tagModified')
        if val is not None:
            tagModified = val
        assert tagModified in [True, False]
        self._modified = tagModified

	ambi = getkw(kw, 'ambient')
	if not ambi is None:
	    if len(ambi)==3 or len(ambi)==4: 
	        self.ambient = OneColor( ambi )
	        GL.glLightModelfv(GL.GL_LIGHT_MODEL_AMBIENT, self.ambient);
	    else: 
	        raise ValueError('length of new color must be 3 or 4') 

	localViewer = getkw(kw, 'localViewer')
	if not localViewer is None:
	    if localViewer in (True,1):
		GL.glLightModelf(GL.GL_LIGHT_MODEL_LOCAL_VIEWER,
			       GL.GL_TRUE);
	    elif localViewer in (False,0):
		GL.glLightModelf(GL.GL_LIGHT_MODEL_LOCAL_VIEWER,
			       GL.GL_FALSE);
	    else: raise AttributeError('localViewer can only be True or False')
	    self.localViewer = localViewer

	twoSide = getkw(kw, 'twoSide')
	if not twoSide is None:
	    if twoSide in (True,1):
		GL.glLightModelf(GL.GL_LIGHT_MODEL_TWO_SIDE,
			       GL.GL_TRUE);
	    elif twoSide in (False,0):
		GL.glLightModelf(GL.GL_LIGHT_MODEL_TWO_SIDE,
			       GL.GL_FALSE);
	    else: raise AttributeError('twoSide can only be True or False')
	    self.twoSide = twoSide

        self.broadcast()
        
	if len(kw):
	    print 'WARNING8: Keyword(s) %s not used' % kw.keys()


    def broadcast(self):
        #print "LightModel.broadcast"
        if self.viewer.rootObject is None: return
#        self.viewer.needsRedraw = 1
        for app in self.applyTo:
            if app and app.winfo_ismapped():  # needed otherwise tests produce seg fault on mesa 
                app.tk.call(app._w, 'makecurrent')
                self.apply()
                app.tkRedraw()
        self.viewer.Redraw()


    def apply(self):
        """setup current lightmodel for current OpenGL context"""

	if self.ambient is not None:
	    GL.glLightModelfv(GL.GL_LIGHT_MODEL_AMBIENT, self.ambient);

        if self.localViewer is True:
            GL.glLightModelf(GL.GL_LIGHT_MODEL_LOCAL_VIEWER, GL.GL_TRUE);
        elif self.localViewer is False:
            GL.glLightModelf(GL.GL_LIGHT_MODEL_LOCAL_VIEWER, GL.GL_FALSE);

        if self.twoSide is True:
            GL.glLightModelf(GL.GL_LIGHT_MODEL_TWO_SIDE, GL.GL_TRUE);
        elif self.twoSide is False:
            GL.glLightModelf(GL.GL_LIGHT_MODEL_TWO_SIDE, GL.GL_FALSE);


    def getState(self):
        """return a dictionary describing this object's state
This dictionary can be passed to the Set method to restore the object's state
"""
        return {'ambient':self.ambient,
                'localViewer':self.localViewer,
                'twoSide':self.twoSide}

        
    def Reset(self):
	"""Restore the default values"""

        # light model state
        self.localViewer = False
        self.twoSide = False

#        from DejaVu2 import preventIntelBug_WhiteTriangles
#        if preventIntelBug_WhiteTriangles:
#            self.ambient = (.4, .4, .4, 1.)
#        else:
#            self.ambient = (.2, .2, .2, 1.)
#        self.Set( ambient=self.ambient, localViewer=self.localViewer, twoSide=self.twoSide )

        self.ambient = (.12, .12, .12, 1.)
        self.Set( localViewer=self.localViewer, twoSide=self.twoSide )
        self._modified = False


    def __init__(self, viewer, **kw):
        self.applyTo = []
        self.viewer = viewer
        self.Reset()  # creates all state attributes
        kw['tagModified'] = False
        apply( LightModel.Set, (self,), kw )
        

    def __repr__(self):
	return '<LightModel localViewer=%d twoSide=%d color=%s)' % \
	     (self.localViewer, self.twoSide, str(self.ambient))



class Light(Transformable):
    """Class for OpenGL light sources"""

    lightIndices = ( GL.GL_LIGHT0, GL.GL_LIGHT1,
		     GL.GL_LIGHT2, GL.GL_LIGHT3,
		     GL.GL_LIGHT4, GL.GL_LIGHT5, 
		     GL.GL_LIGHT6, GL.GL_LIGHT7 )

    def Set(self, **kw):
        """ set light values.  For direction, position, and spot direction,
vals are given in absolute coordinates (independent of camera or
object).  For these three values, the flag is set to 1 when they are
changed.
"""
        #print "Light.Set"
        self.hasBeenCurrent = True # remember the light has been changed

        tagModified = True
        val = getkw(kw, 'tagModified')
        if val is not None:
            tagModified = val
        assert tagModified in [True, False]
        self._modified = tagModified

        self.viewer.currentCamera.Activate()
        
	val = getkw(kw, 'ambient')
	if not val is None:
	    #self.ambient = OneColor( val )
	    #GL.glLightfv(self.num, GL.GL_AMBIENT, self.ambient )
	    if len(val)==3 or len(val)==4: 
	        self.ambient = OneColor( val )
                GL.glLightModelfv(GL.GL_LIGHT_MODEL_AMBIENT, self.ambient)
                GL.glLightfv(self.num, GL.GL_AMBIENT, self.ambient) # needed for mesa 
	    else: 
	        raise ValueError('length of new color must be 3 or 4') 

	val = getkw(kw, 'diffuse')
	if not val is None:
	    self.diffuse = OneColor( val )
	    GL.glLightfv(self.num, GL.GL_DIFFUSE, self.diffuse )

	val = getkw(kw, 'specular')
	if not val is None:
	    self.specular = OneColor( val )
	    GL.glLightfv(self.num, GL.GL_SPECULAR, self.specular )

	val = getkw(kw, 'direction')
	if not val is None:
            val = list(val)
            if len(val)==3: val += [0.]
            assert len(val)==4
	    self.direction = val
	    self.direction[3] = 0.0
            self.dirFlag = 1  # tell the camera to redraw this light
            self.positional = False

	val = getkw(kw, 'position')
	if not val is None:
            val = list(val)
            if len(val)==3: val += [1.]
            assert len(val)==4
	    self.position = val
	    self.position[3] = 1.0
            self.posFlag = 1  # tell the camera to redraw this light
            self.positional = True

	val = getkw(kw, 'spotDirection')
	if not val is None:
            val = list(val)
            if len(val)==3: val += [0.]
            assert len(val)==4
	    self.spotDirection = val
	    self.spotDirection[3] = 0.0
            self.spotFlag = 1  # tell the camera to redraw this light

	val = getkw(kw, 'spotExponent')
	if not val is None:
	    self.spotExponent = float(val)
	    GL.glLightfv(self.num, GL.GL_SPOT_EXPONENT, [self.spotExponent])
            
	val = getkw(kw, 'spotCutoff')
	if not val is None:
            if val > 180.:
                raise ValueError("spotCutoff must be in [0., 90.] or 180.")
	    self.spotCutoff = float( val )
	    GL.glLightfv(self.num, GL.GL_SPOT_CUTOFF, [self.spotCutoff] )

	val = getkw(kw, 'constantAttenuation')
	if not val is None:
	    self.constantAttenuation = float( val )
            if self.constantAttenuation < 0.0:
                raise ValueError("constantAttenuation must be >= 0.0")
	    GL.glLightfv(self.num, GL.GL_CONSTANT_ATTENUATION,
                         [self.constantAttenuation] )

	val = getkw(kw, 'linearAttenuation')
	if not val is None:
	    self.linearAttenuation = float( val )
            if self.linearAttenuation < 0.0:
                raise ValueError("linearAttenuation must be >= 0.0")
	    GL.glLightfv(self.num, GL.GL_LINEAR_ATTENUATION,
                         [self.linearAttenuation] )

	val = getkw(kw, 'quadraticAttenuation')
	if not val is None:
	    self.quadraticAttenuation = float( val )
            if self.quadraticAttenuation < 0.0:
                raise ValueError("quadraticAttenuation must be >= 0.0")
	    GL.glLightfv(self.num, GL.GL_QUADRATIC_ATTENUATION,
                         [self.quadraticAttenuation] )

	val = getkw(kw, 'positional')
	if not val is None:
	    if val is True: self.position[3] = 1.0
	    elif val is False: self.position[3] = 0.0
	    else: raise AttributeError('positional can only be True or False')
	    self.positional = val

	val = getkw(kw, 'enabled')
	if not val is None:
	    if val in (True, 1): GL.glEnable(self.num)
	    elif val in (False, 0): GL.glDisable(self.num)
	    else: raise AttributeError('enabled can only be True or False')
            self.enabled = val
            
	val = getkw(kw, 'visible')
	if not val is None:
	    if val in (True, False): self.visible = val
	    else: raise AttributeError('visible can only be True or False')

	val = getkw(kw, 'lineWidth')
	if not val is None:
	    if val >= 1:
                self.lineWidth = int(val)
	    else: raise AttributeError('lineWidth has to be >= 1')

	val = getkw(kw, 'length')
	if not val is None:
	    if val > 0.0: self.length = float ( val )
	    else: raise AttributeError('length has to be > 0.0')

#	val = getkw(kw, 'antialiased')
#	if not val is None:
#	    if val in (True, False):
#		self.antialiased = val
#	    else: raise ValueError ('antialiased can only be True or False')

	val = getkw(kw, 'rotation')
	if not val is None:
            mat = numpy.reshape(numpy.array(val), (16,)).astype('f')
            self.rotation = mat

	val = getkw(kw, 'translation')
	if not val is None:
            mat = numpy.reshape(numpy.array(val), (3,)).astype('f')
            self.translation = mat

	val = getkw(kw, 'scale')
	if not val is None:
            mat = numpy.reshape(numpy.array(val), (3,)).astype('f')
            self.scale = mat

	val = getkw(kw, 'pivot')
	if not val is None:
            mat = numpy.reshape(numpy.array(val), (3,)).astype('f')
            self.pivot = mat

	if len(kw):
	    print 'WARNING9: Keyword(s) %s not used' % kw.keys()

        #guillaume vareille 9/29/2005 : 
        # was c = self.viewer.cameras[0]
        # was generating alternativly good and wrong rendering when 2 cameras 
        c = self.viewer.currentCamera 
        
        # force light to be update in viewer
        c.Redraw()

        # brodcast to other application that want to know about that light
        # using aftere does not seem to make it better
        #c.after_idle(self.broadcast)
        self.broadcast()
        
        # not needed and causes the master dispaly list to be rebuilt
        #self.viewer.deleteOpenglList()


    def setAmbient(self, val):
        self.Set(ambient=val)


    def setDiffuse(self, val):
        self.Set(diffuse=val)


    def setSpecular(self, val):
        self.Set(specular=val)


    def broadcast(self):
        #print "Light.broadcast"
        for app in self.applyTo:
            # provides a small performance increase but the light is not
            # setup the first time the window comes up
            if app and app.winfo_ismapped():  # needed otherwise tests produce seg fault on mesa 
                app.tk.call(app._w, 'makecurrent')
                self.apply()
                app.tkRedraw()
        
        if self.viewer.materialEditor is not None:
            self.viewer.materialEditor.tk.call(self.viewer.materialEditor._w, 'makecurrent')
            self.apply()
            self.viewer.materialEditor.tkRedraw()
            

    def apply(self):
        """setup this light for current OpenGL context
"""
        #print "Light.apply"

        num = self.num

        if self.ambient is not None:
            GL.glLightfv(num, GL.GL_AMBIENT, self.ambient )

        if self.diffuse is not None:
            GL.glLightfv(num, GL.GL_DIFFUSE, self.diffuse )

        if self.specular is not None:
            GL.glLightfv(num, GL.GL_SPECULAR, self.specular )

        if self.positional is False:
            GL.glLightfv(num, GL.GL_POSITION, self.direction )
        else:
            GL.glLightfv(num, GL.GL_POSITION, self.position )

        if self.spotFlag:
            GL.glLightfv(num, GL.GL_SPOT_DIRECTION,
                         self.spotDirection[:3] )
            GL.glLightfv(num, GL.GL_SPOT_EXPONENT, [self.spotExponent] )
            GL.glLightfv(num, GL.GL_SPOT_CUTOFF, [self.spotCutoff] )
            GL.glLightfv(num, GL.GL_CONSTANT_ATTENUATION,
                         [self.constantAttenuation] )
            GL.glLightfv(num, GL.GL_LINEAR_ATTENUATION,
                         [self.linearAttenuation] )
            GL.glLightfv(num, GL.GL_QUADRATIC_ATTENUATION,
                         [ self.quadraticAttenuation] )
        
        if self.enabled is True:
            GL.glEnable(num)
        else:
            GL.glDisable(num)


    def getState(self):
        """return a dictionary describing this object's state
This dictionary can be passed to the Set method to restore the object's state
"""
        return { 'ambient':self.ambient,
                 'diffuse':self.diffuse,
                 'specular':self.specular,
                 'direction':self.direction,
                 'position':self.position,
		 'spotDirection':self.spotDirection,
		 'spotExponent':self.spotExponent,
		 'spotCutoff':self.spotCutoff,
		 'constantAttenuation':self.constantAttenuation,
		 'linearAttenuation':self.linearAttenuation,
		 'quadraticAttenuation':self.quadraticAttenuation,
		 'positional':self.positional,
		 'enabled':self.enabled,
		 #'antialiased':self.antialiased,
		 'lineWidth':self.lineWidth,
		 'length':self.length,
		 'visible':self.visible,
                 'rotation':list(self.rotation),
                 'translation':list(self.translation),
                 'scale':list(self.scale),
                 'pivot':list(self.pivot) }
        

    def Reset(self):
	"""Restore Light default values"""

	Transformable.__init__(self, self.viewer)
#	self.translation = 10.0
	if self.num == GL.GL_LIGHT0: enabled = True
	else: enabled = False

	self.Set( ambient = [.12, .12, .12, 1.],
		  diffuse = [1., 1., 1., 1.],
		  specular = [.8, .8, .65, 1.],
		  direction = [-1., 1., 1., 1.],
		  position =  [0.0, 0.0, 0.0, 0.0],
		  spotDirection = [0.0, 0.0, -1.0, 0.0],
		  spotExponent = 0.0,
		  spotCutoff = 180.0,
		  constantAttenuation = 1.0,
		  linearAttenuation = 0.0,
		  quadraticAttenuation = 0.0,
		  positional = False,
		  enabled = enabled,
		  #antialiased = False,
		  lineWidth = 2,
		  length = 10.0,
		  visible = False)
	self.hasBeenCurrent = False
        self._modified = False
        

    def __init__(self, num, viewer, **kw):

        self.hasBeenCurrent = False # turn to True if light is modified

        self.applyTo = [viewer.materialEditor]

	self.num = self.lightIndices[num]
    
        self.viewer = viewer
        self.name = 'light '+str(num)
	self.Reset()
	apply( Light.Set, (self,), kw)
        self._modified = False
	self.object = None
        self.dirFlag = 0
        self.posFlag = 0
        self.spotFlag = 0


    def FrameTransform(self, camera=None):
	"""Build the R an RI, the object's frame transformation and inverse"""

        Transformable.FrameTransform(self, camera)
	if self.redirectXform: self.redirectXform.FrameTransform(camera)
	for o in self.copyXform: o.FrameTransform(camera)

	if self.positional is False: return


    def ConcatRotationDir(self, matrix):
	"""Rotate a directional light"""

        self._modified = True
##  	matrix = numpy.reshape( matrix, (4,4) )
##  	newdir = numpy.dot( self.direction, matrix )
	GL.glPushMatrix()
	GL.glLoadIdentity()
	GL.glMultMatrixf(self.Ri)
	GL.glMultMatrixf(matrix)
	GL.glMultMatrixf(self.R)
        mat = numpy.array(GL.glGetDoublev(GL.GL_MODELVIEW_MATRIX)).astype('f')
        #mat = glCleanRotMat(mat).astype('f')
        n = glCleanRotMat(mat).astype('f')
	GL.glPopMatrix()
        #n = numpy.reshape(mat, (4,4))
  	newdir = numpy.dot( self.direction, n )
	self.Set(direction = newdir)
	for o in self.copyXform: o.ConcatRotationDir(matrix)


    def DrawDirectionalLight(self, camera):
	"""Draw a directional light as a line"""

#	if self.antialiased is True:
#	    GL.glEnable(GL.GL_LINE_SMOOTH)
#	    GL.glEnable(GL.GL_BLEND)
#	else:
#	    GL.glDisable(GL.GL_LINE_SMOOTH)
#	    GL.glDisable(GL.GL_BLEND)

	GL.glColor4fv (self.diffuse)
	GL.glLineWidth(self.lineWidth)
	GL.glBegin (GL.GL_LINES)
	GL.glVertex3f ( float(camera.lookAt[0]),float(camera.lookAt[1]),float(camera.lookAt[2]) )
	v = camera.lookAt + self.direction[:3]
	nv = sqrt(numpy.add.reduce(v*v))
	v = (v/nv)*self.length
	GL.glVertex3f ( float(v[0]), float(v[1]), float(v[2]) )
	GL.glEnd()

	GL.glPushMatrix()
	GL.glTranslatef(float(v[0]),float(v[1]),float(v[2]))
	GL.glPushAttrib(GL.GL_ALPHA_TEST);
#	GL.glDisable(GL.GL_LIGHTING);
	GL.glEnable(GL.GL_BLEND)
	GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);  
	pObj = GLU.gluNewQuadric()
	GLU.gluQuadricDrawStyle(pObj, GLU.GLU_FILL)
	GLU.gluQuadricNormals(pObj, GLU.GLU_SMOOTH)
        
	GL.glColor4ub(255,255,255,10);
	GLU.gluDisk(pObj, 0.0, 0.5, 20, 20);
        
	for light_color in range(100):               
		GL.glColor4ub(255-light_color,255-light_color,255-light_color,10)  
		GLU.gluCylinder(pObj, 0.5 - 2.*light_color/255.,0.5 - 2.*(light_color+1)/255. , light_color/255., 20, 20)
		GL.glTranslatef(0.0, 0.0, float(light_color)/255.);
	#GL.glColor4ub(255,255,255,200)  
	#GLU.gluDisk(pObj, 0,0.5 - 2.*(light_color+1)/255. , 20, 20)
	GL.glPopAttrib()
	GL.glPopMatrix()
