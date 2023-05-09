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
# $Header: /mnt/raid/services/cvs/DejaVu2/Materials.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#
# $Id: Materials.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#

import types
import numpy
from opengltk.OpenGL import GL

from DejaVu2 import materialsDef
from DejaVu2 import viewerConst
from DejaVu2.viewerFns import checkKeywords, getkw

propertyNum = {
    'ambient':0, 'AMBIENT':0, 'ambi':0, 'AMBI':0,
    'diffuse':1, 'DIFFUSE':1, 'diff':1, 'DIFF':1,
    'emission':2, 'EMISSION':2, 'emis':2, 'EMIS':2,
    'specular':3, 'SPECULAR':3, 'spec':3, 'SPEC':3, 
    'shininess':4, 'SHININESS':4, 'shini':4, 'SHINI':4,
    'opacity':5, 'OPACITY':5, 'opac':5, 'OPAC':5,
    }



class Materials:
    """Class for material and color properties of an object
We have a numpy array of RGBA values for the ambient, diffuse,
emissive, andspecular component (in this order) and a list of values for
shininess and opcaity. The list of opcaities will be placed as the 4th
column in the diffuse component at rendering time.
  The constructor can be called either with no argument to get a default
material, or with a property Name specifying a material table and an
index into this table (i.e. tropical, 0)."""

    keywords = [
        'ambient',
        'diffuse',
        'emission',
        'specular',
        'shininess',
        'opacity',
        'binding',
        ]

    ambient = AMBIENT = ambi = AMBI = 0
    diffuse = DIFFUSE = diff = DIFF = 1
    emission = EMISSION = emis = EMIS = 2
    specular = SPECUILAR = spec = SPEC = 3
    shininess = SHININESS = shini = SHINI = 4
    opacity = OPACITY = opac = OPAC = 5

    
    def InitDefault(self):
	self.prop[self.ambi] = numpy.array( (( 0.1, 0.1, 0.1, 1.0 ), ), 'f' )
	self.prop[self.diff] = numpy.array( (( 1.0, 1.0, 1.0, 1.0 ), ), 'f'  )
	self.prop[self.emis] = numpy.array( (( 0.0, 0.0, 0.0, 1.0 ), ), 'f'  )
	self.prop[self.spec] = numpy.array( (( 0.4, 0.4, 0.4, 1.0 ), ), 'f'  )
	self.prop[self.shini] = numpy.array( ( 75.0, ), 'f' )
	self.prop[self.opac] = numpy.array( ( 1.0, ), 'f' )
        self._modified = False
        

    def __init__(self, propName=None, propInd=None):
        """prop can be a string fro the properties library"""
        if propName is not None and propInd is not None:
            prop = getattr(materialsDef, propName)[propInd]
            self.prop = [0,0,0,0,0,0]
            for i in (0,1,2,3,4):
                self.prop[i] = numpy.array( prop[i], 'f' )
            self.prop[5] = numpy.array( self.prop[1][:,3] )

        else:
            self.prop = range(6)
            self.InitDefault()
        self.binding = [ viewerConst.OVERALL,]*6

        # any of the AMBIENT, DIFFUSE, EMISSION, and SPECULAR property
        # can be obtained from another one of these properties by multiplying
        # the source property by a scaling factor going from ]0,1].
        # self.getFrom is used by the GetProperty(). This list contains
        # for each property either None or a 2-tuple where the first entry
        # is the source property (string or int) and the second entry is a
        # scaling factor.
        # If the value is None, the requested property is obtained from its
        # own array of values (.prop[property])
        # Else, it is obtained by multiplying the source propertys array
        # of values of the by the scaling factor.
        # for example, by default the ambient component is obtained by
        # multiplying the diffuse component by 0.2 the diffuse component.

        from DejaVu2 import preventIntelBug_WhiteTriangles
        if preventIntelBug_WhiteTriangles:
            self.getFrom = [ None, None, None, None,]
        else:
            # This line was causing diffuse to always be 60% of diffuse hence
            # preventing white color
            #self.getFrom = [ [self.diffuse, .6], [self.diffuse, .6], None, None,]
            self.getFrom = [ [self.diffuse, .6], None, None, None,]


    def GetProperty(self, prop):
        """return a given property and its binding while respecting the
self.getFrom attribute. This is called by RedoDisplayList.
prop can be a string, an integer, or an RGB tuple
"""
        if type(prop) is types.StringType:
            prop = propertyNum[prop]
        assert prop < 4
        if self.getFrom[prop] is None:
            return self.binding[prop], self.prop[prop]
        else:
            sourceProp, sca = self.getFrom[prop]
            if sourceProp is not None:
                if type(sourceProp) is types.StringType:
                    origProp = self.prop[propertyNum[sourceProp]]
                else:
                    origProp = self.prop[sourceProp]
                    
            if type(sourceProp) is types.StringType:
                sourceProp = propertyNum[sourceProp]
                return self.binding[sourceProp], p

            elif type(sourceProp) is types.IntType:
                p = (self.prop[sourceProp]*sca).astype('f')
                if sourceProp==1:
                    p[:,3] = origProp[:,3] # restore opacity
                try:
                    if len(p)==len(self.prop[sourceProp]):
                        return self.binding[sourceProp], p
                    else:
                        return 11, p
                except TypeError:
                    return 10, p

            elif len(sourceProp) == 4: # single RGBA tuple
                try: # duplicate RGB to match length of sca vector
                    sca = numpy.array(sca)
                    sca.shape = (-1,1)
                    p = numpy.array( ((sourceProp),)*len(sca) )*(sca)
                    p[:, 3] = 1
                    return 11, p.astype('f')
                except TypeError: # sca is not a vector
                    p = numpy.array( self.prop[sourceProp])*sca
                    return 10, p.astype('f')

            else:
                sourceProp = numpy.array((1., 1., 0.), 'f')
                p = (self.prop[sourceProp]*sca).astype('f')
                return self.binding[sourceProp], p

        
    def fixOpacity(self):
        """This function should be called after all properties of a material
have been updated. It will create the alpha chanel of the diffuse
component from the opacity property. It will return 1 if there are
alpha values < 1.0 or 0 else"""
        #print "fixOpacity"

        difflen = self.prop[self.diff].shape[0]
        opaclen = self.prop[self.opac].shape[0]
        
        if self.binding[self.opac] == viewerConst.PER_VERTEX and difflen == 1 and opaclen >1:
            # PER_VERTEX bindig && one color ( len(diffuse) == 1) && len(opacity vector) > 1:
            # create diffuse vector which length is equal to len(opacity). Set the alpha
            # chanel of this vector to the values from the opacity vector.
            opac = self.prop[self.opac]
            diff = numpy.ones( (opaclen,4), 'f')
            diff[:,:3] = self.prop[self.diff][0][:3]
            diff[:,3] = opac
            self.prop[self.diff] = diff
            
        elif self.binding[self.opac] == viewerConst.OVERALL or \
                 difflen!=opaclen:
            # OVERALL binding OR  different length of diffuse and opacity vectors:
            # set the alpha channel of the diffuse component to the first value of
            # the opacity vector
            diff = self.prop[self.diff]
            opac = self.prop[self.opac][0]
            alpha = numpy.ones( (diff.shape[0]) )*float(opac)
            diff[:,3] = alpha.astype('f')
        else:

            self.prop[self.diff][:,3] = self.prop[self.opac][:]
#                raise ValueError("diffuse/opacity property shape mismatch")

        alphaMin = numpy.minimum.reduce( self.prop[self.opac] )
        if alphaMin < 1.0: return 1 
        else: return 0
            
##          elif opaclen==1: # duplicate opac[0]
##              ar = numpy.ones( (difflen,), 'f')*self.prop[self.opac][0]
##              self.prop[self.diff][:,3] = ar.astype('f')

##          elif difflen==1: # we have 1 rgb for diffuse but many opacities
##                           # so we duplicate diffuse[rgb][0]
##              ar = numpy.ones( (opaclen,4), 'f')
##              ar[:,3] = ar[:,3] * self.prop[self.diff][:,3]
##              ar[:, 3] = self.prop[self.opac][:]
##              self.prop[self.diff] = ar
##              self.binding[self.diff] = viewerConst.PER_VERTEX
##          else:
##              if self.prop[self.diff].shape[1]==4:
##                  self.prop[self.opac] = self.prop[self.diff][:,3]
##              else:
##                  self.prop[self.opac] = numpy.ones( (1,), 'f')
##                  self.binding[self.opac] = viewerConst.OVERALL
##              #raise ValueError("diffuse/opacity property shape mismatch")
        
##          alphaMin = numpy.minimum.reduce( self.prop[self.opac] )
##          if alphaMin < 1.0: return 1 
##          else: return 0

    
    def AddMaterial(self, values, prop=1):
	"""Add materials to the current set"""

        assert prop in (0,1,2,3,4,5)
        
        values = numpy.array( values, 'f' )
	if prop < self.shini:
	    if values.shape[1] == 3:
		alpha = numpy.ones( (values.shape[0], 1), 'f' )
		values = numpy.concatenate( (values, alpha), 1 )

        self.prop[prop] = numpy.concatenate( (self.prop[prop], values) )


    def SetMaterial(self, values, prop=1, tagModified=True):
	"""Set the materials
WARNING: when back face colors are set, two sided lighting has to be enabled
        
we set RGB values for all properties except for shininess and opacity
If an alpha value are specified they will be ignored
Since IndexedGeomDSPL requires alpha values for all properties
we set them automatically to 1.0 for all properties except for
diffuse. The alpha channel of the diffuse component will be set later
in fixOpacity which should be called after all properties of the
material have been updated.
"""
        if tagModified:
            self._modified = True

        if prop is None: prop = self.diff
        elif type(prop) is types.StringType:
            prop = getattr(self, prop)

        assert prop in (0,1,2,3,4,5)

        values = numpy.array( values, 'f' )
	if prop < self.shini:
            assert len(values.shape)==2 and values.shape[1] in [3,4]
            alpha = numpy.ones( (values.shape[0], 1), 'f' )
            values = numpy.concatenate( (values[:,:3], alpha), 1 )
        else:
            if len(values.shape) != 1:
                values = numpy.array([values], 'f' )
            
        self.prop[prop] = values


    def Set(self, check=1, **kw):

        tagModified = True
        val = getkw(kw, 'tagModified')
        if val is not None:
            tagModified = val
        assert tagModified in [True, False]
        self._modified = tagModified

        if __debug__:
            if check:
                apply( checkKeywords, ('Material',self.keywords), kw)

	val = getkw(kw, 'ambient')
	if not val is None:
            val = numpy.array(val).astype('f')
            self.prop[propertyNum['ambi']] = val

	val = getkw(kw, 'diffuse')
	if not val is None:
            val = numpy.array(val).astype('f')
            self.prop[propertyNum['diff']] = val

	val = getkw(kw, 'emission')
	if not val is None:
            val = numpy.array(val).astype('f')
            self.prop[propertyNum['emis']] = val

	val = getkw(kw, 'specular')
	if not val is None:
            val = numpy.array(val).astype('f')
            self.prop[propertyNum['spec']] = val

	val = getkw(kw, 'shininess')
	if not val is None:
            val = numpy.array(val).astype('f')
            self.prop[propertyNum['shini']] = val

	val = getkw(kw, 'opacity')
	if not val is None:
            val = numpy.array(val).astype('f')
            self.prop[propertyNum['opac']] = val

	val = getkw(kw, 'binding')
	if not val is None:
            val = numpy.array(val).astype('f')
            self.binding[:] = val[:]
            

    def getState(self, includeBinding=True):
        """return a dictionary describing this object's state
This dictionary can be passed to the Set method to restore the object's state
"""
        if includeBinding:
            return {
                'ambient':self.prop[propertyNum['ambi']].tolist(),
                'diffuse':self.prop[propertyNum['diff']].tolist(),
                'emission':self.prop[propertyNum['emis']].tolist(),
                'specular':self.prop[propertyNum['spec']].tolist(),
                'shininess':self.prop[propertyNum['shini']].tolist(),
                'opacity':self.prop[propertyNum['opac']].tolist(),
                'binding': self.binding
                }
        else:
            lDict = {'binding': self.binding}
            if     (self.binding[0] != viewerConst.PER_VERTEX) \
               and (self.binding[0] != viewerConst.PER_PART):
                lDict['ambient'] = self.prop[propertyNum['ambi']].tolist()
            if     (self.binding[1] != viewerConst.PER_VERTEX) \
               and (self.binding[1] != viewerConst.PER_PART):
                lDict['diffuse'] = self.prop[propertyNum['diff']].tolist()
                lDict['opacity'] = self.prop[propertyNum['opac']].tolist()
            if     (self.binding[2] != viewerConst.PER_VERTEX) \
               and (self.binding[2] != viewerConst.PER_PART):
                lDict['emission'] = self.prop[propertyNum['emis']].tolist()
            if     (self.binding[3] != viewerConst.PER_VERTEX) \
               and (self.binding[3] != viewerConst.PER_PART):
                lDict['specular'] = self.prop[propertyNum['spec']].tolist()
            if     (self.binding[4] != viewerConst.PER_VERTEX) \
               and (self.binding[4] != viewerConst.PER_PART):
                lDict['shininess'] = self.prop[propertyNum['shini']].tolist()
            return lDict
