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
# Revision: Guillaume Vareille
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/DejaVu2/datamodel.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#
# $Id: datamodel.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#

"""Classes to store data such as vertices, faces, properties etc..."""

import numpy
import viewerConst
import types

import DejaVu2

def isarray(var):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	return type(var).__name__ == 'array'

class ValuesArray:
    """base class for arrays of data

 ValuesArray is a class used to represent data that can be stored in a
 numpy array. The shape of the array holding the data and the precision
 is set at creation time and only data that matches this shape can be added
 to this array
"""


    def _FirstValues(self, data=None, shape=None, datatype=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Set the initial set of values, always returns a copy"""

	if data is not None:
	    if len(data)==0:
	        self.ashape = (0,) + self.ashape[1:]
		return numpy.zeros( self.ashape, self.type)
	    elif isarray(data):  # data is an numpy array
		if not datatype:
		    if self.type: datatype = self.type
		    else: raise AttributeError('I Need a datatype')

		if data.typecode() == self.type: # right type
		    if shape is None: self.ashape = data.shape
		    if data.shape == self.ashape: return numpy.array(data)
		    else: return numpy.reshape(data, shape )
		else: # bad type
		    data = data.astype(datatype)
		    if shape is None: self.ashape = data.shape
#		    if data.shape == self.ashape: return numpy.array(data)
		    if data.shape == self.ashape: return data
		    else: return numpy.reshape(data, shape)

	    else: # data is not an array
		if datatype and shape:
		    data = numpy.array( data, self.type)
		    return numpy.reshape( data, shape )
		elif datatype:
		    d=numpy.array( data, datatype)
		    self.ashape = d.shape
		    return d
		else:
                    try:
                        d=numpy.array( data,  self.type)
                    except ValueError:
                        raise ValueError('Bad argument')
		    if d.dtype.char in numpy.typecodes['Integer'] or \
		       d.dtype.char in numpy.typecodes['UnsignedInteger']:
			d=numpy.array( data, viewerConst.IPRECISION)
		    else:
			d=numpy.array( data, viewerConst.FPRECISION)
		    self.ashape = d.shape
		    return d

	else: # No data
		if datatype and shape:
		    d = numpy.zeros( (0,), self.type)
		    d.shape = shape
		    return  d
		else:
		    raise AttributeError('I Need a datatype and a shape')


    def __init__(self, name, data=None, shape=None, datatype=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Create an instance, with a name and (shape or data)"""


	assert type(name) == types.StringType

	if shape is not None and data is not None:
	    raise AttributeError("need shape or data attribute")

#	if __debug__:
#	    if data: assert type(data).__name__ == 'array'

	self.name = name
	self.array = None
	self.type = datatype
	self.ashape = shape   # array's shape
	self.array = self._FirstValues(data, shape, datatype)
	self.fixedShape = self.ashape[1:]


    def __len__(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""returns the number of entries"""
	return int(numpy.multiply.reduce(self.array.shape))


    def __repr__(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""List name, number of values and the shape of this array"""

	#d = self.array
	return '<Array %s %d values of type %s shape=%s>' % \
	       (self.name, len(self), self.type, str(self.array.shape) )


    def SetValues(self, values, reshape=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Set values"""

	if values is None: return
	d = self._FirstValues(values)
	if len(self.ashape[1:]) != len(self.fixedShape):
	    raise AttributeError("input array bad shape %s!=%s" % \
				 (str(self.ashape[1:]),
				  str(self.fixedShape)))

	if len(self.array):
	    if not reshape and d.shape[1:] != self.fixedShape:
		raise AttributeError("input array of bad shape")

	self.fixedShape = d.shape[1:]
	self.array = d
	self.ashape = self.array.shape


    def _ShapeValues(self, values):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Check type and shape, coerce if needed"""

	if isarray(values) and values.dtype.char == self.type:
	    if values.shape[1:] == self.fixedShape:
		    return values
	    else:
		    return numpy.reshape( values, ( -1,) + self.fixedShape )
	else:
	    values = numpy.array( values, self.type)
	    return numpy.reshape( values, ( -1,) + self.fixedShape )


    def AddValues(self, values):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Add values"""

	if values is None: return
	if len(self.array)==0:
	    self.SetValues(values)
	    self.fixedShape = self.ashape[1:]
	else:
	    self.array = numpy.concatenate( (self.array,
					       self._ShapeValues(values)) )
	    self.ashape = self.array.shape


class PropertiesArray(ValuesArray):
    """base class for property Arrays"""


    def __init__(self, name, data=None, shape=None, datatype=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

	ValuesArray.__init__(self, name, data, shape, datatype)
	self.Compute = None          # user settable function to compute prop.
	self.args = ()               # arguments to function Compute
	self.computeMode = None      # computation mode
	self.status = viewerConst.UNKNOWN


    def ComputeFunction(self, func, args):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Set the function used to compute this property"""

	assert callable(func)
	assert type(args).__name__=='tuple'
	self.Compute = func
	self.args = args


    def ComputeMode(self, mode):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Set the computation mode for this property"""

	if mode in ( viewerConst.NO_COMPUTATION, viewerConst.AUTO ):
	    self. computeMode = mode
	else:
	    raise AttributeError("Bad property computation mode NO_COMPUTATION,\
 AUTO")


    def PropertyStatus(self, neededNumber):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Check if we have as many property values as we need"""

	if self.array is None or neededNumber == 0:
	    self.status = viewerConst.NONE
	elif len(self) == neededNumber: self.status = viewerConst.SET
	elif len(self) > neededNumber: self.status = viewerConst.TOO_MANY
	elif len(self) < neededNumber: self.status = viewerConst.TOO_FEW


    def AddValues(self, values, neededNumber=-1):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Add values"""
	ValuesArray.AddValues(self, values)
	if neededNumber>0: self.PropertyStatus(neededNumber)


    def GetProperty(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""May be compute and return the values"""

	if self.status < viewerConst.COMPUTED:
	    if self.Compute and self.computeMode==viewerConst.AUTO:
		d = apply(self.Compute, self.args)
		if d is not None:
		    self.SetValues(d)
		    self.status=viewerConst.COMPUTED
		else:
		    self.status=viewerConst.NONE
		    return None
	    else:
		self.status=viewerConst.NONE
		return None

	return self.array


class VectorProperties(PropertiesArray):
    """base class for arrays of vector property values"""


    def __init__(self, name, data=None, shape=None, datatype=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

	PropertiesArray.__init__(self, name, data, shape, datatype)
	if len(self.ashape) < 2:
	    raise AttributeError('VectorProperties require at least a 2D array')

    def __len__(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""returns the number of property vectors"""
	return int(numpy.multiply.reduce(self.array.shape[:-1]))


    def __repr__(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        """List name, number of values and the shape of this array"""

	#d = self.array
	return '<VectorProperties %s %d %d-vectors of type %s shape=%s \
status=%s>' % \
	       (self.name, len(self), self.ashape[-1], self.type,
	       str(self.ashape), self.status )


    def __mul__(self, mat):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Apply the transformation matrix to the vertices"""

	assert self.array.shape[-1]==3
	self.Flatten()
	one = numpy.ones( (self.array.shape[0], 1), viewerConst.FPRECISION )
	c = numpy.concatenate( (self.array, one), 1 )
	self.UnFlatten()
	return numpy.dot( c, numpy.transpose(mat) )[:, :3]


    def Flatten(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Reshape the coordinates to a 2D array of vectors"""
	self.oldshape = self.array.shape
	self.array = numpy.reshape( self.array,  (-1, self.oldshape[-1] ) )


    def UnFlatten(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Reshape the coordinates to initial shape array"""
	self.array = numpy.reshape( self.array, self.oldshape)


    def BoundingBox(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Compute and return the bounding box of the vectors"""

	self.Flatten()
	self.bbmax = numpy.maximum.reduce(self.array)
	self.bbmin = numpy.minimum.reduce(self.array)
	self.UnFlatten()
	return self.bbmin, self.bbmax


class ScalarProperties(PropertiesArray):
    """base class for arrays of scalar property values"""


    def __init__(self, name, data=None, shape=None, datatype=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

	PropertiesArray.__init__(self, name, data, shape, datatype)


    def __repr__(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""List name, number of values and the shape of this array"""

	#d = self.array
	return '<ScalarProperties %s %d values of type %s shape=%s \
status=%d>' % \
	       (self.name, len(self), self.type, str(self.ashape),self.status )


class Set:
    """base class for VertexSet and Face Set"""

    def __init__(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	self.properties = {}


    def AddScalarProperty(self, name, data=None, shape=None, datatype=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""add a scalar property"""

	assert type(name) == types.StringType
	if name in self.properties.keys():
	    raise ValueError('Property name already used')
	p = ScalarProperties(name, data, shape, datatype)
	eval ( 'self.%s = p', name )
	self.properties[name] = p


    def AddVectorProperty(self, name, data=None, shape=None, datatype=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""add a vector property"""

	assert type(name) == types.StringType
	if name in self.properties.keys():
	    raise ValueError('Property name already used')
	p = VectorProperties(name, data, shape, datatype)
	eval ( 'self.%s = p', name )
	self.properties[name] = p


    def ListProperties(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""List names and summary of defined properties"""
	for p in self.properties.keys():
	    print repr(self.properties[p])


class VertexSet(Set):
    """class for vertices"""

    def __init__(self, vert=None, shape=None, norm=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

	Set.__init__(self)
	if vert is not None and len(vert):
	    self.vertices = VectorProperties('vertices', data=vert,
				      datatype=viewerConst.FPRECISION)
	elif shape is not None:
	    self.vertices = VectorProperties( 'vertices', shape=shape,
				       datatype=viewerConst.FPRECISION)
	else:
	    raise AttributeError("I need data or a shape")

	if norm is not None and len(vert):
#	    self.AddVectorProperty('normals', norm,
#				   datatype=viewerConst.FPRECISION)
	    self.normals = VectorProperties('normals', norm,
					    datatype=viewerConst.FPRECISION)
	else:
	    nshape = len(self.vertices.ashape) * (0,)
#	    self.AddVectorProperty('normals', norm,
#				   datatype=viewerConst.FPRECISION)
	    self.normals = VectorProperties('normals', shape=nshape,
					    datatype=viewerConst.FPRECISION)

	if len(self.vertices) > 0: self.vertices.BoundingBox()


    def __repr__(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	return '<%s> with %d elements prop: %s' % \
              (self.__class__, len(self.vertices), str(self.properties.keys()))


    def __len__(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""returns the number of vertices"""
	return len(self.vertices)


    def AddValues(self, values):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Add values"""
	self.vertices.AddValues(values)
	if len(self.vertices) > 0: self.vertices.BoundingBox()


class FaceSet(Set):
    """class for vertices"""

    def __init__(self, faces=None, shape=None, norm=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

	Set.__init__(self)
	if faces:
	    self.faces = VectorProperties('faces', data=faces,
				         datatype=viewerConst.IPRECISION)
	elif shape:
	    self.faces = VectorProperties('faces', shape=shape,
					  datatype=viewerConst.IPRECISION)
	else:
	    raise AttributeError("I need data or a shape")

	if norm:
	    self.AddVectorProperty('normals', norm)
	else:
	    nshape = len(self.faces.ashape) * (0,)
	    self.normals = VectorProperties('normals', shape=nshape,
				   datatype=viewerConst.FPRECISION)


    def __repr__(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	return '<%s> with %d elements prop: %s' % \
              (self.__class__, len(self.faces), str(self.properties.keys()))


    def __len__(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""returns the number of faces"""
	return len(self.faces)

