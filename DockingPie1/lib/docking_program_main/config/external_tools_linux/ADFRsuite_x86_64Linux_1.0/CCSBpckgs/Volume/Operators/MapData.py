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
# Author: A. Omelchenko, M. Sanner
#
# Copyright: M. Sanner TSRI
#
#############################################################################

#
#$Header: /mnt/raid/services/cvs/Volume/Operators/MapData.py,v 1.12.16.1 2017/07/28 01:09:21 annao Exp $
#
#$Id: MapData.py,v 1.12.16.1 2017/07/28 01:09:21 annao Exp $
#
import numpy
import struct
from time import time
from math import log

def isPowerOf2(num):
    if num == 1:
        power = 0
        return 1,  power
    result = 1
    ntemp = num
    power = 0
    while ntemp > 1:
        if ntemp % 2 != 0: result = 0
        ntemp = ntemp/2
        power = power + 1
    if result == 1:
        return 1, power
    power = power + 1
    return 0, power


class MapGridData:
    """An instance of this class can be used to:
- change the data type in a Numeric array
- map the data to a new range of values using a linear or log mapping. This
  operation can be performed in place or create a new array
- pad the array to have dimensions that are powers of 2
"""

    def __call__(self, data, datatype=None, datamap=None, powerOf2=False,
                 fillValue=0.0):
        """ Map(data, datamap=None, datatype=None, powerOf2=False)
    data -- 3D array, source.
    datatype -- The Numeric data type to be used for the destination array
    datamap -- Dictionary that defines the following keys:
                   src_min, src_max: range of values in the source map. None
                   for these values will use the full range
                   dst_min, dst_max: range of values in the destination map
                   map_type: type of mapping. Can be 'linear' or 'log'.
    powerOf2 -- if set to 1, then if the data array dimensions are not
                power of 2 - the returned array will be padded with zeros
                so that its dims are power of 2.
"""
        #if no mapping is required: 
        if datamap is None:
            if datatype:
                if isinstance(data, numpy.ndarray):
                   if data.dtype.char!=datatype: # do the type conversion
                       data = data.astype(datatype)
                else:
                    data = numpy.array( data, datatype )
                
        if not datatype:
            if isinstance(data, numpy.ndarray):
                datatype = data.dtype.char
                
        nx, ny, nz = data.shape        
        # map data to range
        if datamap is not None:
            src_min = datamap['src_min']
            src_max = datamap['src_max']
            dst_min = datamap['dst_min']
            dst_max = datamap['dst_max']
            mtype   = datamap['map_type']

            #nx, ny, nz = data.shape
            arrsize = nx*ny*nz
            #arr_max = numpy.maximum.reduce(data.ravel())
            #arr_min = numpy.minimum.reduce(data.ravel())
            maxif = numpy.maximum.reduce 
            arr_max = maxif(maxif(maxif(data)))
            minif = numpy.minimum.reduce
            arr_min = minif(minif(minif(data)))

            # check src_min
            if src_min is not None:
                assert src_min >= arr_min and src_min < arr_max, \
                       "%f>=%f and %f<=%f failed"%(src_min, arr_min, src_min,
                                                   arr_max)

                if src_max is not None:
                    assert src_min < src_max
            else:
                src_min = arr_min

            # check src_max
            if src_max is not None:
                #assert src_max >= arr_min and src_max < arr_max
                assert src_max >= arr_min and src_max <= arr_max, \
                       "%f>=%f and %f<=%f failed"%(src_max, arr_min, src_max,
                                                   arr_max)
                if src_min is not None:
                    assert src_min < src_max
            else:
                src_max = arr_max

            # check
            assert dst_min < dst_max
            #print type(src_min), type(src_max), type(dst_min), type(dst_max)
            #print "mapping data range %g -- %g to %g -- %g "%(src_min, src_max,
            #                                                  dst_min, dst_max)

            if mtype == 'linear':
                data = self.linearMap( data, src_min, src_max, dst_min,dst_max,
                                       arr_min, arr_max, datatype)
            elif mtype == 'log':
                data = self.logMap( data, src_min, src_max, dst_min, dst_max,
                                    arr_min, arr_max, datatype)

        # pad to make dimensions powers of 2
        if powerOf2:
            nx1, ny1, nz1 = data.shape
            res, power = isPowerOf2(nx1)
            if not res:
                nx1 = 2**power
            res, power = isPowerOf2(ny1)
            if not res:
                ny1 = 2**power
            res, power = isPowerOf2(nz1)
            if not res:
                nz1 = 2**power
            dx, dy, dz = 0, 0, 0
            if nx1 != nx or ny1 != ny or nz1 != nz:
                narr = (numpy.ones( (nx1,ny1,nz1) )*fillValue).astype(datatype)
                narr[:nx, :ny,:nz] = data
                data = narr
##             if nx1 != nx or ny1 != ny or nz1 != nz:
##                 dx = (nx1-nx)
##                 dy = (ny1-ny)
##                 dz = (nz1-nz)
##                 #narr = numpy.zeros( (nx1,ny1,nz1), data.dtype.char)
##                 narr = numpy.zeros( (nx1,ny1,nz1), datatype)
##                 # copy original data
##                 narr[:nx,:ny,:nz] = data[:,:,:]
##                 # duplicate last face in 3 directions
##                 for i in xrange(nx, nx1):
##                     narr[i,:ny,:nz] = data[-1,:,:]
##                 for i in xrange(ny, ny1):
##                     narr[:nx,i,:nz] = data[:,-1,:]
##                 for i in xrange(nz, nz1):
##                     narr[:nx,:ny,i] = data[:,:,-1]
##                 # duplicate edge values
##                 for i in xrange(nx):
##                     carr = numpy.ones( (dy, dz) )*data[i, -1, -1]
##                     #narr[i, ny:, nz:] = carr[:,:].astype(data.dtype.char)
##                     narr[i, ny:, nz:] = carr[:,:].astype(datatype)
##                 for i in xrange(ny):
##                     carr = numpy.ones( (dx, dz) )*data[-1, i, -1]
##                     #narr[nx:, i, nz:] = carr[:,:].astype(data.dtype.char)
##                     narr[nx:, i, nz:] = carr[:,:].astype(datatype)
##                 for i in xrange(nz):
##                     carr = numpy.ones( (dx, dy) )*data[-1, -1, i]
##                     #narr[nx:, ny:, i] = carr[:,:].astype(data.dtype.char)
##                     narr[nx:, ny:, i] = carr[:,:].astype(datatype)
##                 # duplicate far corner arr[-1, -1, -1] in narr[nx:, ny:, nz:]
##                 carr = numpy.ones( (dx, dy, dz) )*data[-1, -1, -1]
##                 #narr[nx:, ny:, nz:] = carr[:,:,:].astype(data.dtype.char)
##                 narr[nx:, ny:, nz:] = carr[:,:,:].astype(datatype)
##                 data = narr
        return data

    
    def linearMap(self, data, data_min, data_max, val_min, val_max,
                  arr_min, arr_max, datatype):
        k2,c2 = self.ScaleMap((val_min, val_max), (data_min, data_max))
        if k2 == 1 and c2 == 0:
            return data
        new_arr = numpy.clip(k2*data+c2, k2*data_min+c2, k2*data_max+c2)
        #return new_arr.astype(data.dtype.char)
        return new_arr.astype(datatype)
    
    def logMap(self, data, data_min, data_max, val_min, val_max,
               arr_min, arr_max, datatype):
        if arr_min < 0:
            diff = abs(arr_min)+1.0
        elif arr_min >= 0 and arr_min < 1.0:
            diff = 1.0
        elif arr_min >= 1.0:
            diff=0
        k1, c1 = self.ScaleMap( (val_min, val_max),
                                (log(arr_min+diff), log(arr_max+diff)) )

        return  (k1*numpy.log10(data+diff)+c1).astype(datatype)



    def ScaleMap(self, val_limit, data_limit, par=1):
        """ Computes coefficients for linear mapping"""
        
        assert data_limit[1]-data_limit[0]
        k = float (val_limit[1]-val_limit[0])
        k = k / (data_limit[1]**par-data_limit[0]**par)
        c = val_limit[1]-k*data_limit[1]**par
        return (k,c)



class MapData :
    
    """ Class for mapping a numeric array of floats to an array of integers.
    A special case where the user can define three mapping intervals for
    the source data array and for the resulting integer array.
    If mapping intervals are specified by supplying data min and max values
    (other than global data minimum and maximum values) and scalar Val_min and
    val_max (in range int_min ... int_limit), then the mapping will proceed as follows:
        [global_data_min, data_min] is mapped to [int_min, val__min],
        [data_min,data_max] is maped to [val _min, val_max],
        [data_max,global_data_max] is mapped to [val_max, int_limit]."""
    
    def __init__(self, int_limit=4095, int_min = 0, int_type = numpy.int16):
        """(int_limit=4095, int_type = numpy.int16)
        int_limit - maximum int value of the resulting integer array;
        int_type - numpy typecode of the integer array. """
        
        self.int_limit = int_limit
        self.arr = None
        self.int_min = int_min
        self.val_max = int_limit
        self.int_type = int_type
    
    def __call__(self, data, data_min=None, data_max=None,
                 val_min=None, val_max=None, map_type='linear', powerOf2=0):
        """ (data, data_min=None, data_max=None, val_min=None, val_max=None,
             map_type = 'linear', powerOf2=0)
        Maps an array of floats to integer values.
        data -- 3D numeric array;
        data_min, data_max -- min and max data values(other than actual
        array minimum/maximum) to map to integer values -
        val_min and val_max - in range (0 ... int_limit);
        map_type -- can be 'linear' or 'log';
        powerOf2 -- if set to 1, then if the data array dimensions are not
                    power of 2 - the returned array will be padded with zeros
                    so that its dims are power of 2.
        """
        # if data_min/max and val_min/max specified then the mapping
        # will proceed as follows:
        # [arr_min, data_min] is mapped to [int_min, val_min],
        # [data_min, data_max] is maped to [val_min, val_max],
        # [data_max, arr_max] is mapped to [val_max, int_limit].
        int_limit = self.int_limit
        int_min = self.int_min
        shape = data.shape
        assert len(shape)==3
        nx, ny, nz = shape
        arrsize = nx*ny*nz
        #arr_max = numpy.maximum.reduce(data.ravel())
        #arr_min = numpy.minimum.reduce(data.ravel())
        maxif = numpy.maximum.reduce 
        arr_max = maxif(maxif(maxif(data)))
        minif = numpy.minimum.reduce
        arr_min = minif(minif(minif(data)))
        #print "min(arr)=%f" % arr_min
        #print "max(arr)=%f" % arr_max
        if val_min != None:
            assert val_min >= 0 and val_min < int_limit
        else:
            val_min = int_min
        if val_max != None:
            assert val_max <= int_limit and val_max > 0
        else:
            val_max = int_limit
        if data_min != None:
            if data_min < arr_min: data_min = arr_min
        else:
            data_min = arr_min
        if data_max != None:
            if data_max > arr_max: data_max = arr_max
        else:
            data_max = arr_max
        print "mapping data_min %4f to val_min %d, data_max %4f to val_max %d"\
              % (data_min, val_min, data_max, val_max)
        if map_type == 'linear':
            k2,c2 = self.ScaleMap((val_min, val_max), (data_min, data_max))
            n_intervals = 3
            if abs(data_min-arr_min) < 0.00001: # data_min==arr_min
                k1,c1 = k2, c2
                n_intervals = n_intervals-1
            else :
                k1, c1 = self.ScaleMap((int_min, val_min),
                                       (arr_min, data_min)) 
            if abs(data_max-arr_max) < 0.00001: # data_max == arr_max
                k3, c3 = k2, c2
                n_intervals = n_intervals-1
            else:
                k3, c3 = self.ScaleMap((val_max, int_limit),
                                       (data_max, arr_max))

            t1 = time()
            #print "n_intervals = ", n_intervals
            if n_intervals == 2:
                if data_max == arr_max:
                    #print "data_max == arr_max"
                    new_arr = numpy.where(numpy.less(data, data_min),
                                            k1*data+c1, k2*data+c2 )
                elif data_min == arr_min:
                    #print "data_min == arr_min"
                    new_arr = numpy.where(numpy.greater_equal(data,
                              data_max), k3*data+c3, k2*data+c2)
            elif n_intervals == 3:
                new_arr1 = numpy.where(numpy.less(data, data_min),
                                        k1*data+c1, k2*data+c2)
                new_arr = numpy.where(numpy.greater_equal(data, data_max),
                                        k3*data+c3, new_arr1)
                del(new_arr1)
            else :
                new_arr = k2*data+c2
            arr = numpy.transpose(new_arr).astype(self.int_type)
            del(new_arr)
            t2 = time()
            print "time to map : ", t2-t1
        elif map_type == 'log':
            if arr_min < 0:
                diff = abs(arr_min)+1.0
            elif arr_min >= 0 and arr_min < 1.0:
                diff = 1.0
            elif arr_min >= 1.0:
                diff=0
            k1, c1 = self.ScaleMap( (int_min, int_limit),
                           (log(arr_min+diff), log(arr_max+diff)) )
            arr=numpy.transpose(k1*numpy.log10(data+diff)+c1).astype(self.int_type)
        self.data_min = data_min
        self.data_max = data_max
        self.val_min = val_min
        self.val_max = val_max    
        if powerOf2:
            nx1, ny1, nz1 = nx, ny, nz
            res, power = isPowerOf2(nx)
            if not res:
                nx1 = 2**power
            res, power = isPowerOf2(ny)
            if not res:
                ny1 = 2**power
            res, power = isPowerOf2(nz)
            if not res:
                nz1 = 2**power
            dx, dy, dz = 0, 0, 0
            if nx1 != nx or ny1 != ny or nz1 != nz:
                #print "new data size: ", nx1,ny1,nz1
                dx = (nx1-nx)/2. ; dy = (ny1-ny)/2. ; dz = (nz1-nz)/2.
                #narr = numpy.zeros((nx1,ny1,nz1), self.int_type)
                #narr[:nx,:ny,:nz] = arr[:,:,:]
                narr = numpy.zeros((nz1,ny1,nx1), self.int_type)
                narr[:nz,:ny,:nx] = arr[:,:,:]
                self.arr = narr
                #arr = numpy.zeros((nx1,ny1,nz1), self.int_type)
                #arr[:nx,:ny,:nz] = new_arr[:,:,:]
                #arr = numpy.transpose(arr).astype(self.int_type)
                #self.arr = arr
                return narr
        self.arr = arr
        return arr


    def ScaleMap(self, val_limit, data_limit, par=1):
        """ Computes coefficients for linear mapping"""
        
        assert data_limit[1]-data_limit[0]
        k=(val_limit[1]-val_limit[0])/(data_limit[1]**par-data_limit[0]**par)
        c=val_limit[1]-k*data_limit[1]**par
        return (k,c)

