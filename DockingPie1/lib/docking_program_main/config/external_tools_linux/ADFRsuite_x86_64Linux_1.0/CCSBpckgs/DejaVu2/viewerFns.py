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
# $Header: /mnt/raid/services/cvs/DejaVu2/viewerFns.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#
# $Id: viewerFns.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#

import numpy
import viewerConst
import types

def checkKeywords(_name, keywords, **kw):
    """test is all kyes in **kw are in list keywords"""

    for key in kw.keys():
        if key not in keywords:
            print 'WARNING: Keyword %s not recognized for %s' % (key,_name)
            #import traceback;traceback.print_stack()
            

def getkw(kw, name):
    """get a dictionary entry and remove it from the dictionary"""

    v = kw.get(name)
    if name in kw.keys(): del kw[name]
    return v


def GetArray(vector, shape=None, precision=viewerConst.FPRECISION ):
    """get a vector and return it in a numpy array"""

    if type(vector).__name__ == 'array':
	if vector.dtype.char==precision:
	    if shape:
		if vector.shape == shape: return vector
		else: return numpy.reshape(vector, shape)
	    else: return vector
	else:
	    vector = numpy.array( vector, precision )
	    if shape:
		if vector.shape == shape: return vector
		else: return numpy.reshape(vector, shape)
	    else: return vector
    else:
	vector = numpy.array( vector, precision )
	if shape:
	    if vector.shape == shape: return vector
	    else: return numpy.reshape(vector, shape)
	else: return vector

    
def read_data(filename, function):
    """Read data from an ASCII file."""

    import string
    result = []
    assert type(filename) == types.StringType
    assert function in (None, float, int)
    print "reading", filename
    f = open(filename, "r")
    while 1:
	line = f.readline();
	if not len(line):
	    break
	if function:
	    datum = map(function, string.split(line))
	    result.append(list(datum))
	else:
	    result.append(line)
    f.close()
    if function == float:
	ar = numpy.array(result).astype(viewerConst.FPRECISION)
    elif function == int:
	ar = numpy.array(result).astype(viewerConst.IPRECISION)
    else:
	ar = result

    return ar


def write_data(filename,data):
    """Write data to an ASCII file."""

    assert type(filename) == types.StringType
    import string
    f = open(filename, "w")
    def writedatum1(datum, f=f):
	line = str(datum)
	line = line[1:-1] + "\n" # slice off the [ , ] characters
	f.write(line)
    def writedatum(datum, f=f):
	import string, numpy
	line = string.join(map(str, numpy.array(datum))," ") + "\n"
	f.write(line)
    map(writedatum,data)
    f.close()

