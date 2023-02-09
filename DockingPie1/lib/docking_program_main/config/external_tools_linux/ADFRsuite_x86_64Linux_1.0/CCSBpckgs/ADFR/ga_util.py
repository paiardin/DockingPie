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
# Copyright: M. Sanner and TSRI 2015
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/ADFR/ga_util.py,v 1.2 2016/12/07 00:38:32 sanner Exp $
#
# $Id: ga_util.py,v 1.2 2016/12/07 00:38:32 sanner Exp $
#
#base definitions for genetic algorithms

from random import random
from math import floor
import numpy

compress=numpy.compress
average=numpy.average

class empty_class: pass

GAError = 'GA Error'

## def nop(x):
##     return x

class FlipCoin:
    def __init__(self, seed):
        from random import Random
        self.random = Random()
        self.setSeed(seed)

    def setSeed(self, seed):
	self.random.seed(seed)

    def flip_coin(self, p):
        return (self.random.random() < p)

def flip_coin(p):
    return (random() < p)

## def shallow_clone(item):
##     new = empty_class()
##     new.__class__ = item.__class__
##     new.__dict__.update(item.__dict__)
##     return new



## NOTE: assume nan or inf never appears
def my_mean(s):
    return numpy.average(s)

def my_std(s):
    return standardDeviation(s)

# a simplified version of scipy.stats.rv._pranv:_choice()
def choice(seq=(0,1)):
    """Return element k with probability 1/len(seq) from non-empty sequence.
    choice(seq=(0,1), size=None)
    """
    lenseq = len(seq)
    if lenseq == 0:
        raise ValueError, '<seq> must not be empty'
    else:
        randomNumber = random()
        #print 'CHOICE random:', randomNumber
        index=randomNumber*lenseq
        return seq[int(floor(index))]
        


## #these are exacly correct, but htey prevent problems with -Inf and Inf
## def my_std(s):
##     a = remove_NaN(s)
##     if len(a) > 1:
##         return stats.std(a)
##     else:
##         return 0.
    #import pdb
    #pdb.set_trace()



## NOTE: THIS IS NAN, INF FRIENDLY    
## def my_mean(s):
##     a = remove_NaN(s)
##     if len(a) > 1:
##         return stats.mean(a)
##     else:
##         return 0.
	
def testflip():
	import time
	b = time.clock()
	for i in range(10000): a = flip_coin(.5)
	e = time.clock()
	print 'rv_flip',e-b
	b = time.clock()
	for i in range(10000): a = flip_coin2(.5)
	e = time.clock()
	print 'wh_flip',e-b
	from rv import random
	b = time.clock()
	for i in range(10000): 
		a = random() < .5
	e = time.clock()
	print 'rv',e-b
	from whrandom import random
	b = time.clock()
	for i in range(10000): 
		a = random() < .5
	e = time.clock()
	print 'wh',e-b




## the following code is taken from MMTK / Scientific.Statistics.__init__.py

def moment(data, order, about=None, theoretical=1):
    data = 1.*numpy.array(data)
    if about is None:
       about = mean(data)
       theoretical = 0
    ln = len(data)-(1-theoretical)
    return 1.*numpy.add.reduce((data-about)**order)/ln

def mean(data):
    "Returns the mean (average value) of |data| (a sequence of numbers)."
    return moment(data, 1, 0)

average = mean

def weightedMean(data, sigma):
    """Weighted mean of a sequence of numbers with given standard deviations.

    |data| is a list of measurements,
    |sigma| a list with corresponding standard deviations.

    Returns weighted mean and corresponding standard deviation.
    """
    #from numpy import array, Float, sqrt, sum

    if len(data) != len(sigma):
        raise ValueError
    data = 1.*numpy.array(data)
    sigma = 1.*numpy.array(sigma)
    nom = sum(data/sigma**2)
    denom = sum(1./sigma**2)
    mean = nom / denom
    sig = sqrt(1./denom)
    return mean, sig

def variance(data):
    "Returns the variance of |data| (a sequence of numbers)."
    return moment(data, 2)

def standardDeviation(data):
    "Returns the standard deviation of |data| (a sequence of numbers)."
    return numpy.sqrt(variance(data))

def median(data):
    "Returns the median of |data| (a sequence of numbers)."
    data = numpy.sort(numpy.array(data))
    l = (len(data)-1)/2.
    return (data[int(numpy.floor(l))]+data[int(numpy.ceil(l))])/2.

def mode(data):
    h = {}
    for n in data:
       try: h[n] = h[n]+1
       except KeyError: h[n] = 1
    a = map(lambda x: (x[1], x[0]), h.items())
    return max(a)[1]

def normalizedMoment(data, order):
    mn = mean(data)
    return moment(data, order, mn)/(moment(data, 2, mn)**order)

def skewness(data):
    "Returns the skewness of |data| (a sequence of numbers)."
    return normalizedMoment(data, 3)

def kurtosis(data):
    "Returns the kurtosis of |data| (a sequence of numbers)."
    return normalizedMoment(data, 4)

##  def chiSquare(data):
##      h = {}
##      for n in data:
##         try: h[n] = h[n]+1
##         except KeyError: h[n] = 1
##      h = numpy.array(h.values())
##      h = h/numpy.add.reduce(h)
##      return moment(h, 2, 1./len(h))

def correlation(data1, data2):
    """Returns the correlation coefficient between |data1| and |data2|,
    which must have the same length."""
    if len(data1) != len(data2):
        raise ValueError, "data series must have equal length"
    data1 = numpy.array(data1)
    data2 = numpy.array(data2)
    data1 = data1 - numpy.add.reduce(data1)/len(data1)
    data2 = data2 - numpy.add.reduce(data2)/len(data2)
    return numpy.add.reduce(data1*data2) / \
           numpy.sqrt(numpy.add.reduce(data1*data1) \
                        * numpy.add.reduce(data2*data2))
