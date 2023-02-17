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
# $Header: /mnt/raid/services/cvs/ADFR/scaling.py,v 1.3 2016/12/07 00:38:32 sanner Exp $
#
# $Id: scaling.py,v 1.3 2016/12/07 00:38:32 sanner Exp $
#

from ga_util import my_std, my_mean, GAError
from numpy import less_equal, choose

# if a score is less the 2 standard deviations below, the average, its score
# is arbitrarily set to zero
class sigma_truncation_scaling:

    def __init__(self,scaling = 2):
	self.scaling = scaling

    def scale(self, pop):
	sc = pop.scores()
	avg = my_mean(sc)
	if len(sc) > 1:
	    dev = my_std(sc)
	else:
            dev = 0

	#print 'SCALING mean %f std %f'%(avg, dev) 
	f = sc - avg + self.scaling * dev
	#print 'SCALED', f
	# document of choose function
	# http://numeric.scipy.org/numpydoc/numpy-9.html#pgfId-36498
	f=choose(less_equal(f,0.),(f,0.))
	# set the fitness
	#print 'FITNESSES',
	for i in range(len(pop)):
	    pop[i]._fitness = f[i]
	    #print '%.3f'%f[i], pop[i]._fitness
	#print
	#print 'fitness for best score', f
	return pop	
