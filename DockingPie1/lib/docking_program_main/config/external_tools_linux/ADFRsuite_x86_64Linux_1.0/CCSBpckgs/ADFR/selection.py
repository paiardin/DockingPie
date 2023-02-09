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
# $Header: /mnt/raid/services/cvs/ADFR/selection.py,v 1.3 2016/12/07 00:38:32 sanner Exp $
#
# $Id: selection.py,v 1.3 2016/12/07 00:38:32 sanner Exp $
#
#genetic algorithm selection routines
#based on galib. 
#exception - these classes only work on the scaled fitness

from ga_util import choice

import numpy

from random import random
rand=random
from math import floor

class selector:
    def update(self,pop):
	pass
    def select(self,pop):
	raise GAError, 'selector.select() must be overridden'
    def clear(self):
		pass

class uniform_selector(selector):
    def select(self,pop,cnt = 1):
	if cnt == 1:
	    return choice(pop)
	res = []
	for i in range(cnt):
	    res.append(choice(pop))
	return res

class rank_selector(selector):
    def select(self,pop,cnt = 1):
	pop.sort()
	studliest = pop[0].fitness()
	tied_for_first = filter(lambda x,y=studliest: x.fitness()==y,pop)
	if cnt == 1:
	    return choice(tied_for_first)
	res = []
	for i in range(cnt):
	    res.append(choice(tied_for_first))
	return res

#scores must all be positive		
class roulette_selector(selector):
    def update(self,pop):
	self.pop = pop[:]
	sz = len(pop)
	if not sz:
	    raise GAError, 'srs_selector - the pop size is 0!'
	f =self.pop.fitnesses()
	f_max = max(f); f_min = min(f)
	if not ( (f_max >= 0 and f_min >= 0) or \
		 (f_max <= 0 and f_min <= 0)):
	    raise GAError, 'srs_selector requires all fitnesses values to be either strictly positive or strictly negative'
	if f_max == f_min: f = ones(shape(f),typecode = Float32)
	self.dart_board = add.accumulate(f / sum(f))
	
    def select(self,pop,cnt = 1):
	returns = []
	for i in range(cnt):
	    dart = rand()
	    idx = 0
	    while dart > self.dart_board[idx]:
		idx += 1
	    returns.append(self.pop[idx])
	if cnt == 1:
	    return returns[0]
	else:
	    return returns
	
    def clear(self): 
	del self.pop

from bisect import bisect

#scores must all be positive		
class srs_selector(selector):


    def update(self, pop):
        # build self.choices, a list of indices in the population used to select for cross and mutations
	sz = len(pop)
	if not sz:
	    raise GAError, 'srs_selector - the pop size is 0!'
        
	f = pop.fitnesses()
	f_max = numpy.max(f); f_min = numpy.min(f)
	if not ( (f_max >= 0. and f_min >= 0.) or 
		 (f_max <= 0. and f_min <= 0.)):
	    raise GAError, 'srs_selector requires all fitnesses values to be either strictly positive or strictly negative - min %f, max %f' %(f_min,f_max)

	f_avg = numpy.sum(f)/sz
	if abs(f_avg) < 1.e-10:
	    e = numpy.ones(numpy.shape(f), 'f')
	else:
	    e = f/f_avg
	self.expected_value = e
	garauntee, chance = divmod(e,1.)
	choices = []
        choiceDict = {}
	for i in xrange(sz):
	    choices = choices + [i] * int(garauntee[i])
            if choiceDict.has_key(i):
                choiceDict[i] += int(garauntee[i])
            else:
                choiceDict[i] = int(garauntee[i])
                
	# now deal with the remainder
        # create a cumulative probability distribution (dart_board)
        # Then generate a random number between 0 and 1 and do a binary search
	sum_tmp = numpy.sum(chance)
	if sum_tmp !=0.0: 
	    for i in range(len(choices),sz):
                cdf = numpy.add.accumulate(chance / sum_tmp)
                idx = bisect(cdf,rand())
	    ## dart_board = add.accumulate(chance / sum_tmp)
	    ## for i in range(len(choices),sz):
	    ## 	dart = rand()
	    ## 	idx = 0
	    ##     while dart > dart_board[idx]: idx = idx + 1
	    ##     choices.append(idx)
                if choiceDict.has_key(idx):
                    choiceDict[idx] += 1
                else:
                    choiceDict[idx] = 1
                
        #if len(choiceDict) < 4:
        #    import pdb
        #    pdb.set_trace()
	self.choices = choices
        #print '  MAX CHOICE', choiceDict


    def select(self, pop, cnt=1):
        """
        Pick cnt individuals out of pop by randomly picking indices in self.choice
        """
        ## MS new implementation
        self._selected = []
        res = []
        lenChoices = len(self.choices)
        for i in range(cnt):
            ind = self.choices[int(floor(random()*lenChoices))]
            # once selected it is removed from the list
            #ind = self.choices.pop([int(floor(random()*lenChoices))])
            self._selected.append(ind)
            res.append(pop[ind])
        if cnt == 1:return res[0]
        else: return res
       
        ## res = []
        ## inds = []
        ## for i in range(cnt):
        ##     ind = choice(self.choices)
        ##     inds.append(ind)
        ##     res.append(pop[ind])

        ## #print "SELECTED", inds
        ## if cnt == 1:
        ##     return res[0]
        ## else:
        ##     return res


    def clear(self): 
        if hasattr(self,'choices'):
            del self.choices		
