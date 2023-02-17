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

## def uniq(l, func=None):
##     """Return a new list with duplicate items removed."""

##     l2 = l[:]	# make a copy
##     d = {}
##     def add_to_dict(value,d=d):
## 	d[`value`] = value
##     map(add_to_dict,l2)
##     l3 = d.values()
##     if len(l2)==len(l3): return(l2)
##     if func: l3.sort(func)
##     else: l3.sort()
##     return l3

def uniq(alist):    # Fastest order preserving
    set = {}
    return [set.setdefault(e,e) for e in alist if e not in set]
 
def uniq3(alist):    # Fastest without order preserving
    set = {}
    map(set.__setitem__, alist, [])
    return set.keys()

"""
from mglutil.util.uniq import uniq, uniq2, uniq3
import time
a=range(100)
b=range(10)
c=a+b

t1=time.time()
for i in range(5000):  x=uniq(c)
print time.time()-t1

t1=time.time()
for i in range(5000):  x=uniq2(c)
print time.time()-t1

t1=time.time()
for i in range(5000):  x=uniq3(c)
print time.time()-t1


>>>
0.865363121033
0.463307857513
0.260641098022
"""
