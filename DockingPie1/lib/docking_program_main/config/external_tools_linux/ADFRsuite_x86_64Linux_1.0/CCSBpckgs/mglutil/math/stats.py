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

import warnings

def stats(values):
    """returns the mimn, max, mean and standard deviation of a list of values"""

    warnings.warn(
        "\n\nWARNING!! This function has been deprecated!!\n \
        Use the stats in Volume/Grid3D.py\n",DeprecationWarning,2)
    
    npts = len(values)
    if npts:
        from math import sqrt
        sum = 0.0
        sumsq = 0.0
        mini = maxi = values[0]
        for v in values:
            sum += v
            sumsq += float(v)*float(v)
            if v<mini:
                mini = v
            if v>maxi:
                maxi = v
        mean = float(sum)/npts
        stdev = sqrt(( sumsq - (sum*sum/float(npts)))/(npts-1))
        return mini, maxi, mean, stdev

    else:
      return (0., 0., 1., 1.)
   
