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

import numpy
from struct import unpack
from  Volume.Grid3D import Grid3DF

"""
The delphi manual specifies the format as (unformatted binary):
character*20
character *10, character *60
real*4: the actual array with the potential map
character*16
real*4: scale (gridpoints per angstrom) and midpoint of the grid
        (in realspace coordinates)

Example: potential.phi
    The map I send you has a grid of 85 in all directions
    The grid centre (43,43,43 as in (85+1)/2)
    The 'midpoint', is at 13.50739       40.60874       85.72395
    The array contains the grid point coordinates (x,y,z) plus the
    electrosatic potential (the 4th number).
"""

class DelphiReaderBin:
    """Read and DelPhi file and return a Grid3D object"""
    
    def read(self, filename, normalize):
        
        sw = "" # used in format to swap bytes if necessary
        self.filename = filename
        
        # open file to read in binary mode
        f = open(filename, 'rb')

        # read the 1024 bytes of the header
        data = f.read(20) # character*20
        data = f.read(10) # character*10
        data = f.read(80) # character*60

        data = f.read() # rest of file
        f.close()

        nbsteps = unpack("i", data[-8:-4])[0]
        scale, midx, midy, midz = unpack("4f", data[-24:-8])
        endstr =  unpack("16c", data[-48:-32])

        values = unpack("%df"%nbsteps**3, data[:-56])

        h = self.header = {'center':[midx, midy, midz], 'scale':scale,
                           'nbSteps':nbsteps}
        
        stepSize = [1./scale, 1./scale, 1./scale]
        step2 = nbsteps/2.
        origin = [ midx-(step2*stepSize[0]),
                   midy-(step2*stepSize[1]),
                   midz-(step2*stepSize[2]) ]
        dataArray = numpy.array(values, 'f')
        dataArray.shape = (nbsteps, nbsteps, nbsteps)
        g = Grid3DF(dataArray , origin, stepSize, h)
        return g


    def describe(self):
        print "CCP4 file: ", self.filename
        print "center   : ", self.header['center']
        print "nbSteps  : ", self.header['nbSteps'] 
        print "scale    : ", self.header['scale']

