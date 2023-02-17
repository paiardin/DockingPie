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
from Volume.Grid3D import Grid3DF

class UHBDReaderASCII:

    def read(self, filename, normalize):
        self.filename = filename

        # open file to read in binary mode
        f = open(filename)
        title = f.readline()
        h = self.header = {'title':title}
        # 
        line = f.readline()
        scale = float(line[:12])
        dum2 = float(line[12:25])
        grdflg, idum2, km, one, km = map( int, line[24:].split() )
        #
        h['grdflg'] =  grdflg
        h['idum2'] = idum2
        h['km'] = km  
        line = f.readline()
        dims = map(int, line[:21].split() )
        scaling = float(line[21:33])
        origin =  (float(line[33:45]), float(line[45:57]), float(line[57:69]))
        h['dims'] =  dims
        h['scaling'] =  scaling
        h['origin'] =  origin

        # skip 2 lines
        line = f.readline()
        line = f.readline()

        # create numpy array to store data
        data = self.data = numpy.zeros( dims, 'f')

        #
        lines = f.readlines()
        f.close()
        
        dimx, dimy, dimz = dims
        dimyz = dimy*dimz
        curline = 0
        curdata = 0
        for i in xrange(dims[0]):
            ii, ij, ik = map( int, lines[curline].split() )
            assert (ii-1)==i
            nblines = dimyz/6
            curline += 1
            for j in xrange(nblines):
                data.flat[curdata:curdata+6] = map( float,
                                                    lines[curline].split() )
                curdata += 6
                curline += 1
            # handle last line
            nbvalleft = dimyz-(nblines*6)
            data.flat[curdata:curdata+nbvalleft] = map( float,
                                                       lines[curline].split())
            curdata += nbvalleft
            curline += 1

        stepSize = h['stepSize'] = (scaling, scaling, scaling)
        volume = Grid3DF(self.data, origin, stepSize, h)
        return volume

    
    def describe(self):
        print "UHBD GridDX file: ", self.filename
        print "nx= ", self.header['dims'][0]
        print "ny= ", self.header['dims'][1]
        print "nz= ", self.header['dims'][2]
        print "xlen= ", self.header['scaling']
        print "ylen= ", self.header['scaling']
        print "zlen= ", self.header['scaling']
        print "min= ", min(self.data.ravel())
        print "max= ", max(self.data.ravel())
        print "origin= ", self.header['origin'] 
        print "stepSize= ", self.header['stepSize'] 


if __name__=='__main__':
    
    from Volume.IO.UHBDGridReader import UHBDReaderASCII
    reader = UHBDReaderASCII()
    vol = reader.read('uhbd_example/mache.50mM.andy.uhbdgrd')
