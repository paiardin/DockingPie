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

from struct import pack
from Volume.Grid3D import Grid3DUC, Grid3DSI, Grid3DF

class VolumeWriterBase:

    pass

class WriteCCP4(VolumeWriterBase):
    """Write a CCP4 file from a Numeric array."""
    
    def write(self, filename, grid3D):
        # open file to read in binary mode
        self.filename = filename
        f = open(filename, 'wb')

        # read the 1024 bytes of the header
        header = grid3D.header
        array = grid3D.data

        # always uss mode=2 (FLOAT values)
        s = pack("4i", array.shape[2], array.shape[1], array.shape[0], 2)
        
        nstart = [0,0,0]
        for i in (0,1,2):
            nstart[i] = round(grid3D.origin[i] / grid3D.stepSize[i])

        if header.has_key('nx'):
            ncstart = round(header['nx']*grid3D.origin[2])
            nrstart = round(header['ny']*grid3D.origin[1])
            nsstart = round(header['nz']*grid3D.origin[0])
            nx = header['nx']
            ny = header['ny']
            nz = header['nz']
        else:
            ncstart = 0
            nrstart = 0
            nsstart = 0
            nx = array.shape[2]
            ny = array.shape[1]
            nz = array.shape[0]

        s += pack("6i", ncstart, nrstart, nsstart, nx, ny, nz)

        acell = header.get('acell', grid3D.stepSize[2]*(array.shape[2]-1))
        bcell = header.get('bcell', grid3D.stepSize[1]*(array.shape[1]-1))
        ccell = header.get('ccell', grid3D.stepSize[0]*(array.shape[0]-1))
        alpha = header.get('alpha', 90.0)
        beta = header.get('beta', 90.0)
        gamma = header.get('gamma', 90.0)
        s += pack("6f", acell, bcell, ccell, alpha, beta, gamma)

        # we always use c-style data ordering
        mapc = 3
        mapr = 2
        maps = 1
        s += pack("3i", mapc, mapr, maps)

        # FIXME we shoudl get real values here if they are not in the header
        amin = header.get('amin', 0.0)
        amax = header.get('amax', 0.0)
        amean = header.get('amean', 0.0)
        s += pack("3f", amin, amax, amean)

        ispg = header.get('ispg', 0)
        nsymbt = header.get('nsymbt', 0)
        lskflg = header.get('lskflg', 0)
        
        s += pack("3i", ispg, nsymbt,lskflg)

        skwmat = header.get('skwmat', (0.0,)*9)
        skwtrn = header.get('skwtrn', (0.0,)*3)

        s += apply( pack, ("12f",) + skwmat + skwtrn )
        
        future_words = header.get('future_words', (0,)*15)
        s += apply( pack, ("15i",) + future_words )

        # we don't use this, but put it here for completion
        s += pack("4c", ' ', ' ', ' ', ' ')

        # Machine stamp
        #machst = unpack(sw+"4c", data[212:216])
        s += pack("4c",  ' ', ' ', ' ', ' ')

        arms = header.get('arms', 0.0)
        s += pack("1f", arms)

        s += pack("i", 0)
        s += apply( pack, ("800c",) + (' ',)*800)
        assert len(s)==1024

        #append an according length of space for the symmetry ops
        s += apply( pack, ("%dc"%nsymbt,) + (' ',)*nsymbt)
        
        f.write('%s'%s)

        if array.flags.contiguous:
            #print 'writing contiguous'
            f.write('%s'% apply ( pack, ("%df"%(array.shape[0]*array.shape[1]*
                                                array.shape[2]),
                                         ) + tuple(array.astype('f').ravel())
                                  )
                    )
        else:
            #print 'writing NON contiguous'
            for i in range(array.shape[0]):
                for j in range(array.shape[1]):
                    for k in range(array.shape[2]):
                        f.write('%f'%array[i][j][k])
        f.close()
import numpy

class WriteRawiv:
    """ Writes volumetric data to .rawiv file format. """

##     def __init__(self):
##         self.file = None
##         self.data = None

    def write_file(self, file, data, dims, stepSize=0, orig = (0., 0., 0.)):
        #data should be a numeric array of either uchar type(Numeric.UnsignedInt8)
        # or float (Numeric.Float32, Numeric.Float16)
        assert len(dims) == 3
        nx, ny, nz = dims

        arraytype = data.dtype.char
        if arraytype == numpy.uint8:
            packType = 'B'
        elif arraytype in (numpy.float32, numpy.float16):
            packType = 'f'
        else:
            print "Error in write_file() , typecode %s not supported"%(arraytype,)
            return

        if not stepSize:
            stepSizeX = (maxX - minX)/(nx-1.)
            stepSizeY = (maxY - minY)/(ny-1.)
            stepSizeZ = (maxZ - minZ)/(nz-1.)
        else:
            try:
                l = len(stepSize)
                assert l == 3
                stepSizeX = stepSize[0]
                stepSizeY = stepSize[1]
                stepSizeZ = stepSize[2]
            except:
                stepSizeX = stepSizeY = stepSizeZ = stepSize

        # coords of the 1st voxel
        minX, minY, minZ = orig

        #coords of the last voxel
        maxX = orig[0] + (nx-1)*stepSizeX
        maxY = orig[1] + (ny-1)*stepSizeY
        maxZ = orig[2] + (nz-1)*stepSizeZ

        # number of vertices in the grid
        size = nx*ny*nz
        # number of cells in the grid
        size1 = (nx-1)*(ny-1)*(nz-1)
        # step size - spacing between one vertex and the next along each dimension.
        st=(minX, minY, minZ, maxX, maxY, maxZ, size, size1, nx, ny, nz,
            orig[0], orig[1], orig[2], stepSizeX, stepSizeY, stepSizeZ)
        
        #open file for writting data:
        of = open(file,"wb")
        #write header:
        of.write(apply(pack, ('>6f5I6f',)+st))
        fmt = ">%d%s"%(size, packType)
        #write data:
        of.write( apply(pack, (fmt,)+tuple(data.ravel())))
        of.close()
    

if __name__=='__main__':
    from volReaders import ReadCCP4
