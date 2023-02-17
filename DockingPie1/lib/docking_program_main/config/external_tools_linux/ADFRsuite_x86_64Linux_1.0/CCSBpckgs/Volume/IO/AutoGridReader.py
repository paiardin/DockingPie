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

class ReadAutoGrid:
    """Read a AutoGrid map file"""
    def read(self, filename, normalize):
        """ Read from AutoGrid map file"""   

        self.SPACING=1.0
        self.CENTER=(0.,0.,0.)
     
        self.header = {'title': 'AutoGrid from %s'%filename}
        
        # Pmv/Grid.py, Class Grid, function ReadAutoGridMap()        
        f = open(filename, 'r')
        try:
            GRID_PARAMETER_FILE = f.readline().split()[1]
        except:
            GRID_PARAMETER_FILE = ''
        self.header['GRID_PARAMETER_FILE'] = GRID_PARAMETER_FILE
        GRID_DATA_FILE = f.readline().split()[1]
        self.header['GRID_DATA_FILE'] = GRID_DATA_FILE
        MACROMOLECULE = f.readline().split()[1]
        self.header['MACROMOLECULE'] = MACROMOLECULE
        
        # spacing
        SPACING = float(f.readline().split()[1])
    
        # number of points and center
        (nx,ny,nz) = f.readline().split()[1:4]
        NELEMENTS = (nx,ny,nz) = (int(nx)+1, int(ny)+1, int(nz)+1)
        (cx,cy,cz) = f.readline().split()[1:4]
        CENTER = ( float(cx),float(cy), float(cz))
        
    
        # read grid points
        points = map( lambda x: float(x), f.readlines())
       
        # data read as z,y,z, swapaxes to make the data x,y,z
        TMPGRIDS = numpy.swapaxes(numpy.reshape( points,(nz,ny,nx)), 0, 2)
        GRIDS = numpy.array(TMPGRIDS).astype('f')
        f.close()
        self.data = GRIDS
        #print "shape***:",self.data.shape
        #print "origin***:",CENTER
        
        origin = (CENTER[0]-(nx/2)*SPACING, CENTER[1]-(ny/2)*SPACING,
                  CENTER[2]-(nz/2)*SPACING)
        stepSize = (SPACING,SPACING,SPACING)
        
        #def __init__(self, data, origin, stepSize, header):
        grid = Grid3DF(self.data, origin, stepSize, self.header)
        
        #print "**", grid.dimensions
       
        return grid
       
