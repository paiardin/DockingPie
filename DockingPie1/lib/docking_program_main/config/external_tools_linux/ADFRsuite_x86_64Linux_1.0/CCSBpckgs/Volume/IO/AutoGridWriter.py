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
from Volume.Grid3D import Grid3D

class WriteAutoGrid:

    def write(self, grid, filename):
        """
        write a Grid3D as an AutoDockGrid
        """
        assert isinstance(grid, Grid3D)

        f = open(filename, 'w')

        ## write header
        ##
        name = grid.header.get('GRID_PARAMETER_FILE', '')
        f.write('GRID_PARAMETER_FILE %s\n'%name)
        name = grid.header.get('GRID_DATA_FILE', '')
        f.write('GRID_DATA_FILE %s\n'%name)
        name = grid.header.get('MACROMOLECULE', '')
        f.write('MACROMOLECULE %s\n'%name)
        spacing = grid.stepSize[0]
        f.write('SPACING %.3f\n'%spacing)
        nx, ny, nz = grid.dimensions
        f.write('NELEMENTS %d %d %d\n'%(nx-1, ny-1, nz-1))
        ox, oy, oz = grid.origin
        center = (ox+(nx/2)*spacing, oy+(ny/2)*spacing, oz+(nz/2)*spacing)
        f.write('CENTER %.3f %.3f %.3f\n'%center)

        ## write the data
        ##
        data = grid.data
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    f.write("%.3f\n"%data[i,j,k])

        f.close()
