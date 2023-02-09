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

def readSurface( name ):
    """Read the files 'name'.vertices and 'name'.triangles and returns
    lists of 6-floats for vertices x,y,z,nx,ny,nz and a list of 3-ints
    for triangles"""
    
    import string
    f = open( name+'.vertices' )
    vdata = f.readlines()
    f.close()

    vdata = map( string.split, vdata )
    vdata = map( lambda x: (float(x[0]), float(x[1]), float(x[2]),
                            float(x[3]), float(x[4]), float(x[5])), vdata )

    f = open( name+'.triangles' )
    tdata = f.readlines()
    f.close()

    tdata = map( string.split, tdata )
    tdata = map( lambda x: (int(x[0]), int(x[1]), int(x[2])), tdata )

    return vdata, tdata
