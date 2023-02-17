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

## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

import os
def readSurface( name ):
    """Read the files 'name'.vertices and 'name'.triangles and returns
    lists of 6-floats for vertices x,y,z,nx,ny,nz and a list of 3-ints
    for triangles"""
    
    import string
    f = open( 'DejaVu'+os.sep+'doc'+os.sep+'Tutorial'+os.sep+name+'.vertices' )
    vdata = f.readlines()
    f.close()

    vdata = map( string.split, vdata )
    vdata = map( lambda x: (float(x[0]), float(x[1]), float(x[2]),
                            float(x[3]), float(x[4]), float(x[5])), vdata )

    f = open( 'DejaVu'+os.sep+'doc'+os.sep+'Tutorial'+os.sep+name+'.triangles' )
    tdata = f.readlines()
    f.close()

    tdata = map( string.split, tdata )
    tdata = map( lambda x: (int(x[0]), int(x[1]), int(x[2])), tdata )

    return vdata, tdata

print 'loading the surface'
v, t = readSurface( 'surface' )

# make a numeric array out of the vertices so we can easily separate vertices
# and normals
import numpy.oldnumeric as Numeric
vn = Numeric.array(v)

print 'getting a viewer'
from DejaVu import Viewer
vi = Viewer()

print 'adding the surface to the viewer'
from DejaVu.IndexedPolygons import IndexedPolygons
srf = IndexedPolygons('myFirstSurface', vertices = vn[:,:3], faces = t)
#                      vnormals = vn[:,3:6], faces = t )
vi.AddObject(srf)

#vertices normals
from DejaVu.Polylines import Polylines
pts = vn.__copy__()
vn[:,3:6] = vn[:,:3]+vn[:,3:6]
pts = Numeric.reshape( vn, (-1,2,3) )

p = Polylines('normals', vertices = pts)
vi.AddObject(p)

#faces normals
#from OpenGL import GL
from geomutils.geomalgorithms import  TriangleNormals
vc = vn[:,:3].__copy__()
#nf = GL.glTriangleNormals(vc, t, 'PER_FACE' )
nf = TriangleNormals(vc, t, 'PER_FACE' )

#face centers
from DejaVu.Spheres import Spheres
pts = Numeric.take(vn[:,:3], t)
cg = Numeric.sum(pts, 1)/3.0
s = Spheres('faceCenters', centers=cg, radii=0.1 )
vi.AddObject(s)

pts = Numeric.concatenate( (cg, cg+nf), 1 )
pts  = Numeric.reshape(pts, (-1,2,3))

pf = Polylines('faceNormals', vertices = pts)
vi.AddObject(pf)

