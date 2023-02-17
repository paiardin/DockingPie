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

import sys; sys.path.insert(0, '.')

#Example Two: Changing materials
#Start this example by making a collection of points, v:
v = 	( (0.,0.,0.), (1.,0.,0.),
	(0.,1.,0.), (1.,1.,0.), 
	(0.,2.,0.), (1.,2.,0.),
	(0.,3.,0.), (1.,3.,0.),
	(0.,4.,0.), (1.,4.,0.),
	(0.,5.,0.), (1.,5.,0.),
	(0.,6.,0.), (1.,6.,0.))

#defining some colors:
RED =(1., 0., 0.)
GREEN = (0., 1., 0.)
BLUE =(0., 0., 1.)

#and collections of colors:
col = ( RED, RED, RED, GREEN, GREEN, GREEN, BLUE, BLUE, BLUE, RED, GREEN, BLUE, RED, GREEN )

col2 = ( RED, RED, RED, RED, RED, RED, RED, GREEN, GREEN, GREEN,
GREEN, GREEN, GREEN, GREEN)

#Define a list to specify the faces of the lines we make later:
ind = (range(14),)

#Start up a viewer :

from DejaVu import Viewer
vi = Viewer()

#and make a line:

from DejaVu.IndexedPolylines import IndexedPolylines
p = IndexedPolylines('testColor', vertices = v, faces = ind,
materials = col)

#and add it to the viewer:
vi.AddObject(p)

#and make another line:
p2 = IndexedPolylines('testColor2', vertices = v, faces = ind,materials
= col2)

#and add it to the viewer:
vi.AddObject(p2)

#With these two objects in the viewer, try changing the current object
#and transforming it.

#Add another line:
norm = ((1.0, 0., 0.0 ),) * 14
pn = IndexedPolylines('testMaterial', vertices = v, faces = ind,
materials = col, vnormals = norm)
vi.AddObject(pn)

#Add another line:
pn2col = IndexedPolylines('testMaterial2', vertices = v, faces
= ind, materials = col2, vnormals = norm)

vi.AddObject(pn2col)

#Finally, try making some rows of spheres colored differently:
from DejaVu.Spheres import Spheres
s1 = Spheres('test', centers = v, radii = (0.4,), materials = col,quality=15)
vi.AddObject(s1)

s2 = Spheres('test', centers = v, radii = (0.4,), materials = col2)
vi.AddObject(s2)

print "end of Example2"