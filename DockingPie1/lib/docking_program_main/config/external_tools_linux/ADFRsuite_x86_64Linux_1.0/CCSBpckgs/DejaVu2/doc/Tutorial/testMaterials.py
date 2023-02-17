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

v = ( (0.,0.,0.),  (1.,0.,0.), 
      (0.,1.,0.),  (1.,1.,0.),
      (0.,2.,0.),  (1.,2.,0.),
      (0.,3.,0.),  (1.,3.,0.),
      (0.,4.,0.),  (1.,4.,0.),
      (0.,5.,0.),  (1.,5.,0.),
      (0.,6.,0.),  (1.,6.,0.))

ind = (range(14),)

RED =   (1., 0., 0.)
GREEN = (0., 1., 0.)
BLUE =  (0., 0., 1.)
col = ( RED, RED, RED, GREEN, GREEN, GREEN, BLUE, BLUE, BLUE,
        RED, GREEN, BLUE, RED, GREEN )

col2 = ( RED, RED, RED, RED, RED, RED, RED,
         GREEN, GREEN, GREEN, GREEN, GREEN, GREEN, GREEN)

from DejaVu.IndexedPolylines import IndexedPolylines
p = IndexedPolylines('testColor', vertices = v, faces = ind, materials = col)

from DejaVu import Viewer
vi = Viewer()
vi.AddObject(p)

p2 = IndexedPolylines('testColor2', vertices = v, faces = ind,
                      materials = col2)
vi.AddObject(p2)

norm = ((1.0, 0., 0.0 ),) * 14
pn = IndexedPolylines('testMaterial', vertices = v, faces = ind,
                      materials = col, vnormals = norm)
vi.AddObject(pn)

pn2col = IndexedPolylines('testMaterial2', vertices = v, faces = ind,
                          materials = col2,
                          vnormals = norm)
vi.AddObject(pn2col)

from DejaVu.Spheres import Spheres
s1 = Spheres('test', centers = v, radii = (0.4,), materials = col)
vi.AddObject(s1)

s2 = Spheres('test', centers = v, radii = (0.4,), materials = col2)
vi.AddObject(s2)

#make OPT=-g
#cp openglutil_num.so /mgl/tools/public_python/1.5.2b1/sgi4DIRIX646/lib/python1.5/site-packages/OpenGL/OpenGL/shared/irix646/openglutil_num.so
