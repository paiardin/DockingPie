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
#Example 1: Getting started

from DejaVu import Viewer

#instantiate our first viewer:
MyViewer = Viewer()

#Now let's display a bunch of lines. To do this we need 3D coordinates:
coords = [ [0.,0.,0.], [1.,1.,0.], [0.,1.,0.], [1.,0.,0.] ]

#an array of indices telling which vertices are to be connected by lines:
#(each entry represents a line between two points whose indices in the coords array are given)
indices = [[0,1], [2,3],[0,2],[2,1],[1,3],[3,0]]

#Alternatively, the lines could be specified where each list represents
#a line through the points whose indices are given. The value -1 terminates
#the line. The first two entries draw the diagonals and the third the box
#itself: indices = [[0,1, -1, 0,0], [2,3,-1,0,0], [0,2,1,3,0]]

#an optional array of materials. Here the tuples represent RGBA values:
materials = ( (1.,1.,1.), (1.,0.,0.), (0.,1.,0.), (0.,0.,1.) )

#We create a geometry object of type IndexedPolylines:
from DejaVu.IndexedPolylines import IndexedPolylines

cross = IndexedPolylines('MyFirstObject')

# add the vertices, lines and materials to that object:
cross.Set(vertices=coords, faces=indices, materials = materials,
          inheritMaterial = 0 )


#add the geometry to the viewer:
MyViewer.AddObject(cross)
MyViewer.Redraw()

#Now the object listbox should have one more line ~MyFirstObject.
#The ~ is used to visualize the hierarchy. By default, AddObject makes the
#new object the child of root.&nbsp; This can be changed programmatically
#In the camera you should see a square with four vertices colored 
#white, red, green and blue. You can use the mouse to transform this object:
###????###

# next, add a clipping plane:
#NOTE:You must make ~MyFirstObject the current object by selecting it in the
#ViewerGUI or by clicking and you should increase the size of "cross" with Shift-Middle Mouse Button 
#before starting this section. ALSO you must select Clip on the Properties panel
#and turn the row1 "on" button so that the clipping plane will interact with
# the current object and the "disp" display button if you want to see the 
#clipping plane. 
cl = MyViewer.clipP[0]

#Make the clippingplane visible
cl.Set(visible=1)
# translate it to the right
cl.ConcatTranslation( (0.5, 0., 0.) )

# activate this cliping plane for our geometry
cross.AddClipPlane(cl)

print "end of Example1"