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

from DejaVu2.Polylines import Polylines

class NormalsViewer:
    """Object that take a DejaVu2 geometry and a viewer and displays the
    geometry's normals in the viewer"""

    def __init__(self, geom, viewer):
        self.geom = geom
        self.normalsGeom = Polylines('normals_for_'+geom.name)
        self.viewer = viewer
        viewer.AddObject(self.normalsGeom, parent=geom)
        self.update()
        
    def update(self):
        vertices = self.geom.getVertices()
        normals = self.geom.getVNormals()
        pts = numpy.concatenate( (vertices, vertices+normals), 1)
        pts = numpy.reshape( pts, (len(vertices),2,-1) ).astype('f')
        self.normalsGeom.Set(vertices = pts) 
        self.viewer.Redraw()
       
