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

from DejaVu.IndexedPolygons import IndexedPolygons
from Volume.Grid3D import Grid3D

class ClipMeshWithMask:
    """Clip method of this class takes a mesh i.e. IndexedPolgons and
selects all vertices which fall onto voxel witha  true value in a mask grid.
It returns a new IndexedPolygons geometry with the triangles for which
3 vertices are selected.
"""
    def __init__(self):
        pass

    def clip(self, mesh, grid):
        assert isinstance(mesh, IndexedPolygons)
        assert isinstance(grid, Grid3D)

        origin = grid.getOriginReal()
        stepSize = grid.getStepSizeReal()
        dx, dy, dz = grid.dimensions

        vertices = mesh.vertexSet.vertices.array
        triangles = mesh.faceSet.faces.array

        # compute the voxel on which each vertex falls
        # array of indiced into grid for the vertices
        vertInd = ((vertices-origin)/stepSize).astype('i')

        # select the vertices on voxels that have a value True
        selVert = []
        vertEquiv = {}
        numVertSel = 0
        nvert = 0
        data = grid.data
        for i,j,k in vertInd:
            if i>=0 and i<dx:
                if j>=0 and j<dy:
                    if k>=0 and k<dz:
                        if data[i,j,k]:
                            selVert.append( vertices[nvert] )
                            vertEquiv[nvert] = numVertSel
                            numVertSel += 1
            nvert += 1
            
        # build a set of faces for which some vertices are selected
        # and keep only selected vertices
        
        selFaces = []
        for i,j,k in triangles:
            nbvs = 0
            v1 = vertEquiv.get(i, None)
            if v1: nbvs +=1
            v2 = vertEquiv.get(j, None)
            if v2: nbvs +=1
            v3 = vertEquiv.get(k, None)
            if v3: nbvs +=1
            if nbvs == 3:
                selFaces.append( (v1,v2,v3) )

        clippedGeom = IndexedPolygons(mesh.name+'_clipped', vertices=selVert,
                                      faces=selFaces)
        return clippedGeom
