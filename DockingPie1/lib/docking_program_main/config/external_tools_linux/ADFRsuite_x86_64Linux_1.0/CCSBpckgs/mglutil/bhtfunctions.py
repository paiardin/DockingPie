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
from bhtree import bhtreelib
#
# $Header: /mnt/raid/services/cvs/python/packages/share1.5/mglutil/bhtfunctions.py,v 1.9.18.1 2017/07/26 22:35:40 annao Exp $
#
# $Id: bhtfunctions.py,v 1.9.18.1 2017/07/26 22:35:40 annao Exp $
#

# ClosePointsDist2: result and dist are empty arrays large enough to contain
# all the points expected to be found. To be safe, they should be the same
# size as the list of coordinates. This function then puts in the results
# array the indices of the close points in the list (as supplied in the array
# ids): the dist array contains the corresponding distances.

BHTree_CUT = 40.0
def findClosestAtoms(obj_verts, atom_coords,
                     cutoff_from=3.5, cutoff_to=BHTree_CUT,
                     instanceMatrices=None):
    """For every vertex in a given set of vertices finds the closest atom.
    Returns an array of corresponding atom indices. Based on bhtree
    library. """
    if not len(obj_verts):
        return []
    natoms = len(atom_coords)
    if instanceMatrices:
        coordv = numpy.ones(natoms *4, "f")
        coordv.shape = (natoms, 4)
        coordv[:,:3] = atom_coords[:]
        new_coords = []
        for m in instanceMatrices:
            new_coords.append(numpy.dot(coordv, \
                               numpy.transpose(m))[:, :3])
        atom_coords = numpy.concatenate(new_coords)
    #print "len atom_coords: ", len(atom_coords)
    bht = bhtreelib.BHtree( atom_coords, None, 10)
    cl_atoms = []
    mdist = cutoff_from
    print "** Bind Geometry to Molecule Info: **"
    print "** looking for closest atoms (cutoff range: %2f-%2f)...   **"%(cutoff_from, cutoff_to)
    cl_atoms = bht.closestPointsArray(obj_verts, mdist)
    while len(cl_atoms) == 0 and mdist <cutoff_to:
        print "** ... no closest atoms found for cutoff = %2f; continue looking ... **"%mdist
        mdist=mdist+0.2
        cl_atoms = bht.closestPointsArray(obj_verts, mdist)
        #print "mdist: ", mdist, "  len cl_atoms: ", len(cl_atoms)
    print "**... done. %d closest atoms found within distance: %2f **"%(len(cl_atoms) , mdist)
    if instanceMatrices:
        if cl_atoms:
            return [x%natoms for x in cl_atoms]
    return cl_atoms


def findNearestAtoms(mol,vertices, **kw):
    """None <- color(mol,vertices2,**kw)
    mol:           reference molecule
    vertices:      list of lists(coordinates): the first three items in each list
                   must be coordinates x,y,z of a point.
                   
    atomIndices is the index of the nearest atom to the vertex, such that
    mol.allAtoms[atomIndices[x]] is the nearest atom to vertices[x]
    vertexIndices is the list of nearest vertices to an atom, such that
    vertexIndices[x] = [vertex1,vertex2,...] are the vertices associated with
    mol.allAtoms[x]
    """
    
    coords = mol.allAtoms.coords
    if not hasattr(mol,'bhtree'):
        print "Building bhtree for ",mol
        ids = numpy.arange(len(coords)).astype('i')
        bhtree = bhtreelib.TBHTree(coords,ids,10,10,9999.0)
        mol.bhtree = bhtree
    
    vertexIndices={}
    atomIndices={}
    for x in range(len(coords)):
        vertexIndices[x+1]=[]

    cutoff=5.
    for x in range(len(vertices)):
        xyz = vertices[x]
        result = numpy.zeros( (len(vertices),) ).astype('i')
        dist = numpy.zeros( (len(vertices),) ).astype('f')
        nb2 = mol.bhtree.ClosePointsDist2(tuple(xyz[:3]), cutoff, result, dist )
        while nb2==0:
            cutoff = cutoff+5.
            nb2 = mol.bhtree.ClosePointsDist2(tuple(xyz[:3]), cutoff, result, dist )
        result = result[:nb2]
        dist = dist[:nb2]
        idx = dist.tolist().index(min(dist))
        atnum = result[idx]+1
        atomIndices[x]=atnum
        vertexIndices[atnum].append(x)
        
    return atomIndices,vertexIndices








