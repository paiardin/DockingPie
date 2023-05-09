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

#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2016
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/PmvApp/extrude.py,v 1.1.4.1 2017/07/13 20:52:52 annao Exp $
#
# $Id: extrude.py,v 1.1.4.1 2017/07/13 20:52:52 annao Exp $
#
from numpy import dot, identity, ones, zeros


class Extruder:
    """
    class to extrude Path instances decsribing shapes along a 3D path
    """

    def __call__(self, path, pathNV, pathBNV, ptsPerRes, atIndices, extrudeVec,
                 sstype, helixShape, sheetShape, coilShape, foff=0):
        # path is the interpolated smooth set of points (4D)
        # pathNV are the vectors in the 2D path plan and orthogonal to the path
        # pathBNV are the binormal vectors (i.e. orthogonal to plane of the
        #         2D sheet)
        # ptsPreRes indicated how many points in path are assigned to each
        #         resdiues
        # atIndices are the atom indices of the atoms used to compute the path
        # extrudeVec is the vector from the first point in the path to the
        #            next one and is used to decide if the cap should be flipped
        # sstype is a list of secondary structure types (one per path point)
        # *Shape     are 2D shapes used for extrusion for the respective sstype
        # foff       is an offset for numbering faces
        
        #assert isinstance(coilShape, Shape2D)
        #assert isinstance(sheetShape, Shape2D)
        #assert isinstance(helixShape, Shape2D)
        self.extrudeVec = extrudeVec
        mat = identity(3)
        extrudedC = []
        extrudedN = []
        # a lookup providing an atom index for each vertex in the cartoon geom
        vertToAtind = []
        # a dictionnary where the key is a residue index and the values is
        # a list of face indices assoviate with that residue in the cartoon geom
        faces4At = {}
        for i in atIndices:
            faces4At[i] = []
        faces = []
        vnum = 0 # number of vertices so far
        fnum = 0 # number of faces so far
        curSst = None
        curAtind = None
        # list of caps faces containing the verices of the face and their
        # indices and where to flip the face or not
        caps = []
        aindex = -1 # residue number
        for i, pt in enumerate(path):
            #rm = min(i/ptsPerRes, len(faces4At)-1)
            atind = atIndices[i]

            if curAtind != atind: # new residue
                aindex += 1
                curAtind = atind
                # find type of ss element and select shape
                try: # try because the last smooth point ends up for nres+1
                    sst = sstype[aindex]
                except IndexError:
                    pass

                if curSst!=sst: # new SS ELEMENT
                    ssElemChanged = True
                    curSst = sst
                    if i > 0:
                        # add ending cap for previous shape
                        nbExtV = len(extrudedC)
                        caps.append( (shape, extrudedC[-shapeLen::stride],
                                      range(nbExtV-shapeLen,nbExtV,stride),
                                      not flipFaces, atIndices[aindex-1]) )
                    # select new shape
                    if sst=='H' or sst=='G' or sst=='I': # helix
                        shape = helixShape
                    elif sst=='B' or sst=='E':
                        shape = sheetShape
                    else:
                        shape = coilShape
                    vertices = shape.contpts
                    normals = shape.contnorm
                    shapeLen = len(vertices)
                    if shape.vertDup:
                        stride = 2
                    else:
                        stride = 1

                    if i>0:
                        # transform new shape at previous path point
                        mat[0, :3] = pathNV[i-1]
                        mat[1, :3] = pathBNV[i-1]
                        mat[2, :3] = path[i-1]
                        tc = dot(vertices, mat)
                        extrudedC.extend(tc[:,:3])
                        vertToAtind.extend([atind]*len(tc))
                        # create cap of first face
                        # FIXME only works for rectangle
                        if normals is not None:
                            tcn = dot(normals, mat)
                        extrudedN.extend(tcn[:,:3])
                        vnum += len(tc)

                        # add cap for begining new shape
                        nbExtV = len(extrudedC)
                        caps.append( (shape, tc[::stride],
                                      range(nbExtV-shapeLen,nbExtV, stride), flipFaces, atind) )
                else:
                    ssElemChanged = False

            # transform shape to reference frame on the path
            mat[0, :3] = pathNV[i]
            mat[1, :3] = pathBNV[i]
            mat[2, :3] = path[i]
            tc = dot(vertices, mat)
            extrudedC.extend(dot(vertices, mat))
            vertToAtind.extend([atind]*len(tc))
            if normals is not None:
                tcn = dot(normals, mat)
                extrudedN.extend(tcn)#[:,:3])
            if i==0:
                # check if normal the face is along extrusion direction
                if shape.vertDup:
                    x1,y1,z1 = tc[0]
                    x2,y2,z2 = tc[2]
                    x3,y3,z3 = tc[4]
                else:
                    x1,y1,z1 = tc[0]
                    x2,y2,z2 = tc[1]
                    x3,y3,z3 = tc[2]
                A = [x2-x1, y2-y1, z2-z1]
                B = [x3-x1, y3-y1, z3-z1]
                cross = [A[1]*B[2]-B[1]*A[2],
                         A[2]*B[0]-B[2]*A[0],
                         A[0]*B[1]-B[0]*A[1]]
                dotp = cross[0]*self.extrudeVec[0]+\
                       cross[1]*self.extrudeVec[1]+\
                       cross[2]*extrudeVec[2]
                if dotp>0:
                    flipFaces = False
                else:
                    flipFaces = True

                # add cap for begining new shape
                nbExtV = len(extrudedC)
                caps.append( (shape, tc[::stride], range(nbExtV-shapeLen,nbExtV,stride),
                              flipFaces, atIndices[0]) )

            if i > 0:
                v0p = vnum-shapeLen
                for fvnum in range(0, shapeLen-stride, stride):
                    if flipFaces:
                        faces.append((vnum+fvnum, vnum+fvnum+stride,
                                      v0p+fvnum+stride, v0p+fvnum))
                    else:
                        faces.append((v0p+fvnum, v0p+fvnum+stride,
                                      vnum+fvnum+stride, vnum+fvnum))
                    faces4At[atind].append(fnum+foff)
                    fnum+=1

                # add last face
                if flipFaces:
                    faces.append((vnum+fvnum+stride, vnum, v0p, v0p+fvnum+stride))
                else:
                    faces.append((v0p+fvnum+stride, v0p, vnum,vnum+fvnum+stride,))
                faces4At[atind].append(fnum+foff)
                fnum+=1
                    
            vnum += shapeLen
        # add last cap
        caps.append( (shape, tc[::stride], range(nbExtV-shapeLen,nbExtV,stride), not flipFaces, atind) )

        vnumNoCaps = vnum
        # add the faces for the caps
        for shape, verts ,vind, flip, atind in caps:
            v, n, f = shape.capGeom(verts, vind, vnum, flip)
            vnum+= len(v)
            extrudedC.extend(v)
            vertToAtind.extend([atind]*len(v))
            extrudedN.extend([n]*len(v))
            faces4At[atind].extend(range(fnum+foff, fnum+foff+len(f)))
            faces.extend(f)
            fnum += len(f)
        return extrudedC, faces, extrudedN, vertToAtind, faces4At
                    
            
