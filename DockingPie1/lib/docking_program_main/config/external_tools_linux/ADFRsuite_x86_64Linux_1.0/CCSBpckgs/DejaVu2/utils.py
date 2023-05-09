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

def RemoveDuplicatedVertices(vertices, faces, vnorms = None):
    """Remove duplicated vertices and re-index the polygonal faces such that
they share vertices
"""
    vl = {}
    vrl = {}
    nvert = 0
    vertList = []
    normList = []
    
        
    for i,v in enumerate(vertices):
        key = '%f%f%f'%tuple(v)
        if not vl.has_key(key):
            vl[key] = nvert
            vrl[i] = nvert
            nvert +=1
            vertList.append(v)
            if vnorms is not None:
                normList.append(vnorms[i])
        else:
            nind = vl[key]
            vrl[i] = nind
            if vnorms is not None:
                vn1 = normList[nind]
                vn2 = vnorms[i]
                normList[nind] = [(vn1[0] +vn2[0])/2,
                                   (vn1[1] +vn2[1])/2,
                                   (vn1[2] +vn2[2])/2 ]

    faceList = []
    for f in faces:
        faceList.append( map( lambda x, l=vrl: vrl[x],  f ) )
    if vnorms is not None:
        return vertList, faceList, normList
    return vertList, faceList

