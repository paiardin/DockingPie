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
## Copyright (c) MGL TSRI 2016
##
################################################################################

from math import sqrt

def triangleArea(p1, p2, p3):
    """Compute the surface area of a triangle.
"""
    x1,y1,z1 = p1
    x2,y2,z2 = p2
    x3,y3,z3 = p3
    dx, dy, dz = x1-x2, y1-y2, z1-z2
    a = sqrt( dx*dx + dy*dy + dz*dz )
    dx, dy, dz = x2-x3, y2-y3, z2-z3
    b = sqrt( dx*dx + dy*dy + dz*dz )
    dx, dy, dz = x1-x3, y1-y3, z1-z3
    c = sqrt( dx*dx + dy*dy + dz*dz )
    s = .5*(a+b+c)
    area = s*(s-a)*(s-b)*(s-c)
    if area <= 0.:
#        print "area = %f for triangles: " % area, p1, p2, p3
        return 0.
    return sqrt(area)

def meshVolume(verts, norm, tri):
    """Compute the Volume of a mesh specified by vertices, their normals, and
indices of triangular faces
"""
    # TEST
    zeronorms = []
    for i, n in enumerate(norm):
        #if n == [0., 0., 0.] or n == (0., 0., 0.):
	if n[0] == 0 and n[1] == 0 and n[2] == 0:
	    #print "normal %d is zero!" % i, n
	    zeronorms.append(i)
    #print "in meshVolume, zeronorms length: ", len(zeronorms), "normals length:", len(norm)
    # Initialize
    volSum = 0.0
    oneThird = 1./3.
    
    # Compute face normals
    trinorm = []
    for t in tri:
        n1 = norm[t[0]]
        n2 = norm[t[1]]
        n3 = norm[t[2]]
	tn = [ (n1[0]+n2[0]+n3[0])*oneThird,
               (n1[1]+n2[1]+n3[1])*oneThird,
               (n1[2]+n2[2]+n3[2])*oneThird ]
	trinorm.append(tn)
    # print trinorm	# TEST

    # Compute volume
    for t,tn in zip(tri, trinorm):
        s1 = verts[t[0]]
        s2 = verts[t[1]]
        s3 = verts[t[2]]
        area = triangleArea(s1,s2,s3)

	g = [ (s1[0]+s2[0]+s3[0])*oneThird,
              (s1[1]+s2[1]+s3[1])*oneThird,
              (s1[2]+s2[2]+s3[2])*oneThird ]
        volSum += (g[0]*tn[0] + g[1]*tn[1] + g[2]*tn[2])*area
    return volSum*oneThird


def findComponents(verts, faces, normals=None, returnOption=0):
    	"""find the components of a geometry.
normals are normals of verts not faces.
returnOptiont: 	component return option
	= 0:		return all components
	= 1:		return all outside surfaces (volume > 0; normals must be given)
The code is based on the Vision node ConnectedComponents.
"""
        fdict = {} 
        vdict = {} #dictionary with key - vertex index,
                   #value - list of face indices in which the vertex is found

        flag1 = True; flag2 = True
        newfaces = [] 
	newverts = []
	if normals is not None:
	    newnorms = []
        while flag2:
            for i, fs in enumerate(faces):
                for v in fs:
                    if not vdict.has_key(v):
                        vdict[v] = [i]
                    else:
                        vdict[v].append(i)
                fdict[i] = fs
            Vco = faces[0][:]
            newfaces1 = []
	    newverts1 = []
	    if normals is not None:
	        newnorms1 = []
            vertinds = {} # keys - vertex indices from the input verts list
                          # values - new vertex indices of current surface
            vcount = 0
            # find a surface
            while flag1:
                _Vco = []
                flag1 = False
                # find all vertices that share the same triangles with the vertices in Vco.
                for vert in Vco:
                    vfs = vdict[vert]
                    for i in vfs:
                        if fdict.has_key(i):
                            flag1 = True
                            fs = fdict.pop(i)

                            fsnew = [] # remapped face (with new vertex idices)
                            for v in fs:
                                if v not in Vco: 
                                    if v not in _Vco:
                                        _Vco.append(v)
                                if not  vertinds.has_key(v):
                                    vertinds[v] = vcount
                                    newverts1.append(verts[v])
				    if normals is not None:
				        newnorms1.append(normals[v])
                                    fsnew.append(vcount) 
                                    vcount = vcount + 1
                                else:
                                    fsnew.append(vertinds[v])
                            newfaces1.append(fsnew) # add found triangle to the list of triangles of current surface

                Vco  = _Vco
            newfaces.append( newfaces1 )
            newverts.append( newverts1 )
	    if normals is not None:
	        newnorms.append( newnorms1 )
		
            if  len(fdict):
                faces = fdict.values()
                fdict = {}
                vdict = {}
                flag1 = True
            else: flag2 = False


	# return all surfaces
	if returnOption == 0:
	    if normals is not None:
	        return newverts, newfaces, newnorms
	    else:
	        return newverts, newfaces
	    
	# return only outside surfaces
	outverts = []
	outfaces = []
	outnorms = []
        for i in range( len(newfaces) ):
	    Nvert = len(outverts)
	    volume = meshVolume(newverts[i], newnorms[i], newfaces[i])
	    print i, len(newverts[i]), len(newfaces[i]), volume	# TEST
	    if volume > 0:
	        outverts += newverts[i]
		outnorms += newnorms[i]
		# update face vertex indices to reflect lengthened vertices
		newfaces_i_mod = []
		for v1, v2, v3 in newfaces[i]:
		    newfaces_i_mod.append( [v1+Nvert, v2+Nvert, v3+Nvert] )
		outfaces += newfaces_i_mod

	return outverts, outfaces, outnorms

