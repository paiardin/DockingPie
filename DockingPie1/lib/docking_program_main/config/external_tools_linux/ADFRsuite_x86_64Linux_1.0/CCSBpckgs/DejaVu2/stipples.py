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

########################################################################
#
# Date: 2016 Authors: Michel Sanner
#
#    sanner@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Michel Sanner and TSRI
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/DejaVu2/stipples.py,v 1.1.4.2 2017/10/05 20:58:56 sanner Exp $
#
# $Id: stipples.py,v 1.1.4.2 2017/10/05 20:58:56 sanner Exp $
#
import numpy

def computeHalfBonds(v, f, c, r=None):
    # compute halfbonds
    fcol = []
    nf = []
    v1 = v.tolist()
    if r is not None:
        r1 = []
    nv = len(v)
    #import pdb; pdb.set_trace()
    num = 0
    for face in f:
        c1 = c[face[0]]
        c2 = c[face[1]]
        ## if abs(numpy.sum(c1-c2))< 0.001: # same color
        ##     nf.append(face)
        ##     fcol.append(c1.tolist())
        ## else:
        p1 = v[face[0]]
        p2 = v[face[1]]
        center = (p1+p2)*.5
        v1.append(center)
        nf.append([face[0],nv])
        fcol.append(c1.tolist())
        nf.append([nv, face[1]])
        fcol.append(c2.tolist())
        nv +=1
        if r is not None:
            r1.append(r[num])
            r1.append(r[num])
        num+=1
        
    #import pdb; pdb.set_trace()
    v = numpy.array(v1)
    faces = nf
    if r is not None:
        return v, nf, fcol, r1
    else:
        return v, nf, fcol
        
def stippleLines(v, f, c, r=None, segLen=0.2, spaceLen=0.15):
    #
    # for a set of indexedpolylines described a vector of 3D points v
    # with colors c and lines drawn between verices indexed in f we compute a
    # new set of vertices, faces and colors representing stippled lines.
    # as follows:
    #
    #  p1--   ----   ----   ----   --p2
    #
    # segLen is te length of a stipple
    # spacelen is the length of space between stipples
    # the stippled line always begins and ends with a half stipple
    # the length of the stipple is adjusted for each line so that the length
    # of the line segment matches the sum of length of stipples and spaces
    #
    if len(f)==0:
        if r is not None:
            return [], [], [], [], []
        else:
            return [], [], [], []
        
    points = v[numpy.array(f)]
    vect = points[:,1]-points[:,0]
    vlength = numpy.linalg.norm(vect, axis=1) # vector length
    vlength = vlength.reshape( (-1,1))
    normalizedV = vect/vlength
    # n segments and n spaces = vlen  n = (vlen+spaceL)/(segL+spaceL)
    numSeg = numpy.round((vlength)/(segLen+spaceLen))
    numSeg = numpy.clip(numSeg, 1,max(numSeg))
    # length after rounding, will be longer or short than vlength
    actualLength = numSeg*(segLen + spaceLen)

    # compute per line segLen to fit in vlength
    diff = vlength-actualLength
    perLineSegLen = segLen+(diff/numSeg)

    verts = []
    faces = []
    fcols = []
    vcols = []
    radii = []
    ct = 0
    for i, p in enumerate(points[:,0]):
        segL = perLineSegLen[i]
        l = segL+spaceLen
        # add first half segment
        verts.append( p )
        vcols.append( c[i] )
        p1 = p + normalizedV[i]*segL*.5
        verts.append( p1 )
        vcols.append( c[i] )
        faces.append( [ct, ct+1] )
        fcols.append( c[i] )
        if r is not None:
            radii.append( r[i] )
        ct+=2
        if numSeg[i]>1:
            for j in range(1, numSeg[i]):
                verts.append( p1 + normalizedV[i]*(j*spaceLen + (j-1)*segL) )
                vcols.append( c[i] )
                verts.append( p1 + normalizedV[i]*j*l )
                vcols.append( c[i] )
                faces.append( [ct, ct+1] )
                fcols.append( c[i] )
                ct += 2
                if r is not None:
                    radii.append( r[i] )
        else:
            j=0 # if numSeg[i]==1, j is undefined or has the value from
                # previous iteration so we need to set it to 0

        # add last half segment
        verts.append( p1 + normalizedV[i]*(j*l+spaceLen) )
        vcols.append( c[i] )
        verts.append( points[i,1] )
        vcols.append( c[i] )
        faces.append( [ct, ct+1] )
        fcols.append( c[i] )
        if r is not None:
            radii.append( r[i] )
        ct += 2
    if r is not None:
        return verts, faces, fcols, vcols, radii
    else:
        return verts, faces, fcols, vcols        
