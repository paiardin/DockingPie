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

from time import time
from mglutil.math.rotax import rotax, rotVectToVect
import numpy

def vdiff(p1, p2):
    # returns p1 - p2
    x1,y1,z1 = p1
    x2,y2,z2 = p2
    return (x1-x2, y1-y2, z1-z2)
    
from math import sqrt
def vnorm(v1):
    x1,y1,z1 = v1
    n = 1./sqrt(x1*x1 + y1*y1 + z1*z1)
    return (x1*n, y1*n, z1*n)

def vlen(v1):
    x1,y1,z1 = v1
    return sqrt(x1*x1 + y1*y1 + z1*z1)

def dot(v1,v2):
    x1,y1,z1 = v1
    x2,y2,z2 = v2
    return (x1*x2) + (y1*y2) + (z1*z2)

def lengthOfProjectedVects(v1, v2):
    """
    return the length of the v1 projected onto v2
    """
    # v1.v2 = |v1|*|v2|*cos(angle)
    # |v1|*cos(angle) = (v1.v2) / |v2|
    x1,y1,z1 = v1
    x2,y2,z2 = v2
    nv2 = sqrt(x2*x2 + y2*y2 + z2*z2)
    s = ( (x1*x2) + (y1*y2) + (z1*z2) ) / nv2
    vpx, vpy, vpz = s*x2, s*y2,s*z2
    return sqrt(( vpx*vpx + vpy*vpy + vpz*vpz))


def area2(a,b,c):
    """
    return 2 * area of triangle abc
    """
    return (a[0]-c[0])*(b[1]-c[1]) - (a[1]-c[1])*(b[0]-c[0])


def insideTriangle(a, b, c, p):
    """
    returns True if P lies within triangle A,B,C
    ABC is counterclockwise
    """
    return area2(a,b,p)>=0 and area2(b,c,p)>=0 and area2(c,a,p)>=0


def triangulate(verts):
    tri = []
    vertsInd = range(len(verts))
    while len(vertsInd) > 2:
        triFound = False
        count = 0
        while not triFound and count < len(vertsInd):
            count += 1
            a, b, c = vertsInd[:3]
            ca = verts[a]
            cb = verts[b]
            cc = verts[c]
            if area2(ca, cb, cc) >= 0:
                for p in vertsInd[3:]:
                    cp = verts[p]
                    if insideTriangle(ca,cb,cc,cp):
                        break
                else:
                    tri.append( (a,c,b) )
                    vertsInd.remove(b)
                    triFound = True
            vertsInd.append(vertsInd.pop(0))
        if count == len(vertsInd):
            print "ERROR Not a simple polygon"
            return []
    return tri

from math import acos, pi

def angleBetweenNormVecsRad(v1, v2):
    """return the angle between vetors v1 and v2. NOTE v1 and v2 have to be
    normalized"""
    # v1.v2 = |v1|*|v2|*cos(angle)
    # angle = acos( v1.v2 / |v1|*|v2|)
    x1,y1,z1 = v1
    x2,y2,z2 = v2
    return acos( (x1*x2) + (y1*y2) +  (z1*z2) )


def angleBetweenVecsDeg(v1, v2):
    """
    return the angle between vetors v1 and v2.
    """
    # v1.v2 = |v1|*|v2|*cos(angle)
    # angle = acos( v1.v2 / |v1|*|v2|)
    x1,y1,z1 = v1
    x2,y2,z2 = v2
    n1 = vlen(v1)
    n2 = vlen(v2)
    arad = acos( (( x1*x2 ) + (y1*y2) + (z1*z2)) / (n1*n2) )
    return 180 * arad /pi

def ccw(poly):
    n = len(poly)
    k = poly.index(min(poly))
    return area2(poly[k-1], poly[k], poly[(k+1)%n]) > 0


def capMesh(vertices, edges, faces, faceEdges, pp, pno, minDist=0.01):
    """
    Finds the section mesh between a mesh and a plane 
	
    vertices: mesh vertices
    edges: mesh edges as (v1,v2)
    faces: mesh faces
    faceEdges: 
    pp: Vector - A point on the plane
    pno: Vector - The cutting plane's normal (normalized)
    minDist: useds to remove points to close in contour
    
    Returns: Mesh - the resulting mesh of the section if any or
                    Boolean - False if no section exists
    """	
    from math import pi
    halfPi = pi*.5
    
##     ppx, ppy, ppz = pp
##     pnox, pnoy, pnoz = pno

##     t1 = time()
##     vertsold = []
##     ed_xsect = {}
##     pr1 = []
##     pr2 = []
##     ang1 = []
##     ang2 = []
##     for ednum, ed in enumerate(edges):
##         # getting a vector from each edge vertices to a point on the plane
##         v1 = vertices[ed[0]]
##         co1 = vdiff( v1, pp)

##         v2 = vertices[ed[1]]
##         co2 = vdiff( v2, pp)

##         # projecting them on the normal vector
##         #proj1 = lengthOfProjectedVects(co1, pno)
##         x1,y1,z1 = co1
##         x2,y2,z2 = pno
##         s = ( (x1*x2) + (y1*y2) + (z1*z2) )
##         vpx, vpy, vpz = s*x2, s*y2,s*z2
##         proj1 = sqrt( (vpx*vpx + vpy*vpy + vpz*vpz) )
##         pr1.append(proj1)
        
##         #proj2 = lengthOfProjectedVects(co2, pno)
##         x1,y1,z1 = co2
##         s = ( (x1*x2) + (y1*y2) + (z1*z2) )
##         vpx, vpy, vpz = s*x2, s*y2,s*z2
##         proj2 = sqrt( (vpx*vpx + vpy*vpy + vpz*vpz) )
##         pr2.append(proj2)

##         if (proj1 != 0):
##             angle1 = angleBetweenVecsDeg(co1, pno)
##         else:
##             print '0 eangle1', ed
##             angle1 = 0

##         if (proj2 != 0):
##             angle2 = angleBetweenVecsDeg(co2, pno)
##         else:
##             angle2 = 0
##             print '0 eangle2', ed

##         ang1.append(angle1)
##         ang2.append(angle2)
##         #Check to see if edge intersects. Also check if edge is coplanar to the
##         #cutting plane (proj1=proj2=0)		
##         if ((proj1 == 0) or (proj2 == 0) or \
##             (angle1 > 90) != (angle2 > 90)) and \
##             (proj1+proj2 > 0) :

##             #edge intersects.
##             proj1 /= proj1+proj2
##             cox = (v2[0]-v1[0])*proj1 + v1[0]
##             coy = (v2[1]-v1[1])*proj1 + v1[1]
##             coz = (v2[2]-v1[2])*proj1 + v1[2]
##             vertsold.append( (cox,coy,coz) )
##             #store a mapping between the new vertices and the mesh's edges
##             ed_xsect[ednum] = len(ed_xsect)
##     print '    time to compute contour', time()-t1, len(vertsold)

    
    t1 = time()
    edges = numpy.array(edges)
    nv1 = numpy.take(vertices, edges[:,0], axis=0)
    nv2 = numpy.take(vertices, edges[:,1], axis=0)
    nco1 = nv1-pp
    ns = numpy.sum( nco1*pno, 1 )
    ns.shape = -1,1
    nvp = pno*ns
    nproj1 = numpy.sqrt(numpy.sum( nvp*nvp, 1 ))

    nco2 = nv2-pp
    ns = numpy.sum( nco2*pno, 1 )
    ns.shape = -1,1
    nvp = pno*ns
    nproj2 = numpy.sqrt(numpy.sum( nvp*nvp, 1 ))

    lpno = vlen(pno)
    nlco1 = numpy.sqrt( numpy.sum( nv1*nv1, 1 ) )
    denom = nlco1*lpno # what id denom is 0. i.e proj=0
    arad1 = numpy.arccos( numpy.sum( nco1*pno, 1) / denom )# * 180 / pi

    nlco2 = numpy.sqrt( numpy.sum( nv2*nv2, 1 ) )
    denom = nlco2*lpno # what id denom is 0. i.e proj=0
    arad2 = numpy.arccos( numpy.sum( nco2*pno, 1) / denom )# * 180 / pi

    verts = []
    contourVerts = {}
    for i in xrange(len(edges)):
        proj1 = nproj1[i]
        proj2 = nproj2[i]
        if ((proj1 == 0) or (proj2 == 0) or \
            (arad1[i] > halfPi) != (arad2[i] > halfPi)) and \
            (proj1+proj2 > 0) :

            #edge intersects.
            proj1 /= proj1+proj2
            v1i, v2i = edges[i]
            v1 = vertices[v1i]
            v2 = vertices[v2i]
            cox = (v2[0]-v1[0])*proj1 + v1[0]
            coy = (v2[1]-v1[1])*proj1 + v1[1]
            coz = (v2[2]-v1[2])*proj1 + v1[2]
            verts.append( (cox,coy,coz) )
            #store a mapping between the new vertices and the mesh's edges
            contourVerts[i] = len(contourVerts)
    print '    time to compute contour1', time()-t1, len(verts)

    if len(verts)==0:
        return [], []
    
    # build list of edges in contour
    t1 = time()
    edgesS = []
    vertToEdges = {} # used to find the 2 edges a vertex belong to
    nbe = 0
    for fe in faceEdges:
        ps = [ contourVerts[key] for key in fe if key in contourVerts]
        if len(ps) == 2:
            v1,v2 = ps

            if v1 in vertToEdges:
                vertToEdges[v1].append(nbe)
            else:
                vertToEdges[v1] = [nbe]

            if v2 in vertToEdges:
                vertToEdges[v2].append(nbe)
            else:
                vertToEdges[v2] = [nbe]
            edgesS.append(tuple(ps))
            nbe += 1
    print '   time to compute edge list', time()-t1

    # compute rotation matrix to cut plane into a Z plane
    rotMat = numpy.array( rotVectToVect(pno, (0,0,1)), 'f')
    vertsR = numpy.dot( verts, rotMat[:3, :3] ).tolist()

    # order vertices along contour
    t0 = time()
    edgeList = range(0,len(edgesS))
    allFaces = []
    allVertices = []
    while len(edgeList):
        t1 = time()
        curE = edgeList[0]
        fstV, curV = edgesS[edgeList[0]]
        ovi = [fstV]
        ov = [vertsR[fstV]]
        ovnr = [verts[fstV]]
        edgeList.remove(curE)
        while curV != fstV:
            ovi.append(curV)
            ov.append(vertsR[curV])
            ovnr.append(verts[curV])
            length = 0
            while length < minDist:
                nextEdge = vertToEdges[curV][0]
                if nextEdge==curE: nextEdge = vertToEdges[curV][1]
                curE = nextEdge
                edgeList.remove(curE)
                v2 = edgesS[curE][0]
                if v2==curV: v2=edgesS[curE][1]
                length = vlen(vdiff(vertsR[ovi[-1]], vertsR[v2]))
                curV = v2
                if curV == fstV: break

        # make sure vertice are counter clockwise
        if not ccw(ov):
            print "reversing"
            ov.reverse()
            ovnr.reverse()
        
        print '           time to order pacth vertices', time()-t1
        #print 'contour length', len(ovi), len(ov)

        # apply rotation to vertices
        #ovr = numpy.dot( ov, rotMat[:3, :3] ).tolist()

        #for v in ovr:
        #    print v[2],
        #print
        #print ovi
        
        # triangulate
        t1 = time()
        facesS = triangulate(ov)
        print '           time to triangulate patch', time()-t1

        allFaces.append(facesS)
        allVertices.append(ovnr)
    print '    time to triangulate all', time()-t0
        
    return allVertices, allFaces


if __name__ == '__main__':
    from DejaVu2.IndexedPolygons import IndexedPolygonsFromFile
    geomS = IndexedPolygonsFromFile('cv', 'mesh')
    faces = geomS.getFaces()
    vertices = geomS.getVertices()
    vnormals = geomS.getVNormals()

    planePoint = (3.,0.,0.)
    planeNormal = (1.,0.,0.)

    vertsC, facesC = cap(vertices, faces, planePoint, planeNormal)

    from DejaVu2 import Viewer
    vi = Viewer()
    
    from DejaVu2.IndexedPolygons import IndexedPolygons
    tet = IndexedPolygons('surfaceMesh', vertices=vertices,
                          faces=faces, vnormals=vnormals,
                          inheritFrontPolyMode=False,
                          frontPolyMode='line',
                          inheritCulling=0, culling='none',
                          inheritShading=0, shading='flat', visible=0)
    vi.AddObject(tet)

    #from DejaVu2.Spheres import Spheres
    #sectionVsph = Spheres('sectionV', vertices=vertsS, radii=(0.02,))
    #vi.AddObject(sectionVsph)
    #from DejaVu2.Points import Points
    #sectionVpts = Points('sectionV', vertices=allvertsS)
    #sectionVpts = Points('sectionV', vertices=[vertsS[14], vertsS[159], vertsS[158]])
    #vi.AddObject(sectionVpts)

    #from DejaVu2.IndexedPolylines import IndexedPolylines
    #line = IndexedPolylines(
    #    'sectionL', vertices=allvertsS, faces=edgesS,
    #    inheritLinewidth=False, Linewidth=3,
    #    inheritMaterial=False, materials=( (0,1,0), ))
    #vi.AddObject(line)

    #from DejaVu2.Polylines import Polylines
    #line = Polylines(
    #    'contour', vertices=(vertsS[:5],),
    #    inheritLinewidth=False, Linewidth=3,
    #    inheritMaterial=False, materials=( (1,1,0), ))
    #vi.AddObject(line)
    i = 0
    for v, f in zip(vertsC, facesC):
        cap = IndexedPolygons('cap%d'%i, vertices=v, faces=f, 
                              inheritFrontPolyMode=False,
                              frontPolyMode='fill',
                              inheritCulling=0, culling='none',
                              inheritShading=0, shading='flat')
        vi.AddObject(cap)
        i += 1


