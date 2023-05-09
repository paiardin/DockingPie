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
# $Header: /mnt/raid/services/cvs/PmvApp/splineFrames.py,v 1.1.4.1 2017/07/13 20:52:52 annao Exp $
#
# $Id: splineFrames.py,v 1.1.4.1 2017/07/13 20:52:52 annao Exp $
#
import numpy
from math import sqrt
from geomutils.geomalgorithms import TriangleNormals

def bspline(v1, v2, v3, t ):
    v4 = [0.,0.,0.,0.]
    frac3 = 0.5 * t*t
    frac1 = 0.5 * (1.-t) * (1.-t)
    frac2 = 1. - (frac1 + frac3)
    for i in (0,1,2,3):
        v4[i] = frac1 * v1[i] + frac2 * v2[i] + frac3 * v3[i]
    return v4


class SplineFrames:
    """
    Class to compute an interpolation spline for helices and approximating
    spline for sheets and coils and a reference frame at each spline point.
    The oxygen atoms are used to control the torsion of plane along the path.
    
    the reference frame is defined by the following 3 unit vectors:
        - parthV (P): a vector pointing from a given spline point to the next
        - pathNV (N): a vector orthiogonal to the path in and in the plane
        - pathBNV (B): a vector orthogonal to the 2 previous ones (i.e. normal
                                                                to the plane)

                             ^  --------------------
                           B | /
                             |/      
                             p1-->  p2-->  p3
                            /    P
                          N/
                           ----------------------
    input:
        caoCoords: coordianates of CA and O atoms for the residues in a chain
        sstype: list secondary structure element types (see prody) default 'C'
        quality: number of points in the spline for each residue
        offset:   amount to offset ctrl points away from CA positions
        off_c:   this parameter has been added to account for DNA/RNA in the
                  coil

    output:
        self.path list of 3D points (quality per residue)
        self.pathV list of vectors pointing from a spline point to the next
        self.parthNV list of vectors orthogonal to
        self.pathBNV: a vector orthogonal to the 2 previous ones (i.e. normal)
        self.edge1
        self.edge2
    """

    def __init__(self, caoCoords=None, quality=5, sstype=None, offset=1.2,
                 off_c=0.5):
        self._nchord = 5
        self._nrib = 2
        self._ribwid = 1.0
        self._cao = None
        self._sstype = None
        self.offset = None
        self._nres = None
        self.path = None
        self.pathV = None
        self.pathNV = None
        
        if caoCoords:
            self.calculate(caoCoords, quality, sstype)
                
    def calculate(self, coords, quality=5, sstype=None, offset=1.2, off_c=0.5, ident=None):
        self._nchord = quality
        self._offset = offset
        self._off_c = off_c
        """
        Generate ctrl points for protein ribbon, based on ideas on
        Carson & Bugg, J.Molec.Graphics 4,121-122 (1986)
    
        Ctrl points for Bspline are generated along a line passing
        through each CA and along the average of the two peptide planes
    
        coords   coordinates of CA and O atoms of all residues in strand
        quality   number of chords/residue
        sstype
        nrib
        offset   amount to offset ctrl points away from CA positions
        off_c   this parameter has been added to account for DNA/RNA in the coil
        ident   is a list of integers uniquely identifying ext point in coords. It will
                be used to build a list of identifier for each point in the smooth path
                tying them back to their control points
        """
        l = len(coords)
        if l<=0: return None
        self._nres = l/2

        if ident is None:
            ident = range(len(coords))
        ## add imaginary CA O before and after real coordinates
        # compute vector CA0-CA1 and create CA O atoms using thsui translation
        vec = [coords[2][0]-coords[0][0],
               coords[2][1]-coords[0][1],
               coords[2][2]-coords[0][2]]
        ca1 = [coords[0][0]-vec[0], coords[0][1]-vec[1], coords[0][2]-vec[2]]
        o1 =  [coords[1][0]-vec[0], coords[1][1]-vec[1], coords[1][2]-vec[2]]
        self._ca1V = vec

        # compute vector CAlast-1-CAlast and create CA O atoms using this translation
        n = len(coords)
        vec = [coords[n-2][0]-coords[n-4][0],
               coords[n-2][1]-coords[n-4][1],
               coords[n-2][2]-coords[n-4][2]]
        ca2 = [coords[n-2][0]+vec[0], coords[n-2][1]+vec[1], coords[n-2][2]+vec[2]]
        o2 =  [coords[n-1][0]+vec[0], coords[n-2][1]+vec[1], coords[n-1][2]+vec[2]]

        coords = numpy.array( [ca1,o1]+ coords.tolist()+[ca2, o2], 'd')
        sstype = [sstype[0]] + sstype.tolist() + [sstype[-1]]
        self._cao = coords
        self._sstype = sstype

        # create control points array
        ctrl = numpy.zeros( (self._nres+1, 2, 4), 'd')

        E = [0.0, 0.0, 0.0]
        F = [0.0, 0.0, 0.0]
        G = [0.0, 0.0, 0.0]
        H = [0.0, 0.0, 0.0]

        if self._nrib > 1: drib = self._ribwid / (self._nrib-1)
        rib2 = (self._nrib+1)*0.5

        # loop over residues
        #for nr in range(self._nres):
        for nr in range(self._nres+1):
            i=2*nr
            if nr<self._nres:
                # A is vector CAi to CAi+1
                #A = coords[i+2] - coords[i]
                A = [coords[i+2][0]-coords[i][0],
                     coords[i+2][1]-coords[i][1],
                     coords[i+2][2]-coords[i][2]]

                # B is vector CAi to Oi
                #B = coords[i+1] - coords[i]
                B = [coords[i+1][0] - coords[i][0],
                     coords[i+1][1] - coords[i][1],
                     coords[i+1][2] - coords[i][2]]

                # C = A x B;  D = C x A
                #C = numpy.array(cross(A,B))
                C = [A[1]*B[2]-B[1]*A[2],
                     A[2]*B[0]-B[2]*A[0],
                     A[0]*B[1]-B[0]*A[1]]

                #D = numpy.array(cross(C,A))
                D = [C[1]*A[2]-A[1]*C[2],
                     C[2]*A[0]-A[2]*C[0],
                     C[0]*A[1]-A[0]*C[1]]

                #D = normalize(D)
                n = 1.0 / sqrt(D[0]*D[0] + D[1]*D[1] + D[2]*D[2])
                D = [D[0]*n, D[1]*n, D[2]*n]

                if i==0:
                    E[:] = D[:] # First peptide, no previous one to average with
                    P = [0.,0.,0.]
                else:
                    # Not first, ribbon cross vector is average of peptide plane
                    # with previous one
                    dot = D[0]*G[0]+D[1]*G[1]+D[2]*G[2]
                    if dot<0.0:
                        B[:] = [-D[0], -D[1], -D[2]]
                    else:
                        B[:] = D[:]
                    E = [G[0]+B[0], G[1]+B[1], G[2]+B[2]]

                    # Offset is along bisector of CA-CA-CA vectors A (H is Ai-1)
                    P = [H[0]-A[0], H[1]-A[1], H[2]-A[2]]
                    #P = normalize(P)
                    n = P[0]*P[0] + P[1]*P[1] + P[2]*P[2]
                    if n == 0.0:
                        p = [0., 0., 0.]
                    else:
                        n = 1.0/sqrt(n)
                        P = [P[0]*n, P[1]*n, P[2]*n]
            else:
                E[:] = G[:]      # Last one, just use last plane
                P = [0.,0.,0.]

            #E = normalize(E) # Normalise vector E
            n = 1.0 / sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2])
            E = [E[0]*n, E[1]*n, E[2]*n]

            if self._sstype is None:
                sst = 'C'
            else:
                sst = self._sstype[nr]
            isHelix = (sst=='H' or sst=='G' or sst=='I')

            # Generate ctrl points
            if not isHelix: # control points for none helical parts
                P = [coords[i][0] + off_c*A[0],
                     coords[i][1] + off_c*A[1],
                     coords[i][2] + off_c*A[2]]
            else: # control point for helices
                  # increasing the offset increases the radius of the helix
                  # increasing the factor multiplying A shifts the per/residue
                  # point down the helix
                  off = 2*offset
                  P = [ off* P[0] + 0.75*A[0],
                        off* P[1] + 0.75*A[1],
                        off* P[2] + 0.75*A[2]]
                  P = [coords[i][0]+P[0], coords[i][1]+P[1], coords[i][2]+P[2]]

            for j in range(self._nrib):
                fr=(float(j+1)-rib2)*drib
                F = [fr*E[0], fr*E[1], fr*E[2]]
                # MS changed.
                ctrl[nr][j][:3] = [P[0]+F[0], P[1]+F[1], P[2]+F[2]]
                ctrl[nr][j][3] = i+2

            # Store things for next residue
            G[:] = E[:]
            H[:] = A[:]

        # end loop over residues
        self.ctrl = ctrl
        self.smooth = smooth = self.ribdrw(ctrl, quality, ident)

        self.edge1 = smooth[0,:,:3]
        self.edge2 = smooth[1,:,:3]
        self.vertices_flat = smooth.reshape( (-1,4) )[:, :3]
        # compute path in middle of the 2D sheet 
        self.path = (smooth[0,:,:3] + smooth[1,:,:3])*0.5
        self.identity = numpy.array(smooth[0,:,3], 'i')
        self.pathV = self.computePathVectors(self.path)
        self.faces2D = self.computeFaces2D(smooth)
        self.pathBNV = self.computeBinormals(smooth, self.faces2D, self.vertices_flat)
        self.pathNV = self.computeNormals(smooth)

    def ribdrw(self, ctrl, nchord, ident, nrib=2 ):
        """
        ctrl: ctrl points 3D
        nchord: how many interpolations per ctrl pt
        ident: ctrl point identifiers
        
        splining from Larry Andrews 7-Nov-1988
        """

        tinc = 1./nchord
        l = (ctrl.shape[0]-1)*nchord
        smooth = numpy.zeros( (nrib, l, 4 ), 'd')
        # calculate spline segments
        for irib in range(nrib):
            # splines go midpoint-to-midpoint
            nb = start = nchord/2
            for pt in range(1,ctrl.shape[0]-1):
                t = 0.0
                for i in range(nchord):
                    smooth[irib][nb] = bspline(ctrl[pt-1][irib], ctrl[pt][irib],
                                               ctrl[pt+1][irib], t )
                    smooth[irib][nb][3] = ident[pt]
                    t = t + tinc
                    nb = nb+1
                    
            # add first start points
            #smooth[irib][nb] = ctrl[-1][irib]
            v = smooth[irib][start]-smooth[irib][start+1]
            for i in range(start):
                smooth[irib][i] = smooth[irib][start]+ v*(start-i)
                smooth[irib][i][3] = ident[0]
                
            # add last points
            v = smooth[irib][nb-1]-smooth[irib][nb-2]
            for i in range(nb, l):
                smooth[irib][i] = smooth[irib][nb-1] + v*(i-nb+1)
                smooth[irib][i][3] = ident[-1]
                
        return smooth#[:, :nb,:]

    def computeFaces2D(self, smooth):
        f = []
        n = smooth.shape[1]
        f = map(lambda x, n1 = n: (x,x+n1,x+n1+1,x+1),
                range(smooth.shape[1]-1))
        faces2D = numpy.array(f).astype('i')
        return faces2D

    def computeBinormals(self, smooth, faces2D, verts2D, fixFlips=True):
        f = numpy.array(faces2D[:,:3])
        binorm = TriangleNormals(verts2D, f, 'PER_VERTEX')
        binorm[smooth.shape[1]-1,:] = binorm[smooth.shape[1]-2,:]
        binorm1 = binorm[:smooth.shape[1]]
        if fixFlips:
            flips = [1]
            flip = 1
            for n in range(len(binorm1)-1):
                dot = binorm1[n][0]*binorm1[n+1][0] + \
                      binorm1[n][1]*binorm1[n+1][1] + \
                      binorm1[n][2]*binorm1[n+1][2]
                if dot < -0.2:
                    flip *= -1
                flips.append(flip)
            return numpy.multiply(binorm1, flip)
        else:
            return binorm1

    def computeNormals(self, smooth):
        #vector in 2D sheet plane orthogonal to spline
        # since width between edge 1 and 2 is 1.0 no normalization needed
        return self.edge2-self.edge1

    def computePathVectors(self, path):
        #vector in 2D sheet from path point to path point
        vect = []
        for i, v in enumerate(path[:-1]):
            diff = path[i+1]-v
            n = 1./sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2])
            vect.append(diff*n)
        vect.append(diff*n)
        return vect

if __name__=='__main__':
    import prody, sys
    # make prody read SS assignement
    prody.confProDy(verbosity='error', auto_secondary=True)

    #
    from time import time
    from MolKit2 import Read
    #mol = Read('1jff.pdb')
    #mol = Read('1crn.pdb')
    mol = Read(sys.argv[1])
    import numpy
    t00 = t0 = time()

    sf = SplineFrames()
    vertSets = []
    for chid in numpy.unique(mol._ag.getChids()):
        caoCoords = mol.select("chid %s name CA O"%chid).getCoords()
        #import pdb; pdb.set_trace()
        sstype = mol.select("protein and chid %s name CA"%chid).getSecstrs()

        print chid, len(caoCoords), caoCoords[0]
        smooth = sf.calculate( caoCoords, sstype=sstype)
        #smooth = ribbon2D( nrib, width, nchords, offset, natoms, caoCoords, sstype )
        vertSets.append( (chid, sf.edge1, sf.edge2, sf.path, sf.pathV, sf.pathNV, sf.pathBNV, sstype, sf._nchord) )

    print 'compute spline frames', time()-t0

    #display = False
    display = True
    if display:
        from PySide import QtGui, QtCore
        app = QtGui.QApplication(sys.argv)

        from DejaVu2.Qt.Viewer import Viewer
        vi = Viewer()
        vi.AddCamera()
        from DejaVu2.Spheres import Spheres
        from DejaVu2.IndexedPolylines import IndexedPolylines
        from DejaVu2.IndexedPolygons import IndexedPolygons
        from extrude import Extruder, Rectangle, Circle
        t0 = time()
        num = 0
        for chid, vs1, vs2, vs3, pv, n, b, sstype, quality in vertSets:
            # add sphere for edges and path
            ## edge1 = Spheres('edge1_%s'%chid, vertices=vs1, radii=(0.15,),
            ##                   inheritMaterial=False, materials=[[1,0,0]])
            ## edge2 = Spheres('edge2_%s'%chid, vertices=vs2, radii=(0.15,),
            ##                 inheritMaterial=False, materials=[[0,1,0]])
            materials = [[1,1,0], [1,0,1], [0,1,1]]
            mat = []
            for i in range(len(caoCoords)/2):
                mat.extend( [materials[i%3]]*quality )
            mat.append([1,1,1])
            path = Spheres('path_%s'%chid, vertices=vs3, radii=(0.1,),
                           inheritMaterial=False, materials=mat)
            
            ## vi.AddObject(edge1)
            ## vi.AddObject(edge2)
            vi.AddObject(path)

            # draw coordinate systems
            ## p1 = numpy.concatenate( (vs3, vs3+pv), 0)
            ## faces = []
            ## length = len(vs3)
            ## for i in range(len(vs3)):
            ##     faces .append( (i, i+length) )
            ## pvs = IndexedPolylines('pv_%s'%chid,vertices=p1, faces=faces,
            ##                        inheritMaterial=False, materials=[[0,1,1]])
            ## #vi.AddObject(pvs)

            ## p1 = numpy.concatenate( (vs3, vs3+n), 0)
            ## normals = IndexedPolylines('IPO_%s'%chid,vertices=p1, faces=faces,
            ##                            inheritMaterial=False, materials=[[1,1,0]])
            ## #vi.AddObject(normals)

            ## p1 = numpy.concatenate( (vs3, vs3+b), 0)
            ## normals = IndexedPolylines('binorm_%s'%chid,vertices=p1, faces=faces,
            ##                            inheritMaterial=False, materials=[[1,0,1]])
            ## #vi.AddObject(normals)

            ## from DejaVu2.Shapes import Rectangle2D, Circle2D

            ## shape1 = Rectangle2D(width=1.2, height=0.2, vertDup=1)
            ## shape2 = Circle2D(radius=0.1)

            ## from PmvApp.extruder import ExtrudeObject, ExtrudeSSElt
            #rect = Rectangle2D(1.5, .2)
            ## fragments = []
            ## curSst = None
            ## curResNum = None
            ## for i in range(len(vs3)):
            ##     resNum = i/quality
            ##     if curResNum != resNum: # new residue
            ##         curResNum = resNum
            ##         try: # try because the last smooth point ends up for nres+1
            ##             sst = sstype[resNum]
            ##         except IndexError:
            ##             pass
            ##         if curSst!=sst: # new SS ELEMENT
            ##             if i > 0:
            ##                 fragments.append( (curSst, matrices) )
            ##             matrices = []
            ##             curSst = sst
            ##     matrices.append( (n[i], b[i], pv[i], vs3[i]) )

            ## for sst, matrices in fragments:
            ##     path = numpy.zeros( (len(matrices), 3) ).astype('f')
            ##     matrix = numpy.zeros( (len(matrices), 3, 3) ).astype('f')
            ##     for i in range(len(matrices)):
            ##         matrix[i][0] = matrices[i][0]
            ##         matrix[i][1] = matrices[i][1]
            ##         matrix[i][2] = matrices[i][3]
            ##         path[i] = matrices[i][2]
            ##         if sst=='H' or sst=='G' or sst=='I':
            ##             shape = shape1
            ##         elif sst=='B' or sst=='E':
            ##             shape = shape1
            ##         else:
            ##             shape = shape2
                        
            ##     extruder1 = ExtrudeObject(path, matrix, shape, cap1=1, cap2=1,
            ##                              larrow=3, warrow=2)
            ##     geom = IndexedPolygons(
            ##         'ss%d'%num, vertices=extruder1.vertices,
            ##         faces=extruder1.faces, vnormals=extruder1.vnormals)
            ##     vi.AddObject(geom)
            ##     num += 1

            ## circ1 = Circle(0.2, quality=10) # for coild
            ## rect = Rectangle(1.5, .2) # for sheets
            ## circ2 = Circle(0.4, quality=10) # for helices
            
            from DejaVu2.Shapes import Rectangle2D, Circle2D

            shape1 = Rectangle2D(width=1.2, height=0.2, vertDup=1)
            shape2 = Circle2D(radius=0.1)
            extruder = Extruder()
            c,f,n = extruder(vs3, pv, n, b, quality, sstype, shape1, shape1, shape2)
            #c,f,vn,fn = extruder(vs3, pv, n, b, quality, sstype, rect, rect, rect)
            #import pdb; pdb.set_trace()
            materials = ([1.,1.,0.],)
            extShapepol = IndexedPolygons('circleEx_%s'%chid,
                                          vertices=c, faces=f, normals=n,
                                          inheritMaterial=False, materials=materials,
                                          inheritCulling=False, culling='none')
            #vi.AddObject(extShapepol)

        # the geometry has sf._res * len(sf.edge1) * len(shape.vertices) vertices
        # 
        print 'extrusion', time()-t0
        sys.exit(app.exec_())
        #execfile('ribbon1.py')
        #pmv.gui().viewer.AddObject(path)
