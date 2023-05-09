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

class SphereMask:
    """Build a discrete grid of voxels representing a sphere
"""

    def __init__(self, radius):
        """Radius is the the sphere radius (integer).
A 3D array of size 2*radius + 1 will contain the discrete sphere
"""
        import numpy
        size = 2*radius + 1
        self.centerIndex = radius
        self.grid = numpy.zeros((size, size, size), numpy.uint8)
        # put ones in voxels inside or on the sphere
        self.Sphere(radius)


    def Sphere( self, radius ):
        "radius is an integeter"
        x = 0
        y = radius
        delta = 2*(1-radius)
        limit = 0
        while y >= limit:
            if delta < 0:
                d = 2*delta + 2*y - 1
                if d > 0:
                    self.Slice(x,y)
                    x +=1
                    y -=1
                    delta = delta + 2*x -2*y +2
                else:
                    x += 1
                    delta = delta + 2*x +1
            elif delta> 0 :
                d = 2*delta - 2*x -1
                if d > 0:
                    self.Slice(x,y)
                    y -= 1
                    delta = delta -2*y +1
                else:
                    self.Slice(x,y)
                    x += 1
                    y -= 1
                    delta = delta + 2*x -2*y + 2
            else:
                self.Slice(x,y)
                x += 1
                y -= 1
                delta = delta +2*x - 2*y +2


    def Slice(self, r, y):
        x = 0
        z = r
        delta = 2*(1-r)
        limit = 0
        self.Track(x,y,z)
        while z >= limit:
            if delta < 0:
                d = 2*delta + 2*z - 1
                if d > 0:
                  x += 1
                  z -= 1
                  delta = delta + 2*x -2*z +2
                  self.Track( x,y,z)
                  self.Track(-x,y,z)
                else:
                    x += 1
                    delta = delta + 2*x + 1
                    self.Track( x,y,z)
                    self.Track(-x,y,z)
            elif delta > 0:
                d = 2*delta - 2*x -1
                if d > 0:
                    z -= 1
                    delta = delta -2*z + 1
                else:
                    x += 1
                    z -= 1
                    delta = delta +2*x -2*z + 2
                    self.Track( x,y,z)
                    self.Track(-x,y,z)
            else:
                x += 1
                z -= 1
                delta = delta + 2*x -2*z + 2
                self.Track( x,y,z)
                self.Track(-x,y,z)


    def Track(self, x,y,z):
        c = self.centerIndex
        for k in range(-z, z+1):
            self.grid[int(c+x), int(c+y), int(c+k)] = 1
        if y != 0:
            for k in range(-z, z+1):
                self.grid[int(c+x), int(c-y), int(c+k)] = 1


class discreteSpheres:
    """build a 3D grid with 1 inside an on the spheres, and 0 outside
Centers are the centers of the spheres and radii a list of radii.
GridSpacing is list of of 3 float indicating the spacing of the grid in
X,Y and Z.  In the case of uneven grid spacing, the radius of the spheres
would become ellipsoinds, and we simply use the sphere that would contain
the ellipsoid.
GridOrigin is a list of 3 float idicating the lower left back corner
of the grid.
grisSize is is a list of 3 integer indicating the number of voxels along
each dimentsion.
"""
    
    def __init__(self, centers, radii, gridSpacing, gridOrigin,
                 gridSize):
        import numpy
        self.grid = numpy.zeros(gridSize, numpy.uint8)
        self.origin = gridOrigin
        self.spacing = gridSpacing
        self.minSpacing = min(gridSpacing)
        self.positions = []
        self.sphereMasks = {}
        #i=0
        for p in centers:
            self.positions.append( [
                int(round((p[0]-gridOrigin[0])/gridSpacing[0])),
                int(round((p[1]-gridOrigin[1])/gridSpacing[1])),
                int(round((p[2]-gridOrigin[2])/gridSpacing[2])) ] )
            #print i, p, self.positions[-1]
            #i += 1
        for r in radii:
            dr = round(r/self.minSpacing)
            if not self.sphereMasks.has_key(dr):
                self.sphereMasks[dr] = SphereMask( dr )

        largeMap = self.grid
        lx, ly, lz = largeMap.shape
        from numpy import logical_or
        i=0
        for p, r in zip(self.positions, radii):
            #print i
            smallMap = self.sphereMasks[round(r/self.minSpacing)].grid
            offx = smallMap.shape[0]/2
            offy = smallMap.shape[0]/2
            offz = smallMap.shape[0]/2
            #print 'off', offx, offy, offz
            px1 = p[0]-offx
            px2 = p[0]+offx+1
            py1 = p[1]-offy
            py2 = p[1]+offy+1
            pz1 = p[2]-offz
            pz2 = p[2]+offz+1
            ox = oy = oz = 0
            ex,ey,ez = smallMap.shape
            lx,ly,lz = largeMap.shape

            #print 'p1', px1, py1, pz1
            #print 'p2', px2, py2, pz2
            #print 'dims', lx, ly, lz

            if px2<0 or py2<0 or pz2<0 or px1>lx or py1>ly or py1>lz:
                print 'discreteSpheres: WARNING: sphere %d outside grid, it is skipped', i
                continue
            if px1<0:
                ox = -px1
                px1 = 0
            if py1<0:
                oy = -py1
                py1 = 0
            if pz1<0:
                oz = -pz1
                pz1 = 0
            if px2>lx:
                ex = lx-px2
                px2 = lx
            if py2>ly:
                ey = ly-py2
                py2 = ly
            if pz2>lz:
                ez = lz-pz2
                pz2 = lx

            #print 'o', ox, oy, oz
            #print 'e', ex, ey, ez
            
            submap = largeMap[px1:px2, py1:py2, pz1:pz2]
            #print 'submap shape', submap.shape
            largeMap[px1:px2, py1:py2, pz1:pz2] = logical_or(
                submap, smallMap[ox:ex, oy:ey, oz:ez]) 
            i += 1


if __name__=='__main__':
    #sm = SphereMask( 3 )
    #largeMap = numpy.zeros((22, 11, 11), numpy.uint8)
    #stencilSphereMap(largeMap, sm.grid, [(5, 5, 5), (15, 5, 5)])

    from MolKit import Read
    mols = Read('1crn.pdb')
    radii = mols[0].defaultRadii()
    maxr= max(radii)
    coords = mols.chains.residues.atoms.coords[:]
    mini = (min(map(lambda x: x[0], coords)), min(map(lambda x: x[1], coords)),
            min(map(lambda x: x[2], coords)))

    maxi = (max(map(lambda x: x[0], coords)), max(map(lambda x: x[1], coords)),
            max(map(lambda x: x[2], coords)))

    padding = 10
    minimum = map( lambda i: int(mini[i] - maxr) - padding, [0,1,2])
    maximum = map( lambda i: int(maxi[i] + maxr + 1) + padding, [0,1,2])
    delta = [maximum[0]-minimum[0], maximum[1]-minimum[1], maximum[2]-minimum[2]]
    #
    dcpk = discreteCPK(coords, radii, 0.5, minimum, [delta[0]*2,delta[1]*2,delta[2]*2])

