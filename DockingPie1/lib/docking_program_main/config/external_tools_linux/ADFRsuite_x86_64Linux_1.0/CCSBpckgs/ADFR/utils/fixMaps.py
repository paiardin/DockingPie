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
# Copyright: M. Sanner and TSRI 2015
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/ADFR/utils/fixMaps.py,v 1.9 2017/06/29 20:31:25 annao Exp $
#
# $Id: fixMaps.py,v 1.9 2017/06/29 20:31:25 annao Exp $
#
import os, time, numpy
import mslib
from time import time
from math import sqrt

from bhtree import bhtreelib
from ADFRcc.adfr import GridMap

def fixMapsFromFiles(rec, tpoints, mapFiles, indent='', myprint=None,
                     msms=None):
    """
    for a list of map files:
      call fixMaps
    """
    maps = []
    for mapFilename in mapFiles:
        _map = GridMap()
        mapsFolder, mapFilename = os.path.split(mapFilename)
        mapType = mapFilename.split('.')[1]
        _map.loadFromMapFile(mapType, mapsFolder, mapFilename)
        maps.append(_map)
    fixmaps(rec, tpoints, maps, indent=indent, myprint=myprint, msms=msms)
    
def getSurface(recAtoms, outsidePoints, indent='', myprint=None, msms=None):
    """
    identify MSMS surface components that is surrounding the outsidePoints
    """
    #import pdb; pdb.set_trace()
    coords = numpy.array(recAtoms.getCoords(), 'f')
    if msms is None:
        t0 = time()
        radii = recAtoms.getData('radius')
        srf = mslib.MSMS(coords=coords, radii=radii)
        # compute reduced surface for all components
        srf.compute_rs(probe_radius=1.5, allComponents=1)
        myprint(indent+'computed receptor suface %.2f (sec)'%(time()-t0))
    else:
        srf = msms
    # find the component closest to the ligand
    # 1 - get vertices for each component
    comp = srf.rsr.fst
    compVerts = []
    allRSv = {}
    rs = []
    while comp:
        rs.append(comp)
        face = comp.ffa
        vd = {}
        while face:
            a, b, c = face._s()
            vd[a] = coords[a]
            vd[b] = coords[b]
            vd[c] = coords[c]
            face = face.nxt
        allRSv.update(vd)
        comp = comp.nxt
        compVerts.append(vd)
    # find smallest distance from outsidePoints atom to RSVertex
    ligAtomsCoords = outsidePoints

    vertInd = allRSv.keys()
    vertInd.sort()
    rsvCoords = []
    for ind in vertInd:
        rsvCoords.append(allRSv[ind])

    bht = bhtreelib.BHtree( rsvCoords, None, 10)
    results = numpy.zeros(5000, 'i')
    dist2 = numpy.zeros(5000, 'f')
    mini = 10000
    minInd = None
    for ligAtCoord in ligAtomsCoords:
        nb = bht.closePointsDist2(tuple(ligAtCoord), 4.0, results, dist2)
        for ind, d2 in zip(results[:nb], dist2[:nb]):
            if d2 < mini:
                mini = d2
                minInd = ind

    minInd = vertInd[minInd]
    # find the components that contain minInd
    comps = []
    for i in range(len(compVerts)):
        if minInd in compVerts[i].keys():
            comps.append(i)

    if len(comps)>1:
        # use the largest one ! . this might not always be right !
        maxi = 0
        for c in comps:
            if len(compVerts[c])> maxi:
                comp = c
                maxi = len(compVerts[c])
    else:
        comp = comps[0]

    myprint( indent+'Using closed surface #%d, receptor atom %s is closest to fill points'%(
        comp,recAtoms._ag[int(minInd)]))

    if rs[comp].ses is None:
        srf.compute_ses(component=comp)
        srf.triangulate(component=comp, density=6.0)

    vf, vi, f = srf.getTriangles()
    verts = vf[:, :3]
    normals = vf[:, 3:6]

    # exclude vertices on edges of analytical surface to avoid
    # singular vertices with bad normals
    bhverts = []
    bhnormals = []
    for i in xrange(len(vf)):
        if vi[i][0] >= 0:
            bhverts.append(verts[i])
            bhnormals.append(normals[i])

    ## # verify that moving receptor atoms are in the cavity
    ## if dockingObject.setting['rmsdRecRef']:
    ##     bhts = bhtreelib.BHtree( bhverts, None, 10)
    ##     Aresults = numpy.zeros(len(bhverts), 'i')
    ##     Adist2 = numpy.zeros(len(bhverts), 'f')
    ##     movingRecAtoms = dockingObject.sortedRecRefAts
    ##     #outsideAtoms = []
    ##     for atom in movingRecAtoms:
    ##         if atom.element=='H':
    ##             atom.outside = None
    ##             continue
    ##         pt = atom.coords
    ##         cut = 2.0
    ##         nb = 0
    ##         while nb==0:
    ##             nb = bhts.closePointsDist2(tuple(pt), cut, Aresults, Adist2)
    ##             if nb == 0:
    ##                 cut += 1.0

    ##         closestSurfInd = Aresults[numpy.argmin(Adist2[:nb])]
    ##         clSP = bhverts[closestSurfInd]
    ##         clN = bhnormals[closestSurfInd]

    ##         # vector for surface to atom
    ##         v = [ pt[0]-clSP[0], pt[1]-clSP[1], pt[2]-clSP[2]]

    ##         # dot product of surface normal with v
    ##         dot = clN[0]*v[0] + clN[1]*v[1] + clN[2]*v[2]
    ##         if dot < 0:
    ##             # implments that only atoms more the 1.0 outside surface
    ##             # make the side chain be rejected
    ##             #n = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])
    ##             #print 'ATOM outside cavity ', atom.full_name(), 'by', n
    ##             #if n>1.0:
    ##             #    outsideAtoms.append(atom)

    ##             # tag atom as outside
    ##             atom.outside = True
    ##         else:
    ##             atom.outside= False
        ## #import pdb
        ## #pdb.set_trace()
        ## from MolKit.protein import ResidueSet
        ## notInPocket = ResidueSet([])
        ## for res in movingRecAtoms.parent.uniq():
        ##     nba = 0
        ##     nbaout = 0
        ##     for a in movingRecAtoms:
        ##         if a.parent != res: continue
        ##         if a.element != 'H':
        ##             nba +=1 # one more atom in this residue
        ##         #print a.name
        ##         if a.outside:
        ##             #print "OUT",a.name
        ##             nbaout +=1
        ##         del a.outside
        ##     #print 'AAAA', res.name, nba, nbaout
        ##     if nbaout==nba: # all are out
        ##         notInPocket.append(res)

        ## if len(notInPocket):
        ##     print 'Side chain atoms for residue are not in pocket', notInPocket.full_name()
        ##     raise ValueError

        #if len(outsideAtoms):
        #    print 'Residue(s) are not in pocket, please remove from flexRec'
        #    for res in AtomSet(outsideAtoms).parent.uniq():
        #        #if len(res.atoms) == len(outsideAtoms):
        #        print "    ", res.name
        #    raise ValueError
    return verts, normals, srf, comp

def fixmaps(recAtoms, tpoints, maps, indent='', myprint=None,
            msms=None):
    """
    for a list of map objects fix the maps
    - recAtoms should be the list of static receptor atoms (ca, cb of moving
    side chains are considered as moving and shoudl not be in recAtoms.
    - tpoints needs to be a list of at least 1 3D point located in the pocket.
    These points are used to select a surface component
    
    """
    #if myprint is None:
    #    myprint = print
    coords = numpy.array(recAtoms.getCoords(), 'f')
    radii = recAtoms._ag._data['radius']
    bht = bhtreelib.BHtree( coords, None, 10)
    result = numpy.zeros( (5000,), 'i' )
    dist2 = numpy.zeros( (5000,), 'f' )
    
    ox, oy, oz = maps[0].getOriginPy()
    sx = sy = sz =  maps[0].getDistBetweenGridPoints()
    nbptx, nbpty, nbptz = maps[0].getNumGridPointsPy()
    #sizeX = (nbptx-1)*sx
    #sizeY = (nbpty-1)*sy
    #sizeZ = (nbptz-1)*sz
    keptPts = []
    cut = 0.0
    nbPoints = nbptx*nbpty*nbptz
    inBoxPtsCounter = 0

    mini = 0.0
    maxi = 1.0
    removedPts = []
    removedijk = []
    keptPts = []
    keptijk = []
    ## put grid point falling inside receptor atomic spheres into removedpts
    ## and others into keptPts

    t0 = time()
    for i in range(nbptx):
        x = ox+i*sx
        for j in range(nbpty):
            y = oy+j*sy
            for k in range(nbptz):
                z = oz+k*sz
                inBoxPtsCounter += 1

                # check if grid point is far enough from closest receptor atom
                nb = bht.closePointsDist2((x,y,z), 1.0, result, dist2)
                keep = True
                for aind, anum in enumerate(result[:nb]):
                    if dist2[aind] < 12.0: # grid point is too close to receptor atom
                        keep = False
                        break

                if keep:
                    keptPts.append( (x, y, z) )
                    keptijk.append( (i, j, k) )
                else:
                    removedPts.append( (x,y,z) )
                    removedijk.append( (i,j,k) )

    #myprint(indent+'removed %d points in %.2f (sec)'%(time()-t0,
    #                                                  len(removedPts)))

    verts, normals, srf, comp = getSurface(
        recAtoms, tpoints, indent=indent, msms=msms, myprint=myprint)
        
    #import pdb; pdb.set_trace()
    srfVerticesBHT = bhtreelib.BHtree(verts, None, 10)
    
    t0 = time()
    keptPts2 = []
    keptijk2 = []
    for pt in keptPts:
        cut = 2.0
        nb = 0
        while nb==0:
            nb = srfVerticesBHT.closePointsDist2(tuple(pt), cut, result, dist2)
            cut += 10.
        vertInd = result[numpy.argmin(dist2[:nb])]
        vx, vy, vz = verts[vertInd]
        n1x = pt[0]-vx
        n1y = pt[1]-vy
        n1z = pt[2]-vz
        nx, ny, nz = normals[vertInd]
        dot = nx*n1x + ny*n1y + nz*n1z
        if dot <= 0.0: # point is inside the surface
            # compute distance to closest surface point
            d2 = n1x*n1x + n1y*n1y + n1z*n1z
            if d2>=0.5625: # if less than 0.75**2 keep it
                removedPts.append(pt)
                removedijk.append( (i,j,k) )
        else:
            keptPts2.append(pt)
            keptijk2.append( (i, j, k) )
            
    myprint(indent+'removed %d points inside the surface in %.2f (sec)'%(
        len(removedPts), time()-t0))
    #import pdb; pdb.set_trace()
    
    t0 = time()
    goodPointsBHT = bhtreelib.BHtree(keptPts2, None, 10)
    result = numpy.zeros( (5000,), 'i' )
    dist2 = numpy.zeros( (5000,), 'f' )

    newData = {}
    oldData = {}
    for amap in maps:
        atype = amap.getMapType()
        newData[atype] = amap.getGridDataPy()
        oldData[atype] = amap.getGridDataPy()

    for pt, ijk in zip(removedPts, removedijk):
        x,y,z = pt
        i,j,k = ijk
        cut = 2.0
        nb = 0
        while nb==0:
            nb = goodPointsBHT.closePointsDist2((x,y,z), cut, result, dist2)
            cut += 6.

        indexOfMinDistVert = result[numpy.argmin(dist2[:nb])]
        vx, vy, vz = keptPts2[indexOfMinDistVert]
        d2 = (x-vx)*(x-vx) + (y-vy)*(y-vy) +(z-vz)*(z-vz)
        im = 0
        for amap1 in maps:
            atype1 = amap1.getMapType()
            #value = values[atype][indexOfMinDistVert]
            value = oldData[atype][keptijk2[indexOfMinDistVert]]
            if d2 < 1.0:
                newData[atype1][i][j][k] = value + sqrt(d2)*100
            else:
                newData[atype1][i][j][k] = value + 100 + d2*1000

    myprint(indent+'Fixed %d maps in %.2f (sec)'%(len(maps), time()-t0))

    # write modified data back into maps and save file
    for amap in maps:
        atype = amap.getMapType()
        amap.setGridData(newData[atype])
        nameWithPath = amap.getMapFilePath()
        folder, name = os.path.split(nameWithPath)
        amap.saveToMapFile(folder, name)

    # use old map reader and writer to save
    ## from Volume.IO.AutoGridWriter import WriteAutoGrid
    ## from Volume.IO.AutoGridReader import ReadAutoGrid
    ## reader = ReadAutoGrid()
    ## writer = WriteAutoGrid()
    ## for amap in maps:
    ##     atype = amap.getMapType()
    ##     path = amap.getMapFilePath()
    ##     name = os.path.split(path)[1].split('.')[0]
    ##     oldmap = reader.read(path, 0)
    ##     oldmap.data[:, :, :] = newData[atype]
    ##     p = os.path.split(path)[0]
    ##     #writer.write(oldmap, '%s_fixed.%s.map'%(os.path.join(p, name) , atype))
    ##     writer.write(oldmap, '%s%s.map'%(os.path.join(p, name) , atype))

    myprint(indent+'Search space reduced by %5.2f percent (%d/%d) and %5.2f (%d/%d) total in %f'%(
        100.-100*len(keptPts2)/float(inBoxPtsCounter), len(keptPts2), inBoxPtsCounter,
        100.-100*len(keptPts2)/float(nbPoints), len(keptPts2), nbPoints, time()-t0))

if __name__=='__main__':
    import sys, numpy
    from glob import glob
    from MolKit2 import Read
    from time import time
    rec = Read(sys.argv[1])
    tpoints = numpy.load(sys.argv[2])
    t0 = time()
    maps = glob(sys.argv[3])
    # assuming no flexible residues
    recAtoms = rec.select()
    fixMapsFromFiles(recAtoms, tpoints, maps)
    print '    TOTAL TIME FOR FIXING MAPS', time()-t0
