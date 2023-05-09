import argparse,os,numpy
from MolKit2 import Read

def spheres2Grid(centers, radii, origin, spacing):
    """Compute a list of (i,j,k) grid point indices covered by spheres
    centered at 'centers' and with radii 'radii' indices on a grid with
    a given orgin and spacing

    [(i,j,k)], [(x,y,z)] = spheres2Grid(centers, radii, origin, spacing)
    """
    from math import ceil
    gpts = [] # list of (i,j,k) grid points covered by molecule
    pts = [] # list of (x,y,z) grid points covered by molecule
    ox, oy, oz = origin
    spacing1 = 1.0/spacing
    used = {}
    for n in xrange(len(centers)):
        x,y,z = centers[n]
        r = radii[n]
        #r=2
        x0 = int((x-ox)*spacing1)
        y0 = int((y-oy)*spacing1)
        z0 = int((z-oz)*spacing1)
        r0 = int(ceil(r*spacing1))
        r2 = r*r
        #import pdb; pdb.set_trace()
        for i in range(x0-r0, x0+1+r0):
            gx = round(ox + i*spacing,3)
            for j in range(y0-r0, y0+1+r0):
                gy = round(oy + j*spacing,3)
                for k in range(z0-r0, z0+1+r0):
                    gz = round(oz + k*spacing,3)
                    d2 = (gx-x)*(gx-x) + (gy-y)*(gy-y) + (gz-z)*(gz-z)
                    if d2<r2 and not used.has_key((i,j,k)):
                        used[(i,j,k)] = True
                        gpts.append([i,j,k])
                        pts.append([gx,gy,gz])
        #import pdb; pdb.set_trace()
    return gpts, pts


def checkMatch(coords1,coords2):
    coord_dict = {}
    coord_dict={}.fromkeys([str(x) for x in coords2])

    count=0
    for coord in coords1:
        if coord_dict.has_key(str(coord)):
            count=count+1
    return len(coords1),len(coords2),float(count)/(len(coords1)+len(coords2)-count)



if __name__=='__main__':
    parser = argparse.ArgumentParser(description='matchPocket', usage="usage: %(prog)s --receptor or --maps [options] filename",
                          version="%prog 0.1")
    parser.add_argument("-l", "--ligand",
                      dest="ligandFile",
                      help="ligand PDBQT file",)
    parser.add_argument("-p", "--pocket",
                      dest="pocketFile",
                      help="pocket PDBQT file",)
    parser.add_argument("-d", "--distance",
                      dest="distance",
                      help="pc for pocketcenter to any ligandatoms, tp for any translational points to ligand center",)
    args = parser.parse_args()
    distanceoption = args.distance
    if args.pocketFile.endswith('npy'):
        pocket = numpy.load(args.pocketFile)
        numpyPocket = True
    elif args.pocketFile.endswith('pdb') or args.pocketFile.endswith('pdbqt'):
        pocket = Read(args.pocketFile)
        numpyPocket = False
    else:
        print 'wrong pocket file format'
        exit()
    ligand = Read(args.ligandFile)
    if distanceoption is None:
        #radiiL=pocket._ag.getRadii()
        if numpyPocket:
            centerL = pocket
            radiiL=[0.9]*len(pocket)
        else:
            centerL=pocket._ag.getCoords()
            radiiL=[0.9]*len(pocket._ag)
        radiiR=ligand._ag.getRadii()
        centerR=ligand._ag.getCoords()
        origin=centerL.tolist()
        #import pdb;pdb.set_trace()
        minx=999;
        miny=999;
        minz=999;
        for pocketgrid in origin:
            if pocketgrid[0]<minx:minx=pocketgrid[0]
            if pocketgrid[1]<miny:miny=pocketgrid[1]
            if pocketgrid[2]<minz:minz=pocketgrid[2]
        mini=[minx,miny,minz]
        ligand_coords,ppts=spheres2Grid(centerR,radiiR,mini,1.0)
        pocket_coords,origin=spheres2Grid(centerL,radiiL,mini,1.0)
        #import pdb; pdb.set_trace()
        coord_dict = {}
        coord_dict={}.fromkeys([str(x) for x in origin])
        count=0
        for coord in ppts:
            if coord_dict.has_key(str(coord)):
                count=count+1
        print len(radiiL),len(ppts),float(count)/(len(ppts)+len(origin)-count)
        #import pdb;pdb.set_trace()
    elif distanceoption == 'pc':
        import pdb;pdb.set_trace()
        from prody import calcCenter
        if numpyPocket:
            numpy.average(pocket,axis=0)
        else:
            pocketcenter=calcCenter(pocket._ag)
        minDist=9999
        for ligatom in ligand._ag:
            import math,pdb
            #pdb.set_trace()
            if ligatom.getType()=='H':
                continue
            coord=ligatom.getCoords()
            distance=math.sqrt(sum( (a - b)**2 for a, b in zip(coord, pocketcenter)))
            if distance<minDist:
                minDist=distance
        print minDist
    elif distanceoption == 'tp':
        from prody import calcCenter
        ligandcenter=calcCenter(ligand._ag)
        minDist=9999
        if numpyPocket:
            pocketPts = pocket
        else:
            pocketPts = pocket._ag.getCoors()
        for tpts in pocketPts:
            import math
            distance=math.sqrt(sum( (a - b)**2 for a, b in zip(tpts, ligandcenter)))
            if distance<minDist:
                minDist=distance
        print minDist


