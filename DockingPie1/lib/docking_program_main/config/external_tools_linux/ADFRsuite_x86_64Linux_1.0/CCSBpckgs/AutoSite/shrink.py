import os,numpy
import argparse
from MolKit2 import Read

def shrinkPocket(arraypoints, finalsize=800, padding=1.0, spacing=1.0, origin=None):
        #centersL = self.centersL
    #import pdb;pdb.set_trace()
    points = []
    points.append([tuple(x) for x in arraypoints])
    lastset=points[0]
    currset=points[0]
    iteration=0
    while iteration<30:
        iteration=iteration+1
        #import pdb;pdb.set_trace()
        for point in lastset:
            counter=0
            neighbors=getNeighbors(point,spacing,dim=3)
            for neighbor in neighbors:
                if neighbor in lastset:
                    #import pdb;pdb.set_trace()                    
                    counter=counter+1
            if counter<=len(neighbors)/2:
                currset.remove(point)
                #removed.append(point)
        if len(currset)<finalsize*1.2:
            break
        lastset=currset
    return currset


def getNeighbors(point,spacing=1.0,dim=2):
    neighbors=[]
    neighbors.append((point[0],point[1]+spacing,point[2]))
    neighbors.append((point[0]+spacing,point[1],point[2]))
    neighbors.append((point[0],point[1],point[2]+spacing))
    neighbors.append((point[0],point[1]-spacing,point[2]))
    neighbors.append((point[0]-spacing,point[1],point[2]))
    neighbors.append((point[0],point[1],point[2]-spacing))
    if dim==1:
        return neighbors
    neighbors.append((point[0]+spacing,point[1]+spacing,point[2]))
    neighbors.append((point[0],point[1]+spacing,point[2]+spacing))
    neighbors.append((point[0]+spacing,point[1],point[2]+spacing))
    neighbors.append((point[0]-spacing,point[1]-spacing,point[2]))
    neighbors.append((point[0],point[1]-spacing,point[2]-spacing))
    neighbors.append((point[0]-spacing,point[1],point[2]-spacing))
    neighbors.append((point[0]+spacing,point[1]-spacing,point[2]))
    neighbors.append((point[0],point[1]+spacing,point[2]-spacing))
    neighbors.append((point[0]+spacing,point[1],point[2]-spacing))
    neighbors.append((point[0]-spacing,point[1]+spacing,point[2]))
    neighbors.append((point[0],point[1]-spacing,point[2]+spacing))
    neighbors.append((point[0]-spacing,point[1],point[2]+spacing))
    neighbors.append((point[0]-spacing,point[1],point[2]+spacing))
    if dim==2:
        return neighbors
    neighbors.append((point[0]-spacing,point[1]+spacing,point[2]-spacing))
    neighbors.append((point[0]-spacing,point[1]-spacing,point[2]-spacing))
    neighbors.append((point[0]-spacing,point[1]+spacing,point[2]+spacing))
    neighbors.append((point[0]-spacing,point[1]-spacing,point[2]+spacing))
    neighbors.append((point[0]+spacing,point[1]+spacing,point[2]-spacing))
    neighbors.append((point[0]+spacing,point[1]-spacing,point[2]-spacing))
    neighbors.append((point[0]+spacing,point[1]+spacing,point[2]+spacing))
    neighbors.append((point[0]+spacing,point[1]-spacing,point[2]+spacing))

    return neighbors

def writePDB(coords, potentials, atypes, filename):
    f = open(filename, 'w')
    for n, (xyz, pot, aty) in enumerate(zip(coords, potentials, atypes)):
        x,y,z = xyz
        f.write('ATOM%7d%4s   FP     1    %8.3f%8.3f%8.3f  1.0%8.3f    0.001 %s\n' %
                (n+1,aty[0], x, y, z, abs(float(pot)), aty[0]))
                        #(n+1,aty[0]+'%s'%n, x, y, z, abs(pot), aty[0]))
    f.close()


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
        x0 = int((x-ox)*spacing1)
        y0 = int((y-oy)*spacing1)
        z0 = int((z-oz)*spacing1)
        r0 = int(ceil(r*spacing1))
        r2 = r*r
        for i in range(x0-r0, x0+1+r0):
            gx = ox + i*spacing
            for j in range(y0-r0, y0+1+r0):
                gy = oy + j*spacing
                for k in range(z0-r0, z0+1+r0):
                    gz = oz + k*spacing
                    d2 = (gx-x)*(gx-x) + (gy-y)*(gy-y) + (gz-z)*(gz-z)
                    if d2<r2 and not used.has_key((i,j,k)):
                        used[(i,j,k)] = True
                        gpts.append((i,j,k))
                        pts.append((gx,gy,gz))

    return gpts, pts



if __name__=='__main__':
    parser = argparse.ArgumentParser(description='matchPocket', usage="usage: %(prog)s --receptor or --maps [options] filename",
                          version="%prog 0.1")
    parser.add_argument("-p", "--pocket",
                      dest="pocketFile",
                      help="ligand PDBQT file",)
    parser.add_argument("-s", "--size",
                      dest="size",
                      help="pocket PDBQT file",)
    args = parser.parse_args()
    
    pocket = Read(args.pocketFile)
    size = int(args.size)
    if size==0:
        size=len(pocket._ag)/5

    #ligand=Read("test.pdb")
    
    radiiR=[0.9]*len(pocket._ag)
    centerR=pocket._ag.getCoords()
    mini=[0,0,0]
    ppts,pocket_coords=spheres2Grid(centerR,radiiR,mini,1.0)
    
    #writePDB(ligand_coords, ['0']*len(ligand_coords),['C']*len(ligand_coords), 'grid.pdb')
    
    
    test=shrinkPocket(pocket_coords,finalsize=size)
    
    
    numpy.save("lll.npy", test)
    
    
    writePDB(test, ['0']*len(test),['C']*len(test), '%s_root.pdb'%(args.pocketFile[:-4]))
    
