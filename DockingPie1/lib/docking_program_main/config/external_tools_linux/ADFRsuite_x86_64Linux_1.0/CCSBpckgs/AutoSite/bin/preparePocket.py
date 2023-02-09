import argparse,os,numpy,json
from MolKit2 import Read
import AutoSite
from AutoSite.clusterNode import clusterNode
from AutoSite.clusterNode import buildJSON
from matchPocket import spheres2Grid

from time import time
from math import sqrt

def load_json(filename):
    with open(filename) as f:
        return json.load(f) 


if __name__=='__main__':
    tt0=time()
    parser = argparse.ArgumentParser(description='preparePocket', usage="usage: %(prog)s -s SIZE",
                          version="%prog 0.1")
    parser.add_argument("-f", "--file",
                      dest="file",
                      help="ligand pdb file",)
    parser.add_argument("-s", "--size",
                      dest="size",
                      help="ligand volume in A^2",)
    parser.add_argument("-p", "--pep",
                      dest="pep",
                      help="using peptide scoring function",)
    args = parser.parse_args()
    
    if args.size:
        size=int(args.size)
    elif args.file:
        ligand = Read(args.file)    
        radiiL=ligand._ag.getRadii()
        centerL=ligand._ag.getCoords()
        origin=centerL.tolist()    
        minx=999;
        miny=999;
        minz=999;
        for pocketgrid in origin:
            if pocketgrid[0]<minx:minx=pocketgrid[0]
            if pocketgrid[1]<miny:miny=pocketgrid[1]
            if pocketgrid[2]<minz:minz=pocketgrid[2]
        mini=[minx,miny,minz]    
        ligand_coords,ppts=spheres2Grid(centerL,radiiL,mini,1.0)
        size=len(ppts)
        #print size
    

    
    obj = load_json('pockettree.txt') 
    tree = buildJSON(obj)
    
    finaloutput=tree.getNodebySize(size)
    #import pdb;pdb.set_trace()
    #finaloutput.sort(key=lambda x: x.score, reverse=True)
    os.system('mkdir pockets')
    for pocket in finaloutput:
        os.system('cp *cl_%03d*pdb ./pockets/%03d.pdb' % (pocket.id, pocket.id))
        if args.pep:
            pocket.score=-pocket.totalE*pocket.buriedness*sqrt(pocket.buriedness)
    finaloutput.sort(key=lambda x: x.score, reverse=True)
    i=1
    for pocket in finaloutput:
        print pocket.score,pocket.id,pocket.totalE,pocket.buriedness,pocket.rg,pocket.size,i
        os.system('mv ./pockets/%03d.pdb ./pockets/ranked_%03d.pdb' % (pocket.id, i))
        os.system('cp *fp_%03d*pdb ./pockets/fp_%03d.pdb' %(pocket.id, i))
        i=i+1


    print "pockets prepareation finished in ",time()-tt0," seconds"
    #import pdb; pdb.set_trace()
