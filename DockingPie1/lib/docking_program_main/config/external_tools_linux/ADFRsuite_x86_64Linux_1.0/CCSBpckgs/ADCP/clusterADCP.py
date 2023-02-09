import os,time
import sys
import numpy
import argparse
from MolKit2 import Read
from mglutil.math.rmsd import RMSDCalculator
from ADFR.utils.cluster import clusterPoses
from prody import writePDB
from MolKit2.AARotamer import AARotamer, CanonicalAARotamers, AARotamerMutator
from prody.measure.contacts import findNeighbors

#Read the energy and Rotamer information for each structures in the pdb
def getEnergy(filename):
    f = open(filename)
    lines = f.readlines()
    f.close()
    if len(lines)==0:
        return []
    totalE=[]
    locE=[]
    extE=[]
    rotamers=[]
    for ln, l in enumerate(lines):
        if not l.startswith('Energy ='):
            continue
        enes=l.split()
        if enes[2] == 'totalE':
            if float(enes[3]) < -1000000 or float(enes[3]) > 1000000:
                return totalE,locE,extE,rotamers
            totalE.append(float(enes[3]))
            locE.append(float(enes[6]))
            #extE.append(float(enes[8]))

            extE.append(0.25*float(enes[3])+0.75*float(enes[8]))
        else:
            totalE.append(float(enes[2]))
            locE.append(float(enes[4]))
            extE.append(float(enes[5]))
        if enes[11] == 'Rotamers:':
            rotstring=''
            for rot in enes[12:]:
                rotstring += str(rot)+" "
            rotamers.append(rotstring)
    return totalE,locE,extE,rotamers

#build the side chain and write out the structure in pdb format.
def buildSC(clusterid,rotamers,receptor=None,outputname=None):
    prot = Read('tmp.pdb')
    mutator = AARotamerMutator()
    i = 0;
    rots = rotamers.split()
    for res in prot._ag.iterResidues():
        resname = res.getResname()
        if resname == 'ALA' or resname == 'GLY':
            i=i+1
            continue
        chid = res.getChid()
        resnum = res.getResnum()
        mutator.mutate(prot, chid, resnum, res.getResname())
        res = prot.select('chid %s and resnum %d'%(chid, resnum))
        rotamer = AARotamer(res, mutator.rotamer.angleDef,mutator.rotamer.angleList)
        res.setCoords(rotamer.getCoordsForRotamer(int(rots[i])))
        i=i+1
    writePDB("%s_ranked_%i.pdb"%(outputname,clusterid+1),prot.select("not deleted and not hydrogen"))
    if receptor is not None:
        pairs = findResPairs(prot._ag.select('not hydrogen and not deleted'),receptor._ag.select("not hydrogen"))
        return pairs, True
    else:
        return True


def findResPairs(refAtoms,receptor):
    pairs = set()
    for pair in findNeighbors(refAtoms,5,receptor):
        pairs.add(str(pair[0].getResnum())+'_'+str(pair[1].getResnum())+'_'+str(pair[1].getChid()))
    return pairs

def clusterPosesInteraction(neighbors, order, cutOff=0.8):
    # cluster a set of solutions in the AutoDock fashion, i.e. use the best
    # solution as a seed and add all solution within cutOff RMSD to this cluster
    # then re-seed the algorithm with the best energy solution not yet clustered
    # until all solution given in "order" are clustered
    #
    # the coordinates are in coords (nsol, natoms, 3)
    # order id a list of indices into coords indicating the subset of
    # solutions to be clustered. These indices are expected to point to
    # solutions sorted by decreasing GA scores
    #
    remainder = order[:]
    clusters = []
    while len(remainder):
        seed = remainder[0]
        #print '%d left to cluster seed=%d'%(len(remainder), seed), remainder
        #import pdb;pdb.set_trace()
        seedPairs = neighbors[seed]
        cluster = [seed]
        notSelected = []
        for i in remainder:
            if i==seed:
                continue
            currPairs = neighbors[i]
            #print '   %d %f'%(i, rmsd)
            overlap = compareContact(seedPairs,currPairs,JC=True)
            if overlap>=cutOff:
                cluster.append(i)
            else:
                notSelected.append(i)
            #print overlap
        #print 'found cluster', cluster, seed
        clusters.append( cluster )
        remainder = notSelected
    return clusters

def compareContact(pairs1,pairs2,JC=False):
    # if JC is true, return the Jaccard coefficient, else, return fraction of contact
    # len1 is the reference
    len1 = len(pairs1)
    len2 = len(pairs2)
    overlap=len(frozenset(pairs1).intersection(pairs2))
    if len1 == 0:
        return 0
    if JC:
        return (overlap+0.0)/(len1+len2-overlap)
    else:
        return (overlap+0.0)/len1

class clusterADCP:

    def __call__(self, **kw):
        
        syst = kw['input']

        if kw['ref']:
            hasRef = True
            ref = Read(kw['ref'])
            refAtoms = ref.select('name C N CA')
            lowestResid = 9999
            for atom in refAtoms:
                if atom.getResnum() < lowestResid:
                    lowestResid = atom.getResnum()
            residDiff = lowestResid - 1
        else:
            hasRef = False


        if kw['rmsd']:
            clusterCutoff = kw['rmsd']
            rmsd = True
        elif kw['nc']:
            if kw['rec'] is None:
                print 'no receptor found, use rmsd instead'
                clusterCutoff = 2.5
                rmsd = True
            else:                
                natContact = True
                rmsd = False
                receptor = Read(kw['rec'])
                clusterCutoff = 0.8
        else:
            clusterCutoff = 2.5
            rmsd = True
        backbone = True        
        modelAtomIndices = []
        models = Read(syst)
        totE,locE,extE,allrotamers=getEnergy(syst)        
        rotamerswithincutoff=[]        
        order=[]
        scores=[]
        bestE = min(extE)
        for i,ext in enumerate(extE):
            if (ext-bestE)<20:
                order.append(i)
                scores.append(ext)

        oorder = numpy.argsort(scores)
        order = numpy.array(order)[oorder]

        if rmsd:
            if hasRef:
                for atom in refAtoms:
                    at = models.select('resnum %d name %s'%(atom.getResnum() - residDiff, atom.getName()))
                    modelAtomIndices.append(at.getIndices()[0])
            else:
                modelAtomIndices.append(models.select('name C N CA').getIndices().tolist())

            rmsdCalc = RMSDCalculator(models._ag._coords[0][modelAtomIndices])
            clusters=clusterPoses(models._ag._coords,order,rmsdCalc,clusterCutoff)

            print "mode |  affinity  | clust. | ref. | clust. | rmsd | best | energy | best |"
            print "     | (kcal/mol) | rmsd   | rmsd |  size  | avg. | rmsd |  avg.  | run  |"
            print "-----+------------+--------+------+--------+------+------+--------+------+"

            eneList = []
            seedList = []
            rsmdsList = []
            bestrmsdRef = 999.
            bestrmsdInd = -1.
            NCRef = 0.
            for i, cl in enumerate(clusters):
                if i >= 100:
                    break;
                rmsdCalc.setRefCoords(models._ag._coords[clusters[0][0]][modelAtomIndices])
                rmsd0 = rmsdCalc.computeRMSD(models._ag._coords[cl[0]][modelAtomIndices])
                if hasRef:
                    rmsdCalc.setRefCoords(refAtoms.getCoords())
                    rmsdRef = rmsdCalc.computeRMSD(models._ag._coords[cl[0]][modelAtomIndices])
                else:
                    rmsdRef = 999.
                rmsds = [rmsdRef]
                if hasRef and rmsdRef < bestrmsdRef:
                    bestrmsdInd = cl[0]
                    bestrmsdRef = rmsdRef
                NoClash = True
                if i < 10 and backbone:
                    writePDB("tmp.pdb",models._ag,cl[0])
                    NoClash = buildSC(i,allrotamers[cl[0]],outputname=kw['input'][:-4])
                    #os.system('rm tmp.pdb')
                ene0 = extE[cl[0]]
                ene = [ene0]
                if len(cl)>1:
                    if hasRef:
                        rmsdCalc.setRefCoords(refAtoms.getCoords())
                    for j in cl[1:]:
                        if hasRef:
                            rmsdRefJ = rmsdCalc.computeRMSD(models._ag._coords[j][modelAtomIndices])
                            rmsds.append(rmsdRefJ)
                            if rmsdRefJ < bestrmsdRef:
                                bestrmsdInd = j
                                bestrmsdRef = rmsdRefJ
                        ene.append(extE[j])
                    print "%4d  %11.1f %7.1f %7.1f  %6d   %5.1f %5.1f %7.1f    %03d "%(
                        i+1, ene0 * 0.59219, rmsd0, rmsdRef, len(cl), numpy.average(rmsds), numpy.min(rmsds), numpy.average(ene) * 0.59219,cl[0])
                else:
                    print "%4d  %11.1f %7.1f %7.1f  %6d      NA      NA    %03d "%(
                        i+1, ene0 * 0.59219, rmsd0, rmsdRef, len(cl),cl[0])    
            if hasRef:
                writePDB("%s_best_%3.1f.pdb"%(kw['input'][:-4],bestrmsdRef),models._ag,bestrmsdInd)


        elif natContact:
            cset = models._ag._coords.shape[0]
            neighbors = []
            start_time = time.time()

            for i in range(cset):       
                pairs = set()
                if i not in order:
                    neighbors.append(pairs)
                    continue
                models._ag.setACSIndex(i)
                for pair in findNeighbors(models._ag.select("name CB or (name CA and resname GLY)"),8,receptor._ag.select("name CB or (name CA and resname GLY)")):
                    #import pdb;pdb.set_trace()
                    pairs.add(str(pair[0].getResnum())+'_'+str(pair[1].getResnum())+'_'+str(pair[1].getChid()))
                neighbors.append(pairs)
            print "finish calculating neighbors for %d poses with %3.1f seconds"%(cset,time.time()-start_time)
            clusters=clusterPosesInteraction(neighbors,order,clusterCutoff)

            if hasRef:
                refpairs = set()
                for pair in findNeighbors(refAtoms,5,receptor._ag.select("not hydrogen")):
                    refpairs.add(str(pair[0].getResnum()-residDiff)+'_'+str(pair[1].getResnum())+'_'+str(pair[1].getChid()))
            print "mode |  affinity  | ref. | clust. | rmsd | energy | best |"
            print "     | (kcal/mol) | fnc  |  size  | stdv |  stdv  | run  |"
            print "-----+------------+------+--------+------+--------+------+"
            eneList = []

            seedList = []
            rsmdsList = []
            bestNCRef = 0.
            bestNCInd = -1.
            NCRef = 0.

            for i, cl in enumerate(clusters):
                NoClash = True
                if backbone and i<100:
                    #prody.writePDB(str(i)+"_bbb.pdb",models._ag,cl[0])
                    writePDB("tmp.pdb",models._ag,cl[0])
                    currpairs, NoClash = buildSC(i,allrotamers[cl[0]],receptor,outputname=kw['input'][:-4])
                    if hasRef:
                        #import pdb;pdb.set_trace()
                        NCRef = compareContact(refpairs,currpairs)
                        if NCRef > bestNCRef:
                            bestNCInd = cl[0]
                            bestNCRef = NCRef
                    ene0 = extE[cl[0]]
                    ene = [ene0]        
                    print "%4d  %11.1f  %7.3f  %6d      NA      NA    %03d "%(i+1, ene0 * 0.59219, NCRef, len(cl),cl[0])
            if hasRef:
                writePDB("%s_best_%3.2f.pdb"%(kw['input'][:-4],bestNCRef),models._ag,bestNCInd)


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='cluster CrankPep output', usage="usage: %(prog)s -i out.pdb -r ref.pdb -c 2.5",version="%prog 0.1")
    parser.add_argument("-i", "--input",dest="input")
    parser.add_argument("-rec", "--rec",dest="rec")
    parser.add_argument("-nc", "--natContacts", type=float,
                       dest="nc", help='native contacts cutoff used in the clustering')
    parser.add_argument("-rmsd", "--rmsd", type=float,
                       dest="rmsd", help='backbone rmsd cutoff used in the clustering')
    parser.add_argument("-ref", "--ref",dest="ref", help='reference peptide structure for calculating rmsd and fnc')
    kw = vars(parser.parse_args())
    runner = clusterADCP()
    runner(**kw)

