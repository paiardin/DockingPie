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
# Author: Pradeep Anand Ravindranath
# Revision: Michel F. Sanner Aug 2016
#
# Copyright: Pradeep Anand Ravindranath and TSRI 2016
#
#########################################################################
#
# $Header: /opt/cvs/AutoSite/bin/AS.py,v 1.11 2017/06/09 23:03:28 zhangyq Exp $
#
# $Id: AS.py,v 1.11 2017/06/09 23:03:28 zhangyq Exp $
#

import argparse, numpy, os, sys, platform, datetime
import AutoSite
from time import time

from MolKit2 import Read
from ADFRcc.adfr import GridMap
parser = argparse.ArgumentParser(description='AutoSite', usage="usage: %(prog)s --receptor or --maps [options] filename",
                      version="%prog 0.1")
parser.add_argument("-r", "--receptor",
                  dest="receptorFile",
                  help="receptor PDBQT file",)
parser.add_argument("-m", "--maps", nargs='+',
                  action="store", # optional because action defaults to "store"
                  dest="MapsFile",
                  help=".zip file with receptor and map",)
parser.add_argument("-o", "--outdir",
                  dest="outdir",
                  help="output directoryreceptor PDBQT file",)
parser.add_argument("-n", "--topn",
                  dest="NoOfOutputs",
                  help="Number of outputs",)
parser.add_argument("--nneighbors",
                  dest="nneigh",
                  help="Number of neighbors",)
parser.add_argument('-bc', '--beginCutoff', dest='bestCutoff', nargs="+", default="default",
                  help="-bc ratio RATIO_to_DEFAULT or -bc user C_CUTOFF O_CUTOFF H_CUTOFF")
parser.add_argument('-steps', '--steps', type=int, dest='steps',
                  help="-steps STEPS")
parser.add_argument('-ligandSize', '--ligandSize', type=int, dest='ligandSize',
                  help="-ligandSize ligandSize")
parser.add_argument('-p', "--pep", dest="pepscore", action="store_true",
                  default=False,
                  help="use peptide scoring function")
parser.add_argument("--numpysave",
                  action="store",
                  dest="NumpySave",
                  help="Save outputs as numpy",)
parser.add_argument("--flexRes",
                  action="store",
                  dest="flexRes",
                  help="get the flexible residues",)
parser.add_argument("--spacing",
                  action="store",
                  dest="spacing",
                  help="grid spacing",)
parser.add_argument("--boxcenter",
                  action="store",
                  dest="boxcenter",
                  help="center of gridbox",)
parser.add_argument("--boxdim",
                  action="store",
                  dest="boxdim",
                  help="dimension of gridbox in Angs",)

def writePDB(coords, potentials, atypes, filename):
    f = open(filename, 'w')
    for n, (xyz, pot, aty) in enumerate(zip(coords, potentials, atypes)):
        x,y,z = xyz
        f.write('ATOM%7d%4s   FP     1    %8.3f%8.3f%8.3f  1.0%8.3f    0.001 %s\n' %
                (n+1,aty[0], x, y, z, abs(pot), aty[0]))
                        #(n+1,aty[0]+'%s'%n, x, y, z, abs(pot), aty[0]))
    f.close()

def usage():
    print "usage: autosite -r <receptor.pdbqt> --spacing x.x -o <output Dir>"
    print "       If you already have C,OA,HD maps:"
    print "            autosite -r <receptor.pdbqt> -m *.map -o <output Dir>"
    
if len(sys.argv)==1:
    usage()
    sys.exit(1)
    
args = parser.parse_args()

receptor = Read(args.receptorFile)
name = os.path.basename(args.receptorFile).split('.')[0]

if args.outdir:
    outdir = args.outdir
else:
    # try to make a folder called after the receptor
    try:
        os.mkdir(name)
        outdir = name
    except OSError:
        raise RuntimeError, "ERROR unable to create folder %s, please use -o option to specify unused folder name for results"%name
    
folder = outdir

def myprint(str):
    sys.stdout.write(str+'\n')
    summaryFP.write(str+'\n')
        
try:
    summaryFP = open(outdir+'_AutoSiteSummary.log', 'w')
except OSError:
    print "ERROR: output file %s already exists, please remove or use the -o command line to provide unused output name"%(outdir+'_AutoSiteSummary.log')
    sys.exit(1)

if args.steps is None:
    revision = 0
else:
    revision = 1
    
t0 = time()
myprint( "#################################################################")
myprint( "# If you used AutoSite in your work, please cite:               #")
myprint( "#                                                               #")
myprint( "# Pradeep Anand Ravindranath and Michel F. Sanner               #")
myprint( "#                                                               #")
myprint( "# AutoSite: an automated approach for pseudoligands             #")
myprint( "# prediction - From ligand binding sites identification to      #")
myprint( "# predicting key ligand atoms                                   #")
myprint( "#                                                               #")
myprint( "# Bioinformatics (2016)                                         #")
myprint( "#                                                               #")
myprint( "# DOI: 10.1093/bioinformatics/btw367                            #")
myprint( "#                                                               #")
myprint( "# Please see http://adfr.scripps.edu/AutoDockFR/autosite.html   #")
myprint( "# for more information                                          #")
myprint( "#################################################################")
myprint( "")
myprint( 'AutoSite v%d.%d.%d started on %s a %s computer'%(
    AutoSite.__version__, revision, AutoSite.__built__, platform.node(), platform.platform()))
myprint( 'Date %s'%datetime.datetime.now().ctime())



coords = receptor._ag.getCoords()
radiiR = receptor._ag.getRadii()

if args.spacing and not args.MapsFile:
    spacing = float(args.spacing)
else:
    spacing = 1.0
if args.nneigh:
    nn = int(args.nneigh)
else:
    nn = 14
if args.boxcenter and args.boxdim:
    center = eval(args.boxcenter)
    size = [int(round(x/spacing)) for x in eval(args.boxdim)]
else:
    mini = numpy.min(coords, 0)
    maxi = numpy.max(coords, 0)
    center = 0.5*(mini+maxi)
    sizef = (maxi-mini) + 2*4
    size = [int(round(x/spacing)) for x in sizef]
#folder = '/tmp'


from AutoSite.compositePoints import CompositePoints
#import pdb;pdb.set_trace()
gc = CompositePoints(receptor, center, size, ['C', 'OA', 'HD'],
            spacing=spacing,
            smooth=0.5,
            flexibleResidues=[],
            folder=folder, atypesOnly=False, fp=True)

if args.MapsFile:
    for mapName in args.MapsFile:
        #import pdb;pdb.set_trace()
        atype = os.path.splitext(os.path.splitext(mapName)[0])[1][1:]
        if atype in ['C','OA', 'HD']:
            myprint('Reading map %s %s'%(atype, mapName))
            gc.mapFiles[atype]=mapName
            if atype == 'C':
                cmap = GridMap()
                cmap.loadFromMapFile('C', '', gc.mapFiles['C'])
                spacing = cmap.getDistBetweenGridPoints()
else:
    myprint('computing maps ....')
    status, msg = gc.run()
    if status!=0:
        myprint("ERROR: running autogrid failed in %s."%gc.folder)
        myprint("  %s"%gc._command)
        if msg.find("no closestH atom was found") >= 0:
            myprint("It seems that hydrogen atoms are missing in the receptor.\n ")
        else:
            myprint("%s. \n "%msg)
        sys.exit(1)

#set cutoff scan mode

if args.bestCutoff[0] == "ratio":
    ratio = float(args.bestCutoff[1])
    cmapFile = gc.mapFiles['C']
    cmap = GridMap()
    cmap.loadFromMapFile('C', '', cmapFile)
    omapFile = gc.mapFiles['OA']
    omap = GridMap()
    omap.loadFromMapFile('OA', '', omapFile)
    hmapFile = gc.mapFiles['HD']
    hmap = GridMap()
    hmap.loadFromMapFile('HD', '', hmapFile)
    cdata = cmap.getGridDataPy()#[::3, ::3, ::3]
    odata = omap.getGridDataPy()
    hdata = hmap.getGridDataPy()
    #cmin = numpy.min(cdata)/ratio
    #omin = numpy.min(odata)/ratio
    #hmin = numpy.min(hdata)/ratio
    cmin=-0.3*ratio
    omin=-0.66*ratio
    hmin=-0.5*ratio
elif args.bestCutoff[0] == "users":
    cmin = float(args.bestCutoff[1])
    omin = float(args.bestCutoff[2])
    hmin = float(args.bestCutoff[3])
else:
    args.bestCutoff[0] == "default"
    cmin = -0.30
    omin = -0.66
    hmin = -0.5
#else:
#    raise ValueError("setBox: ERROR bad cutoff scan mode, got %s"%args.bestCutoff)


if spacing != 1.0:
    scale = int(round(1.0/spacing))
else:
    scale = 1

from AutoSite.utils.clusterTPoints import DensityClustering

if args.steps is not None:
    steps=int(args.steps)
    print "Cutoffs scan starting from", cmin,omin,hmin, " to ", cmin/2,omin/2,hmin/2
    
    tt0=time()
    from AutoSite.compositePoints import CompositePoints
    dcl, headNode = gc.bestCutoffClustering(receptor,spacing=spacing,carbon_cutoff=cmin,oxygen_cutoff=omin,hydrogen_cutoff=hmin,nbSteps=steps)
    print "cutoff scan time:%8.2f seconds"%(time()-tt0)
else:
    print "No cutoff scan and use cutoffs", cmin,omin,hmin
    gc.getASPoints(carbon_cutoff=cmin, oxygen_cutoff=omin, hydrogen_cutoff=hmin)
    myprint('clustering high affinity points ... ')
    dcl = DensityClustering([spacing,spacing,spacing],neighborPts=nn)
    dcl.findClustersD(gc._indices)

#import pdb; pdb.set_trace()
#dcl._clusters = [x for x in dcl._clusters if len(x) >= 50*scale]

myprint('analysing %d clusters ... '%len(dcl._clusters))
from AutoSite.scoreClusters import scoreClusters

clusters, clProp = scoreClusters(receptor, dcl, gc,inflate=True, pepScore = args.pepscore)

if args.steps > 1 and args.ligandSize is not None:
    #import pdb;pdb.set_trace()
    headNode.updateAllNodes(clProp)
    finaloutput=headNode.getNodebySize(args.ligandSize)
    finalclusters = []
    finalclustersProp = []
    for clnode in finaloutput:
        finalclusters.append(clusters[clnode.id-1])
        finalclustersProp.append(clProp[clnode.id-1])


    clustersorted = sorted(finalclusters,key=lambda x:x[4],reverse = True)
    clPropsorted = sorted(finalclustersProp,key=lambda x:x[5], reverse = True)
elif args.steps > 1:
    headNode.updateAllNodes(clProp)
    clustersorted = clusters
    clPropsorted = clProp
    #rebuild the tree for inflated pockets
    #import pdb;pdb.set_trace()
    import json   
    with open('pockettree.txt', 'w') as outfile:
        json.dump(headNode.writeJSON(), outfile)
else:
    clustersorted = sorted(clusters,key=lambda x:x[4],reverse = True)
    clPropsorted = sorted(clProp,key=lambda x:x[5], reverse = True)

if args.NoOfOutputs:
    nbc = min(int(args.NoOfOutputs), len(clustersorted))
else:
    nbc = len(clustersorted)
    
myprint('saving %d clusters and associated feature points in folder %s'%(
    nbc, name))
clN = 1
summary = open('%s/%s_summary.csv'%(outdir,name),'w')
summary.write('%s,%s,%s,%s,%s,%s,%s\n'%('Cluster # | Energy | #points | Radius ofN','e','v','rg','epv','buriedness','v*buriedness^2/rg'))
myprint('\nclust.| Energy| # of |Rad. of | energy |   bns    |score   |Feature Points')
if args.pepscore:
    myprint('number|       |points|gyration|per vol.|buriedness|v*b^2/rg| D | A | Other')
else:
    myprint('number|       |points|gyration|per vol.|buriedness|v*b^1.5 | D | A | Other')
myprint('------+-------+------+--------+--------+----------+--------+---+---+------')


for cl, clp in zip(clustersorted, clPropsorted):
    if args.NoOfOutputs:
        if clN > int(args.NoOfOutputs): continue
    summary.write('%d,%f,%f,%f,%f,%f,%f\n'%(clN,clp[0],clp[1],clp[3],clp[2],clp[4],clp[5]))
    ## indices = cl[0]
    ## c = indices[numpy.where(cl[3]=='C')[0]]
    ## o = indices[numpy.where(cl[3]=='O')[0]]
    ## h = indices[numpy.where(cl[3]=='H')[0]]

    ccoords = cl[1][numpy.where(cl[3]=='C')[0]]
    ocoords = cl[1][numpy.where(cl[3]=='O')[0]]
    hcoords = cl[1][numpy.where(cl[3]=='H')[0]]
    cpot = cl[2][numpy.where(cl[3]=='C')[0]]
    opot = cl[2][numpy.where(cl[3]=='O')[0]]
    hpot = cl[2][numpy.where(cl[3]=='H')[0]]
    #import pdb;pdb.set_trace()
    fpcoords = []
    fppot = []
    fpaty = []
    ## ccoords = [x for i,x in enumerate(cl[1]) if cl[3][i] == 'C']
    ## cpot = [x for i,x in enumerate(cl[2]) if cl[3][i] == 'C']
    ## ocoords = [x for i,x in enumerate(cl[1]) if cl[3][i] == 'O']
    ## opot = [x for i,x in enumerate(cl[2]) if cl[3][i] == 'O']
    ## hcoords = [x for i,x in enumerate(cl[1]) if cl[3][i] == 'H']
    ## hpot = [x for i,x in enumerate(cl[2]) if cl[3][i] == 'H']
    from AutoSite.ASfeaturepoints import featurePts
    fp = featurePts(receptor,ccoords, cpot, ocoords, opot, hcoords, hpot, scale=scale)
    myprint('%5d %8.2f %5d %7.2f %8.2f     %.3f  %8.2f  %3d %3d  %3d'%(
        clN,clp[0],clp[1],clp[3],clp[2],clp[4],clp[5], len(fp.hfp), len(fp.ofp), len(fp.cfp)))
    if args.flexRes:
        fresdict = {}
        for x in fp.bsres:
            xsplit = x.split(':')
            fresdict.setdefault(xsplit[0], [])
            fresdict[xsplit[0]].append(xsplit[1])

        #flexResStr = ''
        flexibleResidues = []
        for k,v in fresdict.iteritems():
            #resstr =  "".join(str(x)+',' for x in v)
            #flexResStr = flexResStr+k+':'+resstr+';'
            reslst = []
            for res in v:
                if res[:3] not in ['ALA','GLY','PRO']:
                    reslst.append((res[:3],res[3:]))
            flexibleResidues.append((k,reslst))
        f = open('%s_fRes_%d.txt'%(name,clN))
        f.write("%s\n"%flexibleResidues)
        f.close()

    fpcoords.extend(fp.cfp); fpcoords.extend(fp.ofp); fpcoords.extend(fp.hfp)
    fppot.extend(fp.cfpp); fppot.extend(fp.ofpp); fppot.extend(fp.hfpp)
    fpaty.extend(['C']*len(fp.cfpp)); fpaty.extend(['O']*len(fp.ofpp)); fpaty.extend(['H']*len(fp.hfpp))
    writePDB(cl[1], cl[2],cl[3], '%s/%s_cl_%03d.pdb'%(outdir,name,clN))
    writePDB(fpcoords, fppot,fpaty, '%s/%s_fp_%03d.pdb'%(outdir,name,clN))
    if args.NumpySave:
        clnpname = '%s/%s_cl_%03d.npy'%(outdir,name,clN)
        fpnpname = '%s/%s_fp_%03d.npy'%(outdir,name,clN)
        numpy.save(clnpname, cl[1])
        numpy.save(fpnpname, fpcoords)
    clN+=1
summary.close()

dt = time()-t0
h,m,s = str(datetime.timedelta(seconds=dt)).split(':')
myprint( 'AutoSite identified and characterized %d clusters in %.2f seconds, i.e. %s hours %s minutes %s seconds '%(len(clustersorted), dt, h, m, s))
summaryFP.close()

