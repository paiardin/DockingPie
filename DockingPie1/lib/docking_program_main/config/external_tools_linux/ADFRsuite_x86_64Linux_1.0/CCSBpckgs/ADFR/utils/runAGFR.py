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
# Copyright: M. Sanner and TSRI 2016
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/ADFR/utils/runAGFR.py,v 1.28.2.8 2017/10/05 20:30:46 annao Exp $
#
# $Id: runAGFR.py,v 1.28.2.8 2017/10/05 20:30:46 annao Exp $
#
import numpy, tempfile, os, sys, shutil, platform, datetime, pickle
from time import time
from glob import glob
from math import ceil

from MolKit2 import Read
from prody.atomic.atom import Atom
from ADFR import checkLigandFile
from ADFR.utils.maps import flexResStr2flexRes
from ADFR.utils.MakeGrids import splitFlexRes
from ADFR.utils.addGradients import addGradientToMaps

from AutoSite.compositePoints import CompositePoints
from AutoSite.utils.clusterTPoints import DensityClustering
from AutoSite.scoreClusters import scoreClusters
from AutoSite.shrink import shrinkPocket

def saveATOMS(mol, filename, selection):
    # save the ATOM or HETATM records for the specified selection 
    # from the original file of mol
    toSave = {}
    for a in selection:
        toSave['%.3f,%.3f,%.3f'%(tuple(a.getCoords()))] = True
    f = open(mol.filename)
    lines = f.readlines()
    f.close()

    f = open(filename, 'w')
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            key = '%.3f,%.3f,%.3f'%(float(line[30:38]),float(line[38:46]),float(line[46:54]))
            if toSave.has_key(key):
                f.write(line)
    f.close()

errorCodes = {0:   """ready to compute grids""",
              100: """please load a receptor""",
              101: """please specify the docking box""",
              102: """box does not overlap receptor.""",
              103: """some flexible residues are outside the box""",
              104: """covalent bond atom(s) are outside the box""",
              105: """please specify a pocket in the "ligand binding pocket" section""",
              106: """none of the pocket points are inside the box""",
              107: """map types are not specified""",
              108: """atoms defining a chi angle not found in receptor""",
              109: """please specify atoms forming the covalent ligand attachment bond""",
              110: """One or more residues cannot be made flexible"""}

actionRecommendedGui = {0: """Compute grids""",
                        100: """load .pdbqt receptor file""",
                        101: """set the docking box.""",
                        102: """change the docking box parameters.""",
                        103: """specify different flexible residues or change the docking box parameters""",
                        104: """ """,
                        105:"""compute pockets (run AutoSite)""",
                        106: """change the docking box parameters """,
                        107: """select map types""",
                        108: """ """,
                        109: """Select covalent bond atoms""",
                        110: """Check flexible residues list"""}                  

actionRecommended = {0: """computeGrids(outFile, flexResStr, spacing, background=False, indent='')""",
                     100: """loadReceptor(filename)""",
                     101: """setBox(mode, padding, spacing)\n ---mode can be:\n ['receptor']: smallest box encompassing the entire receptor (default unless ligand is present).\n ['ligand']  : smallest box encompassing a specified ligand \n ['fill']    : smallest box encompassing fill points\n ['residues', chid, resnamesResnums]: smallest box encompassing the specified residues\n ['user', (cx, cy, cz sx), (sx, sy, sz)]: box centered at (cx, cy, cz) of size sx, sy, sz (units Angstroms) """,
                     102: """ Use setBox(mode, padding, spacing)\n ---mode can be:\n ['receptor']: smallest box encompassing the entire receptor (default unless ligand is present).\n ['ligand']  : smallest box encompassing a specified ligand \n ['fill']    : smallest box encompassing fill points\n ['residues', chid, resnamesResnums]: smallest box encompassing the specified residues\n ['user', (cx, cy, cz sx), (sx, sy, sz)]: box centered at (cx, cy, cz) of size sx, sy, sz (units Angstroms) """,
                     103: """setFlexResidues(flexresStr)\n flexresStr example: 'A:ILE10,GLU34'""",
                     104: """ """,
                     105: """runAutoSite(flexResStr=None, smooth=0.5, cutoff=10)""",
                     106: """setBox(mode, padding, spacing)\n ---mode can be:\n ['receptor']: smallest box encompassing the entire receptor (default unless ligand is present).\n ['ligand']  : smallest box encompassing a specified ligand \n ['fill']    : smallest box encompassing fill points\n ['residues', chid, resnamesResnums]: smallest box encompassing the specified residues\n ['user', (cx, cy, cz sx), (sx, sy, sz)]: box centered at (cx, cy, cz) of size sx, sy, sz (units Angstroms) """,
                     107: """setMapTypes(maptype ).\n---maptype can be  'all' or 'ligand' or [list of atom types]""",
                     108: """ """,
                     109: """ """,
                     110: """Check flexible residues list"""
                     }

class runAGFR:
    """
    class to run AGFR from command line
    """

    def myprint(self, str, newline=True):
        sys.stdout.write(str)
        if self.summaryFP:
            self.summaryFP.write(str)
        if newline:
            sys.stdout.write('\n')
            if self.summaryFP:
                self.summaryFP.write('\n')

    def __init__(self):
        self.summaryFP = None
        self.receptor = None
        self.ligand = None
        self.boxCenter = None
        self.boxLengths = None
        self.boxSize = None #(number of grid points)
        self.spacing = None
        self.pockets = [] # will be a list of fill Points corresponding each
                          # to a docking pocket
                          
        self.atypes = []
        self.data = {}
        self.covalentBond = False
        self.fillPoints = []
        self.flexResAtoms = []
        self.covalentBondToExclude = []
        self.cmdOptions = None
        self.autoSite2 = False # Use atoSite 1.0 by default. True - to use AutoSite2 
        self._wMapEntropy = -0.2 
        self._wMapWeight = 0.6
        self.cutOffValue = None
        
    def loadReceptor(self, filename):
        # make sure receptor is not a ligand and if so load it
        if checkLigandFile(filename):
            self.myprint("ERROR: the file %s contains a torsion tree indicating this is a ligand."%filename)
            raise ValueError("the file %s contains a torsion tree indicating this is a ligand."%filename)
        self.setReceptor(Read(filename))

    def setReceptor(self, mol):
        self.receptor = mol
        self.data['inputReceptor'] = os.path.basename(mol.filename)

    def loadLigand(self, filename):
        # check if a ligand is provided and loaded if so
        if not checkLigandFile(filename):
            self.myprint("ERROR: the file %s does not contain a torsion tree indicating this is a receptor."%filename)
            raise ValueError("the file %s does not contain a torsion tree indicating this is a receptor."%filename)
        self.setLigand(Read(filename))

    def setLigand(self, mol):
        self.ligand = mol
        self.data['inputLigand'] = os.path.basename(mol.filename)

    def loadFillPoints(self, filenames):
        for filename in filenames:
            self.addFillPoints(numpy.load(filename))

    def addFillPoints(self, points):
        self.fillPoints.append(points)

    def setPadding(self, value):
        if self.boxSize is not None:
            x,y,z = self.boxLengths - 2*self.padding
            boxLengths = [x+2*value, y+2*value, z+2*value]
            s = []
            for x in boxLengths:
                n = int(ceil(x/self.spacing)) # requires number of grid points
                if n%2==1:
                    n+=1 # make it even
                s.append(n)
            self.boxSize = s
            self.boxLengths = numpy.array(boxLengths)
            self.data['boxLengths'] = self.boxLengths
        self.padding = value
        self.data['boxPadding'] = value

    def setSpacing(self, spacing):
        self.spacing = spacing
        self.data['spacing'] = spacing
        if self.boxSize is not None:
            s = []
            for x in self.boxLengths:
                n = int(ceil(x/spacing)) # requires number of grid points
                if n%2==1:
                    n+=1 # make it even
                s.append(n)
            self.boxSize = s
            self.data['boxSize'] = self.boxSize

    def setBoxCenter(self, value):
        assert len(value) == 3
        self.boxCenter = numpy.array(value)
        self.data['boxCenter'] = self.boxCenter
        
    def setBox(self, mode, padding, spacing):
        #print "SET BOX: MODE", mode, padding, spacing
        self.spacing = spacing
        userFlag = False
        coords = None
        self.data['boxMode'] = mode
        if mode[0] == "user" and mode[1] in ["receptor","ligand","fill","residues"]:  # User mode with specific center mode
            userFlag = True
            mode[0] = mode[1] # reset mode[0] to "receptor","ligand","fill" or "residues"
                              # to get coords of the atoms and use them to compute box center
            # set sides of the box with user specified values
            self.boxLengths = numpy.array([float(x) for x in mode[2:5]])
        if mode[0]=="receptor":
            coords = self.receptor._ag.getCoords()
        elif mode[0]=="ligand":
            coords = self.ligand._ag.getCoords()
        elif mode[0]=="fill":
            coords = self.fillPoints
        elif mode[0]=="residues":
            flexRes = flexResStr2flexRes(mode[1])
            receptorAtoms, sideChainAtoms = splitFlexRes(self.receptor, flexRes,
                                                         exclude='')
            coords = sideChainAtoms.getCoords()
            #print "RES COORDS", coords
        elif mode[0]=="user":
            # use user specified center of the box 
            self.boxCenter = numpy.array([float(x) for x in mode[1:4]])
            # set sides of the box with user specified values
            self.boxLengths = numpy.array([float(x) for x in mode[4:7]])
            userFlag = True
        else:
            self.myprint( 'ERROR: bad mode expected receptor, ligand, fill, residues, or user, got %s'%mode)
            raise ValueError("setBox: ERROR bad mode expected receptor, ligand, fill, residues, or user, got %s"%mode)
        
        if userFlag: #  User mode with specific center mode and 3 values to set box size
            if coords is not None: # center of the box is computed from the coords of receptor, ligand, fillPoints or flex residues
                mini = numpy.min(coords, 0)
                maxi = numpy.max(coords, 0)
                self.boxCenter = 0.5*(mini+maxi)
            # Compute the size of the box
            self.padding = 0.
            s = []
            for x in self.boxLengths:
                n = int(ceil(x/spacing)) # requires number of grid points
                if n%2==1:
                    n+=1 # make it even
                s.append(n)
            self.boxSize = s
            self.data['boxPadding'] = padding
            self.data['boxCenter'] = self.boxCenter
            self.data['boxLengths'] = self.boxLengths
            self.data['boxSize'] = self.boxSize
            self.data['spacing'] = spacing
        else:
            self.setBoxForCoords(coords, padding, spacing)

    def setBoxForCoords(self, coords, padding=None, spacing=None):
        if padding is None:
            padding = self.padding
        if spacing is None:
            spacing = self.spacing
        mini = numpy.min(coords, 0)
        maxi = numpy.max(coords, 0)
        self.boxCenter = 0.5*(mini+maxi)
        boxLengths = (maxi-mini) + 2*padding # needed length
        # size is th number of grid intervals to cover needed length
        # + 1 to make it grid points
        s = []
        for x in boxLengths:
            n = int(ceil(x/spacing)) # requires number of grid points
            if n%2==1:
                n+=1 # make it even
            s.append(n)
        self.boxSize = s
        self.padding = padding
        self.spacing = spacing
        self.boxLengths = numpy.array([x*spacing for x in self.boxSize])
        
        self.data['boxPadding'] = padding
        self.data['boxCenter'] = self.boxCenter
        self.data['boxLengths'] = self.boxLengths
        self.data['boxSize'] = self.boxSize
        self.data['spacing'] = spacing

    def setCovalentDocking(self, torsionAtIndex,
                           covBondAtIndex, #serial indices of two covalent bond atoms
                           covResStr=None, toRemoveAtoms=None):
        #print "setCovalentDocking:", torsionAtIndex, covBondAtIndex
        id1 = self.receptor.select('serial %d'%covBondAtIndex[0]).getIndices()[0]
        at1 = Atom(self.receptor._ag, id1, 0)
        id2 = self.receptor.select('serial %d'%covBondAtIndex[1]).getIndices()[0]
        at2 = Atom(self.receptor._ag, id2, 0)
        d = {}
        if covResStr is not None:
            for chid, resnums in  flexResStr2flexRes(covResStr):
                d[chid] = [x[1] for x in resnums]
        if not toRemoveAtoms:
            toRemoveAtoms = self.receptor.subTree(at1, at2, limitTo=d)
        id3 = self.receptor.select('serial %d'%torsionAtIndex).getIndices()[0]
        at3 = Atom(self.receptor._ag, id3, 0)
        ## self.data['covalentBondAtom1'] = '%s:%s%d:%s'%(at1.getChid(), at1.getResname(),
        ##                                                at1.getResnum(), at1.getName())        
        ## self.data['covalentBondAtom2'] = '%s:%s%d:%s'%(at2.getChid(), at2.getResname(),
        ##                                                at2.getResnum(), at2.getName())        
        ## self.data['covalentBondTorsionAtom'] = '%s:%s%d:%s (%d)'%(
        ##     at3.getChid(), at3.getResname(), at3.getResnum(), at3.getName(), torsionAtIndex)
        self.data['covalentBondAtom1'] = self.receptor.atomFullName(at1)
        self.data['covalentBondAtom2'] = self.receptor.atomFullName(at2)
        self.data['covalentBondTorsionAtom'] = '%s (%d)'%( self.receptor.atomFullName(at3), torsionAtIndex)
        
        self.data['covalentAtomsCoords'] = at1.getCoords().tolist()+ at2.getCoords().tolist()+ \
                                           at3.getCoords().tolist()
        self.data['covalentBond'] = covBondAtIndex
        self.covalentBond = True
        self.data['covalentRes'] = covResStr
        self.covalentBondToExclude = toRemoveAtoms

    def getAllADatomTypes(self):
        from ADFRcc import getFFParameters
        parameters = getFFParameters()
        self.ADAtomTypes = {}
        #import pdb;pdb.set_trace()
        # MS. March 2016 restrict types to AD4.2 types for now
        AD42atomTypes = {}.fromkeys(
            ['H', 'HD', 'HS', 'C', 'A', 'N', 'NA', 'NS', 'OA', 'OS', 'F',
             'Mg', 'MG', 'P', 'SA', 'S', 'Cl', 'CL', 'Ca', 'CA', 'Mn', 'MN',
             'Fe', 'FE', 'Zn', 'ZN', 'Br', 'BR', 'I', 'Z', 'G', 'GA', 'J', 'Q'])
        for i in range(parameters.numAtomTypes):
            atype = parameters.getAtomTypeByIndex(i).atomTypeName
            if atype in AD42atomTypes:
                self.ADAtomTypes[atype] = True

    def setMapTypes(self, mapTypes):
        if mapTypes=='all':
            self.getAllADatomTypes()
            self.atypes = self.ADAtomTypes.keys()
        elif mapTypes=='ligand':
            if self.ligand is None:
                self.myprint("WARNING: requesting map types from ligand atom types when no ligand was given. Using all map types.")
                #raise ValueError("ERROR: requesting map types from ligand atom types when no ligand was given. ")
                #sys.exit(1)
                self.getAllADatomTypes()
                self.atypes = self.ADAtomTypes.keys()
                mapTypes = 'all'
            else:
                self.atypes = numpy.unique(self.ligand._ag.getData('AD_element'))
        elif isinstance(mapTypes, list):
            self.atypes = mapTypes
        for tt in ["OA", "HD"]:  #make sure these map types are always in the list of selected maptypes.
            # the OA and HD maps will be used to create water map file in the trg file.
            if tt not in self.atypes:
                self.atypes = list(self.atypes)
                self.atypes.append(tt)
        self.data['mapTypes'] = self.atypes
        return mapTypes

    def setFlexResidues(self, flexresStr):
        self.data['flexResStr'] = flexresStr
        #if flexresStr is not None:
        #    receptorAtoms, sideChainAtoms = splitFlexRes(self.receptor, flexResStr2flexRes(flexresStr))
        #    self.flexResAtoms = sideChainAtoms

    ##
    ## grids computing methods based on AutoGrid4
    ##
    def _computeGrids(self, center, size, spacing, atypes, smooth=0.5,
                      flexResStr=None, folder='.', atypesOnly=False,
                      background=False, fp=False, outlev=1):
        # flexResStr is a string like this "A:ILE10,GLU34"
        flexRes = flexResStr2flexRes(flexResStr)
        #print "FLEXRES", flexResStr, flexRes, len(flexRes)
        gc = CompositePoints(
            self.receptor, center, size, atypes, spacing=spacing,
            smooth=smooth, flexibleResidues=flexRes,
            folder=folder, atypesOnly=atypesOnly, fp=fp,
            covalentBondToExclude=self.covalentBondToExclude, outlev=outlev)
        status, msg = gc.run(background=background)

        if status!=0:
            self.myprint("ERROR: running autogrid failed in %s."%gc.folder)
            self.myprint("  %s"%gc._command)
            if msg.find("no closestH atom was found") >= 0:
                self.myprint("It seems that hydrogen atoms are missing in the receptor.\n ")
            else:
                self.myprint("%s. \n "%msg)
            #raise RuntimeError("ERROR: running autogrid failed in %s"%gc.folder)
        # changed headers
        #if flexRes is not None:
        #    gc.addFlexRecHeader('FLEXRES "%s"'%flexResStr)
        return gc, status

    ##
    ## AutoSite computing methods
    ##
    def AutoSiteFill(self, center, length, spacing, flexResStr=None,
                     smooth=0.5, background=False, outlev=1):
        coords = self.receptor._ag.getCoords()
        radiiR = self.receptor._ag.getRadii()
        
        size = [int(ceil(x/spacing)) for x in length]
        #print "SIZE", size
        self.tmpFolder = tempfile.mktemp()
        os.mkdir(self.tmpFolder)
        gc, status = self._computeGrids(
            center, size, atypes=['C','OA','HD'], spacing=spacing,
            smooth=smooth, flexResStr=flexResStr, folder=self.tmpFolder,
            atypesOnly=True, fp=True, background=background, outlev=outlev)
        return gc, status

    def runAutoSite(self, flexResStr=None, smooth=0.5,
                    spacing=1.0, verbose=False, background=False, outlev=1):
        #tpoints = None
        if verbose:
            self.myprint( '\nidentifying pockets using AutoSite ....')
        if flexResStr is None:
            if self.data.has_key('flexResStr'):
                flexResStr = self.data['flexResStr']
        gc, status = self.AutoSiteFill(
            self.boxCenter, self.boxLengths, spacing,
            flexResStr=flexResStr, background=background, smooth=smooth, outlev=outlev)
        # if background is True, the calling programm should wait until the process
        # finishes
        return gc, status

    def afterAutoSite2(self, gc, spacing=1.0, cutoff=10, ligandSize=500, pepScore=False, verbose=False):
        dcl, headNode = gc.bestCutoffClustering(self.receptor,spacing=spacing,carbon_cutoff=-0.36,oxygen_cutoff=-0.792,hydrogen_cutoff=-0.6,nbSteps=4)
        clusters, clProp = scoreClusters(self.receptor, dcl, gc,inflate=True, pepScore = pepScore)
        nbp = numpy.sum(dcl._clen)
        if nbp < 50:
            if verbose:
                self.myprint('    WARNING: found %d pocket(s) adding up to %d fill points.\n'%(len(dcl._clusters),nbp))
                self.myprint('    WARNING: ignoring density clustering and using all %d fill points.\n'%len(gc._indices))
            dcl._clusters = [range(len(gc._indices))]

        headNode.updateAllNodes(clProp)
        finaloutput=headNode.getNodebySize(ligandSize)
        finalclusters = []
        finalclustersProp = []
        for clnode in finaloutput:
            finalclusters.append(clusters[clnode.id-1])
            finalclustersProp.append(clProp[clnode.id-1])


        clustersorted = sorted(finalclusters,key=lambda x:x[4],reverse = True)
        clPropsorted = sorted(finalclustersProp,key=lambda x:x[5], reverse = True)
  
        if verbose:
            self.myprint('    found %d pocket(s)\n'%len(clustersorted))
            self.myprint('    pocket|  energy | # of |Rad. of | energy |   bns    | score  ')
            if pepScore:
                self.myprint('    number|         |points|gyration|per vol.|buriedness|v*b^1.5 ')
            else:
                self.myprint('    number|         |points|gyration|per vol.|buriedness|v*b^2/rg')
            self.myprint('    ------+---------+------+--------+--------+----------+---------')
            n = 0
            for cl, clp in zip(clustersorted, clPropsorted):
                n += 1
                self.myprint('     %4d %9.2f %5d %7.2f   %7.2f    %6.2f    %7.2f'%(
                    n,clp[0],clp[1],clp[3],clp[2],clp[4],clp[5]))
        return  clustersorted, clPropsorted, dcl


    def afterAutoSite(self, gc, spacing=1.0, cutoff=10, verbose=False):
        
        gc.getASPoints()
        # save list of all indices before clustering, in case clustering removes too much
        allIndices = gc._indices[:]
        dcl = DensityClustering([spacing,spacing,spacing], neighborPts=14)
        dcl.findClustersD(gc._indices,  cVolcut=cutoff)
        nbp = numpy.sum(dcl._clen)
        if nbp < 50:
            if verbose:
                self.myprint('    WARNING: found %d pocket(s) adding up to %d fill points.\n'%(len(dcl._clusters),nbp))
                self.myprint('    WARNING: ignoring density clustering and using all %d fill points.\n'%len(gc._indices))
            dcl._clusters = [range(len(gc._indices))]
        clusters, clProp = scoreClusters(self.receptor, dcl, gc)
        clustersorted = sorted(clusters,key=lambda x:x[4],reverse=True)
        clPropsorted = sorted(clProp,key=lambda x:x[5],reverse=True)
        if verbose:
            self.myprint('    found %d pocket(s)\n'%len(clustersorted))
            self.myprint('    pocket|  energy | # of |Rad. of | energy |   bns    | score  ')
            self.myprint('    number|         |points|gyration|per vol.|buriedness|v*b^2/rg')
            self.myprint('    ------+---------+------+--------+--------+----------+---------')
            n = 0
            for cl, clp in zip(clustersorted, clPropsorted):
                n += 1
                self.myprint('     %4d %9.2f %5d %7.2f   %7.2f    %6.2f    %7.2f'%(
                    n,clp[0],clp[1],clp[3],clp[2],clp[4],clp[5]))
        return  clustersorted, clPropsorted, dcl


    def getFillPoints(self, pocketMode, cutoff, clustersorted, verbose=False):
        # This function returns a list of fillPoints for specified poketMode:
        # "all" : fillPoins is a list containing one numpy array of all points
        #         merged from all clusters   ;
        # "best": fillPoins is a list containing one numpy array of points from clustersorted[0]; 
        # "forEach": fillPoints is a list containing numpy arrays of points per each cluster containing more than cutoff number of points;
        # "forTop": fillPoints containes N numpy arrays of points from top N clusters (N = cutoff)
        fillPoints = []
        tpoints = []
        
        if verbose:
            if pocketMode =='all':
            # merge points from all clusters with more than cutoff points
                self.myprint( '    merging clusters ...')
            elif pocketMode =='best':
                self.myprint( '    using best score cluster with %d points'%len(clustersorted[0][1]))
        if pocketMode =='all':
            for cl in clustersorted:
                if len(cl[0])<cutoff:
                    continue
                if self.autoSite2:
                    #Shrink pocket to 1/5 of original size to get better tranlation points 
                    tpoints.extend(shrinkPocket(cl[1],finalsize=len(cl[1])/5))
                else:
                    tpoints.extend(cl[1])
            fillPoints.append(tpoints)
        elif pocketMode == 'best':
            # take coordinates for all points indices in cluster 0
            if self.autoSite2:
                tpoints = shrinkPocket(clustersorted[0][1], finalsize=len(clustersorted[0][1])/5)
            else:
                tpoints = clustersorted[0][1]
            fillPoints.append(tpoints)
        else:
            for cn, cl in enumerate(clustersorted):
                if pocketMode =='forEach' and len(cl[1])<cutoff:
                    continue
                elif pocketMode =='forTop' and cn>=cutoff:
                    continue
                if self.autoSite2:
                    #Shrink pocket to 1/5 of original size to get better tranlation points 
                     fillPoints.append(shrinkPocket(cl[1], finalsize=len(cl[1])/5))
                else:
                    fillPoints.append(cl[1])
        self.tpoints = tpoints
        return fillPoints

    def setFillPoints(self, tpoints):
        assert isinstance(tpoints, (numpy.ndarray, list))
        #self.fillPoints, outside = self.pointsInBox(tpoints)
        self.fillPoints = tpoints

    def pointsInBox(self, pts):
        inside = []
        outside = []
        llx, lly, llz = self.boxCenter - self.boxLengths/2
        urx, ury, urz = self.boxCenter + self.boxLengths/2
        for x,y,z in pts:
            if x>llx and x<urx and y>lly and y<ury and z>llz and z<urz:
                inside.append( (x,y,z) )
            else:
                outside.append( (x,y,z) )
        return inside, outside

    def computeGrids(self, outFile, flexResStr, spacing, background=False,
                     indent='', addGradients=False):
        # create a folder called outfile in which we will compute the maps
        # add TPoints and receptor and zip up as a target object
        outFile = os.path.splitext(outFile)[0]
        destinationFolderPath, destinationFolder = os.path.split(outFile)
        if destinationFolderPath=='':
            destinationFolderPath = '.'
        self.destinationFolderPath = destinationFolderPath
        self.destinationFolder = destinationFolder
        newGridsFolder = os.path.abspath(os.path.join(destinationFolderPath, destinationFolder))
        if os.path.exists(newGridsFolder):
            shutil.rmtree(newGridsFolder)
        try:
            os.makedirs(newGridsFolder)
        except Exception as inst:
            msg = "Failed to make %s folder: %s " % (newGridsFolder, inst.strerror)
            self.myprint (msg)    # the exception instance
            #print inst.filename, inst.message, inst.strerror
            return None, msg
        # wrong to add 1 here because it was already factored in self.boxLengths in setBox
        #size = [int(ceil(x/spacing))+1 for x in self.boxLengths]
        size = self.boxSize
        # compute the grids
        t0 = time()
        #print "IN COMPUTE GRIDS", "box center", self.boxCenter, "size", size, "spacing", spacing
        gc, status = self._computeGrids(self.boxCenter, size, spacing, self.atypes, flexResStr=flexResStr,
                                        background=background, folder=newGridsFolder, outlev=2)
        if status==0 and not background:
            self.myprint(indent+"maps computed in %.2f (sec)"%(time()-t0))
            self.generateTrgFile(gc, newGridsFolder, flexResStr,
                                 addGradients=addGradients)
        return gc, status

    def generateTrgFile(self, gc, gridFolder, flexResStr, indent="", addGradients=False, logFileName=None):
        if len(gc.flexRecAtoms):
            self.myprint(indent+"the following %d flexible receptor atoms did not contribute to the grid calculation:"%len(gc.flexRecAtoms))
            hv = gc.flexRecAtoms.getHierView()
            for res in hv.iterResidues():
                self.myprint(indent+"  %s:%s%d:"% (
                    res.getChid(), res.getResname(), res.getResnum()),
                    newline=False)
                for a in res.iterAtoms():
                    self.myprint("%s,"%a.getName(), newline=False)
                self.myprint("")
                self.data['flexRecFile'] = 'flexRec.pdbqt'
                saveATOMS(self.receptor, os.path.join(gridFolder, 'flexRec.pdbqt'),
                          gc.flexRecAtoms)
        else:
            self.data['flexRecFile'] = ''
        if len(gc.covalentLigAtoms):
            self.myprint(indent+"the following %d covalent ligand atoms did not contribute to the grid calculation:"%len(gc.covalentLigAtoms))
            hv = gc.covalentLigAtoms.getHierView()
            #for a in gc.covalentLigAtoms:
            for res in hv.iterResidues():
                self.myprint(indent+"  %s:%s%d:"% (
                    res.getChid(), res.getResname(), res.getResnum()),
                    newline=False)
                for a in res.iterAtoms():
                    self.myprint("%s,"%a.getName(), newline=False)
                self.myprint("")
            self.data['covalentLigandFile'] = 'covalenLig.pdbqt'
            saveATOMS(self.receptor, os.path.join(gridFolder, 'covalenLig.pdbqt'),
                      gc.covalentLigAtoms)
            self.data['covalentLigandAtomIndices'] = [
                a.getIndex() for a in gc.covalentLigAtoms]
        else:
            self.data['covalentLigandFile'] = ''
            self.data['covalentLigandAtomIndices'] = None
        # fix maps
        mapFiles = glob(os.path.join(gridFolder, '*.map'))
        mtypes = []
        for name in mapFiles:
            mtypes.append(os.path.basename(name).split('.')[-2])
        self.data['mapTypes'] = mtypes
        
        if flexResStr is not None:
            flexStr = flexResStr2flexRes(flexResStr)
            receptorAtoms, sideChainAtoms = splitFlexRes(self.receptor, flexStr)
        else:
            receptorAtoms = self.receptor.select()
        t0 = time()
        if addGradients:
            t0 = time()
            self.myprint(indent+"Adding gradient to maps ...")
            
            addGradientToMaps(mapFiles, mtypes, self.data['spacing'], self.cutOffValue, errorCut=0.01, logFileName=logFileName)
            self.myprint(indent+"done adding gradient to maps %.2f (sec)"%(time()-t0))
        self.data['mapGradients'] = addGradients
        self.data['gradCutOff'] = self.cutOffValue
        
        if "OA" in mtypes and "HD" in mtypes:
            self.makeWmap( os.path.split(mapFiles[0])[0], self.data['spacing'], name="W", ENTROPY=self._wMapEntropy, weight=self._wMapWeight) #os.path.dirname(mapFiles[0]), self.data['spacing'])
            self.data['mapTypes'].append("W")
        # save translation points
        if not self.covalentBond:
            filename = os.path.join(gridFolder, 'translationPoints.npy')
            numpy.save(filename, self.fillPoints)
            self.data['fillPointsFile'] = filename
        self.data['nbFillPoints'] = len(self.fillPoints)
        if not self.data.has_key('date'):
            self.data['date'] = datetime.datetime.now().ctime()
        if not self.data.has_key('platform'):
            self.data['platform'] = platform.platform()
        if not self.data.has_key('node'):
            self.data['node'] = platform.node()
        #add AGFR version number to the datafile
        from Support import version
        self.data["AGFRVERSION"] = version.__version__
        # make sure we have correct box data:
        self.data['boxPadding'] = self.padding
        self.data['boxCenter'] = self.boxCenter
        self.data['boxLengths'] = self.boxLengths
        self.data['boxSize'] = self.boxSize
        self.data['spacing'] = self.spacing
        
        with open(os.path.join(gridFolder, 'data.pkl'), 'wb') as f:
            pickle.dump(self.data, f)
        if self.cmdOptions:
            from ADFR.utils.optParser import makeConfigFile
            cfgfile = os.path.join(gridFolder, self.destinationFolder+".cfg")
            makeConfigFile(self.cmdOptions, cfgfile)

        self.myprint(indent+"making target file %s ..."%(
            self.destinationFolder+'.trg',), newline=False)
        cwd = os.getcwd()
        # find out the path to the gridFolder:
        gridFolderDir = os.path.abspath(os.path.dirname(gridFolder))
        gridFolderName = os.path.basename(gridFolder)
        # we will zip the gridFolder in the gridFolderDir (this is to avoid the error when creating the archive in the
        # cwd on Windows that had no write permission)
        os.chdir(gridFolderDir)
        shutil.make_archive(self.destinationFolder, 'zip', '.',
                            gridFolderName)
        os.chdir(cwd)
        if not os.path.exists(self.destinationFolderPath):
           os.mkdir(self.destinationFolderPath)
        shutil.move(os.path.join(gridFolderDir, self.destinationFolder+'.zip'),  os.path.join(self.destinationFolderPath, self.destinationFolder+'.trg'))
        shutil.rmtree(gridFolder)
        self.myprint(indent+"done.")
        
    def setAutoSiteVersion(self, version, ligandSize=None, pepScore=None):
        self.data["AutoSiteVersion"] = version
        if ligandSize is not None:
            self.data["ligandSize"] = ligandSize
        if pepScore is not None:
            self.data["pepScore"] = pepScore

    def setWmapParams(self, weight=None, ENTROPY=None):
        if ENTROPY is not None:
            self._wMapEntropy = ENTROPY
        if weight is not None:
            self._wMapWeight = weight
        #print "in  setWmapParams:", "weight:", self._wMapWeight, "entropy: ", self._wMapEntropy

    def makeWmap(self, mapsFolder,  spacing, weight=0.6, ENTROPY=-0.2, name="W"):
        """ Combine OA and HD maps into W map.
        """
        from ADFRcc.adfr import GridMap
        OAmap = GridMap()
        OAmap.loadFromMapFile('OA', mapsFolder, 'rigidReceptor.%s.map'%'OA')
        OAdata = OAmap.getGridDataPy()

        HDmap = GridMap()
        HDmap.loadFromMapFile('HD', mapsFolder, 'rigidReceptor.%s.map'%'HD')
        HDdata = HDmap.getGridDataPy()
        wData = self.best(OAdata, HDdata, weight, ENTROPY)
        spacing  = HDmap.getDistBetweenGridPoints()
        from Volume.Grid3D import Grid3DF
        grid = Grid3DF(wData.astype('f'), HDmap.getOriginPy(),
               (spacing, spacing, spacing),
               {'GRID_PARAMETER_FILE':'None',
                'GRID_DATA_FILE':'None',
                'MACROMOLECULE':'None'})
        from Volume.IO.AutoGridWriter import WriteAutoGrid
        writer = WriteAutoGrid()
        writer.write(grid, os.path.join(mapsFolder, 'rigidReceptor.%s.map'%name))
        self.data['wMapEntropy'] = ENTROPY
        self.data['wMapWeight'] = weight
        
    def best(self, OAdata, HDdata, weight, ENTROPY):
        """
        - Return the best value between two map values
        - If any positive values are found, 0 is returned by default
        - by default it should be first=OA, second=HD  
        """
        nx,ny,nz = OAdata.shape
        wData = numpy.zeros(OAdata.shape, 'f')
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    oa = OAdata[i][j][k]
                    hd = HDdata[i][j][k]
                    if oa > 0 or hd > 0: 
                        wData[i][j][k] = ENTROPY
                    else:
                        if oa < hd:
                            wData[i][j][k] = OAdata[i][j][k] * weight
                        else:
                            wData[i][j][k] = HDdata[i][j][k] * weight
        return wData
        

    def saveCmdOptions(self, kw):
        self.cmdOptions = {}
        for opt, val in kw.items():
            if val is not None:
               self.cmdOptions[opt] = val 
   
    def __call__(self, *args, **kw):

        self.myprint( "#################################################################")
        self.myprint( "# If you used AGFR in your work, please cite:                   #")
        self.myprint( "#                                                               #")
        self.myprint( "# P.A. Ravindranath S. Forli, D.S. Goodsell, A.J. Olson and     #")
        self.myprint( "# M.F. Sanner                                                   #")
        self.myprint( "# AutoDockFR: Advances in Protein-Ligand Docking with           #")
        self.myprint( "# Explicitly Specified Binding Site Flexibility                 #")
        self.myprint( "# PLoS Comput Biol 11(12): e1004586                             #")
        self.myprint( "# DOI:10.1371/journal.pcbi.1004586                              #")
        self.myprint( "#                                                               #")
        self.myprint( "# P. Ananad Ravindranath and M.F. Sanner                        #")
        self.myprint( "# AutoSite: an automated approach for pseudoligands prediction  #")
        self.myprint( "# - From ligand binding sites identification to predicting key  #")
        self.myprint( "# ligand atoms                                                  #")
        self.myprint( "# Bioinformatics (2016)                                         #")
        self.myprint( "# DOI:10.1093/bioinformatics/btw367                             #")
        self.myprint( "#                                                               #")
        self.myprint( "# Please see http://adfr.scripps.edu for more information.      #")
        self.myprint( "#################################################################")
        self.myprint( "")
        self.myprint( 'Computing grids on %s a %s computer'%(platform.node(), platform.platform(), ))
        self.myprint( 'Date %s\n'%datetime.datetime.now().ctime())
        self.data['date'] = datetime.datetime.now().ctime()
        self.data['platform'] = platform.platform()
        self.data['node'] = platform.node()
        
        t00 = t0 = time()
        #import pdb; pdb.set_trace()
        # kw ['flexres', 'covalentRes', 'boxMode', 'covalentBondTorsionAtom', 'pocketMode', 'spacing', 'covalentBond', 'ligandFile', 'padding', 'mapTypes', 'smooth', 'pocketCutoff', 'receptorFile', 'outputFile']
        #print "OPTIONS", kw
        self.saveCmdOptions(kw)
        filename = kw.get('outputFile', None)
        if filename:
            filePath = os.path.split(filename)[0]
            filenameBase = os.path.splitext(filename)[0]
            if filePath and not os.path.exists(filePath):
                os.makedirs(filePath)
            self.summaryFP = open(filenameBase+'.log', 'w')
        else:
            self.summaryFP = None

        self.setPadding(kw['padding']) # sets self.padding
        self.setSpacing( kw['spacing']) # sets self.spacing
        self.receptorGradient = kw.get("receptorGradient", True)
        self.cutOffValue = kw.get("recGradVolCut", None)
        if self.receptorGradient:
            if self.cutOffValue is None:
                self.cutOffValue = -1
        # check and get the receptor
        self.myprint( 'loading receptor: %s'%kw['receptorFile'])
        self.loadReceptor(kw['receptorFile'])  # reads receptor molecule from file and assigns Molecule instance to self.receptor 

        t0 = time()
        ligFilename = kw.get('ligandFile', None)
        if ligFilename:
            self.myprint( 'loading ligand: %s\n'%ligFilename)
            self.loadLigand(ligFilename) # reads ligand molecule from file and assigns Molecule instance to self.ligand 

        # check if TPoints are provided and if so load them so that
        pocketMode = kw.get('pocketMode', None)
        boxMode =  kw.get('boxMode', None)
        if pocketMode is None:
            if boxMode is None:
                if ligFilename is None:
                    boxMode = ['receptor']
                    pocketMode = ['best']
                else:
                    boxMode = ['ligand']
                    pocketMode = ['all']
            else:
                if boxMode[0]=='receptor':
                    pocketMode = ['best']
                else:
                    pocketMode = ['all']
                
        self.data['pocketmode'] = pocketMode
        if pocketMode[0]=='user':
            self.myprint( 'loading fill points from: %s\n'% pocketMode[1:])
            self.loadFillPoints(pocketMode[1:])
        # other pocket modes require identifying them with AutoSite
        # which neds the box to be set first
        
        # set box
        if boxMode is None:
           if ligFilename is not None:
               boxMode = ['ligand']
           else:
               boxMode = ['receptor']

        # set the box
        if boxMode[0]!='fill':
            self.setBox(boxMode, self.padding, self.spacing)
            self.myprint( 'set box using %s'% boxMode[0])
            self.myprint( '    Box center: %9.3f %9.3f %9.3f'%tuple(self.boxCenter))
            self.myprint( '    Box length: %9.3f %9.3f %9.3f'%tuple(self.boxLengths))
            spacing = kw['spacing']
            self.myprint( '    Box size  : %9d %9d %9d'%tuple(self.boxSize))
            self.myprint( '    padding   : %9.3f'%self.padding)
            self.myprint( '    spacing   : %9.3f'%self.spacing)
        else:
            self.setBox(['receptor'], self.padding, 1.0)
            
        self.data['flexResStr'] = kw['flexres']
        covalentBond = kw['covalentBond']
        self.data['covalentBond'] = covalentBond
        self.data['covalentRes'] = kw['covalentRes']

        wWeight = kw['waterWeight']
        wEntropy = kw['waterEntropy']
        self.setWmapParams(weight=wWeight, ENTROPY=wEntropy)
        #filename = kw.get('outputFile', None)
        if not filename:
            filename = os.path.basename(self.receptor.filename)
            filenameBase = os.path.splitext(filename)[0]

        if covalentBond is not None:
            s3 = kw.get('covalentBondTorsionAtom', None)
            if s3 is None:
                print 'ERROR: covalentBond specified but no covalentBondTorsionAtoms specified'
                sys.exit(1)
            # set data for 'covalentBondAtom1' , 'covalentBondAtom2', 'covalentBondTorsionAtom', 'covalentAtomsCoords', computes self.covalentBondToExclude.
            self.setCovalentDocking(s3, covalentBond, kw['covalentRes'])
            err = self.checkCovBondParams()
            if len(err):
                msg = "Unable to proceed due to the following error(s):\n"
                for e in err:
                    msg += e[1] + "\n"
                print msg
                sys.exit(1)
            top = 0
        else: # identify pockets
            cutoff = kw['pocketCutoff']
            t0 = time()
            if pocketMode[0] == "user":
                # load TPoints from file
                #tpoints = numpy.load(pocketMode[1])
                pockets = [numpy.load(pocketMode[1])]
                self.tpoints = self.pockets[0]
                self.myprint( '    loading %d fill points from %s'%(pocketMode[1], len(tpoints)))
            else:  # run AutoSite
               err = self.checkFlexResidues( flexResStr=kw['flexres'])
               if len(err):
                   msg = "Unable to proceed due to the following error(s):\n"
                   for e in err:
                       if isinstance (e[1], (tuple, list)):
                           msg += e[1][0] + "\n"
                       else:
                           msg += e[1] + "\n"
                   print msg
                   sys.exit(1)

               gc, process = self.runAutoSite(flexResStr=kw['flexres'],
                        smooth=kw['smooth'], background=False, verbose=True)

               if process!=0:
                    #self.myprint("ERROR: running autogrid failed")
                    raise RuntimeError("ERROR: running autogrid failed")

               #if kw['origAutoSite']:
               asversion = kw.get('autoSiteVersion', 1.0)
               if asversion == 1.1:
                   self.autoSite2 = True
                   self.setAutoSiteVersion("1.1", ligandSize=kw['ligandSize'], pepScore=kw['pepScore']) 
                   #run AutoSite2 for pocket detection
                   self.clustersorted, clPropsorted, dcl = self.afterAutoSite2(gc, ligandSize=kw['ligandSize'], pepScore=kw['pepScore'], verbose=True)
                   
               else: # asversion == 1.0:
                   self.autoSite2 = False
                   self.setAutoSiteVersion("1.0") 
                   #run original AutoSite for pockets
                   self.clustersorted, clPropsorted, dcl = self.afterAutoSite(gc,  verbose=True)
                   
               pockets = self.getFillPoints(pocketMode[0], cutoff, self.clustersorted, verbose=True)
               nfillPoints = sum(map(len, pockets))
               #import pdb;pdb.set_trace()
               ## if kw['origAutoSite'] is False:
               ##     #Shrink pocket to 1/5 of original size to get better tranlation points 
               ##     from AutoSite.shrink import shrinkPocket
               ##     shrinkedPoints = shrinkPocket(pockets[0],finalsize=nfillPoints/5)
               ##     pockets=[]
               ##     pockets.append([numpy.asarray(x) for x in shrinkedPoints])
               ##     nfillPoints = sum(map(len, pockets))
               self.myprint('done. got %s fill Points, in %.2f (sec)'%(nfillPoints, time()-t0))
            if pocketMode[0]=='forEach':
                top = len(self.clustersorted)
            elif pocketMode[0]=='forTop':
                top = cutoff
            else:
                top = 0
                
        # find out which types are needed
        mapTypes = self.setMapTypes(kw['mapTypes'])
        self.myprint( '\nsetting map types using: %s to %s'%(mapTypes, self.atypes))
        t0 = time()
        if top == 0: # we will have a single output file
            size = self.boxSize
            if size is None:
                print "ERROR: docking box is not specified \n"
                sys.exit(1)
            self.myprint('\ncomputing maps for center=(%.3f %.3f %.3f) size=(%.3f %.3f %.3f) dims=(%d %d %d) ...'%(
                self.boxCenter[0], self.boxCenter[1], self.boxCenter[2],
                self.boxLengths[0], self.boxLengths[1], self.boxLengths[2],
                size[0], size[1], size[2]))
            if covalentBond is None:
                self.fillPoints, outside = self.pointsInBox(pockets[0])
                if len(self.fillPoints)==0:
                    self.myprint( 'ERROR: no fill points found inside the docking box, giving up\n')
                    raise RuntimeError("no fill points found inside the docking box")
                else:
                    self.myprint( '    %d points inside the box\n'%len(self.fillPoints))

                if boxMode[0]=='fill':
                    self.myprint( 'set box using fill with %d points'%len(self.fillPoints))
                    self.setBoxForCoords(self.fillPoints, kw['padding'], kw['spacing'])
            
            gc, status = self.computeGrids(filename, kw['flexres'],
                                           kw['spacing'], indent="    ",
                                           addGradients=self.receptorGradient)
            if status !=0:
                #self.myprint('ERROR: AutoGrid failed to run in %s'%gc.folder)
                raise RuntimeError('ERROR: AutoGrid failed to run in %s'%gc.folder)
        else: # we will have several output files
            if pocketMode[0] =='forEach':
                self.myprint( '\ncreating maps for pockets with more than %d points ...'%cutoff)
            elif pocketMode[0] =='forTop':
                self.myprint( '\ncreating maps for %d top rancking pockets ...'%top)


            for n, fp in enumerate(pockets):
                self.fillPoints, outside = self.pointsInBox(fp)
                #if pocketMode[0] =='forEach' and len(self.fillPoints)< cutoff:
                #    continue
                #elif pocketMode[0] =='forTop' and n >= cutoff:
                #    continue
                if boxMode[0]=='fill':
                    self.setBoxForCoords(self.fillPoints, kw['padding'], kw['spacing'])
                filename = filenameBase+'_pocket%03d'%n
                # the box is set and we have valid TPoints, we can compute the grids
                size = self.boxSize#[int(ceil(x/kw['spacing'])) for x in self.boxLengths]
                self.myprint('    computing maps for center=(%.3f %.3f %.3f) size=(%.3f %.3f %.3f) dims=(%d %d %d) ...'%(
                    self.boxCenter[0], self.boxCenter[1], self.boxCenter[2],
                    self.boxLengths[0], self.boxLengths[1], self.boxLengths[2],
                    size[0], size[1], size[2]))
                
                self.myprint('    %d points inside the box\n'%len(self.fillPoints))  
                gc, status = self.computeGrids(
                    filename, kw['flexres'], kw['spacing'],
                    indent="    ")

            if status !=0:
                #self.myprint('ERROR: AutoGrid failed to run in %s'%gc.folder)
                raise RuntimeError('ERROR: AutoGrid failed to run in %s'%gc.folder)
        self.myprint('    done. %.2f (sec)\n'%(time()-t0))
        if self.summaryFP:
            self.summaryFP.close()

    def checkFlexResidues(self, flexResStr):
        # Check if the moving atoms of flex residues are inside the box
        flexresList =  flexResStr2flexRes(flexResStr)
        outsideRes = {}
        err = []
        #import pdb; pdb.set_trace()
        for frchain in flexresList:
            ch = frchain[0]
            for fr in frchain[1]:
                _flres = [[ch, [fr]]]
                flexresAtoms = splitFlexRes(self.receptor, _flres, exclude='CA N C O')[1]
                inside, outside = self.pointsInBox(flexresAtoms.getCoords())
                if len(outside):
                    if not outsideRes.has_key(ch):
                        outsideRes[ch] = ""
                    else:
                        outsideRes[ch]+=", "
                    outsideRes[ch]+="%s%d"%(fr[0], fr[1])
        if len(outsideRes):
            outsideResStr = ""
            for k, v in outsideRes.items():
                outsideResStr+="%s: %s " %(k, v)
            #print "FLEXRES OUTSIDE BOX:", outsideResStr
            err.append((103, "Flexible residue(s) %s outside the box." % outsideResStr))
        #Check that CHI angles atoms can be found by name.
        from ADFR.AARotamers import RotamerLib
        rotlib = RotamerLib()
        for chid, residues in flexresList:
            for resname, resnum in residues:
                try:
                    angleDef, angleList, angleDev = rotlib.get(resname)
                except KeyError as e:
                    errmsg = e.message
                    if not len(errmsg):errmsg = errorCodes[110]
                    err.append((110, [errmsg, (chid, resname, resnum)]))
                    return err
                for i, adef in enumerate(angleDef):
                    if chid == ' ':
                        tatoms = self.receptor.select('chid "%s" resname %s resnum %s name %s %s %s %s'%(
                            (chid, resname, resnum)+tuple(adef[0])))
                    else:
                        tatoms = self.receptor.select('chid %s resname %s resnum %s name %s %s %s %s'%(
                            (chid, resname, resnum)+tuple(adef[0])))
                    if len(tatoms)<4:
                        err.append((108, "Chi angle defining atoms (%s: %s%s) not found in receptor" %(chid, resname, resnum)))
        return err

    def checkCovBondParams(self):
        err = []
        if not self.data.has_key('covalentBond'):
            err.append((109, errorCodes[109]))
        elif not len(self.data['covalentBond']):
            err.append((109, errorCodes[109]))
        else:
            # check that the covalent bond atoms are in Box
            coords = numpy.array(self.data['covalentAtomsCoords']).reshape(3,3)
            inside, outside = self.pointsInBox(coords)
            if len(outside):
                #print "ERROR: %s" % errorCodes[104]
                err.append((104, errorCodes[104]))
        return err

    def checkComputeGrids(self):
        # This method is called by AGFR GUI when any of the parameters exposed
        # by the gui change.

        # The __call__() method (command line agfr) utilizes some of the components
        # (methods) from this function.

        # Check if all conditions are met to compute grids.
        err = []
        # Box
        if self.boxSize is None:
            #print "ERROR: %s" % errorCodes[101]
            return [(101, errorCodes[101])]
        # Receptor
        if self.receptor is None:
            #print "ERROR: %s" % errorCodes[100]
            err.append((100, errorCodes[100]))
        else:
            # Check that at least some atoms of the receptor are in the box:
            coords = self.receptor._ag.getCoords()
            inside, outside = self.pointsInBox(coords)
            if not len(inside):
                #print "ERROR: %s" % errorCodes[102]
                err.append((102, errorCodes[102]))

        # Check if the moving atoms of flex residues are inside the box
        if self.data.has_key('flexResStr') and self.data['flexResStr']:
            _err = self.checkFlexResidues(self.data['flexResStr'])
            if len(_err):
                err.extend(_err)
        # in covalentBond:
        if self.covalentBond == True:
            _err = self.checkCovBondParams()
            if len(_err):
                err.extend(_err)
        else: # not covalent
            # Check that translational point are in the box
            if not len(self.fillPoints):
                #print "ERROR: %s" % errorCodes[105]
                err.append((105, errorCodes[105]))
            else:
                inside, outside = self.pointsInBox(self.fillPoints)
                if not len(inside):
                    #print """ERROR: %s""" % errorCodes[106]
                    err.append((106, errorCodes[106]))
        if not len(self.atypes):
            #print """ERROR: %s""" % errorCodes[107]
            err.append((107, errorCodes[107]))
        if len(err):
            return err
        else:
            return [(0, "Ready to compute maps")]

