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
# $Header: /mnt/raid/services/cvs/ADFR/utils/maps.py,v 1.22.2.1 2017/08/01 00:35:32 annao Exp $
#
# $Id: maps.py,v 1.22.2.1 2017/08/01 00:35:32 annao Exp $
#
import os, sys, tempfile, shutil, pickle
from zipfile import ZipFile
from ADFRcc.adfr import GridMap

class MapsFile:
    """
    proxy for a target file (.trg) produced by AGFR
    """

    def __del__(self):
        # clean up
        if self.isUnzipped():
            shutil.rmtree(self.unzipFolder)

    def printInfo(self, indent=''):
        d = self.getData()
        print indent+'docking target file'
        print indent+'  date       : %s'%d['date']
        print indent+'  node       : %s'%d['node']
        print indent+'  AGFR       : v%s'%d['AGFRVERSION']
        if d.has_key('AutoSiteVersion'):
            print indent+'  AutoSite   : v%s'%d['AutoSiteVersion']
        print
        print indent+'  receptor   : %s'%d['inputReceptor']
        if len(d['flexRecFile']):
            print indent+'     FlexRec : %s in %s'%(d['flexResStr'], d['flexRecFile'])
        else:
            print indent+'    FlexRec  : None'
        if d.has_key('covalentBond') and d['covalentBond'] is not None:
            print indent+'    covBond  : [%s] %s (%d) -- %s (%d)'%(
                d['covalentBondTorsionAtom'],
                d['covalentBondAtom1'], d['covalentBond'][0], 
                d['covalentBondAtom2'], d['covalentBond'][1])
            print indent+'    coords   : [%.3f, %.3f, %.3f] [%.3f, %.3f, %.3f] [%.3f, %.3f, %.3f]'%tuple(d['covalentAtomsCoords'])

            print indent+'    covRes   : %s'%d['covalentRes']
            print indent+'    covFile  : %s'%d['covalentLigandFile']
            print indent+'    ignAtms  : %s'%d['covalentLigandAtomIndices']
        else:
            print indent+'    covBond  : No'
        #if d.has_key('inputLigand'):
        #    print
        #    print '  ligand     : %s'%d['inputLigand']
        print
        print indent+'  box        :  mode %s, padding %.2f'%(d['boxMode'],d['boxPadding'])
        print indent+'    center   : %.3f %.3f %.3f'%tuple(d['boxCenter'])
        print indent+'    length   : %.3f %.3f %.3f'%tuple(d['boxLengths'])
        print indent+'    size     : %.4d %.4d %.4d'%tuple(d['boxSize'])
        print indent+'    spacing  : %.3f'%d['spacing']
        print
        print indent+'  maps       : '
        print indent+'    types    :',
        n = 0
        for at in d['mapTypes']:
            print '%2s '% at,
            if n == 20:
                print '\n'+indent+'              ',
            n +=1
        print
        print indent+'    W map    : weight %.2f entropy %.2f'%(d['wMapWeight'],d['wMapEntropy'])
        if d['mapGradients']:
            print indent+'    gradients: Yes, ',
            if d['gradCutOff']==-1:
                print "kept largest negative cluster"
            else:
                print "kept clusters pockets larger than %d"%d['gradCutOff']
        else:
            print indent+'    gradients: No'
            
        print
        print indent+'  pocketMode : %s'%d['boxMode']
        print indent+'    #fillpts : %d points'%d['nbFillPoints']
        if d.has_key('fillPointsFile'):
            print indent+'    file     : %s'%d['fillPointsFile']
        
    def __init__(self, filename):
        self.filename = filename
        self.unzipFolder = None # temporary folder created to unzip the file 
        self.folder = None # folder inside self.unzipFolder
        self._data = None # will hold the data.pkl content
        self._maps = {} # will contain actual maps as they are read
        
    def isUnzipped(self):
        return self.unzipFolder is not None

    def unzipMaps(self):
        # create folder in temporary space
        tmpFolder = self.unzipFolder = tempfile.mkdtemp()
        zf = ZipFile(self.filename)
        zf.extractall(tmpFolder)
        self.folder = os.path.join(tmpFolder, os.path.splitext(os.path.basename(self.filename))[0])
        self.getData()

    def loadAllMaps(self):
        for typ in self.getMapTypes():
            if self._maps.get(typ, None) is None:
                self.getMap(typ)
        
    def getMap(self, mtype):
        if not self.isUnzipped():
            self.unzipMaps()
        assert mtype in self._data['mapTypes']
        _map = GridMap()
        _map.loadFromMapFile(mtype, self.folder, 'rigidReceptor.%s.map'%mtype)
        self._maps[mtype] = _map
        return _map
    
    def getData(self):
        # get raw data dictionary
        if self._data:
            return self._data
        elif self.isUnzipped():
            with open(os.path.join(self.folder, 'data.pkl')) as f:
                self._data = pickle.load(f)
            return self._data
        else:
            zf = ZipFile(self.filename, 'r')
            dataFile = None
            for name in zf.namelist():
                if os.path.basename(name)=='data.pkl':
                    dataFile = name
                    break
            if dataFile:
                # the replace is needed else on windows we get
                # insecure string pickle. If the file is eactually expanded
                # there \r values are not present
                if os.name=='nt':
                    return pickle.loads(zf.read(dataFile).replace('\r', ''))
                else:
                    return pickle.loads(zf.read(dataFile))
            else:
                return None

    def getMapsFolder(self):
        return self.folder

    def getReceptorFilename(self):
        data = self.getData()
        return data['inputReceptor']

    def getLigandFilename(self):
        data = self.getData()
        return data.get('inputLigand', None)

    def getFlexResStr(self):
        data = self.getData()
        return data.get('flexResStr', None)

    def getFlexRecFile(self):
        data = self.getData()
        return data['flexRecFile']

    def getDate(self):
        data = self.getData()
        return data['date']
    
    def getPlatform(self):
        data = self.getData()
        return data['platform']
    
    def getNode(self):
        data = self.getData()
        return data['node']

    def getCovalentBond(self):
        data = self.getData()
        return data.get('covalentBond', None)
    
    def getCovalentRes(self):
        data = self.getData()
        return data.get('covalentRes', None)
    
    def getCovalentLigandAtomIndices(self):
        data = self.getData()
        return data.get('covalentLigandAtomIndices', None)
    
    def getCovalentBondAtom1(self):
        data = self.getData()
        return data.get('covalentBondAtom1', None)
    
    def getCovalentBondAtom2(self):
        data = self.getData()
        return data.get('covalentBondAtom2', None)

    def getCovalentBondTorsionAtom(self):
        data = self.getData()
        return data.get('covalentBondTorsionAtom', None)
    
    def getCovalentLigandFile(self):
        data = self.getData()
        return data.get('covalentLigandFile', None)
        
    def getMapTypes(self):
        data = self.getData()
        return data['mapTypes']

    def getFillPointsFile(self):
        data = self.getData()
        return data['fillPointsFile']

    def getNumFillPoints(self):
        data = self.getData()
        return data['nbFillPoints']

    def getPocketMode(self):
        data = self.getData()
        return data['pocketmode']

    def getBoxMode(self):
        data = self.getData()
        return data['boxMode']

    def getBoxPadding(self):
        data = self.getData()
        return data['boxPadding']

    def getBoxCenter(self):
        data = self.getData()
        return data['boxCenter']
    
    def getBoxLength(self):
        data = self.getData()
        return data['boxLengths']
    
    def getBoxSize(self):
        data = self.getData()
        return data['boxSize']
    
    def getBoxSpacing(self):
        data = self.getData()
        return data['spacing']
    
    def getAutoSiteVerison():
        data = self.getData()
        return data['AutoSiteVersion']
    
    def getAutoSiteLigandSize():
        #this attribute is for AutoSite2 (version 1.1)
        data = self.getData()
        return data.get('ligandSize', None)
    
    def getAutoSitePepScore():
        #this attribute is for AutoSite2 (version 1.1)
        data = self.getData()
        return data.get('pepScore', None)

    def getAgfrVersion():
        data = self.getData()
        return data['AGFRVERSION']

    def getWaterMapEntropy():
        data = self.getData()
        return data['wMapEntropy']

    def getWaterMapWeight():
        data = self.getData()
        return data['wMapWeight']

    def hasReceptorGradient():
        # returns True if receptor gradient was calculated
        # for the maps, False - otherwise
        data = self.getData()
        return data['']

def flexResStr2flexRes(flexResStr):
    if flexResStr:
        flexRes = []
        for expr in flexResStr.split(';'):
            chid, residues = expr.split(':')
            fr = []
            for res in residues.replace(',', ' ').split():
                resname = res[:3]
                try:
                    resnum = int(res[3:])
                except ValueError:
                    print 'ERROR: invalid syntax\n expect a residue number but got "%s\n"'%res[3:]
                    sys.exit(1)
                fr.append( (resname.upper(), resnum) )           
            flexRes.append( [chid, fr] )
    else:
        flexRes = []
    return flexRes

def getTargetData(filename):
    zf = ZipFile(filename, 'r')
    dataFile = None
    for name in zf.namelist():
        if os.path.basename(name)=='data.pkl':
            dataFile = name
            break
    if dataFile:
        return pickle.loads(zf.read(dataFile))
    else:
        return None

## def printTargetData(d):
##     print 'docking target file'
##     print '  date       : %s'%d['date']
##     print '  node       : %s'%d['node']
##     print '  platform   : %s'%d['platform']
##     print
##     print '  receptor   : %s'%d['inputReceptor']
##     if len(d['flexRecFile']):
##         print '     FlexRec : %s in %s'%(d['flexResStr'], d['flexRecFile'])
##     else:
##         print '    FlexRec  : None'
##     if d['covalentBond'] is not None:
##         print '    covBond  : [%s] %s (%d) -- %s (%d)'%(
##             d['covalentBondTorsionAtom'],
##             d['covalentBondAtom1'], d['covalentBond'][0], 
##             d['covalentBondAtom2'], d['covalentBond'][1])
##         print '    coords   : [%.3f, %.3f, %.3f] [%.3f, %.3f, %.3f] [%.3f, %.3f, %.3f]'%tuple(d['covalentAtomsCoords'])
            
##         print '    covRes   : %s'%d['covalentRes']
##         print '    covFile  : %s'%d['covalentLigandFile']
##     else:
##         print '    covBond  : No'
##     #if d.has_key('inputLigand'):
##     #    print
##     #    print '  ligand     : %s'%d['inputLigand']
##     print
##     print '  box        :  '
##     print '    center   : %.3f %.3f %.3f'%tuple(d['boxCenter'])
##     print '    length   : %.3f %.3f %.3f'%tuple(d['boxLengths'])
##     print '    size     : %.4d %.4d %.4d'%tuple(d['boxSize'])
##     print '    spacing  : %.3f'%d['spacing']
##     print
##     print '  maps       : ',
##     n = 0
##     for at in d['mapTypes']:
##         print '%2s '% at,
##         if n == 30:
##             print '\n                ',
##     print
##     print '  pocketMode : %s'%d['boxMode']
##     print '    #fillpts : %d points'%d['nbFillPoints']
##     print '    file     : %s'%d['fillPointsFile']

## def unzipMaps(filename):
##     """
##         Unzip a map file in a temporary location and return the folder name,
##         receptor name flexible residues descriptor and a status message
##     """
##     #from mglutil.util.packageFilePath import getResourceFolderWithVersion
##     #tmpFolder = os.path.join(getResourceFolderWithVersion(), 'tmp')
##     #if not os.path.exists(tmpFolder):
##     #    os.mkdir(tmpFolder)
##     #elif not os.path.isdir(tmpFolder):
##     #    return None, None, None, "ERROR: %s is not a folder, please remove this file"%tmpFolder
##     import tempfile
##     # create folder in temporary space
##     tmpFolder = tempfile.mktemp()

##     zf = ZipFile(filename)
##     zf.extractall(tmpFolder)
##     mapsFolder = folder = os.path.join(tmpFolder, os.path.splitext(os.path.basename(filename))[0])
##     tpointsFilename = os.path.join(folder, "translationPoints.npy")
##     mapFilesRoot = 'receptor'
##     # find receptor file name i.e PDBQT file which is NOT rigidReceptor.pdbqt
##     for filename in zf.namelist():
##         if filename.endswith('.pdbqt'):
##             filename = os.path.split(filename)[1]
##             recname = os.path.splitext(filename)[0]
##             if recname != 'rigidReceptor':
##                 break
##     receptorFilename = os.path.join(folder, '%s.pdbqt'%recname) 

##     flexResStr = None
##     from ADFRcc.adfr import GridMap
##     #mapTypes = []
##     for mapFileName in zf.namelist():
##         w = mapFileName.split('.')
##         if w[-1]=='map':
##             #mapTypes.append(w[-2])
##             if not (w[-2]=='e' or w[-2]=='d'): # read 'e' map to get center, box size, spacing, and flexres
##                 _map = GridMap()
##                 _f, name = os.path.split(mapFileName)
##                 _map.loadFromMapFile(w[-2], folder.encode('ascii', 'replace'),
##                                      name.encode('ascii', 'replace'))
##                 flexResStr = _map.getFlexRes()
##                 break
##     return mapsFolder, receptorFilename, flexResStr, 'OKAY', _map

if __name__=='__main__':
    mf = MapsFile('4EK3_rec_cmdline.trg')
    # without unzipping the file we can find out the follwing
    assert len(mf._maps) == 0
    assert mf._data == None
    data = mf.getData() # get data without unzipping
    assert len(data) > 0
    mf.printInfo()
    assert mf.getReceptorFilename() is not None
    assert mf.getLigandFilename() is not None
    print mf.getFlexResStr()
    print mf.getFlexRecFile()
    print mf.getDate()
    print mf.getPlatform()
    print mf.getNode()
    print mf.getCovalentBond()
    print mf.getCovalentRes()
    print mf.getCovalentBondAtom1()
    print mf.getCovalentBondAtom2()
    print mf.getCovalentBondTorsionAtom()
    print mf.getCovalentLigandFile()
    mapTypes = mf.getMapTypes()
    print mapTypes
    print mf.getFillPointsFile()
    print mf.getNumFillPoints()
    print mf.getBoxMode()
    print mf.getBoxPadding()
    print mf.getBoxCenter()
    print mf.getBoxLength()
    print mf.getBoxSize()
    print mf.getBoxSpacing()

    assert mf.isUnzipped() == False
    assert mf.getMapsFolder() == None

    # now unzip the file
    mf.unzipMaps()
    assert mf._maps == {}
    assert mf.getMapsFolder() is not None

    emap = mf.getMap('e')
    assert 'e' in mf._maps.keys()
    # need to delete before quiting else we get
    # Exception RuntimeError: 'sys.meta_path must be a list of import hooks' in <bound method MapsFile.__del_

    mf.loadAllMaps()
    del mf
