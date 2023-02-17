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

############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner and TSRI 2015
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/ADFR/utils/MakeGrids.py,v 1.34.2.1 2017/09/12 00:18:05 sanner Exp $
#
# $Id: MakeGrids.py,v 1.34.2.1 2017/09/12 00:18:05 sanner Exp $
#
import platform, os, subprocess, tempfile, numpy, shutil
from ADFR.utils.MakeGpf import ADGPF
from ADFRcc.adfr import GridMap
from MolKit2.selection import Selection
from time import time
from ADFR.utils.maps import flexResStr2flexRes

def splitFlexRes(receptor, flexibleResidues, exclude='N C O'):
    """
    return a list of atoms for the receptor without the flexible side chains
    and a list of atoms for the flexible side chains including Cb and Ca

    flexres is a list of residue names and numbers eg.
        [ ['A', ('GLU', '167'), ('TYR', '628'), ('ARG', '87')],
        ['B', ('GLU', '167')] ]
    """
    selStr = ''
    ct = 0
    for chid, flexRes in flexibleResidues:
        if ct > 0:
            selStr += " or "
        if chid==' ':
            selStr += '(chid "%s" and ('%chid
        else:
            selStr += '(chid %s and ('%chid
        for resname, resnum in flexRes:
            selStr += 'resname %s resnum %s or '%(resname, resnum)

        selStr = selStr[:-3]+ '))'
        if exclude:
            selStr = selStr + ' and not name %s'%exclude # this keeps HN attached to N
        ct += 1    
    sideChainAtoms = receptor.select(selStr)
    # remove H atoms attached to N
    toRemove = []
    for atom in sideChainAtoms:
        if atom.getElement()=='H':
            for natom in atom.iterBonded():
                if natom.isbackbone:
                    toRemove.append(atom.getIndex())

    sideChainAtoms = sideChainAtoms - Selection(receptor._ag, toRemove, '')
    receptorAtoms = receptor.select() - sideChainAtoms
    return receptorAtoms, sideChainAtoms

def findBinary(filename):
    # find the binary in the InstallDir/bin folder (in Linux type machines)
    binary = None
    if (os.name != 'nt'):
        root = os.getenv('ADS_ROOT')
        if root:
            fullName = os.path.abspath(os.path.join(root, "bin", filename))
            if os.path.exists(fullName) and os.access(fullName, os.X_OK):
                binary = fullName
    else: #on Windows the  binary should be in the PATH env variable
        # it is set in the application script (agfr.bat or agfrgui.bat)
        filename = filename+".exe" 
        for path in os.environ['PATH'].split(os.pathsep):
            ff = os.path.join(path, filename)
            if os.path.isfile(ff)and os.access(ff, os.X_OK):
                binary = ff
            break
    return binary

class CalculateAD4Grids:

    def __init__(self, receptor, center, size, atypes, spacing=0.375,
                 smooth=0.5, dielectric=-0.1456, flexibleResidues=[],
                 folder='.', atypesOnly=False, fp=False,
                 covalentBondToExclude=[], outlev=1):
        # covalentBond has to be a list of 2 integers which are r2 r3 indices of receptor atoms
        # covalentRes is a selection string
        self.logExt = '.glg'
        self.paramFile = None
        self.binary = None
        self.setWorkingFolder(folder)
        self.atypes = atypes
        self.mapFiles = {}
        self.outlev = outlev
        
        # set the default binary
        #from mglutil.util.packageFilePath import getBinary
        #adBin = getBinary('autogrid4', 'binaries')
        adBin = findBinary('autogrid4')
        assert(adBin is not None)
        self.setBinary(adBin)
        import ADFR
        # MS use original AutoGrid NOTE: using ADFRparameters.dat changes the min of the OA map :(
        self.paramFile = os.path.join(ADFR.__path__[0], 'Data','AD4.1_bound.dat')

        #if fp == True:
        #    self.paramFile = os.path.join(ADFR.__path__[0], 'Data','AD4.1_bound.dat')
        #else:
        #    self.paramFile = os.path.join(ADFR.__path__[0], 'Data', 'ADFRparameters.dat')

        #print atypesOnly, atypes
        #self.setBinary('/mgl/ms1/people/sanner/python/AD5_adfr/autodockSrc/autogrid/x86_64Linux3/autogrid4NOBHTree')
        #self.setBinary('/mgl/ms1/people/sanner/python/AD5_adfr/autodockSrc/autogrid/x86_64Linux3/autogrid4BHTree')
        self.system_info = platform.uname()
        self.platform = self.system_info[0]
        #if self.platform == 'Linux' or self.platform=='Darwin':
        #   self.platform = 'posix'

        if self.platform == "Windows":
            import ctypes # used by CheckDiskSpace on Windows

        # create receptor with missing side chains if side chains are flexible
        #print 'FR', flexibleResidues
        #import pdb; pdb.set_trace()
        try:
            shutil.copy(receptor.filename, folder)
        except:
            pass
        shutil.copy(self.paramFile, folder)
        if len(flexibleResidues)>0 or len(covalentBondToExclude)>0:
            if len(flexibleResidues)>0:
                receptorAtoms, sideChainAtoms = splitFlexRes(receptor, flexibleResidues)
            else:
                sideChainAtoms = receptor.emptySelection()
            if len(covalentBondToExclude):
                toRemoveAtoms = covalentBondToExclude
            else:
                toRemoveAtoms = receptor.emptySelection()
            #print sideChainAtoms.getNames()
            self.atomsNotInGrid = sideChainAtoms+toRemoveAtoms
            self.flexRecAtoms = sideChainAtoms
            self.covalentLigAtoms = toRemoveAtoms
            self.receptorFilename = self.writeGridReceptorPDBQT(
                receptor, 'rigidReceptor.pdbqt', self.atomsNotInGrid)
            createdReceptor = True
            # update atypes to include types from covalentLigAtoms and receptorFilename
            atypes = list(self.atypes)
            if len(sideChainAtoms):
                atypes += sideChainAtoms.getData('AD_element').tolist()
            if len(toRemoveAtoms):
                atypes += toRemoveAtoms.getData('AD_element').tolist()
            self.atypes = numpy.unique(atypes)
        else:
            self.atomsNotInGrid = receptor.emptySelection()
            self.flexRecAtoms= receptor.emptySelection()
            self.covalentLigAtoms= receptor.emptySelection()
            shutil.copy(receptor.filename, os.path.join(folder, 'rigidReceptor.pdbqt'))
            self.receptorFilename = os.path.join(folder, 'rigidReceptor.pdbqt')
            createdReceptor = False
        
        # create parameter file
        obj=ADGPF(atypesOnly=atypesOnly,
                  paramFile=os.path.basename(self.paramFile),
                  outlev=self.outlev)
        obj.setNpts(size)
        obj.setSpacing(spacing)
        obj.setReceptor_atypes(numpy.unique(receptor._ag.getData("AD_element")))
        obj.setLigand_atypes(self.atypes)
        obj.setReceptor(self.receptorFilename)
        obj.setSmooth(smooth)
        obj.setDielectric(dielectric)
        obj.setGridcenter(center)
        assert obj.isDataValid()[0] == True

        f = open(os.path.join(self.folder, 'rigidReceptor.gpf'), 'w')
        [f.write(line+'\n') for line in obj.getGPFlines()]
        f.close()
        #print "================================================================="
        #for l in obj.getGPFlines():
        #    print l
        self.setParmFile('rigidReceptor.gpf')
        
    #def cleanup(self):
    #    os.remove(os.path.join(self.folder, self.temp_name+'.gpf'))
    #    os.remove(os.path.join(self.folder, self.temp_name+'.glg'))


        
    def writeGridReceptorPDBQT(self, mol, filename, excludedAtoms):
        toExclude = {}
        for a in excludedAtoms:
            toExclude['%.3f,%.3f,%.3f'%(tuple(a.getCoords()))] = True
        f = open(mol.filename)
        lines = f.readlines()
        f.close()
        #filename = os.path.join(self.folder, os.path.split(mol.filename)[1][:-6]+'_NoSC.pdbqt')
        filename = os.path.join(self.folder, filename)
        f = open(filename, 'w')
        #import pdb; pdb.set_trace()
        nbSkip = 0
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                #w = line.split() # cannot use split because when chain id is space the
                # columns are off
                key = '%.3f,%.3f,%.3f'%(float(line[30:38]),float(line[38:46]),float(line[46:54]))
                if not toExclude.has_key(key):
                    f.write(line)
                else:
                    nbSkip += 1
            else:
                f.write(line)
        f.close()
        #print 'removed %d side chain atoms for grid calculation'%nbSkip
        return os.path.join(self.folder, 'rigidReceptor.pdbqt')#filename
    
    def setWorkingFolder(self, path):
        assert os.path.exists(path), 'The directory specified doesn\'t exist or is not accessible:\n\n%s' % path
        #FIXME make sure we have write access
        self.folder = path

    def setBinary(self, binary):
        assert os.path.exists(binary), 'The program specified doesn\'t exist or is not accessible:\n\n%s' % binary
        # FIXME check if executable
        self.binary = binary

    def setParmFile(self, paramFile):
        assert os.path.exists(os.path.join(self.folder, paramFile)), 'The program specified doesn\'t exist or is not accessible:\n\n%s' % paramFile
        self.paramFile = paramFile

    def run(self, background=False):
        assert self.paramFile is not None
        assert self.binary is not None
        assert self.folder is not None

        if self.platform == 'Windows':
            shell=False
        else:
            shell=True

        log_ext = self.logExt
        #queuer = self.adg_queue_parser
        log = os.path.splitext(self.paramFile)[0]+log_ext
        command = '"%s" -p %s -l %s'  % ( self.binary, self.paramFile, log)
        self._command = command
        if self.platform == "Darwin":
            # this is to set DYLIB_LIBRARY_PATH on Darwin so the autogrid4 binary
            # finds lbgomp.1.dylib 
            binFolder = os.path.dirname(self.binary)
            command = 'source "%s/adsenv.sh"; ' % (binFolder,) + command

        #print "COMMAND", command
        cwd = os.getcwd()
        self.logFile = log
        if os.path.exists(self.logFile):
            os.remove(self.logFile)

        self.mapFiles = {}
        for atype in self.atypes:
            self.mapFiles[atype] = '%s.%s.map'%(os.path.splitext(self.receptorFilename)[0], atype)
        self.folder = os.path.abspath(self.folder)
        os.chdir(self.folder)
        if background:
            try:
                p = subprocess.Popen(command,
                                     stdout=subprocess.PIPE , 
                                     stderr=subprocess.PIPE, 
                                     bufsize = 1, cwd=self.folder, shell=shell)
                self.process = p
            finally:
                os.chdir(cwd)
                
            return 0, ""
        else:
            try:
                #print 'CURRENT DIR', os.getcwd()
                #print 'COMMAND', command
                #status = os.system(command)
                pipes = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                std_out, std_err = pipes.communicate()
                #print "OUT:", std_out
                status = pipes.returncode
                #print 'STATUS ============== ', status
                #print "ERR:", std_err, "enderr", len(std_err)
                err=""
                if len(std_err):
                    err = std_err
                    ch1 = std_err.find("ERROR")
                    if ch1 >= 0:
                        err = std_err[ch1+5:].split("\n")[0]
            finally:
                self.process = None
                os.chdir(cwd)
            return status, err

    def getTPoints(self, carbon_cutoff=-0.30):
        cmapFile = self.mapFiles['C']
        cmap = GridMap()
        cmap.loadFromMapFile('C', '', cmapFile)
        cdata = cmap.getGridDataPy()#[::3, ::3, ::3]
        ox, oy, oz = cmap.getOriginPy()
        sx = sy = sz = cmap.getDistBetweenGridPoints()
        goodC = numpy.select( [cdata<carbon_cutoff, cdata>=carbon_cutoff], [cdata,[999]])
        i,j,k = numpy.where(goodC!=999)
        coords = numpy.zeros((len(i),3),'f')
        coords[:,0] = ox + i*sx
        coords[:,1] = oy + j*sy
        coords[:,2] = oz + k*sz
        indices = numpy.zeros((len(i),3),'i')
        indices[:,0] =  i
        indices[:,1] =  j
        indices[:,2] =  k
        self._indices = indices   # list of (i,j,k) in grid
        self._coords = coords     # list of (x,y,z) of selected grid points
        

    def addFlexRecHeader(self, line):
        #mapFilename = self.mapFiles.values()[0]
        #mapFilename = mapFilename[:-5]+'e'+mapFilename[-4:]
        for mapFilename in self.mapFiles.values():
            f = open(mapFilename)
            lines = f.readlines()
            f.close()
            f = open(mapFilename, 'w')
            f.write("%s\n"%line)
            [f.write(l) for l in lines]
            f.close()
        
if __name__=='__main__':
    from MolKit2 import Read
    rec = Read('CDK2/receptors/4EK3_rec.pdbqt')  
    gc = CalculateAD4Grids(rec, (25.834, 27.584, 27.528), (70,70,70), ['C'],
                           flexibleResidues=[
                               ('A', [('ILE','10'), ('VAL', '18'), ('LYS', '33'), ('VAL', '64'),
                                      ('PHE', '80'), ('PHE', '82'), ('GLN', '85'), ('ASP', '86'),
                                      ('LYS', '89'), ('ASN', '132'), ('LEU', '134'),
                                      ('ASP', '145')]
                                )])
    
    #rec = Read('CDK2/receptors/2CCH_rec.pdbqt')
    #gc = CalculateAD4Grids(rec, (25.834, 27.584, 27.528), (70,70,70), ['A', 'C', 'HD', 'N', 'NA', 'OA', 'P'],
    #                       flexRes=[('ILE','10'), ('VAL', '18'), ('LYS', '33'), ('VAL', '64'),
    #                                ('PHE', '80'), ('PHE', '82'), ('GLN', '85'), ('ASP', '86'),
    #                                ('LYS', '89'), ('ASN', '132'), ('LEU', '134'), ('ASP', '145')])
    #gc = CalculateAD4Grids(rec, (25.834, 27.584, 27.528), (70,70,70), ['C'],
    #                       flexRes=[('ILE','10'), ('VAL', '18'), ('LYS', '33'), ('VAL', '64'),
    #                                ('PHE', '80'), ('PHE', '82'), ('GLN', '85'), ('ASP', '86'),
    #                                ('LYS', '89'), ('ASN', '132'), ('LEU', '134'), ('ASP', '145')])
    gc.getTPoints()
    print len(gc._coords)
    
