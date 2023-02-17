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
# Copyright: M. Sanner TSRI 2014
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/PmvApp/fileCmds.py,v 1.8.4.1 2017/07/13 20:55:28 annao Exp $
#
# $Id: fileCmds.py,v 1.8.4.1 2017/07/13 20:55:28 annao Exp $
#

from PmvApp.Pmv import MVCommand
from MolKit2.molecule import MoleculeSet

import sys, os

class MoleculesReader(MVCommand):
    """Command to read molecules.
    filenames --- a list of molecule files,
    addToRecent --- if set to True, adds to the list of application recent files;
    group --- can be None, a MoleculeGroup instance or a name(string) of existing group, if specified , the molecules are added to the group.
    """
    def expandArg0(self, obj):
        if isinstance(obj, list) : return obj
        else: return [obj]

    def doit(self, filename, addToRecent=True, group=None, header=True):

        mol = self.app().readMolecule(filename, addToRecent=addToRecent,
                                      group=group, header=header)
        self.app().pushUndoCmd(
            self.app().deleteMolecule,
            (MoleculeSet('molset', [mol]),), { "redraw":True}
            )
        return mol
    
    ## def undoCmdAfter(self, result, *args, **kw):
    ##     #this is called in VFCmd.afterDoit(). result is a returned value of the doit().
    ##     # (in this case - a list of molecule objects)
    ##     name = ""
    ##     for mol in result:
    ##         name += mol.name +','
    ##     if not hasattr(self.app(), "deleteMolecules"):
    ##         return
    ##     return ([(self.app().deleteMolecules, (result,), {'undoable':0, 'cleanRedo':False, "redraw":True})],
    ##             self.name+ " %s"%name[:-1])

    def checkArguments(self, filenames, addToRecent=True,
                       group=None, header=True):

        assert isinstance(filenames, (list, tuple))
        for name in filenames:
            assert isinstance(name, str), "File names have to be strings"
        assert addToRecent in [0,1,True, False], 'got %s'%addToRecent
        assert header in [0,1,True, False], 'got %s'%addToRecent
        args = (filenames,)
        if group is not None:
            from MolKit2.molecule import MoleculeGroup
            assert isinstance (group ,(MoleculeGroup, str))
        kw = {'addToRecent':addToRecent, 'group':group, 'header':header}
        return args, kw


class ReadAny(MVCommand):
    """Reads Molecule or Python script pmv session, vision network or 3D grids \n
    Package : PmvApp \n
    Module  : fileCmds \n
    Class   : ReadAny \n
    Command : readAny \n
    Synopsis:\n
        None <- readAny([files])
    """        

    def doOneFile(self, file, addToRecent=True):
        
        ext = os.path.splitext(file)[1].lower()
        if ext in ['.py','.rc']:
            if file.find('_pmvnet') > 0 or file.find('_net') > 0:
                vision = self.app().vision
                if vision.isLoader(): vision = vision.loadCommand()
                vision.show()
                vision.ed.loadNetwork(file)
            else:
                self.app().source(file, log=0)
        elif ext == ".psf":
            # Load Pmv session file
            self.app().readFullSession(file)
        elif ext in ['.pdb', '.pdbqt', '.pqr', '.cif', '.mol2']:
            self.app().readMolecule(file, addToRecent=addToRecent)
        else: # not a recognized format
            raise ValueError, "file %s is not a recognized format for %s"%(
                file, self.name)


    def doit(self, files, addToRecent=True):
        for ff in files:
            try:
                self.doOneFile(ff, addToRecent=addToRecent)
            except:
                if self.app().trapExceptions is False:
                    exc_info = sys.exc_info()
                    raise exc_info[1], None, exc_info[2]
                else:
                    msg = 'Error while reading %s' %ff
                    self.app().errorMsg(sys.exc_info(), msg, obj=ff)


    def checkArguments(self, files, addToRecent=True):
        """
        files - list of filenames
        """
        for name in files:
            assert isinstance(name, str), "names in filenames should be strings %s"%repr(name)
        kw = {}
        assert addToRecent in [True, False, 1, 0]
        kw['addToRecent'] =  addToRecent
        return (files,), kw


class ReadPmvSession(MVCommand):
    """Reads Full Pmv Session \n
    Package : PmvApp \n
    Module  : fileCmds \n
    Class   : ReadPmvSession \n
    Command : readPmvSession \n
    Synopsis:\n
        None <- readPmvSession(path)
    """        

    def doit(self, files):
        for name in files:
            try:
                self.app().readFullSession(name)
            except:
                if self.app().trapExceptions is False:
                    exc_info = sys.exc_info()
                    raise exc_info[1], None, exc_info[2]
                else:
                    msg = 'Error while reading %s'%name
                    self.app().errorMsg(sys.exc_info(), msg, obj=file)
                    

    def checkArguments(self, files):
        if not isinstance(files, (list, tuple)):
            files = [files,]
        for file in files:
            assert isinstance(file, str)
        return (files,), {}
        
from mglutil.util.packageFilePath import getCacheFolder

class Fetch(MVCommand):
    """This command reads molecule(s) from the web. /n
    Default format is mmtf. Optional format - pdb.
    """
        
    def expandArg0(self, obj):
        if isinstance(obj, list) : return obj
        else: return [obj]

    def checkArguments(self, pdbID=[], ext="mmtf"):
        assert isinstance(pdbID, (list, tuple))
        for name in pdbID:
            assert isinstance(name, str), "Pdb IDs have to be strings"
        assert ext in ("pdb", "mmtf")
        return (pdbID, ext), {}


    def doit(self, pdbID, ext):
        cacheFolder = getCacheFolder()
        mol, msg = self.app().fetchFile(pdbID, ext, pdbFolder=cacheFolder)
        if mol:
            self.app().pushUndoCmd(
            self.app().deleteMolecule, (MoleculeSet('molset', [mol,]),), { "redraw":True})
        return mol, msg

    def clearPDBCache(self):
        cacheFolder = getCacheFolder()
        dirname = cacheFolder
        if dirname and os.path.exists(dirname):
                filenames = os.listdir(dirname)
                for f in filenames:
                    os.remove(dirname+os.sep+f)

    def checkCache(self, threshold = 1024*1024):
        cacheFolder = getCacheFolder()
        size = self.getCacheSize()
        maxSize = self.app().userpref['PDB Cache Storage (MB)']['value']
        if size > maxSize:
            if not cacheFolder: return
            folder_size = 0
            for (path, dirs, files) in os.walk(cacheFolder):
              for file in files:
                filename = os.path.join(path, file)
                fileSize = os.path.getsize(filename)
                if fileSize > threshold:
                    os.remove(filename)
                else:
                    folder_size += os.path.getsize(filename)
            if (folder_size/(1024*1024.0)) > maxSize:
                self.checkCache(threshold = threshold/2.)
                
    def getCacheSize(self):
        # pick a folder you have ...
        cacheFolder = getCacheFolder()
        if not cacheFolder: return
        folder_size = 0
        for (path, dirs, files) in os.walk(cacheFolder):
          for file in files:
            filename = os.path.join(path, file)
            folder_size += os.path.getsize(filename)
        
        return (folder_size/(1024*1024.0))
                       


class PDBWriter(MVCommand):
    """
    Command to write the given molecule or the given subset of atoms
    of one molecule as PDB file.
    \nPackage : PmvApp
    \nModule  : fileCommands
    \nClass   : PDBWriter 
    \nCommand : writePDB
    \nSynopsis:\n
        None <- writePDB( nodes, filename=None, sort=True,
                          pdbRec=['ATOM', 'HETATM', 'MODRES', 'CONECT'],
                          bondOrigin=('File', 'UserDefined'), ssOrigin=None, **kw)
    \nRequired Arguments:\n    
        nodes --- TreeNodeSet holding the current selection
    \nOptional Arguments:\n
        filename --- name of the PDB file (default=None). If None is given
                  The name of the molecule plus the .pdb extension will be used\n
        sort --- Boolean flag to either sort or not the nodes before writing
                  the PDB file (default=True)\n
        pdbRec --- List of the PDB Record to save in the PDB file.
                  (default: ['ATOM', 'HETATM', 'MODRES', 'CONECT']\n
        bondOrigin --- string or list specifying the origin of the bonds to save
                    as CONECT records in the PDB file. Can be 'all' or a tuple\n

        ssOrigin --- Flag to specify the origin of the secondary structure
                    information to be saved in the HELIX, SHEET and TURN
                    record. Can either be None, File, PROSS or Stride.\n
    """

    def __init__(self):
        MVCommand.__init__(self)

        
    def doit(self, nodes, filename, sort=True, transformed=True,
             pdbRec=['ATOM', 'HETATM', 'MODRES', 'CONECT'],
             bondOrigin=('File', 'UserDefined'), ssOrigin=None):
        
        if transformed:
            oldCoords = {}
            from MolKit2.molecule import Atom
            molecules, atomSets = self.app().getNodesByMolecule(nodes, Atom)
            for mol, atoms in zip(molecules, atomSets):
                coords = atoms.coords
                # save curent coords
                oldCoords[mol] = coords
                # set transformed coords:
                atoms.coords = self.app().getTransformedCoords(mol, coords)
        from MolKit2.pdbWriter import PdbWriter
        writer = PdbWriter()
        try:
            writer.write(filename, nodes, sort=sort, records=pdbRec,
                     bondOrigin=bondOrigin, ssOrigin=ssOrigin)
        except:
            if self.app().trapExceptions is False:
                exc_info = sys.exc_info()
                raise exc_info[1], None, exc_info[2]
            else:
                msg = 'Error while writing %s'%filename
                self.app().errorMsg(sys.exc_info(), msg, obj=filename)
        if transformed:
            for mol, atoms in zip(molecules, atomSets):
                atoms.coords = oldCoords[mol]


    def checkArguments (self, nodes, filename, sort=True, transformed=True,
                 pdbRec=['ATOM', 'HETATM', 'MODRES', 'CONECT'],
                 bondOrigin=('File','UserDefined'), ssOrigin=None):
        """
        \nRequired Argument:\n
        nodes --- TreeNodeSet holding the current selection
        \nOptional Arguments:\n
        filename --- name of the PDB file\n
        sort --- Boolean flag to either sort or not the nodes before writing
                  the PDB file (default=True)\n
        pdbRec --- List of the PDB Record to save in the PDB file.
                  (default: ['ATOM', 'HETATM', 'MODRES', 'CONECT']\n
        bondOrigin --- string or list specifying the origin of the bonds to save
                    as CONECT records in the PDB file. Can be 'all' or a tuple\n

        ssOrigin --- Flag to specify the origin of the secondary structure
                    information to be saved in the HELIX, SHEET and TURN
                    record. Can either be None, File, PROSS or Stride.\n
        transformed --- transformed canbe True or False (default is True)
        """
        assert  isinstance(filename, str), "File names have to be strings"
        assert isinstance (pdbRec , (list, tuple))
        # ??? bondOrigin according to the doc string is a list, tuple or string
        #     how can it be 0 or 1 ???
        if bondOrigin == 0:
            bondOrigin = ('File', 'UserDefined')
        elif bondOrigin == 1:
            bondOrigin = 'all'
        ##
        assert isinstance(bondOrigin, (list, tuple, str))
        assert ssOrigin in [None, "File", "PROSS" , "Stride"]
        assert sort in (True, False, 0, 1)
        args = (nodes, filename)
        kw = {'sort':sort, 'bondOrigin':bondOrigin, 'pdbRec':pdbRec,
              'ssOrigin':ssOrigin, 'transformed':transformed}
        return args, kw




commandClassFromName = {
    'readMolecules' : [MoleculesReader, None],
    #'writePDB':  [PDBWriter, None],
    #'readAny': [ReadAny, None],
    #'readPmvSession' : [ReadPmvSession, None],
    'fetch': [Fetch, None]
    }

    


def initModule(app):
    print "initModule", app
    for cmdName, values in commandClassFromName.items():
        cmdClass, guiInstance = values
        app.addCommand(cmdClass(), cmdName, None)
