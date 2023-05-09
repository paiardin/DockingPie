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

########################################################################
#
# Date: 2015 Authors: Michel Sanner
#
#    sanner@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Michel Sanner and TSRI 2015
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/PmvApp/interactionsCmds.py,v 1.1.4.1 2017/07/13 20:52:52 annao Exp $
#
# $Id: interactionsCmds.py,v 1.1.4.1 2017/07/13 20:52:52 annao Exp $
#
import numpy, prody, os, tempfile, subprocess, platform

from PmvApp.Pmv import MVCommand, formatName, AddAtomsEvent#, AfterAddMoleculeEvent

from mglutil.util.packageFilePath import getBinary

from MolKit2.molecule import Molecule
from MolKit2.selection import Selection, SelectionSet
from MolKit2.Kamaji import KamajiInterface
from MolKit2.protonate import MacroMoleculeProtonator
from MolKit2.obmolprotonator import OBMolProtonator, generateHPdbRecords
from MolKit2.Kamaji import KamajiInterface
from MolKit2.openBabelInterface import ProdyToOBMol, OBMolToPrody

from mglutil.util.io import BufferAsFile

from prody.measure.contacts import findNeighbors

class ProtonateBase(MVCommand):

    def firstArgStringToMol(self, selection):
        
        if isinstance (selection, str):
            self.app().Mols.selectMolKit(selection)
        else:
            return selection
    
class ProtonateWithReduce(ProtonateBase):
    """
    protonate molecule using reduce
    """
    def onAddCmdToApp(self):
        ProtonateBase.onAddCmdToApp(self)
        if platform.uname()[0] == 'Windows':
            self._shell=False
        else:
            self._shell=True
        # get the path to the reduce binary
        self._reducePath = getBinary('reduce', 'binaries')

    def checkArguments(self, molSel, flip=True, redraw=True):
        """
        \nRequired Arguments:\n
        molSel --- MolKit2 selection or selection string\n
        """
        kw = {}
        assert flip in [True, False, 0, 1]
        kw['flip'] = flip
        kw['redraw'] = redraw
        
        assert isinstance(molSel, SelectionSet)

        return (molSel,), kw

    def doit(self, molSel, flip=True, redraw=True):
        ## to protonate this selection

        # FIXME, we have to find out which reduce options to support
        opt = '-build'
        if flip:
            opt = '-build'

        noH = molSel.select('not hydrogen')

        # run reduce
        print "%s %s -"%(self._reducePath, opt)
        proc = subprocess.Popen("%s %s -"%(self._reducePath, opt),
                                stdin=subprocess.PIPE , 
                                stdout=subprocess.PIPE , 
                                stderr=subprocess.PIPE, 
                                bufsize = -1, shell=self._shell)

        # write PDB records to stdin
        prody.writePDBStream(proc.stdin, noH)

        # get PDB records for molecule protonated by reduce
        stdout_value, stderr_value = proc.communicate('through stdin to stdout')
        lines = stdout_value.split('\n')
        if len(lines)==1:
            raise RuntimeError("reduced failed %s"%stderr_value)

        # find added H atoms records and flipped atoms
        Hrec = []
        flips = []
        maxAtNum = max(noH.getSerials())+1
        for line in lines:
            if flip and line[-4:] == 'flip':
                flips.append(line)
            elif (line[0:4]=='ATOM' or line[0:4]=='HETA') and line[77] == 'H':
                Hrec.append(line)

        #if flip:
        #    self.handleFlips(flips)

        # make a molecule with H atoms
        agH = prody.parsePDBStream(BufferAsFile(Hrec))
        agH.setData('atomicNumber', numpy.ones( (len(agH),), 'int'))
        agH.setRadii(numpy.ones( (len(agH),), 'f')*1.2)
        
        # extend the molecule's atom group with hydrogen atoms
        mol = molSel.getAtomGroup().getMolecule()
        lenBefore = len(mol._ag)
        nbBondsBefore = len(mol._ag._bonds)
        #import pdb; pdb.set_trace()
        # add the H atoms
        mol._ag.extend(agH)

        # add the color white and set all H atom colorindices to white
        for name in mol._colors.keys():
            if name=='msms':
               for msname in mol._colors['msms'].keys():
                   mol._colors['msms'][msname] = numpy.concatenate( (mol._colors['msms'][msname], [(1.,1.,1.,1.)]) )
                   cols = mol._ag._data['colorsIndices_msms_%s'%msname]
                   cols[lenBefore:] = len(mol._colors['msms'][msname])-1
            else:
                mol._colors[name] = numpy.concatenate( (mol._colors[name], [(1.,1.,1.,1.)]) )
                cols = mol._ag._data['colorsIndices_%s'%name]
                cols[lenBefore:] = len(mol._colors[name])-1

        # add bonds
        addedAtoms = Selection(mol._ag, range(len(mol._ag)-len(agH),len(mol._ag)), '')

        # find heavy atoms to which H are bound
        atomPairs = findNeighbors(mol.select(), 1.0, addedAtoms)
        bonds = mol._ag._bonds.tolist()
        bo = [mol._ag._bondOrder["%d %d"%(b[0], b[1])] for b in bonds]
        for a1, a2, dist in atomPairs:
            bonds.append( (a1.getIndex(), a2.getIndex()) )
            if bo:
                bo.append(1)
        mol._ag.setBonds(bonds, bo)

        pmv = self.app()

        # handle will replace _ag in mol.geomContainers.atoms and mol.geomContainers.geoms
        event = AddAtomsEvent(molecule=mol, atoms=addedAtoms)
        self.app().eventHandler.dispatchEvent(event)

        # update Lines, CPK, S&B and MSMS geometries by adding the H atoms to the set
        # of atoms displayed for this representation only if the heavy atom is shown
        # in this representation.
        lineAtoms = mol.geomContainer.atoms.get('lines', None)
        if lineAtoms:
            # find h atom for heavy displayed as lines
            atomList = findNeighbors(lineAtoms, 1.0, addedAtoms)
            linH = Selection(mol._ag, [x[1].getIndex() for x in atomList], '')
            mol.geomContainer.atoms['lines'] = lineAtoms | linH
            pmv.displayLines.refreshDisplay(mol)

        sbAtoms = mol.geomContainer.atoms.get('sb', None)
        if sbAtoms:
            # find h atom for heavy displayed as lines
            mol._ag._bondData['sb_cylRadius'][nbBondsBefore:] = 0.2
            mol._ag._data['sb_ballsRadius'][lenBefore:] = 0.3
            atomList = findNeighbors(sbAtoms, 1.0, addedAtoms)
            sbH = Selection(mol._ag, [x[1].getIndex() for x in atomList], '')
            mol.geomContainer.atoms['sb'] = sbAtoms | sbH
            pmv.displaySB.refreshDisplay(mol)

        cpkAtoms = mol.geomContainer.atoms.get('cpk', None)
        if cpkAtoms:
            # find h atom for heavy displayed as lines
            mol._ag._data['cpk_radius'][lenBefore:] = 1.2
            atomList = findNeighbors(cpkAtoms, 1.0, addedAtoms)
            cpkH = Selection(mol._ag, [x[1].getIndex() for x in atomList], '')
            mol.geomContainer.atoms['cpk'] = cpkAtoms | cpkH
            pmv.displayCPK.refreshDisplay(mol)

        if hasattr(mol, '_msmsData'):
            for name in mol._msmsData['msms'].keys():
                pmv.computeMSMS(mol)
                msAtoms = mol.geomContainer.atoms.get('msms_'+name, None)
                if msAtoms:
                    # find h atom for heavy displayed as lines
                    atomList = findNeighbors(msAtoms, 1.0, addedAtoms)
                    msH = Selection(mol._ag, [x[1].getIndex() for x in atomList], '')
                    mol.geomContainer.atoms['msms_'+name] = msAtoms | msH
                    pmv.displayMSMS.refreshDisplay(mol)

        return addedAtoms

         
class ProtonateWithOpenBabel(ProtonateBase):

    def onAddCmdToApp(self):
        if platform.uname()[0] == 'Windows':
            self._shell=False
        else:
            self._shell=True
        # get the path to the reduce binary
        self.OBProtonator = OBMolProtonator()

    def checkArguments(self, molSel, pH=None, force=False, atomList=[],
                       perceiveBondOrder=True, redraw=True):
        """
        \nRequired Arguments:\n
        molSel --- MolKit2 selection or selection string\n

        \nOptional Arguments:\n
        pH          :   protonate the molecule (NOTE set also force to True)
        force       :   remove all previous hydrogens if any
                        otherwise only missing will be added
        atomList    :   protonate only heavy atoms in this list
        perceiveBondOrder: enable bond order perception
        """
        kw = {}
        assert pH is None or pH > 0.0 and pH < 14.0
        assert force in [True, False, 0, 1]
        assert isinstance(atomList, list)
        assert perceiveBondOrder in [True, False, 0, 1]
        kw['pH'] = pH
        kw['force'] = force
        kw['atomList'] = atomList
        kw['perceiveBondOrder'] = perceiveBondOrder
        kw['redraw'] = redraw
        
        assert isinstance(molSel, SelectionSet)

        return (molSel,), kw

    def doit(self, molSel, pH=None, force=False, atomList=[],
             perceiveBondOrder=True, redraw=True):


        ## to protonate this selection using OpenBabel
        obh = []
        # we need to proceed by connected fragments
        l = len(molSel) # i.e. max fragment length
        sel = molSel.copy()
        while len(sel):
            first = molSel.getSerials()[0]
            frag = molSel.select('bonded %d to serial %d'%(l, first))
            obmol = ProdyToOBMol(frag, title='obmol')
            obmolH = self.OBProtonator.addHydrogens(
                obmol, pH=pH, atomList=atomList,
                perceiveBondOrder=perceiveBondOrder)
            obh.extend(generateHPdbRecords(obmolH))
            sel = sel - frag

        ## buid the prody molecule with all H atoms at the end
        # get molecule
        mol = molSel.getAtomGroup().getMolecule()
        # get pdb records for molecule before added H
        #rec = mol.getPDBRecords()
        ag = prody.parsePDBStream(BufferAsFile(obh))
        molH = Molecule(mol.name, ag)
        molH._ag._bonds = None # no need to call buildBondsByDistance
        molH.defaultRadii()
            
        mol = molSel.getAtomGroup().getMolecule()
        mol.addAtoms(molH.select())

        mol.geomContainer.allCoords = mol._ag.getCoords().astype('f')

        Hatoms = mol.select('hydrogen')

        # Create the AddAtomsEvent
        event = AddAtomsEvent(addedAtoms=Hatoms)
        self.app().eventHandler.dispatchEvent(event)

        return Hatoms

class Protonate(MVCommand):
    """- The protonate command will be written to always work on the entire molecule for each selection
  in the selection set passed top the command.
- The command cannot be undone (not easy to implement properly) i.e. residues flipped by reduce cannot be unflipped
- for every molecule it will run reduced which will protonate aminoacids (and nucleic acids ?)
   next it will run Kamaji to:
        - define atom types, bond orders and HBdonnors/acceptors (for the entire molecule, right?)
        - next it should run an OpenBabel-based protonator for ligands and cofactors
   After running this command, atoms will have:
   - types (atomType),
   - bonds will have bondorders save in mol.bondorder[(i,j)]
     where i,j are 0 based atom indices with i<j
   - HB donor acceptor properties ('hbType')
   """

    def onAddCmdToApp(self):
        #self.app().eventHandler.registerListener(DeleteAtomsEvent, self.updateGeom)
        #system_info = platform.uname()
        #platform = system_info[0]
        if platform.uname()[0] == 'Windows':
            self._shell=False
        else:
            self._shell=True
        # get the path to the reduce binary
        self._reducePath = getBinary('reduce', 'binaries')

        # create macro molecule protonator
        self.largeProtonator = MacroMoleculeProtonator()

        # create small moelcule protonator
        self.smallProtonator = SmallMoleculeProtonator()
        
    def updateGeom(self, event):
        pass

    def getLastUsedValues(self, formName='default', **kw):
        """Return dictionary of last used values
"""
        return self.lastUsedValues[formName].copy()

    def firstArgStringToMol(self, selection):
        
        if isinstance (selection, str):
            self.app().Mols.selectMolKit(selection)
        else:
            return selection

    
    def checkArguments(self, molSel, lineWidth=2, negate=False, redraw=True):
        """
        \nRequired Arguments:\n
        molSel --- MolKit2 selection or selection string\n

        \nOptional Arguments:\n
        lineWidth --- int specifying the width of the lines, dots or doted lines
                     representing the selection. (default = 2)\n
        negate --- boolean flag specifying whether or not to negate the
                     current selection. (default = False)\n
        """
        kw = {}
        assert isinstance(lineWidth, (int, float))
        assert lineWidth>0
        assert negate in [True, False, 0, 1]        
        kw['lineWidth'] = lineWidth
        kw['negate'] = negate
        kw['redraw'] = redraw
        
        assert isinstance(molSel, SelectionSet)

        return (molSel,), kw

    def doit(self, molSel, lineWidth=2, negate=False, redraw=True):
        ## to protonate this molecule using reduce

        mol = molSel.getAtomGroup().getMolecule()
        # get PDB reords for molecule protonated by reduce
        lines = self.largeProtonator.addHydrogens(mol)

        ligs = {}
        # loop over ligands cofactors and additives to fix protanation
        # with openbabel
        ligands = mol.select('not protein and not nucleic and not ion')
        obh = []
        if ligands:
            for res in mol.select('not protein and not nucleic and not ion').getHierView().iterResidues():
                resnum = res.getResnum()
                ligs[resnum] = 1
                frag = mol.select('resnum %d'%resnum)
                #stream = Stream()
                #prody.writePDBStream(stream, frag)
                #hrecRes = [l for l in Hrec if int(l[23:26])==resnum]
                #ag = prody.parsePDBStream(BufferAsFile(stream.lines+hrecRes))
                #pmol = Molecule('res', ag)
                #pmol.buildBondsByDistance()
                #obmol = ProdyToOBMol(pmol.select(), title='obmol')
                obmol = ProdyToOBMol(frag, title='obmol')
                p = OBMolProtonator()
                obmolH = p.addHydrogens(obmol, pH=7.4)#, pH=pH, force=force, atomList=atomList,
                                        #perceiveBondOrder=perceiveBondOrder)
                print 'NUMATONS', obmolH.NumAtoms()
                obh.extend(generateHPdbRecords(obmolH))

        #extracted H toms aded by reduce to put them at the end
        Hrec = []
        rec = []
        for line in lines:
            if (line[0:4]=='ATOM' or line[0:4]=='HETA') and line[77] == 'H':
                if not ligs.has_key(int(line[23:26])):
                    Hrec.append(line)
            else:
                rec.append(line) 

        # buid the prody molecule with all H atoms at the end
        ag = prody.parsePDBStream(BufferAsFile(rec+Hrec+obh))
        molH = Molecule(mol.name, ag)
        molH.buildBondsByDistance()
        molH.defaultRadii()
            
        self.app().deleteMolecule(mol)
        # add protonated molecule
        self.app().addMolecule(molH)

        return molH

        
        ## if not hasattr(mol, '_kamaji'):
        ##     mol._kamaji = KamajiInterface(mol)

        ## macroSelections = mol.emptySelection()
        ## for treeName in ['Ligand', 'Other']:
        ##     if getattr(mol._kamaji, 'has%s'%treeName):
        ##        for tree in treeList:
        ##            for selStr in tree.getSelectionString().split("_+_"):
        ##                if len(selStr)==0: continue
        ##                frag = self.mol.select(selStr)
        ##                if len(frag & molSel):
        ##                    print 'sel overlaps with %s %s'%(treeName, selStr)
        ##                    macroSelections = macroSelections | frag

        ## ligSelections = []
        ## for treeName in ['Ligand', 'Other']:
        ##     if getattr(mol._kamaji, 'has%s'%treeName):
        ##        for tree in treeList:
        ##            for selStr in tree.getSelectionString().split("_+_"):
        ##                if len(selStr)==0: continue
        ##                frag = self.mol.select(selStr)
        ##                if len(frag & molSel):
        ##                    print 'sel overlaps with %s %s'%(treeName, selStr)
        ##                    ligSelections.append( inter )

        return molH

    def updateSelections(self, mol, molH):
        # build renumbering map (only necessary if selections exist)
        offset = []
        off = 0
        elements = molH._ag.getElements()
        for i, el in enumerate(elements):
            if el == 'H':
                off += 1
            else:
                offset.append(off)

        ## update existing selections to point to the new molecule
        ##
        for i, sel in enumerate(pmv.curSelection):
            if sel.getAtomGroup() == mol._ag:
                newInd = []
                for ind in sel.getIndices():
                    print ind, offset[ind], ind + offset[ind]
                    newInd.append( ind + offset[ind] )
                pmv.curSelection[i] = Selection(molH._ag, newInd, sel._selstr)
                break

        for name, selSet in pmv.namedSelections.items():
            for i, sel in enumerate(selSet):
                if sel.getAtomGroup() == mol._ag:
                    newInd = []
                    for ind in sel.getIndices():
                        newInd.append( ind + offset[ind] )
                    selSet[i] = Selection(molH._ag, newInd, sel._selstr)
                    break

        # now run Kamaji to identify non amino or nucleic acids
        #interface = KamajiInterface(molH)

        # set atom types, bond orders and HB donor acceptor properties
        #molH._ag.setData('atomType', interface.getAtomTypes())
        #molH._ag.setData('hbType', interface.getHBTypes())
        #molH._ag._bondOrder = interface.getBondOrders()
        # delete unprotonated molecule
        #pmv.deleteMolecule(mol)
        # add protonated molecule
        #pmv.addMolecule(molH)
        

commandClassFromName = {
    #'protonate' : [Protonate, None],
    'protonateWithReduce' : [ProtonateWithReduce, None],
    'protonateWithOpenBabel' : [ProtonateWithOpenBabel, None],

    }


def initModule(viewer, gui=True):
    for cmdName, values in commandClassFromName.items():
        cmdClass, guiInstance = values
        viewer.addCommand(cmdClass(), cmdName, guiInstance)

