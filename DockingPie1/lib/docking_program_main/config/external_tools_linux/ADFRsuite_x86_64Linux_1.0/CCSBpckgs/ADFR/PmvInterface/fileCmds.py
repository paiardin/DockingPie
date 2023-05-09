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
# Copyright: M. Sanner TSRI 2016
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/ADFR/PmvInterface/fileCmds.py,v 1.3 2016/12/07 00:38:32 sanner Exp $
#
# $Id: fileCmds.py,v 1.3 2016/12/07 00:38:32 sanner Exp $
#
import os
from PmvApp.Pmv import MVCommand
from PmvApp.Pmv import AfterAddMoleculeEvent
from MolKit2.selection import SelectionSet
from ADFR.dro import DockingResultsObject

class DockingResultReader(MVCommand):
    """Command to read Docking Result files produced by ADFR (.dro)

    
    addToRecent --- if set to True, adds to the list of applivcation recent files;
    group --- can be None, a MoleculeGroup instance or a name(string) of existing group, if specified , the molecules are added to the group.
    """
    def __init__(self):
        MVCommand.__init__(self)
        self._dros = {}
        
    def checkArguments(self, filenames, addToRecent=True, group=None):

        assert isinstance(filenames, (list, tuple))
        for name in filenames:
            assert isinstance(name, str), "File names have to be strings"
        assert addToRecent in [0,1,True, False], 'got %s'%addToRecent
        args = (filenames,)
        if group is not None:
            from MolKit2.molecule import MoleculeGroup
            assert isinstance (group ,(MoleculeGroup, str))
        kw = {'addToRecent':addToRecent, 'group':group}
        return args, kw

    def expandArg0(self, obj):
        if isinstance(obj, list) : return obj
        else: return [obj]

    def doit(self, filename, addToRecent=True, group=None):
        
        # read dockign result file
        dro = DockingResultsObject()
        dro.load(filename)
        # build ADFR scorer
        dro.makeScorer()

        # delete covlent ligand atoms is any
        clig = dro._rec.select('segment CLIG')
        if clig:
            dro._rec._ag._flags['deleted'][clig.getIndices()] = [True]*len(clig)
        
        # buidl atom sets
        flex = dro._solutions
        rigid = dro._rec

        pmv = self.app()
        # FIXME PmvGUI onNewDisplayMol should be triggered
        name = os.path.splitext(os.path.basename(filename))[0]
        droGroup = pmv.createGroup('%s 1/%d'%(name, dro._solutions.numMols()))
        droGroup._multi = 'conformations'
        droGroup._basename = 'DRO'

        # add ligand molecule
        newmol = dro._solutions
        pmv.addMolecule(newmol, group=droGroup)
        pmv.applyDefaultCommands(newmol, "Molecule")
        event = AfterAddMoleculeEvent(molecule=newmol)
        pmv.eventHandler.dispatchEvent(event)

        # add receptor molecule
        newmol = dro._rec
        pmv.addMolecule(newmol, group=droGroup)
        pmv.applyDefaultCommands(newmol, "Molecule")
        event = AfterAddMoleculeEvent(molecule=newmol)
        pmv.eventHandler.dispatchEvent(event)
        
        if dro._flexRes:
            # add named set for flexible residue atoms
            pmv.addNamedSelection(SelectionSet([dro._FRAtoms],
                                               'flexibleReceptor'),
                                  'flexibleReceptor',
                                  droGroup)

        hbo = pmv.displayHB(dro._solutions, dro._rec)
        
        self._dros[filename] = [dro, hbo]

        return [dro, hbo]

commandClassFromName = {
    'loadDRO': [DockingResultReader, None]
    }
