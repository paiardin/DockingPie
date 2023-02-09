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
# $Header: /mnt/raid/services/cvs/MolKit2/protonate.py,v 1.1.1.1.4.1 2017/07/26 22:03:40 annao Exp $
#
# $Id: protonate.py,v 1.1.1.1.4.1 2017/07/26 22:03:40 annao Exp $
#
import numpy, prody, os, subprocess, platform

from mglutil.util.packageFilePath import getBinary

from MolKit2.molecule import Molecule
from MolKit2.selection import Selection
from MolKit2.Kamaji import KamajiInterface

# helper class to make stdout set of lines look like a file that ProDy can parse
class memFile:
    def __init__(self, lines):
        self.lines = lines
    def readlines(self):
        return self.lines

class MacroMoleculeProtonator:
    """
    Class to protonate macro molecules (protein, amino and nucleic acids)
    using reduce
    """

    def __init__(self):
        system_info = platform.uname()
        if system_info[0] == 'Windows':
            self._shell=False
        else:
            self._shell=True

        self.reducePath = getBinary('reduce', 'binaries')

    def addHydrogens(self, molSel):

        # run reduce
        #print reducePath, '-build', inpFname
        proc = subprocess.Popen("%s -build -"%self.reducePath,
                                stdin=subprocess.PIPE , 
                                stdout=subprocess.PIPE , 
                                stderr=subprocess.PIPE, 
                                bufsize = 1, shell=self._shell)
        prody.writePDBStream(proc.stdin, molSel.select('not hydrogen'))

        stdout_value, stderr_value = proc.communicate('through stdin to stdout')
        output = stdout_value.split('\n')

        # read in the protonated molecule
        if len(output)==1:
            raise RuntimeError("reduced failed %s"%stderr_value)
        else:
            ag = prody.parsePDBStream(memFile(output))
            molH = Molecule('protonated', ag)
            molH.buildBondsByDistance()
            molH.defaultRadii()
        return molH

class SmallMoleculeProtonator:
    """
    Class to protonate small molecules (ligands, co-fators) using OpenBabel
    """

    pass

if __name__=='__main__':
    from MolKit2 import Read
    mol = Read('1crn.pdb')
    p = MacroMoleculeProtonator()
    molH = p.protonate(mol.select())
    
