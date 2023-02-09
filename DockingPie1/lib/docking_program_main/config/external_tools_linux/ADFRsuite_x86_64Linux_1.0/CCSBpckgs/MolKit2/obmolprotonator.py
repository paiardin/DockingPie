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
# Date: 2015 Authors: Stefano Forli
#
#    forli@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Stefano Forli and TSRI 2015
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/MolKit2/obmolprotonator.py,v 1.1.1.1.4.1 2017/07/26 22:03:40 annao Exp $
#
# $Id: obmolprotonator.py,v 1.1.1.1.4.1 2017/07/26 22:03:40 annao Exp $
#
import openbabel as ob



def generateHPdbRecords(mol):
    obconv = ob.OBConversion()
    obconv.SetOutFormat('pdb')
    molString = obconv.WriteString(mol)
    buff = []
    for l in molString.split('\n'):
        if l[0:4] in ['ATOM', 'HETA'] and l[77] == 'H':
            buff.append(l)
    return buff
    #  return [ x for x in molString.split() if x.startswith('ATOM') and x[77] == 'H']


class OBMolProtonator:
    """ """

    def __init__(self):
        """ """

    def addHydrogens(self, mol, pH=None, atomList=[], perceiveBondOrder=True):
        """
            pH          :   protonate the molecule (NOTE set also force to True)
            force       :   remove all previous hydrogens if any
                            otherwise only missing will be added
            atomList    :   protonate only heavy atoms in this list
            perceiveBondOrder: enable bond order perception
        """
        if perceiveBondOrder: 
            mol.AssignSpinMultiplicity(False)
            mol.PerceiveBondOrders()

        if len(atomList):
            for aIdx in atomList:
                atom = mol.GetAtom(aIdx)
                mol.AddHydrogens(atom)
            if perceiveBondOrder:
                mol.PerceiveBondOrders()
            mol.ConnectTheDots()
            molNew = mol
        else:
            molNew = ob.OBMol()
            fragments = mol.Separate()
            for f in fragments:
                f.DeleteHydrogens()
                if pH == None:
                    f.AddHydrogens(False, False)
                else:
                    # from Avogadro hydrogensextension.cpp
                    f.UnsetFlag(ob.OB_PH_CORRECTED_MOL)
                    for a in ob.OBMolAtomIter(f):
                        a.SetFormalCharge(0)
                    f.SetAutomaticFormalCharge(True)
                    f.AddHydrogens(False, True, pH)
            molNew += f                        
        self.fixHydrogenNames(molNew)
        return molNew

    def fixHydrogenNames(self, mol):
        """ assing PDB names to hydrogens that could have been added
            and dont have one.
            It doesn't generate duplicates but it doesn't enforce unique names!
            easy to fix, though
        """
        c = 0
        for res in ob.OBResidueIter(mol):
            for atom in ob.OBResidueAtomIter(res):
                name = res.GetAtomID(atom)
                if name.strip() == '':
                    name = self.eTable.GetSymbol(atom.GetAtomicNum())
                res.SetAtomID(atom, name)
                c+=1

if __name__ == '__main__':
    from sys import argv
    import os
    #from ligProtonator import OBMolProtonator
    import pybel
    ob = pybel.ob

    def OBMolFromFilename(filename):
      obmol = ob.OBMol()
      obconv = ob.OBConversion()
      name, ext = os.path.splitext(filename)
      ext = ext[1:].lower()
      obconv.SetInFormat(ext)
      obconv.ReadFile(obmol, filename)
      return obmol, ext

    # test the code with
    # $ pythonsh ligProtonator.py molecule.ext  [any supported OB format]

    w = ob.OBConversion()
    mol,ext = OBMolFromFilename(argv[1])
    w.SetInAndOutFormats(ext, 'pdb')
    p = OBMolProtonator()


    # single Atom
    mol = p.addHydrogens(mol, atomList=[5,6])
    w.WriteFile(mol, 'test_ATOM-5,6.pdb')    
    print "Result saved as 'test_ATOM-5,6.pdb'"


    # no Ph
    mol = p.addHydrogens(mol)
    w.WriteFile(mol, 'test_NOPH.pdb')    
    print "Result saved as 'test_NOPH.pdb'"
    # pH 7.4 (since the molecule is the same, pH test must be made last, 'cause it alters atype perception)

    mol = p.addHydrogens(mol, pH=7.4)
    w.WriteFile(mol, 'test_pH7.4.pdb')    

"""tested with this input ligand

HETATM    1  O   LIG1    1       1.143   0.492  -0.419  0.00  0.00           O  
HETATM    2  C   LIG1    1       2.074  -0.137   0.178  0.00  0.00           C  
HETATM    3  O   LIG1    1       2.364  -1.365   0.014  0.00  0.00           O  
HETATM    4  C   LIG1    1       2.912   0.646   1.172  0.00  0.00           C  
HETATM    5  C   LIG1    1       2.666   2.009   1.408  0.00  0.00           C  
HETATM    6  C   LIG1    1       3.437   2.730   2.323  0.00  0.00           C  
HETATM    7  C   LIG1    1       4.468   2.099   3.017  0.00  0.00           C  
HETATM    8  C   LIG1    1       4.728   0.748   2.796  0.00  0.00           C  
HETATM    9  C   LIG1    1       3.956   0.028   1.881  0.00  0.00           C  
END

"""
