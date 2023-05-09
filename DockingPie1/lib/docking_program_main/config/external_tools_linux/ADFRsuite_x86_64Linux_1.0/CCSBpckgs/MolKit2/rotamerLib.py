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
# Copyright: M. Sanner TSRI 2015
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/MolKit2/rotamerLib.py,v 1.1.1.1.4.1 2017/07/26 22:03:40 annao Exp $
#
# $Id: rotamerLib.py,v 1.1.1.1.4.1 2017/07/26 22:03:40 annao Exp $
#
import MolKit2, os

class RotamerLib:
    """
This class load the backbone-dependent rotamer library published in
Protein Science, 6, 1661-1681 (1997) by Dunbrack.
ttp://dunbrack.fccc.edu/bbdep/bbind02.May.lib.gz
l
Usage:
     lib = RotamerLib()
     atomNamesForChiAngles = lib.getAngleDef('ARG')
     listOfChiAngles = lib.getAngles('ARG')
     listOfAngleDevistations = lib.getAnglesDev('ARG')
"""

    residueNames = ['ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'HIS', 'ILE',
                    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
                    'TYR', 'VAL']

    
    def __init__(self):
        """Constructor
lib <- RotamerLib()
"""
        # angleDef is a dict with key the amino acid type i.e. ARG and the value
        # is a list of CHI angle definition. Each CHI angle definition contains
        # a list of 4 atoms defining the torsion followed by a list of atoms
        # moved by this CHI angle
        self.angleDef = {}  # [ [A1 A2 A3 A4], [A5, A6 ...] ]

        # angleList is a list of rotameric angles, each rotameric conformation
        # contains a list of N floating point values for N CHI angles. The
        # values are given in degrees
        self.angleList = {}

        # angleList is a list of rotameric angles deviations, each rotameric
        # conformation contains a list of N floating point values for N CHI 
        # angles. The values are given in degrees
        self.angleStddev = {}

        self.loadRotamerLib()


    def loadRotamerLib(self):
        """
load backbone-independent rotamer library May 15, 2002
http://dunbrack.fccc.edu/bbdep/bbind02.May.lib.gz
libFile: bbind02.May.lib, all descriptive lines begin with '#'
defFile: define the chi angles and the atoms that are moved by the chi rotation

Roland L. Dunbrack, Jr., Ph. D.
Institute for Cancer Research, Fox Chase Cancer Center                  
7701 Burholme Avenue, Philadelphia PA 19111
"""
        rotlib = self.angleList
        angdef = self.angleDef
        angledev = self.angleStddev

        ##
        ## parse angle definiton file
        defFile = os.path.join(MolKit2.__path__[0], 'Data',  'rotamer.def')
        input_file = open(defFile, 'r')
        lines = input_file.readlines()

        for line in lines:
            words = line.split()
            if len(words) == 0: continue
            if words[0][0] == '#': continue                
            name = words[0]
            atomsInChiAngle = words[1:5]
            atomsMovedByChiAngle = words[5:]
            if angdef.has_key(name):
                angdef[name].append([atomsInChiAngle, atomsMovedByChiAngle])
            else:
                angdef[name] = [[atomsInChiAngle, atomsMovedByChiAngle]]

        ##
        ## parse angle values and standard deviations list
        libFile = os.path.join(MolKit2.__path__[0], 'Data',
                               'backboneDependentRotamers.lib')
        input_file = open(libFile)
        lines = input_file.readlines()

        for line in lines:
            words = line.split()
            if len(words) ==0: continue
            if words[0][0] == '#': continue                
            name = words[0]     
            # hardwired code for parsing Dunbrack's lib
            number = (len(words) - 11)/2 # number of CHI angles
            chi_angles=[]
            stdev = []
            for i in range(number):
                chi_angles.append(float(words[11 + i*2]))
                stdev.append(float(words[12+ i*2]))
            if rotlib.has_key(name):
                rotlib[name].append(chi_angles)
                angledev[name].append(stdev)                
            else:
                rotlib[name]= [chi_angles]
                angledev[name]= [stdev]


    def getAngleDef(self, residueType):
        """
returns a list of 2 tuples in which the first list contains 4 atom names
defining the torsion angle and the second list contains a list of atoms moved
by this torsion.
For ALA and GLY 2 empty lists are returned
"""
        if residueType.upper() in ['ALA', 'GLY']:
            return [[[],[]]]
        elif residueType.upper() in self.residueNames:
            return self.angleDef[residueType]
        else:
            raise ValueError("%s is not a correct residue type, expect %s"%
                             (residueType, self.residueNames))

    
    def getAngles(self, residueType):
        """
Returns a list of rotameric angles. The length of the list corresponds to the
number of rotamers. Each list contains a list of CHI angles for a given rotamer
"""
        if residueType.upper() in ['ALA', 'GLY']:
            return []
        elif residueType.upper() in self.residueNames:
            return self.angleList[residueType]
        else:
            raise ValueError("%s is not a correct residue type, expect %s"%
                             (residueType, self.residueNames))

        
    def getAnglesDev(self, residueType):
        """
Returns a list of rotameric anglesdeviations.
The length of the list corresponds to the number of rotamers.
Each list contains a list of CHI angles deviations for a given rotamer
"""
        if residueType.upper() in ['ALA', 'GLY']:
            return []
        elif residueType.upper() in self.residueNames:
            return self.angleStddev[residueType]
        else:
            raise ValueError("%s is not a correct residue type, expect %s"%
                             (residueType, self.residueNames))


    def nbRotamers(self, residueType):
        if residueType.upper() in ['ALA', 'GLY']:
            return 0
        elif residueType.upper() in self.residueNames:
            return len(self.angleList[residueType])
        else:
            raise ValueError("%s is not a correct residue type, expect %s"%
                             (residueType, self.residueNames))



if __name__=='__main__':
    lib = RotamerLib()
    for name in lib.residueNames:
        print name, lib.nbRotamers(name)

    print 'GLY', lib.nbRotamers('Gly')
    print 'ALA', lib.nbRotamers('ala')

    try:
        lib.nbRotamers('test')
    except ValueError:
        pass

    for chi in lib.getAngleDef('ASN'):
        print 'def: %s Moving: %s'%(chi[0], chi[1])

    angles = lib.getAngles('ASN')
    dev = lib.getAnglesDev('ASN')
    for a, d in zip(angles, dev):
        print "angles: %7.2f %7.2f dev:  %6.2f %6.2f"%(
            a[0], a[1], d[0], d[1])

    from MolKit2 import Read
    mol = Read('1jff.pdb')
    #hv = mol._ag.getHierView()
    resAtoms = mol.select('chain A resnum 2')
    angles, resname = mol.measureCHIs('A', 2, lib)
    rotind, diff, rotAngles, rotSigma = mol.closestRotamerCHI(angles, resname, lib)
    print 'angles      ', angles
    print 'diff        ', diff
    print 'closest rot', rotind,
    print 'rot Angles  ', rotAngles
    print 'deviation   ', rotSigma
