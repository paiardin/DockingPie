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
#########################################################################
#
# $Header: /mnt/raid/services/cvs/PmvApp/hbCmds.py,v 1.1.4.1 2017/07/13 20:52:52 annao Exp $
#
# $Id: hbCmds.py,v 1.1.4.1 2017/07/13 20:52:52 annao Exp $
#
import numpy
from MolKit2 import hbGeom
from DejaVu2.Cylinders import Cylinders
from DejaVu2.Spheres import Spheres
from DejaVu2.stipples import stippleLines
from PmvApp.displayCmds import DisplayCommand

class HBObj:

    def __init__(self, name, atoms1, atoms2):
        self.name = name
        self.atoms1 = atoms1
        self.atoms2 = atoms2
        self.cylinders = Cylinders('%s_cyl'%name,
                                   visible=0, inheritMaterial=False)
        self.spheres = Spheres('%s_sph'%name,
                               visible=0, inheritMaterial=False)
        self.stippleLength = 0.4
        self.stippleSpace = 0.3
        self.radius = 0.1
        self.colors = [[1,1,0,1]]

class DisplayHB(DisplayCommand):
    """display hydrogen bonds between 2 sets of atoms

    Synopsis:
        None <- displayLines(atoms1, atoms2, width=2,
                             stippleLength=0.2, stippleSpace=0.2)

    arguments:
        atoms1    --- first set of atoms
        atoms2    --- second set of atoms

        width    --- integer > 1 specifying the width of bond representation. (Default:0.2)
        stippleLength --- float specifying the length of the stipple dash (default:0.2)
        stippleSpace  --- float specifying the length of the space between stipples (default:0.2)
        
    """
    
    _argNames = ['width', 'stippleLength', 'stippleSpace']

    def checkArguments(self, nodes1, nodes2, **kw):
        for name in kw.keys():
            if name not in self._argNames:
                raise RuntimeError("%s: unrecognized keyword argument '%s', valid names are %s"%(
                    self.name, name, str(self._argNames)))
        return (nodes1, nodes2), kw

    def __init__(self):
        DisplayCommand.__init__(self)
        self.HBObjects = [] # list of object used to compute and
                            # draw HBonds
        self.molsToHbo = {} # key will be molecules and value is a list of Hbo
        self._atomSetName = "lines"

    def initializeMolForCmd(self, mol):
        """
        """
        return

    def isInitialized(self):
        return True

    def updateModelGlobals(self, hbo, width=None, stippleLength=None,
                           stippleSpace=None, **kw):

        if width is not None:
            assert isinstance(width, int) and width >=1
            hbo.radius = width

        if stippleLength is not None:
            assert isinstance(stippleLength, float) and stippleLength >0.
            hbo.stippleLength = stippleLength

        if stippleSpace is not None:
            assert isinstance(stippleSpace, float) and  stippleSpace>0.
            mol.stippleSpace = stippleSpace

    def doit(self, atoms1, atoms2, **kw):
        ##
        ##  displayLines uses 2 atomsSets 'noBond' and 'lines'
        ##
        if 'donor' not in atoms1._ag._flags.keys():
            print 'ERROR: donors/acceptor flags missing in first atom set'
            return None
        if 'donor' not in atoms2._ag._flags.keys():
            print 'ERROR: donors/acceptor flags missing in second atom set'
            return None

        sel1 = atoms1.select('not deleted and (acceptor or donor)')
        sel2 = atoms2.select('not deleted and (acceptor or donor)')

        if len(sel1)==0:
            print "no donors or acceptors in atom set 1"
            return None
        if len(sel2)==0:
            print "no donors or acceptors in atom set 2"
            return None
        
        hbo = HBObj('HbondObj_%d'%(len(self.HBObjects)+1), sel1, sel2)
        self.HBObjects.append(hbo)
        
        # FIXME name conflicts
        self.app().gui().viewer.AddObject(hbo.cylinders)
        self.app().gui().viewer.AddObject(hbo.spheres)
        
        mol1 = atoms1.getAtomGroup().getMolecule()
        if self.molsToHbo.has_key(mol1):
            self.molsToHbo[mol1].append(hbo)
        else:
            self.molsToHbo[mol1] = [hbo]

        mol2 = atoms2.getAtomGroup().getMolecule()
        if self.molsToHbo.has_key(mol2):
            self.molsToHbo[mol2].append(hbo)
        else:
            self.molsToHbo[mol2] = [hbo]

        self.updateModelGlobals(hbo, **kw)
        self.refreshDisplay(hbo=hbo)
        return hbo
    
    def refreshDisplay(self, mol=None, hbo=None):
        visible = True
        if mol :
            visible = mol.geomContainer.masterGeom.visible
        if mol and self.molsToHbo.get(mol, None) is None and hbo is None:
            return
        if hbo is not None:
            hbos = [hbo]
        else:
            hbos = self.molsToHbo[mol]

        # we should only draw HB if donor and acceptor are visible !
        for hbo in hbos:
            if visible:
                coords, faces = hbGeom(hbo.atoms1, hbo.atoms2)
                col = [0]*len(coords)
                v1, f1, fcol, vcol = stippleLines(
                    numpy.array(coords), faces, col,
                    segLen=hbo.stippleLength, spaceLen=hbo.stippleSpace)

                hbo.cylinders.Set(vertices=v1, faces=f1,
                                  radii=(hbo.radius,), visible=1, 
                                  materials = hbo.colors, tagModified=False) 
                hbo.spheres.Set(vertices=v1, radii=(hbo.radius,),
                                materials = hbo.colors,
                                visible=1, tagModified=False)
            else:
                hbo.cylinders.Set(visible=0)
                hbo.spheres.Set(visible=0)

commandClassFromName = {
    'displayHB' : [DisplayHB, None],
}

