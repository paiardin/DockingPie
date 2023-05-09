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
# Author: Sophie COON, Michel F. SANNER, Anna Omelchenko
#
# Copyright: M. Sanner TSRI 2014
#
#############################################################################

# $Header: /mnt/raid/services/cvs/PmvApp/extrusionCmds.py,v 1.2.4.1 2017/07/13 20:55:28 annao Exp $
#
# $Id: extrusionCmds.py,v 1.2.4.1 2017/07/13 20:55:28 annao Exp $
#

from PmvApp.Pmv import MVCommand
from PmvApp.extruder import Sheet2D
from MolKit2.molecule import AtomSet
from MolKit2.protein import Residue, ResidueSet
from DejaVu2.Spheres import Spheres


class ComputeSheet2DCommand(MVCommand):
    """The ComputeSheet2DCommand class implements methods to compute the sheet2D for each chain in the current selection. Need to specify two control atoms. \n
    Package : PmvApp \n
    Module  : extrusionCmds \n
    Class   : ComputeSheet2DCommand \n
    Command name : computeSheet2D \n
    Synopsis:\n
        None <- ComputeSheet2D(nodes, sheet2DName, ctlAtmName, torsAtmName,
                               buildIsHelix=1, nbrib=2,
                               nbchords=4, width=1.5, offset=1.2) \n
    Required Arguments:\n   
        nodes --- MolKit2 set of nodes describing molecular components \n
        ctlAtmName --- name of first control atom; \n
        torsAtmName --- name of the second atom is used to limit the torsion of the sheet2D. \n
        sheet2DName --- name of Sheet2D , used as a keyword in a dictionary containing all built Sheet2D objects;
    Optional Arguments:\n
        buildIsHelix --- flag, when set to 1 specifies that a helix is defined. \n
        nbrib --- number of ribbons of the sheet2D. \n
        nbchords --- number of points per residue \n
        width  --- width of the sheet2D (float) \n
        offset --- offset of ? (float)  \n

    Required Packages:\n
      MolKit, DejaVu2, mglutil, OpenGL\n
    Examples:\n
      mol = mv.Mols[0]  \n
      mv.extrudeCATrace(mv.getSelection())
    
    """
    def __init__(self):
        MVCommand.__init__(self)

        
    def checkArguments(self, nodes, sheet2DName, ctlAtmName, torsAtmName,
                 buildIsHelix=False, nbrib=2, nbchords=4, width=1.5, offset=1.2
                 ):
        """None <- ComputeSheet2D(nodes, sheet2DName, ctlAtmName, torsAtmName,
                               buildIsHelix=False, nbrib=2, nbchords=4,
                               width=1.5, offset=1.2) \n
        nodes  --- any set of MolKit2.Selection describing molecular components;) \n
        ctlAtmName --- name of first control atom; \n
        torsAtmName --- name of the second atom is used to limit the torsion of the sheet2D. \n
        sheet2DName --- name of Sheet2D , used as a keyword in a dictionary containing all built Sheet2D objects; \n
        buildIsHelix ---- Boolean flag when set to True specifies that a helix is defined. \n
        nbrib --- integer specifying the number of ribbon of the sheet2D.  \n
        nbchords --- integer specifying the number of points per residue \n
        width  --- float specifying the width of the sheet2D  \n
        offset --- float specifying the offset of ?
        """
        if isinstance(nodes, str):
            self.nodeLogString = "'"+nodes+"'"
        nodes = self.app().expandNodes(nodes)
        assert isinstance (sheet2DName, str)
        assert isinstance (ctlAtmName, str)
        assert isinstance (torsAtmName, str)
        assert isinstance(nbchords, int)
        assert isinstance(width, (int, float))
        assert buildIsHelix in [1, 0, True, False]
        kw = {}
        kw['buildIsHelix'] = buildIsHelix
        kw['nbrib'] = nbrib
        kw['nbchords'] = nbchords
        kw['width'] = width
        kw['offset'] = offset
        return (nodes, sheet2DName, ctlAtmName, torsAtmName), kw
        
    def getSheet2DRes(self, chain, ctlAtmName, torsAtmName, buildIsHelix=False):
        isHelix = []
        sheetCoords = []
        sheet2DRes = ResidueSet()
        residues = chain.residues
        for res in residues:
            hasAtm, rCoords = res.getAtmsAndCoords([ctlAtmName, torsAtmName])
            if hasAtm == 0:
                ## MS June 2012: the code below cause 2X2K.pdb and 2IVU.pdb
                ## to only show a short helical part. Instead of stopping
                ## we continue
                continue
##                 if len(sheet2DRes)== 0 or res.atoms[0].hetatm:
##                     # if the residue without CA and O is at the beginning
##                     # go on until you find the right residue
##                     # or res.atoms[0].hetatm added ti fix bug#779
##                     continue
##                 else:
##                     # Stop the sheet2D right there
##                     return  sheet2DRes, sheetCoords, isHelix
                ## end MS June 2012:
            else:
                sheet2DRes.append(res)
                sheetCoords = sheetCoords + rCoords
                if buildIsHelix and hasattr(res, 'secondarystructure'):
                    sr = res.secondarystructure.structureType
                    isHelix.append(sr=='Helix')
                else:
                    isHelix.append(0)
                    
        return sheet2DRes, sheetCoords, isHelix


    def doit(self, nodes, sheet2DName, ctlAtmName, torsAtmName,
             buildIsHelix=False, nbrib=2, nbchords=8, width=1.5, offset=1.2):
        """Compute the sheet2D elements for a molecule.
        1- Get the residues with a CA and O and the coordinates of these
        atoms.
        """
        if len(nodes)==0: return

        # Get the molecules having at least one node in the selection

        # loop over the chains of the molecules. We want to compute
        # the sheet2D for all the molecules not only the chain in the
        # current selection.
        from MolKit2.protein import Chain
        chains = nodes.findType(Chain).uniq()
        for chain in chains:
            # Do not recompute the sheet2D if a entry with the same
            # sheet2DName exists in the dictionary and the sheet2D
            # has been computed with the same chords, ctlAtms...
            if not hasattr(chain, 'sheet2D'):
                chain.sheet2D = {}

            if not chain.sheet2D.has_key(sheet2DName):
                chain.sheet2D[sheet2DName] = Sheet2D()
            else:
                # If one of the parameter of the sheet2D is different then...
                pass

            sheet2DRes, sheetCoords, inHelix = self.getSheet2DRes(
                chain, ctlAtmName, torsAtmName, buildIsHelix)
            
            if sheetCoords is None or len(sheetCoords)<=2:
                chain.sheet2D[sheet2DName] = None
                continue

            s = chain.sheet2D[sheet2DName]
            s.compute(sheetCoords, inHelix, nbchords=nbchords,
                      nbrib=nbrib, width=width, offset=offset )
            s.resInSheet = sheet2DRes



class Nucleic_Acids_properties(MVCommand):
    """The Nucleic_Acids_properties class implements methods for setting
    Nucleic Acids colors and scaling factor. \n
    Package : PmvApp \n
    Module  : extrusionCmds \n
    Class   : Nucleic_Acids_properties \n
    Command name : Nucleic_Acids_properties
    """
    def __init__(self):
        "Constructor for Nucleic_Acids_properties"
        MVCommand.__init__(self)
        self.color_A =  [1,0,0]
        self.color_C = [1,1,0]
        self.color_U = [1,0.5,0]
        self.color_G = [0,0,1]
        self.color_T = [0,1,0]
        self.scale_purine = self.scale_pyrimidine = 1.3
        self.height_purine = self.height_pyrimidine = 0.4
        self.color_backbone = 1


    def checkArguments(self, color_A=[1,0,0], color_G=[0,0,1], color_T=[0,1,0],
                 color_C=[1,1,0], color_U=[1,0.5,0],
                 scale_purine = 1.3, scale_pyrimidine = 1.3, 
                 height_purine = 0.4, height_pyrimidine = 0.4, 
                 color_backbone = 1):
        """None <- Nucleic_Acids_properties(color_A =  [1,0,0], 
                 color_G = [0,0,1], color_T = [0,1,0],
                 color_C = [1,1,0], color_U = [1,0.5,0],
                 scale_purine = 1.3, scale_pyrimidine = 1.3, 
                 height_purine = 0.4, height_pyrimidine = 0.4, 
                 color_backbone = 1)
        """
        assert isinstance(color_A, (list, tuple))
        assert len(color_A) == 3
        assert isinstance(color_G, (list, tuple))
        assert len(color_G) == 3
        assert isinstance(color_T, (list, tuple))
        assert len(color_T) == 3
        assert isinstance(color_C, (list, tuple))
        assert len(color_C) == 3
        assert isinstance(color_U, (list, tuple))
        assert len(color_U) == 3
        assert isinstance (scale_purine, (int, float))
        assert isinstance (scale_pyrimidine, (int, float))
        assert isinstance (height_purine, (int, float))
        assert isinstance (height_pyrimidine, (int, float))
        kw = {}
        kw['color_A'] = color_A
        kw['color_G'] = color_G
        kw['color_C'] = color_C
        kw['color_T'] = color_T
        kw['color_U'] = color_U
        kw['scale_purine'] = scale_purine
        kw['scale_pyrimidine'] = scale_pyrimidine
        kw['height_purine'] = height_purine
        kw['height_pyrimidine'] = height_pyrimidine
        kw['color_backbone'] = color_backbone
        return (), kw
        

    def doit(self, color_A=[1,0,0], color_G=[0,0,1], color_T=[0,1,0],
             color_C=[1,1,0], color_U=[1,0.5,0], 
             scale_purine = 1.3, scale_pyrimidine = 1.3, 
             height_purine = 0.4, height_pyrimidine = 0.4, 
             color_backbone = 1):

        self.color_A = color_A
        self.color_C = color_C
        self.color_U = color_U
        self.color_G = color_G
        self.color_T = color_T
        self.scale_purine = scale_purine
        self.scale_pyrimidine = scale_pyrimidine
        self.height_purine = height_purine
        self.height_pyrimidine = height_pyrimidine

        self.color_backbone = color_backbone

        

        
commandClassFromName = {
    'computeSheet2D' : [ComputeSheet2DCommand, None],
    'Nucleic_Acids_properties' : [Nucleic_Acids_properties, None],
}


def initModule(viewer, gui=False):
    for cmdName, values in commandClassFromName.items():
        cmdClass, guiInstance = values
        viewer.addCommand(cmdClass(), cmdName, guiInstance)
