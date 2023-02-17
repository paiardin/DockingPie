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
# $Header: /mnt/raid/services/cvs/PmvApp/colorCmds.py,v 1.8.4.3 2017/09/07 00:51:28 annao Exp $
#
# $Id: colorCmds.py,v 1.8.4.3 2017/09/07 00:51:28 annao Exp $
#
"""
This Module implements commands to color the current selection different ways.
for example:
    by atoms.
    by residues.
    by chains.
    etc ...
    
"""

import os, sys
from mglutil.util.colorUtil import ToHEX
from PmvApp.colorPalette import ColorPalette, ColorPaletteFunction
from PmvApp.Pmv import DeleteGeomsEvent, AddGeomsEvent, EditGeomsEvent
from mglutil.events import Event

import numpy

from DejaVu2.colorTool import Map, RGBARamp, RedWhiteARamp, WhiteBlueARamp,\
     RedWhiteBlueARamp

from PmvApp.Pmv import MVCommand

from MolKit2.molecule import Molecule, Atom, Residue, Chain
from MolKit2.selection import Selection, SelectionSet
from DejaVu2.colorMap import ColorMap
from opengltk.OpenGL import GL
        

class ColorCommandBase(MVCommand):
    """Base class for Pmv color commands
    """
    _argNames = ['geomsToColor', 'carbonsOnly']

    def checkArguments(self, *args, **kw):
        for name in kw.keys():
            if name not in self._argNames:
                raise RuntimeError("%s: unrecognized keyword argument '%s', valid names are %s"%(
                    self.name, name, str(self._argNames)))
        return args, kw
    
    # virtual function that has to return a list of colors for the specified atoms
    def getColors(self, atoms):
        pass

    def doit(self, atoms, geomsToColor=['all',], carbonsOnly=False):
        """None <--- color(selections, geomsToColor=['all'])"""
        if carbonsOnly:
            atoms = atoms & atoms.select("element C")

        pmv = self.app()
        if 'all' in geomsToColor:
            geomsToColor = pmv.getGeoms(SelectionSet([atoms]))
        elif '*' in geomsToColor:
            geomsToColor = pmv.getGeoms(SelectionSet([atoms]), visibleOnly=False)

        mol = atoms.getAtomGroup().getMolecule()
        gc = mol.geomContainer
        indices = atoms.getIndices()
        bonds = None
        for gName in geomsToColor:
            if gName=='lines' or gName=='noBond':
                col = mol._colors['lines'][mol._ag.getData('colorsIndices_lines').tolist()]
                oldCol = col[indices]
                self.app().pushUndoCmd( self.app().color, (atoms, oldCol, ['lines']), {})
                col[indices,:3] = self.getColors(atoms)
                colList, colIndexList = self.app().indexColors(col)
                mol._colors['lines'] = numpy.array(colList)
                mol._ag.setData('colorsIndices_lines', colIndexList)
                self.app().displayLines.refreshDisplay(mol)

            elif gName=='cpk':
                col = mol._colors['cpk'][mol._ag.getData('colorsIndices_cpk').tolist()]
                oldCol = col[indices]
                self.app().pushUndoCmd( self.app().color, (atoms, oldCol, ['cpk']), {})
                col[indices,:3] = self.getColors(atoms)
                colList, colIndexList = self.app().indexColors(col)
                mol._colors['cpk'] = numpy.array(colList)
                mol._ag.setData('colorsIndices_cpk', colIndexList)
                self.app().displayCPK.refreshDisplay(mol)
                
            elif gName=='sb':
                col = mol._colors['sb_balls'][mol._ag.getData('colorsIndices_sb_balls').tolist()]
                oldCol = col[indices]
                self.app().pushUndoCmd( self.app().color, (atoms, oldCol, ['sb']), {})
                col[indices,:3] = self.getColors(atoms)
                colList, colIndexList = self.app().indexColors(col)
                mol._colors['sb_balls'] = numpy.array(colList)
                mol._ag.setData('colorsIndices_sb_balls', colIndexList)
                mol._colors['sb_cyl'] = numpy.array(colList)
                mol._ag.setData('colorsIndices_sb_cyl', colIndexList)
                self.app().displaySB.refreshDisplay(mol)
            elif gName=='sticks':
                col = mol._colors['sb_cyl'][mol._ag.getData('colorsIndices_sb_cyl').tolist()]
                oldCol = col[indices]
                self.app().pushUndoCmd( self.app().color, (atoms, oldCol, ['sticks']), {})
                col[indices,:3] = self.getColors(atoms)
                colList, colIndexList = self.app().indexColors(col)
                mol._colors['sb_cyl'] = numpy.array(colList)
                mol._ag.setData('colorsIndices_sb_cyl', colIndexList)
                self.app().displaySB.refreshDisplay(mol)
            elif gName=='atomLabels':
                col = mol._colors['atomLabels'][mol._ag.getData('colorsIndices_atomLabels').tolist()]
                oldCol = col[indices]
                self.app().pushUndoCmd( self.app().color, (atoms, oldCol, ['atomLabels']), {})
                col[indices,:3] = self.getColors(atoms)
                colList, colIndexList = self.app().indexColors(col)
                mol._colors['atomLabels'] = numpy.array(colList)
                mol._ag.setData('colorsIndices_atomLabels', colIndexList)
                self.app().labelAtoms.refreshDisplay(mol)

            elif gName=='residueLabels':
                col = mol._colors['residueLabels'][mol._ag.getData('colorsIndices_residueLabels').tolist()]
                oldCol = col[indices]
                self.app().pushUndoCmd( self.app().color, (atoms, oldCol, ['residueLabels']), {})
                col[indices,:3] = self.getColors(atoms)
                colList, colIndexList = self.app().indexColors(col)
                mol._colors['residueLabels'] = numpy.array(colList)
                mol._ag.setData('colorsIndices_residueLabels', colIndexList)
                self.app().labelResidues.refreshDisplay(mol)

            elif gName.startswith("msms_"):
                #molName = gName[5:].split("_surface")[0]
                molName = mol.name
                col = mol._colors['msms']['%s_surface'%molName][mol._ag.getData('colorsIndices_msms_%s_surface'%molName).tolist()]
                oldCol = col[indices]
                self.app().pushUndoCmd( self.app().color, (atoms, oldCol, [gName]), {})
                col[indices,:3] = self.getColors(atoms)
                colList, colIndexList = self.app().indexColors(col)
                mol._colors['msms']['%s_surface'%molName] = numpy.array(colList)
                mol._ag.setData('colorsIndices_msms_%s_surface'%molName, colIndexList)
                self.app().displayMSMS.refreshDisplay(mol)

            elif gName.startswith("cartoon"):
                col = mol._colors['cartoon'][mol._ag.getData('colorsIndices_cartoon').tolist()]
                oldCol = col[indices]
                self.app().pushUndoCmd( self.app().color, (atoms, oldCol, ['cartoon'] ), {})
                col[indices,:3] = self.getColors(atoms)
                colList, colIndexList = self.app().indexColors(col)
                mol._colors['cartoon'] = numpy.array(colList)
                mol._ag.setData('colorsIndices_cartoon', colIndexList)
                self.app().displayCartoon.refreshDisplay(mol)
                
            elif gName.startswith("coarseMolSurf"):
                print 'NOT YET'
                ## geom = gc.geoms[gName]
                ## oldCol = self.getAtomColors(geom).copy()
                ## self.app().pushUndoCmd( self.app().color, (atoms, oldCol, [gName]), {})
                ## key = (gName+"colors").replace("-", "")
                ## col = atoms.getAtomGroup().getData(key)
                ## col[indices] = allColors #elementColors[atoms.getData('atomicNumber')]
                ## surfAtomInds = gc.boundGeom[gName]['atoms'].getIndices()
                ## surfCol = col[surfAtomInds] # colors of surface atoms
                ## cl_atoms = gc.boundGeom[gName]['cl_atoms']
                ## material = numpy.take(surfCol, cl_atoms, axis=0).astype('f')
                ## gc.geoms[gName].Set(materials=material, redo=1, inheritMaterial=False,
                ##                     #tagModified=False,
                ##                     transparent='implicit')
                ## mol._ag.setData(key, col)
            else:
                print 'ERROR: Color %s not implemented'%gName

class ColorCommand(ColorCommandBase):
    """The ColorCommand provide a command for coloring a user specified set of atoms with user specified set of colors.

    Synopsis:
      None <--- color(atoms, colors, geomsToColor=['all'])

    Required Arguments:

      nodes --- any set of MolKit2.Selection describing molecular components

      colors --- list of rgb or [rgbs] tuple of the same length as atoms

    Optional Arguments:

      geomsToColor --- list of the name of geometries to color default is ['all']
                       meaning all graphical representations of atoms (including
                       hidden ones). Use '*' for all visible representations.

    Package : PmvApp
    Module  : colorCmds
    Class   : ColorCommand
    Command : color
    """

    def doit(self, selection, colors, geomsToColor=['all',]):
        """None <--- color(selection, geomsToColor, colors)
        colors --- is an array of colors for the selection (len(selection) == len(colors))
        """
        assert len(selection)==len(colors)
        self._colors = colors
        ColorCommandBase.doit(self, selection, geomsToColor)
        
    def getColors(self, atoms):
        return self._colors[:, :3]
    

from .pmvPalettes import elementColors, sequentialColors, rainbow256Colors

class ColorByAtomType(ColorCommandBase):
    """The colorByAtomType command allows the user to color the given geometry representing the given nodes using the atomtype coloring scheme where:N :Blue ; C :Gray ; O : Red ; S : Yellow ; H : Cyan; P: Magenta;UNK:green.
    \nPackage : PmvApp
    \nModule  : colorCmds
    \nClass   :  ColorByAtomType
    \nCommand : colorbyAtomType
    \nDescription:\n
    This coloring scheme gives some information on the atomic composition of
    the given nodes.\n
    \nSynopsis:\n
    None <- colorByAtomType(nodes, geomsToColor=['all'])\n
    nodes       : any set of MolKit2.Selection describing molecular components\n
    geomsToColor: list of the name of geometries to color default is 'all'\n
    Keywords: color, atom type\n
    """

    def getColors(self, atoms):
        return elementColors[atoms.getData('atomicNumber')]
    

class ColorByMolecule(ColorCommandBase):
    """The colorByMolecules command allows the user to color the given geometries representing the given selection by molecules. A different color is assigned to each molecule. \n
    Package : PmvApp \n
    Module  : colorCmds \n
    Class   : ColorByMolecule \n
    Command : colorByMolecules \n
    Synopsis:\n
      None <- colorByMolecules(selection, geomsToColor=['all'], carbonsOnly=False)\n
      selection --- any set of MolKit2.Selection describing molecular components\n
      geomsToColor --- list of the name of geometries to color default is ['all']\n
      carbonsOnly --- flag (True, False) 
                        When this flag is set to True only carbon atoms \n
                        will be assigned color.
      Keywords --- color, chain\n
    """

    def getColors(self, sel):
        mol = sel.getAtomGroup().getMolecule()
        return sequentialColors[self.app().Mols.index(mol) % len(sequentialColors)]
        

class ColorByChain(ColorByMolecule):
    """The colorByChain command allows the user to color the given geometries representing the given selection by chain. A different color is assigned to each chain.
    \nPackage : PmvApp
    \nModule  : colorCmds
    \nClass   : ColorByChain
    \nCommand : colorByChains
    \nSynopsis:\n
      None <- colorByChains(selection, geomsToColor=['all'], carbonsOnly=False)\n
      selection --- any set of MolKit2.Selection describing molecular components\n
      geomsToColor --- list of the name of geometries to color default is 'all'\n
      Keywords --- color, chain\n
    """

    def getColors(self, sel):
        mol = sel.getAtomGroup().getMolecule()
        #colors = []
        colors = numpy.array([[1.,1.,1.]]*len(mol._ag), 'f')
        chnum = 0
        for chain in mol._ag.iterChains():
            selChain = sel & chain
            #colors.extend( [sequentialColors[chnum%len(sequentialColors)].tolist()]*len(selChain))
            chinds = selChain.getIndices()
            colors[chinds] = sequentialColors[chnum%len(sequentialColors)]
            chnum += 1
        return colors[sel.getIndices()]

class ColorRainbow(ColorByMolecule):
    """The RainbowColor command colors molecules using a rainbow color map.\n
    Package : PmvApp \n
    Module  : colorCmds \n
    Class   : ColorRainbow \n
    Command : colorRainbow \n
    Synopsis:\n
    None <- colorRainbow(selection, geomsToColor=['all']) \n
    selection --- atom selection \n
    geomsToColor --- list of geometries (names) to be colored
    """  
    def getColors(self, sel):
        mol = sel.getAtomGroup().getMolecule()
        indices = sel.getIndices()
        colors = Map(indices, rainbow256Colors, 0, len(mol._ag))
        return colors
                
class ColorRainbowChain(ColorByMolecule):
    """The RainbowColorChain command colors molecules using a rainbow color map per chain.\n
    Package : PmvApp \n
    Module  : colorCmds \n
    Class   : ColorRainbow \n
    Command : colorRainbow \n
    Synopsis:\n
    None <- colorRainbow(selection, geomsToColor=['all']) \n
    selection --- atom selection \n
    geomsToColor --- list of geometries (names) to be colored
    """  
    def getColors(self, sel):
        mol = sel.getAtomGroup().getMolecule()
        #colors = []
        chnum = 0
        # we build an array of clors spanning all chains
        molCol = numpy.ones( (mol._ag.numAtoms(), 3), 'f')
        for chain in mol._ag.iterChains():
            indices = range(chain.numAtoms())
            mini = min(chain.getIndices())
            colors = Map(indices, rainbow256Colors, 0, len(indices))
            # write rainbown colors for the atoms of this chain into the array
            molCol[chain.getIndices()] = colors
            #selChain = sel & chain
            #colors.extend( colors[selChain.getIndices()].tolist() )
        return molCol[sel.getIndices()]
                
class CustomColor(ColorCommand):
    
    def doit(self, selection, color, geomsToColor=['all',], carbonsOnly=False):
        self._color = color
        ColorCommandBase.doit(self, selection, geomsToColor, carbonsOnly)

    def getColors(self, atoms):
        return [self._color]

class ColorBySecondaryStructure(ColorCommandBase):
    """Command to color selected geometry by secondary structure
    element.
    Package : PmvApp \n
    Module  : colorCmds \n
    Class   : ColorBySecondaryStructure \n
    Command : colorBySS \n
    Synopsis:\n
    None <- colorBySS(selection, geomsToColor=['all'], carbonsOnly=False)\n
    selection --- any set of MolKit2.Selection describing molecular components\n
    geomsToColor --- list of the name of geometries to color default is ['all']\n
    carbonsOnly --- flag (True, False) 
                    When this flag is set to True only carbon atoms \n
                    will be assigned color.
    """
    ssColorMap = {
                  "H": [0.937, 0.0, 0.5], # Helix
                  "G": [0.937, 0.0, 0.5], # Helix
                  "I": [0.937, 0.0, 0.5], # Helix
                  "B": [1.0,   1.0, 0.0], # Sheet
                  "E": [1.0,   1.0, 0.0], # Sheet
                  "C": [1.0,   1.0, 1.0], # Coil 
                  "S": [1.0,   1.0, 1.0], # Coil 
                  "T": [1.0,   1.0, 1.0], # Coil 
                  }
                  
    def getColors(self, sel):
        inds = sel.getIndices()
        ss = sel.getAtomGroup().getSecstrs()
        if ss is None:
            mol = sel.getAtomGroup().getMolecule()
            mol.assignSecondaryStructureWithPross()
            ss = sel.getAtomGroup().getSecstrs()
        coil = self.ssColorMap["C"]
        return numpy.array([self.ssColorMap.get(ch, coil) for ch in ss[inds]])
        

class MoleculesRainbow(ColorCommandBase):

    def getColors(self, atoms):
        numOfMols = len(self.app().Mols)
        mol = atoms.getAtomGroup().getMolecule()
        color = Map([self.app().Mols.index(mol)], rainbow256Colors, mini=0, maxi=numOfMols)
        return color[0]
                

class SetOpacity(MVCommand):
    """Sets opacities to specified geometries
    Package : PmvApp \n
    Module  : colorCmds \n
    Class   : SetOpacity \n
    Command : setOpacity \n
    Description:\n
    Synopsis:\n
      None <--- setOpacity(selection, opacities, geoms=['all']\n
    Required Arguments:\n 
      selection --- any set of MolKit2.Selection describing molecular components\n
    Optional Arguments:\n  
      opacities --- a float or a list of opacity values\n
      geoms--- list of names of geometries, default is ['all']\n
    """

    def makeOpacityEvent(self, *args, **kw):
        event = SetOpacityEvent(args, kw)
        self.app().eventHandler.dispatchEvent(event)
        
    def initializeMolForCmd(self, patoms):
        self.initializedFor[patoms.getMolecule()] = True

    def getAvailableGeoms(self, selections, showUndisplay=0):
        """Method to build a dictionary containing all the geometries
        available in the scene."""
        #build alist of visible geometries
        geomsAvailable = {} # use dict to avoid duplication
        for selection in selections:
            mol = selection.getAtomGroup().getMolecule()
            geomC = mol.geomContainer
            for geomName, geom in geomC.geoms.items():
                if geom.visible and len(geom.vertexSet.vertices):
                    g = geom.parent
                    while g.visible and g.parent:
                        g = g.parent
                    if g.name=='root':
                        #if geomName in ['noBond', 'singleBonds', 'doubleBonds', 'tripleBonds', 'aromaticBonds']:
                        #    geomName = 'lines'
                        geomsAvailable[geomName] = True
        return geomsAvailable.keys()


    def getAtomOpacities(self, geom):
        from DejaVu2.viewerConst import OVERALL, PER_VERTEX
        if geom.materials[GL.GL_FRONT].binding[1]==OVERALL:
            col = geom.materials[GL.GL_FRONT].prop[1][3]
            return numpy.resize( col, (len(geom.vertexSet.vertices), 3) )
        elif geom.materials[GL.GL_FRONT].binding[1]==PER_VERTEX:
            return geom.materials[GL.GL_FRONT].prop[1][:,3]
        else:
            raise
        
    def checkArguments(self, selection, opacities, geoms=['all']):
        """
           selection --- single selection i.e. subset of one molecule
           geomsToColor --- name of geometry objects to color
           opacities --- list of opacities. has to match the len of selection"""
        if len(selection)==0:
            raise ValueError, '%s: No molecular fragment found for %s'%(self.name, str(selection))
        if hasattr(opacities, "__len__"):
            assert len(selection) == len(opacities)
        assert isinstance(geoms, (list, tuple))
        geoms = [x for x in geoms if x not in [' ', '']]
        if 'all' in geoms:
            geoms = self.getAvailableGeoms(selection)
        elif '*' in geoms:
            geoms = self.getAvailableGeoms(selection, showUndisplay=1)
        if "sticksAndBalls" in geoms:
             geoms.remove("sticksAndBalls")
             geoms.extend(["sticks", "balls"])
        if "lines" in geoms:
             geoms.remove("lines")
             geoms.extend(['noBond', 'singleBonds']) #, 'doubleBonds', 'tripleBonds', 'aromaticBonds'])
        mol = selection.getAtomGroup().getMolecule()
        for geomName in geoms:
            assert geomName in mol.geomContainer.geoms.keys()
        return (selection, opacities), {'geoms':geoms}


    def doit(self, selection, opacities, geomsToColor=['all']):
        """None <--- color(selection, geomsToColor, colors)"""

        mol = selection.getAtomGroup().getMolecule()
        for geomName in geomsToColor: 
            geom = mol.geomContainer.geoms[geomName]
            newOpacities = self.getAtomOpacities(geom)
            self.app().pushUndoCmd( self.app().setOpacity, (selection,  newOpacities.copy()), {'geoms':[geomName]})
            indices = selection.getIndices()
            try:
                len(opacities)
            except:
                opacities = numpy.ones((len(newOpacities),), 'f') * opacities
            newOpacities[indices] = opacities[indices]
            geom.Set(opacity = newOpacities, inheritMaterial=False, transparent=True)
        if self.createEvents:
            event = EditGeomsEvent("opacity", [selection, geomsToColor, opacities, self.name[5:11]])
            self.app().eventHandler.dispatchEvent(event)
        self.makeOpacityEvent(selection, opacities=opacities, geoms=geomsToColor)

        
    #def getName(self, args, kw, full=False):
    #    return MVCommand.getName(self, args, kw, full=True)


class ColorFromPalette(ColorCommand):#, MVAtomICOM):
    """The ColorFromPalette class is the base class from which all the color commands using a colorPalette implemented for PMV will derive. \n
    Package : PmvApp \n
    Module  : colorCmds \n
    Class   : ColorFromPalette \n
    Description:\n
    It implements the general functionalities needed to retrieve the colors given a palette and a set of nodes.\n
     Synopsis:\n
     None <- colorFromPalette(nodes, geomsToColor=['all'])\n
     Required Arguments:\n      
     nodes---TreeNodeSet holding the current selection\n
     geomsToColor---list of names of geometries to color, default is ['all']\n
    
    """
    # THE USER SHOULD BE ABLE TO CREATE A COLORPALETTE AND USE IT TO COLOR THE
    # SELECTED NODES WITH IT.
    
    def __init__(self):
        ColorCommand.__init__(self)
        #MVAtomICOM.__init__(self)
        #self.flag = self.flag | self.objArgOnly

    
    def undoCmdBefore(self, nodes, geomsToColor=['all',]):
        # these commands do not require the color argument since colors are
        # gotten from a palette
        # we still can use the ColorCommand.undoCmdBefore method by simply
        # passing None for the color argument
        return ColorCommand.undoCmdBefore(self, nodes, None, geomsToColor)
        

    def doit(self, nodes, geomsToColor=['all',]):
        # these commands do not require the color argument since colors are
        # gotten from a palette
        # we still can use the ColorCommand.undoCmdBefore but first we get
        # the colors. This also insures that the colors are not put inside the
        # log string for these commands

        colors = self.getColors(nodes)
        ColorCommand.doit(self, nodes, colors=colors, geomsToColor=geomsToColor)

            
    def onAddCmdToApp(self):
        # these commands use a color command to undo their effect
        # so we make sure it is loaded and we place its name into
        # undoCmdsString
        from mglutil.util.defaultPalettes import ChooseColor, \
             ChooseColorSortedKeys
        self.palette = ColorPalette(
            'Color palette', ChooseColor, readonly=0, info='Color palette',
            sortedkeys=ChooseColorSortedKeys )
            
        if not self.app().commands.has_key('color'):
            self.app().lazyLoad('colorCmds', commands=['color'], package='PmvApp')
        self.undoCmdsString= self.app().color.name
        
        
    def getColors(self, nodes):
        return self.palette.lookup( nodes )


    def checkArguments(self, nodes, geomsToColor=['all']):
        """
           nodes --- TreeNodeSet holding the current selection \n
           geomsToColor --- list of of geometry names to color,
                           default is ['all'] """
        #print "checkArguments:", self.name
        #assert nodes
        #if isinstance (nodes, str):
        #   self.nodeLogString = "'"+nodes+"'"
        #self.molSet = mols
        #self.atmSet = atoms
        if not nodes:
            raise ValueError, '%s: No molecular fragment found for %s'%(self.name, str(nodes))
        assert isinstance(geomsToColor, (list, tuple))
        geomsToColor = [x for x in geomsToColor if x not in [' ', '']]
        if 'all' in geomsToColor:
            geomsToColor = self.getAvailableGeoms(nodes)
        elif '*' in geomsToColor:
            geomsToColor = self.getAvailableGeoms(nodes, showUndisplay=1)
        assert len(geomsToColor)
        return (nodes,), {'geomsToColor':geomsToColor}


class ColorByDG(ColorFromPalette):
    """The colorByDG command allows the user to color the given geometries representing the given nodes using David Goodsell's coloring
    scheme. \n
    Package : PmvApp \n
    Module  : colorCmds \n
    Class   : ColorByDG \n
    Command : colorAtomsUsingDG\n
    Synopsis:\n
    None <--- colorByDG(nodes, geomsToColor=['all'])\n
    Arguments:\n
    nodes --- any set of MolKit2.Selection describing molecular components\n
    geomsToColor --- list of the name of geometries to color default is 'all'\n
    Keywords --- color, David Goodsell's coloring scheme\n
    """
    
    def __init__(self):
        ColorFromPalette.__init__(self)

        self.DGatomIds=['ASPOD1','ASPOD2','GLUOE1','GLUOE2', 'SERHG',
                        'THRHG1','TYROH','TYRHH',
                        'LYSNZ','LYSHZ1','LYSHZ2','LYSHZ3','ARGNE','ARGNH1','ARGNH2',
                        'ARGHH11','ARGHH12','ARGHH21','ARGHH22','ARGHE','GLNHE21',
                        'GLNHE22','GLNHE2',
                        'ASNHD2','ASNHD21', 'ASNHD22','HISHD1','HISHE2' ,
                        'CYSHG', 'HN']

        
    def onAddCmdToApp(self):
        from PmvApp.pmvPalettes import DavidGoodsell, DavidGoodsellSortedKeys
        c = 'Color palette for coloring using David Goodsells colors'

        self.palette = ColorPaletteFunction(
            'DavidGoodsell', DavidGoodsell, readonly=0, info=c,
            sortedkeys=DavidGoodsellSortedKeys, lookupFunction=self.lookupFunc)
        if not self.app().commands.has_key('color'):
            self.app().lazyLoad('colorCmds', commands=['color'], package='PmvApp')
        self.undoCmdsString = self.app().color.name


    def lookupFunc(self, atom):
        assert isinstance(atom, Atom)
        if atom.name in ['HN']:
            atom.atomId = atom.name
        else:
            atom.atomId=atom.parent.type+atom.name
        if atom.atomId not in self.DGatomIds: 
            atom.atomId=atom.element
        return atom.atomId



class ColorByResidueType(ColorFromPalette):
    """The colorByResidueType command allows the user to color the given geometries representing the given nodes using the Rasmol coloring scheme. \n
    Package : PmvApp \n
    Module  : colorCmds \n
    Class   : ColorByResidueType \n
    Command : colorByResidueType \n
    where:\n
    ASP, GLU    bright red       CYS, MET       yellow\n
    LYS, ARG    blue             SER, THR       orange\n
    PHE, TYR    mid blue         ASN, GLN       cyan\n
    GLY         light grey       LEU, VAL, ILE  green\n
    ALA         dark grey        TRP            pink\n
    HIS         pale blue        PRO            flesh\n
    Synopsis:\n
      None <- colorByResidueType(nodes, geomsToColor=['all'])\n
      nodes --- any set of MolKit2.Selection describing molecular components.\n
      geomsToColor --- list of the name of geometries to color default is 'all'\n
    Keywords --- color, Rasmol, residue type\n
    """
    
    def onAddCmdToApp(self):
        from PmvApp.pmvPalettes import RasmolAmino, RasmolAminoSortedKeys
        c = 'Color palette for Rasmol like residues types'
        self.palette = ColorPalette(
            'RasmolAmino', RasmolAmino, readonly=0, info=c,
            sortedkeys = RasmolAminoSortedKeys, lookupMember='type')
        if not self.app().commands.has_key('color'):
            self.app().lazyLoad('colorCmds', commands=['color'], package='PmvApp')
        self.undoCmdsString = self.app().color.name


    def getColors(self, nodes):
        return self.palette.lookup( nodes.findType(Residue) )
    


class ColorShapely(ColorFromPalette):
    """The colorByShapely command allows the user to color the given geometries representing the given nodes using the Shapely coloring scheme where each residue has a different color. (For more information please refer to the pmv tutorial). \n
    Package : PmvApp \n
    Module  : colorCmds \n
    Class   : ColorShapely \n
    Command : colorResiduesUsingShapely\n
    Synopsis:\n
      None <- colorResiduesUsingShapely(nodes, geomsToColor=['all'])\n
      nodes --- any set of MolKit2.Selection describing molecular components\n
      geomsToColor --- list of the name of geometries to color default is ['all']\n
      Keywords --- color, shapely, residue type\n
    """
    def onAddCmdToApp(self):
        from PmvApp.pmvPalettes import Shapely

        self.palette = ColorPalette('Shapely', Shapely, readonly=0,
            info='Color palette for shapely residues types',
                                    lookupMember='type')
        if not self.app().commands.has_key('color'):
            self.app().lazyLoad('colorCmds', commands=['color'], package='PmvApp',
                                   )
            #colorCmd = self.app().color(loadOnly=True)
        self.undoCmdsString = self.app().color.name


    def getColors(self, nodes):
        if not nodes:
            raise ValueError, "%s: Missing molecular fragment to color" % self.name
        return self.palette.lookup( nodes.findType(Residue) )



##
## FIXME ramp should be a colorMap object with editor and stuff
## it should be possible to pass it as an argument. if none is given a default
## a colorMap object with a default RGBRamp would be used
## This is also true for Palettes.
##
from DejaVu2.colorTool import Map

class ColorFromRamp(ColorFromPalette):
    """The colorFromRamp class implements the functionality to color the given geometries representing the given nodes using a colorMap created from the Ramp.
    \nPackage : PmvApp
    \nModule  : colorCmds
    \nClass   : ColorFromRamp
    \nSynopsis:\n
      None <- colorFromRamp(nodes, geomsToColor=['all'])\n
      nodes --- any set of MolKit2.Selection describing molecular components\n
      geomsToColor --- list of the name of geometries to color default is ['all']\n
      Keywords --- color, ramp\n
    """
 
    def __init__(self):
        ColorFromPalette.__init__(self)
        #self.flag = self.flag | self.objArgOnly
        from DejaVu2.colorTool import RGBRamp, Map
        self.ramp = RGBRamp()

    


class ColorByInstance(ColorFromPalette):
    """Command to color the current selection by instance using a Rainbow palette.
    \nPackage : PmvApp
    \nModule  : colorCmds
    \nClass   : ColorByInstance
    \nCommand : colorByInstance
    \nSynopsis:\n
      None <- colorByInstance(nodes, geomsToColor=['all'])\n
      nodes --- any set of MolKit2.Selection describing molecular components\n
      geomsToColor --- list of the name of geometries to color default is ['all']\n
    """

    def onAddCmdToApp(self):
        if not self.app().commands.has_key('color'):
            self.app().lazyLoad('colorCmds', commands=['color'], package='PmvApp')
        self.undoCmdsString= self.app().color.name
        from mglutil.util.defaultPalettes import Rainbow, RainbowSortedKey
        c = 'Color palette molecule number'
        c = ""
        self.palette = ColorPaletteFunction(
            'Rainbow', Rainbow, readonly=0, info=c,
            lookupFunction = lambda x, length=len(RainbowSortedKey): \
            x%length, sortedkeys=RainbowSortedKey)


    def onAddObjectToViewer(self, obj):
        self.objectState[obj] = {'onAddObjectCalled':True}
        obj.number = self.app().Mols.index(obj)


    def doit(self, nodes, geomsToColor=['all']):
        # nodes is an atom set
        for m in self.molSet:
            geomc = m.geomContainer
            #moreGeoms = []
            #if 'lines' in geomsToColor:
            #    moreGeoms = ['bonded', 'bondorder','nobnds']
            #for g in geomsToColor+moreGeoms:
            for g in geomsToColor:
                try:
                    ge = geomc.geoms[g]
                    colors = self.palette.lookup(range(len(ge.instanceMatricesFortran)))
                    ge.Set(materials=colors, inheritMaterial=0)#, tagModified=False)
                    ge.SetForChildren(inheritMaterial=True, recursive=True)
                except:
                    if self.app().trapExceptions is False:
                        exc_info = sys.exc_info()
                        raise exc_info[1], None, exc_info[2]
                    else:
                        msg = '%s: Error while coloring %s for molecule %s'%(self.name, g, m.name,)
                        self.app().errorMsg(sys.exc_info(), msg, obj=self.atmSet)
                #ColorCommand.doit(self, m, colors, geomsToColor)
            self.app()._executionReport.addSuccess('%s for molecule %s success'% (self.name, m.name))


class ColorByProperties(ColorCommand):
    """
    Command to color the current selection according to the integer
    or float properties, or by defining a function.
    \nPackage : PmvApp
    \nModule  : colorCmds
    \nClass   : ColorByProperties
    \nCommand : colorByProperty
    \nSynopsis:\n
      None <- colorByProperty(nodes, geomsToColor, property,colormap='rgb256')\n
      nodes --- any set of MolKit2.Selection describing molecular components\n
      geomsToColor --- list of the name of geometries to color default is 'all'\n
      property ---  property name of type integer or float or property defined by a function returning a list of float or int.\n
      colormap ---  either a string representing a colormap or a DejaVu2.ColorMap instance.\n
    """

    levelOrder = {'Atom':0 , 
                  'Residue':1,
                  'Chain':2,
                  'Molecule':3 }

    def __init__(self):
        ColorCommand.__init__(self)
        #self.flag = self.flag & 0
        self.level = Atom  # class at this level (i.e. Molecule, Residue, Atom)
        

    def onAddCmdToApp(self):
        # the following commands use a color command to undo their effect
        # so we make sure it is loaded and we place its name into
        # undoCmdsString
        if not self.app().commands.has_key('color'):
            self.app().lazyLoad('colorCmds', commands=['color'], package='PmvApp')
        if not self.app().commands.has_key('saveSet'):
            self.app().lazyLoad('selectionCmds', commands=['saveSet'], package='PmvApp')
            #colorCmd = self.app().color(loadOnly=True)
        self.app().loadColormap()
        self.undoCmdsString = self.app().color.name
        self.molDict = {'Molecule':Molecule,
                        'Atom':Atom, 'Residue':Residue, 'Chain':Chain}

        self.leveloption={}
        for name in ['Atom', 'Residue', 'Molecule', 'Chain']:
            col = self.app().levelColors[name]
            bg = ToHEX((col[0]/1.5,col[1]/1.5,col[2]/1.5))
            ag = ToHEX(col)
            self.leveloption[name]={'bg':bg,'activebackground':ag,
                                    'borderwidth':3}
        self.propValues = None
        self.level = "Molecule"
        self.propertyLevel = self.level


    def getPropValues(self, nodes, prop, propertyLevel=None):
        #print "getPropValues", self, nodes, prop, propertyLevel
        try:
            if propertyLevel is not None:
                lNodesInLevel = nodes.findType(self.molDict[propertyLevel])
                self.propValues = getattr(lNodesInLevel, prop)
            else:
                self.propValues = getattr(nodes, prop)
        except:
            from Pmv import formatName
            msg= "%s.getPropValues: nodes(%s) do not have property %s" % (self.name, formatName(nodes.buildRepr(), 60), prop)
            raise RuntimeError(msg)


    def undoCmdBefore(self, nodes, geomsToColor,  property,
                        propertyLevel=None, colormap='rgb256',
                        mini=None, maxi=None, carbonsOnly=False):
        # nodes is an atom set
        if carbonsOnly:
            nodes= AtomSet([x for x in nodes if x.element=='C'])
        return ColorCommand.undoCmdBefore(self, nodes, colors=None, geomsToColor=geomsToColor)


    def doit(self, nodes, geomsToColor, property, propertyLevel=None,
             colormap='rgb256', mini=None, maxi=None, carbonsOnly=False):

        #print 'VVVV', len(nodes), property, propertyLevel
        #nodes = nodes.findType(self.app().selectionLevel)
        #print 'CCCCCC', len(nodes)
        if "sticksAndBalls" in geomsToColor:
            geomsToColor.remove("sticksAndBalls")
            geomsToColor.extend(["sticks", "balls"])
        if isinstance (colormap, str):
            colormap = self.app().colorMaps[colormap]
        # Get the list of values corresponding the the chosen property
        # if not already done ?
        self.getPropValues(nodes, property, propertyLevel)
        # build the color ramp.
        selectioncol = colormap.Map(self.propValues, mini=mini, maxi=maxi)
        # Call the colorProp method
        self.colorProp(nodes, geomsToColor, selectioncol,
                       propertyLevel, carbonsOnly )
        
        ## insideIntervalnodes = []
        ## for i in range(len(nodes)):
        ##     if self.propValues[i] >= mini and self.propValues[i] <= maxi:
        ##         insideIntervalnodes.append(nodes[i])
        ## if nodes[0].__class__.__name__.endswith('Atom'):
        ##     lSet = AtomSet(insideIntervalnodes)
        ## elif nodes[0].__class__.__name__.endswith('Residue'):
        ##     lSet = ResidueSet(insideIntervalnodes)
        ## elif nodes[0].__class__.__name__.endswith('Chain'):
        ##     lSet = ChainSet(insideIntervalnodes)
        ## elif nodes[0].__class__.__name__.endswith('Molecule'):
        ##     lSet = MoleculeSet(insideIntervalnodes)
        ## elif nodes[0].__class__.__name__.endswith('Protein'):
        ##     lSet = ProteinSet(insideIntervalnodes)

##         if 'ColorByProperties' in self.app().sets.keys():
##             self.app().sets.pop('ColorByProperties')
##         self.app().saveSet(lSet, 'ColorByProperties', log=False,
##                 comments="""Last set created by colorByProperties,
## contains only nodes with chosen property between mini and maxi.
## This set is ovewritten each time ColorByProperties is called.
## """)
        #geomEditEventss
        
        event = EditGeomsEvent("color", [nodes,[geomsToColor, selectioncol, self.name[5:11]]])
        self.app().eventHandler.dispatchEvent(event)
        #return lSet


    def colorProp(self, nodes, geomsToColor, selectioncol, propertyLevel='Atom', carbonsOnly=False):
        # nodes is an atom set
        if (propertyLevel is not None) and \
           self.levelOrder[propertyLevel] < self.levelOrder[self.level]:
            nodes = nodes.findType(self.molDict[propertyLevel])     

        # loop over the node and assign the right color to the atoms.
        deltaOpac = 0.0
        for gName in geomsToColor:
            for i, n in enumerate(nodes):
                if not carbonsOnly or n.chemElem == "C":
                    #if n.colors.has_key(gName):
                    n.colors[gName] = tuple(selectioncol[i][:3])
                    if n.opacities.has_key(gName):
                        newOpac = selectioncol[i][3]
                        oldOpac = n.opacities[gName]
                        deltaOpac = deltaOpac + (oldOpac-newOpac)
                        n.opacities[gName] = newOpac
        for mol in self.molSet:
            try:
                updatedGeomsToColor = []
                for gName in geomsToColor:
                    if not mol.geomContainer.geoms.has_key(gName): continue
                    geom = mol.geomContainer.geoms[gName]
                    if geom.children != []:
                        # get geom Name:
                        childrenNames = [x.name for x in geom.children]
                        updatedGeomsToColor = updatedGeomsToColor + childrenNames
                        for childGeom in geom.children:
                            childGeom.Set(inheritMaterial=0, redo=0, tagModified=False)
                    else:
                        updatedGeomsToColor.append(gName)
                        geom.Set(inheritMaterial=0, redo=0, tagModified=False)

                updateOpac = (deltaOpac!=0.0)
                mol.geomContainer.updateColors(updatedGeomsToColor,
                                               updateOpacity=updateOpac)
                self.app()._executionReport.addSuccess('%s: colored molecule %s successfully'% (self.name, mol.name))
            except:
                if self.app().trapExceptions is False:
                    exc_info = sys.exc_info()
                    raise exc_info[1], None, exc_info[2]
                else:
                    msg = '%s: Error in colorProp() for molecule %s'%(self.name, mol.name)
                    self.app().errorMsg(sys.exc_info(), msg, obj=self.atmSet)
            

    def checkArguments(self, nodes, geomsToColor, property,
                       propertyLevel='Atom', colormap='rgb256',
                       mini=None, maxi=None, carbonsOnly=False):
        """None <- colorByProperty(nodes, geomsToColor, property,colormap='rgb256', **kw)
        \nnode --- TreeNodeSet holding the current selection
        \ngeomsToColor --- the list of the name geometries to be colored
        \nproperty ---   property name of type integer or float or property defined by a function returning a list of float or int.
        \ncolormap--- either a string representing a colormap or a DejaVu2.ColorMap instance.
        """
        assert nodes
        if isinstance (nodes, str):
            self.nodeLogString = "'"+nodes+"'"
        mols, atms = self.getNodes(nodes)
        self.molSet = mols
        self.atmSet = atms
        assert isinstance(geomsToColor, (list, tuple))
        geomsToColor = [x for x in  geomsToColor if x not in [' ', '']]
        if 'all' in geomsToColor:
            geomsToColor = self.getAvailableGeoms(mols)
        if '*' in geomsToColor:
            geomsToColor = self.getAvailableGeoms(mols, showUndisplay=1)
        assert len(geomsToColor)
        assert propertyLevel in ['Atom', 'Residue', 'Molecule', 'Chain']
        assert isinstance(property, str)
        if mini is not None:
            assert isinstance (mini, (int, float))
        if maxi is not None:
            assert isinstance (maxi, (int, float))
            assert maxi >= mini
        if isinstance (colormap, str):
            assert self.app().colorMaps.has_key(colormap)
        else:
            assert isinstance(colormap, ColorMap)

        kw= {}
        kw['colormap'] = colormap
        kw['mini'] = mini
        kw['maxi'] = maxi
        kw['propertyLevel'] = propertyLevel
        return (atms, geomsToColor, property), kw



class ColorByExpression(ColorByProperties):
    """The colorByExpression command allows the user to color the given geometries representing the given nodes evaluated by  python function or lambda function. \n
    Package : PmvApp \n
    Module  : colorCmds \n
    Class   : ColorByExpression \n
    Command : colorByExpression \n
    Synopsis:\n
     None <- colorByExpression(nodes, geomsToColor, function, colormap='rgb256', min=None, max=None) \n
     nodes --- TreeNodeSet holding the current selection \n
     geomsToColor --- the list of the name geometries to be colored \n
     function --- python function or lambda function that will be evaluated with the given nodes \n
     colormap ---  can either be a string which is the name of a loaded colormap or a DejaVu2.colorMap.ColorMap instance.
  """  
    
    # comments for map function definition window
    mapLabel = """Define a function to be applied
on each node of the current selection: """

    # code example for map function definition window
    mapText = '\
#This is a demo function for this widget.\n\
def foo(object):\n\
\tif hasattr(object, "number"):\n\
\t\treturn object._uniqIndex\n\
\telse:\n\
\t\treturn 0\n\
\n'

    # comments for function to operate on selection
    funcLabel = """Define a function to be applied to
the current selection:"""
    
    # code example for function to operate on selection
    funcText ='#def foo(selection):\n#\tvalues = []\n#\t#loop on the current selection\n#\tfor i in xrange(len(selection)):\n#\t\t#build a list of values to color the current selection.\n#\t\tif selection[i].number > 20:\n#\t\t\tvalues.append(selection[i].number*2)\n#\t\telse:\n#\t\t\tvalues.append(selection[i].number)\n#\t# this list of values is then returned.\n#\treturn values\n'

    def __init__(self):
        ColorCommand.__init__(self)
        #self.flag = self.flag & 0

    def onAddCmdToApp(self):
        ColorByProperties.onAddCmdToApp(self)
        self.evalFlag = 0
        self.propValues=None
        

        
    def undoCmdBefore(self, nodes, geomsToColor,
                        function='lambda x: x._uniqIndex',
                        colormap='rgb256'):
        return ColorCommand.undoCmdBefore(self, nodes, colors=None, geomsToColor=geomsToColor)


    def checkArguments(self, nodes,  geomsToColor=['all'],
                 function='lambda x: x._uniqIndex',
                 colormap='rgb256'):
        """
        nodes --- TreeNodeSet holding the current selection \n
        geomsToColor --- the list of the name geometries to be colored \n
        function --- python function or lambda function that will be evaluated with the given nodes \n
        colormap ---  can either be a string which is the name of a loaded colormap or a DejaVu2.colorMap.ColorMap instance.
        """
        assert nodes
        if isinstance (nodes, str):
            self.nodeLogString = "'"+nodes+"'"
        mols, atms = self.getNodes(nodes)
        self.molSet = mols
        self.atmSet = atms
        assert isinstance(geomsToColor, (list, tuple))
        geomsToColor = [x for x in  geomsToColor if x not in [' ', '']]
        if 'all' in geomsToColor:
            geomsToColor = self.getAvailableGeoms(mols)
        if '*' in geomsToColor:
            geomsToColor = self.getAvailableGeoms(mols, showUndisplay=1)
        assert len(geomsToColor)
        if isinstance (colormap, str):
            assert self.app().colorMaps.has_key(colormap)
        else:
            assert isinstance(colormap, ColorMap)
        assert isinstance (function, str)
        kw = {}
        kw['function'] = function
        kw['colormap'] = colormap
        return (atms, geomsToColor), kw
    

    def getPropValues(self, nodes, function):
        func = evalString(function)
        self.propValues = func(nodes)

        
    def doit(self, nodes,  geomsToColor,
             function='lambda x: x._uniqIndex', colormap='rgb256'):
        #nodes is an atom set
        if "sticksAndBalls" in geomsToColor:
            geomsToColor.remove("sticksAndBalls")
            geomsToColor.extend(["sticks", "balls"])
        if isinstance(colormap, str):
            colormap = self.app().colorMaps[colormap]
        # get the values
        self.getPropValues(nodes, function)

        # Get the color corresponding the values
        selectioncol = colormap.Map(self.propValues, colormap.mini,
                                    colormap.maxi)
        
        self.colorProp(nodes, geomsToColor, selectioncol)
        #geomEditEventss
        event = EditGeomsEvent("color", [nodes,[geomsToColor, selectioncol, self.name[5:11]]])
        self.app().eventHandler.dispatchEvent(event)



## class RainbowColor(ColorFromPalette):
##     """The RainbowColor command colors molecules using a rainbow color map.\n
##     Package : PmvApp \n
##     Module  : colorCmds \n
##     Class   : RainbowColor \n
##     Command : colorRainbow \n
##     Synopsis:\n
##     None <- colorRainbow(nodes, geomsToColor=['all']) \n
##     nodes --- TreeNodeSet holding the current selection \n
##     geomsToColor --- list of geometries (names) to be colored
##     """  
    

##     def doit(self, nodes, geomsToColor=['all']):
##         # nodes is an atom set
##         nodes._uniqNumber = range(1, len(nodes)+1)  
##         self.app().colorByProperty(nodes, geomsToColor, '_uniqNumber',
##                                 propertyLevel='Atom', colormap='rgb256',
##                                 mini=1., maxi=len(nodes), 
##                                 setupUndo=False)
        


## class RainbowColorByChain(ColorFromPalette):
##     """The RainbowColor command colors molecule's chains using a rainbow color map.\n
##     Package : PmvApp \n
##     Module  : colorCmds \n
##     Class   : RainbowColor \n
##     Command : colorRainbowByChain\n
##     Synopsis:\n
##      None <- colorRainbowByChain(nodes, geomsToColor=['all']) \n
##      nodes --- TreeNodeSet holding the current selection \n
##      geomsToColor --- list of geometries (names) to be colored
##     """  
##     def doit(self, nodes, geomsToColor=['all']):
##         # nodes is an atomset
##         #from time import time
##         #t1 = time()
##         chains = nodes.findType(Chain)
##         if nodes == nodes.top.uniq().allAtoms:
##             chains = chains.uniq()
##             for chain in chains:
##                 atoms = chain.findType(Atom)
##                 atoms._uniqNumber = range(1, len(atoms)+1)
##                 self.app().colorByProperty(atoms, geomsToColor, '_uniqNumber',
##                              propertyLevel='Atom', colormap='rgb256',
##                              mini=1., maxi=len(atoms), 
##                              setupUndo=False)
##         else:
##             dd = {}.fromkeys(chains.uniq(), AtomSet())
##             for i, ch in enumerate(chains): dd[ch].append(nodes[i])
##             for chain, atoms in dd.items():
##                 atoms._uniqNumber = range(1, len(atoms)+1)
##                 self.app().colorByProperty(atoms, geomsToColor, '_uniqNumber',
##                                     propertyLevel='Atom', colormap='rgb256',
##                                     mini=1., maxi=len(atoms), 
##                                     setupUndo=False)
##         #print "done RainbowColorByChain in :", time() -t1



class ColorByLineGeometryColor(ColorFromPalette):
    """ Colors selected geometries by the Lines color scheme.\n
    Package : PmvApp \n
    Module  : colorCmds \n
    Class   : ColorByLineGeometryColor \n
    Command : colorByLinesColor\n
    Synopsis:\n
     None <- colorByLinesColor(nodes, geomsToColor=['all']) \n
     nodes --- TreeNodeSet holding the current selection \n
     geomsToColor --- list of geometries (names) to be colored"""
            
    def doit (self, nodes, geomsToColor=['all']):
        # nodes is an atom set
        if 'lines' in geomsToColor:
            geomsToColor.remove('lines')
            if not len(geomsToColor): return
        if "sticksAndBalls" in geomsToColor:
            geomsToColor.remove("sticksAndBalls")
            geomsToColor.extend(["sticks", "balls"])
        for geom in geomsToColor:
            for a in nodes:
                if a.colors.has_key(geom):
                    del a.colors[geom]
            if geom == 'sticks' or geom == 'balls':
                self.app().displaySB(nodes)
            elif geom == 'cpk':
                self.app().displayCPK(nodes)
            elif geom == "MSMS-MOL":
                self.app().computeMSMS(nodes)
            #CAballs CAsticks
            elif geom == "CAballs" or geom == "CAsticks":
                self.app().displayBackboneTrace(nodes)
            elif geom == "secondarystructure":
                self.app().displayExtrudedSS(nodes)
            elif geom.find("Labels")>0:
                ColorFromPalette.doit(self, nodes, [geom,])
                

    def getColors(self, nodes):
        return nodes.colors.get('lines', [(1.,1.,1.),])


    
commandClassFromName = {
    'color' : [ColorCommand,  None],
    'colorByAtomType' : [ColorByAtomType, None],
    'colorByMolecules' : [ColorByMolecule, None],
    'colorByChains' : [ColorByChain, None],
    'colorRainbow' : [ColorRainbow, None],
    'colorRainbowChain' : [ColorRainbowChain, None],
    'customColor': [CustomColor, None],
    'setOpacity': [SetOpacity, None],
    'moleculesRainbow': [MoleculesRainbow, None],
    'colorBySS': [ColorBySecondaryStructure, None]
    #'colorByResidueType' : [ColorByResidueType, None],
    #'colorAtomsUsingDG' : [ColorByDG, None],
    #'colorResiduesUsingShapely' : [ColorShapely, None],
    #'colorByInstance' : [ColorByInstance, None],
    #'colorByProperty' : [ColorByProperties, None],
    #'colorByExpression' : [ColorByExpression, None],
    #'colorByLinesColor' :[ColorByLineGeometryColor, None]
    }


def initModule(viewer):
    for cmdName, values in commandClassFromName.items():
        cmdClass, guiInstance = values
        viewer.addCommand(cmdClass(), cmdName, guiInstance)
