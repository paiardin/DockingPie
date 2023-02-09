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
# Author: Sophie I. COON, Michel F. SANNER, Anna Omelchenko
#
# Copyright: M. Sanner TSRI 2014
#
#############################################################################
#
# $Header: /mnt/raid/services/cvs/PmvApp/secondaryStructureCmds.py,v 1.6.4.1 2017/07/13 20:55:28 annao Exp $
#
# $Id: secondaryStructureCmds.py,v 1.6.4.1 2017/07/13 20:55:28 annao Exp $
#



import numpy

from opengltk.OpenGL import GL
from DejaVu2.IndexedPolygons import IndexedPolygons
from DejaVu2.Shapes import Shape2D, Triangle2D, Circle2D, Rectangle2D,\
     Square2D, Ellipse2D

from MolKit2.tree import TreeNodeSet
from MolKit2.molecule import Atom, AtomSet, Molecule, MoleculeSet
from MolKit2.protein import Protein, Residue, Chain, ResidueSet, ProteinSet
from MolKit2.protein import SecondaryStructure, SecondaryStructureSet, \
     Helix, Strand, Turn, Coil
from NucleicBases import Add_Nucleic_Bases

from PmvApp.Pmv import MVCommand
from PmvApp.extruder import Sheet2D, ExtrudeSSElt, ExtrudeNA
from PmvApp.displayCmds import DisplayCommand
from PmvApp.colorCmds import ColorFromPalette
from PmvApp.colorPalette import ColorPalette
from PmvApp.Pmv import AfterDeleteAtomsEvent
from PmvApp.Pmv import DeleteGeomsEvent, AddGeomsEvent, EditGeomsEvent

from AppFramework.App import RemoveGeometryEvent, AddGeometryEvent



class ComputeSecondaryStructureCommand(MVCommand):
    """The computeSecondaryStructure command gets the information on the secondary structure of each molecule from the current selection. This information is then used to create objects describing the various secondary structure elements. \n
    Package : PmvApp \n
    Module  : secondaryStructureCmds \n
    Class   : ComputeSecondaryStructureCommand \n
    Command name : computeSecondaryStructure \n
    Description:\n 
    The SS element object belonging to a chain are then grouped into sets.\n
    A new level is added in the 4 level hierarchy...\n
    The information is taken from the file when available or using stride when\n
    available. This command can be used as an interactive command. \n
    Synopsis:\n
        None <--- ComputeSS(nodes, molMode={}) \n
    Required Arguments:\n   
        nodes --- any set for MolKit2.Selection describing molecular components \n
    Optional Arguments:\n
        molmode --- dictionary key: molecule name, value : 'From File' or 'From Stride'. \n
    Required Packages:\n
      MolKit, DejaVu2, mglutil, OpenGL, ViewerFramework \n
    Known bugs:\n
      None \n
    Examples:\n
      mol = mv.Mols[0] \n
      mv.computeSecondaryStructure(mol)
    """
    

    def __init__(self):
        MVCommand.__init__(self)
        #self.flag = self.flag | self.objArgOnly


    def onRemoveObjectFromViewer(self, obj):
        """
        Method to delete sets created by this command.
        This is done to prevent circular references and memory leaks.
        """
        if self.app().undoableDelete__: return
        if not hasattr(obj, 'chains'): return
        for c in obj.chains:
            if not hasattr(c, 'secondarystructureset') :
                continue
            if c.secondarystructureset is None:
                delattr(c.secondarystructureset)
            else:
                # Cleaning up all the circular references created in this
                # command
                while len(c.secondarystructureset)!=0:
                    if hasattr(c.secondarystructureset[0], 'exElt'):
                        delattr(c.secondarystructureset[0], 'exElt')
                    delattr(c.secondarystructureset[0], 'children')
                    delattr(c.secondarystructureset[0], 'residues')
                    delattr(c.secondarystructureset[0], 'start')
                    delattr(c.secondarystructureset[0], 'end')
                    delattr(c.secondarystructureset[0], 'parent')
                    delattr(c.secondarystructureset[0], 'chain')
                    delattr(c.secondarystructureset[0], 'top')
                    del(c.secondarystructureset[0])


    def onAddCmdToApp(self):
        # Try to import stride and set the flag to 1 if succeeds and to 0 if
        # does not
        try:
            import stride
            self.haveStride = 1
        except:
            self.haveStride = 0

        # Load the dependent commands if not already loaded in the
        # application
        if not self.app().commands.has_key('saveSet'):
            self.app().lazyLoad('selectionCmds', commands=['saveSet'],
                             package='PmvApp')

    def checkArguments(self, nodes, molModes=None):
        """None <--- computeSecondaryStructure(nodes, molModes = None) \n
        nodes --- TreeNodeSet holding the current selection. \n
        moldMode --- dictionary {name of the protein: 'From File', 'From PROSS', or 'From Stride'},\n
                  'From File' to get the information from the file,\n
                  'From Pross' to use the pross code to assing SS\n
                  'From Stride' to use stride (requires stride to be installed).
        """
        if isinstance(nodes, str):
            self.nodeLogString = "'"+nodes+"'"
        nodes = self.app().expandNodes(nodes)
        assert len(nodes)
        if molModes:
            assert isinstance(molModes, dict)
            for mode in molModes.values():
                if mode:
                    assert mode in ['From File', 'From Pross', 'From Stride']
        kw = {}
        kw['molModes'] = molModes
        return (nodes,), kw
        

    def doit(self, nodes, molModes=None):
        molecules, nodeSets = self.app().getNodesByMolecule(nodes)
        # Loop over the molecules
        for mol in molecules:
            try:
                for c in mol.chains:
                    c.ribbonType()

                # Determine what mode to use to get the information
                if not molModes:
                    if mol.hasSS:
                        # Information has already been computed then
                        # continue
                        continue
                    else:
                        # Find out the possibilities and set the mode to
                        # one of them.
                        if mol.parser:                    
                            if mol.parser.hasSsDataInFile():
                                # File contains information the mode will be 'From
                                # File'
                                mode = 'From File'
                            else:
                                mode = 'From Pross'
                        elif self.haveStride:
                            # Stride is available on the platform
                            # but no info in the file then stride will be used
                            mode='From Stride'
                        else:
                            mode='From Pross'
                else:
                    #print molModes
                    # a mode to get the information has been specified
                    # for the given molecules
                    if molModes and not molModes.has_key(mol.name):
                        # if the mode has not been specified for a molecule
                        # print a message and continue
                        raise RuntimeError, '%s: No mode has been specified for %s'% (self.name, mol.name)
                    else:
                        # Set the mode to the given value.
                        mode = molModes[mol.name]

                # if this mode has already been used pass.
                if mode in mol.hasSS: continue

                # if secondarystructure have been computed once using another
                # mode need to clean up first
                elif mol.hasSS != []:
                    self.clean(mol)

                # Then compute using the new given mode.
                #if mode is None:
                #    mol.secondaryStructureFromFile()

                if mode == 'From File':
                    # If both modes available try file first if fails use stride
                    # instead.
                    #if not mol.parser.hasSsDataInFile():
                    #    # GIVE FEEDBACK TO THE USER !!!!
                    #    self.warningMsg("WARNING: "+mol.name + \
                    #                    ".pdb does not contain Secondary \
                    #                    Structure information.")
                    #    continue
                    #else:
                        mol.secondaryStructureFromFile()
                        #self.savesets(mol)

                elif mode == 'From Stride':
                    if not self.haveStride:
                        raise RuntimeError , "%s: Stride is not available on \
                        this computer to compute the secondary structure of " %(self.name, mol.name+".pdb")
                    else:
                        mol.secondaryStructureFromStride()
                        #self.savesets(mol)

                elif mode == 'From Pross':
                    mol.secondaryStructureFromPross()
                    #self.savesets(mol)

                #self.savesets(mol)
                self.app()._executionReport.addSuccess('Computed secondary structure for molecule %s successfully'% mol.name, obj=mol)
            except:
                msg = 'Error while computing secondary structure for molecule %s'%mol.name
                self.app().errorMsg(sys.exc_info(), msg, obj=mol)


    def savesets(self, mol):
        for c in mol.chains:
            if not hasattr(c, 'secondarystructureset'): continue
            for ss in c.secondarystructureset:
                name = "%s%s"%(ss.name, ss.chain.id)
                if ss.residues:  #Bugfix for  #1033
## MS calling the command slows this down a lot
## but the sets are needed to extrude, so we add sets 'by hand'
## side effect: no vision nodes for these sets
#                    self.app().saveSet(ss.residues, mol.name+':'+name[-1]
#                                    +':'+name[:-1],
#                                    '%s-%s' %(ss.residues[0].name,
#                                              ss.residues[-1].name),
#                                    )
                    name = mol.name+':'+name[-1] +':'+name[:-1]
                    ss.residues.comments = '%s-%s'%(ss.residues[0].name,
                                                    ss.residues[-1].name)
                    self.app().sets.add(name, ss.residues)

    
    def clean(self, mol):
        """
        This method is called when getting the secondary structure information
        using stride after having from file and vice versa. It is used to
        delete all the secondary structure objects and attributes created
        previously."""
        # Compute secondary structure creates the following:
        # - Secondary structure elements
        # - Save the secondary structure elements residues as a set
        # - new mol attribute hasSS which is a list
        #from PmvApp.selectionCommands import sets__
        molName = mol.name
        mol.hasSS = []
        # Here need to check if the secondary structure element have
        # been extruded or not. If yes need to clean up that as well.
        if hasattr(mol, '_ExtrudeSecondaryStructureCommand__hasSSGeom')\
        and mol._ExtrudeSecondaryStructureCommand__hasSSGeom:
            self.app().extrudeSecondaryStructure.clean(mol)

        for chain in mol.chains:
            # delete the secondarystructureset
            if not hasattr(chain, 'secondarystructureset'): continue
            for ss in chain.secondarystructureset:
                name = "%s%s"%(ss.name, ss.chain.id)
                setName = mol.name+':'+name[-1]+':'+name[:-1]
                if self.app().sets.has_key(setName):
                    del self.app().sets[setName]
                #del sets__[mol.name+':'+name[-1]+':'+name[:-1]]
            delattr(chain, 'secondarystructureset')
            
            # delete the secondarystructure attribute of the residues when
            # existing.
            resTest = [delattr(x, 'secondarystructure') for x in chain.residues if hasattr(x, 'secondarystructure')]
            
        # call the onRemoveObjectFromViewer for the mol.
        self.app().undoableDelete__ = False
        self.onRemoveObjectFromViewer(mol)
        del self.app().undoableDelete__

        # Also need to clean up the sheet2D information.
        for c in mol.chains:
            if hasattr(c, 'sheet2D') and c.sheet2D.has_key('ssSheet2D'):
                del c.sheet2D['ssSheet2D']

        


class ExtrudeSecondaryStructureCommand(MVCommand):
    """The ExtrudeCommand allows the user to represent the secondary structure elements by extruding 2D geometries along a 3D path.To execute this command use the entry 'extrude Secondary Structure' under the 'Compute' menu in the menu bar. 
   The panel that appears lets the user choose the 2D shapes for the extrusion. The entry 'default' in the listChooser lets do a traditional ribbon representation.nbchords represents the number of points in the path3D corresponding to one residue. The higher this parameter is the smoother the extruded geometries will look.gapBeg allows the user to introduce a gap of gapBeg points the extruded geometrie before each residue.gapEnd allows the user to introduce a gap of gapEnd points the extruded geometrie after each residue.The value of this two parameters depend on the value of the nbchords parameter and on each other's value.Once you clique OK on this panel another panel appears to let the user caracterize the chosen 2D geometry.Once the user chose all the parameters an ExtrudeSSElt object is created for each secondary structure element. The geometries associated to each secondary structure element are then updated with the new vertices and faces.Finally the displaySSCommand is executed.This command has the objArgsOnly flag. \n
   Package : PmvApp \n
   Module  : secondaryStructureCommands \n
   Class   : ExtrudeSecondaryStructureCommand \n
   Command name : extrudeSecondaryStructure \n
   Synopsis:\n
            None <--- extrudeSecondaryStructure(nodes, shape1=None, shape2=None,frontcap=1, endcap=True, arrow=True, nbchords=8, gapBeg=False,gapEnd=False, larrow=2, display=True) \n
   Required Arguments:\n    
        nodes ---  TreeNodeSet holding the current selection(mv.getSelection()) \n
   Optional Arguments:\n    
        shape1 &
        shape2 --- DejaVu2.Shapes.Shape2D objects. shape1 will be used to
                  represent the helix and strand, shape2 to represent coils and
                  turns. \n
        frontcap &
        endcap   --- Boolean flag when set to True a  cap will be added to the
                   geom either at the front or at the end \n
        arrow  --- Boolean flag when set to True an arow will be added to the
                   geometry representing the strand.\n
        nbchords --- Nb of points per residues in the smooth array \n
        gapBeg&  
        gapEnd  --- defines gap at the beginning or the end of each residue. \n
        larrow  --- lenth of the arrow if arrow boolean flag set to 1 \n
        display  --- Boolean flag when set to True the displaySecondaryStructure
                   is called automatically
    """

    def __init__(self):
        MVCommand.__init__(self)
        #self.flag = self.flag | self.objArgOnly


    def pickedVerticesToAtoms(self, geom, vertInd):
        """
        This function gets called when a picking or drag select event has
        happened. It gets called with a geometry and the list of vertex
        indices of that geometry that have been picked.
        This function is in charge of turning these indices into an AtomSet
        This function takes the following arguments:
        geom   : geometry picked, instance of a class derived from DejaVu2.Geom
                 (IndexedPolygons, IndexedPolylines.....)
        vertInd: list of integer representing the indices of the picked
                 vertices in the given geometry geom. 
        """

        # this function gets called when a picking or drag select event has
        # happened. It gets called with a geometry and the list of vertex
        # indices of that geometry that have been selected.
        # This function is in charge of turning these indices into an AtomSet
        ss = geom.SS
        l = []
        for vi in vertInd:
            resInd = ss.exElt.getResIndexFromExtrudeVertex( vi )
            try:
                l.append(ss.children[int(resInd)].atoms.get('CA')[0])
            except:
                l.append(ss.children[int(resInd)].atoms[0])
        return AtomSet( AtomSet( l ) )


    def atomPropToVertices(self, geom, residues, propName, propIndex=None):
        """Function called to compute the array of properties"""
        if residues is None or len(residues)==0 : return None
        propVect = []
        if not propIndex is None:
            propIndex = 'secondarystructure'
        for r in residues:
            try:
                prop = getattr(r.atoms.get('CA')[0], propName)
            except  IndexError:
                prop = getattr(r.atoms[0], propName)
            if not propIndex is None:
                propVect.append(prop.get(propIndex, prop['lines']))
            else:
                propVect.append(prop)
        geom.SS.exElt.setResProperties(propVect, propName, residues)
        properties = geom.SS.exElt.getExtrudeProperties( residues, propName )
        return properties


    def onAddObjectToViewer(self, obj):
        self.objectState[obj] = {'onAddObjectCalled':True}
        # private flag to specify whether or not the geometries for the SS
        # have been created.
        obj.__hasSSGeom = 0
        if self.app().commands.has_key('dashboard'):
            self.app().dashboard.resetColPercent(obj, '_showRibbonStatus')


    def createGeometries(self, obj):
        if obj.__hasSSGeom :
            return
        from DejaVu2.Geom import Geom
        geomC = obj.geomContainer

        if not geomC.geoms.has_key('secondarystructure'):
            
            t = Geom('secondarystructure', shape=(0,0), protected=True)
            geomC.addGeom( t, parent=geomC.masterGeom, redo=0 ) 
        else:
            t = geomC.geoms['secondarystructure']
        
        for a in obj.allAtoms:
            a.colors['secondarystructure']=(1.,1.,1.)
            a.opacities['secondarystructure']=1.

        for c in obj.chains:
            if not hasattr(c, 'secondarystructureset'):
                continue
            for ss in c.secondarystructureset:
                name = "%s%s"%(ss.name, ss.chain.id)
                g = IndexedPolygons(name, visible=0, pickableVertices=1, protected=True,)
                if self.app().userpref['Sharp Color Boundaries for MSMS']['value'] == 'blur':
                    g.Set(inheritSharpColorBoundaries=False, sharpColorBoundaries=False,)
                #g.RenderMode(GL.GL_FILL, face=GL.GL_FRONT, redo=0)
                #g.Set(frontPolyMode=GL.GL_FILL,redo=0)

                g.SS = ss
                
                geomC.atomPropToVertices[name] = self.atomPropToVertices
                geomC.geomPickToAtoms[name] = self.pickedVerticesToAtoms
                geomC.geomPickToBonds[name] = None
                geomC.addGeom(g, parent=t, redo=0 )
                self.managedGeometries.append(g)
                #geomC.addGeom(g,self,parent=t, redo=0 )
                geomC.atoms[name] = ResidueSet()
                
        obj.__hasSSGeom = 1


    def onAddCmdToApp(self):
        self.app().lazyLoad("secondaryStructureCmds",
               commands = ['computeSecondaryStructure', 'displayExtrudedSS'],
                         package="PmvApp")

        self.app().lazyLoad("extrusionCmds",
               commands=['computeSheet2D', "Nucleic_Acids_properties"], package="PmvApp")


    def clean(self, obj):
        if not hasattr(obj, 'chains'): return
        for c in obj.chains:
            if hasattr(c, 'residuesInSS'):
                delattr(c, 'residuesInSS')

            if not hasattr(c, 'secondarystructureset'):
                continue
            for ss in c.secondarystructureset:
                # Have to remove specifically geoms.SS and geoms.mol
                # from the geomContainer and the viewer
                g = obj.geomContainer.geoms[ss.name+c.id]
                del(g.SS)
                del(g.mol)
                g.Set(visible=0, tagModified=False)
                g.protected = False
                event = RemoveGeometryEvent(g)
                self.app.eventHandler.dispatchEvent(event)
                # the application's GUI should add a listener for this event,
                # and the method to:
                #self.app().GUI.VIEWER.RemoveObject(g)
                del obj.geomContainer.geoms[ss.name+c.id]
                del obj.geomContainer.atoms[ss.name+c.id]
        obj.__hasSSGeom=0


    def onRemoveObjectFromViewer(self, obj):
        if self.objectState.has_key(obj):
            self.objectState.pop(obj)
        if self.app().undoableDelete__:
            return
        if not hasattr(obj, 'chains'): return
        for c in obj.chains:
            if hasattr(c, 'residuesInSS'):
                delattr(c, 'residuesInSS')

            if not hasattr(c, 'secondarystructureset'):
                continue
            for ss in c.secondarystructureset:
                # Have to remove specifically geoms.SS and geoms.mol
                # from the geomContainer and the viewer
                g = obj.geomContainer.geoms[ss.name+c.id]
                del(g.SS)
                del(g.mol)
                g.Set(visible=0, tagModified=False)
                g.protected = False
                event = RemoveGeometryEvent(g)
                self.app.eventHandler.dispatchEvent(event)
                # the application's GUI should add a listener for this event,
                # and the method to:
                #self.app().GUI.VIEWER.RemoveObject(g)
                del obj.geomContainer.geoms[ss.name+c.id]
                del obj.geomContainer.atoms[ss.name+c.id]
        obj.__hasSSGeom=0
        


    def checkArguments(self, nodes, shape1=None, shape2=None, frontcap=True,
                       endcap=True, arrow=True, nbchords=8, gapBeg=0, gapEnd=0,
                       larrow=2, display=True, width=1.2, height=0.2, radius=0.1,
                       updateNucleicAcidsPropertiesGUI=False,
                       only=True, negate=False):
        """Required Arguments:\n
        nodes ---  TreeNodeSet holding the current selection
                   (mv.getSelection()) \n
        Optional Arguments:\n
        shape1 &
        shape2 --- DejaVu2.Shapes.Shape2D objects. shape1 will be used to \n
                  represent the helix and strand, shape2 to represent coils and\n
                  turns.\n
        frontcap &
        endcap --- Boolean flag when set to True a  cap will be added to the \n
                   geom either at the front or at the end \n
        arrow --- Boolean flag when set to True an arow will be added to the \n
                   geometry representing the strand. \n
        nbchords --- Nb of points per residues in the smooth array \n
        gapBeg&  
        gapEnd --- defines gap at the beginning or the end of each residue. \n
        larrow  --- length of the arrow if arrow boolean flag set to 1 \n
        display --- Boolean flag when set to True the displaySecondaryStructure
                   is called automatically
        width, height, radius --- if shape1 is not specified, these parameters \n
               are used to create shape1 (Rectangle2D(withd, height)) \n
               and shape2 (Circle2D(radius)) .          
        """

        if isinstance (nodes, str):
            self.nodeLogString = "'"+nodes+"'"
        nodes = self.app().expandNodes(nodes)

        assert isinstance(nbchords, int)
        assert gapEnd<=len(nodes)
        assert gapBeg<=len(nodes)
        if shape1:
            assert isinstance(shape1 , Shape2D)
        if shape2:
            assert isinstance(shape2 , Shape2D)
        assert frontcap in (True,False, 1, 0)
        assert endcap in (True, False, 1, 0)
        assert arrow in (True, False, 1, 0)
        assert display in (True, False, 1, 0)
        assert isinstance (larrow, (int, float))
        assert isinstance (width, (int, float))
        assert isinstance (height, (int, float))
        assert isinstance (radius, (int, float))
        kw = {}
        kw['shape1'] = shape1
        kw['shape2'] = shape2
        kw['frontcap'] = frontcap
        kw['endcap'] = endcap
        kw['arrow'] = arrow
        kw['nbchords'] = nbchords
        kw['gapBeg'] = gapBeg
        kw['gapEnd'] = gapEnd
        kw['larrow'] = larrow
        kw['display'] = display
        kw['width'] = width
        kw['height'] = height
        kw['radius'] = radius
        kw['updateNucleicAcidsPropertiesGUI'] = updateNucleicAcidsPropertiesGUI
        kw['only'] = only
        kw['negate'] = negate
        #print "kw.has_key('only')=", kw.has_key('only')
        #print kw.get('only', 'no_value')
        return (nodes,), kw
        

    def doit(self, nodes, shape1=None, shape2=None, frontcap=True, endcap=True,
             arrow=True, nbchords=8, gapBeg=0, gapEnd=0, larrow = 2,
             display=True, width=1.2, height=0.2, radius=0.1,
             updateNucleicAcidsPropertiesGUI=False,
             only=True, negate=False):
        """ nodes, shape1, shape2=None, frontcap=True, endcap=True, arrow=True,
        nbchords=8, gapBeg=0, gapEnd=1, display=True"""
        #print "2: kw.has_key('only')=", kw.has_key('only'), ':', 
        #print kw.get('only', 'no_value')
        shape1o = shape1
        shape2o = shape2
        molecules, residueSets = self.app().getNodesByMolecule(nodes, Residue)
        if shape1 is None:
            shape1 = Rectangle2D(width=width, height=height, vertDup=1)
            shape2 = Circle2D(radius=radius)

        # highlight selection
        selMols, selResidues = self.app().getNodesByMolecule(self.app().activeSelection.get(),
                                                          Residue)
        molSelectedResiduesDict = dict( zip( selMols, selResidues) )
        # Create a sheet2 object.
        for mol, residues in map(None, molecules, residueSets):
            try:
                if not self.objectState.has_key(mol):
                    self.onAddObjectToViewer(mol)
                if not mol.hasSS:
                    # Compute the secondarystructure if not there
                    self.app().computeSecondaryStructure(mol)
                if not hasattr(mol,'__hasSSGeom') or not mol.__hasSSGeom:
                    # Need here to change
                    self.createGeometries(mol)
                reswithss = residues.get(lambda x:
                                         hasattr(x, 'secondarystructure'))
                if reswithss is None:
                    raise RuntimeError, "%s: no secondary structure in specified nodes of molecule "% (self.name, mol.name)

                selectionSS = reswithss.secondarystructure.uniq()
                chains = residues.parent.uniq()

                # highlight selection
                if molSelectedResiduesDict.has_key(mol) and len(molSelectedResiduesDict[mol]) > 0:
                    lHighlight = True
                else:
                    lHighlight = False

                for i in range(len(chains)):
                    chain = chains[i]
                    chain.ssExtrusionParams = { # used to save session
                        'shape1' : shape1o,
                        'shape2' : shape2o,
                        'frontcap' : frontcap,
                        'endcap' : endcap,
                        'arrow' : arrow,
                        'nbchords' : nbchords,
                        'gapBeg' : gapBeg,
                        'gapEnd' : gapEnd,
                        'larrow' : larrow
                        }
                    newsheet = 0
                    if not hasattr(chain, 'sheet2D'):
                        chain.sheet2D = {}
                    if not hasattr(chain,'secondarystructureset'):
                        self.app().warningMsg('%s: no secondary structure set for chain %s in molecule %s'%(self.name, chain.id, mol.name))
                        chain.sheet2D['ssSheet2D'] = None
                        continue

                    ssSet = chain.secondarystructureset

                    # 1- Check if the sheet2D for a secondary structure has been
                    # computed already.
                    if chain.sheet2D.has_key('ssSheet2D'):
                        if chain.sheet2D['ssSheet2D'] is None:
                            newsheet = 0
                            continue
                        elif chain.sheet2D['ssSheet2D'].chords != nbchords:
                            rt = chain.ribbonType()
                            if rt=='NA':
                                ExtrudeNA(chain)
                                newsheet = 1
                            elif rt=='AA':
                                self.app().computeSheet2D(chain, 'ssSheet2D',
                                   'CA','O', buildIsHelix=1,
                                   nbchords=nbchords)
                                newsheet = 1
                            else:
                                newsheet = 0
                        else:
                            newsheet = 0

                    elif not chain.sheet2D.has_key('ssSheet2D'):
                        rt = chain.ribbonType()
                        if rt=='NA':
                            ExtrudeNA(chain)
                            newsheet = 1
                        elif rt=='AA':
                            self.app().computeSheet2D(chain, 'ssSheet2D',
                                                   'CA', 'O',buildIsHelix=1,
                                                   nbchords=nbchords)
                            newsheet = 1
                        else:
                            newsheet = 0

                    if newsheet:
                        sd = chain.sheet2D['ssSheet2D']
                        # then create a pointer to the sheet2D for each secondary structures.
                        ssSet.sheet2D = sd
                        if sd is None : continue
                    # Do the extrusion ONLY for the ss having a residue in the
                    # selection
                    removeSS =[]
                    #from PmvApp.selectionCommands import sets__
                    for SS in ssSet:
                        # test here if all the residues of the sselt are
                        # in the residue set used
                        # to compute the sheet2D. if not remove the ss.
                        if SS.sheet2D is None:
                            continue
                        #if filter(lambda x, rs = SS.sheet2D.resInSheet:
                        #          not x in rs, SS.residues):
                        if [x for x in SS.residues if not x in SS.sheet2D.resInSheet]:
                            self.app().warningMsg("%s: Removing %s from secondary structure set(molecule %s). One or more residues doesn't have CA and O"%(self.name, SS.name, mol.name))
                            # remove the SS from the set and etc....
                            #delattr(SS.residues, 'secondarystructure')
                            #ssSet.remove(SS)
                            removeSS.append(SS)
                            name = "%s%s"%(SS.name, SS.chain.id)
                            setName = mol.name+':'+name[-1]+':'+name[:-1]
                            if self.app().sets.has_key(setName):
                                del self.app().sets[setName]
                            #del sets__[mol.name+':'+name[-1]+':'+name[:-1]]
                            g = mol.geomContainer.geoms[name]
                            g.protected = False
                            event = RemoveGeometryEvent(g)
                            self.app.eventHandler.dispatchEvent(event)
                            # the application's GUI should add a listener for
                            # this event, and the method to:
                            # self.app().GUI.VIEWER.RemoveObject(g)
                            continue
                        name = "%s%s"%(SS.name, SS.chain.id)

                        if not SS in selectionSS:
                            continue
                        if isinstance(SS, Strand):
                            arrowf = arrow
                        else:
                            arrowf = 0
                        if not shape2 is None:
                            if SS.__class__.__name__ in ['Strand', 'Helix']:
                                SS.exElt = ExtrudeSSElt(
                                    SS, shape1, gapEnd , gapBeg, frontcap, endcap,
                                    arrowf,larrow)
                            elif SS.__class__.__name__ in ['Coil', 'Turn']:
                                rt = chain.ribbonType()
                                if rt=='NA':
                                    NAp = self.app().Nucleic_Acids_properties
                                    if NAp.isLoader():
                                        NAp = NAp.loadCommand()
                                    sc = max(NAp.scale_pyrimidine, NAp.scale_purine)
                                    #shape2 = Circle2D(radius=sc/2.5)
                                    shape3 = Circle2D(radius=sc/5.)
                                    SS.exElt = ExtrudeSSElt(
                                        SS, shape3, gapEnd, gapBeg, frontcap,
                                        endcap, arrowf)
                                elif rt=='AA':
                                    SS.exElt = ExtrudeSSElt(
                                        SS, shape2, gapEnd, gapBeg, frontcap,
                                        endcap, arrowf)

                        else:
                            SS.exElt = ExtrudeSSElt(SS, shape1, gapEnd , gapBeg,
                                                    frontcap, endcap, arrowf,
                                                    larrow)
                        resfaces, resfacesDict = SS.exElt.getExtrudeResidues(SS.residues)
                        g = mol.geomContainer.geoms[name]
    ##                     # MS triangulate faces
    ##                     trifaces = []
    ##                     for f in resfaces:
    ##                         trifaces.append( (f[0],f[1],f[3]) )
    ##                         if f[2]!=f[3]:
    ##                             trifaces.append( (f[1],f[2],f[3]) )

                        # highlight selection
                        g.resfacesDict = resfacesDict
                        highlight = []
                        if lHighlight is True:# and chain in residueSet :
                            highlight = [0]*len(SS.exElt.vertices)
                            for lResidue in molSelectedResiduesDict[mol]:
                                if resfacesDict.has_key(lResidue):
                                    for lFace in resfacesDict[lResidue]:
                                        for lVertexIndex in lFace:
                                            highlight[int(lVertexIndex)] = 1

                        g.Set(vertices=SS.exElt.vertices,
                              highlight=highlight,
                              faces = resfaces,
    ##                          faces=trifaces,
                              vnormals=SS.exElt.vnormals, redo=0,
                              tagModified=False)

                        if chain.ribbonType()=='NA':
                            geom_bases = Add_Nucleic_Bases(
                                g, self.app().Nucleic_Acids_properties)
                            event = AddGeometryEvent(geom_bases, parent=g, redo=False)
                            self.app.eventHandler.dispatchEvent(event)
                            # the GUI of the application should create this event
                            # listenter with the method that will: 
                            #self.app().GUI.VIEWER.AddObject(geom_bases, parent=g)
                            if geom_bases not in g.children: 
                                g.children.append(geom_bases)
                                geom_bases.parent = g
                                #geom_bases.fullName = g.fullName+'|'+geom_bases.name
                    for SS in removeSS:
                        delattr(SS.residues, 'secondarystructure')
                        ssSet.remove(SS)
                        
                self.app()._executionReport.addSuccess('Extruded SS for molecule %s successfully'%
                    mol.name, obj=residues)

            except:
                msg = 'Error while displaying lines for molecule %s'%mol.name
                self.app().errorMsg(sys.exc_info(), msg, obj=residues)
        if display:
            kw = {'only':True, 'negate':negate}
            #print "calling displayExtrudedSS with ", kw
            self.app().displayExtrudedSS(*(nodes,),  **kw)
#            gg = g.viewer.FindObjectByName('root|1dwb_0|secondarystructure|Coil1L')
#            print 'AFTER DISPLAY', gg
            #self.app().displayExtrudedSS(nodes)

        event = EditGeomsEvent(
            'SSextrude', [nodes,[shape1, shape2, frontcap, endcap,
                                 arrow, nbchords, gapBeg, gapEnd, larrow,
                                 display, updateNucleicAcidsPropertiesGUI]])
        self.app().eventHandler.dispatchEvent(event)

        


class ExtrudeSecondaryStructureCommandUnic(ExtrudeSecondaryStructureCommand):
    """The ExtrudeCommand allows the user to represent the secondary structure 
	elements by extruding 2D geometries along a 3D path
    Package : PmvApp \n
    Module  : secondaryStructureCmds \n
    Class   : ExtrudeSecondaryStructureCommand \n
    Command name : extrudeSecondaryStructure \n
    Synopsis:\n
            None <--- extrudeSecondaryStructure(nodes, shape1=None, shape2=None,frontcap=1, endcap=True, arrow=True, nbchords=8, gapBeg=False,gapEnd=False, larrow=2, display=True) \n
    Required Arguments:\n    
        nodes ---  TreeNodeSet holding the current selection(mv.getSelection()) \n
    Optional Arguments:\n    
        shape1 &
        shape2 --- DejaVu2.Shapes.Shape2D objects. shape1 will be used to \n
                  represent the helix and strand, shape2 to represent coils and \n  
                  turns. \n
        frontcap &
        endcap   --- Boolean flag when set to True a  cap will be added to the \n
                   geom either at the front or at the end \n
        arrow  --- Boolean flag when set to True an arow will be added to the \n
                   geometry representing the strand. \n
        nbchords --- Nb of points per residues in the smooth array \n
        gapBeg&  
        gapEnd  --- defines gap at the beginning or the end of each residue. \n
        larrow  --- lenth of the arrow if arrow boolean flag set to 1 \n
        display  --- Boolean flag when set to True the displaySecondaryStructure
                   is called automatically
    """

    def __init__(self):
        ExtrudeSecondaryStructureCommand.__init__(self)

    def createGeometries(self, obj):
        if hasattr(obj,'__hasSSGeom') :
            return
        from DejaVu2.Geom import Geom
        geomC = obj.geomContainer
        if not geomC.geoms.has_key('SS'):
            
            t = Geom('SS', shape=(0,0), protected=True)
            geomC.addGeom( t, parent=geomC.masterGeom, redo=0 ) 
        else:
            t = geomC.geoms['SS']
        
        for a in obj.allAtoms:
            a.colors['SS']=(1.,1.,1.)
            a.opacities['SS']=1.
        for c in obj.chains:
            #if not hasattr(c, 'secondarystructureset'):
            #    continue
            #for ss in c.secondarystructureset:
                name = "SS%s"%(c.id)
                g = IndexedPolygons(name, visible=0, pickableVertices=1, protected=True,)
                if self.app().userpref['Sharp Color Boundaries for MSMS']['value'] == 'blur':
                    g.Set(inheritSharpColorBoundaries=False, sharpColorBoundaries=False,)
                g.Set(frontPolyMode=GL.GL_FILL)

                #g.SS = ss
                
                geomC.atomPropToVertices[name] = self.atomPropToVertices
                geomC.geomPickToAtoms[name] = self.pickedVerticesToAtoms
                geomC.geomPickToBonds[name] = None
                geomC.addGeom(g, parent=t, redo=0 )
                self.managedGeometries.append(g)
                #geomC.addGeom(g,self,parent=t, redo=0 )
                geomC.atoms[name] = ResidueSet()
                atoms = c.findType(Atom)
                for a in atoms:
                    a.colors[name]=(1.,1.,1.)
                    a.opacities[name]=1.								
                for ss in c.secondarystructureset:
                    sname = "%s%s"%(ss.name, ss.chain.id)
                    geomC.atoms[sname] = ResidueSet()               
        obj.__hasSSGeom = 1


    def doit(self, nodes, shape1=None, shape2=None, frontcap=True, endcap=True,
             arrow=True, nbchords=8, gapBeg=0, gapEnd=0, larrow = 2,
             display=True, width=1.2, height=0.2, radius=0.1,
             updateNucleicAcidsPropertiesGUI=False, only=True, negate=False):
        """ nodes, shape1, shape2=None, frontcap=True, endcap=True, arrow=True,
        nbchords=8, gapBeg=0, gapEnd=1, display=True"""
        #print "2: kw.has_key('only')=", kw.has_key('only'), ':', 
        #print kw.get('only', 'no_value')

        molecules, residueSets=self.app().getNodesByMolecule(nodes, Residue)
        if len(molecules)==0: return
        if shape1 is None:
            shape1 = Rectangle2D(width=width, height=height, vertDup=1)
            shape2 = Circle2D(radius=radius)

        # highlight selection
        selMols, selResidues = self.app().getNodesByMolecule(self.app().activeSelection.get(), Residue)
        molSelectedResiduesDict = dict( zip( selMols, selResidues) )

        # Create a sheet2 object.
        for mol, residues in map(None, molecules, residueSets):
            try:
                if not mol.hasSS:
                    # Compute the secondarystructure if not there
                    self.app().computeSecondaryStructure(mol)
                if not hasattr(mol,'__hasSSGeom') or not mol.__hasSSGeom:
                    # Need here to change
                    self.createGeometries(mol)
                reswithss = residues.get(lambda x:
                                         hasattr(x, 'secondarystructure'))
                if reswithss is None:
                    raise RuntimeError, '%s: no secondary structure in specified nodes for molecule %s' % (self.name, mol)

                selectionSS = reswithss.secondarystructure.uniq()
                chains = residues.parent.uniq()

                # highlight selection
                if molSelectedResiduesDict.has_key(mol) and len(molSelectedResiduesDict[mol]) > 0:
                    lHighlight = True
                else:
                    lHighlight = False

                for i in range(len(chains)):
                    chain = chains[i]
                    newsheet = 0
                    if not hasattr(chain, 'sheet2D'):
                        chain.sheet2D = {}

                    if not hasattr(chain,'secondarystructureset'):
                        self.app().warningMsg('%s:no secondary structure set for chain: %s in molecule %s'%(self.name, chain.id, mol.name))
                        chain.sheet2D['ssSheet2D'] = None
                        continue

                    ssSet = chain.secondarystructureset
                    # 1- Check if the sheet2D for a secondary structure has been
                    # computed already.
                    if chain.sheet2D.has_key('ssSheet2D'):
                        if chain.sheet2D['ssSheet2D'] is None:
                            newsheet = 0
                            continue
                        elif chain.sheet2D['ssSheet2D'].chords != nbchords:
                            rt = chain.ribbonType()
                            if rt=='NA':
                                ExtrudeNA(chain)
                                newsheet = 1
                            elif rt=='AA':
                                self.app().computeSheet2D(chain, 'ssSheet2D',
                                   'CA','O', buildIsHelix=1,
                                   nbchords=nbchords)
                                newsheet = 1
                            else:
                                newsheet = 0
                        else:
                            newsheet = 0

                    elif not chain.sheet2D.has_key('ssSheet2D'):
                        rt = chain.ribbonType()
                        if rt=='NA':
                            ExtrudeNA(chain)
                            newsheet = 1
                        elif rt=='AA':
                            self.app().computeSheet2D(chain, 'ssSheet2D',
                                                   'CA', 'O',buildIsHelix=1,
                                                   nbchords=nbchords)
                            newsheet = 1
                        else:
                            newsheet = 0
                    if newsheet:
                        sd = chain.sheet2D['ssSheet2D']
                        # then create a pointer to the sheet2D for each secondary structures.
                        ssSet.sheet2D = sd
                        if sd is None : continue
                    # Do the extrusion ONLY for the ss having a residue in the
                    # selection
                    removeSS =[]
                    faces=[]
                    vertices=[]	
                    normals=[]							
                    #from PmvApp.selectionCommands import sets__
                    name = "SS"+chain.id
                    g = mol.geomContainer.geoms[name]
                    for SS in ssSet:
                        # test here if all the residues of the sselt are
                        # in the residue set used
                        # to compute the sheet2D. if not remove the ss.
                        if SS.sheet2D is None:
                            continue
                        if [x for x in SS.residues if not x in SS.sheet2D.resInSheet]:
                            self.app().warningMsg("%s: Removing %s from secondary structure set (%s). One or more residues doesn't have CA and O"%(self.name, SS.name, mol.name))

                            # remove the SS from the set and etc....
                            #delattr(SS.residues, 'secondarystructure')
                            #ssSet.remove(SS)
                            removeSS.append(SS)
                            #name = "%s%s"%(SS.name, SS.chain.id)
                            #del self.app().sets[mol.name+':'+name[-1]+':'+name[:-1]]
                            #del sets__[mol.name+':'+name[-1]+':'+name[:-1]]
                            #g = mol.geomContainer.geoms[name]
                            #g.protected = False
                            #if self.app().hasGui:self.app().GUI.VIEWER.RemoveObject(g)
                            continue
                        name = "%s%s"%(SS.name, SS.chain.id)

                        if not SS in selectionSS:
                            continue
                        if isinstance(SS, Strand):
                            arrowf = arrow
                        else:
                            arrowf = 0
                        if not shape2 is None:
                            if SS.__class__.__name__ in ['Strand', 'Helix']:
                                SS.exElt = ExtrudeSSElt(SS,shape1,gapEnd ,
                                                        gapBeg, frontcap, endcap,
                                                        arrowf,larrow)
                            elif SS.__class__.__name__ in ['Coil', 'Turn']:
                                if chain.ribbonType()=='NA':
                                    NAp = self.app().Nucleic_Acids_properties
                                    sc = max(NAp.scale_pyrimidine, NAp.scale_purine)
                                    shape2 = Circle2D(radius=sc/2.5)                                    
                                SS.exElt = ExtrudeSSElt(SS,shape2, gapEnd,
                                                        gapBeg, frontcap, endcap,
                                                        arrowf)

                        else:
                            SS.exElt = ExtrudeSSElt(SS, shape1, gapEnd , gapBeg,
                                                    frontcap, endcap, arrowf,
                                                    larrow)
                        resfaces, resfacesDict = SS.exElt.getExtrudeResidues(SS.residues)
                        #g = mol.geomContainer.geoms[name]
    ##                     # MS triangulate faces
    ##                     trifaces = []
    ##                     for f in resfaces:
    ##                         trifaces.append( (f[0],f[1],f[3]) )
    ##                         if f[2]!=f[3]:
    ##                             trifaces.append( (f[1],f[2],f[3]) )

                        # highlight selection
                        g.resfacesDict = resfacesDict
                        highlight = []
                        if lHighlight is True:# and chain in residueSet :
                            highlight = [0]*len(SS.exElt.vertices)
                            for lResidue in molSelectedResiduesDict[mol]:
                                if resfacesDict.has_key(lResidue):
                                    for lFace in resfacesDict[lResidue]:
                                        for lVertexIndex in lFace:
                                            highlight[int(lVertexIndex)] = 1

                        faces.extend(numpy.array(resfaces)+len(vertices))
                        vertices.extend(SS.exElt.vertices)
                        normals.extend(SS.exElt.vnormals)										
                        if chain.ribbonType()=='NA':
                            geom_bases = Add_Nucleic_Bases(g, 
                                               self.app().Nucleic_Acids_properties)
                            event = AddGeometryEvent(geom_bases, parent=g, redo=False)
                            self.app.eventHandler.dispatchEvent(event)
                            if geom_bases not in g.children:
                                g.children.append(geom_bases)
                                geom_bases.parent = g
                                #geom_bases.fullName = g.fullName+'|'+geom_bases.name
                    g.Set(vertices=vertices, highlight=highlight,
                          faces=faces,  vnormals=normals, redo=0,
                          tagModified=False)                       

                    for SS in removeSS:
                        delattr(SS.residues, 'secondarystructure')
                        ssSet.remove(SS)

                    atoms = chain.findType(Atom)
                    self.app().bindGeomToMolecularFragment(g, atoms)
                self.app()._executionReport.addSuccess('extruded SS for molecule %s successfully'%
                        mol.name, obj=residues)
            except:
                msg = 'Error while extruding SS for molecule %s'%mol.name
                self.app().errorMsg(sys.exc_info(), msg, obj=residues)
        if display:
            kw = {}
            if kw.get('only', 0):
                kw['only'] = 1
            kw['negate'] = negate
            #print "calling displayExtrudedSS with ", kw
            apply(self.app().displayExtrudedSS,(nodes,),  kw)


class DisplayExtrudedSSCommand(DisplayCommand):

    """ The DisplaySSCommand displays the geometries representing the secondary structure elements of the current selection.To execute this command use the 'Display Secondary Structure' entry under the 'Display' menu in the menu bar. \n
    Package : PmvApp \n
    Module  : secondaryStructureCmds \n
    Class   : DisplayExtrudedSSCommand \n
    Command name : displaySecondaryStructure \n
    Synopsis:\n
        None <- displaySecondaryStructure(nodes, only=False, negate=False) \n
    Required Arguments:\n
        nodes --- TreeNodeSet holding the current selection \n
    Optional Arguments:\n
        only --- allows the user to display only the current selection when set to 1 \n
        negate --- allows to undisplay the current selection when set to 1. \n
    This command is undoable.
    """
                                       
    def onAddCmdToApp(self):
        self.app().lazyLoad("secondaryStructureCmds",
            commands=['computeSecondaryStructure', 'extrudeSecondaryStructure'],
                         package="PmvApp")
            
        self.app().lazyLoad("extrusionCmds",
                            commands=['computeSheet2D'], package="PmvApp")

        self.app().eventHandler.registerListener(AfterDeleteAtomsEvent, self.handleAfterDeleteAtoms)


    def handleAfterDeleteAtoms(self, event):
        """Function to update geometry objects created by this command
upon atom deletion.
\nevent --- instance of a VFEvent object
"""
        # split event.objects into atoms sets per molecule
        molecules, ats = self.app().getNodesByMolecule(event.objects)
        # loop over molecules to update geometry objects
        for mol, atomSet in zip(molecules, ats):
            #if no backbone atoms are deleted ss does not change
            if len(atomSet.get("backbone")) == 0:
                continue
            
            for atom in atomSet.get("CA"):
                atom.parent.hasCA = False
                atom.parent.CAatom = None
            for atom in atomSet.get("O"):
                atom.parent.hasO = False
                atom.parent.Oatom = None
            
            if not  mol.geomContainer.geoms.has_key('secondarystructure'): continue
            kw = mol.geomContainer.geoms['secondarystructure'].kw.copy()
            #cleanup SS atomset in geomcontainer             
            self.app().computeSecondaryStructure.clean(mol)
            
            #this for molparser.hasSsDataInFile to trurn false
            while 'HELIX' in mol.parser.keys:
                mol.parser.keys.remove('HELIX')
            while 'SHEET' in mol.parser.keys:
                mol.parser.keys.remove('SHEET')
            while 'TURN' in mol.parser.keys:
                mol.parser.keys.remove('TURN')
                
            for chain in mol.chains:
                chain.ribbonType(noCache=True)
            
            #mol.hasSS = False
            mol.secondaryStructureFromPross()
            mol.__hasSSGeom=0

            self.app().extrudeSecondaryStructure(mol, display=0)
            self(mol, **kw)

##     def undoCmdBefore(self, nodes, only=False, negate=False, **kw):
##         if len(nodes)==0 : return
##         #molecules = nodes.top.uniq()
##         molecules, residueSets = self.getNodes(nodes)
##         #for mol, res in map(None, molecules, residueSets):
##         negateCmds = []
##         for mol in molecules:
##             resWithSS = mol.findType(Residue).get(
##                 lambda x:hasattr(x,'secondarystructure'))
##             if resWithSS is None:
##                 continue
##             SSinMol = resWithSS.secondarystructure.uniq()
##             #resWithSS = res.get(lambda x: hasattr(x,'secondarystructure'))
##             #SSinSel = resWithSS.secondarystructure.uniq()
##             #mol.geomContainer.atoms['secondarystructure']=resWithSS.atoms
##             set = ResidueSet()
##             if mol.geomContainer.geoms.has_key('SS'):
##                 for ch in mol.chains :
##                     set = set + mol.geomContainer.atoms['SS'+ch.id].parent
##                     if len(set)==0: # nothing is displayed
##                         negateCmds.append((self, (mol,), {'negate':True, 'redraw':True},))

##                     else:
##                         negateCmds.append((self, (set,), {'only':True, 'redraw':True}) )
##             else :
##                 for ss in SSinMol:
##                     set = set + mol.geomContainer.atoms[ss.name+ss.chain.id]
##                     if len(set)==0: # nothing is displayed
##                         negateCmds.append((self, (mol,), {'negate':True, 'redraw':True}))
##                     else:
##                         negateCmds.append((self,  (set,), {'only':True, 'redraw':True}))
##         if len(negateCmds):
##             return (negateCmds, self.name)


    def doit(self, nodes, only=False, negate=False, redraw=True):
        """ displays the secondary structure for the selected treenodes """

        #print self.name, "in display with only=", only, " and negate=", negate
        ###############################################################
        def drawResidues(SS, res, only, negate, uniq=False):
            mol = SS.chain.parent
            name = '%s%s'%(SS.name, SS.chain.id)
            _set = mol.geomContainer.atoms[name]
            inres = [x for x in _set if not x in res]
            if len(inres) == 0:
                # res and _set are the same
                if negate:
                    _set = ResidueSet()
                    setOff = res
                    setOn = None
                else:
                    _set = res
                    setOff = None
                    setOn = res
            else:
                # if negate, remove current res from displayed _set
                if negate :
                    setOff = res
                    setOn = None
                    _set = _set - res

                else: # if only, replace displayed _set with current res
                    if only:
                        setOff = _set - res
                        setOn = res
                        _set = res
                    else:
                        _set = res.union(_set)
                        setOff = None
                        setOn = _set

##             # if negate, remove current res from displayed _set
##             if negate :
##                 _set = _set - res

##             else: # if only, replace displayed _set with current res
##                 if only:
##                     _set = res
##                 else:
##                     _set = res.union(_set)
##             # now, update the geometries:
## 	    if len(_set)==0:
##                 mol.geomContainer.geoms[name].Set(visible=0, tagModified=False)
##                 mol.geomContainer.atoms[name] = ResidueSet()
##                 return

            #the rest is done only if there are some residues           
            
            mol.geomContainer.setAtomsForGeom(name, _set)
            #print _set
            if not hasattr(SS, 'exElt'):
                return setOn, setOff

            if isinstance(SS, Coil):
                gapBefore = SS.gapBefore
                gapAfter = SS.gapAfter
            else:
                gapBefore = gapAfter = False
                
            resfaces, resfacesDict = SS.exElt.getExtrudeResidues(
                _set, gapBefore, gapAfter)
            
##          # MS triangulate faces
##             trifaces = []
##             for f in resfaces:
##                 trifaces.append( (f[0],f[1],f[3]) )
##                 if f[2]!=f[3]:
##                     trifaces.append( (f[1],f[2],f[3]) )
##             g.Set(faces=trifaces, vnormals=SS.exElt.vnormals,
            if uniq:
                return setOn, setOff, resfaces, SS.exElt.vnormals, \
                       SS.exElt.vertices            
            g = mol.geomContainer.geoms[name]
            col = mol.geomContainer.getGeomColor(name)
            g.Set(faces=resfaces, vnormals=SS.exElt.vnormals,
                  visible=1, materials=col, inheritMaterial=False,
                  tagModified=False)
            if SS.chain.ribbonType()=='NA':
                faces = []
                colors = [] 
                for residue in _set:
                    faces.extend(residue._base_faces)
                    colors.extend(residue._coil_colors)
                    try:
                        atm = residue.atoms.get('CA')[0]
                    except:
                        atm = residue.atoms[0]
                    if len(residue._coil_colors) :
                        atm.colors["secondarystructure"] = residue._coil_colors[0]#??
                    else:
                        atm.colors["secondarystructure"] = residue._coil_colors
                g.children[0].Set(faces=faces)
                if colors:
                    g.Set(materials=colors, inheritMaterial=False)
                    
                if self.app().Nucleic_Acids_properties.color_backbone:
                    g.Set(inheritMaterial=False)
                    
                else:
                    g.Set(inheritMaterial=True)
                mol.geomContainer.atoms['Bases'] = ResidueSet()
                #mol.geomContainer.atoms[name] = ResidueSet()

            return setOn, setOff

###############################################################

        molecules, residueSets = self.app().getNodesByMolecule(nodes, Residue)
        setOn = ResidueSet([])
        setOff = ResidueSet([])
        for mol, residues in map(None, molecules, residueSets):
            try:
                if not mol.hasSS:
                    self.app().computeSecondaryStructure(mol)
                    self.app().extrudeSecondaryStructure(mol, display=0)
                reswithss = residues.get(lambda x:
                                         hasattr(x,'secondarystructure'))
                if reswithss is None:
                    raise RuntimeError, '%s: no secondary structure in specified nodes of %s molecule' % (self.name, mol.name)

                SSInSel = reswithss.secondarystructure.uniq()
                chainsInSel = residues.parent.uniq()
                for c in mol.chains:
                    if not hasattr(c, 'secondarystructureset'):
                        continue
                    if not hasattr(c, 'sheet2D'):
                        self.app().warningMsg("%s: chain '%s'(%s) does not have a sheet2D computed"%(self.name, c.full_name(), mol.name))
                        continue
                    elif (c.sheet2D.has_key('ssSheet2D') and \
                          c.sheet2D['ssSheet2D'] is None): continue
                    if mol.geomContainer.geoms.has_key('SS') :
                        faces=[]
                        vertices=[]	
                        normals=[]							
                        name = "SS"+c.id
                        g = mol.geomContainer.geoms[name]
                    SS, resInSS = self.getResiduesBySS(residues, c)
                    for s in xrange(len(c.secondarystructureset)):
                        ss = c.secondarystructureset[s]
                        res = resInSS[s]
                        if ss in SSInSel and not hasattr(ss, 'exElt') \
                           and negate == 0:
                            self.app().extrudeSecondaryStructure(res, display=0)

                        if mol.geomContainer.geoms.has_key('SS') :
                            son, sof, f, n, v = drawResidues(ss, res, only , negate ,uniq=True)
                            faces.extend(numpy.array(f)+len(v))
                            vertices.extend(v)
                            normals.extend(n)
                        else :
                            son, sof = drawResidues(ss, res, only , negate )

                        if son: setOn += son
                        if sof: setOff += sof

                    if mol.geomContainer.geoms.has_key('SS') :
                        g.Set(visible=not negate)
                kw = {"only":only, "negate":negate}
                if mol.geomContainer.geoms.has_key('secondarystructure'):
                    mol.geomContainer.geoms['secondarystructure'].kw = kw
                self.app()._executionReport.addSuccess('displayed secondary structure for molecule %s successfully'%
                    mol.name, obj=residues)
            except:
                msg = 'Error while displaying secondary structure for molecule %s'%mol.name
                self.app().errorMsg(sys.exc_info(), msg, obj=residues)
        event = EditGeomsEvent('SSdisplay', [nodes,[only, negate,redraw]],
                               setOn=setOn, setOff=setOff)
        self.app().eventHandler.dispatchEvent(event)

					
    def checkArguments(self, nodes, only=False, negate=False, redraw=True):
        """ Required Arguments:\n
            nodes  ---  TreeNodeSet holding the current selection \n
            Optional Arguments:\n
            only ---  flag when set to 1 only the current selection will be displayed as secondarystructures \n
            negate ---  flag when set to 1 undisplay the current selection"""
        if isinstance(nodes, str):
            self.nodeLogString = "'"+nodes+"'"
        nodes = self.app().expandNodes(nodes)

        kw = {}
        assert only in [True, False, 1, 0]
        assert negate in [True, False, 1, 0]
        kw['only'] = only
        kw['negate'] = negate
        kw['redraw']=redraw
        return (nodes,),kw
    

    def getResiduesBySS(self, residues, chain):
        resWithSS = residues.get(lambda x: hasattr(x, 'secondarystructure'))
        residuesInSS = []
        for ss in chain.secondarystructureset :
            res = resWithSS.get(lambda x, ss=ss:x.secondarystructure==ss)
            if res is None:
                res = ResidueSet()
            residuesInSS.append(res)
        return chain.secondarystructureset, residuesInSS


class UndisplayExtrudedSSCommand(DisplayCommand):
    """ UndisplaySSCommand is an interactive command to undisplay part of
    the molecule when represented as extruded secondary structure. \n
    Package : PmvApp \n
    Module  : secondaryStructureCmds \n
    Class   : UndisplayExtrudedSSCommand \n
    Command name : undisplaySecondaryStructure \n
    Synopsis:\n
         None <--- undisplaySecondaryStructure(nodes, **kw) \n
    Required Arguments:\n    
         nodes --- TreeNodeSet holding the current selection
    """
    def onAddCmdToApp(self):
        if not self.app().commands.has_key('displayExtrudedSS'):
            self.app().lazyLoad('secondaryStructureCmds',
                             commands=['displayExtrudedSS'], package='PmvApp')

        
    def checkArguments(self, nodes, **kw):
        """ nodes --- TreeNodeSet holding the current selection
        """
        if isinstance(nodes, str):
            self.nodeLogString = "'"+nodes+"'"
        nodes = self.app().expandNodes(nodes)

        kw = {'negate':1}
        return (nodes,), kw


    def doit(self, nodes, **kw):
        self.app().displayExtrudedSS(nodes, **kw)



class RibbonCommand(MVCommand):
    """ The RibbonCommand is a shortcut to visualize a traditional Ribbon
    representation of the current selection. It first executes getSSCommand
    then the extrudeSSCommand with the default values for all the parameters.
    This command is undoable. \n
    Package : PmvApp \n
    Module  : secondaryStructureCmds \n
    Class   : RibbonCommand \n
    Command name : ribbon \n
    Synopsis:\n    
        None <- ribbon(nodes, only=False, negate=False) \n
    Required Arguments:\n    
        nodes ---  TreeNodeSet holding the current selection \n
    Optional Arguments:\n   
        only --- flag when set to 1 only the current selection
                  will be displayed \n
        negate --- flag when set to 1 undisplay the current selection
    """
    def __init__(self):
        MVCommand.__init__(self)
        #self.flag = self.flag | self.objArgOnly
        #self.flag = self.flag | self.negateKw

   
    def undoCmdBefore(self, *args, **kw):
        return None # this prevents this function from trying to create a negation
                    # command since it only calls other VFCommands who can negate themselves
    
    def onAddCmdToApp(self):
        self.app().lazyLoad("secondaryStructureCmds",
           commands=['extrudeSecondaryStructure', 'displayExtrudedSS',
            'computeSecondaryStructure'], package="PmvApp")


    def checkArguments(self, nodes, shape1=None, shape2=None, frontcap=True,
                       endcap=True, arrow=True, nbchords=8, gapBeg=0, gapEnd=0,
                       larrow=2, display=True,  redraw=True, width=1.2,
                       height=0.2, radius=0.1, updateNucleicAcidsPropertiesGUI=False,
                       only=True, negate=False):
        """Required Arguments:\n
        nodes ---  TreeNodeSet holding the current selection
                   (mv.getSelection()) \n
        Optional Arguments (will be passed to extrudeSecondaryStructure() ):\n
        shape1 &
        shape2 --- DejaVu2.Shapes.Shape2D objects. shape1 will be used to \n
                  represent the helix and strand, shape2 to represent coils and\n
                  turns.\n
        frontcap &
        endcap --- Boolean flag when set to True a  cap will be added to the \n
                   geom either at the front or at the end \n
        arrow --- Boolean flag when set to True an arow will be added to the \n
                   geometry representing the strand. \n
        nbchords --- Nb of points per residues in the smooth array \n
        gapBeg&  
        gapEnd --- defines gap at the beginning or the end of each residue. \n
        larrow  --- length of the arrow if arrow boolean flag set to 1 \n
        display --- Boolean flag when set to True the displaySecondaryStructure
                   is called automatically \n
        only --- flag when set to 1 only the current selection
                  will be displayed \n
        negate ---  flag when set to 1 undisplay the current selection \n
        width, height, radius --- if shape1 is not specified, these parameters \n
               are used to create shape1 (Rectangle2D(withd, height)) \n
               and shape2 (Circle2D(radius)) 
                   """

        if isinstance (nodes, str):
            self.nodeLogString = "'"+nodes+"'"
        nodes = self.app().expandNodes(nodes)

        assert isinstance(nbchords, int)
        assert gapEnd<=len(nodes)
        assert gapBeg<=len(nodes)
        if shape1:
            assert isinstance(shape1 , Shape2D)
        if shape2:
            assert isinstance(shape2 , Shape2D)
        assert frontcap in (True,False, 1, 0)
        assert endcap in (True, False, 1, 0)
        assert arrow in (True, False, 1, 0)
        assert display in (True, False, 1, 0)
        assert isinstance (larrow, (int, float))
        assert isinstance (width, (int, float))
        assert isinstance (height, (int, float))
        assert isinstance (radius, (int, float))
        kw = {}
        kw['shape1'] = shape1
        kw['shape2'] = shape2
        kw['frontcap'] = frontcap
        kw['endcap'] = endcap
        kw['arrow'] = arrow
        kw['nbchords'] = nbchords
        kw['gapBeg'] = gapBeg
        kw['gapEnd'] = gapEnd
        kw['larrow'] = larrow
        kw['display'] = display
        kw['width'] = width
        kw['height'] = height
        kw['radius'] = radius
        kw['updateNucleicAcidsPropertiesGUI'] = updateNucleicAcidsPropertiesGUI
        kw['only'] = only
        kw['negate'] = negate
        kw['redraw']=redraw
        #print "kw.has_key('only')=", kw.has_key('only')
        #print kw.get('only', 'no_value')
        return (nodes,), kw

        
    def doit(self, nodes, only=False, negate=False, redraw=True, **kw):
        self.app().computeSecondaryStructure( nodes)
        kw.update({'only':only, 'negate':negate})
        #print self.name, "doit:", kw
        self.app().extrudeSecondaryStructure( nodes, **kw)


class ColorBySSElementType(ColorFromPalette):
    """Command to color the given geometry by secondary structure
    element. (Rasmol color code) \n
    Package : PmvApp \n
    Module  : secondaryStructureCmds \n
    Class   : ColorBySSElementType
     """
    def onAddCmdToApp(self):
        from PmvApp.pmvPalettes import SecondaryStructureType
        c = 'Color palette for secondary structure element type:'
        self.palette = ColorPalette(
            'SecondaryStructureType', SecondaryStructureType,
            info=c, lookupMember = 'structureType')

        if not self.app().commands.has_key('color'):
            self.app().lazyLoad('colorCmds', ['color'], 'PmvApp')
        self.undoCmdsString= self.app().color.name


    def getColors(self, nodes):
        res = nodes.findType(Residue)
        resWithSS = res.get(lambda x: hasattr(x, 'secondarystructure'))
        if resWithSS is None: return None, None
        return resWithSS, self.palette.lookup(resWithSS.secondarystructure)


    def getNodes(self, nodes, returnNodes=False):
        """expand nodes argument into a list of atoms and a list of
        molecules."""
        nodes = self.app().expandNodes(nodes)
        res = nodes.findType(Residue).uniq()
        resWithSS = res.get(lambda x: hasattr(x,'secondarystructure'))
        if resWithSS is None or len(resWithSS)==0:
            atoms = AtomSet()
            molecules = ProteinSet()
        else:
            atoms = resWithSS.atoms
            molecules= resWithSS.top.uniq()
        if returnNodes:
            return molecules, atoms, nodes
        else:
            return molecules, atoms

            
    def doit(self, nodes, geomsToColor):
        # this command do not require the color argument since colors are
        # gotten from a palette
        # we still can use the ColorCommand.undoCmdBefore but first we get
        # the colors. This also insures that the colors are not put inside the
        # command's log string
        #print self.name , "geomsToColor", geomsToColor
        #nodes is AtomSet 
        resWithSS, colors = self.getColors(nodes)
        if colors is None: return
        for g in geomsToColor:
            if len(colors)==1 or len(colors)!=len(nodes):
                for a in nodes:
                    a.colors[g] = tuple( colors[0] )
            else:
                for a, c in map(None, nodes, colors):
                    a.colors[g] = tuple(c)

        updatedGeomsToColor = []
        for mol in self.molSet:
            try:
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

                mol.geomContainer.updateColors(updatedGeomsToColor)
                self.app()._executionReport.addSuccess('%s: colored molecule %s successfully'% (self.name, mol.name))
            except:
                msg = 'Error while coloring secondary structure for molecule %s'% mol.name
                self.app().errorMsg(sys.exc_info(), msg, obj=self.atmSet)
        #geomEditEventss
        event = EditGeomsEvent("color", [nodes,[geomsToColor, colors, self.name[5:11]]])
        self.app().eventHandler.dispatchEvent(event)



commandClassFromName = {
    'computeSecondaryStructure' : [ComputeSecondaryStructureCommand, None],
    'extrudeSecondaryStructure' : [ExtrudeSecondaryStructureCommand, None],
    'extrudeSecondaryStructureUnic' : [ExtrudeSecondaryStructureCommandUnic, None],
    'displayExtrudedSS' : [DisplayExtrudedSSCommand, None],
    'colorBySecondaryStructure' : [ColorBySSElementType, None],
    'undisplayExtrudedSS' : [UndisplayExtrudedSSCommand,  None],
    'ribbon' : [RibbonCommand, None],
}

def initModule(viewer, gui=False):
    for cmdName, values in commandClassFromName.items():
        cmdClass, guiInstance = values
        viewer.addCommand(cmdClass(), cmdName, guiInstance)


