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

# -*- coding: utf-8 -*-
"""
Created on Mon Jul  7 14:15:36 2014

@author: ludo
"""
import numpy as np

from PmvApp.displayCmds import DisplayCommand
#from PmvApp.Pmv import DeleteAtomsEvent, AddAtomsEvent, EditAtomsEvent,\
#     ShowMoleculesEvent, DeleteGeomsEvent, AddGeomsEvent, EditGeomsEvent
import pdb
from PmvApp.Pmv import EditGeomsEvent
#from MolKit2.molecule import Atom, AtomSet, Molecule, MoleculeSet, BondSet
from MolKit2.molecule import Atom, Molecule, Residue, Chain
#from MolKit2.protein import Protein, Residue, Chain
#
from opengltk.OpenGL import GL

class DisplayHyperBalls(DisplayCommand):
    """ The displayHyperBalls command allows the user to display/undisplay the given nodes using a CPK representation, where each atom is represented with a sphere. A scale factor and the quality of the spheres are user
    defined parameters.
    \nPackage : PmvApp
    \nModule  : displayCommands
    \nClass   : displayHyperBalls
    \nCommand : displayHyperBalls
    \nSynopsis:\n
        None <- displayHyperBalls(nodes, only=False, negate=False, 
                           scaleFactor=1.0, cpkRad=0.0, quality=0,  
                           shrink=0.1, bScale=0.05, )\n
        nodes --- any set of MolKit2.Selection describing molecular components\n
        only --- Boolean flag specifying whether or not to only display
                     the current selection\n
        negate --- Boolean flag specifying whether or not to undisplay
                     the current selection\n
        scaleFactor --- specifies a scale factor that will be applied to the atom
                     radii to compute the spheres radii. (default is 1.0)\n
        quality  --- specifies the quality of the spheres. (default 10)\n
        cpkRad  --- offset , used to comute spheres radii
        propertyName --- if specified,  the property is used to compute the spheres \n
                        radii (Sphere Radii = cpkRad + property * scaleFactor), \n
        propertyLevel --- can be one of the following: 'Atom', 'Residue', 'Chain', 'Molecule'.
        keywords --- display, CPK, space filling, undisplay.\n
    """

        
    def onAddCmdToApp(self,):
        DisplayCommand.onAddCmdToApp(self)

        if not self.app().commands.has_key('assignAtomsRadii'):
            self.app().lazyLoad('editCmds', commands=['assignAtomsRadii'], package='PmvApp')

        from mglutil.util.colorUtil import ToHEX
        self.molDict = {'Molecule':Molecule,
                        'Atom':Atom, 'Residue':Residue, 'Chain':Chain}
        self.nameDict = {Molecule:'Molecule', Atom:'Atom', Residue:'Residue',
                         Chain:'Chain'}

        self.leveloption={}
        for name in ['Atom', 'Residue', 'Molecule', 'Chain']:
            col = self.app().levelColors[name]
            bg = ToHEX((col[0]/1.5,col[1]/1.5,col[2]/1.5))
            ag = ToHEX(col)
            self.leveloption[name]={#'bg':bg,
                                    'activebackground':bg, 'selectcolor':ag ,
                                    'borderwidth':3 , 'width':15}
        self.propValues = None
        self.getVal = 0
        #self.propertyLevel = self.nameDict[self.app().selectionLevel]
        self.propertyLevel = "Molecule"
        self.options={"CPK":[0.01,0.2,0.3],
                      "VDW":[0.01,0.01,0.99],
                      "LIC":[0.01,0.26,0.26],
                      "HBL":[0.3,0.4,0.4]}

        

    def onAddObjectToViewer(self, obj):
        """Adds the hb geometry and the cpk Atomset to the object's
        geometry container
        """
        print "onAddObjectToViewer"
        self.objectState[obj] = {'onAddObjectCalled':True}
        if self.app().commands.has_key('dashboard'):
            self.app().dashboard.resetColPercent(obj, '_showCPKStatus')

        obj.allAtoms.cpkScale = 1.0
        obj.allAtoms.cpkRad = 0.0
        geomC = obj.geomContainer
        # CPK spheres
        from DejaVu2.hpBalls import hyperBalls
        #from pyQutemol.qutemolSpheres import QuteMolSpheres
        g = hyperBalls("hpballs" ,visible=0, protected=True)
        geomC.addGeom(g, parent=geomC.masterGeom, redo=0)
        self.managedGeometries.append(g)
        geomC.geomPickToBonds['hpballs'] = None
        for atm in obj.allAtoms:
            atm.colors['hpballs'] = (1.0, 1.0, 1.0)
            atm.opacities['hpballs'] = 1.0
#        geomC.masterGeom.isScalable =1 
#        g.initialiseCamera()
#        self.vf.GUI.VIEWER.useMasterDpyList=0
        try :
            print "try changing useMasterDpyList"
            self.app().gui().viewer.useMasterDpyList=0
            print "ok ",self.app().gui().viewer.useMasterDpyList
        except :
            print "cant access viewer camera"

    def initializeMolForCmd(self, patoms):
        """Adds the hb geometry and the cpk Atom Group to the object's
        geometry container
        """
        obj = patoms.getMolecule()
        self.initializedFor[patoms.getMolecule()] = True
        #if self.app().commands.has_key('dashboard'):
        #    self.app().dashboard.resetColPercent(obj, '_showCPKStatus')
#        obj.cpkScale = 1.0
#        obj.cpkRad = 0.0
        geomC = obj.geomContainer
        #obj.setcpkScale = self.setCpkScale
        #obj.getcpkScale = self.getCpkScale
        #obj.setcpkScale(obj, 1.0)
        obj._cpkScale = np.ones( (len(patoms),), 'f')
        obj._cpkRad = np.zeros( (len(patoms),), 'f')
        
        #obj.setcpkRad = self.setCpkRad
        #obj.getcpkRad = self.getCpkRad
        #obj.setcpkRad(obj, 0.0)

        # CPK spheres
        if not geomC.geoms.has_key('hpballs'):
#            g = Spheres( "cpk", quality=0, visible=0, protected=True)
#            #faces=[[x] for x in xrange(len(geomC.allCoords))])
#            geomC.addGeom(g, parent=geomC.masterGeom, redo=0)
#            g.vertexSet.vertices.array = geomC.allCoords   
#            self.managedGeometries.append(g)
#            geomC.geomPickToBonds['cpk'] = None
#            g._firstDisplay = True # used in doit to use line colors for first display
            from DejaVu2.hpBalls import hyperBalls
            #from pyQutemol.qutemolSpheres import QuteMolSpheres
            g = hyperBalls("hpballs" ,visible=0, protected=True)
            geomC.addGeom(g, parent=geomC.masterGeom, redo=0)
            self.managedGeometries.append(g)
            geomC.geomPickToBonds['hpballs'] = None
#            obj._cpkRad = np.zeros( (len(patoms),), 'f')
#            for atm in obj.allAtoms:
#                atm.colors['hpballs'] = (1.0, 1.0, 1.0)
#                atm.opacities['hpballs'] = 1.0
            g._firstDisplay = True #
            #all coords : pmv.Mols[0].geomContainer.allCoords
        try :
            print "try changing useMasterDpyList"
            self.app().gui().viewer.useMasterDpyList=0
            print "ok ",self.app().gui().viewer.useMasterDpyList
        except :
            print "cant access viewer camera"

    def onRemoveObjectFromViewer(self, obj):
        if self.objectState.has_key(obj):
            self.objectState.pop(obj)    


#    def undoCmdBefore(self, nodes, only=False, negate=False, scaleFactor=0.05,
#                        cpkRad=0.0, quality=0, shrink=0.1, bScale=0.05,byproperty = False,
#                        propertyName=None, propertyLevel='Molecule',
#                        setScale=True, unitedRadii=True, redraw=True ):
#
#        #print "in undoCmdBefore: only= ", only, "negate=", negate, "scaleFactor=", scaleFactor, "cpkRad=", cpkRad, "quality=", quality, "property = ", propertyName, "propertyLevel = ", propertyLevel
#
#        kw = self.getLastUsedValues()
#        geomSet = []
#        for mol in self.app().Mols:
#            if not self.objectState.has_key(mol):
#                self.onAddObjectToViewer(mol)
#            geomSet = geomSet + mol.geomContainer.atoms['hpballs']
#        #print "geomset:", len(geomSet)
#        if len(geomSet) == 0:
#            # If nothing was displayed before, the undo of a display
#            # command is to undisplay what you just displayed 
#            kw['negate'] = True
#            #print "kw:", kw
#            return ( [(self, (nodes,), kw)], self.name)
#        else:
#            # If the current was already displayed and only is False and negate
#            # is False then the undo is to display this set using the last used
#            # values.
#            
#            # If something was already displayed, the undo of a display
#            # command is to display ONLY what was displayed before.
#            #kw['only'] = True
#            #print "kw:", kw
#            return ( [(self, (geomSet,), kw)], self.name)
#


    def doit(self, molSel, only=False, negate=False, scaleFactor=0.05,
             cpkRad=0.0, quality=0,  shrink=0.1, bScale=0.05, byproperty = False, propertyName=None,
             propertyLevel='Molecule', setScale=True, unitedRadii=True, redraw=True):

        #print "in  DOIT: ", "scaleFactor=", scaleFactor, "cpkRad=", cpkRad, 'bScale=', bScale, 'shrink=',shrink
#        scaleFactor, "cpkRad=", cpkRad, "quality=", quality, 
#        "shrink =" , shrink, "bScale ", bScale,"property=", propertyName, 
#        "propertyLevel = ", propertyLevel

        ##############################################################
        def drawAtoms( mol, sel, only, negate, scaleFactor, cpkRad, quality,
                       propertyName = None, propertyLevel="Molecule",
                       setScale=True):
            """
            handle all the atoms in one molecule
            """
            gc = mol.geomContainer
            ggeoms = gc.geoms            
            if setScale:
               setInds =  sel.getIndices()
               mol._cpkRad[setInds] = cpkRad
               mol._cpkScale[setInds] = scaleFactor                
#                atms.cpkScale = scaleFactor
#                atms.cpkRad = cpkRad
            _set = gc.atoms['hpballs']
            
            if len(sel) == len(mol._ag):
                if negate:
                    _set = mol.emptySelection()
                else:
                    _set = mol.select()
            else:
                ##if negate, remove current atms from displayed _set
                if negate:
                    _set = _set - sel

                else:
                    _set = _set | sel

            # at this point _set is what will still be displayed as CPK after this cmd
            gc.setAtomsForGeom('hpballs', _set)

            g = mol.geomContainer.geoms['hpballs']
            if len(_set)==0: # nothing is diplayed as CPK anymore for this molecule
                g.Set(visible=0, tagModified=False) # hide the geom
            else: # CPK is still displayed for this molecule
                #_set.sort()   # EXPENSIVE...
                # The following assumes that the lines geometry has been added:
#                colors = [x.colors.get('hpballs', x.colors['lines']) for x in _set]
#                # check if all colors are the same
#                cd = {}.fromkeys(['%f%f%f'%tuple(c) for c in colors])
#                if len(cd)==1:
#                    colors = [colors[0]]
                
                setInds =  _set.getIndices()
                cpkSF = np.array( mol._cpkScale)[setInds]
#                rad = _set.radius
                rad =  mol._ag.getRadii()[setInds]
#                if propertyName:
#                    propvals = self.getPropValues(_set, propertyName, propertyLevel)
#                    rad = propvals
#                cpkRad = numpy.array(_set.cpkRad)                
#                sphRad = cpkRad + cpkSF*numpy.array(rad)
#                sphRad = sphRad.tolist()
                #this is done in the shader, should be removed
                atmRad = np.array(rad)[setInds]#_set.radius)
                cpkRad = np.array( mol._cpkRad)[setInds]
                sphRad = cpkRad + cpkSF*atmRad
#                sphRad = cpkRad + atmRad#cpkSF
                sphRad = (sphRad/10.0).tolist()
                #setInds =  _set.getIndices()
#                bonds, atnobnd = _set.bonds
                bonds = _set.getBonds()
                #pdb.set_trace()
                #indices = map(lambda x: (x.atom1._bndIndex_,
                #                         x.atom2._bndIndex_), bonds)
                #coords is gc.allCoords[indices]
                #faces=[[x] for x in setInds]
                #pdb.set_trace()
                if g._firstDisplay: # use line colors for first display
                    g._firstDisplay = False
                    g.Set(vertices=gc.allCoords[setInds],bonds=bonds[1], inheritMaterial=False, 
                      materials=ggeoms['singleBonds'].materials[GL.GL_FRONT].prop[1], radii=sphRad, visible=1,
                      quality=quality,aScale=scaleFactor, shrink=shrink, 
                      bScale=bScale,tagModified=False)
                else :
                    g.Set(vertices=gc.allCoords[setInds],bonds=bonds[1], inheritMaterial=False,
                       radii=sphRad, visible=1,
                      quality=quality,aScale=scaleFactor, shrink=shrink, 
                      bScale=bScale,tagModified=False)
               # highlight selection
#                selMols, selAtms = self.app().getNodesByMolecule(self.app().activeSelection.get())
#                lMolSelectedAtmsDict = dict( zip( selMols, selAtms) )
#                lSelectedAtoms = lMolSelectedAtmsDict.get(mol, None)
#                if lSelectedAtoms is not None:
#                    lVertexClosestAtomSet = mol.geomContainer.atoms['hpballs']
#                    if len(lVertexClosestAtomSet) > 0:
#                        lVertexClosestAtomSetDict = dict(zip(lVertexClosestAtomSet,
#                                                             range(len(lVertexClosestAtomSet))))
#                        highlight = [0] * len(lVertexClosestAtomSet)
#                        for i in range(len(lSelectedAtoms)):
#                            lIndex = lVertexClosestAtomSetDict.get(lSelectedAtoms[i], None)
#                            if lIndex is not None:
#                                highlight[lIndex] = 1
#                        g.Set(highlight=highlight)
#            return setOn, setOff
        
        ##################################################################
        setOn = setOff = []        
        #molecules, atmSets = self.getNodes(nodes)
        mol = molSel.getAtomGroup().getMolecule()
        isInitialized = self.initializedFor.get(molSel._molecule(), False)
        if not isInitialized:
            self.initializeMolForCmd(molSel.getAtomGroup())        
#        molecules, atmSets = self.app().getNodesByMolecule(nodes, Atom)
#        setOn = AtomSet([])
#        setOff = AtomSet([])
#        for mol in molecules:
#            if not self.objectState.has_key(mol):
#                self.onAddObjectToViewer(mol)
#            hasRadii = False
#            try:
#                radii = mol.allAtoms.radius
#                hasRadii = True
#            except:
#                pass
#            #print "united:", unitedRadii, "mol.unitedRadii:", mol.unitedRadii
#            if not hasRadii or mol.unitedRadii != unitedRadii:
#                #print "defaultRadii"
#                mol.unitedRadii = unitedRadii
#                mol.defaultRadii(united=unitedRadii, overwrite=True)


        #try:
        #    radii = molecules.allAtoms.radius
        #except:
        #    self.app().assignAtomsRadii(molecules, 
        #                             overwrite=False, topCommand=False)
        if byproperty:
            drawAtoms(
                mol, molSel, only, negate, scaleFactor, cpkRad,
                quality, propertyName, propertyLevel,setScale=setScale)
        else:
            drawAtoms(mol, molSel, only, negate, 
                                 scaleFactor, cpkRad, quality,
                                 setScale=setScale)       
#        for mol, atms in map(None, molecules, atmSets):
#            if byproperty:
#                son, sof = drawAtoms(
#                    mol, atms, only, negate, scaleFactor, cpkRad,
#                    quality, propertyName, propertyLevel,setScale=setScale)
#            else:
#                son, sof = drawAtoms(mol, atms, only, negate, 
#                                     scaleFactor, cpkRad, quality,
#                                     setScale=setScale)
#            if son: setOn += son
#            if sof: setOff += sof
        self.app().pushUndoCmd(
            self.app().displayHyperBalls, (molSel,),
            {'negate': not negate, 'redraw':redraw,
             'scaleFactor':scaleFactor , 'cpkRad':cpkRad, 'quality':quality,
             'byproperty':byproperty, 'propertyName':propertyName,
             'propertyLevel':propertyLevel, 'setScale':setScale, 'unitedRadii':unitedRadii})
#        if only and len(molecules)!=len(self.app().Mols):
#            mols = self.app().Mols - molecules
#            for mol in mols:
#                only = False
#                negate = True
#                if byproperty:
#                    son, sof = drawAtoms(
#                        mol, mol.allAtoms, only, negate, scaleFactor, cpkRad,
#                        quality, propertyName, propertyLevel, setScale)
#                else:
#                    son, sof = drawAtoms(mol, mol.allAtoms, only, negate, 
#                                         scaleFactor, cpkRad, quality,
#                                         setScale=setScale)
#                if son: setOn += son
#                if sof: setOff += sof
                
        if self.createEvents:
            only=False
            event = EditGeomsEvent(
            'hpballs', [molSel,[only, negate,scaleFactor,cpkRad,quality,shrink, bScale,
                           byproperty,propertyName, propertyLevel,redraw]],
                                    setOn=setOn, setOff=setOff)
            self.app().eventHandler.dispatchEvent(event)

        
    def checkArguments(self, nodes, only=False, negate=False, 
                 scaleFactor=0.05,  cpkRad=0.0, quality=0,  shrink=0.1, bScale=0.05,
                 byproperty=False,
                 propertyName=None, propertyLevel='Molecule',
                 setScale=True, unitedRadii=True, redraw=True):
        """
           nodes --- TreeNodeSet holding the current selection
           \nonly --- Boolean flag: when True only the current selection will be
                    displayed
           \nnegate --- Boolean flag: when True undisplay the current selection
           \nscaleFactor --- value used to compute the radii of the sphere representing the atom (Sphere Radii = cpkRad + atom.radius * scaleFactor (default is 1.0)
           \ncpkRad --- value used to compute the radii of the sphere representing the atom (Sphere Radii = cpkRad + atom.radius * scaleFactor (default 0.0). This arg
           \ncorresponds to the 'Offset' on the gui input form.
           \nquality --- sphere quality default value is 10
           \npropertyName --- if specified,  the property is used to compute the sphere radii (Sphere Radii = cpkRad + property * scaleFactor),
           \npropertyLevel --- can be one of the following: 'Atom', 'Residue', 'Chain', 'Molecule'.
           \nsetScale --- when true atm.scaleFactor and atm.cpkRad are set"""

        #print "in  checkArguments: ", "scaleFactor=", scaleFactor, "cpkRad=", cpkRad, 'bScale=', bScale, 'shrink=',shrink
        kw = {}
        assert isinstance(scaleFactor, (int , float, long))
        assert isinstance(cpkRad, (int , float, long))
        assert isinstance(quality, int)
        assert isinstance(shrink, (int , float, long))
        assert isinstance(bScale, (int , float, long))

        assert (quality>=0)
        assert only in [True, False, 0, 1]
        assert negate in [True, False, 0, 1]
        assert byproperty in [True, False, 0, 1]
        assert setScale in [True, False, 0, 1]
        assert unitedRadii in [True, False, 0, 1]        
        assert propertyLevel in ['Atom', 'Residue', 'Chain', 'Molecule']
        if propertyName is not None:
            if not isinstance(propertyName, str):
                assert isinstance (propertyName, (list, tuple))
                propertyName = propertyName[0]
                assert isinstance(propertyName, str)
            byproperty = True
        if isinstance (nodes, str):
            self.nodeLogString = "'" + nodes +"'"
        #nodes = self.app().expandNodes(nodes)
        #assert nodes            
        kw['only'] = only
        kw['negate'] = negate
        kw['scaleFactor'] = scaleFactor
        kw['cpkRad'] = cpkRad
        kw['quality'] = quality
        kw['shrink'] = shrink
        kw['bScale'] = bScale
        kw['setScale']= setScale
        kw['unitedRadii'] = unitedRadii
        kw['propertyName'] = propertyName
        kw['propertyLevel']= propertyLevel
        kw['byproperty'] = byproperty

        #self.lastUsedValues['default']['byproperty'] = byproperty
        #self.lastUsedValues['default']['propertyName'] = propertyName
        #self.lastUsedValues['default']['unitedRadii'] = unitedRadii
        kw['redraw'] = redraw
        return (nodes,), kw


    def getPropValues(self, nodes, prop, propertyLevel=None):
        try:
            if propertyLevel is not None:
                #lNodesInLevel = self.app().activeSelection.get().findType(self.molDict[propertyLevel])
                lNodesInLevel = nodes.findType(self.molDict[propertyLevel])     
                num = len(lNodesInLevel.findType(Atom))
                propValues = getattr(lNodesInLevel, prop)
            else:
                propValues = getattr(nodes, prop)
        except Exception, error:
            from mglutil.errors import MGLError
            msg= "DisplayHyperBalls.getPropValues: nodes do not have the selected property %s" % prop
            MGLError(error, msg, raiseException=True)
            propValues=None
        return propValues


    def getProperties(self, selection):

        properties = [x[0] for x in selection[0].__dict__.items() if (
            x[0][0]!='_' and isinstance(x[1], (int, float)))]
        #print "getProperties:", properties
        return properties


        
class UndisplayHyperBalls(DisplayCommand):
    """ The undisplayCPK command is a picking command allowing the user to undisplay the CPK geometry representing the picked nodes when used as a picking command or the given nodes when called from the Python shell.
    \nPackage : PmvApp
    \nModule  : displayCmds
    \nClass   : UndisplayHyperBalls
    \nCommand : undisplayHyperBalls
    \nSynopsis:\n
        None <- undisplayHyperBalls(nodes, **kw)\n
        nodes : any set of MolKit2.Selection describing molecular components\n
        keywords: undisplay, HyperBalls\n
    """

    def onAddCmdToApp(self):
        DisplayCommand.onAddCmdToApp(self)
        if not self.app().commands.has_key('displayHyperBalls'):
            self.app().lazyLoad('displayHyperBallsCmds', commands=['displayHyperBalls'], package='pyHyperballs')


    def checkArguments(self, nodes, **kw):
        kw['negate' ] = True
        return self.app().displayHyperBalls.checkArguments(nodes, **kw)


    def doit(self, nodes, **kw):
        return self.app().displayHyperBalls(nodes, **kw)
        
commandClassFromName = {
    'displayHyperBalls': [DisplayHyperBalls, None],
    'undisplayHyperBalls' : [UndisplayHyperBalls,  None],
    }


def initModule(viewer, gui=True):
    for cmdName, values in commandClassFromName.items():
        cmdClass, guiInstance = values
        viewer.addCommand(cmdClass(), cmdName, guiInstance)
        
#to add in PmvApp
#pmv.lazyLoad('displayHyperBallsCmds', commands=[
#    'displayHyperBalls', 'undisplayHyperBalls'],
#     package='pyHyperballs')
#
