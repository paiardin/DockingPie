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
# Author: Michel F. SANNER, , Anna Omelchenko
#
# Copyright: M. Sanner TSRI 2014
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/PmvApp/editCmds.py,v 1.2.4.2 2017/10/12 00:26:18 annao Exp $
#
# $Id: editCmds.py,v 1.2.4.2 2017/10/12 00:26:18 annao Exp $
#
#

from PmvApp.Pmv import MVCommand, AfterAddMoleculeEvent, RefreshDisplayEvent

from MolKit2.selection import Selection, SelectionSet

class AssignAtomsRadiiCommand(MVCommand):
    """This command adds radii to all atoms loaded in the application. Only default radii for now.
    \nPackage:PmvApp
    \nModule :editCmds
    \nClass:AssignAtomsRadiiCommand
    \nCommand:assignAtomsRadii
    \nSynopsis:\n
    None <- mv.assignAtomsRadii(nodes, united=1, overwrite=1,**kw)\n
    \nRequired Arguments:\n
    nodes --- TreeNodeSet holding the current selection
    \nOptional Arguments:\n
    \nunited   --- (default=1) flag to specify whether or not to consider
            hydrogen atoms. When hydrogen are there the atom radii
            is smaller. 

    overwrite ---(default=1) flag to specify whether or not to overwrite
            existing radii information.\n

    
    """

    def __init__(self):
        MVCommand.__init__(self)
        #self.flag = self.flag | self.objArgOnly


    def doit(self, nodes, united=True, overwrite=False):
        #nodes = self.vf.expandNodes(nodes)
        if not nodes: return
        for sel in nodes:
            mol = sel.getAtomGroup().getMolecule()
            # Reassign the radii if overwrite is True
            if overwrite is True:
                mol.unitedRadii = united
                mol.defaultRadii(united=united, overwrite=overwrite)
            # Reassign the radii if different.
            elif mol.unitedRadii != united:
                mol.unitedRadii = united
                mol.defaultRadii(united=united, overwrite=overwrite)
        

    def checkArguments(self, nodes, united=True, overwrite=False, **kw):
        """ None <- mv.assignAtomsRadii(nodes, united=True, overwrite=False,**kw)
        \nRequired Arguments:\n
         nodes --- TreeNodeSet holding the current selection 

         \nOptional Arguments:\n
          united --- (default=True) Boolean flag to specify whether or not
                    to consider hydrogen atoms. When hydrogen are there the
                    atom radii is smaller.\n 

          overwrite --- (default=True) Boolean flag to specify whether or not to overwrite
                    existing radii information.\n

        """
        if isinstance(nodes, str):
            self.nodeLogString = "'"+nodes+"'"
        nodes = self.app().expandNodes(nodes)
        if not len(nodes): return (), {}
        kw = {}
        kw['united'] = united
        kw['overwrite'] = overwrite
        return (nodes, ), kw


class TypeAtomsAndBonds(MVCommand):
    """This command assigns atom types and bond order. after this command is run on a molecule
    each atom will have atomType, hbType and the atom Group will have _bondOrders. Note the Hydrogen
    bond might be missing if the molecule is not protonated.
    \nPackage:PmvApp
    \nModule :editCmds
    \nClass:TypeAtomsAndBonds
    \nCommand:typeAtomsAndBonds
    \nSynopsis:\n
    None <- mv.typeAtomsAndBonds(nodes,**kw)\n
    \nRequired Arguments:\n
    nodes --- TreeNodeSet holding the current selection
    """

    def __init__(self):
        MVCommand.__init__(self)
        #self.flag = self.flag | self.objArgOnly

    def doit(self, selection, redraw=True):
        mol = selection.getAtomGroup().getMolecule()
        mol.typeAtomsAndBonds()        
        event = RefreshDisplayEvent(molecule=mol)
        self.app().eventHandler.dispatchEvent(event)

    def checkArguments(self, nodes, redraw=True):
        """ None <- mv.assignAtomsRadii(nodes, united=True, overwrite=False,**kw)
        \nRequired Arguments:\n
         nodes --- TreeNodeSet holding the current selection 
        """
        #if isinstance(nodes, str):
        #    self.nodeLogString = "'"+nodes+"'"
        #nodes = self.app().expandNodes(nodes)
        kw = {}
        kw['redraw'] = redraw
        assert isinstance(nodes, SelectionSet)
        return (nodes, ), kw

from PySide import QtCore, QtGui
from time import time, sleep
import openbabel as ob
import sys, os, numpy
import pybel # else the plugins dont; work :(

class MinimizerProgress(object):

    def __init__(self, mol, step, energy, coords, status):
        self._time = time()
        self.mol = mol
        self.coords = coords
        self.energy = energy
        self.status = status
        self.step = step
        
class MySignal(QtCore.QObject):
        sig = QtCore.Signal(MinimizerProgress)

## class Minimizer(QtCore.QThread):

##     #progress = QtCore.Signal(MinimizerProgress) # create a custom signal we can subscribe

##     def __init__(self, pmol, nsteps, upFreq=10):
##         # calling super is needed to avoid
##         # AttributeError: 'PySide.QtCore.pyqtSignal' object has no attribute 'connect'
##         super(Minimizer, self).__init__()
##         self.mol = pmol
##         self.nsteps = nsteps
##         self.upFreq = upFreq
##         from MolKit2.openBabelInterface import ProdyToOBMol, OBMolToPrody
##         self.obmol2 = ProdyToOBMol(pmol.select())
##         forcefield="mmff94"
##         self.ff = pybel._forcefields[forcefield]
##         success = self.ff.Setup(self.obmol2)
##         if success is False:
##             raise ValueError, 'ffSetup failed'
##         print 'E:', self.ff.Energy(0)
##         self.stepSignal = MySignal()
##         self.endSignal = MySignal()
        
##     def run(self):
##         #obj = MinimizerProgress(self.mol, -2, self.ff.Energy(0), None, None)
##         #self.progress.emit(obj)
##         #self.signal.sig.emit('after init')
##         #self.progress.emit('before init')
##         self.ff.SteepestDescentInitialize(self.nsteps,1e-5)
##         #obj = MinimizerProgress(self.mol, -1, self.ff.Energy(0), None, None)
##         #self.progress.emit(obj)
##         #self.signal.sig.emit('after init')
##         nbs = 0
##         while nbs < self.nsteps:
##             status = self.ff.SteepestDescentTakeNSteps(self.upFreq)
##             self.ff.GetCoordinates(self.obmol2)
##             coords = numpy.zeros( (self.obmol2.NumAtoms(), 3), 'f')
##             sys.stdout.write('step - %d\n'%nbs)
##             sys.stdout.flush()
##             for j, atom in enumerate(ob.OBMolAtomIter(self.obmol2)):
##                 n = self.obmol2.obToProdyIndex[j]
##                 coords[n][0] = atom.GetX()
##                 coords[n][1] = atom.GetY()
##                 coords[n][2] = atom.GetZ()
##             nbs += self.upFreq
##             self.stepSignal.sig.emit('step %d'%nbs)

##             #obj = MinimizerProgress(self.mol, nbs, self.ff.Energy(0),
##             #                        coords, status)
##             #self.progress.emit(obj)
##         self.endSignal.sig.emit('end')
##         #print 'MINIMIZER THREAD dieing'
##         #self.terminate() # cause FATAL python error and Abort

class SpawnOBmin(QtCore.QThread):

    def __init__(self, pmol, nsteps, upFreq=10):
        # calling super is needed to avoid
        # AttributeError: 'PySide.QtCore.pyqtSignal' object has no attribute 'connect'
        super(SpawnOBmin, self).__init__()
        self.mol = pmol
        self.nsteps = nsteps
        self.upFreq = upFreq
        self.obminsh = os.path.join(os.path.split(sys.executable)[0],
                                     'pythonsh') + ' obmin.py'
        from MolKit2.openBabelInterface import ProdyToOBMol, OBMolToPrody
        self.obmol = ProdyToOBMol(pmol.select())
        self.stepSignal = MySignal()
        self.endSignal = MySignal()

    def run(self):
        import platform, subprocess

        try:
            from Queue import Queue, Empty
        except ImportError:
            from queue import Queue, Empty  # python 3.x

        from threading  import Thread

        if platform.uname()[0] == 'Windows':
            _shell=False
        else:
            _shell=True

        proc = subprocess.Popen(self.obminsh,
                                stdin=subprocess.PIPE , 
                                stdout=subprocess.PIPE , 
                                stderr=subprocess.PIPE, 
                                bufsize=-1, shell=_shell)

        # write mini params
        proc.stdin.write("minimizer steepest\n")
        proc.stdin.write("nbSteps %d\n"%self.nsteps)
        proc.stdin.write("upFreq %d\n"%self.upFreq)
        proc.stdin.write("EndParameters\n")
        
        # now write molecule
        obconv = ob.OBConversion()
        obconv.SetOutFormat("mol2")
        molStr = obconv.WriteString(self.obmol)
        proc.stdin.write(molStr)

        # close stdin to for parsing
        proc.stdin.close()

        #print proc.stdout.readlines()
        #print proc.stderr.readlines()
        
        ## while True:
        ##     if proc.poll() == None:
        ##         #print 'Reading obmin.stdout'
        ##         line = proc.stdout.readline()
        ##         print line
        ##     else:
        ##         #self.endSignal.sig.emit('end')
        ##         print proc.stderr.readlines()
        ##         print 'done'
        ##         break
        
        length = len(self.mol._ag)        
        coords = numpy.zeros( (length,3), 'f')
        _c = numpy.zeros( (3*length,), 'f')
        while True:
            if proc.poll() == None:
                #print 'Reading obmin.stdout'
                line = proc.stdout.readline()
                if len(line)==0: continue
                w = line.split()
                if w[0]=='STEP':
                    nbs = int(w[1])
                    ene = float(w[3])
                elif w[0]=='COORDS':
                    _c[:] = w[1:]
                    j = 0
                    for i in range(length):
                        n = self.obmol.obToProdyIndex[i]
                        coords[n][0] = _c[j]
                        coords[n][1] = _c[j+1]
                        coords[n][2] = _c[j+2]
                        j += 3
                    obj = MinimizerProgress(self.mol, nbs, ene, coords, True)
                    self.stepSignal.sig.emit(obj)
            else:
                #self.endSignal.sig.emit('end')
                break
       
class MinimizeMolecule(MVCommand):
    """This command

    \nPackage:PmvApp
    \nModule :editCmds
    \nClass:MinimizeMolecule
    \nCommand:minimize
    \nSynopsis:\n
    None <- mv.minimize(nodes,**kw)\n
    \nRequired Arguments:\n
    nodes --- TreeNodeSet
    """

    def __init__(self):
        MVCommand.__init__(self)
        #self.flag = self.flag | self.objArgOnly

    def endcb(self, ev):
        print 'FOGO', ev
        return
    
    def stepcb(self, ev):
        #print 'FOGO', ev
        #return
        if ev.coords is not None:
            print 'CB Step: %5d, Energy: %f'%(ev.step, ev.energy)
            ev.mol.geomContainer.allCoords[:] = ev.coords
            ev.mol._ag._setCoords(ev.coords)
            self.app().displayLines(ev.mol)
            self.app().gui().viewer.OneRedraw()

    def doit(self, selection, redraw=True):
        mol = selection.getAtomGroup().getMolecule()
        #self.minimizer = Minimizer(mol, 100, 10)
        self.minimizer = SpawnOBmin(mol, 10000, 10)
        self.minimizer.stepSignal.sig.connect(self.stepcb)#, QtCore.Qt.QueuedConnection)
        self.minimizer.endSignal.sig.connect(self.endcb)#, QtCore.Qt.QueuedConnection)
        if not self.minimizer.isRunning():
            self.minimizer.start()

    def checkArguments(self, nodes, redraw=True):
        """ None <- mv.assignAtomsRadii(nodes, united=True, overwrite=False,**kw)
        \nRequired Arguments:\n
         nodes --- TreeNodeSet holding the current selection 
        """
        #if isinstance(nodes, str):
        #    self.nodeLogString = "'"+nodes+"'"
        #nodes = self.app().expandNodes(nodes)
        kw = {}
        kw['redraw'] = redraw
        assert isinstance(nodes, SelectionSet)
        return (nodes, ), kw

        
class BiologicalUnit(MVCommand):
    """This command

    \nPackage:PmvApp
    \nModule :editCmds
    \nClass:BiologicalUnit
    \nCommand:biologicalUnit
    \nSynopsis:\n
    MoleculeSet <- mv.BiologicalUnit(nodes,**kw)\n
    \nRequired Arguments:\n
    nodes --- SelectionSet
    """

    def doit(self, selection, redraw=True):
        mol = selection.getAtomGroup().getMolecule()
        bioUnits = mol.biologicalUnit()
        app = self.app()
        for bm in bioUnits:
            try:
                app.addMolecule(bm)
                app.applyDefaultCommands(bm, "Molecule")
                event = AfterAddMoleculeEvent(molecule=bm)
                app.eventHandler.dispatchEvent(event)
            except:
                ## if app.trapExceptions is False:
                ##     exc_info = sys.exc_info()
                ##     raise exc_info[1], None, exc_info[2]
                ## else:
                msg = 'Error adding biological unit %s'%(bm.name)
                app.errorMsg(sys.exc_info(), msg, obj=None)
                app._executionReport.finalize()
        return bioUnits

    def checkArguments(self, nodes, redraw=True):
        """ None <- mv.assignAtomsRadii(nodes, united=True, overwrite=False,**kw)
        \nRequired Arguments:\n
         nodes --- TreeNodeSet holding the current selection 
        """
        #if isinstance(nodes, str):
        #    self.nodeLogString = "'"+nodes+"'"
        #nodes = self.app().expandNodes(nodes)
        kw = {}
        kw['redraw'] = redraw
        assert isinstance(nodes, SelectionSet)
        return (nodes, ), kw


commandClassFromName = {
    'assignAtomsRadii' : [AssignAtomsRadiiCommand,  None],
    'typeAtomsAndBonds' : [TypeAtomsAndBonds,  None],
    'minimize' : [MinimizeMolecule,  None],
    'biologicalUnit' : [BiologicalUnit,  None],
    }


def initModule(viewer, gui=True):
    for cmdName, values in commandClassFromName.items():
        cmdClass, guiInstance = values
        viewer.addCommand(cmdClass(), cmdName, guiInstance)
