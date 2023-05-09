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
# Copyright: M. Sanner and TSRI 2015
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/ADFR/GUI/ADGridGUI.py,v 1.13 2016/12/07 00:38:32 sanner Exp $
#
# $Id: ADGridGUI.py,v 1.13 2016/12/07 00:38:32 sanner Exp $
#

import sys, weakref, os, numpy, string
from PySide import QtCore, QtGui

from ADFR.GUI import ICONPATH
from MolKit2.selection import Selection

class FlexResList(QtGui.QDialog):
    frstrSignal = QtCore.Signal(str)

    def __init__(self, parent=None):
        super(FlexResList, self).__init__(parent)
        self.setWindowTitle("select flexible resdiues")

        self.allButton = QtGui.QCheckBox("All")
        self.allButton.clicked.connect(self.selectAll)
        self.resListWidget = QtGui.QListWidget()
        self.resListWidget.itemClicked.connect(self.clickedRes)

    def fill(self, resList):
        self.resListWidget.clear()
        # order residues
        chids = numpy.unique([x.split(':')[0] for x in resList.keys()])
        chids.sort()
        num = 0
        for chid in chids:
            resnums = []
            names = []
            on = []
            for res, chk in resList.items():
                chidr, resstr = res.split(':')
                if chidr==chid:
                    resnums.append(int(resstr[3:]))
                    names.append( res )
                    on.append(chk)
            order = numpy.argsort(resnums)
            for n in order:
                item = QtGui.QListWidgetItem(names[n])
                if on[n]:
                    item.setCheckState(QtCore.Qt.Checked)
                else:
                    item.setCheckState(QtCore.Qt.Unchecked)
                self.resListWidget.addItem(item)
        
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept)
        
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.allButton)
        layout.addWidget(self.resListWidget)
        layout.addWidget(self.buttonBox)
        # Set dialog layout
        self.setLayout(layout)
        
    def selectAll(self):
        if self.allButton.isChecked():
            self.resListWidget.setDisabled(True)
        else:
            self.resListWidget.setDisabled(False)
            
    def clickedRes(self, item):
        d = {}
        #import pdb; pdb.set_trace()
        for i in range(self.resListWidget.count()):
            item = self.resListWidget.item(i)
            if item.checkState()==QtCore.Qt.Checked:
                chid, resstr = item.text().split(':')
                if not d.has_key(chid):
                    d[chid] = [resstr]
                else:
                    d[chid].append(resstr)
        chids = d.keys()
        chids.sort()
        s = ''
        for chid in chids:
            s += '%c:'%chid
            for res in d[chid]:
                s += '%s '%res
            s = s[:-1]
            s += ';'
        s = s[:-1]
        self.frstrSignal.emit(s)

class MyQLineEdit(QtGui.QLineEdit):
    recall = QtCore.Signal(str)

    def keyPressEvent(self, e):
        key = e.key()
        if key == QtCore.Qt.Key_Up:
            self.recall.emit('up')
        elif key == QtCore.Qt.Key_Down:
            self.recall.emit('down')
        else:
           QtGui.QLineEdit.keyPressEvent(self, e)

class ADGridGUI(QtGui.QWidget):

    cornerChangedSignal = QtCore.Signal()
    flexResChanged = QtCore.Signal()

    def __init__(self, app, spacing, smooth=0.5, parent=None, geom=None,
                 receptor=None, ligand=None):
        if geom:
            from DejaVu2.Box import NiceBox
            assert isinstance(geom, NiceBox)
        super(ADGridGUI, self).__init__(parent)
        self.app = weakref.ref(app)
        self.geom = geom
        self.receptor = receptor
        self.ligand = ligand
        self._flexRes = [] # list of prody residues for flexible residues
        self._flexResStr = '' # string in FlexRes entry widget
        self._flexResAtoms = None # MolKit2 selection with all atoms in flexible residues
        self._flexResSCAtoms = None # MolKit2 selection with all moving side chain atoms in flexible residues
        self._flexRecDescr = [] # [ (chid:[resnum1, resnum2]), ]
        self.cmdHistory = []
        self.historyPointer = 0
        self._ptsChanged = False
        self._szChanged = False
        self._spChanged = False
        self._initialized = False
        self._suspended = False
        self._baseSize = None
        self.LL = None
        self.UR = None
        
        self.mainLayout = QtGui.QVBoxLayout(self)
        self.groupBox = QtGui.QGroupBox()
        self.mainLayout.addWidget(self.groupBox)
        self.groupBox.setTitle(self.tr("docking box"))
        layout = QtGui.QVBoxLayout(self.groupBox)
        self.gridLayout = QtGui.QGridLayout()
        layout.addLayout(self.gridLayout)
        self.gridLayout.addWidget(QtGui.QLabel("x:"), 0, 1, 1, 1, QtCore.Qt.AlignCenter)
        self.gridLayout.addWidget(QtGui.QLabel("y:"), 0, 2, 1, 1, QtCore.Qt.AlignCenter)
        self.gridLayout.addWidget(QtGui.QLabel("z:"), 0, 3, 1, 1, QtCore.Qt.AlignCenter)
        self.gridLayout.addWidget(QtGui.QLabel("center:"), 1, 0, 1, 1, QtCore.Qt.AlignRight)
        self.gridLayout.addWidget(QtGui.QLabel("size (A):"), 2, 0, 1, 1, QtCore.Qt.AlignRight)
        self.gridLayout.addWidget(QtGui.QLabel("intervals:"), 3, 0, 1, 1, QtCore.Qt.AlignRight)
        self.gridLayout.addWidget(QtGui.QLabel("spacing:"), 4, 0, 1, 1, QtCore.Qt.AlignRight)
        self.gridLayout.addWidget(QtGui.QLabel("smoothing:"), 4, 2, 1, 1, QtCore.Qt.AlignRight)
        self.gridLayout.addWidget(QtGui.QLabel("cmd:"), 5, 0, 1, 1, QtCore.Qt.AlignRight)
        hlayout = QtGui.QHBoxLayout()
        layout.addLayout(hlayout)

        # addd button to set box to cover recetor or ligand with padding
        w = QtGui.QPushButton()
        w.clicked.connect(self.setGridFullReceptor)
        hlayout.addWidget(w)
        w.setText(" grid on rec.")
        w.setDisabled(True)
        self.gridForRec = w
        
        w = QtGui.QPushButton()
        w.clicked.connect(self.setGridFullLigand)
        hlayout.addWidget(w)
        w.setText(" grid on lig.")
        w.setDisabled(True)
        self.gridForLig = w

        w = QtGui.QLabel('padding:')
        w.setAlignment(QtCore.Qt.AlignRight)
        hlayout.addWidget(w)
        
        w = QtGui.QDoubleSpinBox()
        w.setDecimals(3)
        w.setValue(4.0)
        w.valueChanged.connect(self.paddingChanged)
        hlayout.addWidget(w)
        self.gridPaddingWidget = w
        
        self.centerSpinBoxes = []
        for i in range(3):
            widget = QtGui.QDoubleSpinBox()
            widget.valueChanged.connect(self.centerChanged)
            widget.setDecimals(3)
            widget.setMinimum(-9999.)
            self.centerSpinBoxes.append(widget)
            self.gridLayout.addWidget(widget, 1, i+1, 1, 1)
        
        self.sizefSpinBoxes = []
        for i in range(3):
            widget = QtGui.QDoubleSpinBox()
            widget.valueChanged.connect(self.sizeChanged)
            widget.setDecimals(3)
            widget.setMinimum(spacing)
            self.sizefSpinBoxes.append(widget)
            self.gridLayout.addWidget(widget, 2, i+1, 1, 1)
            
        self.pointsSpinBoxes = []
        for i in range(3):
            widget = QtGui.QSpinBox()
            widget.setRange(2, 2048)
            widget.valueChanged.connect(self.pointsChanged)
            self.pointsSpinBoxes.append(widget)
            self.gridLayout.addWidget(widget, 3, i+1, 1, 1)

        widget = QtGui.QDoubleSpinBox()
        widget.setDecimals(3)
        widget.setValue(spacing)
        widget.setMinimum(0.00001)
        widget.setSingleStep(0.025)
        widget.valueChanged.connect(self.spacingChanged)
        self.spacingWidget = widget
        self.gridLayout.addWidget(widget, 4, 1, 1, 1)
 
        widget = QtGui.QDoubleSpinBox()
        widget.setDecimals(3)
        widget.setValue(smooth)
        widget.setMinimum(0.00001)
        self.smoothWidget = widget
        self.gridLayout.addWidget(widget, 4, 3, 1, 1)
        
        self.cmdWidget = MyQLineEdit()
        self.cmdWidget.recall.connect(self.handleRecall)
        self.gridLayout.addWidget(self.cmdWidget, 5, 1, 1, 3)
        self.cmdWidget.returnPressed.connect(self.executeCmd)

        self.groupBox1 = QtGui.QGroupBox()
        self.mainLayout.addWidget(self.groupBox1)
        self.groupBox1.setTitle(self.tr("flexible receptor sidechains"))
        self.gridLayout1 = QtGui.QGridLayout(self.groupBox1)
        #self.gridLayout1.addWidget(QtGui.QLabel("flexRes:"), 0, 0, 1, 1)#, QtCore.Qt.AlignRight)
        self.flexResWidget = MyQLineEdit()
        self.gridLayout1.addWidget(self.flexResWidget, 0, 1, 1, 1)
        self.flexResWidget.returnPressed.connect(self.setFlexRes)
        self.flexResWidget.setDisabled(True)
        self.flexResWidget.setStyleSheet("QLineEdit{background: grey;}")

        self.setFlexRexButton = QtGui.QPushButton(QtGui.QIcon(os.path.join(ICONPATH, 'cube.png')),'')
        self.setFlexRexButton.clicked.connect(self.showFlexResChooser)
        
        #self.flexResButtonMenu = QtGui.QMenu()
        #w.setMenu(self.flexResButtonMenu)        
        self.gridLayout1.addWidget(self.setFlexRexButton, 0, 2, 1, 1)

        self._initialized = True
        if receptor is not None:
            self.setReceptor(receptor)

        if ligand is not None:
            self.setLigand(ligand)
 
        self.spacingChanged(spacing)
        self.ADatomTypes()

        def _flexResChanged(frstr):
            self.flexResWidget.setText(frstr)
            self.setFlexRes()
            
        self.flexResDialog = FlexResList()
        self.flexResDialog.frstrSignal.connect(_flexResChanged)

    def getFlexResDialogInput(self):
        resNames = []
        df = {}
        for res in self._flexRes:
            df['%c:%s%d'%(res.getChid(), res.getResname(), res.getResnum())] = True

        d = {}
        for atom in self.selectCAInBox():
            key = '%c:%s%d'%(atom.getChid(), atom.getResname(), atom.getResnum())
            if df.has_key(key):
                d[key] = True
            else:
                d[key] = False
        return d
    
    def showFlexResChooser(self):
        self.flexResDialog.fill(self.getFlexResDialogInput())
        #self.dialog = FlexResList(d)

        self.flexResDialog.show()
        #self.dialog.exec_()
        #dialog._raise()
        
    def selectCAInBox(self):
        if self.receptor is not None:
            serials = ' '.join([str(x) for x in self.selectAtomsInBox().getSerials()])
            selstr = 'same residue as (serial %s) and ca'%serials
            #selstr = 'ca and x>%f and x<%f and y>%f and y<%f and z>%f and z<%f'%(
            #    self.LL[0], self.UR[0], self.LL[1], self.UR[1], self.LL[2], self.UR[2])
            inBox = self.receptor.select(selstr)
            if inBox is not None:
                return inBox
            #    return self.removeNH(inBox)

    def isBoxOnReceptor(self):
        return not self.selectAtomsInBox() is None
    
    def selectAtomsInBox(self):
        if self.receptor is not None:
            selstr = '(not bb) and (not name CB) and (not resname PRO ALA GLY) and '
            selstr += 'x>%f and x<%f and y>%f and y<%f and z>%f and z<%f'%(
                self.LL[0], self.UR[0], self.LL[1], self.UR[1], self.LL[2], self.UR[2])
            inBox = self.receptor.select(selstr)
            if inBox is not None:
                return self.removeNH(inBox)
            
    def setGridFullReceptor(self, silent=False):
        self._suspended = True
        coords = self.receptor._ag.getCoords()
        mini = numpy.min(coords, 0)
        maxi = numpy.max(coords, 0)
        center = 0.5*(mini+maxi)
        lengths = (maxi-mini) + self.gridPaddingWidget.value()
        for i in range(3):
            self.centerSpinBoxes[i].setValue(center[i])
        spacing = self.spacingWidget.value()

        if spacing < 1.0 and not silent:
            msgBox = QtGui.QMessageBox(self)
            msgBox.setText("The grid spacing is currently %.3f.\n Large maps covering the entire receptor are usually coputed using larger spacing values, typically 1.0."%spacing)
            msgBox.setInformativeText("Do you want to use a spacing of 1.0 Angstroms?")
            msgBox.setStandardButtons(QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
            msgBox.setDefaultButton(QtGui.QMessageBox.Yes)
            ret =  msgBox.exec_()
            if ret == QtGui.QMessageBox.Yes:
                self.spacingWidget.setValue(1.0)

        for i in range(3):
            self.sizefSpinBoxes[i].setValue(lengths[i]+self.gridPaddingWidget.value())

        self._suspended = False
        self.spacingChanged(self.spacingWidget.value())
        self._baseSize = [x.value()-self.gridPaddingWidget.value() for x in self.sizefSpinBoxes]

    def setGridFullLigand(self):
        self._suspended = True
        coords = self.ligand._ag.getCoords()
        mini = numpy.min(coords, 0)
        maxi = numpy.max(coords, 0)
        center = 0.5*(mini+maxi)
        lengths = (maxi-mini) + self.gridPaddingWidget.value()
        for i in range(3):
            self.centerSpinBoxes[i].setValue(center[i])
        for i in range(3):
            if i ==2:
                self._suspended = False
            self.sizefSpinBoxes[i].setValue(lengths[i]+self.gridPaddingWidget.value())
        self._baseSize = [x.value()-self.gridPaddingWidget.value() for x in self.sizefSpinBoxes]
    
    def ADatomTypes(self):
        from ADFRcc import getFFParameters
        parameters = getFFParameters()
        self.ADAtomTypes = {}
        for i in range(parameters.numAtomTypes):
            self.ADAtomTypes[parameters.getAtomTypeByIndex(i).atomTypeName] = True
        #print len(self.ADAtomTypes), self.ADAtomTypes
        
    def setReceptor(self, mol):
        from MolKit2.molecule import Molecule
        assert isinstance(mol, Molecule)
        self.receptor = mol
        self.flexResWidget.setDisabled(False)
        self.gridForRec.setDisabled(False)
        self.flexResWidget.setStyleSheet("QLineEdit{background: white;}")
        
    def setLigand(self, mol):
        from MolKit2.molecule import Molecule
        assert isinstance(mol, Molecule)
        self.ligand = mol
        self.gridForLig.setDisabled(False)

    def computeCorners(self):
        cx, cy, cz = [x.value() for x in self.centerSpinBoxes]
        sx, sy, sz = [x.value() for x in self.sizefSpinBoxes]
        self.LL = [cx - (sx*.5), cy - (sy*.5), cz - (sz*.5)]
        self.UR = [cx + (sx*.5), cy + (sy*.5), cz + (sz*.5)]
        if self.receptor and not self._suspended:
            self.cornerChangedSignal.emit()
            if self.flexResDialog.isVisible():
                self.flexResDialog.fill(self.getFlexResDialogInput())
            
    def paddingChanged(self, value):
        if self._baseSize:
            self._suspended = True
            for i in range(3):
                self.sizefSpinBoxes[i].setValue( self._baseSize[i] + value )
            self._suspended = False
            self.computeCorners()
            
    def spacingChanged(self, value):
        #print 'SPACING'
        if value > 0 and self._initialized:
            self._suspended = True
            for i in range(3):
                self._ptsChanged = False
                self._szChanged = False
                self._spChanged = True
                self.sizefSpinBoxes[i].setSingleStep(value)
                self.sizefSpinBoxes[i].setMinimum(value)
                self.sizefSpinBoxes[i].setValue(self.pointsSpinBoxes[i].value() * value)
            self._suspended = False
            self.computeCorners()
        
    def pointsChanged(self, value):
        #print 'POINTS'
        if not self._initialized: return
        spacing = self.spacingWidget.value()
        if not self._szChanged and not self._spChanged:
            for i in range(3):
                self.sizefSpinBoxes[i].setValue( self.pointsSpinBoxes[i].value() * spacing )
        self.computeCorners()
        self._ptsChanged = True
        self._szChanged = False
        self._spChanged = False

    def sizeChanged(self, value=None):
        #print 'SIZE', not self._ptsChanged and not self._spChanged
        #import traceback
        #traceback.print_stack()
        if not self._initialized: return
        spacing = self.spacingWidget.value()
        if not self._ptsChanged and not self._spChanged:
            for i in range(3):
                self.pointsSpinBoxes[i].setValue( int(round(self.sizefSpinBoxes[i].value() / spacing) ))
        self.computeCorners()
        self._ptsChanged = False
        self._szChanged = True
        self._spChanged = False

        if self.geom:
            self.geom.setSides( self.sizefSpinBoxes[0].value(),
                                self.sizefSpinBoxes[1].value(),
                                self.sizefSpinBoxes[2].value() )
        
    def centerChanged(self, value=None):
        if not self._initialized: return
        self.computeCorners()
        if self.geom:
            self.geom.setCenter( self.centerSpinBoxes[0].value(),
                                 self.centerSpinBoxes[1].value(),
                                 self.centerSpinBoxes[2].value() )

    def handleRecall(self, what):
        if what=='up':
            if self.historyPointer > 0:
                self.historyPointer -= 1
                self.cmdWidget.clear()
                self.cmdWidget.insert(self.cmdHistory[self.historyPointer])
        elif what=='down':
            if self.historyPointer < len(self.cmdHistory):
                self.historyPointer += 1
                self.cmdWidget.clear()
                if self.historyPointer < len(self.cmdHistory):
                    self.cmdWidget.insert(self.cmdHistory[self.historyPointer])

    def _getIntOtFloat(self, s):
        try:
            value = int(s)
            return 'int', value
        except ValueError:
            try:
                value = float(s) 
                return 'float', value
            except ValueError:
                return None, None

    def removeNH(self, atoms):
        # remove H atoms attached to N atoms
        toRemove = []
        for atom in atoms:
            if atom.getElement()=='H':
                for natom in atom.iterBonded():
                    if natom.isbackbone:
                        toRemove.append(atom.getIndex())
        #for resnum in atoms.getResnums():
            #natoms = self.receptor.select( 'resnum %s and name N'%(resnum,))
            #if natoms:
                #for n in natoms:
                    #for atom in n.iterBonded():
                        #if atom.getElement()=='H':
                            #toRemove.append(atom.getIndex())
        atoms = atoms - Selection(self.receptor._ag, toRemove, '')
        return atoms
    
    def setFlexRes(self):
        ## expect string of type chid1:Resnum1,Resnum2,Resnum3 ;  chid2:Resnum1,Resnum2,Resnum3
        ## if ":" is missing chid will match all chains
        ## residues can also be specified using ResnameResnum and can beseparated by space or ,
        self._flexRes = []
        #self._flexResStr = [] # list of (
        self._flexResAtoms = None # MolKit2 selection with all atoms in flexible residues
        self._flexResSCAtoms = None # MolKit2 selection with all moving side chain atoms in flexible resid

        frstr = self.flexResWidget.text().encode('ascii', 'replace')
        if len(frstr)==0:
            self.flexResChanged.emit()
            return
        selStr = ""
        errors = []
        hv = self.receptor._ag.getHierView()
        #import pdb; pdb.set_trace()
        flexResAtoms = [] # list of atom indices for all atoms in FlexRes
        for expr in frstr.split(';'): #Loop over expressions delimited by ;
            dum = expr.split(':')
            if len(dum)==1:
                chids = numpy.unique(self.receptor._ag.getChids())
                residues = dum[0]
            else:
                chids = [dum[0]]
                residues = dum[1]
                selStr += "%s:"%dum[0]
            for res in residues.replace(',', ' ').split():
                if res[0] in string.digits: # number only
                    resname = None
                    resnum = int(res)
                else: # name and number
                    resname = res[:3]
                    try:
                        resnum = int(res[3:])
                    except ValueError:
                        msgBox = QtGui.QMessageBox(self)
                        msgBox.setText("ERROR: invalid syntax\n expect a residue number but got %s\nValid syntax examples:\n\t ILE10 VAL34\n\t 10,34\n\t A:10 34 ; B: 6 45\n"%res[3:])
                        msgBox.exec_()
                        self.flexResWidget.setText('')
                        return
                    
                for chid in chids:
                    res = hv.getResidue(chid, resnum)
                    if res is None:
                        errors.append('Residue %d not found in chain %s'%(
                            resnum, chid))
                    else:
                        flexResAtoms.extend(res.getIndices())
                        if resname:
                            if res.getResname()==resname:
                                self._flexRes.append(res)
                                #self._flexResStr.append( (res.getResname(), res.getResnum()) )
                                selStr += "%s%s,"%(resname,resnum)
                            else:
                                errors.append('Residue %d is %s which does not match specified %s'%(
                                    resnum, res.getResname(), resname))
                        else:
                            self._flexRes.append(res)
                            #self._flexResStr.append( (res.getResname(), res.getResnum()) )
                            selStr += "%s%s,"%(res.getResname(),resnum)
            selStr = selStr[:-1] # remove trailing ","
            selStr += ';'
        selStr = selStr[:-1] # remove trailing ";"
        self._flexResAtoms = Selection(self.receptor._ag, flexResAtoms, '')
        self._flexResSCAtoms = self.removeNH(self._flexResAtoms.select('not bb and not ca'))

        # check that grid covers flexRec atoms
        flexRecInside, outsideAtoms = self.atomInBox(self._flexResSCAtoms)
        if not flexRecInside:
            outSideResnames = ''
            resnums = outsideAtoms.getResnums()
            chids = outsideAtoms.getChids()
            done = {}
            outSideRes = []
            for chid, resnum in zip(chids, resnums):
                key = '%s:%d'%(chid, resnum)
                if done.has_key(key): continue
                done[key] = True
                res = hv.getResidue(chid, resnum)
                outSideRes.append(res)
                outSideResnames += '%s%d '%(res.getResname(), res.getResnum())
            msgBox = QtGui.QMessageBox(self)
            msgBox.setText('atoms for residues %s are outside the grid\nWould you like change the grid to cover these atoms ? If you choose NO this(ese) residue(s) will bre removed from the flexible residues list'%outSideResnames)
            msgBox.setStandardButtons(QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
            ret =  msgBox.exec_()

            if ret == QtGui.QMessageBox.Yes: # change grid to cover flexres atoms
                coords = self._flexResSCAtoms.getCoords()
                mini = numpy.min(coords, 0)
                maxi = numpy.max(coords, 0)
                center = [x.value() for x in self.centerSpinBoxes]
                size = [x.value() for x in self.sizefSpinBoxes]
                corner1 = center[0]-size[0]*.5, center[1]-size[1]*.5, center[2]-size[2]*.5
                corner2 = center[0]+size[0]*.5, center[1]+size[1]*.5, center[2]+size[2]*.5
                nc1 = numpy.min( (mini, corner1), 0) # new corner1
                nc2 = numpy.max( (maxi, corner2), 0) # new corne2
                nc = 0.5*(nc1 + nc2) # new center
                for i in range(3):
                    self.centerSpinBoxes[i].setValue(nc[i])
                    self.sizefSpinBoxes[i].setValue(nc2[i]-nc1[i]+self.gridPaddingWidget.value())
                    self._baseSize = [x.value() for x in self.sizefSpinBoxes]
            else: # remove offending residues from flexRes list
                for res in outSideRes:
                    self._flexRes.remove(res)
                selStr, descr = self.buildFlexRec(self._flexRes)
                #import pdb; pdb.set_trace()
                self.flexResWidget.setText(selStr)
                self.flexResDialog.fill(self.getFlexResDialogInput())
                self.setFlexRes()
        self.flexResChanged.emit()

        if len(errors):
            msgBox = QtGui.QMessageBox(self)
            msgBox.setText('\n'.join(errors))
            msgBox.setStandardButtons(QtGui.QMessageBox.Ok)
            ret =  msgBox.exec_()
            
        #print self._flexRes
        #print self._flexResStr
        #print selStr
        self.flexResWidget.setText(selStr)

    def buildFlexRec(self, resList):
        # use self._flexRes to build self._flexResStr and self._flexRecDescr
            chains = {}
            for res in resList:
                chid = res.getChid()
                if chains.has_key(chid):
                    chains[chid].append(res)
                else:
                    chains[chid] = [res]
            chids = chains.keys()
            chids.sort()
            selStr = ''
            descr = []
            for chid in chids:
                selStr+="%c:"%chid
                resl = []
                for res in chains[chid]:
                    resname = res.getResname()
                    resnum = res.getResnum()
                    resl.append((resname, resnum))
                    selStr+="%s%d,"%(resname, resnum)
                descr. append( (chid, resl) )
                selStr = selStr[:-1]
                selStr += ';'
            selStr = selStr[:-1]
            return selStr, descr
        
    def atomInBox(self, atoms):
        # check is atoms are in the grid box, If not identify atoms outside
        center = [x.value() for x in self.centerSpinBoxes]
        size = [x.value() for x in self.sizefSpinBoxes]
        p1x, p1y, p1z = center[0]-size[0]*.5, center[1]-size[1]*.5, center[2]-size[2]*.5
        p2x, p2y, p2z = center[0]+size[0]*.5, center[1]+size[1]*.5, center[2]+size[2]*.5
        coords = atoms.getCoords()
        mini = numpy.min(coords, 0)
        maxi = numpy.max(coords, 0)
        if mini[0]<p1x or mini[1]<p1y or mini[2]<p1z or \
            maxi[0]>p2x or maxi[1]>p2y or maxi[2]>p2z:
            outside = atoms.select('x<%f or x>%f or y<%f or y>%f or z<%f or z>%f'%(
                p1x, p2x, p1y, p2y, p1z, p2z))
            return False, outside
        return True, None
    
    def executeCmd(self):
        cmd = self.cmdWidget.text()
        self.cmdHistory.append(cmd)
        self.historyPointer += 1
        self.cmdWidget.clear()
        cmds = cmd.split(';')
        from time import time
        for cmd in cmds:
            t0 = time
            words = cmd.split()
            if words[0].lower()=='center':
                if len(words)==2 and words[1]=='flexres':
                    # FIXME implement this
                    center = [0,0,0] 
                    for i in range(3):
                        self.centerSpinBoxes[i].setValue(0.)

                elif len(words)==4:
                    for i in range(3):
                        self.centerSpinBoxes[i].setValue(float(words[i+1]))
                    self.centerChanged()
                else:
                    print 'Syntax Error: center x y z or center flexres'

            elif words[0].lower()=='size':
                size = None
                if len(words)==2:
                    vtype, size = self._getIntOtFloat(words[1])
                    if vtype=='int': # number of points
                        size = size*self.spacing 
                        self.size = [size, size, size]
                    elif vtype=='float': # size in angstroms
                        size = [size, size, size]                
                elif len(words)==4:
                    vtype, size = self._getIntOtFloat(words[1])
                    if vtype=='int': # number of points
                        nx = size
                        ny = int(words[2])
                        nz = int(words[3])
                        if 2*(nx/2)-nx==0: #event number
                            nx += 1
                        if 2*(ny/2)-ny==0: #event number
                            ny += 1
                        if 2*(nz/2)-nz==0: #event number
                            nz += 1
                        s = self.spacingWidget.value()
                        size = [nx*s, ny*s, nz*s]
                    elif vtype=='float': # size in angstroms
                        s = self.spacingWidget.value()
                        sx = size
                        sy = float(words[2])
                        sz = float(words[3])
                        nx = round(sx/s)
                        ny = round(sy/s)
                        nz = round(sz/s)
                        size = [nx*s, ny*s, nz*s]
                    else:
                        print 'Syntax Error: size d or size dx dy dz where d is int for grid points or float for angstroms'
                else:
                    print 'Syntax Error: size d or size dx dy dz where d is int for grid points or float for angstroms'
                if size is not None:
                    for i in range(3):
                        self.sizefSpinBoxes[i].setValue(size[i])
                    self._baseSize = [x.value() for x in self.sizefSpinBoxes]
                    self.sizeChanged()
            print 'Ran', cmd, 'in', time()-t0

           
if __name__=='__main__':
    import sys
    app = QtGui.QApplication(sys.argv)

    from DejaVu2.Qt.Viewer import Viewer
    window = QtGui.QWidget()
    layout = QtGui.QGridLayout(window)
    vi = Viewer()
    camera = vi.AddCamera()

    from DejaVu2.Box import NiceBox
    box = NiceBox('grid')
    box.addToViewer(vi)
    
    widget = ADGridGUI(vi, 0.375, geom=box)

    from MolKit2 import Read
    mol = Read('Astex/receptorsuniq-85/1n1m_rec.pdbqt')
    widget.setReceptor(mol)
    
    #widget.spinBoxSX.setValue(3)
    #widget.spinBoxSY.setValue(3)
    #widget.spinBoxSZ.setValue(3)
    
    layout.addWidget(widget, 0, 0, 1, 1)
    layout.addWidget(camera, 0, 1, 1, 1)
    
    #print widget.receptorContent.keys()
    #k,v = widget.receptorContent.items()[0]
    #print v.keys()
    window.resize(500, 400)
    window.show()
    app.exec_()
