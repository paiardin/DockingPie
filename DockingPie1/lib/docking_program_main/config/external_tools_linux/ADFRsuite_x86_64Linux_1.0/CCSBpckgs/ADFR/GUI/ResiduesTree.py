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
# Copyright: M. Sanner and TSRI 2016
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/ADFR/GUI/ResiduesTree.py,v 1.16.2.2 2017/08/02 00:27:34 annao Exp $
#
# $Id: ResiduesTree.py,v 1.16.2.2 2017/08/02 00:27:34 annao Exp $
#
from PySide import QtCore, QtGui
from ADFR.GUI import ICONPATH
from prody.measure.contacts import findNeighbors
import os, numpy

from MolKit2.selection import Selection

class ResiduesTree(QtGui.QWidget):

    resChanged = QtCore.Signal()

    def __init__(self, rec, lig=None, checkedItems=None, parent=None):
        # checkedItems is a a dict where the key is 'chid:resnum' and the
        # value is True or False, only residues with a key in checked are shown
        super(ResiduesTree, self).__init__(parent)
        self._suspend = False
        self._resToItem = {}
        self._resCA = []
        self._CAindices = set([])
        self.receptor = rec
        self.ligand = None
        self.checkedItems = checkedItems
        self.buildUI()
        if lig:
            self.setLigand(lig)

    def onCheckGroup(self):
        if self.groupBox1.isChecked():
            self.selectClose()
        
    def buildUI(self):
        self.groupBox1 = QtGui.QGroupBox(self)
        self.groupBox1.setTitle(self.tr("residues within"))
        self.groupBox1.setCheckable(True)
        self.groupBox1.setChecked(False)
        self.groupBox1.clicked.connect(self.onCheckGroup)
        self.distanceWidget = w = QtGui.QDoubleSpinBox()
        w.setToolTip("select residues with atoms within this distance of ligand atoms")
        w.setMinimum(0.1)
        w.setSingleStep(0.2)
        w.setDecimals(2)
        w.setValue(3.0)
        self.groupBox1.setDisabled(True)
        w.valueChanged.connect(self.selectClose)

        self.buildTree()

        layout = QtGui.QVBoxLayout()
        layout1 = QtGui.QHBoxLayout()
        self.groupBox1.setLayout(layout1)
        layout.addWidget(self.groupBox1)
        layout1.addWidget(self.distanceWidget)
        layout1.addWidget(QtGui.QLabel(' of ligand'))
        
        layout.addWidget(self.treeWidget)
        self.setLayout(layout)
        #treeWidget.itemClicked is called when we click on box or label
        # treeWidget.itemClicked.connect(self.onClick)

        # treeWidget.itemChanged is only called when box toggles
        self.treeWidget.itemChanged.connect(self.onClick)

    def buildTree(self):
        self.treeWidget = treeWidget = QtGui.QTreeWidget()
        hv = self.receptor.getHierView()
        chIcon = QtGui.QIcon(os.path.join(ICONPATH, 'chain.png'))
        resIcon = QtGui.QIcon(os.path.join(ICONPATH, 'residue.png'))
        _CAindices = []
        for chain in hv.iterChains():
            chid = chain.getChid()
            chItem = QtGui.QTreeWidgetItem(treeWidget.invisibleRootItem(), chid)
            chItem.setExpanded(True)
            chItem.setIcon(0,chIcon)
            
            for res in chain.iterResidues():
                resnum = res.getResnum()
                key1 = '%s:%d'%(chid,resnum)
                if self.checkedItems and key1 not in self.checkedItems:
                    continue
                resname = res.getResname()
                resItem = QtGui.QTreeWidgetItem(chItem )
                resItem.setText(0, '%s%d'%(resname,resnum ))
                resItem.setIcon(0, resIcon)
                resItem._pmvObj = res
                if self.checkedItems:
                    if self.checkedItems['%s:%d'%(chid,resnum)]:
                        #sel = res.select('ca')
                        #if sel is not None:
                        sel = res.select("chain %s resnum %d"%(chid,resnum))
                        #self._CAindices.append(sel.getIndices()[0])
                        _CAindices.extend(sel.getIndices())
                        resItem.setCheckState(0, QtCore.Qt.Checked)
                    else:
                        resItem.setCheckState(0, QtCore.Qt.Unchecked)
                else:
                    resItem.setCheckState(0, QtCore.Qt.Unchecked)
                key = '%s:%s%d'%(chid,resname,resnum)
                self._CAindices = set(_CAindices)
                self._resToItem[key] = resItem

    def selectClose(self):
        #if self.ligand is None:
        #    return

        atoms = findNeighbors(self.ligand._ag,
                              self.distanceWidget.value(),
                              self.receptor)
        #import pdb; pdb.set_trace()
        if len(atoms):
            self._suspend = True
            indices = numpy.unique([x[1].getSerial() for x in atoms])
            resnums = self.receptor._ag._data['resnum']
            chids = self.receptor._ag._data['chain']
            if self.checkedItems:
                indicesInBox = []
                for ind in indices:
                    if self.checkedItems.has_key('%s:%d'%(chids[ind],
                                                          resnums[ind])):
                        indicesInBox.append(ind)
                indices = indicesInBox
            selStr = '(same residue as serial %s) and ca'%' '.join([str(x) for x in indices])
            self._resCA = resCA = self.receptor.select(selStr)
            #for atom in resCA:
            #    print '%s%d '%(atom.getResname(), atom.getResnum())
            # uncheck all
            for item in self._resToItem.values():
                item.setCheckState(0, QtCore.Qt.Unchecked)
            # check close
            for i, ca in enumerate(resCA):
                key = '%s:%s%d'%(ca.getChid(), ca.getResname(), ca.getResnum())
                self._resToItem[key].setCheckState(0, QtCore.Qt.Checked)
                if i==0:
                    # make first check residue visible
                    self.treeWidget.scrollToItem(self._resToItem[key])
            #print len(atoms), len(resCA), time()-t0
            #self._CAindices = resCA.getIndices().tolist()
            self._CAindices = set(resCA.getIndices())
            self._suspend = False
        else:
            self._resCA = []
            self._CAindices = set([])
        self.resChanged.emit()

    def getAtoms(self):
        if len(self._CAindices)==0:
            return None
        else:
            return self.receptor.select(
                'same residue as index %s'%' '.join(
                    [str(x) for x in self._CAindices]))
        
    def setLigand(self, ligand):
        self.ligand = ligand
        self.groupBox1.setDisabled(False)

    def onClick(self, item, column):
        if self._suspend: return
        res = item._pmvObj
        if item.checkState(column) == QtCore.Qt.Checked:
            #self._CAindices.append(ca.getIndices()[0])
            self._CAindices.update(res.getIndices())
        else:
            #self._CAindices.remove(ca.getIndices()[0])
            self._CAindices = self._CAindices - set(res.getIndices())
        self.resChanged.emit()

class ResiduesTreeDialog(QtGui.QDialog):
    closedSignal = QtCore.Signal()

    def __init__(self, rec, lig=None, checkedItems=None, title='No Name',
                 parent=None):
        super(ResiduesTreeDialog, self).__init__(parent)
        self.setWindowTitle(title)
        self.resTreeWidget = ResiduesTree(rec, lig, checkedItems, parent)
        self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        self.buttonBox = QtGui.QDialogButtonBox()
        self.buttonBox.accepted.connect(self.accept)

        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.resTreeWidget)
        layout.addWidget(self.buttonBox)
        # Set dialog layout
        self.setLayout(layout)
        #self.setWindowFlags( self.windowFlags() & ~QtCore.Qt.WindowCloseButtonHint)

    def reject(self):
        QtGui.QDialog.reject(self)
        self.closeEvent()
        
    def closeEvent(self, evnt=None):
        self.closedSignal.emit()
        
if __name__=='__main__':
    import sys

    from MolKit2 import Read
    rec = Read('4EK3_rec.pdbqt')
    lig = Read('4EK4_lig.pdbqt')
    app = QtGui.QApplication(sys.argv)
    from time import time
    t0 = time()
    widget = ResiduesTree(rec, lig)
    #widget.setLigand(lig)
    #def p(chid, resname, checked):
    #    print chid, resname, checked
    #widget.resClickedSignal.connect(p)
    
    #print time()-t0
    widget.resize(150,300)
    widget.show()
    sys.exit(app.exec_())
