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
# Date: 2000 Author: Ruth Huey Michel F. SANNER
#
#    sanner@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Michel Sanner and TSRI
#
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/PmvApp/GUI/Qt/coarseMolecularSurface.py,v 1.1.4.1 2017/07/13 20:46:08 annao Exp $
#
# $Id: coarseMolecularSurface.py,v 1.1.4.1 2017/07/13 20:46:08 annao Exp $
#

import sys
import string
import urllib
from PySide import QtCore, QtGui

surfNames = ['MSMSSurface', 'CoarseMS']
resNames = ['very smooth', 'medium smooth', 'smooth', 'custom']
isoValueNames = ['fast approximation', 'precise value', 'custom']

class Form(QtGui.QDialog):

    def __init__(self, parent=None):
        super(Form, self).__init__(parent)
        topLayout = QtGui.QVBoxLayout()
        formLayout = QtGui.QFormLayout()

        self.setWindowTitle("Compute Coarse Molecular Surface")
        #layout = QGridLayout(self) #immediate cbox;surflabel;
        #widgets:
        w = self.immediateCheckBox = QtGui.QCheckBox()
        w.toggled.connect(self.somethingChanged)
        formLayout.addRow(self.tr("&Immediate:"), self.immediateCheckBox)

        w = self.surfComboBox = QtGui.QComboBox(self)
        for name in surfNames:
            self.surfComboBox.addItem(name)      
        w.currentIndexChanged.connect(self.somethingChanged)
        formLayout.addRow(self.tr("&Surface name:"), self.surfComboBox)

        w = self.checkBox2 = QtGui.QCheckBox()
        w.toggled.connect(self.somethingChanged)
        formLayout.addRow(self.tr("&Per Molecule:"), self.checkBox2)

        w = self.gridSizeSpinBox = QtGui.QSpinBox()
        w.valueChanged.connect(self.somethingChanged)
        formLayout.addRow(self.tr("Grid size: "), self.gridSizeSpinBox)

        w = self.paddingSpinBox = QtGui.QSpinBox()
        w.valueChanged.connect(self.somethingChanged)
        formLayout.addRow(self.tr("padding: "), self.paddingSpinBox)

        w = self.surfResNamesComboBox = QtGui.QComboBox()
        for name in resNames:
            self.surfResNamesComboBox.addItem(name)
        w.currentIndexChanged.connect(self.somethingChanged)
        formLayout.addRow(self.tr("&Surface resolution: "), self.surfResNamesComboBox)

        w = self.surfResSpinBox = QtGui.QDoubleSpinBox()
        w.valueChanged.connect(self.somethingChanged)
        formLayout.addRow(self.surfResSpinBox)
        #formLayout.addRow(self.tr("&Surface resolution: "), self.surfResSpinBox)

        w = self.isoValueNameComboBox = QtGui.QComboBox()        
        for name in isoValueNames:
            self.isoValueNameComboBox.addItem(name)      
        w.currentIndexChanged.connect(self.somethingChanged)
        formLayout.addRow(self.tr("&Isocontour values:"), self.isoValueNameComboBox)

        w = self.checkBox3 = QtGui.QCheckBox()
        w.toggled.connect(self.somethingChanged)
        formLayout.addRow(self.tr("Bind Surface to molecule"), self.checkBox3)

        w = self.checkBox4 = QtGui.QCheckBox()
        w.toggled.connect(self.somethingChanged)
        formLayout.addRow(self.tr("Check surface components"), self.checkBox4)

        self.buttons = QtGui.QDialogButtonBox(
            QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel,
            QtCore.Qt.Horizontal, self)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)

        topLayout.addLayout(formLayout)
        topLayout.addWidget(self.buttons)
        self.setLayout(topLayout)

    def getPramsFromGui(self):
        #gridSize, bindGeom, isovalue, immediate, padding,
        #perMol, nodes, resolution, log
        gridSize = self.gridSizeSpinBox.value()
        bindGeom = self.checkBox3.isChecked()
        isovalue = isoValueNames[self.isoValueNameComboBox.currentIndex()]
        immediate = self.immediateCheckBox.isChecked()
        padding = self.paddingSpinBox.value()
        perMol = self.checkBox2.isChecked()
        surfResType = self.surfResNamesComboBox.currentText()
        surfRes = self.surfResSpinBox.value()
        surfResName = self.surfResNamesComboBox.currentText()
        bind_surface_to_molecule = self.checkBox3.isChecked()
        check_surface_components = self.checkBox4.isChecked()
        print 'run coarseMol with permol=%d'%perMol
        dict = {}
        dict['gridSize'] = 32
        dict['bindGeom'] = 1
        dict['immediate'] = False
        dict['padding'] = 0.0
        dict['perMol'] = 1
        dict['isovalue'] = 'fast approximation'
        #dict['nodes'] = 
        dict['resolution'] = -0.3
        dict['log'] = 1
        return {'perMol':perMol, 'immediate': immediate, 'gridSize':gridSize, 'bindGeom':bindGeom,'isovalue':isovalue,'padding':padding,}
    
    def issuePmvCommand(self, **args):
        print "issue pmv command here"
        for k,v  in args:
            print k, '--', v
        self.pmv.computeCoarseMolecularSurface(**args)
        
    def somethingChanged(self):
        if self.immediateCheckBox.isChecked():
            # we are in immediate mode we need to apply the change
            self.issuePmvCommand(**self.getPramsFromGui())
        else:
            # not immediate mode, do nothing
            return
        
    def accept(self):
        print 'you clicked accept', self.immediateCheckBox.isChecked()
        QtGui.QDialog.accept(self)
        
        ## surfNameMessage = QLabel("Select/type surface name:") 
        ## surfComboBox  = QComboBox(self) #@@
        ## for name in surfNames:
        ##     surfComboBox.addItem(name)      
        ## checkBox2 = QCheckBox("&Per Molecule") 

        ## gridSizeLabel = QLabel("Grid size: ")
        ## gridSizeSpinBox = QSpinBox()
        ## paddingLabel = QLabel("padding:")
        ## paddingSpinBox = QDoubleSpinBox()

        ## gridSizeSpinBox = QSpinBox()  
        ## paddingSpinBox = QDoubleSpinBox()
        ## surfResLabel = QLabel("Surface resolution:")
        ## surfResComboBox = QComboBox(self)
        ## for name in resNames:
        ##     ind = resNames.index(name)
        ##     surfResComboBox.insertItem(ind, name)
        ## surfResSpinBox = QDoubleSpinBox()
        ## isoContourLabel = QLabel("Isocontour values:")
        ## isoValComboBox = QComboBox(self)
        ## for name in isoValueNames:
        ##     ind = isoValueNames.index(name)
        ##     isoValComboBox.insertItem(ind+1, name)
        ## #isoValSpinBox = QDoubleSpinBox()
        ## checkBox3 = QCheckBox("Bind Surface to molecule")
        ## checkBox4 = QCheckBox("Check surface components")
        ## buttonBox = QDialogButtonBox(Qt.Horizontal)

        #buttonBox.addButton(QDialogButtonBox.Ok)
        #buttonBox.addButton(QDialogButtonBox.Cancel) 
        ## #buttonLayout = QHBoxLayout()
        ## #buttonLayout.addStretch()
        ## #buttonLayout.addWidget(buttonBox)
        ## #buttonLayout.addWidget(cancelButton)
        ## #layout2 = QGridLayout()
        ## #checkBox1 = QCheckBox("&Immediate")
        ## #layout.addWidget(checkBox1, 0,0)
        ## #nameLabel = QLabel("Select/type surface name:")
        ## #layout.addWidget(nameLabel, 1, 0)
        ## cb_layout = QGridLayout(self)   
        ## cb_layout.addWidget(checkBox1, 0, 0)        # checkBox1 'Immediate'
        ## cb_layout.addWidget(surfNameMessage, 1, 0)        # text 'Select/type surface name' 
        ## cb_layout.addWidget(surfComboBox, 2,0)      # QComboBox for surfNames
        ## cb_layout.addWidget(checkBox2, 3, 0)        # Per Molecule
        ## cb_layout.addWidget(gridSizeLabel, 4, 0)    # gridSizeLabel 
        ## cb_layout.addWidget(gridSizeSpinBox, 4, 1) # gridSizeSpinBox  
        ## cb_layout.addWidget(paddingLabel, 5, 0)     # paddingLabel 
        ## cb_layout.addWidget(paddingSpinBox, 5, 1)   # paddingDoubleSpinBox
        ## cb_layout.addWidget(surfResLabel, 6, 0)     # surfResLabel
        ## cb_layout.addWidget(surfResComboBox, 7,0)   # surfResComboBox
        ## cb_layout.addWidget(surfResSpinBox, 8,0)    # surfResSpinBox
        ## cb_layout.addWidget(isoContourLabel, 9, 0)  # isoContourLabel
        ## cb_layout.addWidget(isoValComboBox, 10, 0)  # isoValComboBox
        ## #cb_layout.addWidget(isoValSpinBox, 11, 0)   # isoValSpinBox
        ## cb_layout.addWidget(checkBox3, 12, 0)       # Bind Surface to molecule
        ## cb_layout.addWidget(checkBox4, 13, 0)       # Check surface components
        ## cb_layout.addWidget(buttonBox, 14, 0)        # Compute/Dismiss


        ## #layout.addLayout(buttonLayout, 2, 0,1,3)


        ## self.setLayout(cb_layout)

        #global_layout = QBoxLayout()
        #global_layout.addLayout(QBoxLayout())
        #widgets:
      

        #surfComboBox  = QComboBox(self) #@@
        #for name in surfNames:
        #    surfComboBox.addItem(name)
        #checkBox2 = QCheckBox("&Per Molecule") 

        #paddingLabel = QLabel("padding: ")
        ##paddingSpinBox = QDoubleSpinBox()
        #paddingLineEdit = QLineEdit("padding")
        #surfResLabel = QLabel("Surface resolution:")
        #surfResComboBox = QComboBox(self)
        #for name in resNames:
        #    ind = resNames.index(name)
        #    surfResComboBox.insertItem(ind, name)
        #surfResSpinBox = QDoubleSpinBox()
        #isoContourLabel = QLabel("Isocontour values:")
        #isoValComboBox = QComboBox(self)
        #for name in isoValueNames:
        #    ind = isoValueNames.index(name)
        #    isoValComboBox.insertItem(ind+1, name)
        #isoValSpinBox = QDoubleSpinBox()
        #checkBox3 = QCheckBox("Bind Surface to molecule")
        #checkBox4 = QCheckBox("Check surface components")
        #computeButton = QPushButton("&Compute")
        #dismissButton = QPushButton("&Dismiss")
        #buttonBox = QDialogButtonBox(Qt.Horizontal)
        #buttonBox.addButton(QDialogButtonBox.Ok)
        #buttonBox.addButton(QDialogButtonBox.Cancel)
        #buttonBox = QDialogButtonBox(QDialogButtonBox.Apply|QDialogButtonBox.Cancel)
 
        #cb_layout = QGridLayout(self)   
        #cb_layout.addWidget(checkBox1, 0, 0)        # checkBox1 'Immediate'
        #cb_layout.addWidget(nameLabel, 1, 0)        # text 'Select/type surface name' 
        #cb_layout.addWidget(surfComboBox, 2,0)      # QComboBox for surfNames
        #cb_layout.addWidget(checkBox2, 3, 0)        # Per Molecule
        #cb_layout.addWidget(gridSizeLabel, 4, 0)    # gridSizeLabel 
        #cb_layout.addWidget(gridSizeLineEdit, 4, 1) # gridSizeLineEdit  
        #cb_layout.addWidget(paddingLabel, 5, 0)     # paddingLabel 
        #cb_layout.addWidget(paddingLineEdit, 5, 1)   # paddingSpinBox
        #cb_layout.addWidget(surfResLabel, 6, 0)     # surfResLabel
        #cb_layout.addWidget(surfResComboBox, 7,0)   # surfResComboBox
        #cb_layout.addWidget(surfResSpinBox, 8,0)    # surfResSpinBox
        #cb_layout.addWidget(isoContourLabel, 9, 0)  # isoContourLabel
        #cb_layout.addWidget(isoValComboBox, 10, 0)  # isoValComboBox
        #cb_layout.addWidget(isoValSpinBox, 11, 0)   # isoValSpinBox
        #cb_layout.addWidget(checkBox3, 12, 0)       # Bind Surface to molecule
        #cb_layout.addWidget(checkBox4, 13, 0)       # Check surface components
        #cb_layout.addWidget(buttonBox, 14, 0)        # Compute/Dismiss

        #cb_layout.addWidget(self.checkBox1, 0,0)
        #layout = QGridLayout()
        
        #nameText = "Select/type surface name:"
        #self.nameLabel.setText(nameText)
        #layout.addWidget(self.checkBox1, 1,1)
        #layout.addWidget(self.nameLabel,2,1)
        #self.setWindowTitle("Coarse MSMS Surface")

if __name__== '__main__':
    app = QtGui.QApplication(sys.argv)
    form = Form()
    # non modal dialog
    form.show()
    form.raise_()
    form.activateWindow()
    #modal dialog
    #value = form.exec_()
    app.exec_()

