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

from PySide import QtCore, QtGui
from mglutil.util.callback import CallbackFunction

class ReduceParamsGUI(QtGui.QWidget):
    def __init__(self, name, parent=None):
        super(ReduceParamsGUI, self).__init__(parent)

        self.mainLayout = gLayout = QtGui.QGridLayout()
        self.flipWidget = QtGui.QCheckBox("Flip")
        self.mainLayout.addWidget(self.flipWidget, 0, 0, 1, 1)
        self.flipWidget.stateChanged.connect(self.changeFlip)
        self.setLayout(gLayout)
        
    def changeFlip(self, value):
        print 'flip', value

    def getParams(self):
        return {'flip':self.flipWidget.isChecked()}
        
class OBmolAddHParamsGUI(QtGui.QWidget):
    def __init__(self, name, parent=None):
        super(OBmolAddHParamsGUI, self).__init__(parent)

        self.mainLayout = gLayout = QtGui.QGridLayout()
        self.phCheckBox = QtGui.QCheckBox("pH 0-14")
        self.phCheckBox.setChecked(True) 
        self.mainLayout.addWidget(self.phCheckBox, 0, 0, 1, 1)
        self.phCheckBox.stateChanged.connect(self.enablePh)
        
        self.phSpinBox = QtGui.QDoubleSpinBox()
        self.phSpinBox.setRange(0.0, 14.0)
        self.phSpinBox.setSingleStep(0.5)
        self.phSpinBox.setValue(7.4)
        self.mainLayout.addWidget(self.phSpinBox, 0, 1, 1, 1)
        
        self.forceCheckBox = QtGui.QCheckBox("force")
        self.forceCheckBox.setChecked(False)
        self.mainLayout.addWidget(self.forceCheckBox, 1, 0, 1, 1)

        self.boCheckBox = QtGui.QCheckBox("perceive Bond Order")
        self.boCheckBox.setChecked(True)
        self.mainLayout.addWidget(self.boCheckBox, 1, 1, 1, 1)

        altomListLabel = QtGui.QLabel("Atom list")
        self.mainLayout.addWidget(altomListLabel, 2, 0, 1, 1)
        self.atomListWidget = QtGui.QLineEdit()
        self.mainLayout.addWidget(self.atomListWidget, 2, 1, 1, 1)
        
        self.setLayout(gLayout)
        
    def enablePh(self, checked):
        self.phSpinBox.setEnabled(checked)
        

    def getParams(self):
        ph = None
        if self.phCheckBox.isChecked():
            ph = self.phSpinBox.value()
        return {'pH':ph, 'force': self.forceCheckBox.isChecked(), 'atomList':[],
                'perceiveBondOrder':self.boCheckBox.isChecked()}


class AddHydrogenGUI(QtGui.QDialog):
    def __init__(self, parent=None):
        super(AddHydrogenGUI, self).__init__(parent)
        self.buildUI()

    def buildUI(self):
        radio1 = QtGui.QRadioButton("Reduce")
        radio1.toggled.connect(self.changeMethod)
        #radio1.toggled.connect(CallbackFunction(self.changeMethod, 'reduce'))
        radio2 = QtGui.QRadioButton("OpenBabel")
        #radio2.toggled.connect(CallbackFunction(self.changeMethod, 'openBabel'))
        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(radio1)
        hbox.addWidget(radio2)

        groupBox = QtGui.QGroupBox("Method used for adding hydrogen atoms")
        groupBox.setChecked(False) 
        groupBox.setLayout(hbox)

        self.paramPanelsWidget = QtGui.QStackedWidget()

        w = ReduceParamsGUI('reduce')
        self.reducedPrams = w
        self.paramPanelsWidget.addWidget(w)

        w = OBmolAddHParamsGUI('OBaddH')
        self.obPrams = w
        self.paramPanelsWidget.addWidget(w)

        vertLayout = QtGui.QVBoxLayout()
        vertLayout.addWidget(groupBox)
        vertLayout.addWidget(self.paramPanelsWidget, 1)

        # OK and Cancel buttons
        self.buttons = QtGui.QDialogButtonBox(
            QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel,
            QtCore.Qt.Horizontal, self)
        vertLayout.addWidget(self.buttons)

        self.setLayout(vertLayout)
        radio1.setChecked(True)
        self.radio1 = radio1
        
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)

    def changeMethod(self, checked):
        if checked:
            self.paramPanelsWidget.setCurrentIndex(0)
        else:
            self.paramPanelsWidget.setCurrentIndex(1)

    def getAddHOpt(self):
        if self.radio1.isChecked():
            return 'reduce', self.reducedPrams.getParams()
        else:
            return 'OB', self.obPrams.getParams()
        

    @staticmethod
    def getAddHPrams(molname, parent = None):
        dialog = AddHydrogenGUI(parent)
        dialog.setWindowTitle(QtGui.QApplication.translate("Dialog", "Adding Hydrogen for %s" % molname, None, QtGui.QApplication.UnicodeUTF8))

        result = dialog.exec_()
        opts = dialog.getAddHOpt()
        return (opts, result == QtGui.QDialog.Accepted)
    
if __name__=="__main__":
    import sys
    from optparse import OptionParser
    app = QtGui.QApplication(sys.argv)
    opt, ok = AddHydrogenGUI.getAddHPrams()

    print opt, ok
