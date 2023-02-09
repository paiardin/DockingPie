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

# $Header:
#
# $Id:

from PySide import QtGui, QtCore

class SetDistanceDialog(QtGui.QDialog):
    def __init__(self, parent=None, callBack=None):
        super(SetDistanceDialog, self).__init__(parent)
        self.callBack = callBack
        self.buildUI()

    def buildUI(self):
        self.groupBox = QtGui.QGroupBox()
        groupLayout = QtGui.QHBoxLayout()
        distLabel = QtGui.QLabel("Set distance:")    
        self.distBox = QtGui.QDoubleSpinBox()
        self.distBox.setSingleStep(0.5)
        self.distBox.setMinimum(0)
        if self.callBack:
            self.distBox.valueChanged.connect(self.callBack)
        self.distBox.setValue(5.0)
        groupLayout.addWidget(distLabel)
        groupLayout.addWidget(self.distBox)
        self.groupBox.setLayout(groupLayout)
        vertLayout = QtGui.QVBoxLayout()
        vertLayout.addWidget(self.groupBox)
        
        # OK and Cancel buttons
        self.buttons = QtGui.QDialogButtonBox(
            QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel,
            QtCore.Qt.Horizontal, self)
        vertLayout.addWidget(self.buttons)
        self.setLayout(vertLayout)

        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)

    def getParams(self):
        return {"distance": self.distBox.value()}

    @staticmethod
    def getDialogParams(parent = None, title=None, callBack=None):
        dialog = SetDistanceDialog(parent, callBack=callBack)
        if not title:
            title = "Set distance"
        dialog.setWindowTitle(QtGui.QApplication.translate("Dialog", title, None, QtGui.QApplication.UnicodeUTF8))
        dialog.setWindowModality(QtCore.Qt.WindowModal) # WindowModal
        result = dialog.exec_()
        opts = dialog.getParams()
        return (opts, result == QtGui.QDialog.Accepted)



if __name__=="__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    opt, ok = SetDistanceDialog.getDialogParams()
    print opt, ok
