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
# Author: Pradeep Anand Ravindranath
#
# Copyright: Pradeep Anand Ravindranath and TSRI 2016
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/AutoSite/GUI/AutoSiteGUI.py,v 1.2 2016/12/07 00:38:33 sanner Exp $
#
# $Id: AutoSiteGUI.py,v 1.2 2016/12/07 00:38:33 sanner Exp $
#
import os, sys, weakref, tempfile, shutil, numpy, string
from time import time
from PySide import QtCore, QtGui
from MolKit2.selection import Selection, SelectionSet
from ADFR.GUI import ICONPATH
from mglutil.util.packageFilePath import getResourceFolderWithVersion
from PmvApp.Pmv import MolApp

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

class BoxParameters(QtGui.QGroupBox):

    def __init__(self, app, parent=None):
        super(BoxParameters, self).__init__(parent)
        self.setTitle('box parameters')
        self.app = weakref.ref(app)
        self.cmdHistory = []
        self.historyPointer = 0

        gridLayout = QtGui.QGridLayout()
        gridLayout.addWidget(QtGui.QLabel("x:"), 0, 1, 1, 1, QtCore.Qt.AlignCenter)
        gridLayout.addWidget(QtGui.QLabel("y:"), 0, 2, 1, 1, QtCore.Qt.AlignCenter)
        gridLayout.addWidget(QtGui.QLabel("z:"), 0, 3, 1, 1, QtCore.Qt.AlignCenter)
        gridLayout.addWidget(QtGui.QLabel("center:"), 1, 0, 1, 1, QtCore.Qt.AlignRight)
        gridLayout.addWidget(QtGui.QLabel("size:"), 2, 0, 1, 1, QtCore.Qt.AlignRight)
        gridLayout.addWidget(QtGui.QLabel("spacing:"), 3, 0, 1, 1, QtCore.Qt.AlignRight)
        gridLayout.addWidget(QtGui.QLabel("smoothing:"), 3, 2, 1, 1, QtCore.Qt.AlignRight)
        gridLayout.addWidget(QtGui.QLabel("cmd:"), 4, 0, 1, 1, QtCore.Qt.AlignRight)

        self._initialized = False
        self.centerSpinBoxes = []
        values = app.boxGeom.center
        for i in range(3):
            widget = QtGui.QDoubleSpinBox()
            widget.valueChanged.connect(self.centerChanged)
            widget.setDecimals(3)
            widget.setMinimum(-9999.)
            widget.setValue(values[i])
            self.centerSpinBoxes.append(widget)
            gridLayout.addWidget(widget, 1, i+1, 1, 1)

        self.sizeSpinBoxes = []
        values = app.boxGeom.sides
        for i in range(3):
            widget = QtGui.QDoubleSpinBox()
            widget.valueChanged.connect(self.sizeChanged)
            widget.setDecimals(3)
            widget.setMinimum(app._spacing)
            widget.setValue(values[i])
            self.sizeSpinBoxes.append(widget)
            gridLayout.addWidget(widget, 2, i+1, 1, 1)

        self.spacingWidget = widget = QtGui.QDoubleSpinBox()
        widget.setDecimals(3)
        widget.setValue(app._spacing)
        widget.setMinimum(0.00001)
        widget.setSingleStep(0.025)
        widget.valueChanged.connect(self.spacingChanged)
        gridLayout.addWidget(widget, 3, 1, 1, 1)
 
        self.smoothWidget = widget = QtGui.QDoubleSpinBox()
        widget.setDecimals(3)
        widget.setValue(app._smooth)
        widget.setMinimum(0.00001)
        widget.valueChanged.connect(self.smoothChanged)
        widget.setSingleStep(0.1)
        gridLayout.addWidget(widget, 3, 3, 1, 1)
        
        self.cmdWidget = MyQLineEdit()
        self.cmdWidget.recall.connect(self.handleRecall)
        gridLayout.addWidget(self.cmdWidget, 4, 1, 1, 3)
        self.cmdWidget.returnPressed.connect(self.executeCmd)

        self.setLayout(gridLayout)
        self._initialized = True

    def sizeChanged(self, value=None):
        if not self._initialized: return
        app = self.app()
        padding = app.paramsWidget.gridPaddingWidget.value()
        lengths = [x.value() for x in self.sizeSpinBoxes]
        app._baseSize[:] = [x-2*padding for x in lengths]
        app.boxGeom.setSides( *lengths)
        app.onBoxChange()

    def centerChanged(self, value=None):
        if not self._initialized: return
        app = self.app()
        center = [x.value() for x in self.centerSpinBoxes]
        app.boxGeom.setCenter( *center)
        app.onBoxChange()

    def spacingChanged(self, value):
        app = self.app()
        app._spacing = value
        app.onBoxChange()

    def smoothChanged(self, value):
        app = self.app()
        app._smooth = value

    ##
    ## methods for cmd entry
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

    def executeCmd(self):
        cmd = self.cmdWidget.text()
        self.cmdHistory.append(cmd)
        self.historyPointer += 1
        self.cmdWidget.clear()
        cmds = cmd.split(';')
        app = self.app()
        for cmd in cmds:
            words = cmd.split()
            if words[0].lower()=='center'[:len(words[0])]:
                if len(words)==4:
                    size = [float(x) for x in words[1:4]]
                    for i in range(3):
                        self.centerSpinBoxes[i].setValue(size[i])
                else:
                    msgBox = QtGui.QMessageBox(self)
                    msgBox.setText("ERROR: Bad syntax. expected center x y z")
                    msgBox.exec_()

            elif words[0].lower()=='size'[:len(words[0])]:
                size = None
                if len(words)==2:
                    vtype, size = self._getIntOtFloat(words[1])
                    size = [size, size, size]                
                elif len(words)==4:
                    s = self.spacingWidget.value()
                    sx = float(words[1])
                    sy = float(words[2])
                    sz = float(words[3])
                    nx = round(sx/s)
                    ny = round(sy/s)
                    nz = round(sz/s)
                    size = [nx*s, ny*s, nz*s]
                else:
                    msgBox = QtGui.QMessageBox(self)
                    msgBox.setText("ERROR: bad syntax. expected size x or size x y z")
                    msgBox.exec_()

                if size is not None:
                    for i in range(3):
                        self.sizeSpinBoxes[i].setValue(size[i])

            elif words[0].lower()=='dimensions'[:len(words[0])]:
                size = None
                if len(words)==2:
                    vtype, size = self._getIntOtFloat(words[1])
                    size = size*self.spacingWidget.value() 
                    size = [size, size, size]
                elif len(words)==4:
                    vtype, size = self._getIntOtFloat(words[1])
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
                else:
                    msgBox = QtGui.QMessageBox(self)
                    msgBox.setText("ERROR: bad syntax. expected nx, ny nz as int for nb grid points")
                    msgBox.exec_()

                if size is not None:
                    for i in range(3):
                        self.sizeSpinBoxes[i].setValue(size[i])

class readLogThread(QtCore.QThread):
    progress = QtCore.Signal(str) # create a custom signal we can subscribe
                                  # to to emit update commands
    progressBar = QtCore.Signal(float)

    def __init__(self, process, logFile, parent=None):
        super(readLogThread,self).__init__(parent)
        self.exiting = False
        self.process = process
        self.logFile = logFile
        self.nbLogLines = 0

    def run(self):
        while True:
            #print 'IN RUN', self.process.poll() is None
            if self.process.poll() == None:
                self.msleep(500)
                #print 'reading STDOUT'
                #for line in iter(self.process.stdout.readline, ""):
                #    self.progress.emit('out: '+line[:-1])
                #print 'reading STDERR'
                #for line in iter(self.process.stderr.readline, ""):
                #    self.progress.emit('err: '+line[:-1])
                #self.msleep(1)
                #print 'reading LOG'
                try:
                    f = open(self.logFile)
                    #self.progressBar.emit(1)
                except IOError:
                    continue
                lines = f.readlines()
                f.close()
                if len(lines)>self.nbLogLines:
                    #print 'NBL', self.nbLogLines
                    text = ''.join(lines[self.nbLogLines:])
                    line = lines[-1]
                    if len(line)>26 and line[26]=='%':
                        #print 'FUGU', float(line.split()[2][:-1])/100.
                        self.progressBar.emit( float(line.split()[2][:-1]) )
                    else:
                        self.progress.emit( 0. )
                    self.progress.emit(text)
                    self.nbLogLines = len(lines)
            else:
                self.progressBar.emit(1.0)
                break
        self.progress.emit('readLogThread THREAD DONE')


class MapsFolderDialog(QtGui.QDialog):

    def __init__(self, name, folder=None, parent=None):
        super(MapsFolderDialog, self).__init__(parent)
        self.setWindowTitle("select folder for saving maps file")

        labelName = QtGui.QLabel("name:")
        self.nameEntry = QtGui.QLineEdit()
        self.nameEntry.setText(name)
        
        label = QtGui.QLabel('Enter folder for saving maps')
        browseButton = QtGui.QPushButton('browse')
        browseButton.clicked.connect(self.browseForFolder)
        
        self.folderPath = QtGui.QLineEdit()
        self.folderPath.returnPressed.connect(self.setPath)
        self.folderPath.setMinimumWidth(400)
        if folder is None:
            self.folderPath.setText(os.getcwd())
        else:
            self.folderPath.setText(folder)

        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok |
                                                QtGui.QDialogButtonBox.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        #self.buttonBox.buttons()[0].setDisabled(True)
        
        layouth1 = QtGui.QHBoxLayout()
        layouth1.addWidget(labelName)
        layouth1.addWidget(self.nameEntry)

        layouth = QtGui.QHBoxLayout()
        layouth.addWidget(browseButton)
        layouth.addWidget(self.folderPath)
        
        layout = QtGui.QVBoxLayout()
        layout.addWidget(label)
        layout.addLayout(layouth1)
        layout.addLayout(layouth)
        layout.addWidget(self.buttonBox)
        self.setLayout(layout)

    def browseForFolder(self):
        dialog = QtGui.QFileDialog(self)
        dialog.setFileMode(QtGui.QFileDialog.Directory)
        dialog.setOption(QtGui.QFileDialog.ShowDirsOnly)
        if dialog.exec_():
            fileNames = dialog.selectedFiles()
            self.folderPath.setText(fileNames[0])
            self.setPath()

    def setPath(self):
        folder = self.folderPath.text()
        if not os.path.exists(folder):
            msgBox = QtGui.QMessageBox(self)
            msgBox.setText("ERROR: Folder %s does not exist"%folder)
            msgBox.exec_()
            self.folderPath.setText('')
            return
        if not os.path.isdir(folder):
            msgBox = QtGui.QMessageBox(self)
            msgBox.setText("ERROR: %s is not a folder"%folder)
            msgBox.exec_()
            self.folderPath.setText('')
            return
        self.buttonBox.buttons()[0].setDisabled(False)
        
class AtomTypeSelector(QtGui.QDialog):

    def __init__(self, atomList, parent=None):
        super(AtomTypeSelector, self).__init__(parent)
        self.setWindowTitle("select atom types")

        self.allButtonOn = QtGui.QCheckBox("select all")
        self.allButtonOn.clicked.connect(self.selectAll)

        self.allButtonOff = QtGui.QCheckBox("deselect all")
        self.allButtonOff.clicked.connect(self.deselectAll)

        self.atypeList = QtGui.QListWidget()
        self.items = []
        for atype, checked in atomList.items():
            item = QtGui.QListWidgetItem(atype)
            self.items.append(item)
            if checked:
                item.setCheckState(QtCore.Qt.Checked)
            else:
                item.setCheckState(QtCore.Qt.Unchecked)
            self.atypeList.addItem(item)
                
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept)

        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.allButtonOn)
        layout.addWidget(self.allButtonOff)
        layout.addWidget(self.atypeList)
        layout.addWidget(self.buttonBox)
        # Set dialog layout
        self.setLayout(layout)

    def selectAll(self):
        for item in self.items:
            item.setCheckState(QtCore.Qt.Checked)
            
    def deselectAll(self):
        for item in self.items:
            item.setCheckState(QtCore.Qt.Unchecked)

    def getTypesString(self):
        s = ''
        for item in self.items:
            if item.checkState()==QtCore.Qt.Checked:
                s += '%s '%item.text()
        return s[:-1]

class ADFRGridMapParametersWidget(QtGui.QWidget):
    
    def __init__(self, parent=None):

        super(ADFRGridMapParametersWidget, self).__init__(parent)
        self._tpok = False
        self._gridok = True
        self.carbon_cutoff = -0.3
        
        ## create group for input
        ##
        self.groupBox1 = QtGui.QGroupBox(self)
        self.groupBox1.setTitle(self.tr("receptor"))
        layout = QtGui.QHBoxLayout()
        self.groupBox1.setLayout(layout)
        w = QtGui.QPushButton("PDBQT ...")
        layout.addWidget(w)
        self.loadRecButton = w

        #w = QtGui.QPushButton("Maps (.zip) ...")
        #layout.addWidget(w)
        #self.loadMapsButton = w
        
        ## create group for ligand
        ##
        self.groupBox2 = QtGui.QGroupBox(self)
        self.groupBox2.setTitle(self.tr("[ligand]"))
        layout = QtGui.QHBoxLayout()
        self.groupBox2.setLayout(layout)
        w = QtGui.QPushButton("PDBQT...")
        layout.addWidget(w)
        self.loadLigButton = w

        ## create group for placing the box
        ##
        self.groupBox3 = QtGui.QGroupBox(self)
        self.groupBox3.setTitle(self.tr("Grid box"))
        layout1 = QtGui.QVBoxLayout(self.groupBox3) # VBox because we might add manual grid widget below later
        self.groupBox3Layout = layout1
        layout = QtGui.QGridLayout() # grid layout to hold bar and padding
        layout1.addLayout(layout)
        
        self.toolBar1 = QtGui.QToolBar()
        self.toolBar1.setOrientation(QtCore.Qt.Horizontal)
        layout.addWidget(self.toolBar1, 0, 0, 2, 1)
        
        self.ASGridAct = act = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'center_grid_AutoSite.png')),
            'AutoSite', self)
        act.setStatusTip('place grid based on AutoSite prediction')
        #act.setDisabled(True)
        #self.toolBar1.addAction(act)
        
        self.recGridAct = act = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'center_grid_REC.png')),
            'grid entire receptor', self)
        act.setStatusTip('make the grid cover the receptor')
        act.setDisabled(True)
        self.toolBar1.addAction(act)

        self.ligGridAct = act = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'center_grid_LIG.png')),
            'grid entire ligand', self)
        act.setStatusTip('make the grid cover the ligand')
        #act.setCheckable(True)
        act.setDisabled(True)
        self.toolBar1.addAction(act)

        self.TPGridAct = act = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'center_grid_TRANS.png')),
            'grid entire translational points', self)
        act.setStatusTip('make the grid cover selected translational points')
        #act.setCheckable(True)
        act.setDisabled(True)
        self.toolBar1.addAction(act)

        self.residuesGrdAct = act = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'center_grid_RES.png')),
            'set grid manually', self)
        act.setStatusTip('set grid around residues')
        act.setCheckable(True)
        act.setDisabled(True)
        self.toolBar1.addAction(act)

        self.manualGrdAct = act = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'center_grid_HAND.png')),
            'set grid manually', self)
        act.setStatusTip('control grid parameters manually')
        act.setCheckable(True)
        act.setDisabled(True)
        self.toolBar1.addAction(act)

        # add grid padding widget
        w = QtGui.QLabel('padding:')
        w.setAlignment(QtCore.Qt.AlignLeft)
        layout.addWidget(w, 0, 1, 1, 1)

        self.gridPaddingWidget = w = QtGui.QDoubleSpinBox()
        w.setDecimals(3)
        w.setValue(4.0)
        layout.addWidget(w, 1, 1, 1, 1)

        #act.triggered.connect(self.gridGUI.setGridFullLigand)
        act.setDisabled(True)
        self.toolBar1.addAction(act)
        self.groupBox3.setDisabled(True)

        self.groupBox4 = QtGui.QGroupBox(self)
        self.groupBox4.setTitle(self.tr("AutoSite"))
        layout = QtGui.QVBoxLayout(self.groupBox4)
        self.ASGridAct = w = QtGui.QPushButton("Identify binding sites")
        layout.addWidget(w)


        ## create group for flexible residues
        ##
        self.groupBox5 = QtGui.QGroupBox(self)
        self.groupBox5.setTitle(self.tr("[Binding Site residues]"))
        layout = QtGui.QGridLayout(self.groupBox5)

        w = self.flexResWidget = MyQLineEdit()
        layout.addWidget(self.flexResWidget, 0, 0, 1, 1)
        #self.flexResWidget.setStyleSheet("QLineEdit{background: grey;}")

        self.setFlexResButton = QtGui.QPushButton(QtGui.QIcon(
            os.path.join(ICONPATH, 'cube.png')),'')
        layout.addWidget(self.setFlexResButton, 0, 1, 1, 1)
        self.setFlexResButton.setCheckable(True)
        self.groupBox4.setDisabled(True)

        self.setFlexResGrdButton = QtGui.QPushButton(QtGui.QIcon(
            os.path.join(ICONPATH, 'center_grid_RES.png')),'')
        layout.addWidget(self.setFlexResGrdButton, 0, 2, 1, 1)
        self.groupBox4.setDisabled(True)

        ## create group for translational points
        ##
        self.groupBox6 = QtGui.QGroupBox(self)
        self.groupBox6.setTitle(self.tr("translation points"))
        layout = QtGui.QVBoxLayout(self.groupBox6)
        self.computeTPointsButton = w = QtGui.QPushButton("compute focused fills")
        layout.addWidget(w)

        self.clustersWidget = w = QtGui.QTableWidget(1,1,self.groupBox6)
        layout.addWidget(w)

        self.featurePtsButton = w = QtGui.QPushButton("compute feature points")
        layout.addWidget(w)
        
        self.groupBox5.setDisabled(True)

        ## create group for computing affinity maps
        ##
        self.groupBox7 = QtGui.QGroupBox(self)
        self.groupBox7.setTitle(self.tr("Save"))
        layout = QtGui.QGridLayout(self.groupBox7)



        hlayout = QtGui.QHBoxLayout()
        self.toolBar2 = QtGui.QToolBar()

        self.TPOKAct = act = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'cancel24.png')),
            'translation points are needed', self)
        act.setStatusTip('Translation points are needed')
        self.toolBar2.addAction(act)

        self.gridOKAct = act = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'ok24.png')),
            'Docking box covers moving atoms', self)
        act.setStatusTip('docking box covers moving receptor atoms')
        act.setDisabled(True)
        self.toolBar2.addAction(act)
        if sys.platform=='darwin':
            self.toolBar2.setMaximumWidth(100)
        elif sys.platform=='linux2':
            self.toolBar2.setMaximumWidth(70)
        elif sys.platform=='win32':
            self.toolBar2.setMaximumWidth(80)
        else:
            self.toolBar2.setMaximumWidth(100)

        hlayout.addWidget(self.toolBar2)

        w = self.computeGridsButton = QtGui.QPushButton("generate maps ... ")
        w.setDisabled(True)
        hlayout.addWidget(w)
        layout.addLayout(hlayout, 2, 0, 1, 3)

        self.progressBar = QtGui.QProgressBar()
        layout.addWidget(self.progressBar, 3, 0, 1, 3)
        self.progressBar.setMinimum(0)
        #self.progressBar.setMaximum(100)
        self.progressBar.setValue(0)

        mainLayout = QtGui.QGridLayout()
        mainLayout.setContentsMargins(2,2,2,2)
        mainLayout.setSpacing(2)
        mainLayout.addWidget(self.groupBox1, 0, 0, 1, 1)
        mainLayout.addWidget(self.groupBox2, 0, 1, 1, 1)
        mainLayout.addWidget(self.groupBox3, 1, 0, 1, 2)
        mainLayout.addWidget(self.groupBox4, 2, 0, 1, 2)
        mainLayout.addWidget(self.groupBox5, 3, 0, 1, 2)
        mainLayout.addWidget(self.groupBox6, 4, 0, 1, 2)
        mainLayout.addWidget(self.groupBox7, 5, 0, 1, 2)
        self.setLayout(mainLayout)

    def configureGenerateMaps(self):
        if self._tpok and self._gridok:
            self.computeGridsButton.setDisabled(False)
        else:
            self.computeGridsButton.setDisabled(True)

    def disableGridCheck(self, value):
        self.gridOKAct.setDisabled(value)
        
    def handleTPointsSignal(self, value):
        if value:
            self.TPGridAct.setDisabled(False)
            #self.TPOKAct.setDisabled(False)
            self.TPOKAct.setStatusTip('translation points are OKAY')
            self.TPOKAct.setText('translation points are OKAY')
            self.TPOKAct.setIcon(QtGui.QIcon(os.path.join(ICONPATH, 'ok24.png')))
            self._tpok = True
        else:
            self.TPGridAct.setDisabled(True)
            #self.TPOKAct.setDisabled(True)
            self.TPOKAct.setStatusTip('translation points needed')
            self.TPOKAct.setText('translation points needed')
            self.TPOKAct.setIcon(QtGui.QIcon(os.path.join(ICONPATH, 'cancel24.png')))
            self._tpok = False
        self.configureGenerateMaps()
        
    def handleGridOKSignal(self, value):
        if value:
            self.gridOKAct.setIcon(QtGui.QIcon(os.path.join(ICONPATH, 'ok24.png')))
            self.gridOKAct.setText('docking box covers moving receptor atoms')
            self.gridOKAct.setStatusTip('docking box covers moving receptor atoms')
            self._gridok = True
        else:
            self._gridok = False
            self.gridOKAct.setIcon(QtGui.QIcon(os.path.join(ICONPATH, 'cancel24.png')))
            self.gridOKAct.setText('docking box does NOT covers moving receptor atoms')
            self.gridOKAct.setStatusTip('docking box does NOT covers moving receptor atoms')
        self.configureGenerateMaps()

class GridGUI(QtGui.QWidget):

    jobDone = QtCore.Signal()
    flexResChanged = QtCore.Signal()
    TPointsOK = QtCore.Signal(bool)
    gridOK = QtCore.Signal(bool)
    
    def __init__(self, parent=None, pmv=None, pmvViewer=None, eventHandler=None):
        super(GridGUI, self).__init__(parent)

        self._labDisplayMode = 0 # 0: hidden, 1: residues in box, 2: flexRes
        self._LL = [0, 0, 0] # lower left box corner
        self._UR = [1, 1, 1] # upper right box corner
        self._baseSize = [0,0,0]
        self._smooth = 0.5
        self._spacing = 0.375
        self._clusterItems = []
        self._noFlexResParse = False # set to True to set entry without parsing string
        self._flexRes = [] # list of prody residues for flexible residues
        self._flexResSCAtoms = None # MolKit2 selection with all moving side chain atoms in flexible residues
        self._flexResAtoms = []

        self.jobDone.connect(self.finishUp)
        self.flexResChanged.connect(self.onFlexResChanged)
        
        self.getAllADatomTypes()
        self.receptor = None # will be he prody molecule for the receptor
        self.ligand = None # will be he prody molecule for the ligand

        if pmv is None:
            self.createPmvApp()
        else:
            assert isinstance(Pmv, MolApp)
            self.pmv = pmv

        self.buildUI(pmvViewer=pmvViewer, parent=parent)

        self.tmpFolder = os.path.join(getResourceFolderWithVersion(), 'tmp')
        if not os.path.exists(self.tmpFolder):
            os.mkdir(self.tmpFolder)
        elif not os.path.isdir(self.tmpFolder):
            msgBox = QtGui.QMessageBox(self)
            msgBox.setText("ERROR: %s is not a folder, please remove this file"%self.tmpFolder)
            msgBox.setStandardButtons(QtGui.QMessageBox.Ok)
            msgBox.exec_()
            return

        def _flexResChanged(frstr):
            self.paramsWidget.flexResWidget.setText(frstr)
            self.setFlexRes()

    def getAllADatomTypes(self):
        from ADFRcc.adfr import Parameters
        parameters = Parameters.getParameters()
        self.ADAtomTypes = {}
        for i in range(parameters.numAtomTypes):
            self.ADAtomTypes[parameters.getAtomTypeByIndex(i).atomTypeName] = True

    def createPmvApp(self, eventHandler=None):
        pmv = MolApp(eventHandler=eventHandler)
        pmv.trapExceptions = False
        pmv.lazyLoad('bondsCmds', package='PmvApp')
        pmv.lazyLoad('fileCmds', package='PmvApp')
        pmv.lazyLoad('displayCmds', package='PmvApp')
        pmv.lazyLoad('editCmds', package='PmvApp')
        pmv.displayLines.loadCommand()
        pmv.lazyLoad("colorCmds", package="PmvApp")
        pmv.color.loadCommand()
        pmv.lazyLoad("selectionCmds", package="PmvApp")
        pmv.lazyLoad('deleteCmds', package='PmvApp')
        pmv.lazyLoad('labelCmds', package='PmvApp')

        pmv.lazyLoad('msmsCmds', package='PmvApp')
        pmv.lazyLoad('displayHyperBallsCmds', package='PmvApp')
        pmv.lazyLoad('interactionsCmds', package='PmvApp')
        pmv.lazyLoad('coarseMolecularSurfaceCmds', package='PmvApp')

        pmv.setOnAddObjectCmd('Molecule', [pmv.displayLines, pmv.colorByAtomType])
        self.pmv = pmv

    def setGridVisible(self, value):
        # value is 0 for unchecked and 2 for checked for checkbox
        # not(value==0) make it work for 0, 1, 2, False, True
        self.boxGeom.master.Set(visible = not(value==0))
        for c in self.boxGeom.master.children:
            if c.name=='faces':
                c.Set(visible = 0)
            else:
                c.Set(visible = not(value==0))

    def buildUI(self, pmvViewer=None, parent=None):

        w = self.paramsWidget = ADFRGridMapParametersWidget(parent)
        w.loadRecButton.clicked.connect(self.getReceptorFilename)
        #w.loadMapsButton.clicked.connect(self.getReceptorMapsFilename)
        w.loadLigButton.clicked.connect(self.getLigandFilename)
        #w.editAtypesButton.clicked.connect(self.pickATypes)
        w.computeGridsButton.clicked.connect(self.computeGrids)
        w.ASGridAct.clicked.connect(self.setGridAutoSite)
        w.computeTPointsButton.clicked.connect(self.computeTPoints)
        w.clustersWidget.itemClicked.connect(self.showHideCluster)
        w.featurePtsButton.clicked.connect(self.computeFPoints)
        w.flexResWidget.returnPressed.connect(self.setFlexRes)
        w.setFlexResButton.clicked.connect(self.showFlexResChooser)
        w.setFlexResGrdButton.clicked.connect(self.setGridFlexRes)
        w.gridPaddingWidget.valueChanged.connect(self.paddingChanged)
        #w.compAllButton.clicked.connect(self.toggleComputeAll)
        

        w.recGridAct.triggered.connect(self.setGridFullReceptor)
        w.ligGridAct.triggered.connect(self.setGridFullLigand)
        w.TPGridAct.triggered.connect(self.setGridClusterPoints)
        w.residuesGrdAct.triggered.connect(self.showResiduesGridControls)
        w.manualGrdAct.triggered.connect(self.showManualGridControls)
        
        if sys.platform=='darwin':
            self.paramsWidget.setMaximumWidth(400)
        elif sys.platform=='linux2':
            self.paramsWidget.setMaximumWidth(300)
        elif sys.platform=='win32':
            self.paramsWidget.setMaximumWidth(400)
        else:
            self.paramsWidget.setMaximumWidth(400)
        self.TPointsOK.connect(w.handleTPointsSignal)
        self.gridOK.connect(w.handleGridOKSignal)
        
        ## create the vertical tool bar
        ##
        self.toolBarV = QtGui.QToolBar()
        self.toolBarV.setOrientation(QtCore.Qt.Vertical)

        act = self.focusBoxAction = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'crosshair.png')),
            'Center On Grid', self)
        act.setStatusTip('focus on docking box')
        act.triggered.connect(self.focusOnGrid)
        self.toolBarV.addAction(act)

        self.toolBarV.addSeparator()

        act = self.showReceptorAction = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'receptor2_NO.png')),
            'Show Receptor', self)
        act.setStatusTip('show/hide receptor')
        act.triggered.connect(self.showHideReceptor)
        act.setCheckable(True)
        act.setDisabled(True)
        self.toolBarV.addAction(act)
        
        act = self.showLigandAction = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'ligand64_NO.png')),
            'Show Receptor', self)
        act.setStatusTip('show/hide ligand')
        act.triggered.connect(self.showHideLigand)
        act.setCheckable(True)
        act.setDisabled(True)
        self.toolBarV.addAction(act)
        
        act = self.showGridBoxAction = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'cube.png')),
            'Show Receptor', self)
        act.setStatusTip('show/hide labels docking box')
        act.triggered.connect(self.showHideGridBox)
        act.setCheckable(True)
        act.setChecked(True)
        self.toolBarV.addAction(act)
        
        act = self.showResLabAction = QtGui.QAction(
            QtGui.QIcon(os.path.join(ICONPATH, 'labels.png')),
            'Show Receptor', self)
        act.setStatusTip('show/hide labels ')
        act.triggered.connect(self.showHideResidueLabels)
        act.setCheckable(True)
        act.setDisabled(True)
        self.toolBarV.addAction(act)
        
        
        from PmvApp.GUI.Qt.PmvGUI import PmvViewer
        if pmvViewer is None:
            ## create the Pmv Viewer
            self.viewer = PmvViewer(self.pmv, master=self)
        else:
            assert isisntance(pmvViewer, PmvViewer)
            self.Viewer = viewer

        from DejaVu2.Box import NiceBox
        self.boxGeom = NiceBox('gridOutline')
        self.boxGeom.addToViewer(self.viewer)
        self.setGridVisible(True)

        from DejaVu2.Points import Points
        self.TPoints = Points('tpoints', visible=0, inheritMaterial=False,
                              materials=[(1,0,0)], inheritPointWidth=False,
                              pointWidth=4.)
        self.viewer.AddObject(self.TPoints)

        # place widgets
        mainLayout = QtGui.QHBoxLayout()
        mainLayout.addWidget(self.paramsWidget)
        mainLayout.addWidget(self.toolBarV)
        mainLayout.addWidget(self.viewer.cameras[0])

        self.setLayout(mainLayout)

    def getReceptorMapsFilename(self):
        filename, selfilter = QtGui.QFileDialog().getOpenFileName(
            self, self.tr("Read Receptor Maps"), '',
            self.tr("zip Files (*.zip);; All files (*)"))
        
        if filename:
            self.getReceptorMaps(filename.encode('ascii', 'replace'))

    def getReceptorMaps(self, filename):
        from zipfile import ZipFile
        from ADFRcc.adfr import GridMap
        zf = ZipFile(filename)
        zf.extractall(self.tmpFolder)
        folder = os.path.join(self.tmpFolder, os.path.splitext(os.path.basename(filename))[0])
        # find receptor file name i.e PDBQT file which is NOT receptor.pdbqt
        for filename in zf.namelist():
            if filename.endswith('.pdbqt'):
                filename = os.path.split(filename)[1]
                recname = os.path.splitext(filename)[0]
                if recname != 'receptor':
                    break

        self.loadReceptor(os.path.join(folder, '%s.pdbqt'%recname))

        # Find out all map types and load 'e' map for map info 
        atomTypes = ''
        nb = 0
        flexResStr = None
        _mapLoaded = False
        for mapFileName in zf.namelist():
            w = mapFileName.split('.')
            if w[-1]=='map':
                atomTypes += '%s '%w[-2]
                nb += 1
                if not _mapLoaded and not (w[-2]=='e' or w[-2]=='d'): # read 'e' map to get center, box size, spacing, and flexres
                    _map = GridMap()
                    _f, name = os.path.split(mapFileName)
                    _map.loadFromMapFile(w[-2], folder, name)
                    center = _map.getCenterPy()
                    nx, ny, nz = _map.getNumGridPointsPy()
                    self._spacing = spacing = _map.getDistBetweenGridPoints()
                    sx = spacing*(nx-1)
                    sy = spacing*(ny-1)
                    sz = spacing*(nz-1)
                    print 'ASDASda', mapFileName, center, sx, sy, sz
                    flexResStr = _map.getFlexRes()
                    _mapLoaded = True

        #print 'SADASDA', nb, atomTypes, 'F', flexResStr, 'F', center, sx, sy, sz, spacing
        # set box
        self.boxGeom.setCenter( *center)
        self.boxGeom.setSides( sx, sy, sz)
        self._baseSize[:] = (sx, sy, sz)
        self.paramsWidget.gridPaddingWidget.setValue(0.0)

        # set flexRes
        if flexResStr:
            self.paramsWidget.flexResWidget.setText(flexResStr)
            self.setFlexRes()
            
        # set atom types
        if nb==len(self.ADAtomTypes): # all types
            self.paramsWidget.atypesLabel.setText('')
            self.paramsWidget.compAllButton.setChecked(True)
        else:
            self.paramsWidget.compLigandAtypes.setChecked(True)
            self.paramsWidget.atypesLabel.setText(atomTypes[:-1])

        # set TPoints
            tp = numpy.load(os.path.join(folder, 'translationPoints.npy'))
            self._clusters = [[[], tp, [], []]]
            self._allClusterColors = [[0, 1, 0]]

            clustersWidget = self.paramsWidget.clustersWidget
            clustersWidget.clear()
            self._clusterItems = []
            newItem = QtGui.QListWidgetItem()
            self._clusterItems = [newItem]
            newItem.setText('All Clusters')
            newItem.setCheckState(QtCore.Qt.Checked)
            clustersWidget.insertItem(0, newItem)

            self.TPoints.Set(visible=True, vertices=tp, materials=[[0, 1, 0]])
            self.TPointsOK.emit(True)

    def getReceptorFilename(self):
        filename, selfilter = QtGui.QFileDialog.getOpenFileName(
            self, self.tr("Read Receptor"), '',
            self.tr("PDBQT Files (*.pdbqt);; All files (*)"))
        if filename:
            self.loadReceptor(filename.encode('ascii', 'replace'))
        
    def loadReceptor(self, filename):
        if self.receptor:
            self.pmv.deleteMolecule(self.receptor)
        mols = self.pmv.readMolecules(filename)
        self.setReceptor(mols[0][0])
        #import pdb; pdb.set_trace()
        self.pmv.customColor(mols[0][0].select('element C'), [(0.,1.,1.)], geomsToColor=['lines'])
        self.viewer.Reset_cb()
        self.viewer.Normalize_cb()
        self.viewer.Center_cb()

    def getLigandFilename(self):
        filename, selfilter = QtGui.QFileDialog.getOpenFileName(
            self, self.tr("Read Ligand"), '',
            self.tr("PDBQT Files (*.pdbqt);; All files (*)"))
        if filename:
            self.loadLigand(filename.encode('ascii', 'replace'))

    def loadLigand(self, filename):
        if self.ligand:
            self.pmv.deleteMolecule(self.ligand)
        mols = self.pmv.readMolecules(filename)
        self.setLigand(mols[0][0])
        self.pmv.displaySticksAndBalls(mols[0][0])
        self.pmv.customColor(mols[0][0].select('element C'), [(1.,1.,0.)])
        self.viewer.Reset_cb()
        self.viewer.Normalize_cb()
        self.viewer.Center_cb()
       
    def setReceptor(self, mol):
        self.receptor = mol
        self.showReceptorAction.setDisabled(False)
        self.showReceptorAction.setChecked(True)
        self.showResLabAction.setDisabled(False)
        #self.defaultReceptorBox()
        self.setGridFullReceptor()
        self.paramsWidget.groupBox3.setDisabled(False)
        self.paramsWidget.ASGridAct.setDisabled(False)
        self.paramsWidget.recGridAct.setDisabled(False)
        self.paramsWidget.manualGrdAct.setDisabled(False)
        self.paramsWidget.residuesGrdAct.setDisabled(False)
        self.paramsWidget.groupBox4.setDisabled(False)
        self.paramsWidget.groupBox5.setDisabled(False)
        self.paramsWidget.groupBox6.setDisabled(False)
        self.pmv.computeMSMS(self.receptor, pRadius=1.5, density=6.0)
        
    def setLigand(self, mol):
        self.ligand = mol
        atypes = numpy.unique(mol._ag.getData('AD_element'))
        if 'C' not in atypes:
            atypes.insert(0, 'C')
        self.paramsWidget.groupBox3.setDisabled(False)
        self.showLigandAction.setDisabled(False)
        self.showLigandAction.setChecked(True)
        self.paramsWidget.ligGridAct.setDisabled(False)
        self.paramsWidget.editAtypesButton.setDisabled(False)
        self.paramsWidget.atypesLabel.setText(' '.join(atypes))
        self.paramsWidget.compLigandAtypes.setChecked(True)

    def updateDisplay(self):
        # handle lines of receptor
        if self.showReceptorAction.isChecked():
            resWithAtomsInBox = self.getResWithAtomsInBox()
            gc = self.receptor.geomContainer
            opac = numpy.ones( len(gc.allCoords) )*0.3
            if resWithAtomsInBox:
                opac[resWithAtomsInBox.getIndices()] = 1.
            gc.geoms['singleBonds'].Set(
                opacity=opac, transparent=True, inheritMaterial=False, polyFace='front')
            self.pmv.displayLines(self.receptor)
        else:
            if self.receptor:
                #self.pmv.displayLines(self.receptor, negate=True)
                self.pmv.displayMSMS(self.receptor, negate=True)

        # handle residue labels
        if self._labDisplayMode==0 and self.receptor: # hide all residue labels
            self.pmv.labelResidues(self.receptor, negate=True)
        elif self._labDisplayMode==1: # label residues in box
            atoms =  self.selectCAInBox()
            if atoms:
                self.pmv.labelResidues(SelectionSet([atoms]))
        elif self._labDisplayMode==2: # label flexres
            #import pdb; pdb.set_trace()
            if self.receptor:
                self.pmv.labelResidues(self.receptor, negate=True)
            if self._flexResAtoms:
                self.pmv.labelResidues(SelectionSet([self._flexResAtoms]))
            
    def checkBoxCoversFlexRes(self):
        # check the ligand atoms are inside the box
        if not self._flexResAtoms:
            self.paramsWidget.gridOKAct.setIcon(QtGui.QIcon(os.path.join(ICONPATH, 'ok24.png')))
            #self.paramsWidget.disableGridCheck(True)
        else:
            #self.paramsWidget.disableGridCheck(False)
            boxselstr = 'x>%f and x<%f and y>%f and y<%f and z>%f and z<%f'%(
                self._LL[0], self._UR[0], self._LL[1],
                self._UR[1], self._LL[2], self._UR[2])
            if len(self._flexResAtoms) == len(self._flexResAtoms.select(boxselstr)):
                self.paramsWidget.gridOKAct.setIcon(QtGui.QIcon(os.path.join(ICONPATH, 'ok24.png')))
            else:
                self.paramsWidget.gridOKAct.setIcon(QtGui.QIcon(os.path.join(ICONPATH, 'cancel24.png')))
            
    def onBoxChange(self, invalidateClusterPoints=True):
        cx, cy, cz = self.boxGeom.center
        sx, sy, sz = self.boxGeom.sides
        self._LL = [cx - (sx*.5), cy - (sy*.5), cz - (sz*.5)]
        self._UR = [cx + (sx*.5), cy + (sy*.5), cz + (sz*.5)]

        if invalidateClusterPoints:
            self.paramsWidget.clustersWidget.clear()
            self.TPoints.Set(visible=0)
            self.paramsWidget.computeGridsButton.setDisabled(True)
            self.paramsWidget.TPGridAct.setDisabled(True)
            self._clusterItems = []

        self.checkBoxCoversFlexRes()
        if invalidateClusterPoints:
            self.TPointsOK.emit(False)
        self.updateDisplay()
        
    def showHideReceptor(self):
        self.updateDisplay()

    def showHideLigand(self):
        self.pmv.showMolecules(self.ligand,
                               show=self.showLigandAction.isChecked())
    def showHideGridBox(self):
        self.setGridVisible(self.showGridBoxAction.isChecked())
        self.TPoints.Set(visible=self.showGridBoxAction.isChecked())
        
    def showHideResidueLabels(self):
        self._labDisplayMode = (self._labDisplayMode+1)%3
        if self._labDisplayMode==2 and self._flexResSCAtoms is None:
            self._labDisplayMode=0
        self.updateDisplay()
        
    def focusOnGrid(self):
        oldRecVis = oldLigVis = None
        if self.receptor:
            oldRecVis = self.receptor.geomContainer.geoms['master'].visible
            self.receptor.geomContainer.geoms['master'].Set(visible=False)
        if self.ligand:
            oldLigVis = self.ligand.geomContainer.geoms['master'].visible
            self.ligand.geomContainer.geoms['master'].Set(visible=False)
        oldGridVis = self.showGridBoxAction.isChecked()
        if not oldGridVis:
            self.setGridVisible(True)
        rot = self.viewer.rootObject.rotation
        self.viewer.Reset_cb()
        self.viewer.Normalize_cb()
        self.viewer.Center_cb()
        self.viewer.rootObject.Set(rotation=rot)
        self.setGridVisible(oldGridVis)
        if oldRecVis is not None:
            self.receptor.geomContainer.geoms['master'].Set(visible=oldRecVis)
        if oldLigVis is not None:
            self.ligand.geomContainer.geoms['master'].Set(visible=oldLigVis)

    def selectCAInBox(self):
        # get list of CA atoms for residues with moving side chain atoms in box
        if self.receptor is not None:
            serials = ' '.join([str(x) for x in self.selectAtomsInBox().getSerials()])
            selstr = 'same residue as (serial %s) and ca'%serials
            inBox = self.receptor.select(selstr)
            if inBox is not None:
                return inBox

    def selectAtomsInBox(self):
        # select moving side chains atoms in the docking box
        if self.receptor is not None:
            selstr = '(not bb) and (not name CB) and (not resname PRO ALA GLY) and '
            selstr += 'x>%f and x<%f and y>%f and y<%f and z>%f and z<%f'%(
                self._LL[0], self._UR[0], self._LL[1],
                self._UR[1], self._LL[2], self._UR[2])
            inBox = self.receptor.select(selstr)
            if inBox is not None:
                return self.removeNH(inBox)

    def removeNH(self, atoms):
        # remove H atoms attached to N atoms
        toRemove = []
        for atom in atoms:
            if atom.getElement()=='H':
                for natom in atom.iterBonded():
                    if natom.isbackbone:
                        toRemove.append(atom.getIndex())
        atoms = atoms - Selection(self.receptor._ag, toRemove, '')
        return atoms

    def getResWithAtomsInBox(self):
        # select all moving receptor side chain atoms inside the box
        # then expand to full residues
        inBoxAtoms = self.selectAtomsInBox()
        if inBoxAtoms:
            serials = ' '.join([str(x) for x in inBoxAtoms.getSerials()])
            selStr = 'same residue as (serial %s)'%serials
            resWithAtomsInBox = self.receptor.select(selStr)
            return resWithAtomsInBox

    def pickATypes(self):
        self.paramsWidget.compLigandAtypes.setChecked(True)
        self.paramsWidget.atypesLabel.setDisabled(False)
        d = {}
        atstr = self.paramsWidget.atypesLabel.text().split()
        for at in self.ADAtomTypes.keys():
            if at in atstr:
                d[at] = True
            else:
                d[at] = False
        dialog = AtomTypeSelector(d, self)
        result = dialog.exec_()
        self.paramsWidget.atypesLabel.setText( dialog.getTypesString() )

    def getFlexResDialogInput(self):
        # convert selection string into dict with key 'Chid:ResnameResnum' and
        # value True is this residue in 
        resNames = []
        df = {}
        for res in self._flexRes:
            df['%c:%s%d'%(res.getChid(), res.getResname(), res.getResnum())] = True

        d = {}
        for atom in self.selectCAInBox():
            key = '%s:%s%d'%(atom.getChid(), atom.getResname(), atom.getResnum())
            if df.has_key(key):
                d[key] = True
            else:
                d[key] = False
        return d

    def flexResSelStrFromCA(self, caAtoms):
        # build flexRes selection string from a list of CA atoms
        d = {}
        for atom in caAtoms:
            chid = atom.getChid()
            resnum = atom.getResnum()
            resname = atom.getResname()
            if not d.has_key(chid):
                d[chid] = ['%s%d'%(atom.getResname(),atom.getResnum())]
            else:
                d[chid].append('%s%d'%(atom.getResname(), atom.getResnum()))
        chids = d.keys()
        chids.sort()
        s = ''
        descr = []
        for chid in chids:
            s += '%c:'%chid
            resl = []
            for res in d[chid]:
                s += '%s,'%res
                resl.append((res[:3], res[3:]))
            descr. append( (chid, resl) )
            s = s[:-1]
            s += ';'
        s = s[:-1]
        return s, descr
    
    def showFlexResChooser(self):
        def setFR():
            caAtoms = self._flexResDialog.resTreeWidget.getAtoms()
            if caAtoms:
                allResAtoms = self.receptor.select('same residue as index %s'%' '.join(
                    [str(x) for x in caAtoms.getIndices()]))
                self._flexResAtoms = allResAtoms
                self._flexResSCAtoms = self.removeNH(allResAtoms.select('not bb and not ca'))
                selStr, flexResDescr = self.flexResSelStrFromCA(caAtoms)
            else:
                self._flexResAtoms = []
                self._flexResSCAtoms = []
                selStr = ''
            self._noFlexResParse = True
            self.paramsWidget.flexResWidget.setText(selStr)
            self._noFlexResParse = False
            self.flexResChanged.emit()
            
        def accept():
            self.paramsWidget.flexResWidget.setDisabled(False)
            self.paramsWidget.setFlexResGrdButton.setDisabled(False)
            self.paramsWidget.setFlexResButton.setChecked(False)

        if self.paramsWidget.setFlexResButton.isChecked():
            self.paramsWidget.flexResWidget.setDisabled(True)
            self.paramsWidget.setFlexResGrdButton.setDisabled(True)
            from ADFR.GUI.ResiduesTree import ResiduesTreeDialog
            caInbox = self.selectCAInBox()
            checked = {}
            for ca in caInbox:
                checked['%s:%d'%(ca.getChid(),ca.getResnum())] = False
            if self._flexResSCAtoms:
                selstr, descr = self.flexResSelStrFromCA(self._flexResSCAtoms)
                for chid, residues in descr:
                    for resname, resnum in residues:
                        checked['%s:%s'%(chid,resnum)] = True
            self._flexResDialog = ResiduesTreeDialog(caInbox, self.ligand, checked)
            self._flexResDialog.resTreeWidget.resChanged.connect(setFR)
            self._flexResDialog.accepted.connect(accept)
            self._flexResDialog.show()
        else:
            self.paramsWidget.flexResWidget.setDisabled(False)
            self.paramsWidget.setFlexResGrdButton.setDisabled(False)
            self._flexResDialog.accept()
            
    ## def defaultReceptorBox(self):
    ##     coords = self.receptor._ag.getCoords()
    ##     mini = numpy.min(coords, 0)
    ##     maxi = numpy.max(coords, 0)
    ##     lengths = (maxi-mini) + 2*self.paramsWidget.gridPaddingWidget.value()
    ##     center = mini + 0.5*lengths
    ##     self._baseSize[:] = lengths - 2*self.paramsWidget.gridPaddingWidget.value()
    ##     self.boxGeom.setCenter( *center)
    ##     self.boxGeom.setSides( *lengths)
    ##     self.onBoxChange()

    def showResiduesGridControls(self):
        def setBox():
            atoms = self._resControlWidget.resTreeWidget.getAtoms()
            if atoms:
                coords = atoms.getCoords()
                mini = numpy.min(coords, 0)
                maxi = numpy.max(coords, 0)
                lengths = (maxi-mini) + 2*self.paramsWidget.gridPaddingWidget.value()
                center = mini + 0.5*lengths
                self._baseSize[:] = lengths - 2*self.paramsWidget.gridPaddingWidget.value()
                self.boxGeom.setCenter( *center)
                self.boxGeom.setSides( *lengths)

                #nc1 = numpy.min( (mini, self._LL), 0) # new corner1
                #nc2 = numpy.max( (maxi, self._UR), 0) # new corne2
                #nc = 0.5*(nc1 + nc2) # new center
                #nl = (nc2 - nc1)
                #self._baseSize[:] = nl - 2*self.paramsWidget.gridPaddingWidget.value()
                #self.boxGeom.setCenter( *nc)
                #self.boxGeom.setSides( *nl)
                self.onBoxChange()

        def accept():
            self.paramsWidget.toolBar1.setDisabled(False)
            self.paramsWidget.residuesGrdAct.setChecked(False)

        if self.paramsWidget.manualGrdAct.isChecked():
            self.paramsWidget.manualGrdAct.setChecked(False)
            self.showManualGridControls()

        if self.paramsWidget.residuesGrdAct.isChecked():
            from ADFR.GUI.ResiduesTree import ResiduesTreeDialog
            self._resControlWidget = ResiduesTreeDialog(self.receptor.select(), 
                                                        self.ligand)
            self._resControlWidget.resTreeWidget.resChanged.connect(setBox)
            self._resControlWidget.accepted.connect(accept)
            self._resControlWidget.show()
            self.paramsWidget.toolBar1.setDisabled(True)
        else:
            self.paramsWidget.toolBar1.setDisabled(False)
            self._resControlWidget.accept()
            
    def showManualGridControls(self):
        if self.paramsWidget.manualGrdAct.isChecked():
            self._manualControlWidget = BoxParameters(self)
            self.paramsWidget.groupBox3Layout.addWidget(self._manualControlWidget)
        else:
            self._manualControlWidget.deleteLater()

    def setGridAutoSite(self):
        temp_name = next(tempfile._get_candidate_names())
        folder = os.path.join(self.tmpFolder, temp_name)
        os.mkdir(folder)
        spacing = 1.0
        coords = self.receptor._ag.getCoords()
        mini = numpy.min(coords, 0)
        maxi = numpy.max(coords, 0)
        center = 0.5*(mini+maxi)
        sizef = (maxi-mini) + 2*self.paramsWidget.gridPaddingWidget.value()
        size = [int(round(x/spacing)) for x in sizef]
        self.AutoSiteFill(center = center,size = size, spacing=1.0, folder=folder)

        
    def setGridFlexRes(self):
        coords = self._flexResAtoms.getCoords()
        mini = numpy.min(coords, 0)
        maxi = numpy.max(coords, 0)
        #center = 0.5*(mini+maxi)
        #lengths = (maxi-mini) + 2*self.paramsWidget.gridPaddingWidget.value()

        padding = self.paramsWidget.gridPaddingWidget.value()
        nc1 = numpy.min( (mini-padding, numpy.array(self._LL)), 0) # new corner1
        nc2 = numpy.max( (maxi+padding, numpy.array(self._UR)), 0) # new corne2
        nc = 0.5*(nc1 + nc2) # new center
        nl = (nc2 - nc1)
        self._baseSize[:] = nl - 2*self.paramsWidget.gridPaddingWidget.value()
        self.boxGeom.setCenter( *nc)
        self.boxGeom.setSides( *nl)
        self.onBoxChange()

    def setGridFullReceptor(self):
        if self.paramsWidget.manualGrdAct.isChecked():
            self.paramsWidget.manualGrdAct.setChecked(False)
            self.showManualGridControls()
        coords = self.receptor._ag.getCoords()
        mini = numpy.min(coords, 0)
        maxi = numpy.max(coords, 0)
        center = 0.5*(mini+maxi)
        lengths = (maxi-mini) + 2*self.paramsWidget.gridPaddingWidget.value()
        self._baseSize[:] = lengths - 2*self.paramsWidget.gridPaddingWidget.value()
        self.boxGeom.setCenter( *center)
        self.boxGeom.setSides( *lengths)
        self.onBoxChange()
        
    def setGridFullLigand(self):
        if self.paramsWidget.manualGrdAct.isChecked():
            self.paramsWidget.manualGrdAct.setChecked(False)
            self.showManualGridControls()
        coords = self.ligand._ag.getCoords()
        mini = numpy.min(coords, 0)
        maxi = numpy.max(coords, 0)
        center = 0.5*(mini+maxi)
        lengths = (maxi-mini) + 2*self.paramsWidget.gridPaddingWidget.value()
        self._baseSize[:] = lengths - 2*self.paramsWidget.gridPaddingWidget.value()
        self.boxGeom.setCenter( *center)
        self.boxGeom.setSides( *lengths)
        self.onBoxChange()

    def setGridClusterPoints(self):
        if self.paramsWidget.manualGrdAct.isChecked():
            self.paramsWidget.manualGrdAct.setChecked(False)
            self.showManualGridControls()
        if self.paramsWidget.residuesGrdAct.isChecked():
            self.paramsWidget.residuesGrdAct.setChecked(False)
            self.showResiduesGridControls()
        coords = self.getClusterPoints()
        mini = numpy.min(coords, 0)
        maxi = numpy.max(coords, 0)
        center = 0.5*(mini+maxi)
        lengths = (maxi-mini) + 2*self.paramsWidget.gridPaddingWidget.value()
        self._baseSize[:] = lengths - 2*self.paramsWidget.gridPaddingWidget.value()
        self.boxGeom.setCenter( *center)
        self.boxGeom.setSides( *lengths)
        self.onBoxChange(invalidateClusterPoints=False)

    def paddingChanged(self, value):
        x,y,z = self._baseSize
        self.boxGeom.setSides( x+2*value, y+2*value, z+2*value )
        self.onBoxChange()

    def setFlexRes(self):
        ## expect string of type chid1:Resnum1,Resnum2,Resnum3 ;  chid2:Resnum1,Resnum2,Resnum3
        ## if ":" is missing chid will match all chains
        ## residues can also be specified using ResnameResnum and can beseparated by space or ,
        if self._noFlexResParse:
            return
        self._flexRes = []
        #self._flexResStr = [] # list of (
        self._flexResAtoms = None # MolKit2 selection with all atoms in flexible residues
        self._flexResSCAtoms = None # MolKit2 selection with all moving side chain atoms in flexible resid

        frstr = self.paramsWidget.flexResWidget.text().encode('ascii', 'replace')
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
                    try:
                        resnum = int(res)
                    except ValueError:
                        errors.append('ERROR: invalid syntax\n expect a residue number but got "%s"\n'%res)
                        continue
                else: # name and number
                    resname = res[:3]
                    try:
                        resnum = int(res[3:])
                    except ValueError:
                        errors.append('ERROR: invalid syntax\n expect a residue number but got "%s\n"'%res[3:])
                        continue
                    
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
        if self._flexResAtoms:
            self._flexResSCAtoms = self.removeNH(self._flexResAtoms.select('not bb and not ca'))
        #print "FRA", len(self._flexResAtoms), len(self._flexResSCAtoms), self._flexResSCAtoms.getResnames()
        # check that grid covers flexRec atoms
        #flexRecInside, outsideAtoms = self.atomInBox(self._flexResSCAtoms)
        #if not flexRecInside:
        #    outSideResnames = ''
        #    resnums = outsideAtoms.getResnums()
        #    chids = outsideAtoms.getChids()
        #    done = {}
        #    outSideRes = []
        #    for chid, resnum in zip(chids, resnums):
        #        key = '%s:%d'%(chid, resnum)
        #        if done.has_key(key): continue
        #        done[key] = True
        #        res = hv.getResidue(chid, resnum)
        #        outSideRes.append(res)
        #        outSideResnames += '%s%d '%(res.getResname(), res.getResnum())
        #    msgBox = QtGui.QMessageBox(self)
        #    msgBox.setText('atoms for residues %s are outside the grid\nWould you like change the grid to cover these atoms ? If you choose NO this(ese) residue(s) will bre removed from the flexible residues list'%outSideResnames)
        #    msgBox.setStandardButtons(QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        #    ret =  msgBox.exec_()

        #    if ret == QtGui.QMessageBox.Yes: # change grid to cover flexres atoms
        #        coords = self._flexResSCAtoms.getCoords()
        #        mini = numpy.min(coords, 0)
        #        maxi = numpy.max(coords, 0)
        #        center = [x.value() for x in self.centerSpinBoxes]
        #        size = [x.value() for x in self.sizefSpinBoxes]
        #        corner1 = center[0]-size[0]*.5, center[1]-size[1]*.5, center[2]-size[2]*.5
        #        corner2 = center[0]+size[0]*.5, center[1]+size[1]*.5, center[2]+size[2]*.5
        #        nc1 = numpy.min( (mini, corner1), 0) # new corner1
        #        nc2 = numpy.max( (maxi, corner2), 0) # new corne2
        #        nc = 0.5*(nc1 + nc2) # new center
        #        for i in range(3):
        #            self.centerSpinBoxes[i].setValue(nc[i])
        #            self.sizefSpinBoxes[i].setValue(nc2[i]-nc1[i]+2*self.gridPaddingWidget.value())
        #            self._baseSize = [x.value() for x in self.sizefSpinBoxes]
        #    else: # remove offending residues from flexRes list
        #        for res in outSideRes:
        #            self._flexRes.remove(res)
        #        selStr, descr = self.buildFlexRec(self._flexRes)
        #        #import pdb; pdb.set_trace()
        #        self.paramsWidget.flexResWidget.setText(selStr)
        #        self.flexResDialog.fill(self.getFlexResDialogInput())
        #        self.setFlexRes()
        self.flexResChanged.emit()
        
        if len(errors):
            msgBox = QtGui.QMessageBox(self)
            msgBox.setText('\n'.join(errors))
            msgBox.setStandardButtons(QtGui.QMessageBox.Ok)
            ret =  msgBox.exec_()
            
        #print self._flexRes
        #print self._flexResStr
        #print selStr
        # string might have been expanded to add residues names
        self.paramsWidget.flexResWidget.setText(selStr)

    def atomInBox(self, atoms):
        # check is atoms are in the grid box, If not identify atoms outside
        center = self.boxGeom.center
        size = self.boxGeom.sides
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

    def toggleComputeAll(self):
        if self.paramsWidget.compAllButton.isChecked():
            #self.atypesEntry.setDisabled(True)
            self.paramsWidget.atypesLabel.setDisabled(True)
        else:
            #self.atypesEntry.setDisabled(False)
            self.paramsWidget.atypesLabel.setDisabled(True)

#######################################################################END
    
            
    ## def setAtypes(self):
    ##     #atstring = self.atypesEntry.text()
    ##     atstring = self.paramsWidget.atypesLabel.text()
    ##     atstring = atstring.replace(',', ' ')
    ##     atypes = []
    ##     for atype in atstring.split():
    ##         if not self.ADAtomTypes.has_key(atype):
    ##             msgBox = QtGui.QMessageBox(self)
    ##             msgBox.setText("ERROR: %s is not a valid AutoDock5 atom type"%atype)
    ##             msgBox.setStandardButtons(QtGui.QMessageBox.Ok)
    ##             ret =  msgBox.exec_()
    ##         else:
    ##             atypes.append(atype)
    ##     #self.atypesEntry.setText(' '.join(atypes))
    ##     self.paramsWidget.atypesLabel.setText(' '.join(atypes))
    ##     self.compLigandAtypes.setChecked(True)

    def onFlexResChanged(self):

            #self.paramsWidget.disableGridCheck(True)
        self.pmv.displaySticksAndBalls(self.receptor, negate=True)
        if self._flexResSCAtoms:
            self.pmv.displaySticksAndBalls(self._flexResAtoms)
            self.pmv.customColor(self._flexResAtoms.select('element C'), [(1.,0.5,0.)],
                                 geomsToColor=['sticks', 'balls'])
            gc = self.receptor.geomContainer
            # FIXME for some reason the stickas and balls get the transparency of the lines
            gc.geoms['sticks'].Set(transparent=False)
            gc.geoms['balls'].Set(transparent=False)
            #gc.geoms['sticks'].Set(opacity=0.8)
            #gc.geoms['balls'].Set(opacity=0.8)
        self.updateDisplay()

    def AutoSiteFill(self,center = None,size = None, spacing=None, folder=None, background=True ):
        from AutoSite.fillBuriedness import Buriedness
        coords = self.receptor._ag.getCoords()
        radiiR = self.receptor._ag.getRadii()
        if center is None and size is None:
            gc = self._computeGrids(atypes=['C','OA','HD'], spacing=spacing, useFlexRec=False,
                                folder=folder, background=background, atypesOnly=False, fp=True)
        else:
            gc = self._computeGrids(center, size, atypes=['C','OA','HD'], spacing=spacing, useFlexRec=False,
                                folder=folder, atypesOnly=False, fp=True)
        ###gc.getTPoints()
        gc.getASPoints()
        #print 'Computed %d TPoints in '%len(gc._coords)
        #import pdb; pdb.set_trace()
        from AutoSite.utils.clusterTPoints import DensityClustering
        dcl = DensityClustering([spacing,spacing,spacing])
        dcl.findClustersD(gc._indices)

        self.colorMap = [(0,1,0), (1,0.5,0), (1,1,0), (0,1,1), (1,0,1), ]
        colors = [None]*len(gc._indices)
        nbLargeCl = 0
        self._clusters = [[gc._indices, gc._coords,gc._potential,gc._atype]]
        
        
        smallClustersCoords = []
        smallClustersInds = []
        smallClustersPots = []
        smallClustersAtypes = []
        lclusters = []
        clProp = []
        for length, clinds in zip(dcl._clen, dcl._clusters):
            if length<50:
                color = (1,0,0)
                smallClustersInds.extend( gc._indices[clinds] )
                smallClustersCoords.extend( gc._coords[clinds] )
                smallClustersPots.extend( gc._potential[clinds] )
                smallClustersAtypes.extend( gc._atype[clinds] )
                
            else:
                color = self.colorMap[nbLargeCl%len(self.colorMap)]
                nbLargeCl += 1
                dmetric = Buriedness(gc._coords[clinds], gc._atype[clinds], coords.tolist(),radiiR )
                clc = gc._coords[clinds]
                cx,cy,cz = numpy.mean(clc,axis=0)
                dst2 = 0.0
                for px,py,pz in clc:
                    dst2 += (px-cx)*(px-cx)+(py-cy)*(py-cy)+(pz-cz)*(pz-cz)
                complen = len(clc)
                Rg2 = (dst2/complen)
                from math import sqrt
                Rg = sqrt(Rg2)
                fburied = dmetric.NumericalBurriedness()
                metric = (len(clc)*fburied*fburied)/Rg
                lclusters.append( [gc._indices[clinds], gc._coords[clinds], gc._potential[clinds], gc._atype[clinds], metric] )
                clProp.append([metric, len(clc),fburied, Rg])
                
            for ind in clinds:
                colors[ind] = color
        if len(lclusters) !=0 :
            tmpclsort = sorted(lclusters,key=lambda x:x[4],reverse = True)
            self._clusters.extend([x[0],x[1],x[2],x[3]] for x in tmpclsort)
        self._clusters.append( [smallClustersInds, smallClustersCoords, smallClustersPots, smallClustersAtypes])
        self._allClusterColors = colors
        clustersWidget = self.paramsWidget.clustersWidget
        clustersWidget.clear()
        self._clusterItems = []
        data = ['All', 'Volume', 'Rg','Buriedness','metric']
                #import pdb;pdb.set_trace()
        for j in range (len(data)):
            if j !=0:
                clustersWidget.insertColumn(j)
            newItem = QtGui.QTableWidgetItem(str(data[j]))
            clustersWidget.setItem(0,j, newItem)
            if j == 0:
                self._clusterItems.append(newItem)
                newItem.setCheckState(QtCore.Qt.Unchecked)

        #newItem = QtGui.QTableWidgetItem()
        #newItem.setText('All')
        #newItem.setCheckState(QtCore.Qt.Checked)
        #self._clusterItems.append(newItem)
        #clustersWidget.insertRow(1)
        #clustersWidget.setItem(0,0, newItem)
        i = 1
        for inds, coords,pots,atypes in self._clusters[1:]:
            clustersWidget.insertRow(i)
            if i < len(self._clusters)-1:
                data = [i, len(inds), clProp[i-1][3], clProp[i-1][2], clProp[i-1][0]]
            else:
                data = ['Small cluster (<50)', len(inds), '-','-','-']
                #import pdb;pdb.set_trace()
            for j in range (len(data)):
                #if j !=0 and i ==1:
                #    clustersWidget.insertColumn(j)
                newItem = QtGui.QTableWidgetItem(str(data[j]))
                clustersWidget.setItem(i,j, newItem)
                if j == 0:
                    self._clusterItems.append(newItem)
                    newItem.setCheckState(QtCore.Qt.Unchecked)
            i += 1
        #clustersWidget.show()
        self.TPoints.Set(visible=True, vertices=gc._coords, materials=colors)
        self.TPointsOK.emit(True)

        shutil.rmtree(folder)


    def computeTPoints(self):
        # create scratch folder name
        temp_name = next(tempfile._get_candidate_names())
        folder = os.path.join(self.tmpFolder, temp_name)
        os.mkdir(folder)
        self.AutoSiteFill(center = None,size = None, spacing=0.75, folder=folder, background=False )

    def computeFPoints(self):
        import numpy
        self.TPoints.Set(visible=0)
        #self.TPointsOK.emit(False)
        clusterPoints = self.getClusterPoints()
        clusterPotentials = self.getClusterPotentials()
        clusterAtomtypes = self.getClusterAtypes()
        if len(clusterPoints)==0:
            msgBox = QtGui.QMessageBox(self)
            msgBox.setText("ERROR: no translation points selected")
            msgBox.exec_()
            return
        ccoords = [x for i,x in enumerate(clusterPoints) if clusterAtomtypes[i] == 'C']
        cpot = [x for i,x in enumerate(clusterPotentials) if clusterAtomtypes[i] == 'C']
        ocoords = [x for i,x in enumerate(clusterPoints) if clusterAtomtypes[i] == 'O']
        opot = [x for i,x in enumerate(clusterPotentials) if clusterAtomtypes[i] == 'O']
        hcoords = [x for i,x in enumerate(clusterPoints) if clusterAtomtypes[i] == 'H']
        hpot = [x for i,x in enumerate(clusterPotentials) if clusterAtomtypes[i] == 'H']
        from AutoSite.ASfeaturepoints import featurePts
        fp = featurePts(self.receptor,ccoords, cpot, ocoords, opot, hcoords, hpot)
        from DejaVu2.Spheres import Spheres
        self.TSpheres1 = Spheres('tsph1', visible=0, inheritMaterial=False,\
                              materials=[(0,1,0)], inheritPointWidth=False,\
                              pointWidth=4.)
        self.TSpheres1.Set(visible=True, vertices=numpy.array(fp.cfp),materials=[(0,1,0)])
        self.viewer.AddObject(self.TSpheres1)
        self.TSpheres2 = Spheres('tsph2', visible=0, inheritMaterial=False,\
                              materials=[(0,0,1)], inheritPointWidth=False,\
                              pointWidth=4.)
        self.TSpheres2.Set(visible=True, vertices=numpy.array(fp.ofp),materials=[(1,0,0)])
        self.viewer.AddObject(self.TSpheres2)
        self.TSpheres3 = Spheres('tsph3', visible=0, inheritMaterial=False,\
                              materials=[(1,1,1)], inheritPointWidth=False,\
                              pointWidth=4.)
        self.TSpheres3.Set(visible=True, vertices=numpy.array(fp.hfp),materials=[(1,1,1)])
        self.viewer.AddObject(self.TSpheres3)
        
        ## flexResStr = ''
        ## for x in fp.bsres:
        ##     flexResStr = flexResStr+x+';'
        #print flexResStr
        ## import pdb;pdb.set_trace()
        fresdict = {}
        for x in fp.bsres:
            xsplit = x.split(':')
            fresdict.setdefault(xsplit[0], [])
            fresdict[xsplit[0]].append(xsplit[1])

        flexResStr = ''
        for k,v in fresdict.iteritems():
            resstr =  "".join(str(x)+',' for x in v)
            flexResStr = flexResStr+k+':'+resstr+';'
            

        print flexResStr
        self.paramsWidget.flexResWidget.setText(flexResStr)
        self.setFlexRes()
        #self.flexResChanged.emit() 

        
        
    def showHideCluster(self, item):
        clustersWidget = self.paramsWidget.clustersWidget
        row = clustersWidget.row(item)
        #import pdb; pdb.set_trace()
        if row==0: # all
            if item.checkState() == QtCore.Qt.CheckState.Unchecked: # we clicked to uncheck
                self.TPoints.Set(visible=False)
                self.TPointsOK.emit(False)
            else: # we clicked to check ALL clusters
                for n in range(1, clustersWidget.rowCount()):
                    clustersWidget.item(n,0).setCheckState(QtCore.Qt.Unchecked)
                self.TPoints.Set(visible=True, vertices=self._clusters[0][1],
                                 materials=self._allClusterColors)
                self.TPointsOK.emit(True)
        else:
            clustersWidget.item(0,0).setCheckState(QtCore.Qt.Unchecked)
            coords = []
            colors = []
            for n in range(1, clustersWidget.rowCount()):
                if clustersWidget.item(n,0).checkState()==QtCore.Qt.Checked:
                    c = self._clusters[n][1]
                    coords.extend( c )
                    colors.extend( [self.colorMap[(n-1)%len(self.colorMap)]]*len(c) )
            if coords:
                self.TPoints.Set(visible=True, vertices=coords, materials=colors)
                self.TPointsOK.emit(True)

            else:
                self.TPoints.Set(visible=False)
                self.TPointsOK.emit(False)
                #self.TSpheres1.Set(visible=False)
                #self.TSpheres2.Set(visible=False)
                #self.TSpheres3.Set(visible=False)
                self.paramsWidget.flexResWidget.setText('')
                self.setFlexRes()

    def isBoxOnReceptor(self):
        if not self.gridGUI.isBoxOnReceptor():
            msgBox = QtGui.QMessageBox(self)
            msgBox.setText("WARNING: Grid box does not overlap with receptor")
            msgBox.exec_()
            return False
        else:
            return True

    def getClusterPoints(self):
        item0 = self._clusterItems[0]
        if item0.checkState() == QtCore.Qt.Checked: # all
            clusterPoints = self._clusters[0][1] # cluster points coordinates
        else:
            clusterPoints = [] # cluster points coordinates
            i = 1
            #import pdb; pdb.set_trace()
            for item in self._clusterItems[1:]:
                if item.checkState() == QtCore.Qt.Checked:
                    #clusterPoints.extend( self._clusters[i][1].tolist() )
                    clusterPoints.extend( self._clusters[i][1] )
                i += 1
        return clusterPoints

    def getClusterPotentials(self):
        item0 = self._clusterItems[0]
        if item0.checkState() == QtCore.Qt.Checked: # all
            clusterPotentials = self._clusters[0][2] # cluster points potentials
        else:
            clusterPotentials = [] # cluster points coordinates
            i = 1
            #import pdb; pdb.set_trace()
            for item in self._clusterItems[1:]:
                if item.checkState() == QtCore.Qt.Checked:
                    #clusterPoints.extend( self._clusters[i][1].tolist() )
                    clusterPotentials.extend( self._clusters[i][2] )
                i += 1
        return clusterPotentials

    def getClusterAtypes(self):
        item0 = self._clusterItems[0]
        if item0.checkState() == QtCore.Qt.Checked: # all
            clusterAtomtypes = self._clusters[0][3] # cluster points atom types
        else:
            clusterAtomtypes = [] # cluster points coordinates
            i = 1
            #import pdb; pdb.set_trace()
            for item in self._clusterItems[1:]:
                if item.checkState() == QtCore.Qt.Checked:
                    #clusterPoints.extend( self._clusters[i][1].tolist() )
                    clusterAtomtypes.extend( self._clusters[i][3] )
                i += 1
        return clusterAtomtypes
    
    def computeGrids(self):
        recName = os.path.splitext(os.path.split(self.receptor.filename)[1])[0]

        # ask for destination folder
        dialog = MapsFolderDialog(recName, parent=self)
        destinationFolder = None
        if dialog.exec_():
            destinationFolderPath = dialog.folderPath.text()
            destinationFolder = dialog.nameEntry.text()
        else:
            return
        if destinationFolderPath is None:
            return

        self.destinationFolder = destinationFolder
        self.destinationFolderPath = destinationFolderPath
        
        #import pdb; pdb.set_trace()
        # create scratch folder name
        folder = os.path.join(self.tmpFolder, destinationFolder)

        if os.path.exists(folder):
            shutil.rmtree(folder)
        os.mkdir(folder)

        # save translation points
        clusterPoints = self.getClusterPoints()
        if len(clusterPoints)==0:
            msgBox = QtGui.QMessageBox(self)
            msgBox.setText("ERROR: no translation points selected")
            msgBox.exec_()
            return

        if self.paramsWidget.compAllButton.isChecked():
            atypes = self.ADAtomTypes.keys()
        else:
            atypes = [x.encode('ascii', 'replace') for x in self.paramsWidget.atypesLabel.text().split()]

        # write translation points
        filename = os.path.join(folder, 'translationPoints.npy')
        numpy.save(filename, clusterPoints)

        background = True
        self.paramsWidget.computeGridsButton.setDisabled(True)
        gc = self._computeGrids(atypes=atypes, useFlexRec=True, folder=folder,
                                background=background)

        if not background:
            self.makeTarFile()

    def makeTarFile(self):
        # create tar file
        try:
            cwd = os.getcwd()
            os.chdir(self.tmpFolder)
            shutil.make_archive(self.destinationFolder, 'zip', '.',
                                self.destinationFolder)
            shutil.copy(self.destinationFolder+'.zip', self.destinationFolderPath)
        finally:
            os.chdir(cwd)
            self.paramsWidget.computeGridsButton.setDisabled(False)
        # copy tar to user-specified location

    def _computeGrids(self, center=None, size=None, atypes=None, spacing=None,
                      useFlexRec=False, folder='.',
                      background=False, atypesOnly=False, fp = False):
        #from ADFR.utils.MakeGrids import CalculateAD4Grids
        if center is None:
            center = self.boxGeom.center[:]
        if spacing is None:
            spacing = self._spacing
            if size is None:
                size = [int(round(x/spacing)) for x in self.boxGeom.sides]
        else:
            if size is None:
                size = [int(round(x/spacing)) for x in self.boxGeom.sides]

        if atypes is None:
            atypes = [x.encode('ascii', 'replace') for x in self.paramsWidget.atypesLabel.text().split()]
        #import pdb;pdb.set_trace()
        flexResDescr=[]
        selStr = ''
        if useFlexRec:
            #import pdb; pdb.set_trace()
            if self._flexResAtoms:
                selStr, flexResDescr = self.flexResSelStrFromCA(
                    self._flexResAtoms.select('ca'))
        ## gc = CalculateAD4Grids(self.receptor, center, size, atypes,
        ##                        spacing=spacing,
        ##                        smooth=self._smooth,
        ##                        flexibleResidues=flexResDescr,
        ##                        folder=folder, atypesOnly=atypesOnly)
        from AutoSite.compositePoints import CompositePoints
        gc = CompositePoints(self.receptor, center, size, atypes,
                               spacing=spacing,
                               smooth=self._smooth,
                               flexibleResidues=flexResDescr,
                               folder=folder, atypesOnly=atypesOnly, fp = fp)
        t0 = time()
        status = gc.run(background=background)

        if not background:
            if status==0:
                print 'Computed grids for %s in '%atypes, time()-t0
            else:
                msgBox = QtGui.QMessageBox(self)
                msgBox.setText("ERROR: running autogrid failed")
                msgBox.exec_()
            # changed headers
            if useFlexRec and len(selStr):
                gc.addFlexRecHeader('FLEXRES %s'%selStr)
        else:
            logFile = os.path.join(folder, gc.temp_name+'.glg')
            logReader = readLogThread(status, logFile)        
            logReader.progress.connect(self.updateText, QtCore.Qt.QueuedConnection)
            logReader.progressBar.connect(self.updateProgress, QtCore.Qt.QueuedConnection)
            if not logReader.isRunning():
                logReader.start()
            self.logReader = logReader
            if useFlexRec and len(selStr):
                self._selStr = selStr
            else:
                self._selStr = None
            self._gc = gc
        return gc

    def finishUp(self):
        if self._selStr:
            self._gc.addFlexRecHeader('FLEXRES %s'%self._selStr)
        self.makeTarFile()
        
    def updateText(self,text):
        #print '++++++++++++++++++++++++++++++++++++++++++'
        #print len(text)
        #self.text.append(text)
        if text=='readLogThread THREAD DONE':
            self.paramsWidget.progressBar.setValue(0)
            self.jobDone.emit()
            
    def updateProgress(self, value):
        self.paramsWidget.progressBar.setValue(value)
        
if __name__=='__main__':
    import sys

    from optparse import OptionParser
    parser = OptionParser(usage="usage: %prog ligand mapFolder [options] filename",
                      version="%prog 0.1")
    parser.add_option("-r", "--receptor",
                      action="store", # optional because action defaults to "store"
                      dest="receptorFile",
                      help="receptor PDBQT file",)
    parser.add_option("-l", "--ligand",
                      action="store", # optional because action defaults to "store"
                      dest="ligandFile",
                      help="ligand PDBQT file",)
    parser.add_option("-m", "--maps",
                      action="store", # optional because action defaults to "store"
                      dest="receptorMaps",
                      help="zip file containing receptor maps and receptor structure",)

    (options, args) = parser.parse_args()
    app = QtGui.QApplication(sys.argv)
    widget = GridGUI()

    if len(sys.argv)==2:
        widget.pmv.readMolecules(sys.argv[1])
    widget.resize(1000,400)
    widget.show()
    if options.receptorFile:
        widget.loadReceptor(options.receptorFile)
    if options.ligandFile:
        widget.loadLigand(options.ligandFile)
    if options.receptorMaps:
        widget.getReceptorMaps(options.receptorMaps)
    
    #widget.viewer.OneRedraw()
    widget.raise_()
    timer = QtCore.QTimer()
    timer.singleShot(200, widget.viewer.OneRedraw)
    timer.start(1)

    sys.exit(app.exec_())
   
