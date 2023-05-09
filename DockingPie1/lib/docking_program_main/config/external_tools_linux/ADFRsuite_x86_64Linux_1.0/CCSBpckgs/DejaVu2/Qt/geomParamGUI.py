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

#
# $Header: /mnt/raid/services/cvs/DejaVu2/Qt/geomParamGUI.py,v 1.1.4.1 2017/07/13 22:25:38 annao Exp $
#
# $Id: geomParamGUI.py,v 1.1.4.1 2017/07/13 22:25:38 annao Exp $
from PySide import QtCore, QtGui
from mglutil.util.callback import CallbackFunction
import os

class BondOrderWidgets(QtGui.QWidget):
    toggleBO = QtCore.Signal(bool)
    doubleBondSeparation = QtCore.Signal(float)
    tripleBondSeparation = QtCore.Signal(float)
    aromaticValueChanged =  QtCore.Signal(float)
    
    def __init__(self, parent=None, labels={"db": "double bond separation:", "tb":"triple bond separation:", "aromatic":"aromatic line width:"}, values={"db":(0.12, [0.1, 10.0]), "tb":(0.2, [0.1, 10.0]), "aromatic":(1, [1, 10])}):
        super(BondOrderWidgets, self).__init__()
        self.parent = parent
        self.immediate = True
        self.layout = QtGui.QVBoxLayout()
        #display BO widgets
        aromaticWidth, aromRange = values["aromatic"]
        dbSep, dbRange = values["db"]
        tbSep, tbRange = values["tb"]
        self.displayBO = QtGui.QGroupBox("bond order")
        self.displayBO.setCheckable(True)
        self.displayBO.setChecked(True)
        displayBOLayout = QtGui.QGridLayout()
        # Double Bond 
        lab1 = QtGui.QLabel(labels["db"])
        self.dbSepSpinBox = dbSepSpinBox = QtGui.QDoubleSpinBox()
        if dbSep == "auto":
            dbSepSpinBox.setSpecialValueText("auto")
            dbRange[0]=0.0
            dbSep = 0.0
        dbSepSpinBox.setRange(*dbRange)#(0.1, 10.0)
        dbSepSpinBox.setSingleStep(0.01)
        dbSepSpinBox.setValue(dbSep)
        dbSepSpinBox.valueChanged[float].connect(self.setDBSeparation_cb)
        # Triple Bond
        lab2 = QtGui.QLabel(labels["tb"])
        self.tbSepSpinBox = tbSepSpinBox = QtGui.QDoubleSpinBox()
        if tbSep == "auto":
            tbSepSpinBox.setSpecialValueText("auto")
            tbRange[0]=0.
            tbSep = 0.
        tbSepSpinBox.setRange(*tbRange)#(0.1, 10.0)
        tbSepSpinBox.setSingleStep(0.01)
        tbSepSpinBox.setValue(tbSep)
        #tbSepSpinBox.setEnabled(displayBO)
        tbSepSpinBox.valueChanged[float].connect(self.setTBSeparation_cb)
        # Aromatic
        lab3 =  QtGui.QLabel(labels["aromatic"])
        if type(aromaticWidth) == int:
            self.aromLineWidthSpinBox =  aromLineWidthSpinBox = QtGui.QSpinBox()
            aromLineWidthSpinBox.setSingleStep(1)
        else:
            self.aromLineWidthSpinBox =  aromLineWidthSpinBox = QtGui.QDoubleSpinBox()
            aromLineWidthSpinBox.setSingleStep(0.01)
        aromLineWidthSpinBox.setRange(*aromRange)#(1., 10.)
        aromLineWidthSpinBox.setValue(aromaticWidth)
        aromLineWidthSpinBox.valueChanged.connect(self.setAromatic_cb)
        #add widgets to the layout
        displayBOLayout.addWidget(lab1, 0, 0)
        displayBOLayout.addWidget(dbSepSpinBox, 0, 1)
        displayBOLayout.addWidget(lab2, 1, 0)
        displayBOLayout.addWidget(tbSepSpinBox, 1, 1)
        displayBOLayout.addWidget(lab3, 2, 0)
        displayBOLayout.addWidget(aromLineWidthSpinBox, 2, 1)     
        self.displayBO.toggled.connect(self.toggleDisplayBO_cb)
        self.displayBO.setLayout(displayBOLayout)
        self.displayBO.adjustSize() # this is to fix "qDrawShadeRect: Invalid parameters" warning (???) (A.O.)
        self.layout.setContentsMargins(0,0,0,0)
        self.layout.addWidget(self.displayBO)
        self.setLayout(self.layout)

    def toggleDisplayBO_cb(self, val):
        #print "toggleDisplayBO_cb", val
        if self.immediate:
            self.toggleBO.emit(val)

    def setDBSeparation_cb(self, val):
        #print "setDBSeparation_cb", val, type(val)
        if self.immediate:
            self.doubleBondSeparation.emit(val)
        
    def setTBSeparation_cb(self, val):
        if self.immediate:
            self.tripleBondSeparation.emit(val)

    def setAromatic_cb(self, val):
        if self.immediate:
            #print "setAromatic_cb", val, type(val)
            self.aromaticValueChanged.emit(val)

    def getWidgetsValues(self):
        kw = {}
        kw['displayBondOrder'] = self.displayBO.isChecked()
        kw['doubleBondSep'] = self.dbSepSpinBox.value()
        kw['tripleBondSep'] = self.tbSepSpinBox.value()
        kw['aromaticLinewidth'] = self.aromLineWidthSpinBox.value()
        return kw

    def setWidgetsValues(self, kw):
        self.displayBO.setChecked(kw.get('displayBondOrder', True))
        self.tbSepSpinBox.setValue(kw.get('tripleBondSep',0.2))
        self.dbSepSpinBox.setValue(kw.get('doubleBondSep',0.12))
        self.aromLineWidthSpinBox.setValue(kw.get('aromaticLinewidth',1))

class LinesModelGlobalParams(QtGui.QWidget):
    """User input form for setting the global parameters of
    lines geometry such as:
    width, stipple length and space, bond order parameters(
    double & triple bond separation and aromatic linewidth), soft color boundaries. """
    
    # SIGNALS
    lineWidthChanged=QtCore.Signal(int)
    stippleLengthChanged=QtCore.Signal(float)
    stippleSpaceChanged=QtCore.Signal(float)
    softColorBoundariesToggled=QtCore.Signal(bool)
    # + BondOrderWidgets signals

    def __init__(self, geoms=None, parent=None):
        super(LinesModelGlobalParams, self).__init__()
        self.geoms = geoms
        self.parent= parent
        self.addWidgets()
        self.immediate = True

    def addWidgets(self):
        self.propLayout = propLayout = QtGui.QVBoxLayout()
        # get initial values of line properties
        lineWidth = 2
        stippleSpace = 0.2
        stippleLength = 0.2
        stippleLines = False
        #line width widgets
        lineWidthLayout = QtGui.QHBoxLayout()
        # linewidth spinbox 
        linewidthLabel = QtGui.QLabel("line width:")
        self.linewidthSpinBox = linewidthSpinBox = QtGui.QSpinBox()
        linewidthSpinBox.setRange(1, 10)
        linewidthSpinBox.setSingleStep(1)
        linewidthSpinBox.setValue(lineWidth)
        linewidthSpinBox.valueChanged[int].connect(self.setLineWidth_cb)

        lineWidthLayout.addWidget(linewidthLabel)
        lineWidthLayout.addWidget(self.linewidthSpinBox)
        # stipple lines length and space group with 2 spinboxes
        self.stippleGroup = QtGui.QGroupBox("stipple")
        stippleLinesLayout = QtGui.QGridLayout()
        lab1 = QtGui.QLabel("length:")
        self.stippleLengthSpinBox = stippleLengthSpinBox = QtGui.QDoubleSpinBox()
        stippleLengthSpinBox.setRange(0.1, 10.0)
        stippleLengthSpinBox.setSingleStep(0.01)
        stippleLengthSpinBox.setValue(stippleLength)
        lab2 = QtGui.QLabel("space:")
        self.stippleSpaceSpinBox = stippleSpaceSpinBox = QtGui.QDoubleSpinBox()
        stippleSpaceSpinBox.setRange(0.1, 10.0)
        stippleSpaceSpinBox.setSingleStep(0.1)
        stippleSpaceSpinBox.setValue(stippleSpace)
        stippleLengthSpinBox.valueChanged[float].connect(self.setStippleLength_cb)
        stippleSpaceSpinBox.valueChanged[float].connect(self.setStippleSpace_cb)
        stippleLinesLayout.addWidget(lab1, 0, 0)
        stippleLinesLayout.addWidget(stippleLengthSpinBox, 0, 1)
        stippleLinesLayout.addWidget(lab2, 1, 0)
        stippleLinesLayout.addWidget(stippleSpaceSpinBox, 1, 1)
        self.stippleGroup.setLayout(stippleLinesLayout)
        self.stippleGroup.adjustSize() # this is to prevent "qDrawShadeRect: Invalid parameters" warning (???)
        # bonorder widgets group 
        self.bondOrderGroup = BondOrderWidgets()
        
        propLayout.addLayout(lineWidthLayout)
        self.softBoundaries = QtGui.QCheckBox("soft color boundaries")
        self.softBoundaries.toggled.connect(self.colorBoundaries_cb)
        propLayout.addWidget(self.stippleGroup)
        propLayout.addWidget(self.bondOrderGroup)#.displayBO)
        propLayout.addWidget(self.softBoundaries)
        #propLayout.addWidget(self.stippleGroup)
        self.setLayout(propLayout)

    #CALLBACKS
    def setLineWidth_cb(self, val):
        if self.immediate:
            self.lineWidthChanged.emit(val)

    def setStippleLength_cb(self, val):
        if self.immediate:
            self.stippleLengthChanged.emit(val)

    def setStippleSpace_cb(self, val):
        if self.immediate:
            self.stippleSpaceChanged.emit(val)

    def saveWidgetsValues(self):
        # Form widgets values
        initval = self.initialFormValues={}
        initval['linewidth'] = self.linewidthSpinBox.value()
        initval.update(self.bondOrderGroup.getWidgetsValues())
        initval['softColorBoundaries'] = self.softBoundaries.isChecked()

    def restoreWidgetsValues(self):
        # restore initial form values
        self.immediate = False
        kw = self.initialFormValues
        bo = {}
        bo['displayBondOrder'] = kw['displayBondOrder']
        bo['doubleBondSep'] = kw['doubleBondSep']
        bo['tripleBondSep'] = kw['tripleBondSep']
        bo['aromaticLinewidth'] = kw['aromaticLinewidth'] 
        self.linewidthSpinBox.setValue(kw.get('linewidth', 2))
        self.stippleLengthSpinBox.setValue(kw.get('stippleLength', 0.2))
        self.stippleSpaceSpinBox.setValue(kw.get('stippleSpace', 0.2))
        self.softBoundaries.setChecked(kw['softColorBoundaries'])
        self.immediate = True
        self.bondOrderGroup.immediate = False
        self.bondOrderGroup.setWidgetsValues(bo)
        self.bondOrderGroup.immediate = True

    def colorBoundaries_cb(self, val):
        if self.immediate:
            self.softColorBoundariesToggled.emit(val)


from mglutil.util.packageFilePath import findFilePath        
PMVICONPATH = findFilePath('Icons', 'PmvApp.GUI')

class StrokeWidgets(QtGui.QWidget):
    #SIGNALS
    toggleStroke = QtCore.Signal(bool, int)
    backFaceColorModeChanged = QtCore.Signal(str)
    backFaceColorChanged = QtCore.Signal(float, float, float, float, bool, bool)
    lineWidthChanged = QtCore.Signal(int)
    
    def __init__(self, gui, geoms=None, parent=None):
        super(StrokeWidgets, self).__init__()
        self.parent = parent
        self.geoms = geoms
        self.parent= parent
        self.addWidgets()
        self.atomColors = None
        self.immediate = True

    def addWidgets(self):
        self.propLayout = propLayout = QtGui.QVBoxLayout()
        self.strokeGroup = QtGui.QGroupBox("Stroke")
        self.strokeGroup.setCheckable(True)
        self.strokeGroup.setChecked(False)
        self.strokeGroup.toggled.connect(self.setStroke_cb)
        strokeLayout = QtGui.QGridLayout()
        self.strokeGroup.setLayout(strokeLayout)
        # line width label and spinbox
        widthLab = QtGui.QLabel("width:")
        self.widthBox = QtGui.QSpinBox()
        self.widthBox.setRange(1, 10)
        self.widthBox.setValue(2)
        self.widthBox.valueChanged[int].connect(self.setLineWidth_cb)
        strokeLayout.addWidget(widthLab, 0, 0)
        strokeLayout.addWidget(self.widthBox, 0, 1)
        strokeLayout.setSpacing(10)

        # COLOR
        popupButton = QtGui.QPushButton("Color:")
        menu = QtGui.QMenu(self)
        colorsIcon = QtGui.QIcon(os.path.join(PMVICONPATH, 'colorChooser24.png'))
        custom = QtGui.QAction(colorsIcon, "Custom", self)
        custom.setIconVisibleInMenu(True)
        self.connect(custom, QtCore.SIGNAL("triggered()"), self.showCustomColorDialog_cb)
        menu.addAction(custom)
        action = menu.addAction("Same as atom")
        cb = CallbackFunction(self.backFaceColorModeChanged.emit, "sameAsFront")
        self.connect(action, QtCore.SIGNAL('triggered()'), cb)
        popupButton.setMenu(menu)
        strokeLayout.addWidget(popupButton, 1, 0)
        
        self.strokeGroup.adjustSize()
        propLayout.setContentsMargins(0,0,0,0)
        propLayout.addWidget(self.strokeGroup)
        self.setLayout(propLayout)
            
    def setStroke_cb(self, val=None):
        if not self.immediate: return
        width = self.widthBox.value()
        self.toggleStroke.emit(val, width)

    def setLineWidth_cb(self, val):
        if self.immediate:
            self.lineWidthChanged.emit(val)

    def showCustomColorDialog_cb(self):
        self.colorDialog = colorDialog = QtGui.QColorDialog(QtCore.Qt.green, self)
        colorDialog.currentColorChanged.connect(self.colorChanged_cb)
        colorDialog.finished.connect(self.colorDialogClosed_cb)
        colorDialog.colorSelected.connect(self.setCustomColor_cb)
        colorDialog.open()
        
    def colorDialogClosed_cb(self, val):
        #val == 0 (Cancel)  val == 1 (OK)
        if val == 0:
            # need to restore saved colors
            self.backFaceColorChanged.emit(0., 0., 0., 0., False, True)
            
    def colorChanged_cb(self, val):
        # save color of the geometry (if it has not been saved)
        # set the new color
        #print "colorChanged_cb", val
        rgb = val.getRgbF()
        self.backFaceColorChanged.emit(rgb[0], rgb[1], rgb[2], rgb[3], True, False)
        
    def setCustomColor_cb(self, val):
        #print "setCustomColor_cb", val
        rgb = val.getRgbF()
        self.backFaceColorChanged.emit(rgb[0], rgb[1], rgb[2], rgb[3], False, False)

    def saveWidgetsValues(self):
        kw = self.initialFormValues = {}
        kw['stroke'] = self.strokeGroup.isChecked()
        kw['width'] = self.widthBox.value()

    def restoreWidgetsValues(self):
        if len(self.initialFormValues):
            self.immediate = False
            kw = self.initialFormValues
            self.strokeGroup.setChecked(kw['stroke'])
            self.widthBox.setValue(kw['width'])
            self.immediate = True

class SBModelGlobalParams(QtGui.QWidget):
    """User input form for setting the global parameters of
    stick & balls geometry: quality, bond order parameters
   (double & triple bond separation and aromatic arck width)stroke """

    #SIGNALS
    qualityChanged = QtCore.Signal(int)
    softColorBoundariesToggled=QtCore.Signal(bool)
    # +  BondOrderWidgets signals
    # + Stroke widget signals

    def __init__(self, gui, geoms=[], parent=None):
        super(SBModelGlobalParams, self).__init__()
        self.geoms = geoms
        #self.parent = parent
        self.immediate = True
        self.gui = gui
        self.addWidgets()

    def addWidgets(self):
        self.layout = propLayout = QtGui.QVBoxLayout()
        # Quality
        qualLayout = QtGui.QHBoxLayout()
        qualLabel = QtGui.QLabel("quality:")
        self.qualBox  = QtGui.QSpinBox()
        self.qualBox.setRange(0, 5)
        self.qualBox.setSingleStep(1)
        self.qualBox.valueChanged.connect(self.setQuality_cb)
        qualLayout.addWidget(qualLabel)
        qualLayout.addWidget(self.qualBox)
        propLayout.addLayout(qualLayout)

        #Bond Order Widgets
        bolabels={"db": "double bond separation:", "tb":"triple bond separation:", "aromatic":"aromatic sphere rad:"}
        vals = {"db":("auto", [0.01, 10.0]), "tb":("auto", [0.01, 10.0]), "aromatic":( 0.15, [0., 10.])}
        self.bondOrderGroup = BondOrderWidgets(labels=bolabels, values=vals)
        propLayout.addWidget(self.bondOrderGroup)#.displayBO)
        # Stroke widget
        self.strokeWidget = StrokeWidgets(self.gui)
        propLayout.addWidget(self.strokeWidget)
        self.softBoundaries = QtGui.QCheckBox("soft color boundaries")
        self.softBoundaries.toggled.connect(self.colorBoundaries_cb)
        propLayout.addWidget(self.softBoundaries)
        self.setLayout(propLayout)
        #if self.parent:
        #    self.parent.setLayout(propLayout)

    def setQuality_cb(self, val):
        if not self.immediate: return
        self.qualityChanged.emit(val)

    def saveWidgetsValues(self):
        # Form widgets values
        initval = self.initialFormValues={}
        initval['softColorBoundaries'] = self.softBoundaries.isChecked()
        initval.update(self.bondOrderGroup.getWidgetsValues())
        self.strokeWidget.saveWidgetsValues()

    def restoreWidgetsValues(self):
        # restore initial form values
        self.immediate = False
        kw = self.initialFormValues
        bo = {}
        bo['displayBondOrder'] = kw['displayBondOrder']
        bo['doubleBondSep'] = kw['doubleBondSep']
        bo['tripleBondSep'] = kw['tripleBondSep']
        bo['aromaticLinewidth'] = kw['aromaticLinewidth']
        self.softBoundaries.setChecked(kw['softColorBoundaries'])
        self.immediate = True
        self.bondOrderGroup.immediate = False
        self.bondOrderGroup.setWidgetsValues(bo)
        self.bondOrderGroup.immediate = True
        self.strokeWidget.restoreWidgetsValues()

    def colorBoundaries_cb(self, val):
        if self.immediate:
            self.softColorBoundariesToggled.emit(val)

class MsmsSurfaceParams(QtGui.QWidget):
    """ """

    #SIGNALS
    qualityChanged = QtCore.Signal(int)
    probeRadiusChanged = QtCore.Signal(float)
    computeSES = QtCore.Signal(float, int)
    computeSAS = QtCore.Signal(int)
    softColorBoundariesToggled=QtCore.Signal(bool)
    #+ Stroke widget signals

    def __init__(self, gui, geoms=[], parent=None):
        super(MsmsSurfaceParams, self).__init__()
        self.geoms = geoms
        self.immediate = True
        self.gui = gui
        self.addWidgets()

    def addWidgets(self):
        self.layout = propLayout = QtGui.QVBoxLayout()
        surfTypeBox = QtGui.QGroupBox("Surface Type")
        self.sesRadioB = radio1 = QtGui.QRadioButton("Solvent Excluded")
        self.sasRadioB = radio2 = QtGui.QRadioButton("Solvent Accessible")
        radio1.setChecked(True)
        radio1.toggled.connect(self.solventExcluded_cb)
        radio2.toggled.connect(self.solventAccessible_cb)
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(radio1)
        vbox.addWidget(radio2)
        #vbox.addStretch(1)
        surfTypeBox.setLayout(vbox)
        surfTypeBox.adjustSize()
        propLayout.setContentsMargins(0,0,0,0)
        propLayout.addWidget(surfTypeBox)
        #self.stackWidget = QtGui.QStackedWidget()
        # Prob Rad
        probeRadLayout = QtGui.QHBoxLayout()
        probeRadLabel = QtGui.QLabel("Prob Radius:")
        self.probeRadBox = QtGui.QDoubleSpinBox()
        self.probeRadBox.setRange(0., 10.)
        self.probeRadBox.setSingleStep(0.1)
        self.probeRadBox.setValue(1.5)
        self.probeRadBox.valueChanged.connect(self.setProbRad_cb)
        probeRadLayout.addWidget(probeRadLabel)
        probeRadLayout.addWidget(self.probeRadBox)
        propLayout.addLayout(probeRadLayout)
        # Quality
        qualLayout = QtGui.QHBoxLayout()
        qualLabel = QtGui.QLabel("quality:")
        self.qualBox  = QtGui.QSpinBox()
        self.qualBox.setValue(6)
        self.qualBox.setRange(0, 10)
        self.qualBox.setSingleStep(1)
        
        self.qualBox.valueChanged.connect(self.setQuality_cb)
        qualLayout.addWidget(qualLabel)
        qualLayout.addWidget(self.qualBox)
        propLayout.addLayout(qualLayout)
        # Stroke
        self.strokeWidget = StrokeWidgets(self.gui)
        propLayout.addWidget(self.strokeWidget)
        # Soft color boundaries checkbutton
        self.softBoundaries = QtGui.QCheckBox("soft color boundaries")
        self.softBoundaries.toggled.connect(self.colorBoundaries_cb)
        self.softBoundaries.setChecked(True)
        propLayout.addWidget(self.softBoundaries)
        propLayout.addStretch(1)
        self.setLayout(propLayout)

    def solventExcluded_cb(self, val):
        #print "solventExcluded_cb", val
        self.probeRadBox.setEnabled(val)
        if val and self.immediate:
            rad = self.probeRadBox.value()
            qual = self.qualBox.value()
            self.computeSES.emit(rad, qual)

    def solventAccessible_cb(self, val):
        #print "solventAccessible_cb", val
        if val and self.immediate:
            qual = self.qualBox.value()
            self.computeSAS.emit(qual)

    def setProbRad_cb(self, val):
        if self.immediate:
            self.probeRadiusChanged.emit(val)
    
    def setQuality_cb(self, val):
        if self.immediate:
            self.qualityChanged.emit(val)

    def restoreWidgetsValues(self):
        self.immediate = False
        kw = self.initialFormValues
        self.strokeWidget.restoreWidgetsValues()
        self.probeRadBox.setValue(kw['probeRadius'])
        self.qualBox.setValue(kw['quality'])
        if kw['sesChecked']:
            self.sesRadioB.setChecked(True)
        else:
            self.sasRadioB.setChecked(True)
        self.softBoundaries.setChecked(kw['softColorBoundaries'])
        self.immediate = True

    def saveWidgetsValues(self):
        # Form widgets values
        initval = self.initialFormValues={}
        self.strokeWidget.saveWidgetsValues()
        initval['probeRadius'] = self.probeRadBox.value()
        initval['quality'] = self.qualBox.value()
        initval['sesChecked'] = self.sesRadioB.isChecked()
        initval['softColorBoundaries'] = self.softBoundaries.isChecked()

    def colorBoundaries_cb(self, val):
        if self.immediate:
            self.softColorBoundariesToggled.emit(val)


class CartoonModelGlobalParams(QtGui.QWidget):
    """ """

    #SIGNALS
    qualityChanged = QtCore.Signal(int)
    shapeChanged = QtCore.Signal(str)
    scaleChanged = QtCore.Signal(float)
    
    def __init__(self, gui, geoms=[], parent=None):
        super(CartoonModelGlobalParams, self).__init__()
        self.geoms = geoms
        self.immediate = True
        self.gui = gui
        self.qualityDict = {"low":3, "medium":10, "high":20}
        self.addWidgets()

    def addWidgets(self):
        self.layout = propLayout = QtGui.QVBoxLayout()
        #quality combo box
        qualLayout = QtGui.QHBoxLayout()
        qualLabel = QtGui.QLabel("quality")
        self.qualityBox = QtGui.QComboBox()
        for item in ["low", "medium", "high"]:
           self.qualityBox.addItem(item)
        self.qualityBox.setCurrentIndex(0)
        self.qualityBox.currentIndexChanged.connect(self.qualityChanged_cb)
        qualLayout.addWidget(qualLabel)
        qualLayout.addWidget(self.qualityBox)
        self.layout.addLayout(qualLayout)
        # Shape group with "default" and "cylinder" radiobuttons
        shapeBox = QtGui.QGroupBox("Shape")
        self.defaultRadioB = radio1 = QtGui.QRadioButton("default")
        self.cylinderRadioB = radio2 = QtGui.QRadioButton("cylinder")
        radio1.setChecked(True)
        radio1.toggled.connect(self.defaultShape_cb)
        radio2.toggled.connect(self.cylinderShape_cb)
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(radio1)
        vbox.addWidget(radio2)
        vbox.addStretch(1)
        shapeBox.setLayout(vbox)
        shapeBox.adjustSize()
        propLayout.addWidget(shapeBox)
        # scale spin box
        scaleLayout = QtGui.QHBoxLayout()
        scaleLabel = QtGui.QLabel("scale")
        self.scaleSpinBox = scaleSpinBox = QtGui.QDoubleSpinBox()
        scaleSpinBox.setRange(0.1, 5.0)
        scaleSpinBox.setSingleStep(0.1)
        scaleSpinBox.setValue(1.0)
        scaleSpinBox.valueChanged.connect(self.setScale_cb)
        scaleLayout.addWidget(scaleLabel)
        scaleLayout.addWidget(scaleSpinBox)
        propLayout.addLayout(scaleLayout)
        # Stroke
        self.strokeWidget = StrokeWidgets(self.gui)
        propLayout.addWidget(self.strokeWidget)
        self.setLayout(propLayout)

    def defaultShape_cb(self, val):
        if val and self.immediate:
            self.shapeChanged.emit("default")
        
    def cylinderShape_cb(self, val):
        if val and self.immediate:
            self.shapeChanged.emit("cylinder")

    def qualityChanged_cb(self, val):
        if not self.immediate: return
        val = self.qualityBox.currentText()
        self.qualityChanged.emit(self.qualityDict[str(val)])

    def setScale_cb(self, val):
        if self.immediate:
            self.scaleChanged.emit(val)
            
    def saveWidgetsValues(self):
        # Form widgets values
        initval = self.initialFormValues={}
        initval['quality'] = str(self.qualityBox.currentText())
        initval['scale'] = self.scaleSpinBox.value()
        initval['defaultShape'] = self.defaultRadioB.isChecked()
        self.strokeWidget.saveWidgetsValues()

    def restoreWidgetsValues(self):
        self.immediate = False
        kw = self.initialFormValues
        ind = self.qualityBox.findText(kw['quality'])
        self.qualityBox.setCurrentIndex(ind)
        if kw['defaultShape']:
            self.defaultRadioB.setChecked(True)
        else:
            self.cylinderRadioB.setChecked(True)
        self.scaleSpinBox.setValue(kw['scale'])
        self.strokeWidget.restoreWidgetsValues()
        self.immediate = True
