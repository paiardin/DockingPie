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
# $Header: /mnt/raid/services/cvs/PmvApp/GUI/Qt/advancedParamGui.py,v 1.1.4.2 2017/10/17 01:17:49 annao Exp $
#
# $Id: advancedParamGui.py,v 1.1.4.2 2017/10/17 01:17:49 annao Exp $

from PySide import QtCore, QtGui
from mglutil.util.callback import CallbackFunction
import weakref, os
import numpy
from mglutil.util.packageFilePath import findFilePath
PMVICONPATH = findFilePath('Icons', 'PmvApp.GUI')
from mglutil.util.callback import CallbackFunction

class ColorsWidget(QtGui.QWidget):

    beforeColorChanged = QtCore.Signal()

    def __init__(self, gui, geomNames=[], parent=None):
        super(ColorsWidget, self).__init__(parent)
        #self.layout = layout = QtGui.QGridLayout()
        self.immediate = True
        self.savedColors = []
        self.app = app = gui.selection[0].getAtomGroup().getMolecule().app()
        self.gui = gui
        try:
            len(geomNames)
        except:
            geomNames = [geomNames]
        self.geomNames=geomNames
        self.colorCmds = {'atom type':app.colorByAtomType,
                          'molecules': app.colorByMolecules,
                          'chains': app.colorByChains,
                          'N to C rainbow': app.colorRainbow,
                          'rainbow chain': app.colorRainbowChain,
                          'sec structure': app.colorBySS,
                          'custom':app.customColor}
        self.colorCmdsNames = ['atom type', 'molecules', 'chains',
                          'N to C rainbow', 'rainbow chain',
                          'sec structure']#, 'custom']

        #self.layout = layout = QtGui.QVBoxLayout()
        self.layout = layout = QtGui.QHBoxLayout()
        self.COnly = QtGui.QCheckBox("'C' only")
        self.COnly.setChecked(False)
        popupButton = QtGui.QPushButton("Color by:")
        self.menu = menu = QtGui.QMenu(self)
        self.colorsIcon = QtGui.QIcon(os.path.join(PMVICONPATH, 'colorChooser24.png'))
        self.buildMenu()
        popupButton.setMenu(menu)
        layout.addWidget(popupButton)
        
        self.COnly.setToolTip('Apply coloring to Carbon atoms only')
        self.COnly.toggled.connect(self.buildMenu)
        layout.addWidget(self.COnly)
        self.setLayout(layout)

    def buildMenu(self):
        cOnly = self.COnly.isChecked()
        self.menu.clear()
        for name in  self.colorCmdsNames:
            if name == "atom type" and cOnly == True:
                continue
            cb = CallbackFunction(self.color_cb, self.geomNames, self.colorCmds[name])
            action = self.menu.addAction(name)
            self.connect(action, QtCore.SIGNAL('triggered()'), cb)
        colorsIcon = QtGui.QIcon(os.path.join(PMVICONPATH, 'colorChooser24.png'))
        custom = QtGui.QAction(colorsIcon, "custom", self)
        custom.setIconVisibleInMenu(True)
        self.connect(custom, QtCore.SIGNAL("triggered()"), self.custom_cb)
        self.menu.addAction(custom)
        
        
    def color_cb(self, geomNames, cmd):
        self.beforeColorChanged.emit()
        cOnly = self.COnly.isChecked()
        if geomNames[0] == "msms":
            for sel in self.gui.selection:
                names = []
                for name in self.gui.getSurfName(sel):
                    names.append("msms_%s"%name)
                cmd(sel, names,  carbonsOnly=cOnly)
        else:
            for sel in self.gui.selection:
                cmd(sel, geomNames,  carbonsOnly=cOnly)

    def custom_cb(self):
        self.colorDialog = colorDialog = QtGui.QColorDialog(QtCore.Qt.green, self)
        colorDialog.currentColorChanged.connect(self.colorChanged_cb)
        colorDialog.finished.connect(self.colorDialogClosed_cb)
        colorDialog.colorSelected.connect(self.setCustomColor_cb)
        self.savedColors = []
        colorDialog.open()
        
    def colorDialogClosed_cb(self, val):
        #val == 0 (Cancel)  val == 1 (OK)
        if val == 0:
            if len(self.savedColors):
                if hasattr(self.gui, "restoreGeomColor"):
                    self.gui.restoreGeomColor(self.savedColors)
            
    def colorChanged_cb(self, val):
        if not len(self.savedColors):
            if hasattr (self.gui, "saveCurrentGeomColor"):
                self.savedColors = self.gui.saveCurrentGeomColor()
        self.setCustomColor_cb(val)
        
    def setCustomColor_cb(self, val):
        self.beforeColorChanged.emit()
        #print "setCustomColor_cb", val
        if hasattr (self.gui, "setGeomColor"):
            rgb = val.getRgbF()
            self.gui.setGeomColor(rgb, carbonsOnly=self.COnly.isChecked())
        else:
            print "Custom color for " , self.gui, " is not implemented" 

class OpacityWidget(QtGui.QWidget):
    opacityChanged = QtCore.Signal(float)
    def __init__(self, parent=None):
        super(OpacityWidget, self).__init__(parent)
        self.layout = layout = QtGui.QGridLayout()#QHBoxLayout()
        self.immediate = True 
        lab1 = QtGui.QLabel("Opacity:")
        self.opacityBox = QtGui.QDoubleSpinBox()
        self.opacityBox.setRange(0.0, 1.0)
        self.opacityBox.setSingleStep(0.1)
        self.opacityBox.setValue(1.0)
        layout.addWidget(lab1, 0, 0)
        #layout.addWidget(self.opacityBox, alignment=QtCore.Qt.AlignLeft)
        layout.addWidget(self.opacityBox, 0, 1)# alignment=QtCore.Qt.AlignLeft)
        self.opacityBox.valueChanged.connect(self.setOpacity_cb)
        layout.setContentsMargins(0,0,0,0)
        self.setLayout(layout)
        self.saveWidgetsValues()

    def setOpacity_cb(self, val):
        if self.immediate:
            self.opacityChanged.emit(val)

    def saveWidgetsValues(self):
        # Form widgets values
        self.initialFormValues={'opacity': self.opacityBox.value()}

    def restoreWidgetsValues(self):
        # restore initial form values
        self.immediate = False
        self.opacityBox.setValue(self.initialFormValues['opacity'])
        self.immediate = True

class CustomizeWidget(QtGui.QWidget):
    
    def __init__(self, name, molSel, selName, perAtomTabs=["Properties", "Colors"],
                 perMoleculeTabs=[], parent=None):
        super(CustomizeWidget, self).__init__(parent)
        self.geomName = name
        self.parent = parent
        self.gui = None # will be a weakref to the main gui
        self.mainLayout = QtGui.QGridLayout()
        self.tabWidget = tw = QtGui.QTabWidget()
        self.selection = molSel
        for sel in self.selection:
            mol =  sel.getAtomGroup().getMolecule()
            if not mol._renderingProp.has_key(self.geomName):
                print "molecule %s has no global properties for geom name %s" %(mol.name, self.geomName)
        tw.addTab(QtGui.QWidget(), selName)
        tw.addTab(QtGui.QWidget(), "Entire Molecule")
        self.mainLayout.addWidget(tw)
        self.setLayout(self.mainLayout)
        tw.currentChanged.connect(self.onSetCurrentTab_cb)
        # create tab widget inside selName (per atom properties) tab
        self.perAtomPropLayout = QtGui.QVBoxLayout()
        self.perAtomPropTabWidget = QtGui.QTabWidget()
        for prop in perAtomTabs:
            self.perAtomPropTabWidget.addTab(QtGui.QWidget(), prop)
        self.perAtomPropLayout.addWidget(self.perAtomPropTabWidget)
        self.tabWidget.widget(0).setLayout(self.perAtomPropLayout)
        self.perAtomPropTabWidget.currentChanged.connect(self.perAtomPropTabChanged_cb)
        
        # create tab widget for global properties (entire molecule)
        entireMolTab = self.tabWidget.widget(1)
        self.globalPropTabWidget = None
        if len(perMoleculeTabs):
            self.globalPropLayout = QtGui.QVBoxLayout()
            self.globalPropTabWidget = QtGui.QTabWidget()
            for prop in perMoleculeTabs:
                self.globalPropTabWidget.addTab(QtGui.QWidget(), prop)
            self.globalPropTabWidget.currentChanged.connect(self.globPropTabChanged_cb)
            self.globalPropLayout.addWidget(self.globalPropTabWidget)
            entireMolTab.setLayout(self.globalPropLayout)
        self.modelGlobalWidgets = None
        self.colorWidgets = None
        self.opacityWidget = None
        self.geoms = []
        QtGui.QWidget.connect(self.parent, QtCore.SIGNAL("immediateStatusChanged(bool)"), self.setImmediateStatus_cb)

    def updateSelection(self, selection, selName):
        self.selection = selection
        self.tabWidget.setTabText(0, selName)
        
    def onSetCurrentTab_cb(self, index):
        #print "onSetCurrentTab", index
        #applyB = self.parent.buttons.button(self.parent.buttons.Apply)
        if index == 1:
            self.addPropPerMoleculeWidgets()
            #if applyB.isEnabled():
            #    applyB.setEnabled(False)
            #self.parent.immediateCheckBox.setEnabled(False)
        #else:
            #self.parent.immediateCheckBox.setEnabled(True)
            #val = self.parent.immediateCheckBox.isChecked()
            #applyB.setEnabled(not val)

    def perAtomPropTabChanged_cb(self, index):
        #print "perAtomPropTabChanged_cb", index
        #applyB = self.parent.buttons.button(self.parent.buttons.Apply)
        if self.perAtomPropTabWidget.tabText(index) == "Colors":
            self.addColorsWidgets()
            #if applyB.isEnabled():
                #applyB.setEnabled(False)
                #self.parent.immediateCheckBox.setEnabled(False)
        #else:
        #    self.parent.immediateCheckBox.setEnabled(True)
        #    val = self.parent.immediateCheckBox.isChecked()
        #    applyB.setEnabled(not val)

    def globPropTabChanged_cb(self, val=None):
        pass
            
    def setImmediateStatus_cb(self, val):
        self.immediate = val

    def saveInitValues(self):
        self.initialValues = {}
    
    def addPropPerMoleculeWidgets(self):
        pass

    def addColorsWidgets(self):
        pass
    
    def setOpacity(self, val):
        for sel in self.selection:
            mol =  sel.getAtomGroup().getMolecule()
            inds = sel.getIndices()
            if not self.initialValues.has_key(sel):
                self.initialValues[sel] = {}
            dataLabel = 'opacity_%s'%self.geomName
            if mol._ag.getData(dataLabel) is None:
                opacity = numpy.ones(len(mol._ag), 'f')
                if not self.initialValues[sel].has_key("opacity"):
                    self.initialValues[sel]["opacity"]=opacity.copy()
                opacity[inds] = val
                mol._ag.setData(dataLabel, opacity)
            else:
                if not self.initialValues[sel].has_key("opacity"):
                    self.initialValues[sel]["opacity"] = mol._ag._data[dataLabel].copy() 
                mol._ag._data[dataLabel][inds] = val
            self.displayCmd.refreshDisplay(mol)

    def beforeGeomColorChanged(self):
        self.saveCurrentGeomColor()

    def saveCurrentGeomColor(self):
        geomcol = []
        for sel in self.selection:
            mol =  sel.getAtomGroup().getMolecule()
            col = mol._colors[self.geomName].copy()
            inds = mol._ag.getData('colorsIndices_%s'%self.geomName).copy()
            geomcol.append([sel, col, inds])
            if not self.initialValues.has_key(sel):
                self.initialValues[sel] = {}
            if not self.initialValues[sel].has_key('colors'):
                self.initialValues[sel]['colors'] = [col, inds.copy()]
        return geomcol

    def restoreGeomColor(self, geomcol):
        #print "restoreGeomColor", self.geomName
        for sel, color, inds in geomcol:
            mol =  sel.getAtomGroup().getMolecule()
            mol._ag.setData('colorsIndices_%s'%self.geomName, inds)
            mol._colors[self.geomName] = color
            self.displayCmd.refreshDisplay(mol)
            
    def setGeomColor(self, rgb, carbonsOnly=False):
        #print "setGeomColor:", rgb
        for sel in self.selection:
            if carbonsOnly:
                sel = sel & sel.select("element C")
            mol = sel.getAtomGroup().getMolecule()
            inds = sel.getIndices()
            col = mol._colors[self.geomName]
            mol._ag._data['colorsIndices_%s'%self.geomName][inds] = len(col)
            mol._colors[self.geomName] = numpy.concatenate((col, numpy.array([rgb])))
            self.displayCmd.refreshDisplay(mol)

    def cancel(self):
        # restore initial values of the molecular selection and the widgets
        pass

    def getPropValues(self):
        pass

class IndexedGeomCustomizeWidget(CustomizeWidget):
    
    def __init__(self, name, molSel, selName, perAtomTabs=["Properties", "Colors"],
                 perMoleculeTabs=[], parent=None):
        super(IndexedGeomCustomizeWidget, self).__init__(name, molSel, selName, perAtomTabs=perAtomTabs, perMoleculeTabs=perMoleculeTabs, parent=parent)
        self.backFaceColor=[]
        self.initialValues = {}
        ### the following attributes should be set by the class that inherits from this class
        self.displayCmd = None 
        ###

    def updateGlobalProp(self, kw, geomName=None):
        if not geomName: geomName = self.geomName
        for sel in self.selection:
            mol =  sel.getAtomGroup().getMolecule()
            prop = mol._renderingProp[geomName]
            if not self.initialValues.has_key(sel):
                    self.initialValues[sel] = {}
            for propName, value in kw.items():
                if not self.initialValues[sel].has_key(propName):
                    try:
                        if isinstance(prop[propName], str) and self.displayCmd._renderToPmv.has_key(prop[propName]):
                            self.initialValues[sel][propName]= self.displayCmd._renderToPmv[prop[propName]]
                        else:
                            self.initialValues[sel][propName] = prop[propName]
                    except:
                        print "No rendering properties %s for geometry %s" %(propName, geomName)  
            self.displayCmd.updateModelGlobals(mol, prop, **kw)
            self.displayCmd.refreshDisplay(mol)

    def setLineWidth(self, val):
        kw = {"linewidth":val}
        self.updateGlobalProp(kw)

    def setStroke(self, val, strokeWidth):
        # callback of the Stroke widget (sets Stroke on/off)
        if val:
            kw = {'culling':'none', 'frontRendering':'solid',
                  'backRendering':'mesh', 'linewidth':strokeWidth}
            #Check if the initial backFaceColor has been saved. If it has not,
            # then the color has not been modified --> we set the initial stroke color
            # to black.
            setColor = True
            sel = self.selection[0]
            if self.initialValues.has_key(sel):
                if self.initialValues[sel].has_key("backfaceColor"):
                    setColor = False
            if setColor:
                # initial stroke color is black
                kw["backfaceColor"]=[0.0, 0.0, 0.0, 1.0]
        else:
            kw = {'culling':'back', 'frontRendering':'solid',
                  'backRendering':'solid', 'linewidth':strokeWidth}
        self.updateGlobalProp(kw)

    def setBackFaceMode(self, mode=""):
        # this is connected to the Stroke widget's "Color"-->"Same as atom" menu callback 
        mode = str(mode)
        if mode == "sameAsFront":
            kw = {"backfaceColor":mode}
            self.updateGlobalProp(kw)

    def setBackFaceColor(self, r, g, b, a, saveCurrentColor=False, restoreColor=False):
        # callback of the Stroke widget's  custom color dialog.
        #print "setBackFaceColor:", "saveCurrentColor:", saveCurrentColor, "restoreColor", restoreColor
        if restoreColor:
            if len(self.backFaceColor):
                #print "restoring backface color"
                for sel, color in  self.backFaceColor:
                    mol =  sel.getAtomGroup().getMolecule()
                    prop = mol._renderingProp[self.geomName]
                    self.displayCmd.updateModelGlobals(mol, prop,
                              backfaceColor=color)
                    self.displayCmd.refreshDisplay(mol)
            self.backFaceColor=[]
            return
        if saveCurrentColor and len(self.backFaceColor) == 0:
            #print "saving back color"
            for sel in self.selection:
                mol =  sel.getAtomGroup().getMolecule()
                self.backFaceColor.append([sel, mol._renderingProp[self.geomName]['backfaceColor']])
        elif saveCurrentColor==False and restoreColor==False:
            self.backFaceColor = []
        kw = {"backfaceColor":[r,g,b,a]}
        #print "updateGlobalProp: backfaceColor"
        self.updateGlobalProp(kw)

    def setQuality(self, val):
        kw = {"quality":val}
        self.updateGlobalProp(kw)

    def setSoftColorBoundaries(self, val):
        kw = {"sharpColorBoundaries": not val}
        self.updateGlobalProp(kw)

    
class LinesCustomizeWidget(CustomizeWidget):
    
    def __init__(self,name, molSel, selName, parent=None):
        super(LinesCustomizeWidget, self).__init__(name, molSel, selName, parent=parent)
        self.displayCmd = self.selection[0].getAtomGroup().getMolecule().app().displayLines
        self.addPropPerAtomWidgets()
        self.modelGlobalWidgets = None
        self.colorWidgets = None
        self.saveInitValues()
        self.immediate = True
        
    def  addPropPerAtomWidgets(self):
        """Adds widgets to the Properties per Atom tab"""
        # Properties tab:
        selTab = self.perAtomPropTabWidget.widget(0)
        self.propLayout = propLayout = QtGui.QVBoxLayout()
        # get initial values of line properties
        self.stippleLineBox = QtGui.QCheckBox("stipple lines")
        self.stippleLineBox.setChecked(False)
        self.stippleLineBox.toggled.connect(self.toggleStippleLines_cb)
        propLayout.addWidget(self.stippleLineBox)
        propLayout.addStretch(1)
        selTab.setLayout(propLayout)

    def toggleStippleLines_cb(self, val):
        #if not self.parent.immediateCheckBox.isChecked():
        if not self.immediate:
            return
        for sel in self.selection:
            mol =  sel.getAtomGroup().getMolecule()
            mol._ag._flags['stippled'][sel.getIndices()] = val
            mol.app().displayLines.refreshDisplay(mol)

    def saveInitValues(self):
        # line properties of the selection
        self.initialValues = {}
        self.initialFormValues={'stippleLines':self.stippleLineBox.isChecked()}
        for sel in self.selection:
            mol =  sel.getAtomGroup().getMolecule()
            self.initialValues[sel] = {'stippleLines': mol._ag._flags['stippled'][sel.getIndices()][:]}
        if self.opacityWidget:
            self.opacityWidget.saveWidgetsValues()
        if self.modelGlobalWidgets:
            self.modelGlobalWidgets.saveWidgetsValues()

    def getPropValues(self):
        ind = self.tabWidget.currentIndex()
        if ind <0 : return
        kw = {}
        if ind == 0 : # first tab (properties per atom)
            stippleLines = self.stippleLineBox.isChecked()
            kw['stippleLines'] = stippleLines
                
    def cancel(self):
        # restore initial values of the geometry's parameters  and the widgets
        stippleLines = self.initialFormValues.get('stippleLines')
        if stippleLines is not None:
            self.immediate = False
            self.stippleLineBox.setChecked(stippleLines)
            self.immediate = True
        if self.modelGlobalWidgets:
            self.modelGlobalWidgets.restoreWidgetsValues()
        if self.opacityWidget:
            self.opacityWidget.restoreWidgetsValues()
        # restore initial values of the molecule
        for sel, molkw in self.initialValues.items():
            self.setGeomPropValues(molkw, sel)

    def setGeomPropValues(self, propValues, sel):
        kw = {}
        kw.update(propValues)
        stippleLines = kw.get('stippleLines', None)
        mol =  sel.getAtomGroup().getMolecule()
        if stippleLines is not None:
            kw.pop('stippleLines')
            mol._ag._flags['stippled'][sel.getIndices()] = stippleLines
        opacity = kw.get('opacity', None)
        if opacity is not None:
            kw.pop("opacity")
            dataLabel = 'opacity_%s'%self.geomName
            mol._ag.setData(dataLabel, opacity)
        if kw.has_key("colors"):
            colors, inds = kw.pop("colors")
            mol._ag.setData('colorsIndices_lines', inds)
            mol._colors['lines'] = colors
        if len(kw):
            mol.app().displayLines.updateModelGlobals(mol, **kw)
        mol.app().displayLines.refreshDisplay(mol)

    def addPropPerMoleculeWidgets(self):
        from DejaVu2.Qt.geomParamGUI import LinesModelGlobalParams
        if not self.modelGlobalWidgets:
            self.modelGlobalWidgets = glw = LinesModelGlobalParams()
            layout = QtGui.QVBoxLayout()
            layout.addWidget(glw)
            self.tabWidget.widget(1).setLayout(layout)
            #connect globalpropwidgets signals:
            bo = glw.bondOrderGroup
            bo.toggleBO.connect(self.toggleDisplayBO)
            bo.doubleBondSeparation.connect(self.setDBSeparation)
            bo.tripleBondSeparation.connect(self.setTBSeparation)
            bo.aromaticValueChanged.connect(self.setAromaticLineWidth)
            glw.lineWidthChanged.connect(self.setLineWidth)
            glw.stippleLengthChanged.connect(self.setStippleLength)
            glw.stippleSpaceChanged.connect(self.setStippleSpace)
            glw.softColorBoundariesToggled.connect(self.setSoftColorBoundaries)
            self.modelGlobalWidgets.saveWidgetsValues()

    def updateGlobalProp(self, propName, value):
        for sel in self.selection:
            mol =  sel.getAtomGroup().getMolecule()
            if not self.initialValues.has_key(mol):
                self.initialValues[mol] = {}
            if not self.initialValues[mol].has_key(propName):
                self.initialValues[mol][propName] = mol._renderingProp['lines'][propName]
            self.displayCmd.updateModelGlobals(mol, **{propName:value})
            self.displayCmd.refreshDisplay(mol)

    def setSoftColorBoundaries(self, val):
        self.updateGlobalProp("sharpColorBoundaries", not val)

    def setLineWidth(self, val):
        self.updateGlobalProp("linewidth", val)
        
    def toggleDisplayBO(self, val):
        #print "toggleDisplayBO:", val
        self.updateGlobalProp("displayBondOrder", val)
                
    def setStippleLength(self, val):
        #print "Stipple Length:", val
        self.updateGlobalProp("stippleLength", val)

    def setStippleSpace(self, val):
        #print "Stipple Space:", val
        self.updateGlobalProp("stippleSpace", val)

    def setDBSeparation(self, val):
        #print "DBSeparation:", val
        self.updateGlobalProp("doubleBondSep", val)

    def setTBSeparation(self, val):
        #print "setTBSeparation:", val
        self.updateGlobalProp("tripleBondSep", val)

    def setAromaticLineWidth(self, val):
        #print "setAromaticLineWidth", val, type(val)
        self.updateGlobalProp('aromaticLinewidth', int(val))

    def addColorsWidgets(self):
        if self.colorWidgets: return
        self.colorWidgets = ColorsWidget(self, geomNames=["lines"])
        self.colorWidgets.beforeColorChanged.connect(self.beforeGeomColorChanged)
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.colorWidgets)
        self.opacityWidget = OpacityWidget()
        layout.addWidget(self.opacityWidget)
        self.opacityWidget.opacityChanged.connect(self.setOpacity)
        layout.addStretch(1)
        self.perAtomPropTabWidget.widget(1).setLayout(layout)


class RadiusWidget(QtGui.QWidget):
    radiusChanged = QtCore.Signal(float, str)
    
    def __init__(self, propList, parent=None, defaultSize=1.0):
        super(RadiusWidget, self).__init__()
        self.immediate = True
        self.radLayout = radLayout = QtGui.QVBoxLayout()
        sizeLayout =  QtGui.QHBoxLayout()
        self.propList = propList
        self.defaultSize = defaultSize
        lab = QtGui.QLabel("Size:")
        self.sizeSpinBox = sizeSpinBox = QtGui.QDoubleSpinBox()
        sizeSpinBox.setRange(0.0, 10.0)
        sizeSpinBox.setSingleStep(0.1)
        sizeSpinBox.setValue(self.defaultSize)
        sizeSpinBox.valueChanged.connect(self.radSizeChanged_cb)
        self.byPropGroup =  QtGui.QGroupBox("scale by property:")
        self.byPropGroup.setCheckable(True)
        self.byPropGroup.toggled.connect(self.byPropGroupToggled_cb)
        vbox = QtGui.QVBoxLayout()
        self.propComboBox = QtGui.QComboBox()
        for item in propList:
           self.propComboBox.addItem(item)
        self.propComboBox.currentIndexChanged.connect(self.propertyChanged_cb)
        vbox.addWidget(self.propComboBox)
        self.byPropGroup.setLayout(vbox)
        sizeLayout.addWidget(lab)
        sizeLayout.addWidget(sizeSpinBox)
        radLayout.addLayout(sizeLayout)
        radLayout.addWidget(self.byPropGroup)
        if parent:
            parent.addLayout(radLayout)
        self.initialFormValues = self.getWidgetsValues()

    def getWidgetsValues(self):
        return {'radius': self.sizeSpinBox.value(),
                'scaleByProp': self.byPropGroup.isChecked(),
                'property': str(self.propComboBox.currentText())}

    def setWidgetsValues(self, kw):
        self.immediate = False
        self.sizeSpinBox.setValue(kw['radius'])
        self.byPropGroup.setChecked(kw['scaleByProp']),
        ind = self.propComboBox.findText(kw['property'])
        self.propComboBox.setCurrentIndex(ind)
        self.immediate = True

    def propertyChanged_cb(self, val):
        #print "propertyChanged_cb", val
        if not self.immediate: return
        prop = str(self.propComboBox.currentText())
        size = self.sizeSpinBox.value()
        self.radiusChanged.emit(size, str(prop))

    def radSizeChanged_cb(self, rad):
        if not self.immediate: return
        if self.byPropGroup.isChecked():
            prop = str(self.propComboBox.currentText())
            self.radiusChanged.emit(rad, prop)
        else:
            self.radiusChanged.emit(rad, "")

    def byPropGroupToggled_cb(self, val):
        #print "byPropGroupToggled_cb:", val
        if not self.immediate: return
        if val:
            prop = str(self.propComboBox.currentText())
            rad = self.sizeSpinBox.value()
            self.radiusChanged.emit(rad, prop)
        else:
            rad = self.sizeSpinBox.value()
            self.radiusChanged.emit(rad, "")


class CpkCustomizeWidget(IndexedGeomCustomizeWidget):
    
    def __init__(self,name, molSel, selName, parent=None):
        super(CpkCustomizeWidget, self).__init__(name, molSel, selName, parent=parent)
        self.displayCmd = self.selection[0].getAtomGroup().getMolecule().app().displayCPK
        self.addPropPerAtomWidgets()
        self.modelGlobalWidgets = None
        self.colorWidgets = None
        self.backFaceColor = []
        self.saveInitValues()
        self.immediate = True
        
    def  addPropPerAtomWidgets(self):
        """Adds widgets to the Properties per Atom tab"""
        # Properties tab:
        selTab = self.perAtomPropTabWidget.widget(0)
        self.propLayout = propLayout = QtGui.QVBoxLayout()
        prop = ['vdw']
        for sel in self.selection:
            if 'charge' not in prop and sel._ag.getCharges() is not None:
                prop.append('charge')
            if 'anisotropic factor' not in prop and sel._ag.getAnisous() is not None:
                prop.append('anisotropic factor') 
            if 'temperature factor' not in prop and sel._ag.getBetas() is not None:
                prop.append('temperature factor')
        self.radiusWidget = RadiusWidget(prop, propLayout)
        self.radiusWidget.radiusChanged.connect(self.setCPKRadius)
        propLayout.addStretch(1)
        selTab.setLayout(propLayout)
        
    def saveInitValues(self):
        self.initialValues = {}
        self.initialFormValues= self.radiusWidget.getWidgetsValues()
        if self.modelGlobalWidgets:
            self.modelGlobalWidgets.saveWidgetsValues()
        if self.opacityWidget:
            self.opacityWidget.saveWidgetsValues()

    def setCPKRadius(self, size, prop=""):
        for sel in self.selection:
            if len(prop):
                rad = self.computeCPKRadius(size, str(prop), sel)
            else:
                rad = size
            mol =  sel.getAtomGroup().getMolecule()
            # save initial cpk_radius value (if it has not been saved yet) 
            if not self.initialValues.has_key(sel):
                self.initialValues[sel] = {}
            if not self.initialValues[sel].has_key("cpk_radius"):
                cpkRad = mol._ag.getData("cpk_radius")
                self.initialValues[sel]['cpk_radius'] = cpkRad[sel.getIndices()]
            mol.app().displayCPK.updateModel(sel, radii=rad)
            mol.app().displayCPK.refreshDisplay(mol)
        
    def computeCPKRadius(self, size, prop, sel):
        mol = sel.getAtomGroup().getMolecule()
        inds = sel.getIndices()
        radii = None
        if prop == "vdw":
            radii = mol._ag.getRadii()
        elif prop == 'charge':   
            radii = mol._ag.getCharges()
        elif prop == 'anisotropic factor':
            radii = mol._ag.getAnisous()
        elif prop == 'temperature factor':
            radii = mol._ag.getBetas()
        if radii is not None:
            radii = radii*size
            return radii[inds]

    def setQuality(self, val):
        kw = {"quality":val}
        self.updateGlobalProp(kw)
    
    def getPropValues(self):
        ind = self.tabWidget.currentIndex()
        if ind <0 : return
        if ind == 0 : # Properties tab
           kw = {}
           vals = self.radiusWidget.getWidgetsValues()
           size = vals['radius']
           if vals['scaleByProp']:
               prop = str(vals['property'])
               self.setCPKRadius(size, prop)
           else:
               self.setCPKRadius(size, "")

    def  addPropPerMoleculeWidgets(self):
        from DejaVu2.Qt.geomParamGUI import StrokeWidgets
        if not self.modelGlobalWidgets:
            self.modelGlobalWidgets = glw =StrokeWidgets(self)
            layout = QtGui.QVBoxLayout()
            layout.addWidget(glw)
            self.tabWidget.widget(1).setLayout(layout)
            glw.toggleStroke.connect(self.setStroke)
            glw.backFaceColorModeChanged.connect(self.setBackFaceMode)
            glw.backFaceColorChanged.connect(self.setBackFaceColor)
            glw.lineWidthChanged.connect(self.setLineWidth)
            layout.addStretch(1)
            self.modelGlobalWidgets.saveWidgetsValues()

    def addColorsWidgets(self):
        if self.colorWidgets: return
        self.colorWidgets = ColorsWidget(self, geomNames=["cpk"])
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.colorWidgets)
        self.colorWidgets.beforeColorChanged.connect(self.beforeGeomColorChanged)
        self.opacityWidget = OpacityWidget()
        layout.addWidget(self.opacityWidget)
        self.opacityWidget.opacityChanged.connect(self.setOpacity)
        layout.addStretch(1)
        self.perAtomPropTabWidget.widget(1).setLayout(layout)

    def setImmediateStatus_cb(self, val):
        self.immediate = val
        self.radiusWidget.immediate = val
        #if self.modelGlobalWidgets:
         #   self.modelGlobalWidgets.immediate = val

    def cancel(self):
        # restore initial values of the geometry's parameters  and the widgets
        if self.modelGlobalWidgets:
            self.modelGlobalWidgets.restoreWidgetsValues()
        if self.opacityWidget:
            self.opacityWidget.restoreWidgetsValues()

        for sel in self.selection:
            mol =  sel.getAtomGroup().getMolecule()
            kw = self.initialValues.get(sel, {})
            display = False
            if len(kw): display = True
            if kw.has_key('cpk_radius'):
                rad = kw.pop('cpk_radius')
                mol.app().displayCPK.updateModel(sel, radii=rad)
            if kw.has_key('opacity'):
                opacity = kw.pop("opacity")
                dataLabel = 'opacity_%s'%self.geomName
                mol._ag.setData(dataLabel, opacity)
            if kw.has_key("colors"):
                colors, inds = kw.pop("colors")
                mol._ag.setData('colorsIndices_cpk', inds)
                mol._colors['cpk'] = colors
            if len(kw):
                prop = mol._renderingProp['cpk']
                mol.app().displayCPK.updateModelGlobals(mol, prop, **kw)
            if display:    
                mol.app().displayCPK.refreshDisplay(mol)
        self.radiusWidget.setWidgetsValues(self.initialFormValues)

class SBCustomizeWidget(IndexedGeomCustomizeWidget):
    
    def __init__(self,name, molSel, selName, parent=None):
        super(SBCustomizeWidget, self).__init__(name, molSel, selName, parent=parent)
        self.displayCmd = self.selection[0].getAtomGroup().getMolecule().app().displaySB
        self.geomName = 'sb'
        self.addPropPerAtomWidgets()
        self.modelGlobalWidgets = None
        self.colorWidgets = None
        self.immediate = True #self.parent.immediateCheckBox.isChecked()
        self.backFaceColor = []
        self.sticksColor = []
        self.saveInitValues()

    def  addPropPerAtomWidgets(self):
        """Adds widgets to the Properties per Atom tab"""
        # Properties tab
        selTab = self.perAtomPropTabWidget.widget(0)
        # Balls widgets
        self.propLayout = propLayout = QtGui.QVBoxLayout()
        prop = ['vdw']
        for sel in self.selection:
            if 'charge' not in prop and sel._ag.getCharges() is not None:
                prop.append('charge')
            if 'anisotropic factor' not in prop and sel._ag.getAnisous() is not None:
                prop.append('anisotropic factor') 
            if 'temperature factor' not in prop and sel._ag.getBetas() is not None:
                prop.append('temperature factor')
        self.defaultBallRad = 0.3
        self.ballRadWidget = RadiusWidget(prop, defaultSize=self.defaultBallRad)
        self.ballRadWidget.byPropGroup.setChecked(False)
        self.ballRadWidget.sizeSpinBox.setSingleStep(0.05)
        self.defaultBallRad
        self.ballRadWidget.radiusChanged.connect(self.setBallRadius)
        self.ballGroup = QtGui.QGroupBox("Balls")
        self.ballGroup.setLayout(self.ballRadWidget.radLayout)
        # Cylinders widgets
        self.cylGroup = QtGui.QGroupBox("Sticks")
        cylLayout = QtGui.QHBoxLayout()
        lab = QtGui.QLabel("Radius:")
        self.cylRadSpinBox = QtGui.QDoubleSpinBox()
        self.cylRadSpinBox.setRange(0.0, 10.0)
        self.cylRadSpinBox.setSingleStep(0.05)
        self.cylRadSpinBox.setValue(0.2)
        self.cylRadSpinBox.valueChanged.connect(self.setCylRadius_cb)
        cylLayout.addWidget(lab)
        cylLayout.addWidget(self.cylRadSpinBox)
        self.cylGroup.setLayout(cylLayout)
        propLayout.addWidget(self.ballGroup)
        propLayout.addWidget(self.cylGroup)
	selTab.setLayout(propLayout)

    def saveInitValues(self):
	self.initialValues = {}
	self.initialFormValues= self.ballRadWidget.getWidgetsValues()
	self.initialFormValues['cylRadius'] = self.cylRadSpinBox.value()
        if self.modelGlobalWidgets:
	    self.modelGlobalWidgets.saveWidgetsValues()
        if self.opacityWidget:
            self.opacityWidget.saveWidgetsValues()

    def setBallRadius(self, size, prop=""):
        #import pdb; pdb.set_trace()
        cylRad = self.cylRadSpinBox.value()
        for sel in self.selection:
	    mol =  sel.getAtomGroup().getMolecule()
	    if not self.initialValues.has_key(sel):
                self.initialValues[sel] = {}
            if not self.initialValues[sel].has_key('balls_radius'):
		brad = mol._ag.getData("sb_ballsRadius")
		if brad is not None:
		    self.initialValues[sel]['balls_radius'] = brad[sel.getIndices()]
	    if len(prop):
                rad = self.computeBallRadius(size, str(prop), sel)
            else:
                rad = size

            mol.app().displaySB.updateModel(sel, ballsRadii=rad, cylRadii=cylRad)
            mol.app().displaySB.refreshDisplay(mol)
            
    def computeBallRadius(self, size, prop, sel):
        mol = sel.getAtomGroup().getMolecule()
        inds = sel.getIndices()
        if prop == "vdw":
            radii = mol._ag.getRadii()
        elif prop == 'charge':   
            radii = mol._ag.getCharges()
        elif prop == 'anisotropic factor':
            radii = mol._ag.getAnisous()
        elif prop == 'temperature factor':
            radii = mol._ag.getBetas()
        if radii is not None:
            radii = radii*size
            return radii[inds]

    def setCylRadius_cb(self, rad):
        if not self.immediate: return
        self.setCylRadius(rad)

    def setCylRadius(self, rad):
        for sel in self.selection:
	    mol =  sel.getAtomGroup().getMolecule()
	    if not self.initialValues.has_key(sel):
                self.initialValues[sel] = {}
            if not self.initialValues[sel].has_key('cyl_radius'):
	        crad = mol._ag._bondData.get("sb_cylRadius")
		if crad is not None:
		    bonds = sel.getBonds()
		    bnums = [mol._ag._bondIndex['%d %d'%(b[0],b[1])] for b in bonds[1]+bonds[2]+bonds[3]]
		    self.initialValues[sel]['cyl_radius'] = crad[bnums]

            ballRad = self.getBallRadFromWidgets(sel)
            mol.app().displaySB.updateModel(sel, cylRadii=rad, ballsRadii=ballRad)
            mol.app().displaySB.refreshDisplay(mol)

    def getPropValues(self):
        ind = self.tabWidget.currentIndex()
        if ind <0 : return
        if ind == 0 : # Properties tab
           vals = self.ballRadWidget.getWidgetsValues()
           size = vals['radius']
           if vals['scaleByProp']:
               prop = str(vals['property'])
               self.setBallRadius(size, prop)
           else:
               self.setBallRadius(size, "")
           cylRad = self.cylRadSpinBox.value()
           self.setCylRadius(cylRad)
                         
    def getBallRadFromWidgets(self, sel):
        vals = self.ballRadWidget.getWidgetsValues()
        size = vals['radius']
        if vals['scaleByProp']:
            prop = str(vals['property'])                  
            return self.computeBallRadius(size, str(prop), sel)
        else:
            return size

    def addPropPerMoleculeWidgets(self):
        from DejaVu2.Qt.geomParamGUI import SBModelGlobalParams
        if not self.modelGlobalWidgets:
            self.modelGlobalWidgets =SBModelGlobalParams(self)
            layout = QtGui.QVBoxLayout()
            layout.addWidget(self.modelGlobalWidgets)
            self.tabWidget.widget(1).setLayout(layout)
            self.modelGlobalWidgets.qualityChanged.connect(self.setQuality)
            self.modelGlobalWidgets.softColorBoundariesToggled.connect(self.setSoftColorBoundaries)
            bo = self.modelGlobalWidgets.bondOrderGroup
            bo.toggleBO.connect(self.toggleDisplayBO)
            bo.doubleBondSeparation.connect(self.setDBSeparation)
            bo.tripleBondSeparation.connect(self.setTBSeparation)
            bo.aromaticValueChanged.connect(self.setAromaticSphereRad)
            stw = self.modelGlobalWidgets.strokeWidget
            stw.toggleStroke.connect(self.setStroke)
            stw.backFaceColorModeChanged.connect(self.setBackFaceMode)
            stw.backFaceColorChanged.connect(self.setBackFaceColor)
            stw.lineWidthChanged.connect(self.setLineWidth)
            layout.addStretch(1)
            self.modelGlobalWidgets.saveWidgetsValues()

    def toggleDisplayBO(self, val):
        #print "toggleDisplayBO:", val
        self.updateGlobalProp({'displayBondOrder':val})

    def setDBSeparation(self, val):
        if val == 0.0: val = 'auto'
        self.updateGlobalProp({'doubleBondSep': val})

    def setTBSeparation(self, val):
        if val == 0.0: val = 'auto'
        self.updateGlobalProp({'tripleBondSep': val})

    def setAromaticSphereRad(self, rad):
        #print "setAromaticSphereRad", rad
        self.updateGlobalProp({'aromaticSphRad':rad})

    def addColorsWidgets(self):
        if self.colorWidgets: return
        layout = QtGui.QVBoxLayout()
        self.ballColorGroup = QtGui.QGroupBox("Sticks and Balls")
        self.colorWidgets = ColorsWidget(self, geomNames=['sb'])
        self.colorWidgets.beforeColorChanged.connect(self.beforeGeomColorChanged)
        ballLayout = QtGui.QVBoxLayout()
        ballLayout.addWidget(self.colorWidgets)
        self.ballColorGroup.setLayout(ballLayout)
        layout.addWidget(self.ballColorGroup)
        self.sticksColorGroup = QtGui.QGroupBox("Sticks")
        self.sticksColorGroup.setCheckable(True)
        self.sticksColorGroup.setChecked(False)
        #self.sticksColorGroup.toggled.connect()
        sticksColorLayout = QtGui.QHBoxLayout()
        stickColorIcon = QtGui.QIcon(os.path.join(PMVICONPATH, 'colorChooser24.png'))
        stickColorB = btn = QtGui.QPushButton(stickColorIcon, "")
        iconSize = QtCore.QSize(32,32)
        btn.setIconSize(iconSize)
        btn.setMaximumSize(QtCore.QSize(35, 35))
        btn.clicked.connect(self.showSticksCustomColorDialog_cb)
        btn.setToolTip("choose custom color")
        sticksColorLayout.addWidget(btn)
        self.SticksCOnly = QtGui.QCheckBox("'C' only")
        self.SticksCOnly.setChecked(False)
        self.SticksCOnly.setToolTip('Apply coloring to Carbon atoms only')
        sticksColorLayout.addWidget(self.SticksCOnly)
        
        self.sticksColorGroup.setLayout(sticksColorLayout)
        layout.addWidget(self.sticksColorGroup)
        self.opacityWidget = OpacityWidget()
        layout.addWidget(self.opacityWidget)
        self.opacityWidget.opacityChanged.connect(self.setOpacity)
        layout.addStretch(1)
        self.perAtomPropTabWidget.widget(1).setLayout(layout)

    def saveCurrentGeomColor(self):
        geomcol = []
        for sel in self.selection:
            mol =  sel.getAtomGroup().getMolecule()
            bcol = mol._colors['sb_balls'].copy()
            binds = mol._ag.getData('colorsIndices_sb_balls').copy()
            ccol = mol._colors['sb_cyl'].copy()
            cinds = mol._ag.getData('colorsIndices_sb_cyl').copy()
            geomcol.append([sel, bcol, binds, ccol, cinds])
            if not self.initialValues.has_key(sel):
                self.initialValues[sel] = {}
            if not self.initialValues[sel].has_key('colors'):
                self.initialValues[sel]['colors'] = [bcol, binds.copy(), ccol, cinds.copy()]
        return geomcol

    def restoreGeomColor(self, geomcol):
        for sel, bcolor, binds, ccolor, cinds in geomcol:
            mol =  sel.getAtomGroup().getMolecule()
            mol._ag.setData('colorsIndices_sb_balls', binds)
            mol._colors['sb_balls'] = bcolor
            mol._ag.setData('colorsIndices_sb_cyl', cinds)
            mol._colors['sb_cyl'] = ccolor
            mol.app().displaySB.refreshDisplay(mol)
            
    def setGeomColor(self, rgb, carbonsOnly=False):
        for sel in self.selection:
            if carbonsOnly:
                sel = sel & sel.select("element C")
            mol = sel.getAtomGroup().getMolecule()
            inds = sel.getIndices()
            bcol = mol._colors['sb_balls']
            mol._ag._data['colorsIndices_sb_balls'][inds] = len(bcol)
            mol._colors['sb_balls'] = numpy.concatenate((bcol, numpy.array([rgb])))
            ccol = mol._colors['sb_cyl']
            mol._ag._data['colorsIndices_sb_cyl'][inds] = len(ccol)
            mol._colors['sb_cyl'] = numpy.concatenate((ccol, numpy.array([rgb])))
            mol.app().displaySB.refreshDisplay(mol)

    def showSticksCustomColorDialog_cb(self):
        self.sticksColorDialog = colorDialog = QtGui.QColorDialog(QtCore.Qt.green, self)
        colorDialog.currentColorChanged.connect(self.sticksColorChanged_cb)
        colorDialog.finished.connect(self.sticksColorDialogClosed_cb)
        colorDialog.colorSelected.connect(self.setSticksColor_cb)
        colorDialog.open()
        
    def sticksColorDialogClosed_cb(self, val):
        #val == 0 (Cancel)  val == 1 (OK)
        if val == 0:
            # need to restore saved colors
            if len(self.sticksColor):
                #print "restoring sticks color"
                for sel, color, inds in  self.sticksColor:
                    mol =  sel.getAtomGroup().getMolecule()
                    mol._ag.setData('colorsIndices_sb_cyl', inds)
                    mol._colors['sb_cyl'] = color
                    mol.app().displaySB.refreshDisplay(mol)
            self.sticksColor=[]
            
    def sticksColorChanged_cb(self, val):
        # save color of the geometry (if it has not been saved)
        # set the new color
        rgb = val.getRgbF()
        if len(self.sticksColor) == 0:
            #print "saving sticks color"
            for sel in self.selection:
                mol =  sel.getAtomGroup().getMolecule()
                self.sticksColor.append([sel, mol._colors['sb_cyl'].copy(), mol._ag.getData('colorsIndices_sb_cyl').copy()])
        self.setSticksColor(rgb, carbonsOnly=self.SticksCOnly.isChecked())
        
    def setSticksColor_cb(self, val):
        #Called after Color dialog OK button is clicked
        #print "setCustomColor_cb", val
        rgb = val.getRgbF()
        self.sticksColor = []
        self.setSticksColor(rgb, carbonsOnly=self.SticksCOnly.isChecked())

    def setSticksColor(self, rgb, carbonsOnly=False):
        for sel in self.selection:
            if carbonsOnly:
                sel = sel & sel.select("element C")
            mol = sel.getAtomGroup().getMolecule()
            inds = sel.getIndices()
            ccol = mol._colors['sb_cyl']
            mol._ag._data['colorsIndices_sb_cyl'][inds] = len(ccol)
            mol._colors['sb_cyl'] = numpy.concatenate((ccol, numpy.array([rgb])))
            mol.app().displaySB.refreshDisplay(mol)
        
    def cancel(self):
        # restore initial values of the geometry's parameters and the widgets
        self.immediate = False
        if self.modelGlobalWidgets:
            self.modelGlobalWidgets.restoreWidgetsValues()
        if self.opacityWidget:
            self.opacityWidget.restoreWidgetsValues()
        for sel in self.selection:
            kw = self.initialValues.get(sel, {})
            display = len(kw)>0
            if not display: continue
            brad = kw.get('balls_radius')
            crad = kw.get('cyl_radius')
            mol =  sel.getAtomGroup().getMolecule()
            if brad is not None or crad is not None:
                mol.app().displaySB.updateModel(sel, ballsRadii=brad, cylRadii=crad)
            if kw.has_key('opacity'):
                opacity = kw.pop("opacity")
                dataLabel = 'opacity_%s'%self.geomName
                mol._ag.setData(dataLabel, opacity)
            if kw.has_key("colors"):
                bcol, binds, ccol, cinds = kw.pop("colors")
                mol._ag.setData('colorsIndices_sb_balls', binds)
                mol._colors['sb_balls'] = bcol
                mol._ag.setData('colorsIndices_sb_cyl', cinds)
                mol._colors['sb_cyl'] = ccol
            if self.modelGlobalWidgets:
                if brad is not None:
                    kw.pop('balls_radius')
                if crad is not None:
                    kw.pop('cyl_radius')
		if len(kw):
		    mol.app().displaySB.updateModelGlobals(mol, mol._renderingProp['sb'], **kw)
	    if display:
	        mol.app().displaySB.refreshDisplay(mol)
        self.ballRadWidget.setWidgetsValues(self.initialFormValues)
        self.cylRadSpinBox.setValue(self.initialFormValues.get('cylRadius', 0.2))
        self.immediate = True

    def setImmediateStatus_cb(self, val):
        self.immediate = val
        self.ballRadWidget.immediate = val
        #if self.modelGlobalWidgets:
        #    self.modelGlobalWidgets.immediate = val

class SurfaceTableWidget(QtGui.QWidget):

    def __init__(self, gui, selection, surfName, parent=None):
        super(SurfaceTableWidget, self).__init__(parent)
        self.immediate = True
        self.header = ["name", "show", "flip", "vol", "area"]
        ncolumns = len(self.header)
        self.surfName = surfName
        self.molSel = selection
        mol = selection.getAtomGroup().getMolecule()
        surf = mol._msmsData['msms'][surfName][mol._ag.getACSIndex()]
        self.gui = gui
        self.maxVolume = None
        self.minVolume = None
        self.layout = QtGui.QVBoxLayout()
        hbox1 = QtGui.QHBoxLayout()
        self.allVisibleB = QtGui.QPushButton("all visible")
        self.allVisibleB.clicked.connect(self.allVisible_cb)
        self.invertB = QtGui.QPushButton("invert")
        self.invertB.clicked.connect(self.invert_cb)
        hbox1.addWidget(self.allVisibleB)
        hbox1.addWidget(self.invertB)
        self.layout.addLayout(hbox1)
        # Table
        nrows = surf.sesr.nb
        self.table = table = QtGui.QTableWidget(nrows, ncolumns, self)
        table.setHorizontalHeaderLabels(self.header)
        table.horizontalHeader().setResizeMode(ncolumns-1, QtGui.QHeaderView.Stretch)
        table.setColumnWidth(0, 90)
        table.setColumnWidth(1, 40)
        table.setColumnWidth(2, 40)
        table.setColumnWidth(3, 70)
        table.setColumnWidth(4, 70)
        self.layout.addWidget(table)
        table.horizontalHeader().setSortIndicator(0, QtCore.Qt.AscendingOrder)
        
        # Cavities check button
        self.cavitiesCB = QtGui.QCheckBox("Cavities")
        self.cavitiesCB.setChecked(mol.app().userpref['Compute cavities by default']['value']=='yes')
        self.cavitiesCB.toggled.connect(self.cavities_cb)
        self.layout.addWidget(self.cavitiesCB, alignment=QtCore.Qt.AlignLeft)
        # Volume Filter
        volFilterLayout = QtGui.QGridLayout()
        self.volFilterGroup = QtGui.QGroupBox("volume filter")
        self.volFilterGroup.setCheckable(True)
        self.volFilterGroup.setChecked(False)
        self.volFilterGroup.setLayout(volFilterLayout)
        self.volFilterGroup.toggled.connect(self.volumeFilterToggled_cb)
        self.volFilterGroup.adjustSize()
        minLab = QtGui.QLabel("min")
        self.volMinB = QtGui.QDoubleSpinBox()
        maxLab = QtGui.QLabel("max")
        self.volMaxB = QtGui.QDoubleSpinBox()
        volFilterLayout.addWidget(minLab, 0,0)
        volFilterLayout.addWidget(self.volMinB, 0, 1)
        volFilterLayout.addWidget(maxLab, 1, 0)
        volFilterLayout.addWidget(self.volMaxB, 1, 1)
        self.layout.addWidget(self.volFilterGroup)
        self.volFilterGroup.adjustSize()
        self.layout.setContentsMargins(0,0,0,0)
        table.itemChanged.connect(self.itemChanged_cb)
        self.setLayout(self.layout)
        self.fillTable() # this method computes self.maxVolume and self.minVolume
        # and sets the range and the value of the spin boxes
        self.volMinB.valueChanged.connect(self.setMinVolume_cb)
        self.volMaxB.valueChanged.connect(self.setMaxVolume_cb)
        self.saveInitialValues()

    def saveInitialValues(self):
        self.initValues = {'cavities':self.cavitiesCB.isChecked(),
                           'volume filter':self.volFilterGroup.isChecked(),
                           'minvolume':self.volMinB.value(),
                           'maxvolume':self.volMaxB.value(),
                           }
        
    def fillTable(self):
        mol = self.molSel.getAtomGroup().getMolecule()
        surf = mol._msmsData['msms'][self.surfName][mol._ag.getACSIndex()]
        s = surf.sesr.fst
        volumes = self.volumes = []
        names = self.names = []
        areas = self.areas = []
        while s:
            volume = s.n_ses_volume
            if volume > 0:
                names.append("outer surf")
            else:
                names.append("cavity")
            volumes.append(round(abs(volume), 2))
            areas.append(s.n_ses_area)
            s = s.nxt
        # sort list volumes in descending order
        self.sortedVolIndices = inds = numpy.argsort(volumes)[::-1]
        maxVolume = volumes[inds[0]]
        minVolume = volumes[inds[-1]]
        self.immediate = False
        if minVolume == maxVolume:
            self.volFilterGroup.setChecked(False)
        self.fillTable_()
        self.volMaxB.setRange(0, maxVolume+100)
        self.volMinB.setRange(0, maxVolume+100)
        self.volMaxB.setValue(self.maxVolume)
        self.volMinB.setValue(self.minVolume)
        self.immediate = True
        #print "maxval:" , self.maxVolume, "minval:", self.minVolume
        #if not self.table.isSortingEnabled():
        #    self.table.setSortingEnabled(True)

    def fillTable_(self):
        self.immediate = False
        count = 0
        rowcount = 0
        self.tableVisibleCheckItems = []
        self.tableFlipCheckItems = []
        mol = self.molSel.getAtomGroup().getMolecule()
        prop = mol._renderingProp['msms'][self.surfName]
        minval = round(self.volMinB.value(), 2)
        maxval = round(self.volMaxB.value(), 2)
        minVolume = maxVolume = 0
        for count, ind in enumerate(self.sortedVolIndices):
            volume = self.volumes[ind]
            text = "%s %d"%(self.names[ind], count+1)
            if self.volFilterGroup.isChecked() and len(self.volumes) > 1:
                if volume > maxval: continue
                elif volume < minval: break
            nrows = self.table.rowCount()
            #print "addVolumeRow", "nrows", nrows,  "rowcount", rowcount
            if rowcount > nrows-1:
                #print "inserting row", rowcount
                self.table.insertRow(rowcount)
            #column name
            nameItem = QtGui.QTableWidgetItem(text)
            self.table.setItem(rowcount, 0, nameItem)
            #column show
            show = prop['compVisible'][ind] 
            checkItem = QtGui.QTableWidgetItem()
            checkItem.setFlags(QtCore.Qt.ItemIsUserCheckable|
                               QtCore.Qt.ItemIsEnabled)
            if show:
                checkItem.setCheckState(QtCore.Qt.Checked)
            else:
                checkItem.setCheckState(QtCore.Qt.Unchecked)                
            self.table.setItem(rowcount, 1, checkItem)
            self.tableVisibleCheckItems.append([checkItem, ind])
            #column  flip
            flip = prop['flip'][ind]
            flipItem = QtGui.QTableWidgetItem()
            flipItem.setFlags(QtCore.Qt.ItemIsUserCheckable|
                              QtCore.Qt.ItemIsEnabled)
            if flip:
                flipItem.setCheckState(QtCore.Qt.Checked)
            else:
                flipItem.setCheckState(QtCore.Qt.Unchecked)
            self.table.setItem(rowcount, 2, flipItem)
            self.tableFlipCheckItems.append([flipItem, ind])
            #column vol            
            volumeItem = QtGui.QTableWidgetItem("%.2f"%volume)
            self.table.setItem(rowcount, 3, volumeItem)
            #column area
            areaItem = QtGui.QTableWidgetItem("%.2f" % self.areas[ind])
            self.table.setItem(rowcount, 4, areaItem)
            if  rowcount == 0:
                minVolume = maxVolume = volume
            else:
                if volume < minVolume: minVolume = volume
                elif volume > maxVolume: maxVolume = volume
            rowcount += 1
        self.minVolume = minVolume
        self.maxVolume = maxVolume
        self.immediate = True

    def itemChanged_cb(self, item):
        if not self.immediate: return
        column = item.column()
        if column  == 1 or column  == 2:
            checked = item.checkState() == QtCore.Qt.Checked
            name = self.table.item(item.row(), 0).text()
            nitem = int(name.split(" ")[-1])-1
            volind = self.sortedVolIndices[nitem]
            #print "Item changed", n, item.row(), item.column()
            mol = self.molSel.getAtomGroup().getMolecule()
            if column == 1: #"show"
                #if not self.initValues.has_key('visible'):
                #    self.initValues['visible']=mol._renderingProp['msms'][self.surfName]['compVisible'][:]
                mol._renderingProp['msms'][self.surfName]['compVisible'][volind] = checked
                self.gui.app.displayMSMS.refreshDisplay(mol, names=[self.surfName])
            else: #"flip"
                #if not self.initValues.has_key('flipped'):
                #    self.initValues['flipped']=mol._renderingProp['msms'][self.surfName]['flip'][:]
                mol._renderingProp['msms'][self.surfName]['flip'][volind] = checked
                self.gui.app.displayMSMS.refreshDisplay(mol, names=[self.surfName])
            
    def allVisible_cb(self):
        self.immediate = False
        mol = self.molSel.getAtomGroup().getMolecule()
        #if not self.initValues.has_key('visible'):
        #    self.initValues['visible']=mol._renderingProp['msms'][self.surfName]['compVisible'][:]
        for checkB, ind in self.tableVisibleCheckItems:
            checkB.setCheckState(QtCore.Qt.Checked)
            mol._renderingProp['msms'][self.surfName]['compVisible'][ind] = True
        self.immediate = True
        self.gui.app.displayMSMS.refreshDisplay(mol, names=[self.surfName])

    def invert_cb(self):
        self.immediate = False
        mol = self.molSel.getAtomGroup().getMolecule()
        #if not self.initValues.has_key('visible'):
        #    self.initValues['visible']=mol._renderingProp['msms'][self.surfName]['compVisible'][:]
        for checkB, ind in self.tableVisibleCheckItems:
            if checkB.checkState() == QtCore.Qt.Checked:
                checkB.setCheckState(QtCore.Qt.Unchecked)
                checked = False
            else: 
                checkB.setCheckState(QtCore.Qt.Checked)
                checked = True
            mol._renderingProp['msms'][self.surfName]['compVisible'][ind] = checked
        self.immediate = True
        self.gui.app.displayMSMS.refreshDisplay(mol, names=[self.surfName])

    def cavities_cb(self, val):
        #print "cavities_cb:", val, "selection", self.molSel
        mol = self.molSel.getAtomGroup().getMolecule()
        self.gui.app.computeMSMS(mol.select(), surfName=self.surfName, cavities=val)
        self.gui.app.displayMSMS(self.molSel, surfNames=[self.surfName])
        self.table.clearContents()
        if self.table.rowCount()> 1:
            self.table.setRowCount(1)
        self.fillTable()

    def setMinVolume_cb(self, val):
        if not self.immediate: return
        #print "setMinVolume_cb", val,
        if self.minVolume != val:
            self.setMinMaxVolume(val, self.maxVolume)
            self.minVolume = val
        #print "minVolume = ",  self.minVolume

    def setMaxVolume_cb(self, val):
        if not self.immediate: return
        #print "setMaxVolume_cb", val,
        if self.maxVolume != val:
            self.setMinMaxVolume(self.minVolume, val)
            self.maxVolume = val
        #print "maxVolume = ",  self.maxVolume

    def setMinMaxVolume(self, minval, maxval):
        mol = self.molSel.getAtomGroup().getMolecule()
        prop = mol._renderingProp['msms'][self.surfName]
        #import pdb; pdb.set_trace()
        for ind in self.sortedVolIndices:
            volume = self.volumes[ind]
            if volume >= minval and volume <= maxval:
                prop['compVisible'][ind] = True
            else:
                prop['compVisible'][ind] = False
        self.table.clearContents()
        if self.table.rowCount()> 1:
            self.table.setRowCount(1)
        self.fillTable_()
        self.gui.app.displayMSMS.refreshDisplay(mol, names=[self.surfName])


    def volumeFilterToggled_cb(self, val):
        if not self.immediate: return
        update = False
        oldmaxVol = maxVol = self.volMaxB.value()
        oldminVol = minVol = self.volMinB.value()
        if not val: # unchecked
            maxVol = self.volumes[self.sortedVolIndices[0]]
            minVol = self.volumes[self.sortedVolIndices[-1]]
        #if maxVol != self.maxVolume:
        if abs(maxVol-self.maxVolume) > 0.01:
            update = True
        if abs(minVol-self.minVolume) > 0.01:
            update = True
        if update:
            mol = self.molSel.getAtomGroup().getMolecule()
            prop = mol._renderingProp['msms'][self.surfName]
            for ind in self.sortedVolIndices:
                volume = self.volumes[ind]
                if volume > oldmaxVol or volume < oldminVol:
                    if val: prop['compVisible'][ind] = False
                    else: prop['compVisible'][ind] = True
            self.table.clearContents()
            if self.table.rowCount()> 1:
                self.table.setRowCount(1)
            self.fillTable_()
            self.gui.app.displayMSMS.refreshDisplay(mol, names=[self.surfName])

    def setNewSurface(self, selection, surfName):
        #print "setNewSurface:", len(self.molSel), len(selection), surfName
        self.surfName = surfName
        self.molSel = selection
        self.immediate = False
        self.volFilterGroup.setChecked(False)
        self.updateTable()
        self.cavitiesCB.setChecked(len(self.volumes) > 1)
        self.immediate = True
        self.saveInitialValues()

    def updateTable(self):
        self.table.clearContents()
        if self.table.rowCount()> 1:
            self.table.setRowCount(1)
        self.fillTable()

    def cancel(self):
        # cancel all the changes
        # Cavities check button
        cavities = self.initValues['cavities']
        if cavities != self.cavitiesCB.isChecked():
            self.cavitiesCB.setChecked(cavities)
        # volume filter group
        self.immediate = True
        if cavities:
            volfilter = self.initValues['volume filter']
            if volfilter != self.volFilterGroup.isChecked():
                self.immediate = False
                # set initial max and min values
                self.volMinB.setValue(self.initValues['minvolume'])
                self.volMaxB.setValue(self.initValues['maxvolume'])
                self.immediate = True
                self.volFilterGroup.setChecked(volfilter)
            else:
                if volfilter == True:
                    maxVol = self.volMaxB.value()
                    minVol = self.volMinB.value()
                    if abs(maxVol-self.initValues['maxvolume']) > 0.01:
                        self.volMaxB.setValue(self.initValues['maxvolume'])
                    if abs(minVol-self.initValues['minvolume']) > 0.01:
                        self.volMinB.setValue(self.initValues['minvolume'])
            #display = False
            # if the checkbuttons in the table were used self.initValues will have
            # the 'visible' or|and 'flipped' keys. 
            #if self.initValues.has_key('visible'):
            #    mol = self.molSel.getAtomGroup().getMolecule()
            #    visible = self.initValues.pop('visible')
            #    mol._renderingProp['msms'][self.surfName]['compVisible'] = visible
            #    self.immediate = False
            #    for checkB, ind in self.tableVisibleCheckItems:
            #        checkB.setCheckState(QtCore.Qt.Checked) if visible[ind] else checkB.setCheckState(QtCore.Qt.Unchecked)
            #    self.immediate = True
            #    display = True
            #if self.initValues.has_key('flipped'):
            #    mol = self.molSel.getAtomGroup().getMolecule()
            #    flipped = self.initValues.pop('flipped')
            #    mol._renderingProp['msms'][self.surfName]['flip'] = flipped
            #    self.immediate = False
            #    for checkB, ind in self.tableFlipCheckItems:
            #        checkB.setCheckState(QtCore.Qt.Checked) if flipped[ind] else checkB.setCheckState(QtCore.Qt.Unchecked)
            #    self.immediate = True
            #    display = True
            #    if display:
            #        self.gui.app.displayMSMS.refreshDisplay(mol, names=[self.surfName])

            
class MsmsCustomizeWidget(IndexedGeomCustomizeWidget):
    
    def __init__(self,name, molSel, selName, perAtomTabs=["Colors"], perMoleculeTabs=["Surface parameters", "Surfaces"], parent=None):
        super(MsmsCustomizeWidget, self).__init__(name, molSel, selName, perAtomTabs=perAtomTabs, perMoleculeTabs=perMoleculeTabs, parent=parent)
        self.app = self.selection[0].getAtomGroup().getMolecule().app()
        self.displayCmd = self.app.displayMSMS
        self.geomName = "msms"
        self.surfaces = {}
        self.modelGlobalWidgets = None
        self.colorWidgets = None
        self.addColorsWidgets()
        self.surfaceTableWidget = None
        self.backFaceColor = []
        self.saveInitValues()
        self.immediate = True

    def updateSelection(self, selection, selName):
        self.selection = selection
        self.tabWidget.setTabText(0, selName)
        if self.surfaceTableWidget:
            names = []
            for sel in self.selection:
                names.append([sel, self.getSurfName(sel)])
                if len(names) > 1:
                    self.surfPageWidget.setCurrentIndex(0)
                else:
                    if self.surfaceTableWidget == "Not available":
                        self.surfaceTableWidget = SurfaceTableWidget(self, names[0][0], names[0][1][0])
                        self.surfPageWidget.addWidget(self.surfaceTableWidget)
                    else:
                        self.surfaceTableWidget.setNewSurface(names[0][0], names[0][1][0])
                    self.surfPageWidget.setCurrentIndex(1)

    def getSurfName(self, sel):
        mol = sel.getAtomGroup().getMolecule()
        gc = mol.geomContainer
        names = []
        for name in gc.atoms.keys():
            if name.find("msms") >= 0:
                if gc.geoms.has_key(name): 
                    names.append(gc.geoms[name]._pmvName)
        return names
    
    def saveInitValues(self):
        self.initialValues = {}
        if self.modelGlobalWidgets:
            self.modelGlobalWidgets.saveWidgetsValues()
        if self.opacityWidget:
            self.opacityWidget.saveWidgetsValues()

    def addPropPerMoleculeWidgets(self):
        from DejaVu2.Qt.geomParamGUI import MsmsSurfaceParams
        if not self.modelGlobalWidgets:
            self.modelGlobalWidgets = glw = MsmsSurfaceParams(self)
            layout = QtGui.QVBoxLayout()
            layout.addWidget(glw)
            self.globalPropTabWidget.widget(0).setLayout(layout)
            glw.qualityChanged.connect(self.setQuality)
            glw.probeRadiusChanged.connect(self.setProbeRadius)
            glw.computeSES.connect(self.computeSES)
            glw.computeSAS.connect(self.computeSAS)
            glw.softColorBoundariesToggled.connect(self.setSoftColorBoundaries)
            glw.strokeWidget.toggleStroke.connect(self.setStroke)
            glw.strokeWidget.backFaceColorModeChanged.connect(self.setBackFaceMode)
            glw.strokeWidget.backFaceColorChanged.connect(self.setBackFaceColor)
            glw.strokeWidget.lineWidthChanged.connect(self.setLineWidth)
            layout.addStretch(1)
            self.modelGlobalWidgets.saveWidgetsValues()
            
    def globPropTabChanged_cb(self, index):
        #applyB = self.parent.buttons.button(self.parent.buttons.Apply)
        if self.globalPropTabWidget.tabText(index) == "Surfaces":
            if not self.surfaceTableWidget:
                self.addSurfaceTableWidget()
            else:
                self.surfaceTableWidget.updateTable()
            #if applyB.isEnabled():
                #applyB.setEnabled(False)
                #self.parent.immediateCheckBox.setEnabled(False)
        #else:
            #self.parent.immediateCheckBox.setEnabled(True)
            #val = self.parent.immediateCheckBox.isChecked()
            #applyB.setEnabled(not val)

    def setProbeRadius(self, probeRad):
        for sel in self.selection:
            names = self.getSurfName(sel)
            mol =  sel.getAtomGroup().getMolecule()
            if not self.initialValues.has_key(sel):
                self.initialValues[sel]={}
            for name in names:
                if not self.initialValues[sel].has_key(name):
                    self.initialValues[sel][name] = {}
                if not self.initialValues[sel][name].has_key('probeRadius'):
                    self.initialValues[sel][name]['probeRadius'] = mol._msmsData['params'][name]['probeRadius']
                atoms = mol._msmsData['atoms'][name]
                self.app.computeMSMS(atoms, surfName=name, probeRadius=probeRad)
            self.app.displayMSMS(sel, surfNames=names)

    def setQuality(self, dens):
        for sel in self.selection:
            if not self.initialValues.has_key(sel):
                self.initialValues[sel]={}
            names = self.getSurfName(sel)
            mol =  sel.getAtomGroup().getMolecule()
            for name in names:
                if not self.initialValues[sel].has_key(name):
                    self.initialValues[sel][name] = {}
                if not self.initialValues[sel][name].has_key('quality'):
                    self.initialValues[sel][name]['quality'] = mol._msmsData['params'][name]['density']    
                atoms = mol._msmsData['atoms'][name]
                self.app.computeMSMS(atoms, surfName=name, density=dens)
            self.app.displayMSMS(sel, surfNames=names)

    def updateGlobalProp(self, kw):
        for sel in self.selection:
            names = self.getSurfName(sel)
            mol =  sel.getAtomGroup().getMolecule()
            if not self.initialValues.has_key(sel):
                self.initialValues[sel]={}
            for name in names:
                if not self.initialValues[sel].has_key(name):
                    self.initialValues[sel][name] = {}
                prop = mol._renderingProp['msms'][name]
                for propName, value in kw.items():
                    if not self.initialValues[sel][name].has_key(propName):
                        try:
                            if isinstance(prop[propName], str) and self.displayCmd._renderToPmv.has_key(prop[propName]):
                                self.initialValues[sel][name][propName]= self.displayCmd._renderToPmv[prop[propName]]
                            else:
                                self.initialValues[sel][name][propName] = prop[propName]
                        except:
                            print "No rendering properties %s for geometry %s" %(propName, geomName)
                self.app.displayMSMS.updateModelGlobals(mol, name, prop, **kw)
            self.app.displayMSMS.refreshDisplay(mol, names)

    def setPropInitValues(self, sel, name):
        if not self.initialValues[sel].has_key(name):
            self.initialValues[sel][name] = {}
        if not self.initialValues[sel][name].has_key('radii'):
            mol = sel.getAtomGroup().getMolecule()
            radii = mol._msmsData['params'][name]['radii']
            if radii is None:
                radii = sel.getRadii().copy()
            self.initialValues[sel][name]['radii'] = radii
            if not self.initialValues[sel][name].has_key('probeRadius'):
                self.initialValues[sel][name]['probeRadius'] = mol._msmsData['params'][name]['probeRadius']
            if not self.initialValues[sel][name].has_key('quality'):
                self.initialValues[sel][name]['quality'] = mol._msmsData['params'][name]['density']

    def computeSES(self, probeRad, qual):
        #import pdb; pdb.set_trace()
        for sel in self.selection:
            if not self.initialValues.has_key(sel):
                self.initialValues[sel]={}
            mol =  sel.getAtomGroup().getMolecule()
            names = self.getSurfName(sel)
            #print "computeSES:", sel, names
            for name in names:
                self.setPropInitValues(sel, name)
                atoms = mol._msmsData['atoms'][name]
                radii = mol._ag.getRadii()
                self.app.computeMSMS(atoms, surfName=name, radii=radii, probeRadius=probeRad, density=qual)
            self.app.displayMSMS(sel, surfNames=names)

    def computeSAS(self, qual):
        #import pdb; pdb.set_trace()
        for sel in self.selection:
            if not self.initialValues.has_key(sel):
                self.initialValues[sel]={}
            mol =  sel.getAtomGroup().getMolecule()
            names = self.getSurfName(sel)
            for name in names:
                self.setPropInitValues(sel, name)
                if not self.initialValues[sel].has_key(name):
                   self.initialValues[sel][name] = {}
                if not self.initialValues[sel][name].has_key('radii'):
                   self.initialValues[sel][name]['radii'] = mol._msmsData['params'][name]['radii']
                radii = mol._ag.getRadii() + 1.5
                atoms = mol._msmsData['atoms'][name]
                self.app.computeMSMS(atoms, radii=radii, probeRadius=0.1, density=qual)
            self.app.displayMSMS(sel, surfNames=names)

    def setStroke(self, val, strokeWidth):
        # callback of the Stroke widget (sets Stroke on/off)
        if val:
            kw = {'culling':'none', 'frontRendering':'solid',
                  'backRendering':'mesh', 'linewidth':strokeWidth}
            #Check if the initial backFaceColor has been saved. If it has not,
            # then the color has not been modified --> we set the initial stroke color
            # to black.
            sel = self.selection[0]
            setColor = True
            if self.initialValues.has_key(sel):
                name = self.getSurfName(sel)[0]
                if self.initialValues[sel].has_key(name):
                    if self.initialValues[sel][name].has_key("backfaceColor"):
                        setColor = False
            if setColor:
                # initial stroke color is black
                kw["backfaceColor"]=[0.0, 0.0, 0.0, 1.0]
        else:
            kw = {'culling':'back', 'frontRendering':'solid',
                  'backRendering':'solid', 'linewidth':strokeWidth}
        self.updateGlobalProp(kw)

    def setBackFaceColor(self, r, g, b, a, saveCurrentColor=False, restoreColor=False):
        #print "setBackFaceColor:", "saveCurrentColor:", saveCurrentColor, "restoreColor", restoreColor
        if restoreColor:
            if len(self.backFaceColor):
                #print "restoring backface color"
                for cl in  self.backFaceColor:
                    sel = cl[0]
                    mol =  sel.getAtomGroup().getMolecule()
                    names = []
                    for name, color in cl[1:]:
                        prop = mol._renderingProp['msms'][name]
                        self.displayCmd.updateModelGlobals(mol, name, prop,
                                                           backfaceColor=color)
                        names.append(name)
                    self.displayCmd.refreshDisplay(mol, names)
            self.backFaceColor=[]
            return
        if saveCurrentColor and len(self.backFaceColor) == 0:
            #print "saving back color"
            for sel in self.selection:
                mol =  sel.getAtomGroup().getMolecule()
                names = self.getSurfName(sel)
                cl = [sel]
                for name in names:
                    cl.append([name, mol._renderingProp['msms'][name]['backfaceColor']])
                self.backFaceColor.append(cl)
        elif saveCurrentColor==False and restoreColor==False:
            self.backFaceColor = []
        kw = {"backfaceColor":[r,g,b,a]}
        #print "updateGlobalProp: backfaceColor"
        self.updateGlobalProp(kw)

    def addColorsWidgets(self):
        if self.colorWidgets: return
        names = ["msms"]
        self.colorWidgets = ColorsWidget(self, geomNames=names)
        self.colorWidgets.beforeColorChanged.connect(self.beforeGeomColorChanged)
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.colorWidgets)
        self.opacityWidget = OpacityWidget()
        layout.addWidget(self.opacityWidget)
        self.opacityWidget.opacityChanged.connect(self.setOpacity)
        layout.addStretch(1)
        self.perAtomPropTabWidget.widget(0).setLayout(layout)

    def saveCurrentGeomColor(self):
        geomcol = []
        for sel in self.selection:
            mol =  sel.getAtomGroup().getMolecule()
            names = self.getSurfName(sel)
            item = [sel]
            if not self.initialValues.has_key(sel):
                self.initialValues[sel] = {}
            for name in names:
                col = mol._colors['msms'][name].copy()
                inds = mol._ag.getData('colorsIndices_msms_%s'%name)
                item.append([name, col, inds])
                if not self.initialValues[sel].has_key(name):
                    self.initialValues[sel][name] = {}
                if not self.initialValues[sel][name].has_key('colors'):
                    self.initialValues[sel][name]['colors'] = [col, inds.copy()]
            geomcol.append(item)
        return geomcol

    def restoreGeomColor(self, geomcol):
        for item in geomcol:
            sel = item[0]
            names = []
            mol =  sel.getAtomGroup().getMolecule()
            for name, color, inds in item[1:]:
                mol._ag.setData('colorsIndices_msms_%s'%name, inds)
                mol._colors['msms'][name] = color
                names.append(name)
            mol.app().displayMSMS.refreshDisplay(mol, names)
            
    def setGeomColor(self, rgb, carbonsOnly=False):
        #print "setGeomColor:", rgb
        for sel in self.selection:
            if carbonsOnly:
                sel = sel & sel.select("element C")
            inds = sel.getIndices()
            mol = sel.getAtomGroup().getMolecule()
            names = self.getSurfName(sel)
            for name in names:
                col = mol._colors['msms'][name]
                mol._ag._data['colorsIndices_msms_%s'%name][inds] = len(col)
                mol._colors['msms'][name] = numpy.concatenate((col, numpy.array([rgb])))
            mol.app().displayMSMS.refreshDisplay(mol, names)
    
    def setOpacity(self, val):
        for sel in self.selection:
            if not self.initialValues.has_key(sel):
                self.initialValues[sel] = {}
            mol =  sel.getAtomGroup().getMolecule()
            inds = sel.getIndices()
            names = self.getSurfName(sel)
            for name in names:
                if not self.initialValues[sel].has_key(name):
                    self.initialValues[sel][name] = {}
                dataLabel = 'opacity_msms_%s'%name
                if mol._ag.getData(dataLabel) is None:
                    opacity = numpy.ones(len(mol._ag), 'f')
                    if not self.initialValues[sel][name].has_key("opacity"):
                        self.initialValues[sel][name]["opacity"]=opacity.copy()
                    opacity[inds] = val
                    mol._ag.setData(dataLabel, opacity)
                else:
                    if not self.initialValues[sel][name].has_key("opacity"):
                        self.initialValues[sel][name]["opacity"] = mol._ag._data[dataLabel].copy()
                    mol._ag._data[dataLabel][inds] = val
            self.displayCmd.refreshDisplay(mol, names)

    def addSurfaceTableWidget(self):
        if self.surfaceTableWidget: return
        tab = self.globalPropTabWidget.widget(1)
        layout = QtGui.QVBoxLayout()
        self.surfPageWidget = QtGui.QStackedWidget()
        lab = QtGui.QLabel("This panel is not available for multiple molecules.\nTo use this panel select a single molecule.")
        self.surfPageWidget.addWidget(lab)
        names = []
        for sel in self.selection:
            names.append([sel, self.getSurfName(sel)])
        if len(names) > 1:
            self.surfaceTableWidget = "Not available"
            self.surfPageWidget.setCurrentIndex(0)
        else:
            sel = names[0][0]
            surfName = names[0][1][0]
            mol = sel.getAtomGroup().getMolecule()
            self.surfaceTableWidget = SurfaceTableWidget(self, sel, surfName)
            self.surfPageWidget.addWidget(self.surfaceTableWidget)
            self.surfPageWidget.setCurrentIndex(1)
        layout.addWidget(self.surfPageWidget)
        tab.setLayout(layout)
        
    def cancel(self):
        # restore initial values of the geometry's parameters and the widgets
        if self.modelGlobalWidgets:
            self.modelGlobalWidgets.restoreWidgetsValues()
        if self.opacityWidget:
            self.opacityWidget.restoreWidgetsValues()
        for sel in self.selection:
            if not self.initialValues.get(sel): continue
            mol =  sel.getAtomGroup().getMolecule()
            names = self.getSurfName(sel)
            for name in names:
                kw = self.initialValues[sel].get(name, {})
                if len(kw):
                    prop = mol._renderingProp['msms'][name]
                    probeRadius = None
                    dens = None
                    radii = None
                    recompute = False
                    if kw.has_key('probeRadius'):
                       probeRadius = kw.pop('probeRadius')
                       recompute = True
                    if kw.has_key('quality'):
                        dens = kw.pop('quality')
                        recompute = True
                    if kw.has_key('radii'):
                        radii = kw.pop('radii')
                        recompute = True
                    if recompute:
                        atoms = mol._msmsData['atoms'][name]
                        self.app.computeMSMS(atoms, surfName=name, 
                                             probeRadius=probeRadius,
                                             radii=radii, density=dens)
                    if kw.has_key('opacity'):
                        opacity = kw.pop("opacity")
                        dataLabel = 'opacity_msms_%s'%name
                        mol._ag.setData(dataLabel, opacity)
                    if kw.has_key("colors"):
                        col, inds = kw.pop("colors")
                        mol._ag.setData('colorsIndices_msms_%s'%name, inds)
                        mol._colors['msms'][name] = col
                    if len(kw):
                        self.app.displayMSMS.updateModelGlobals(mol, name, prop,
**kw)
                    self.app.displayMSMS(sel, surfNames=[name])
        if self.surfaceTableWidget:
            self.surfaceTableWidget.cancel()

    def getPropValues(self):
        pass

    
from DejaVu2.Shapes import Rectangle2D, Circle2D, Ellipse2D
class  CartoonCustomizeWidget(IndexedGeomCustomizeWidget):

    def __init__(self,name, molSel, selName, perAtomTabs=["Colors"], parent=None):
        super(CartoonCustomizeWidget, self).__init__(name, molSel, selName, perAtomTabs=perAtomTabs, parent=parent)
        self.app = self.selection[0].getAtomGroup().getMolecule().app()
        self.displayCmd = self.app.displayCartoon
        self.geomName="cartoon"
        self.modelGlobalWidgets = None
        self.colorWidgets = None
        self.addColorsWidgets()
        self.backFaceColor = []
        self.saveInitValues()
        self.immediate = True

    def saveInitValues(self):
        #print "saveInitValues", self
        self.initialValues = {}
        self.shape = "default"
        self.quality = 10
        self.scale = 1.0
        if self.modelGlobalWidgets:
            self.modelGlobalWidgets.saveWidgetsValues()
            #default = self.modelGlobalWidgets.initialFormValues['defaultShape']
            #if default: self.shape = "default"
            #else: self.shape = "cylinder"
        if self.opacityWidget:
            self.opacityWidget.saveWidgetsValues()

    def addPropPerMoleculeWidgets(self):
        from DejaVu2.Qt.geomParamGUI import CartoonModelGlobalParams
        if not self.modelGlobalWidgets:
            self.modelGlobalWidgets = glw = CartoonModelGlobalParams(self)
            layout = QtGui.QVBoxLayout()
            layout.addWidget(glw)
            self.tabWidget.widget(1).setLayout(layout)
            glw.qualityChanged.connect(self.setQuality)
            glw.shapeChanged.connect(self.setShape)
            glw.scaleChanged.connect(self.setScale)
            glw.strokeWidget.toggleStroke.connect(self.setStroke)
            glw.strokeWidget.backFaceColorModeChanged.connect(self.setBackFaceMode)
            glw.strokeWidget.backFaceColorChanged.connect(self.setBackFaceColor)
            glw.strokeWidget.lineWidthChanged.connect(self.setLineWidth)
            layout.addStretch(1)
            self.modelGlobalWidgets.saveWidgetsValues()

    def setPropInitValues(self, sel):
        if not self.initialValues.has_key(sel):
             self.initialValues[sel] = {}
        if not self.initialValues[sel].has_key('shape'):
            mol = sel.getAtomGroup().getMolecule()
            kw = mol._renderingProp.get('cartoon', None)
            if kw is None: return
            shapeDict = {}
            shapeDict['helixShape'] = kw.get('helixShape', None)
            shapeDict['sheetShape'] = kw.get('sheetShape', None)
            shapeDict['coilShape'] = kw.get('coilShape', None)
            shapeDict['shapeNABB'] = kw.get('shapeNABB', None)
            self.initialValues[sel]['shape'] = shapeDict
            
    def setQuality(self, quality):
        scale = self.scale
        for sel in self.selection:
            self.setPropInitValues(sel)
            if self.shape == "default":
                self.app.computeCartoon(sel, quality=quality, scale=scale)
            else:
               helixShape = Circle2D(radius=0.2*scale, quality=quality)
               sheetShape = Circle2D(radius=0.2*scale, quality=quality)
               coilShape = Circle2D(radius=0.2*scale, quality=quality)
               shapeNABB = Circle2D(radius=0.2*scale, quality=quality)
               self.app.computeCartoon(sel, helixShape=helixShape, sheetShape=sheetShape, coilShape=coilShape, shapeNABB=shapeNABB, quality=quality)
            self.displayCmd(sel)
        self.quality = quality

    def setShape(self, shape="default"):
        assert shape in ["default", "cylinder"]
        quality = self.quality
        scale = self.scale
        for sel in self.selection:
            self.setPropInitValues(sel)
            if shape == "default":
                self.app.computeCartoon(sel, quality=quality, scale=scale)
            else:
               helixShape = Circle2D(radius=0.2*scale, quality=quality)
               sheetShape = Circle2D(radius=0.2*scale, quality=quality)
               coilShape = Circle2D(radius=0.2*scale, quality=quality)
               shapeNABB = Circle2D(radius=0.2*scale, quality=quality)
               self.app.computeCartoon(sel, helixShape=helixShape, sheetShape=sheetShape, coilShape=coilShape, shapeNABB=shapeNABB, quality=quality)
            self.displayCmd(sel)
        self.shape = shape

    def setScale(self, scale):
        quality = self.quality
        shape = self.shape
        for sel in self.selection:
            self.setPropInitValues(sel)
            if shape == "default":
                self.app.computeCartoon(sel, quality=quality, scale=scale)
            else:
               helixShape = Circle2D(radius=0.2*scale, quality=quality)
               sheetShape = Circle2D(radius=0.2*scale, quality=quality)
               coilShape = Circle2D(radius=0.2*scale, quality=quality)
               shapeNABB = Circle2D(radius=0.2*scale, quality=quality)
               self.app.computeCartoon(sel, helixShape=helixShape, sheetShape=sheetShape, coilShape=coilShape, shapeNABB=shapeNABB, quality=quality)
            self.displayCmd(sel)
        self.scale = scale

    def addColorsWidgets(self):
        if self.colorWidgets: return
        names = ["cartoon"]
        self.colorWidgets = ColorsWidget(self, geomNames=names)
        self.colorWidgets.beforeColorChanged.connect(self.beforeGeomColorChanged)
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.colorWidgets)
        #self.opacityWidget = OpacityWidget()
        #layout.addWidget(self.opacityWidget)
        #self.opacityWidget.opacityChanged.connect(self.setOpacity)
        layout.addStretch(1)
        self.perAtomPropTabWidget.widget(0).setLayout(layout)

    def cancel(self):
        # restore initial values of the geometry's parameters and the widgets
        if self.modelGlobalWidgets:
            self.modelGlobalWidgets.restoreWidgetsValues()
        #import pdb; pdb.set_trace()
        for sel in self.selection:
            mol = sel.getAtomGroup().getMolecule()
            initVals = self.initialValues.get(sel, {})
            if not len(initVals): continue
            if initVals.has_key("shape"):
                shapeDict = initVals.get('shape', None)
                rendProp = mol._renderingProp['cartoon']
                if shapeDict is not None:
                    helixShape = shapeDict['helixShape']
                    sheetShape = shapeDict['sheetShape']
                    coilShape = shapeDict['coilShape']
                    shapeNABB = shapeDict['shapeNABB']
                else:
                    helixShape = rendProp.get('helixShape')
                    sheetShape = rendProp.get('sheetShape')
                    coilShape  = rendProp.get('coilShape')
                    shapeNABB = rendProp.get('shapeNABB')
                self.app.computeCartoon(sel, helixShape=helixShape, sheetShape=sheetShape, coilShape=coilShape, shapeNABB=shapeNABB)
                initVals.pop('shape')
            if initVals.has_key("colors"):
                colors, inds = initVals.pop("colors")
                mol._ag.setData('colorsIndices_cartoon', inds)
                mol._colors['cartoon'] = colors
            if len(initVals):
                self.displayCmd.updateModelGlobals(mol, **initVals)
            self.displayCmd.refreshDisplay(mol)

class AtomLabelsCustomoizeWidget(CustomizeWidget):
    
    def __init__(self,name, molSel, selName, perAtomTabs=["Colors"], parent=None):
        super(AtomLabelsCustomoizeWidget, self).__init__(name, molSel, selName, parent=parent, perAtomTabs=perAtomTabs)
        self.displayCmd = self.selection[0].getAtomGroup().getMolecule().app().labelAtoms
        self.addPropPerAtomWidgets()
        self.modelGlobalWidgets = None
        self.colorWidgets = None
        self.addColorsWidgets()
        self.saveInitValues()
        #print "label geom name:", self.geomName

    def addPropPerAtomWidgets(self):
        pass

    def addColorsWidgets(self):
        if self.colorWidgets: return
        names = [self.geomName]
        self.colorWidgets = ColorsWidget(self, geomNames=names)
        self.colorWidgets.beforeColorChanged.connect(self.beforeGeomColorChanged)
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.colorWidgets)
        #self.opacityWidget = OpacityWidget()
        #layout.addWidget(self.opacityWidget)
        #self.opacityWidget.opacityChanged.connect(self.setOpacity)
        layout.addStretch(1)
        self.perAtomPropTabWidget.widget(0).setLayout(layout)

    def cancel(self):
        for sel, kw in self.initialValues.items():
            mol =  sel.getAtomGroup().getMolecule()
            if kw.has_key("colors"):
                colors, inds = kw.pop("colors")
                mol._ag.setData('colorsIndices_%s'%self.geomName, inds)
                mol._colors[self.geomName] = colors
                self.displayCmd.refreshDisplay(mol)

class ResidueLabelsCustomoizeWidget(AtomLabelsCustomoizeWidget):
    
    def __init__(self,name, molSel, selName, perAtomTabs=["Colors"], parent=None):
        super(ResidueLabelsCustomoizeWidget, self).__init__(name, molSel, selName, parent=parent, perAtomTabs=perAtomTabs)
        self.displayCmd = self.selection[0].getAtomGroup().getMolecule().app().labelResidues

        
class CustomizeRepresentationsDialog(QtGui.QDialog):
    """
    This class creates a dialog to customize properties of available geometries
    for selected PMV object.
    """
    tabWidgetClass = {"lines": LinesCustomizeWidget,
                      "cpk": CpkCustomizeWidget,
                      "sb": SBCustomizeWidget,
                      "msms": MsmsCustomizeWidget,
                      "cartoon": CartoonCustomizeWidget,
                      #"CoarseMS": CoarseCustomizeWidget,
                      "atomLabels": AtomLabelsCustomoizeWidget,
                      "residueLabels": ResidueLabelsCustomoizeWidget,
                      }

    def __init__(self, selectionSet, selectionName, app, parent=None):
        super(CustomizeRepresentationsDialog, self).__init__(parent)
        self.selSet = selectionSet # pmv selectionSet
        self.selName = selectionName
        self.app = weakref.ref(app)
        self.buildUI(parent)
        
    def buildUI(self, parent):
        # list item on the left, controling apage widget on the right
        # containing ParamTabWidgets corresponding to a given geometry
        lw = self.contentsWidget = QtGui.QListWidget()
        self.contentsWidget.setMaximumWidth(110)
        self.contentsWidget.setMinimumWidth(110)

        #self.immediateCheckBox = QtGui.QCheckBox("Immediate")
        #self.immediateCheckBox.setChecked(True)
        #self.immediateCheckBox.toggled.connect(self.immediate_cb)
        
        pw = self.pagesWidget = QtGui.QStackedWidget()
        self.pageWidgetsDict = {}
        self.createTabWidgets()
        
        #self.paramsLayout = QtGui.QGridLayout()
        #pw.setLayout(self.paramsLayout)

        self.contentsWidget.setCurrentRow(0)
        self.contentsWidget.currentItemChanged.connect(self.contentsWidgetItemChanged_cb)
        buttonBox = QtGui.QGroupBox()
        buttonBoxLayout = QtGui.QGridLayout()
        self.buttons = QtGui.QDialogButtonBox(
            #QtGui.QDialogButtonBox.Ok| QtGui.QDialogButtonBox.Apply |QtGui.QDialogButtonBox.Cancel,
            QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok,
            QtCore.Qt.Horizontal, self)
        #buttonBoxLayout.addWidget(self.immediateCheckBox, 0, 0, 1, 2)
        buttonBoxLayout.addWidget(self.buttons, 0, 2)
        buttonBox.setLayout(buttonBoxLayout)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        #self.buttons.button(self.buttons.Apply).setDefault(False)
        #self.buttons.button(self.buttons.Apply).clicked.connect(self.apply_cb)
        self.buttons.button(self.buttons.Cancel).clicked.connect(self.cancel_cb)
        self.buttons.button(self.buttons.Ok).clicked.connect(self.ok_cb)
        #self.buttons.button(self.buttons.Apply).setEnabled(False)
        self.buttons.button(self.buttons.Ok).setDefault(False)
        
        self.hlayout = QtGui.QHBoxLayout()
        self.setLayout(self.hlayout)
        self.hlayout.addWidget(lw)
        self.vlayout = QtGui.QVBoxLayout()
        self.hlayout.addLayout(self.vlayout)
        self.vlayout.addWidget(pw)
        self.vlayout.addWidget(buttonBox)

        self.saveInitValues()

    def getGeoms(self, obj, visibleOnly=True):
        # return a dictionnary in with the key tells which geometries
        # are available for this selectionSet obj
        geoms = {}
        moleculesAndSelections = [[x.getAtomGroup().getMolecule(),x] for x in obj]
        for mol, sele in moleculesAndSelections:
            gc = mol.geomContainer
            for gname, ats in gc.atoms.items():
                name = None
                if gc.geoms.has_key(gname):
                    if visibleOnly and not gc.geoms[gname].visible:
                        continue
                    if gc.displayedAs([gname], sele, 'fast'):
                        if gname.find("msms") == 0:
                            name = "msms"
                        elif gname.startswith("chain_") and gc.geoms[gname].parent.name == "cartoon":
                            name = "cartoon"
                        else:
                            name = gname
                if name:
                    if not geoms.has_key(name):
                        geoms[name] = []
                    geoms[name].append(sele)
        return geoms 

    def createTabWidgets(self):
        pw = self.pagesWidget
        geoms = self.getGeoms(self.selSet)#self.app().getGeoms(self.selSet)
        #print "createTabWidgets names:", geoms.keys()
        surfName = None
        for name, selection in geoms.items():
            widgetClass = self.tabWidgetClass.get(name, None)
            if widgetClass:
                if not self.pageWidgetsDict.has_key(name):
                    #add a new widget for geometry name
                    widget = widgetClass(name, selection, self.selName, parent=self)
                    self.pageWidgetsDict[name] = widget #append(widget)
                    pw.addWidget(widget)
                    QtGui.QListWidgetItem(name, self.contentsWidget)
                else:
                    # widget for geometry exists, set a new selection
                    # for the widget
                    widget = self.pageWidgetsDict[name]
                    widget.updateSelection(selection, self.selName)
                    items = self.contentsWidget.findItems(name, QtCore.Qt.MatchExactly)
                    if not(len(items)):
                        # name is not in the list box, add it to the list of  available geometries
                        QtGui.QListWidgetItem(name, self.contentsWidget)
        for name in self.pageWidgetsDict.keys():
            if name not in geoms.keys():
                # we have created a widget for geometry name ,
                # but this geometry is not in the current selection.
                # Remove the name from the list chooser.
                txt = str(self.contentsWidget.currentItem().text())
                items = self.contentsWidget.findItems(name, QtCore.Qt.MatchExactly)
                if len(items):
                    row = self.contentsWidget.row(items[0])
                    self.contentsWidget.takeItem(row)
                if name == txt and self.contentsWidget.count()>0:
                    self.contentsWidget.setCurrentRow(0)

    def contentsWidgetItemChanged_cb(self, item):
        #import pdb; pdb.set_trace()
        name = str(item.text())
        w = self.pageWidgetsDict.get(name, None)
        if w:
          self.pagesWidget.addWidget(w)
          self.pagesWidget.setCurrentWidget(w)
        
    def apply_cb(self):
        """Callback of the Apply button"""
        #print 'APPLY'
        self.pagesWidget.currentWidget().getPropValues()

    def cancel_cb(self):
        #self.pagesWidget.currentWidget().cancel()
	for widget in self.pageWidgetsDict.values():
	    widget.cancel()
        QtGui.QDialog.reject(self)

    def ok_cb(self):
        self.pagesWidget.currentWidget().getPropValues()
        QtGui.QDialog.accept(self)

    def setSelection(self, selSet):
        self.selSet = selSet
        for widget in self.pageWidgetsDict.values():
            widget.selection = selSet
        
    def immediate_cb(self, val=None):
        """Callback of the immediate check button. Disables/enables the Apply pushbutton"""
        applyB = self.buttons.button(self.buttons.Apply)
        applyB.setEnabled(not val)
        self.emit(QtCore.SIGNAL("immediateStatusChanged(bool)"), val)

    def updateGUI(self, selSet, selName):
        self.selSet=selSet
        self.selName=selName
        self.createTabWidgets()
        self.saveInitValues()

    def saveInitValues(self):
        for widget in self.pageWidgetsDict.values():
            widget.saveInitValues()

    def keyPressEvent(self, evt):
        if evt.key() == QtCore.Qt.Key_Enter or evt.key() == QtCore.Qt.Key_Return:
            return
        QtGui.QDialog.keyPressEvent(self, evt)

if __name__=="__main__":
    import sys
    from optparse import OptionParser
    app = QtGui.QApplication(sys.argv)
    opt, ok = AdvancedParamGUI.getParams()

    print opt, ok
