#
# $Header: /mnt/raid/services/cvs/ADFR/GUI/trgMapsGUI.py,v 1.1.2.11 2018/01/16 01:07:56 annao Exp $
#
# $Id: trgMapsGUI.py,v 1.1.2.11 2018/01/16 01:07:56 annao Exp $
#

from PySide import QtCore, QtGui
import numpy
import sys, os

from ADFR.utils.maps import MapsFile
from UTpackages.UTisocontour import isocontour
isocontour.setVerboseLevel(0)
from DejaVu2.IndexedPolygons import IndexedPolygons
from mglutil.util.callback import CallbackFunction
from ADFRcc.adfr import GridMap
from DejaVu2.Box import NiceBox
from DejaVu2.Qt.colorpicker import ColorPicker
from DejaVu2.Qt.QColorMapEditor import HLSAColorMapEditor
from DejaVu2.colorMap import ColorMap
from DejaVu2.colorTool import colorsForValueRanges, RGBARamp
from DejaVu2.Textured2DArray import textured2DArray

from mglutil.util.packageFilePath import findFilePath
PMVICONPATH = findFilePath('Icons', 'PmvApp.GUI')
ADFRICONPATH = findFilePath('Icons', 'ADFR.GUI')

from PmvApp.pmvPalettes import AtomElements, elementColors
ElementLookup = {'A': 'C', 'OA': 'O', 'NA': 'N',
                'NX': 'N', 'SA': 'S', 'HD': 'H'}

class TrgMapsGui(QtGui.QDialog):
    closedSignal = QtCore.Signal()
    
    def __init__(self, parent=None, app=None):
        super(TrgMapsGui, self).__init__(parent)
        self.setWindowTitle("Trg Map GUI")
        self.app = None
        if app:
            self.app = app#weakref.ref(app)
        self.viewer = None
        self.trg_data = {}
        self.iso_data = {}
        self.default_maptypes = ['C', 'OA', 'HD']
        self.isoGeomPropWidget = None
        self.orthoGeomPropWidgets = {}
        self.orthoGeomPropWidget = None
        self.numIsoItems = 0
        self.numOrthoItems = 0
        # build GUI

        # Tree widget
        self.vlayout1 = QtGui.QVBoxLayout()
        tw = self.treeWidget = QtGui.QTreeWidget()
        tw.setColumnCount(1)
        tw.setHeaderLabels(['Trg objects',])
        tw.itemExpanded.connect(self.onItemExpanded)
        tw.itemSelectionChanged.connect(self.onItemClick)
        tw.itemDoubleClicked.connect(self.onItemDoubleClick)
        # to build menu on right click on item:
        tw.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        tw.customContextMenuRequested.connect(self.buildTreeMenu)
        #tw.setStyleSheet("QTreeWidget { background: lightGray}")
        self.vlayout1.addWidget(tw)
        
        # geom properties widget
        self.vlayout2 = QtGui.QVBoxLayout()
        pw = self.propWidget = QtGui.QStackedWidget()
        self.stackedWidgets = {}
        label = QtGui.QLabel("Geom Properties")
        self.vlayout2.addWidget(label)
        self.vlayout2.addWidget(pw)
        #self.vlayout2.addStretch()
        self.vlayout2.setContentsMargins(0, 0, 0, 0)
        self.buttonBox = QtGui.QGroupBox()
        buttonBoxLayout = QtGui.QHBoxLayout()
        self.addButton = QtGui.QPushButton()
        #self.addButton.setIcon(QtGui.qApp.style().standardIcon(QtGui.QStyle.SP_FileDialogStart))
        self.addButton.setIcon(QtGui.QIcon(os.path.join(PMVICONPATH, 'open.png')))
        self.addButton.clicked.connect(self.addTrg)
        self.addButton.setToolTip("Open trg file")
        self.removeButton = QtGui.QPushButton()
        self.removeButton.setIcon(QtGui.QIcon(os.path.join(ADFRICONPATH, "cross.png")))
        self.removeButton.clicked.connect(self.deleteTrgItem)
        self.removeButton.setToolTip("Remove trg object")
        #self.visibleButton = QtGui.QPushButton()#"Show/hide")
        #self.visibleButton.setIcon(QtGui.QIcon(os.path.join(ADFRICONPATH, "eye-icon.png")))
        #self.visibleButton.setToolTip("show/hide trg object")

        #self.visibleButton.clicked.connect(self.showHideTrg)
        self.boxButton = QtGui.QPushButton()
        self.boxButton.setIcon(QtGui.QIcon(os.path.join(ADFRICONPATH, "cube1.png")))
        self.boxButton.setCheckable(True)
        self.boxButton.setChecked(False)
        self.boxButton.clicked.connect(self.showHideBox)
        self.boxButton.setToolTip("show/hide box")
        buttonBoxLayout.addWidget(self.addButton)
        buttonBoxLayout.addWidget(self.removeButton)
        #buttonBoxLayout.addWidget(self.visibleButton)
        buttonBoxLayout.addWidget(self.boxButton)
        self.buttonBox.setLayout(buttonBoxLayout)
        self.vlayout1.addWidget(self.buttonBox)
        
        self.mainlayout = QtGui.QHBoxLayout()
        self.splitter = splitter = QtGui.QSplitter()
        self.mainlayout.addWidget(splitter)
        self.leftPaneWidget = QtGui.QWidget()
        self.leftPaneWidget.setLayout(self.vlayout1)
        self.rightPaneWidget = QtGui.QWidget()
        self.rightPaneWidget.setLayout(self.vlayout2)
        splitter.addWidget(self.leftPaneWidget)
        splitter.addWidget(self.rightPaneWidget)
        #self.mainlayout.addLayout(self.vlayout1)
        #self.mainlayout.addLayout(self.vlayout2)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 3)
        self.setLayout(self.mainlayout)
        
    def buildTreeMenu(self, point):
        # right-click on a tree item pops up a menu 
        item = self.treeWidget.itemAt(point)
        #print "ITEM", item.text(0)
        parent = item.parent()
        menu = QtGui.QMenu(self)
        obj = item._object
        #vis = item._visible
        vis = False
        # visStr will appear in the menu
        visStr = "Show"
        if list(item.foreground(0).color().getRgb()[:3]) == [0,0,0]:
            #item is currently visible (the font is black)
            vis = True
            visStr = "Hide"
        if isinstance(obj, MapsFile):
            action = QtGui.QAction("Add/remove map types", self)
            action.triggered.connect(CallbackFunction(self.mapTypeSelector_cb, item))
            menu.addAction(action)
            action = QtGui.QAction("%s trg"%visStr, self)
            menu.addAction(action)
            action.triggered.connect(CallbackFunction(self.showHideTrg, item, not vis))
            #import pdb; pdb.set_trace()
        else:
            if isinstance(obj, GridMap):
                nch = item.childCount()
                if nch > 0:
                    action = QtGui.QAction("%s map"%visStr, self)
                    menu.addAction(action)
                    action.triggered.connect(CallbackFunction(self.showHideMap, item, not vis))   # isocontour menu action
                # orthoslice menu action
                directions = ["X", "Y", "Z"]
                #find all ortho_ items in the tree (if any of 'ortho_X' or 'ortho_Y' or 'orto_Z' exists exclude it from drop down menu)
                isoPos = nch
                nisoItems = 0
                orthoPos = nch
                for i in range(nch):
                    itchild = item.child(i)
                    txt = str(itchild.text(0))
                    #print "menu item:", txt
                    if txt.startswith("ortho"):
                        dd = txt.split("_")[1]
                        if dd in directions:
                            directions.remove(dd)
                        orthoPos = i+1
                    elif txt.startswith("iso"):
                        nisoItems += 1
                        isoPos = i+1
                action = QtGui.QAction("Add isocontour", self)
                action.triggered.connect(CallbackFunction(self.addIsocontour, item, isoPos, nisoItems+1))
                menu.addAction(action)
                if len(directions):
                    submenu = menu.addMenu("Add orthoslice")
                    for txt in directions:
                        action = QtGui.QAction(txt, self)
                        action.triggered.connect(CallbackFunction(self.addOrthosice, item, txt, orthoPos))
                        submenu.addAction(action)
                        
            elif isinstance(obj, IndexedPolygons): # iso or ortho geometry
                action = QtGui.QAction(visStr, self)
                menu.addAction(action)
                action.triggered.connect(CallbackFunction(self.showHideGeomObject, item, not vis))
        menu.popup(QtGui.QCursor.pos())
        #menu.popup(self.treeWidget.viewport().mapToGlobal(point))
        
    def onItemDoubleClick(self, item, column, expandCollapse=True):
        #show/hide item
        if expandCollapse:
            if item.isExpanded():
                self.treeWidget.collapseItem(item)
            else:
                self.treeWidget.expandItem(item)
        if not hasattr(item, "_object"): return
        visible = False
        if list(item.foreground(0).color().getRgb()[:3]) == [0,0,0]:
            visible = True
        if isinstance(item._object, MapsFile): #trgitem
            self.showHideTrg(item, not visible)
        elif isinstance(item._object,  GridMap): # map item
            self.showHideMap(item, not visible)
        elif isinstance(item._object, IndexedPolygons): #iso_# item
            self.showHideGeomObject(item, not visible)

    def setItemBold(self, item):
        font = item.font(0)
        font.setBold(True)
        item.setFont(0, font)
        
    def openTrgFile(self):
        # get trg file from file dialog
        filename, selfilter = QtGui.QFileDialog().getOpenFileName(
            self, self.tr("Open trg file"), '',
            self.tr("trg Files (*.trg);; All files (*)"))
        return str(filename)

    def addTrg(self, filename=None, boxGeom=None):
        # add trg item to the tree widget
        if not filename:
            filename = self.openTrgFile()
        if not filename: return            
        trgname = os.path.splitext(os.path.basename(filename))[0]
        item = QtGui.QTreeWidgetItem(self.treeWidget, [trgname])
        self.setItemBold(item)
        item.dummyChild = QtGui.QTreeWidgetItem(item)
        item.dummyChild._object = None
        if self.trg_data.has_key(trgname) and self.trg_data[trgname]['filename']==os.path.abspath(filename):

            _item = self.treeWidget.findItems(trgname, QtCore.Qt.MatchExactly, column=0)
            if len(_item) :
                self.treeWidget.setCurrentItem(_item[0], 0)

                reply = QtGui.QMessageBox.question(self, "Trg file exists",
                         "Trg  %s is in the gui. Replace it?"% trgname,
                          QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
                if reply == QtGui.QMessageBox.Yes:
                    self.deleteTrgItem(ask=False)
                else:
                    return
        data = self.getTargetData(filename.encode('ascii', 'replace'))
        
        self.trg_data[trgname] = {'filename':os.path.abspath(filename),
                                  'mapTypes':data['mapTypes']}
        mf = MapsFile(filename)
        item._object = mf
        item._visible = True
        if self.viewer:
            if not boxGeom:
                boxGeom = NiceBox('%s_box'%trgname)
                boxGeom.addToViewer(self.viewer)
                center = data['boxCenter']
                boxGeom.setCenter(*center)
                sides = data['boxLengths']
                boxGeom.setSides(*sides)
                boxGeom.edges.Set(materials=[[0., 0.6, 1.,]], inheritMaterial=0)
                boxGeom.corners.Set(materials=[[0., 0.6, 1.,]], inheritMaterial=0)
            self.trg_data[trgname]['boxGeom'] = boxGeom
        #print "DATA", data

    def setBoxVisible(self, boxGeom, val=None):
        if val is None:
            val = not boxGeom.master.visible
        boxGeom.master.Set(visible=val)
        for c in boxGeom.master.children:
            if c.name=='faces':
                c.Set(visible = 0)
            else:
                c.Set(visible = val)
        

    def onItemExpanded (self, item):
        # called when item's plus icon is clicked
        #print 'Expand XXXXXX', item, obj, hasattr(item, 'dummyChild')
        obj = item._object
        if isinstance (obj, MapsFile):
            if hasattr(item, 'dummyChild'):
                item.removeChild(item.dummyChild)
                del item.dummyChild
                if item.childCount() == 0:
                    obj = item._object
                    trg_name = item.text(0)
                    trg_maptypes = self.trg_data[trg_name]['mapTypes']
                    items = []
                    for mt in self.default_maptypes:
                        if mt in trg_maptypes: 
                            items.append(self.addMap(mt, item))

    def addMap(self, mapType, parent):
        # add map item to the tree widget
        item = QtGui.QTreeWidgetItem(parent, [mapType])
        from ADFRcc.adfr import GridMap
        item._object = GridMap() # create a "dummy" object without reading any
        # maps yet
        trgName = str(parent.text(0))
        if not self.iso_data.has_key(trgName):
            self.iso_data[trgName] = {}
        item._visible = True
        self.setItemBold(item)
        return item
    
    # Isocontour GUI
    def addIsocontour(self, parent, pos, isonum):
        # add isocontour item to the tree widget
        mapObj = parent._object
        mapType = str(parent.text(0)) 
        trgname = str(parent.parent().text(0))
        if not mapObj.getMapType(): #  this is a "dummy" object. The map file has not been read. Set the map type and read the map file:
            mfObj = parent.parent()._object
            mapObj =  mfObj.getMap(mapType)
            parent._object = mapObj
            mapObj._gridData = mapObj.getGridDataPy().astype(numpy.float32)
            self.trg_data[trgname][mapType] = mapObj #GridMap() instance
        # Compute isocontour and get the geometry 
        self.computeIsocountour(trgname, mapType, mapObj)
        data = mapObj._gridData
        minval = float(data.min())
        maxval = 0. # setting this to 0 (real data max values can be very large, the user can change this value in the GUI)
        isoval = 0.5*minval
        if minval == 0:
            maxval = float(data.max())
            isoval = maxval
        if isoval > maxval: isoval = maxval
        #print "addIsocontour: isoval", isoval, type(isoval) 
        g = self.getIsocontourGeom(trgname, mapType, isoval, maxval, geomName="isocontour_%s_%s_%d"%(trgname, mapType, isonum))
        g._isoValue = isoval
        # add iso item to the tree under the parent at index pos
        item = QtGui.QTreeWidgetItem()
        item.setText(0, "iso_%d"%isonum)
        parent.insertChild(pos,item)
        self.treeWidget.setCurrentItem(item, 0)
        item._object = g
        item._visible = True
        self.numIsoItems += 1
        self.setItemBold(item)
        if self.viewer:
            self.addIsoGeomToViewer(trgname, mapType, g, parent=None, addBox=False)
            boxGeom = self.trg_data[trgname]['boxGeom']
            if not self.boxButton.isChecked():
                self.boxButton.setChecked(True)
            self.setBoxVisible(boxGeom, True)
        # create or update  geom parameters widgets
        if not self.isoGeomPropWidget:
            self.isoGeomPropWidget = IsoGeomProp(data=data, geom=g, app=self)
            self.propWidget.addWidget(self.isoGeomPropWidget)
            self.stackedWidgets["iso"] = self.isoGeomPropWidget
        else:
            self.isoGeomPropWidget.setEnabled(g.visible)
            self.isoGeomPropWidget.setGeom(g)
            self.isoGeomPropWidget.updateIsoData(data, maxval=maxval, isovalue=isoval)
        self.propWidget.setCurrentWidget(self.isoGeomPropWidget)
        self.treeWidget.expandItem(parent)
        #if mapitem is not visible(grayed out), make it visible:
        if list(parent.foreground(0).color().getRgb()[:3]) != [0,0,0]:
            parent.setForeground(0, QtGui.QColor(0,0,0)) #make it black
            parent._visible = True
        return item

    def computeIsocountour(self, trgname, mapType, mapObj):
        data = mapObj.getGridDataPy().astype(numpy.float32)
        spacing = mapObj.getDistBetweenGridPoints()
        newgrid3D = numpy.ascontiguousarray(
            numpy.reshape( numpy.transpose(data),
                           (1, 1)+tuple(data.shape) ) , 'f')
        origin = mapObj.getOriginPy().astype('f')
        stepsize = numpy.array([spacing, spacing, spacing]).astype('f')
        #print "########## Computing Isocontour. mapType : %s, spacing %f, data size:" %(mapType, spacing), data.shape, origin
        dataSet = isocontour.newDatasetRegFloat3D(newgrid3D, origin, 
                                                            stepsize)
        iso_data = self.iso_data[trgname]
        if iso_data.has_key(mapType):
            isocontour.delDatasetReg(iso_data[mapType][0])
        # have to keep newgrid3D (iso_data structure points to newgrid3D ) 
        iso_data[mapType] = [dataSet, newgrid3D]

    def getIsocontourGeom(self, trgname, mapType, isoval, maxval, geomName, material=None):
        # Returns IndexedPolygons geometry for mapType and specified iso value
        # (if the isocontour has been computed for this name).
        # Returns None if there is no isodata for the mapType.
        #import pdb; pdb.set_trace()
        g = None
        vert, norm, tri = self.getIsoData(trgname, mapType, isoval)
        #print tri[:10], norm[:10]
        #tmp = tri[:,0].copy()
        #tri[:,0] = tri[:,2]
        #tri[:,2] = tmp
        #norm = -norm
        #print tri[:10], norm
        if vert is not None:
            #print " ########### creating geometry", geomName
            g = IndexedPolygons(geomName, inheritMaterial=False, inheritLineWidth=0, inheritFrontPolyMode=False)
            g.Set(inheritSharpColorBoundaries=False, sharpColorBoundaries=True)
            g.Set(vertices=vert,vnormals=norm,faces=tri, inheritMaterial=False, inheritFrontPolyMode=False)
            g._isoValue = isoval
            g._maxIsoValue = maxval
            if material:
                g.inheritMaterial = False
                g.Set(materials=material)
            else:
                el = ElementLookup.get(mapType, mapType)
                col = AtomElements.get(el, None)
                if not col is None:
                    g.inheritMaterial = False
                    g.Set(materials=[col])
            #if vert is not None:
            #    g.sortPoly()
        return g

    def setIsoValue(self, geom, isoval):
        names = geom.name.split("_") # should be [trgname, mapType, index]
        mapType = names[-2]
        trgnames = names[1:-2]
        if len(trgnames)>1 : # trgname contained "_"
            from string import join
            trgname = join(trgnames, "_")
        else: trgname = trgnames[0]
            
        vert, norm, tri = self.getIsoData(trgname, mapType, isoval)
        #tmp = tri[:,0].copy()
        #tri[:,0] = tri[:,2]
        #tri[:,2] = tmp
        #norm = -norm
        if vert is not None:
            geom.Set(vertices=vert,vnormals=norm,faces=tri, inheritMaterial=False, inheritFrontPolyMode=False)

    def getIsoData(self, trgname, mapType, isoval):
        iso_data = self.iso_data[trgname]
        if not iso_data.has_key(mapType):
            print "No map data is available for %s"% mapType
            return None, None, None
        dataSet = iso_data[mapType][0]
        isoc = isocontour.getContour3d(dataSet, 0, 0, isoval,
                                           isocontour.NO_COLOR_VARIABLE)
        col = numpy.zeros((isoc.nvert)).astype('f')
        vert = numpy.zeros((isoc.nvert,3)).astype('f')
        norm = numpy.zeros((isoc.nvert,3)).astype('f')
        tri = numpy.zeros((isoc.ntri,3)).astype('i')
        isocontour.getContour3dData(isoc, vert, norm, col, tri, 0)
        tmp = tri[:,0].copy()
        tri[:,0] = tri[:,2]
        tri[:,2] = tmp
        norm = -norm
        return vert, norm , tri

    # ortho slice gui
    def addOrthosice(self, parent, direction, pos):
        mapObj = parent._object
        mapType = str(parent.text(0))
        trgname = str(parent.parent().text(0))
        if not mapObj.getMapType(): #  this is a "dummy" object. The map file has not been read. Set the map type and read the map file:
            mfObj = parent.parent()._object
            mapObj =  mfObj.getMap(mapType)
            parent._object = mapObj
            self.trg_data[trgname][mapType] = mapObj

        g = textured2DArray('OrthoSlice_%s_%s_%s'%(direction, trgname, mapType),
                                       inheritLighting=False,
                                       lighting=False)
        slice, verts, minval, maxval = self.get2DOrthoSlice(mapObj, direction, 0)
        g.Set(vertices=verts, array=slice)
        g._trgname = trgname
        g._maptype = mapType
        g._orthoSliceNum = {"X":0, "Y":0, "Z":0}
        #g._orthoValRanges = {"X":[], "Y":[], "Z":[]}
        # get number of slices in specified direction
        nslices = mapObj._gridData.shape[{'X':0, 'Y':1, 'Z':2}.get(direction)]
        rangeList = [minval, maxval]
        w_key = trgname + " "+mapType+" " + direction
        if not self.orthoGeomPropWidgets.get(w_key, None):
            # create orthoslice GUI
            w = OrthoSliceGui(callback=self.setOrthoColorMap, direction=direction, nslices=nslices, valrange=rangeList, geom=g)
            self.propWidget.addWidget(w)
            self.orthoGeomPropWidget = w
            self.orthoGeomPropWidgets[w_key] = w
            ramp = w.getColorRamp()
            self.setOrthoColorMap(geom=g, rgbaRamp=w.getColorRamp(),
                                  rangeList=rangeList)
        else:
            w = self.orthoGeomPropWidgets[w_key]
            w.setGeom(g)
        self.propWidget.setCurrentWidget(w)

        # add ortho item to the tree under the parent at index pos
        #item = QtGui.QTreeWidgetItem(parent, ["ortho_%s"%direction])
        item = QtGui.QTreeWidgetItem()
        item.setText(0, "ortho_%s"%direction)
        parent.insertChild(pos,item)
        self.treeWidget.setCurrentItem(item, 0)
        item._object = g
        item._visible = True
        item._orthowidget = w
        self.numOrthoItems += 1
        self.setItemBold(item)
        self.treeWidget.expandItem(parent)
        #if mapitem is not visible(grayed out), make it visible:
        if list(parent.foreground(0).color().getRgb()[:3]) != [0,0,0]:
            parent.setForeground(0, QtGui.QColor(0,0,0)) #make it black
            parent._visible = True
        if self.viewer:
            self.viewer.AddObject(g, parent=None)
        #return item
        
    def get2DOrthoSlice(self, mapObj, axis, sliceNum):
        """Return a 2D array corresponding to an axis aligned slice orthogonal
        to axis at index along this axis.  We also return the real vertices
        corresponding the corner of the quad for this slice.
        axis cand be 'x', 'y' or 'z'
        sliceNum varies from 0 to max(extent(axis))
        """

        if hasattr(mapObj, "_gridData"):
            data = mapObj._gridData
        else:
            data = mapObj.getGridDataPy().astype(numpy.float32)
            mapObj._gridData = data
        spacing = mapObj.getDistBetweenGridPoints()
        stepsize = numpy.array([spacing, spacing, spacing]).astype('f')
        origin = mapObj.getOriginPy().astype('f')
        minval = float(data.min())
        if minval < 0:
            maxval = -minval
        else:
            maxval = float(data.max())
        # this following code is taken from Volume/Grid3D.py            
        if axis=='X': axisInd = 0
        elif axis=='Y': axisInd = 1
        elif axis=='Z': axisInd = 2
        dims = data.shape
        dx, dy, dz = (dims[0]-1,dims[1]-1,dims[2]-1)
        sx, sy, sz = stepsize
        ox, oy, oz = origin
        ex, ey, ez = (ox+sx*dx, oy+sy*dy, oz+sz*dz)
        if axis=='X':
            slice = numpy.array(data[sliceNum,:,:])
            x = ox + (sx*0.5)+(sliceNum*sx)
            vertices = [ [x,oy,oz], [x,ey,oz], [x,ey,ez], [x,oy,ez] ]

        elif axis=='Y':
            slice = numpy.array(data[:,sliceNum,:])
            y = oy + (sy*0.5)+(sliceNum*sy)
            #vertices = [ [ox,y,oz], [ox,y,ez], [ex,y,ez], [ex,y,oz] ]
            vertices = [ [ox,y,oz], [ex,y,oz], [ex,y,ez], [ox,y,ez] ]

        elif axis=='Z':
            slice = numpy.array(data[:,:,sliceNum])
            z = oz + (sz*0.5)+(sliceNum*sz)
            vertices = [ [ox,oy,z], [ex,oy,z], [ex,ey,z], [ox,ey,z] ]

        return slice, vertices, minval, maxval

    def setOrthoColorMap(self, geom=None, direction=None, rgbaRamp=None, rangeList=None, slicen=None):
        #rangeList = [minval, maxval]
        # callback of the OrthoSlice widget
        #print "!!!!!! setOrthoColorMap, geom", geom.name, "min-max:", rangeList, "slice:", slicen, "ramp len:",
        #if rgbaRamp is not None: print len(rgbaRamp)
        
        if not geom:
            geom = self.orthoGeomPropWidget.geom
        if geom == None: return
        if direction is None:
            direction = self.orthoGeomPropWidget.direction
        if slicen is not None:
            mapObj = self.trg_data[geom._trgname][geom._maptype]
            slice, verts, minval, maxval = self.get2DOrthoSlice(mapObj, direction, slicen)
            geom.Set(vertices=verts, array=slice)
        colormap = geom.colormap
        if rgbaRamp is not None:
            if rangeList is not None:
                minval = rangeList[0]
                maxval = rangeList[1]
            else:
                minval = 'not passed'
                maxval = 'not passed'
            colormap.configure(ramp=rgbaRamp, mini=minval, maxi=maxval)
            geom.Set(colormap=colormap)
        else:
            if rangeList is not None:
                colormap.configure(mini=rangeList[0], maxi=rangeList[1])
                geom.Set(colormap=colormap)
        #print "IN CALLBACK min and max", geom.colormap.mini, geom.colormap.maxi
    
    # general GUI events 
    def closeEvent(self, evnt=None):
        self.closedSignal.emit()
        
    def isItemVisible(self, item):
        pass

    def onItemClick(self):
        # callback of the tree widget item
        item = self.treeWidget.currentItem()
        disable = True
        if hasattr(item, "_object"): 
            if isinstance(item._object, MapsFile): #trgitem
                disable = False
                if self.isoGeomPropWidget:
                    self.isoGeomPropWidget.setEnabled(False)
                if self.orthoGeomPropWidget:
                    self.orthoGeomPropWidget.setEnabled(False)
            elif isinstance(item._object,  GridMap): # map item
                if self.isoGeomPropWidget:
                    self.isoGeomPropWidget.setEnabled(False)
                if self.orthoGeomPropWidget:
                    self.orthoGeomPropWidget.setEnabled(False)
            elif isinstance(item._object, IndexedPolygons): #iso or ortho item
                mapitem = item.parent()
                data = mapitem._object._gridData #getGridDataPy().astype(numpy.float32)
                geom = item._object
                if isinstance(item._object, textured2DArray): # ortho
                    w = item._orthowidget
                    newDir = item.text(0).split("_")[1]
                    self.orthoGeomPropWidget = w 
                    w.setEnabled(geom.visible)
                    self.propWidget.setCurrentWidget(w)
                else: # iso_## item 
                    self.propWidget.setCurrentWidget(self.isoGeomPropWidget)
                    self.isoGeomPropWidget.setEnabled(geom.visible)
                    self.isoGeomPropWidget.updateIsoData(data, maxval=geom._maxIsoValue, isovalue=geom._isoValue)
                    self.isoGeomPropWidget.setGeom(geom)
        self.removeButton.setDisabled(disable)
        #self.visibleButton.setDisabled(disable)

    def mapTypeSelector_cb(self, item):
        # pop up a widget to select map types
        trg_name = item.text(0)
        trg_maptypes = self.trg_data[trg_name]['mapTypes']
        mtDict = {}.fromkeys(trg_maptypes, False)
        added_maptypes = []
        for i in range(item.childCount()):
            ch = item.child(i)
            if ch._object:
                mt = ch.text(0)
                mtDict[mt] = True
                added_maptypes.append(mt)
        # create dialog with checkbuttons to select map types
        dialog = MapTypeChooser(mtDict, self)
        result = dialog.exec_()
        if not result: return
        atypesstr = dialog.getTypesString() 
        mapTypes = atypesstr.split()
        mapTypes.sort()
        # dictionary containing maptypes to add to the tree (with True value)
        newDict = {}.fromkeys(mapTypes, True)
        for mt in added_maptypes:
            if newDict.has_key(mt):
                newDict.pop(mt)
            else:
                newDict[mt] = False
        for mt, val in newDict.items():
            if val:
                self.addMap(mt, item)
            else:
                self.removeMap(str(mt), item)
        self.treeWidget.expandItem(item)
        #import pdb; pdb.set_trace()

    def addIsoGeomToViewer(self, trgname, mapName, g, parent=None, addBox=False):
        vi = self.viewer
        if addBox:
            from DejaVu2.Box import Box
            mapObj = self.maps[mapName]
            box = Box('BoundingBox')
            pt1 = mapObj.getOriginPy().astype('f')
            dims = mapObj._gridData.shape #getGridDataPy().shape
            stepsize = mapObj.getDistBetweenGridPoints()
            pt2 = [pt1[0]+(stepsize*(dims[0]-1)),
                   pt1[1]+(stepsize*(dims[1]-1)),
                   pt1[2]+(stepsize*(dims[2]-1))]
            #print "Box corners", pt1, pt2
            box.Set(cornerPoints=(pt1, pt2))
            vi.AddObject(box,parent=parent)

        vi.AddObject(g, parent = parent)
        vi.Redraw()

    def removeMap(self, mapType, trgItem):
        for i in range(trgItem.childCount()):
            mapItem = trgItem.child(i)
            if not mapItem: continue
            if mapItem.text(0) == mapType:
                mapObj = mapItem._object
                #print "removing map", mapType
                if mapObj:
                    del mapObj
                for i in range(mapItem.childCount()):
                    geomItem = mapItem.child(i)
                    if hasattr(geomItem, "_object"):
                        if isinstance(geomItem._object, textured2DArray):
                            w = geomItem._orthowidget
                            self.propWidget.removeWidget(w)
                            self.numOrthoItems -=1
                        else:
                            self.numIsoItems -= 1
                        self.viewer.RemoveObject(geomItem._object)
                        self.viewer.Redraw()
                trgItem.removeChild(mapItem)
                trgName = trgItem.text(0)
                if self.trg_data[trgName].has_key(mapType):
                    self.trg_data[trgName].pop(mapType)
                if self.iso_data.has_key(trgName):
                    if self.iso_data[trgName].has_key(mapType):
                        self.iso_data[trgName].pop(mapType)
                if self.numIsoItems == 0 and self.isoGeomPropWidget:
                    self.propWidget.removeWidget(self.isoGeomPropWidget)
                    self.isoGeomPropWidget = None
                if self.numOrthoItems == 0 and self.orthoGeomPropWidget:
                    self.orthoGeomPropWidget = None
                return
        
    def getTargetData(self, filename):
        #get some info from the target file without unzipping it                        
        from zipfile import ZipFile
        import pickle
        zf = ZipFile(filename, 'r')
        dataFile = None
        for name in zf.namelist():
            if os.path.basename(name)=='data.pkl':
                dataFile = name
                break
        if dataFile:
            return pickle.loads(zf.read(dataFile))
        else:
            return None

    def deleteTrgItem(self, ask=True):
        # callback of the "cross" button.
        # If a trg item is selected in the tree widget- remove it.
        item = self.treeWidget.currentItem()
        if item and hasattr(item, "_object"):
            if isinstance(item._object, MapsFile):
                # this is a trg item , pop up "confirm delete" dialog
                trgName = item.text(0)
                if ask:
                    reply = QtGui.QMessageBox.question(self, "Delete",
                                                       "do you want to delete %s ?"% trgName,
                                                       QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
                else:
                    reply = QtGui.QMessageBox.Yes
                if reply == QtGui.QMessageBox.Yes:
                    itemInd = self.treeWidget.indexOfTopLevelItem(item)
                    #print "Removing item %d %s from the tree" % (itemInd, trgName)
                    #import pdb; pdb.set_trace()
                    for i in range(item.childCount()):
                        mapItem = item.child(i)
                        mapObj = mapItem._object
                        if mapObj:
                            del mapObj
                        if self.viewer:
                            for i in range(mapItem.childCount()):
                                geomItem = mapItem.child(i)
                                if hasattr(geomItem, "_object"):
                                    if isinstance(geomItem._object, textured2DArray):
                                        w = geomItem._orthowidget
                                        self.propWidget.removeWidget(w)
                                        self.numOrthoItems -=1
                                    else:
                                        self.numIsoItems -= 1
                                    self.viewer.RemoveObject(geomItem._object)
                                    self.viewer.Redraw()
                    del item._object
                    trgFile = self.trg_data[trgName]['filename']
                    self.trg_data.pop(trgName)
                    if self.iso_data.has_key(trgName):self.iso_data.pop(trgName) 
                    self.treeWidget.takeTopLevelItem(itemInd)
                    if self.numIsoItems == 0 and self.isoGeomPropWidget:
                        self.propWidget.removeWidget(self.isoGeomPropWidget)
                        self.isoGeomPropWidget = None
                    if self.numOrthoItems == 0 and self.orthoGeomPropWidget:
                        self.orthoGeomPropWidget = None
                    # When the TrgMapGui is closed and reopened again the object tree is
                    # updated from the calling GUI (self.app) list of trg filenames.
                    # We do not want the deleted item appear when the gui is reopened,
                    # so we need to remove the trg filename from this list.
                    if self.app and hasattr(self.app, 'trgFiles'):
                        if trgFile in self.app.trgFiles:
                            self.app.trgFiles.remove(trgFile)

    def showHideTrg(self, item, show):
        if not item: return
        trgName = item.text(0)
        geoms = []
        mapItems = []
        item._visible = show
        #import pdb; pdb.set_trace()
        for i in range(item.childCount()):
            mapItem = item.child(i)
            if not hasattr (mapItem, "_visible"):
                continue
            if mapItem._visible:
                self.showHideMapGeoms(mapItem, show=show)
        boxGeom = self.trg_data[trgName]['boxGeom']
        if show:
           item.setForeground(0, QtGui.QColor(0,0,0)) # make text black
        else:
           item.setForeground(0, QtGui.QColor(125, 125, 125)) # make text grey
        self.setBoxVisible(boxGeom, show)


    def showHideMap(self, mapItem, show):
        if not mapItem:return
        if mapItem.childCount() == 0: return
        mapItem._visible = show
        self.showHideMapGeoms(mapItem, show=show)
        #see if the parent becomes visible or not
        self.updateTreeVisible(mapItem)
                
    def showHideGeomObject(self, item, show):
        #show/hide isocontour or orthoslice geometry
        if not item:return
        geom = item._object
        geom.Set(visible=show)
        item._visible = show
        if show:
            item.setForeground(0, QtGui.QColor(0,0,0)) # make text black
        else:
           item.setForeground(0, QtGui.QColor(125, 125, 125)) # make text grey
        if isinstance(geom, textured2DArray):
            self.orthoGeomPropWidget.setEnabled(show)
        else:
            self.isoGeomPropWidget.setEnabled(show)
        self.updateTreeVisible(item)

    def updateTreeVisible(self, item):
        parent = item.parent()
        while parent:
            show = False
            for i in range(parent.childCount()):
                item = parent.child(i)
                if item._visible:
                    show = True
                    break
            if show:
                parent.setForeground(0, QtGui.QColor(0,0,0)) # make text black
            else:
                parent.setForeground(0, QtGui.QColor(125, 125, 125)) # make text grey
            parent._visible = show
            parent = parent.parent()

        
    def getMapGeoms(self, mapItem):
        # return a list of all geometries under mapItem in the tree
        # and a visible flag. Visible is set to True if at least one geometry in
        # the list is visible
        visible = False
        geoms = []
        for i in range(mapItem.childCount()):
            geomItem = mapItem.child(i)
            if hasattr(geomItem, "_object"): 
                if geomItem._object.visible:
                    visible = True
                geoms.append(geomItem._object)
        return geoms, visible

    def showHideMapGeoms(self, mapItem, show=True):
        trgItem = mapItem.parent()
        mapVisible = mapItem._visible # tag of the map item in the tree
        mapShow = False
        ngeoms = 0
        for i in range(mapItem.childCount()):
            geomItem = mapItem.child(i)
            if hasattr(geomItem, "_object"):
                itemVisible = geomItem._visible #tag of the iso item int tree 
                geom = geomItem._object
                ngeoms += 1
                if show:
                    if itemVisible and not geom.visible:
                        geom.Set(visible=True)
                        geomItem.setForeground(0, QtGui.QColor(0,0,0)) # make text black
                        mapShow = True
                else:
                    if itemVisible and geom.visible:
                        geom.Set(visible=False)
                        geomItem.setForeground(0, QtGui.QColor(125, 125, 125)) # make text grey
        if ngeoms == 0: mapShow = show
        if mapShow:
            #set the text of hte mapItem black
            mapItem.setForeground(0, QtGui.QColor(0,0,0))
        else: # grey
          mapItem.setForeground(0, QtGui.QColor(125, 125, 125))
        return mapShow
    
    def showHideBox(self):
        item = self.treeWidget.currentItem()
        parent = item.parent()
        while parent is not None:
            item = parent
            parent = item.parent()
        trgname = item.text(0)
        boxGeom = self.trg_data[trgname].get('boxGeom', None)
        val = self.boxButton.isChecked()
        if boxGeom:
            # toggle visibility of the box
            self.setBoxVisible(boxGeom, val)

class MapTypeChooser(QtGui.QDialog):

    def __init__(self, mapList, parent=None):
        super(MapTypeChooser, self).__init__(parent)
        self.setWindowTitle("maps shown in tree")
        layout =  QtGui.QVBoxLayout()
        layout1 =  QtGui.QHBoxLayout()
        self.allButtonOn = QtGui.QCheckBox("select all")
        self.allButtonOn.clicked.connect(self.selectAll)

        self.allButtonOff = QtGui.QCheckBox("deselect all")
        self.allButtonOff.clicked.connect(self.deselectAll)
        layout1.addWidget(self.allButtonOn)
        layout1.addWidget(self.allButtonOff)
        layout.addLayout(layout1)
        self.checkbuttonGroup = QtGui.QGroupBox()
        self.checkbuttonGroup.setCheckable(False)
        buttonGroupLayout = QtGui.QGridLayout()
        self.mapList = mapList
        ncol = 8
        row = 0
        col = 0
        self.items = []
        mtypes = mapList.keys()
        mtypes.sort()
        for count, mtype in enumerate(mtypes):
            row = count/ncol
            bb = QtGui.QCheckBox(mtype)
            bb.setChecked(mapList[mtype])
            self.items.append(bb)
            buttonGroupLayout.addWidget(bb, row, col)
            bb.clicked.connect(self.itemClicked_cb)
            col += 1
            if col == ncol: col = 0
        self.checkbuttonGroup.setLayout( buttonGroupLayout)
        layout.addWidget(self.checkbuttonGroup)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.cancel_cb)
        self.okbutton = self.buttonBox.button(QtGui.QDialogButtonBox.Ok)
        layout.addWidget(self.buttonBox)
        self.setLayout(layout)

    def selectAll(self):
        for item in self.items:
            item.setCheckState(QtCore.Qt.Checked)
        self.okbutton.setDisabled(False)

    def deselectAll(self):
        for item in self.items:
            item.setCheckState(QtCore.Qt.Unchecked)
        self.okbutton.setDisabled(True)
        
    def itemClicked_cb(self):
        nitems = len(self.items)
        count = 0
        for item in self.items:
            if item.checkState()==QtCore.Qt.Unchecked:
                count+=1
        if count == nitems: #all unchecked, disable "OK" button
            self.okbutton.setDisabled(True)
        else:
            self.okbutton.setDisabled(False)

    def getTypesString(self):
        s = ''
        for item in self.items:
            if item.checkState()==QtCore.Qt.Checked:
                s += '%s '%item.text()
        return s[:-1]
    
    def cancel_cb(self):
        n = 0
        for atype, checked in self.mapList.items():
            item = self.items[n]
            if checked:
                item.setCheckState(QtCore.Qt.Checked)
            else:
                item.setCheckState(QtCore.Qt.Unchecked)
            n+=1
        self.reject()
        

        
import math
    
class Histogram(QtGui.QWidget):
    valChanged = QtCore.Signal(float)
    
    def __init__(self, parent=None, data=None, isovalue=0.0):
        super(Histogram, self).__init__()
        self.rwidth = self.minwidth = 220
        self.rheight = 60
        self.setMinimumWidth(self.rwidth+20)
        self.setMinimumHeight(self.rheight+15)
        self.parent = parent
        self.data = None
        self.histData = None
        self.rx = 10
        self.ry = 15
        self.drawLine = False
        self.pos = None
        self.setMouseTracking(True)
        self.minVal = self.maxVal = 0.0
        self.isovalue = isovalue
        if data is not None:
            self.setData(data, redraw=False)
        
    def setData(self, data, minval=None, maxval=0., isovalue=None, redraw=True):
        #set histogram data
        if minval == None:
            minval = data.min()
        if maxval == None:
            maxval == data.max()
        if minval == maxval == 0:
            maxval = data.max()
        self.data = data
        self.histData = numpy.histogram(data.flat,bins=self.rwidth, range=[minval, maxval])
        self.minVal = minval
        self.maxVal = maxval
        if isovalue is not None:
            self.setIsoVal(isovalue, redraw=False)
            #self.isovalue = isovalue
        if redraw:
            self.update()

    def getValFromPos(self, x):
        dd = (self.maxVal - self.minVal)/self.rwidth
        return (x - self.rect.left())*dd + self.minVal

    def setMaxValue(self, val, redraw=True):
        self.histData = numpy.histogram(self.data.flat,bins=self.rwidth, range=[self.minVal, val])
        self.maxVal = val
        if redraw:
            self.update()
        
    def mousePressEvent(self, event):
        x = event.x()
        y = event.y()
        if self.isInsideRect(x,y):
            self.pos = [x, y]
            self.isovalue = self.getValFromPos(x)
            self.update()
            self.drawLine = True
        else:
            self.drawLine = False

    def mouseReleaseEvent(self, event):
        self.drawLine = False

    def mouseMoveEvent(self, event):
        if not self.drawLine: return
        #distance_from_center = round(((event.y() - 250)**2 + (event.x() - 500)**2)**0.5)
        #self.label.setText('Coordinates: ( %d : %d )' % (event.x(), event.y()) + "Distance from center: " + str(distance_from_center))
        x = event.x()
        y = event.y()
        if self.isInsideRect(x,y):
            self.pos = [x, y]
            self.isovalue = self.getValFromPos(x)
            self.update()
            self.valChanged.emit(self.isovalue)

    def isInsideRect(self, x, y):
        xLeft = self.rect.left()
        xRight = self.rect.right()
        yTop = self.rect.top()
        yBottom = self.rect.bottom()
        #print "isInsideRect, X:", x, xLeft, xRight, "Y:", y, yTop, yBottom 
        xinside = yinside = False
        if x >= xLeft and x <= xRight:
            xinside = True
        if y >= yTop and y <= yBottom:
            yinside = True
        if xinside and yinside : return True
        else: return False

    def setIsoVal(self, val, redraw=True):
        #compute the position of the current value line based on given val
        # redraw the histogram
        if val > self.maxVal or val < self.minVal: return
        y = self.rect.top()
        dd = (self.maxVal - self.minVal)/self.rwidth
        x = int((val-self.minVal)/dd + self.rect.left())
        self.pos = [x, y]
        self.isovalue = val
        if redraw:
            self.update()
            self.valChanged.emit(self.isovalue)
            
    def paintEvent(self, event):
        # Redraws the widget
        #print "PAINT EVENT", event
        if self.parent:
            # this is to resize the histogram rectangle when the parent widget's size  changes.
            w = self.parent.width() - self.rx*2 - 4
            if w > self.minwidth:
                if self.rwidth != w:
                    #resizing the histogram
                    self.rwidth = w
                    #print "RESIZING:", self.rwidth
                    if self.pos:
                        # get new coordinates for the vertical "iso value" line:
                        self.setIsoVal(self.isovalue, redraw=False)
            else:
                self.rwidth = self.minwidth
        painter = QtGui.QPainter()
        painter.begin(self)
        pen = QtGui.QPen()
        pen.setColor(QtCore.Qt.black)
        pen.setWidth(1)
        painter.setPen(pen)
        
        brush = QtGui.QBrush(QtCore.Qt.SolidPattern)
        brush.setColor("#DDDDDD")
        painter.setBrush(brush)
        self.rect = rect = QtCore.QRect(self.rx, self.ry, self.rwidth, self.rheight)        
        painter.drawRect(rect)
        xLeft = rect.left()
        xRight = rect.right()
        yTop = rect.top()
        yBottom = rect.bottom()
        height = rect.height()
        width = rect.width()
        #print "parent width", self.parent.width(), "rwidth", self.rwidth, "xr", xRight, "xl", xLeft
        # ---- Histogram itself ----------------------------------------------------
        if self.histData is None:
            painter.end()                    
            return
        _bins = self.histData[0]
        #import pdb; pdb.set_trace()
        nbBins = len(_bins)#_bins.size()

        # Find maximum height in bins unit
        heightMax = max(_bins)
        _heightMax = 1
        #print "heightMax", heightMax
        #Avoid giggling graph: do not update heightmax if variation <5%
        if abs(_heightMax-heightMax)/float(_heightMax) > 0.05:
            _heightMax = heightMax

        # Scale histogram from bins unit to pixels unit
        # handle upscaling and downscaling in a different way
        myPolygon = QtGui.QPolygon()
        linearGradient = QtGui.QLinearGradient(0, 0, 0, height)
        pen.setStyle(QtCore.Qt.SolidLine)

        if( nbBins < width ):
            wScale = width/float(nbBins)
            hScale = height/float(_heightMax)
            hScaleLog = height/math.log(float(_heightMax))

            pen.setColor("#00aaee")
            painter.setPen(pen)
            linearGradient.setColorAt(0.2, QtCore.Qt.white)
            linearGradient.setColorAt(1.0, "#00aaee")
            painter.setBrush(linearGradient)

            #bins
            pen.setColor("#016790")
            painter.setPen(pen)
            linearGradient.setColorAt(0.2, QtCore.Qt.white)
            linearGradient.setColorAt(1.0, "#016790")
            painter.setBrush(linearGradient)
            myPolygon.clear()
            myPolygon << QtCore.QPoint(xRight, yBottom) << QtCore.QPoint(xLeft, yBottom)
            for i in range(nbBins):
                myPolygon << QtCore.QPoint(xLeft+wScale*i, yTop+hScale*(_heightMax-_bins[i]))
            painter.drawPolygon(myPolygon)
        else :
            wScale = float(nbBins-1)/(width-1)
            hScale = height/float(_heightMax)
            hScaleLog = height/math.log(float(_heightMax))
            #bins
            pen.setColor("#016790")
            painter.setPen(pen)
            linearGradient.setColorAt(0.2, QtCore.Qt.white)
            linearGradient.setColorAt(1.0, "#016790")
            painter.setBrush(linearGradient)

            myPolygon.clear()
            myPolygon << QtCore.QPoint(xRight, yBottom) << QtCore.QPoint(xLeft, yBottom)
            for i in range(int(width)):
                n = int(wScale*i)
                myPolygon << QtCore.QPoint(xLeft+i, yTop+hScale*(_heightMax-_bins[n]))
            painter.drawPolygon(myPolygon)
        if not self.pos:
            if self.isovalue is not None:
                self.setIsoVal(self.isovalue, redraw=False)
            else:
                self.pos = [xRight, yBottom]
        x = self.pos[0]
        # add text with corresponding isovalue on top of the histogram rectangle
        font = QtGui.QFont()
        font.setPixelSize(7)
        painter.setFont(font)
        painter.drawText(x-5, yTop-15, 45, 15, QtCore.Qt.AlignLeft, self.tr("%.2f"%self.isovalue))
        # draw a vertical line inside the rectangle at the position of the mouse
        # (if the rectangle was clicked and/or the mouse was moved) 
        pen.setColor("#FF0000")
        pen.setWidth(2)
        painter.setPen(pen)
        painter.drawLine(x, yTop, x, yBottom)
        #self.valChanged.emit(self.isovalue)
        painter.end()

class PopupWireOptions(QtGui.QWidget):
    valChanged = QtCore.Signal(int)
    
    def __init__(self, point, parent = None, lineWidth=2):    
        super(PopupWireOptions, self).__init__(parent)
        self.immediate = True
        self.mainLayout = layout = QtGui.QHBoxLayout(self)
        label = QtGui.QLabel("line width")
        lineWBox = QtGui.QSpinBox()
        lineWBox.setRange(1, 10)
        lineWBox.setValue(lineWidth)
        lineWBox.valueChanged.connect(self.setLineWidth)
        layout.addWidget(label)
        layout.addWidget(lineWBox)
        layout.setContentsMargins(2, 0, 0, 0)
        self.setLayout(layout)
        self.adjustSize()
        self.setWindowFlags(QtCore.Qt.Popup)

        # by default, a widget will be placed from its top-left corner, so
        # we need to move it to the left based on the widgets width
        self.move(point)

    def setLineWidth(self, val):
        self.valChanged.emit(val)

class PopupOutlineOptions(PopupWireOptions):
    def __init__(self,  point, outline, parent = None):
        super(PopupOutlineOptions, self).__init__(point, parent, outline.lineWidth)
        colorsIcon = QtGui.QIcon(os.path.join(PMVICONPATH, 'colorChooser24.png'))
        colorsB = QtGui.QPushButton()
        colorsB.setIcon(colorsIcon)
        colorsB.clicked.connect(self.outlineColorDialog_cb)
        self.mainLayout.addWidget(colorsB)
        self.adjustSize()
        self.currentColor = None
        self.outline = outline

    def outlineColorDialog_cb(self):
        color = self.outline.color
        qcolor = QtGui.QColor()
        qcolor.setRgbF(color[0], color[1], color[2])
        qcolor.setAlphaF(color[3])
        self.colorDialog = colorDialog = QtGui.QColorDialog(qcolor)
        colorDialog.setOption(QtGui.QColorDialog.ShowAlphaChannel)
        self.currentColor = None
        colorDialog.currentColorChanged.connect(self.outlineColorChanged_cb)
        colorDialog.finished.connect(self.outlineColorDialogClosed_cb)
        
        colorDialog.colorSelected.connect(self.setOutlineColor_cb)
        colorDialog.open()

    def outlineColorChanged_cb(self, val):
        if self.currentColor is None:
            #save current front polygons color
            self.currentColor = self.outline.color[:]
        self.setOutlineColor_cb(val)
        
    def outlineColorDialogClosed_cb(self, val):
        if val == 0: # Cancel selected color
            # need to restore saved colors
            if self.currentColor is not None:
                self.outline.Set(color=self.currentColor)
                self.outline.geom.viewer.Redraw()
                self.currentColor = None
                
    def setOutlineColor_cb(self, val):
        # color dialog's callback
        rgb = val.getRgbF()
        alpha = val.alphaF()
        self.outline.Set(color=[rgb[0], rgb[1], rgb[2], alpha])
        self.outline.geom.viewer.Redraw()

class GeomAppearanceWidget(QtGui.QWidget):
    def __init__(self, parent=None, geom=None):
        super(GeomAppearanceWidget, self).__init__()
        self.mainLayout=QtGui.QVBoxLayout()
        self.mainLayout.setContentsMargins(0,0,0,0)
        # two toolbars with action widgets (icons) for setting geometry front and back
        # drawing mode and colors.
        self.ftoolBar = ftoolBar = QtGui.QToolBar()
        self.btoolBar = btoolBar = QtGui.QToolBar()
        self.frontCB = QtGui.QCheckBox("f")
        self.backCB = QtGui.QCheckBox("b")

        ## self.fopacityBox = QtGui.QDoubleSpinBox()
        ## self.fopacityBox.setRange(0.0, 1.0)
        ## self.fopacityBox.setSingleStep(0.1)
        ## self.fopacityBox.setFrame(False)
        ## self.fopacityBox.setToolTip("front polygons opacity")
        
        ## self.bopacityBox = QtGui.QDoubleSpinBox()
        ## self.bopacityBox.setRange(0.0, 1.0)
        ## self.bopacityBox.setSingleStep(0.1)
        ## self.bopacityBox.setFrame(False)
        ## self.bopacityBox.setToolTip("back polygons opacity")

        self.setGeomValues = True
        if geom is not None:
            self.setGeom(geom) # this method sets the "f" and "b" checkbuttons and the values for the opacity spinboxes
        else: self.geom = None

        self.frontCB.stateChanged.connect(self.cull_cb)
        self.frontCB.setToolTip('Check to show front facing polygons')
        self.backCB.stateChanged.connect(self.cull_cb)
        self.backCB.setToolTip('Check to show back facing polygons')
        
        pointsIcon = QtGui.QIcon(os.path.join(PMVICONPATH, 'point24.png'))
        wireIcon = QtGui.QIcon(os.path.join(PMVICONPATH, 'wire24.png'))
        shadedIcon = QtGui.QIcon(os.path.join(PMVICONPATH, 'smooth24.png'))
        outlineIcon = QtGui.QIcon(os.path.join(PMVICONPATH, 'outline24.png'))
        colorsIcon = QtGui.QIcon(os.path.join(PMVICONPATH, 'colorChooser24.png'))
        actions = []
        # front polygons action icons (set drawing mode, colors)
        act = self.frontPointsB = QtGui.QAction(pointsIcon, "draw points for front polygons", self)# "points")
        actions.append((act, 'front', 'points'))
        act = self.frontWireB = QtGui.QAction(wireIcon, "draw wires for front polygons",self)#"wire")
        actions.append((act, 'front', 'wire'))
        act = self.frontShadedB = QtGui.QAction(shadedIcon, "draw shaded polygons for front polygons", self)# "shaded")
        actions.append((act, 'front','gouraud'))        
        act = self.frontOutlinedB = QtGui.QAction(outlineIcon, "draw shaded polygons with outline for front polygons", self)
        actions.append((act, 'front', 'outlined'))
        act = self.frontColorB = QtGui.QAction(colorsIcon, "color front polygons", self)
        actions.append((act, 'front', 'colors'))

        # back polygons action icons
        act = self.backPointsB = QtGui.QAction(pointsIcon, "draw points for back polygons", self)# "points")
        actions.append((act, 'back', 'points'))
        act = self.backWireB = QtGui.QAction(wireIcon, "draw wires for back polygons", self)#"wire")
        actions.append((act, 'back', 'wire'))
        act = self.backShadedB = QtGui.QAction(shadedIcon, "draw shaded polygons for back polygons", self)# "shaded")
        actions.append((act, 'back', 'gouraud'))        
        act = self.backOutlinedB = QtGui.QAction(outlineIcon, "draw shaded polygons with outline for back polygons", self)
        actions.append((act, 'back', 'outlined'))
        act = self.backColorB = QtGui.QAction(colorsIcon, "color back polygons", self)
        actions.append((act, 'back', 'colors')) 
        for act, face, mode in actions:
            act.triggered.connect(CallbackFunction(self.setMode_cb, mode, face))

        #self.fopacityBox.valueChanged.connect(CallbackFunction(self.setOpacity_cb, "front"))
        #self.bopacityBox.valueChanged.connect(CallbackFunction(self.setOpacity_cb, "back"))
        # add widgets and actions to the toolbars
        # front polygons toolbar:
        ftoolBar.addWidget(self.frontCB)
        # add "spacer" widget to align all buttons in the toolbar to the right 
        spacerW1 = QtGui.QWidget(self)
        spacerW1.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Preferred)
        spacerW1.setVisible(True)
        ftoolBar.addWidget(spacerW1)
        
        ftoolBar.addAction(self.frontPointsB)
        ftoolBar.addAction(self.frontWireB)
        ftoolBar.addAction(self.frontShadedB)
        ftoolBar.addAction(self.frontOutlinedB)
        ftoolBar.addAction(self.frontColorB)
        #ftoolBar.addWidget(self.fopacityBox)
        # allow the right-click on the tool bar action widgets to set the line width of the polygons:
        # front wire
        w = ftoolBar.widgetForAction(self.frontWireB)
        w.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        w.customContextMenuRequested.connect(CallbackFunction(self.actionRightClick_cb, 'front', 'wire', w))
        #front outlined
        w = ftoolBar.widgetForAction(self.frontOutlinedB)
        w.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        w.customContextMenuRequested.connect(CallbackFunction(self.actionRightClick_cb, 'front', 'outlined',w))
        # back polygons toolbar:
        btoolBar.addWidget(self.backCB)
        spacerW2 = QtGui.QWidget(self)
        spacerW2.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Preferred)
        spacerW2.setVisible(True)
        btoolBar.addWidget(spacerW2)
        btoolBar.addAction(self.backPointsB)
        btoolBar.addAction(self.backWireB)
        
        btoolBar.addAction(self.backShadedB)
        btoolBar.addAction(self.backOutlinedB)
        
        btoolBar.addAction(self.backColorB)
        #btoolBar.addWidget(self.bopacityBox)
        # right-click on actions:
        # back wire
        w = btoolBar.widgetForAction(self.backWireB)
        w.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        w.customContextMenuRequested.connect(CallbackFunction(self.actionRightClick_cb, 'back', 'wire', w))
        #back outlined
        w = btoolBar.widgetForAction(self.backOutlinedB)
        w.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        w.customContextMenuRequested.connect(CallbackFunction(self.actionRightClick_cb, 'back', 'outlined', w))
        # make icon size smaller
        #btoolBar.setIconSize(QtCore.QSize(22,22))
        #ftoolBar.setIconSize(QtCore.QSize(22,22))
        #print "tool bar icon size:", btoolBar.iconSize()
        self.mainLayout.addWidget(self.ftoolBar)
        self.mainLayout.addWidget(self.btoolBar)
        self.setLayout(self.mainLayout)

    def cull_cb(self, status):
        # callback of the "f" and "b" check boxes in the toolbars
        if not self.setGeomValues: return
        cf = not self.frontCB.isChecked()
        cb = not self.backCB.isChecked()
        if cf and cb:
            self.culling = 'front_and_back'
        elif cf:
            self.culling = 'front'
        elif cb:
            self.culling = 'back'
        else:
            self.culling = 'none'
        self.geom.Set(culling=self.culling)
        #print self.geom.culling
        self.geom.viewer.deleteOpenglList()
        self.geom.viewer.Redraw()

    def actionRightClick_cb(self, face, mode, actwidget, point):
        # callback of the right-click on the 'wire' and 'outlined' action icons
        # in the front and back toolbars
        # point is QPoint instance with coordinates of the mouse pointer  
        # map that point as a global position
        global_point = actwidget.mapToGlobal(point)
        wpos = global_point + QtCore.QPoint(0, actwidget.height()/2)
        #print "handle_right_click", face, mode, point, global_point, wpos
        lw = self.geom.lineWidth
        w = None
        if mode == "outlined":
            if hasattr (self.geom, "outline"):
                w = PopupOutlineOptions(wpos, self.geom.outline, parent=self)
        else:
            w = PopupWireOptions(wpos, parent=self, lineWidth=lw)
        if w:
            w.valChanged.connect(CallbackFunction(self.setLineWidth, mode))
            w.show()

    def setLineWidth(self, mode, val):
        #print "setting line width", mode, val
        if self.geom:
            if mode == "wire":
                self.geom.Set(lineWidth=val, inheritLineWidth=False)
            elif mode == "outlined":
                if hasattr (self.geom, "outline"):
                    self.geom.outline.Set(lineWidth=val)
                    self.geom.viewer.Redraw()
                    
    def setMode_cb(self, mode, face):
        # left-click on the icons in the toolbar. Sets the drawing mode and colors
        # for front and back polygons
        if mode=='colors':
            color = self.getCurrentColor(face)
            if len(color) == 1:
                qcolor = QtGui.QColor()
                qcolor.setRgbF(color[0][0], color[0][1], color[0][2])
            else:
                qcolor = QtCore.Qt.green
            self.colorDialog = colorDialog = QtGui.QColorDialog(qcolor)
            colorDialog.setOption(QtGui.QColorDialog.ShowAlphaChannel)
            colorDialog.currentColorChanged.connect(CallbackFunction(self.colorChanged_cb, face))
            self.currentFcolor = None
            self.currentBcolor = None
            colorDialog.finished.connect(CallbackFunction(self.colorDialogClosed_cb, face))

            colorDialog.colorSelected.connect(CallbackFunction(self.setColor, face))
            colorDialog.open()
        else:
            self.renderMode = mode
            geom = self.geom
            
            if mode=='points':
                frontPolyMode='point'
                if face=='front':
                    geom.Set(inheritFrontPolyMode=False,
                             frontPolyMode='point', outline=False,
                             culling=self.culling)
                else:
                    geom.Set(inheritBackPolyMode=False,
                             backPolyMode='point', outline=False,
                             culling=self.culling)
            elif mode=='wire':
                if face=='front':
                    geom.Set(inheritFrontPolyMode=False, frontPolyMode='line',
                             outline=False, polyFace=face, culling=self.culling)
                else:
                    geom.Set(inheritBackPolyMode=False, backPolyMode='line',
                             outline=False, culling=self.culling)
            elif mode=='gouraud':
                if face=='front':
                    geom.Set(inheritFrontPolyMode=False, frontPolyMode='fill',
                             inheritShading=False, shading='smooth',
                             outline=False, culling=self.culling)
                else:
                    geom.Set(inheritBackPolyMode=False, backPolyMode='fill',
                             inheritShading=False, shading='smooth',
                             outline=False, culling=self.culling)
            elif mode=='outlined':
                if face=='front':
                    geom.Set(inheritFrontPolyMode=False, frontPolyMode='fill',
                             inheritShading=False, shading='smooth',
                             outline=True, culling=self.culling)
                else:
                    geom.Set(inheritBackPolyMode=False, backPolyMode='fill',
                             inheritShading=False, shading='smooth',
                             outline=True, culling=self.culling)
            
    def setOpacity_cb(self, face, val):
        """Callback of the opacity spin boxes"""
        if not self.setGeomValues : return
        if face == "front":
            self.geom.Set(opacity=float(val), inheritMaterial=False, transparent=True)
        else:
            self.geom.Set(opacity=float(val), inheritMaterial=False, transparent=True,
                     polyFace='back')
        
    def colorChanged_cb(self, face, val):
        # color dialog's callback (when a new color is selected)
        if face == "front":
            if self.currentFcolor is None:
                #save current front polygons color
                self.currentFcolor = self.getCurrentColor(face)
        else:            
            if self.currentBcolor is None:
                #save current back polygons color
                self.currentBcolor = self.getCurrentColor(face)
        self.setColor(face, val)

    def getCurrentColor(self, face):
        # return current rgba of the geometry (self.geom)
        geom = self.geom
        from DejaVu2.viewerConst import OVERALL, PER_VERTEX
        from opengltk.OpenGL import GL
        if face == "front":
            if geom.materials[GL.GL_FRONT].binding[1]==OVERALL:
                color = numpy.copy(geom.materials[GL.GL_FRONT].prop[1])
            elif geom.materials[GL.GL_FRONT].binding[1]==PER_VERTEX:
                color = numpy.copy(geom.materials[GL.GL_FRONT].prop[1])
        else:            
            if geom.materials[GL.GL_BACK].binding[1]==OVERALL:
                color = numpy.copy(geom.materials[GL.GL_BACK].prop[1])
            elif geom.materials[GL.GL_BACK].binding[1]==PER_VERTEX:
                color = numpy.copy(geom.materials[GL.GL_BACK].prop[1])
        return color

    def setColor(self, face, val):
        # color dialog's callback
        rgb = val.getRgbF()
        alpha = val.alphaF()
        self.geom.frontAndBack=False
        self.geom.Set(materials=[[rgb[0], rgb[1], rgb[2], alpha]], redo=1,
                      tagModified=False, polyFace=face,
                      transparent='implicit', inheritMaterial=0)
            
    def colorDialogClosed_cb(self, face, val):
        #called when color dialog is closed. val == 0 (Cancel)  val == 1 (OK)
        if val == 0: # Cancel selected color
            # need to restore saved colors
            if face == "front":
                if self.currentFcolor is not None:
                    self.geom.Set(materials=self.currentFcolor, tagModified=False,
                                  transparent='implicit', polyFace="front", inheritMaterial=0, redo=1)
                    self.currentFcolor = None
            else:
                if self.currentBcolor is not None:
                    self.geom.Set(materials=self.currentBcolor, tagModified=False,
                                  transparent='implicit', polyFace="back", inheritMaterial=0, redo=1)
                    self.currentBcolor = None

    def setGeom(self, geom):
        # set self.geom. Update "f" & "b" checkbuttons and opacity value boxes.
        self.geom = geom
        from DejaVu2.viewerConst import CULLINGS_REV
        self.culling = CULLINGS_REV[self.geom.culling]
        self.setGeomValues = False
        if self.culling == "front":
            self.frontCB.setChecked(False)
            self.backCB.setChecked(True)
        elif self.culling == "back":
            self.frontCB.setChecked(True)
            self.backCB.setChecked(False)
        elif self.culling == "front_and_back":
            self.frontCB.setChecked(False)
            self.backCB.setChecked(False)
        else: # 'none'
            self.frontCB.setChecked(True)
            self.backCB.setChecked(True)
        #fopac = self.getCurrentColor("front")[0][3]
        #self.fopacityBox.setValue(fopac)
        #bopac = self.getCurrentColor("back")[0][3]
        #self.bopacityBox.setValue(bopac)
        self.setGeomValues = True

class IsoGeomProp(QtGui.QWidget):
    def __init__(self, parent=None, data=None, geom=None, app=None):
        super(IsoGeomProp, self).__init__(parent)
        mainLayout = QtGui.QVBoxLayout()
        self.isoValGroup = QtGui.QGroupBox("iso value")
        self.isoValGroup.setCheckable(False)
        vbox = QtGui.QGridLayout()
        self.geom = geom
        isovalue = 0.0
        if hasattr (geom, "_isoValue"):
            isovalue = geom._isoValue
        self.app = app # TrgMapsGui instance
        # histogram widget
        self.histW = hist = Histogram(parent=self.isoValGroup, data=data, isovalue=isovalue)
        hist.valChanged.connect(self.setGeomIsoValue_cb)
        isolayout = QtGui.QVBoxLayout()
        isolayout.addWidget(hist)
        vbox.addLayout(isolayout,0,0)
        # add a label for min value of the histogram
        valsLayout = QtGui.QHBoxLayout()
        self.minLabel = minLabel = QtGui.QLabel("%.2f"%hist.minVal)
        valLabel = QtGui.QLabel("val:")
        # entry (edit) widget to set isovalue
        self.isovalEdit = isovalEdit = QtGui.QLineEdit()
        isovalEdit.returnPressed.connect(self.setGeomIsoValue_textcb)
        metrics = QtGui.QFontMetrics(isovalEdit.font())
        isovalEdit.setFixedWidth(metrics.width("000000"))
        isovalEdit.setText("%.2f" % isovalue)
        self.maxEdit = maxEdit = QtGui.QLineEdit()
        #self.maxEdit.setInputMask("9")
        maxEdit.setFixedWidth(metrics.width("000000"))
        maxEdit.setText("%.2f"%hist.maxVal)
        maxEdit.returnPressed.connect(self.setMaxIsoVal)
        maxLabel = QtGui.QLabel("max:")
        valsLayout.addWidget(minLabel)
        valsLayout.addWidget(valLabel)
        valsLayout.addWidget(isovalEdit)
        valsLayout.addStretch(1)
        valsLayout.addWidget(maxLabel, stretch=1, alignment=QtCore.Qt.AlignRight)
        valsLayout.addWidget(maxEdit, stretch=1, alignment=QtCore.Qt.AlignRight)
        vbox.addLayout(valsLayout, 1,0)
        vbox.setContentsMargins(0,0,0,0)
        self.isoValGroup.setLayout(vbox)
        mainLayout.addWidget(self.isoValGroup)

        self.appearanceW = GeomAppearanceWidget(self, geom=geom)
        mainLayout.addWidget(self.appearanceW)
        
        self.setLayout(mainLayout)
        #print "line edit width:", self.isovalEdit.width()

    def setGeomIsoValue_textcb(self):
        # callback of the "Enter" or "return" key press in the "val:" lineEdit widget
        if not self.geom: return
        val = self.isovalEdit.text()
        try:
            val = float(val)
        except:
            self.isovalEdit.setText("%.2f" % self.geom._isoValue)
            return

        self.histW.setIsoVal(val)
        #self.setGeomIsoValue(val)

    def setGeomIsoValue_cb(self, val):
        # this is called when the histogramm isovalue changes
        self.setGeomIsoValue(val)
        self.isovalEdit.setText("%.2f" % val)
        
    def setGeomIsoValue(self, val):
        # calls  the application setIsoValue method
        #print "setGeomIsoValue:", val
        if not self.app: return
        if not self.geom: return
        self.app.setIsoValue(self.geom, val)
        self.geom._isoValue = val

    def updateIsoData(self, data, minval=None, maxval=None, isovalue=None):
        self.histW.pos = None # this will prevent from drawing the current iso value vertical line
        self.histW.setData(data, minval, maxval, isovalue=isovalue)
        _minval = self.histW.minVal
        _maxval = self.histW.maxVal
        if _maxval is not None:
            self.maxEdit.setText("%.2f"%_maxval)
        if _minval is not None:
            self.minLabel.setText("%.2f"%_minval)
        if isovalue is not None:
            self.isovalEdit.setText("%.2f" % isovalue)

    def setGeom(self, geom, value=None):
        self.geom = geom
        if value is not None:
            self.histW.setIsoVal(value)
        self.appearanceW.setGeom(geom)

    def setMaxIsoVal(self):
        txt = self.maxEdit.text()
        val = float(txt)
        self.histW.setMaxValue(val)
        self.geom._maxIsoValue = val
        gval = self.geom._isoValue
        if gval < self.histW.minVal:
            gval = self.histW.minVal
            self.geom._isoValue = gval
        elif gval > self.histW.maxVal:
            gval = self.histW.maxVal
            self.geom._isoValue = gval
        self.histW.setIsoVal(gval)

class OrthoSliceGui(QtGui.QWidget):
    
    def __init__(self, parent=None, callback=None, geom=None, direction="X", nslices=None, valrange=None):
        super(OrthoSliceGui, self).__init__(parent)
        self.callBack = callback
        self.geom = geom
        self.direction = direction
        #mainLayout = QtGui.QHBoxLayout()
        mainLayout = QtGui.QVBoxLayout()
        mainLayout.setContentsMargins(0,0,0,0)
        slLayout = QtGui.QHBoxLayout()
        self.sliceLabel = QtGui.QLabel(direction)
        slLayout.addWidget(self.sliceLabel)
        self.sliceW = QtGui.QSlider()
        self.sliceW.setOrientation(QtCore.Qt.Horizontal)
        self.sliceW.setMinimum(0)
        if nslices is not None:
            self.sliceW.setMaximum(nslices-1)
        self.sliceW.valueChanged.connect(self.setSlice_cb)
        self.sliceW.sliderReleased.connect(self.sliderReleased_cb)
        slLayout.addWidget(self.sliceW)
        self.immediate = True
        self.colorpicker = HLSAColorMapEditor(chanelSelector='combo')
        self.initRgbaRamp()
        self.colorpicker.mapChange.connect(self.color_cb)
        colLayout = QtGui.QVBoxLayout()
        colLayout.setContentsMargins(0,0,0,0)
        colLayout.addWidget(self.colorpicker)
        minMaxlayOut = QtGui.QHBoxLayout()
        minlab = QtGui.QLabel("min:")
        self.minLE = QtGui.QLineEdit()
        self.minLE.setMinimumSize(60, 20)
        #self.minLE.setMinimumHeight(20)
        maxlab = QtGui.QLabel("max:")
        self.maxLE = QtGui.QLineEdit()
        #self.maxLE.setSizePolicy(QtGui.QSizePolicy.Ignored, QtGui.QSizePolicy.Preferred)
        self.maxLE.setMinimumSize(60, 20)
        #self.maxLE.setMinimumHeight(20)
        if valrange is not None:
            self.setMinMax(valrange[0], valrange[1])
        self.maxLE.returnPressed.connect(self.maxVal_cb)
        self.minLE.returnPressed.connect(self.minVal_cb)
        minMaxlayOut.addWidget(minlab)
        minMaxlayOut.addWidget(self.minLE)
        minMaxlayOut.addWidget(maxlab)
        minMaxlayOut.addWidget(self.maxLE)
        mainLayout.addLayout(colLayout)
        colLayout.addLayout(minMaxlayOut)
        mainLayout.addLayout(slLayout)
        self.setLayout(mainLayout)
        self.values = None
        self.valid = True

    def initRgbaRamp(self, ncolors=64, alpha=0.5):
        ramp = RGBARamp(ncolors)
        ramp[:,3]=alpha
        self.colorpicker.setRampFromRGBA(ramp)

    def setMaxSlice(self, val):
        self.immediate = False
        self.sliceW.setMaximum(val)
        self.immediate = True
        
    def setSlice(self, val):
        self.immediate = False # this is to prevent the callback
        self.sliceW.setValue(val)
        self.immediate = True

    def sliderReleased_cb(self):
        self.values = None

    def setSlice_cb(self, val):
        # callback of the slider
        if not self.immediate: return
        #print "SLIDER val:", val
        self.callBack(geom=self.geom, slicen=val, direction=self.direction)
        
    def getSliceVal(self):
        return self.sliceW.value()
        
    # Callbacks of the min and max lineEdit widgets (on "Return" press)
    def maxVal_cb(self):
        try:
            val = float(self.maxLE.text())
        except:
            self.setValid(self.maxLE, False)
            return
        minval = self.minLE.text()
        try:
            minval = float(minval)
        except:
            self.setValid(self.minLE, False)
            return
        if minval > val:
            self.setValid(self.maxLE, False)
            return
        self.setValid(self.maxLE)
        if self.valid and self.callBack:
            self.callBack(geom=self.geom, rangeList=[minval, val], direction=self.direction)
        #self.maxValChanged.emit(val)
        
    def minVal_cb(self):
        try:
            val = float(self.minLE.text())
        except:
            self.setValid(self.minLE, False)
            return
        maxval = self.maxLE.text()
        try:
            maxval = float(maxval)
        except:
            self.setValid(self.maxLE, False)
            return
        if val > maxval:
            self.setValid(self.minLE, False)
            return
        self.setValid(self.minLE)
        if self.valid and self.callBack:
            self.callBack(geom=self.geom, rangeList=[val, maxval], direction=self.direction)
        #self.minValChanged.emit(val)

    def setValid(self, lineEdit=None, valid=True):
        if not lineEdit: lineEdit=self.minLE
        #print "setValid", valid
        if valid:
            lineEdit.setStyleSheet('''QLineEdit{color:black}''')
        else:
            lineEdit.setStyleSheet('''QLineEdit{color:red}''')
        self.valid = valid

    def setMinMax(self, minval, maxval):
        # sets the text inside LineEdit widgets
        try:
            minval = float(minval)
        except: return
        try:
            maxval = float(maxval)
        except:  return
        self.maxLE.setText("%.2f"%maxval)
        self.minLE.setText("%.2f"%minval)

    def color_cb(self):
        #callback of the coloreditor widget
        if self.callBack:
            self.callBack(geom=self.geom, rgbaRamp=self.getColorRamp(),direction=self.direction)

    def setGeom(self, geom):
        self.geom = geom
        
    def updateWidgets(self):
        rgbaRamp = self.geom.colormap.ramp
        #print "setGeom for", geom.name
        #print colorRamp
        #print
        #self.colorpicker.setRampFromRGBA(rgbaRamp)
        self.immediate = False
        ramp = [QtGui.QColor.fromRgbF(r,g,b,a).getHsl() for r,g,b,a, in rgbaRamp]
        self.colorpicker.setRampFromHSLA(ramp)
        maxi = self.geom.colormap.maxi
        mini = self.geom.colormap.mini
        setMinMax = True
        if mini is not None and maxi is not None:
            try:
                assert (mini <= maxi)
            except:
               setMinMax = False
            if setMinMax:
                self.setValid(self.minLE, True)
                self.minLE.setText("%.2f"%mini)
                self.setValid(self.maxLE, True)
                self.maxLE.setText("%.2f"%maxi)
        if hasattr (self.geom, "_orthoSliceNum"):
            sliderVal = self.geom._orthoSliceNum[self.direction]
        else:
            sliderVal = 0
        self.sliceW.setValue(sliderVal)
        self.immediate = True
                
    def setDirection(self, direction):
        self.direction = direction
        self.sliceLabel.setText(direction)

    def getValues(self):
        # returns values of the line edit widgets and coloreditor ramp
        if not self.valid: return []
        minval = self.minLE.text()
        if not len(minval): return []
        maxval = self.maxLE.text()
        if not len(maxval): return []
        return [float(minval), float(maxval), self.getColorRamp()]

    def getColorRamp(self):
        ramp = self.colorpicker.getRGBRamp()
        #print "COLOR RAMP", ramp
        return ramp #self.colorpicker.getRGBRamp()
        

if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    mapw  =  TrgMapsGui()
    mapw.addTrg("ligPocket.trg")
    mapw.show()
    sys.exit(app.exec_())

#else:
#    vi = pmv.gui().viewer
#    mapw  =  TrgMapsGui(app=pmv)
#    mapw.viewer = vi
#    mapw.addTrg("ligPocket.trg")
#    mapw.show()
    
    

