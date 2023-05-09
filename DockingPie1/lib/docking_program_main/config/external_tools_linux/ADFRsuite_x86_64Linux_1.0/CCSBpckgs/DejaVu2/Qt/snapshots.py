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

from PySide import QtGui, QtCore

from DejaVu2.states import setRendering, getRendering, getOrientation
from DejaVu2.moveGeom import MoveGeom
from mglutil.util.callback import CallbackFunction

from os import path
import os
import time
import weakref



def _getRendering(viewer, checkAnimatable=False):
    state = getRendering(viewer, checkAnimatable)
    if state.has_key("camera"):
        state.pop("camera")
    if viewer.currentCamera.ssao: # ambient occlusion is on
        state["camera"] = {"near":2.0}
    else:
        state["camera"]= {'near': viewer.currentCamera.near}
    return state

iconw = 50
iconh = 50

### The following ImageQt class is from PIL library. (PIL is using PyQt4)
##
# (Internal) Turns an RGB color into a Qt compatible color integer.
# We cannot use qRgb directly for this, since it returns a long
# integer, but setColorTable can only deal with integers.

def rgb(r, g, b):
    # use qRgb to pack the colors, and then turn the resulting long
    # into a negative integer with the same bitpattern.
    return (QtGui.qRgb(r, g, b) & 0xffffff) - 0x1000000

class ImageQt:
    
    ##
    # A reference to a QImage object.

    image = None

    def __init__(self, im):

        data = None
        colortable = None

        if im.mode == "1":
            format = QtGui.QImage.Format_Mono
        elif im.mode == "L":
            format = QtGui.QImage.Format_Indexed8
            colortable = []
            for i in range(256):
                colortable.append(rgb(i, i, i))
        elif im.mode == "P":
            format = QtGui.QImage.Format_Indexed8
            colortable = []
            palette = im.getpalette()
            for i in range(0, len(palette), 3):
                colortable.append(rgb(*palette[i:i+3]))
        elif im.mode == "RGB":
            data = im.tostring("raw", "BGRX")
            format = QtGui.QImage.Format_RGB32
        elif im.mode == "RGBA":
            try:
                data = im.tostring("raw", "BGRA")
            except SystemError:
                # workaround for earlier versions
                r, g, b, a = im.split()
                from PIL import Image
                im = Image.merge("RGBA", (b, g, r, a))
            format = QtGui.QImage.Format_ARGB32
        else:
            raise ValueError("unsupported image mode %r" % im.mode)

        # must keep a reference, or Qt will crash!
        self.data = data or im.tostring()

        self.image = QtGui.QImage(self.data, im.size[0], im.size[1], format)

        if colortable:
            self.image.setColorTable(colortable)

    def __getattr__(self, attr):
        return getattr(self.image, attr)


class Snapshot:

    def __init__(self, viewer, name, image=None):
        self.viewer = viewer
        self.object = obj = viewer.rootObject
        self.moveGeom = mg = MoveGeom(obj)
        self.transformation = mg.getTransformation(obj)
        self.cameraAttr = ca = mg.getCameraAttr()
        self.rendering = _getRendering(viewer, checkAnimatable=True)
        self.image=image
        self.name = name

    def interpolate(self, transf=None, camAttr=None, time=0.5):
        if transf is None:
            transf = self.moveGeom.getTransformation(self.object)
        if camAttr is None:
            camAttr = self.moveGeom.getCameraAttr()
        #setRendering(self.viewer, self.rendering)
        self.moveGeom.interpolate(time, self.object, transf=[{self.object:transf},
                                                             {self.object:self.transformation}],
                           camattr=[camAttr, self.cameraAttr])
        setRendering(self.viewer, self.rendering)
        
        
        
class SnapshotsGUI(QtGui.QWidget):
    
    def __init__(self, app, parent=None):
        super(SnapshotsGUI, self).__init__(parent)
        self.mainLayout =  layout = QtGui.QVBoxLayout()
        self.snapshotsLW = QtGui.QListWidget()
        layout.addWidget(self.snapshotsLW)
        w = 210
        h = 100
        self.snapshotsLW.setMinimumSize(QtCore.QSize(w, h))

        #self.label = QtGui.QLabel(self)
        #self.label.setMinimumSize(iconw, iconh)
        #layout.addWidget(self.label)
        
        group = QtGui.QGroupBox("tools")
        groupLayout = QtGui.QGridLayout()
        addButton = QtGui.QPushButton('Add Snapshot')
        addButton.clicked.connect(self.addSnapShot)
        groupLayout.addWidget(addButton, 0, 0)
        removeButton = QtGui.QPushButton('Remove Snapshot')
        removeButton.clicked.connect(self.removeSnapShot)
        groupLayout.addWidget(removeButton, 0, 1)

        #up-down buttons
        ## self.upButton = QtGui.QPushButton()
        ## self.upButton.setIcon(self.style().standardIcon(QtGui.QStyle.SP_ArrowUp))
        ## self.downButton = QtGui.QPushButton()
        ## self.downButton.setIcon(self.style().standardIcon(QtGui.QStyle.SP_ArrowDown))
        ## groupLayout.addWidget(self.downButton, 0,2)        
        ## groupLayout.addWidget(self.upButton, 0,3)

        group.setLayout(groupLayout)
        layout.addWidget(group)

        self.setLayout(layout)
        self.count = 0
        self.snapshotsLW.itemClicked.connect(self.run_cb)
        
        self.app = app
        self.viewerName = app.pmv.name
        self.snapshots = {}
        self.currentRow = None
        self.interpTime = 0.5

    def run_cb(self, item=None):
        print "run_cb item:", item, type(item)
        #print "item 2:", self.listWidget.currentItem()
        snapshotName = str(item.text())
        snapshot = self.snapshots[snapshotName]
        snapshot.interpolate(time=self.interpTime)

    def setInterpolationTime(self, time):
        self.interpTime = time

    def addSnapShot(self):
        snapshot = self.recordOrient()
        im = self.getIconImage()
        qtIm = ImageQt(im)
        txt = "Snapshot %s" % str(self.count)
        #self.qtImages[txt] = qtIm
        snapshot.image = qtIm
        pix = QtGui.QPixmap.fromImage(qtIm.image)
        item = QtGui.QIcon(pix)
        
        self.snapshotsLW.insertItem(self.count, QtGui.QListWidgetItem(item, txt))
        self.snapshotsLW.setIconSize(QtCore.QSize(iconw, iconh))
        #self.label.setPixmap(pix)
        self.count+=1

    def getIconImage(self):
        from PIL import Image, ImageChops
        vi = self.app.viewer
        vi.OneRedraw()
        im = vi.currentCamera.GrabFrontBuffer()
        imc = None
        def autocrop(im, bgcolor):
            if im.mode != "RGB":
                im = im.convert("RGB")
            #print "image size:",  im.size, bgcolor
            bg = Image.new("RGB", im.size, bgcolor)
            diff = ImageChops.difference(im, bg)
            bbox = diff.getbbox()
            print "image bbox:", bbox
            if bbox:
                return im.crop(bbox)
            return None # no contents

        #imc = autocrop(im, tuple([int(round(x*255)) for x in c.backgroundColor[:3]]))
        if imc is None:
            print "imc None"
            imc = im
        print "image size 2:", imc.size 
        return imc.resize((iconw, iconh), Image.ANTIALIAS)


    def recordOrient(self, event=None):
        """
        build default orientation transition (left clicks)
        """
        obj = self.app.viewer.rootObject
        snName = self.checkName("Snapshot %d"%self.count)
        snapshot = Snapshot(self.app.viewer, snName) 
        self.snapshots[snName] =  snapshot
        return snapshot
        
    def checkName(self, name):
        """check if the name exists in the self.snapshots or in the sequence player.
        If exists - create unique name"""
        allnames = self.snapshots.keys()
        if name in allnames:
            i = 1
            while(name in allnames):
                name = name+"_%d"%i
                i = i+1
        return name

    def removeSnapShot(self):
        item = self.snapshotsLW.currentItem()
        if not item: return
        txt = str(item.text())
        print "removing current item",  txt
        row = self.snapshotsLW.currentRow()
        item = self.snapshotsLW.takeItem(row)
        if self.snapshots.has_key(txt):
            self.snapshots.pop(txt)


class SnapshotDialog(QtGui.QDialog):
    def __init__(self, app, parent=None):
        super(SnapshotDialog, self).__init__(parent)
        vertLayout = QtGui.QVBoxLayout()
        w = SnapshotsGUI(app)
        vertLayout.addWidget(w)
        # OK and Cancel buttons
        self.buttons = QtGui.QDialogButtonBox(
            #QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel,
            QtGui.QDialogButtonBox.Ok,
            QtCore.Qt.Horizontal, self)
        vertLayout.addWidget(self.buttons)
        self.setLayout(vertLayout)
        self.buttons.accepted.connect(self.accept)
        #self.buttons.rejected.connect(self.reject)
        
    @staticmethod
    def addSnapshots(app, parent=None):
        #import pdb; pdb.set_trace()      
        dialog = SnapshotDialog(app)
        dialog.setWindowModality(QtCore.Qt.WindowModal) # WindowModal
        result = dialog.exec_()
        #dialog.show()
        return (result == QtGui.QDialog.Accepted)

        


if __name__ == '__main__':
  
    import sys
  
    app = QtGui.QApplication(sys.argv)
    ok = SnapshotDialog.addSnapshots(None)
    ## window = Snapshots()
    ## window.show()
    ## print "size2:", window.snapshotsLW.contentsSize()
    ## sys.exit(app.exec_())
        
## app = QtGui.QApplication([])

## window = QtGui.QWidget()
## layout = QtGui.QVBoxLayout(window)
## list = QtGui.QListWidget()
## items = [str(i) for i in range(100)]
## list.addItems(items)
## #list.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
## #list.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
## #list.setFixedSize(list.sizeHintForColumn(0) + 2 * list.frameWidth(), list.sizeHintForRow(0) * list.count() + 2 * list.frameWidth())
## layout.addWidget(list)
## window.show()
## app.exec_()
