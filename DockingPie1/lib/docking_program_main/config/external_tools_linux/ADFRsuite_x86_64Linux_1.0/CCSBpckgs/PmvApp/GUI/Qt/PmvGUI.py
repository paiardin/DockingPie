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
# Date: 2014 Authors: Michel Sanner
#
#    sanner@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Michel Sanner and TSRI
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/PmvApp/GUI/Qt/PmvGUI.py,v 1.59.4.6 2017/10/26 23:12:16 annao Exp $
#
# $Id: PmvGUI.py,v 1.59.4.6 2017/10/26 23:12:16 annao Exp $
#
import os, sys, weakref, numpy, platform
from time import time

from PySide import QtGui, QtCore

from MolKit2.molecule import Atom, Molecule, MoleculeSet, Residue, Chain
from MolKit2.selection import SelectionSet, selectionTypes

from mglutil.util.callback import CallbackFunction
from mglutil.util.uniq import uniq
from mglutil.util.packageFilePath import findFilePath
PMVICONPATH = findFilePath('Icons', 'PmvApp.GUI')

from DejaVu2.Spheres import Spheres
from DejaVu2.Cylinders import Cylinders
from DejaVu2.IndexedPolygons import IndexedPolygons

#from PmvApp.Pmv import Selection, numOfSelectedVerticesToSelectTriangle
from PmvApp.Pmv import PmvSelection, numOfSelectedVerticesToSelectTriangle
from MolKit2.selection import Selection, SelectionSet
from PmvApp.group import Group

use_ipython_shell=False

try:
    from IPython.qt.console.rich_ipython_widget import RichIPythonWidget
    from IPython.qt.inprocess import QtInProcessKernelManager
    from IPython.lib import guisupport
    use_ipython_shell=True
except ImportError:
    print "Warning: failed to import IPython.qt module. Cannot use IPython shell widget,\nusing PyShell instead." 

def busy_cursor(function):
    def new_function(self):
        QtGui.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        try:
            function(self)
        except Exception as e:
            raise e
            print("Error {}".format(e.args[0]))
        finally:
            QtGui.QApplication.restoreOverrideCursor()
    return new_function

class GetNewName(QtGui.QInputDialog):
    # dialog to get one or more new names for dashboard entries
    
    # class to get unique class names
    def __init__(self, objItems, pmvGUI, text="New name for"):
	self.pmvGUI = pmvGUI
        self.objItems = objItems
	QtGui.QInputDialog.__init__(self, pmvGUI)
	self.setInputMode(QtGui.QInputDialog.TextInput)
        if len(objItems)==1:
            obj0, item0 = objItems[0]
            if item0 is None:
                name = obj0.name
            else:
                name = item0.text(0)
            self.setLabelText('%s %s'%(text, name))
	else:
            itemNames = ''
            for obj, item in objItems:
                if item:
                    itemNames += item.text(0)
                else:
                    itemNames += obj.name
                itemNames += ', '
            itemNames = itemNames[:-2]
            if len(itemNames)>63:
                itemNames = itemNames[:30] + "..." + itemNames[-30:]
            self.setLabelText('%s %s'%(text, itemNames))
            
    def done(self, result):
        # validate names
	if result==QtGui.QDialog.Accepted:
	    names = self.textValue().encode('ascii', 'ignore').split(',')
            if len(names)==1:
                names = names*len(self.objItems)
            elif len(names)!=len(self.objItems):
		self.pmvGUI.statusBar().showMessage(self.tr(
			"for %d items to be renames provide 1 or %d comma separated names.try again..."))
                return
            else:
                # make sure the name is unique
                if len(names)!=len(uniq(names)):
                    self.pmvGUI.statusBar().showMessage(self.tr(
			"Selection names have to be unique.try again..."))
                    return

            # at this point we should have good names
            object0, item0 = self.objItems[0]
            if isinstance(object0, PmvSelection):
                # makes sure the names do not exist yet in selections
                for name in names:
                    # FIXME: we could check against names in list of children of invisible root ?
                    if name in self.pmvGUI.pmv.namedSelections.keys() or name=='Current Selection' or \
                           name in self.pmvGUI.pmv.Mols.name:
                        self.pmvGUI.statusBar.showMessage(self.tr(
                            "Invalid name. Selection names must be unique (%s) ..try again..."%name))
                QtGui.QInputDialog.done(self, result)
	    else:
		QtGui.QInputDialog.done(self, result)
	else:
	    QtGui.QInputDialog.done(self, result)

from PySide import QtCore
import urllib2                   
from xml.dom import minidom

class FetchGUI(QtGui.QDialog):

    def __init__(self, parent=None):
        super(FetchGUI, self).__init__(parent)
        hLayout = QtGui.QHBoxLayout()
        self.pdbidEditW = QtGui.QLineEdit("")
        self.searchButton = QtGui.QPushButton()
        self.searchButton.setIcon(QtGui.QIcon(os.path.join(PMVICONPATH, "search-icon.png")))
        self.searchButton.clicked.connect(self.showTableDialog)
        self.searchButton.setToolTip("Open a keyword search dialog")
        hLayout.addWidget(self.pdbidEditW)
        hLayout.addWidget(self.searchButton)
        layout = QtGui.QFormLayout()
        layout.addRow("Pdb Id:", hLayout)
        self.formats=["mmtf", "pdb"]
        self.formatBox = QtGui.QComboBox()
        self.formatBox.addItems(self.formats)
        layout.addRow("Format:", self.formatBox)

        self.buttons = QtGui.QDialogButtonBox(QtCore.Qt.Horizontal, self)
        fetchButton = QtGui.QPushButton(self.tr("&Fetch"))
        self.buttons.addButton(fetchButton, QtGui.QDialogButtonBox.AcceptRole)
        self.buttons.addButton(QtGui.QDialogButtonBox.Close)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        fetchButton.setDefault(False)
        layout.addRow("", self.buttons)
        self.setLayout(layout)
        self.pdbidEditW.setFocus()
        
    def accept(self):
        if not len(self.pdbidEditW.text()):
            flags = QtGui.QMessageBox.StandardButton.Ok
            msg = "Pdb Id is not specified"
            response = QtGui.QMessageBox.warning(self, "Warning!", msg, flags)
        else:
            QtGui.QDialog.accept(self)
            
    def showTableDialog(self):
        # opens a second dialog with an entry for keyword search and a 
        # QTableWidget that will contain the results of the search.
        dialog = TableWidgetDialog("Search PDB", parent=self)
        result = dialog.exec_()
        selected = dialog.table.selectedIndexes()
        #print "TableWidgetDialog:", result
        if result == QtGui.QDialog.Accepted:
            # update the entry with selected PdbIds 
            pdbId = ""
            for item in selected:
                column = item.column()
                if column == 0:
                    pdbId += str(item.data()) + " "
            self.pdbidEditW.setText(pdbId)
            
    @staticmethod
    def getName(parent = None):
        dialog = FetchGUI(parent=parent)
        dialog.setWindowTitle(QtGui.QApplication.translate("Dialog", "Fetch molecule" , None, QtGui.QApplication.UnicodeUTF8))
        #dialog.setWindowModality(QtCore.Qt.WindowModal)
        result = dialog.exec_()
        #print "FetchGUI result:", result
        molname = dialog.pdbidEditW.text()
        molformat = dialog.formatBox.currentText()
        return ([str(molname), str(molformat)], result == QtGui.QDialog.Accepted)

import re

class TableWidgetDialog(QtGui.QDialog):
    searchResult = QtCore.Signal(str)
    searchEnd = QtCore.Signal()

    uniprotACIDpattern = re.compile("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")

    url = 'http://www.rcsb.org/pdb/rest/search'
    txtSearch = """<orgPdbQuery>
    
<queryType>org.pdb.query.simple.AdvancedKeywordQuery</queryType>

<description>Text Search for: %s</description>

<keywords>%s</keywords>

</orgPdbQuery>"""

    uniprotACSearch = """<orgPdbQuery>

<queryType>org.pdb.query.simple.UpAccessionIdQuery</queryType>

<description>Simple query for a list of Uniprot Accession IDs: %s</description>

<accessionIdList>%s</accessionIdList>

</orgPdbQuery>"""

    def __init__(self, title, parent=None):
        super(TableWidgetDialog, self).__init__(parent)
        self.dataList = []
        self.header = header = ['ID',  '#Res', 'Method', 'Reso', 'Date', 'Title']
        ncolumns = len(header)
        
        self.setMinimumSize(800, 400)
        self.setWindowTitle(title)

        hLayout = QtGui.QHBoxLayout()
        self.keywordEditW =  QtGui.QLineEdit()
        #self.keywordEditW.setMaximumWidth(150)
        searchLabel = QtGui.QLabel("Keyword:")
        self.searchButton = QtGui.QPushButton("Search")
        self.searchButton.clicked.connect(self.fillTableDialog)
        self.keywordEditW.returnPressed.connect(self.fillTableDialog)
        hLayout.addWidget(searchLabel)
        hLayout.addWidget(self.keywordEditW)
        hLayout.addWidget(self.searchButton)
        
        self.table = table = QtGui.QTableWidget(1, ncolumns, self)
        table.setHorizontalHeaderLabels(header)
        table.horizontalHeader().setResizeMode(ncolumns-1, QtGui.QHeaderView.Stretch)
        #table.resizeColumnsToContents()
        table.setColumnWidth(0, 45)
        table.setColumnWidth(1, 45)
        table.setColumnWidth(2, 80)
        table.setColumnWidth(3, 50)
        table.setColumnWidth(4, 80)
        table.cellClicked.connect(self.cellClicked_cb)
        lfont = QtGui.QFont()
        lfont.setPointSize(8)
        table.setFont(lfont)
        layout = QtGui.QVBoxLayout()
        layout.addLayout(hLayout)
        layout.addWidget(table)

        # enable sorting
        table.horizontalHeader().setSortIndicator(0, QtCore.Qt.AscendingOrder)
        # setSortingEnabled() calls the sort() of the table_model in the default Discending order.
        #setSortIndicator() is used to change the initial sort order to Ascending.
        #table.setSortingEnabled(True)
        
        #add progress bar , "close", "accept" buttons
        buttonsLayout = QtGui.QHBoxLayout()
        self.progress = QtGui.QProgressBar(self)
        self.numIdFoundLabel = QtGui.QLabel("")
        self.buttons = QtGui.QDialogButtonBox(QtCore.Qt.Horizontal, self)
        self.acceptButton = acceptButton = QtGui.QPushButton(self.tr("&Use selected PDB id"))
        self.buttons.addButton(acceptButton, QtGui.QDialogButtonBox.AcceptRole)
        self.buttons.addButton(QtGui.QDialogButtonBox.Close)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        acceptButton.setDefault(False)
        acceptButton.setEnabled(False)
        buttonsLayout.addWidget(self.progress)
        buttonsLayout.addWidget(self.numIdFoundLabel)
        buttonsLayout.addWidget(self.buttons)

        #self.buttons.button(self.buttons.Close).setDefault(True)
        #layout.addWidget(self.buttons)
        layout.addLayout(buttonsLayout)
        self.setLayout(layout)
        self.searchResult.connect(self.fillTable)

    def keyPressEvent(self, evt):
        if evt.key() == QtCore.Qt.Key_Enter or evt.key() == QtCore.Qt.Key_Return:
            return
        QtGui.QDialog.keyPressEvent(self, evt)
        
    def fillTableDialog(self):
        # display search results in the table
        keyword = str(self.keywordEditW.text())
        self.table.clearContents()
        if self.table.rowCount()> 1:
            self.table.setRowCount(1)
        if len(keyword):
            import thread
            thread.start_new_thread(self.keywordSearch, (keyword,) )
        
    def keywordSearch(self, keyword):
        match = self.uniprotACIDpattern.match(keyword.strip())
        if match is None:
            req = urllib2.Request(self.url, data=self.txtSearch%(keyword, keyword))
        else:
            req = urllib2.Request(self.url, data=self.uniprotACSearch%(keyword, keyword))
            
        f = urllib2.urlopen(req)
        self._result = f.read()
        #print 'time search', time()-t0
        pdbIds = self._result.replace('\n', ',')[:-1] 
        url = 'http://www.rcsb.org/pdb/rest/customReport.xml?pdbids=%s&customReportColumns=structureId,residueCount,experimentalTechnique,resolution,releaseDate,structureTitle&service=wsfile&format=csv '%pdbIds
        try:
            response = urllib2.urlopen(url)
            result = response.read()
            #self.searchResult.emit( result.replace('"','') ) # this calls self.fillTable
            self.searchResult.emit( result) # this calls self.fillTable
        except urllib2.HTTPError:
            self.searchResult.emit( 'No match for: "%s"'%keyword )
    
    def fillTable(self, result):
        if result.startswith('No match for:'):
            return
        lines = result.split('\n')[1:-1]
        nlines = len(lines)
        self.numIdFoundLabel.setText("Found %d entities"%nlines)
        for n, line in enumerate(lines):
            ncolumns = len(self.header)
            for j, txt in enumerate(line.split('",')):
                newitem = QtGui.QTableWidgetItem(txt.replace('"',''))
                self.table.setItem(self.table.rowCount()-1, j, newitem)
                if j == ncolumns-1:
                    font = QtGui.QFont()
                    font.setUnderline(True)
                    #font.setPointSize(8)
                    newitem.setFont(font)
                    newitem.setForeground(QtGui.QBrush(QtGui.QColor("blue")))
                    if n < nlines-1:
                        self.table.insertRow(self.table.rowCount())
            self.progress.setValue(int(100*n/nlines))
        self.progress.setValue(0)
        if not self.table.isSortingEnabled():
            self.table.setSortingEnabled(True)
              
    def cellClicked_cb(self, row, column):
        #enable "Use selected PDB IDs " button if any pdbid is selected, disable otherwise 
        selecdedPdb = 0
        selected = self.table.selectedIndexes()
        enabled = self.acceptButton.isEnabled()
        for item in selected:
            if item.column() == 0:
                if item.data() is not None:
                    selecdedPdb+=1
        if selecdedPdb > 0 and not enabled:
            self.acceptButton.setEnabled(True)
        elif selecdedPdb == 0 and enabled:
            self.acceptButton.setEnabled(False)
        # open link in the browser if a cell in the last column is clicked    
        if column == len(self.header)-1:
            import webbrowser
            pdbid = str(self.table.item(row, 0).text())
            url = "http://www.rcsb.org/pdb/explore/explore.do?structureId=%s" %pdbid
            webbrowser.open(url, new=0, autoraise=True)

class GetCLoneSelectionName(GetNewName):
    # dialog to get one or more new names for dashboard entries
    
    def done(self, result):
        # validate names
	if result==QtGui.QDialog.Accepted:
	    names = self.textValue().encode('ascii', 'ignore').split(',')
            if len(names) > 1:
		self.pmvGUI.statusBar.showMessage(self.tr(
			"Only one name is allowed.Try again..."))
                return

            # at this point we should have good names
            object0, item0 = self.objItems[0]
            if isinstance(object0, PmvSelection):
                ## # makes sure the names do not exist yet in selections
                ## for name in names:
                ##     # FIXME: we could check against names in list of children of invisible root ?
                ##     if name in [mol.name for mol in self.pmvGUI.pmv.Mols]:
                ##         self.pmvGUI.statusBar.showMessage(self.tr(
                ##             "Invalid name. Molecule name must be unique (%s) ..try again..."%name))
                ##         return
                QtGui.QInputDialog.done(self, result)
	    else:
		QtGui.QInputDialog.done(self, result)
	else:
	    QtGui.QInputDialog.done(self, result)

def contentFont():
        font = QtGui.QFont()
        font.setStyleStrategy(QtGui.QFont.PreferAntialias)
 
        #if sys.platform == 'darwin':
        #    font.setPixelSize(9)
        #    font.setFamily('Arial')
        #else:
        #    font.setPixelSize(11)
        #    font.setFamily('Arial')
        font.setPixelSize(11)
        font.setFamily('Arial')
 
        return font

class GridGroup:

    def __init__(self, name, children=None):
        self.name = name
        if children is None:
            children= []
        self.children = children


from DejaVu2.Qt.Viewer import Viewer


class Worker(QtCore.QThread):
    done = QtCore.Signal()

    def __init__(self, func, args, kw):
        QtCore.QThread.__init__(self)
        self.func = func
        self.args = args
        self.kw = kw
        
    #A QThread is run by calling it's start() function, which calls this run()
    #function in it's own "thread". 
    def run(self):
        self.func(*self.args, **self.kw)
        self.done.emit()
        
class SaveSessionDialog(QtGui.QDialog):

    def __init__(self, func, args=[], kw={}, estimatedSize=None, parent=None):
        super(SaveSessionDialog, self).__init__(parent)
        self.filename = args[0]
        self.estimatedSize = estimatedSize
        self.func = func
        self.args = args
        self.kw = kw
        self.layout = QtGui.QVBoxLayout()
        self.progressBar = QtGui.QProgressBar()
        self.progressBar.minimum = 0
        self.progressBar.maximum = 0
        self.layout.addWidget(self.progressBar)
        self.fileSizeLabel = QtGui.QLabel(self)
        self.fileSizeLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Sunken)
        self.fileSizeLabel.setText("file size 0 Mb")
        self.layout.addWidget(self.fileSizeLabel)

        #self.pushButton = QtGui.QPushButton()
        #self.pushButton.setText("Cancel")
        #self.pushButton.released.connect(self.cancel)
        #self.layout.addWidget(self.pushButton)

        self.setLayout(self.layout)
        self.buttons = QtGui.QDialogButtonBox(QtCore.Qt.Horizontal, self)

        self.worker = Worker(func, args, kw)
        self.worker.done.connect(self.accept)
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.setProgress)
        self.timer.start(0.1)
        self.worker.start()

    #def cancel(self):
        #print 'rejecting'
        #val = self.reject()
        #print 'got val\n stoping timer'
        #self.timer.stop()
        #print 'done \n terminating thread'
        #self.worker.terminate()
        #self.worker.wait()
        #print 'done'

        #return self.reject()
        #return QtGui.QDialog.Rejected
    
    def setProgress(self):
        #print 'SET PROGRSS'
        filesize = os.path.getsize(self.filename)
        sizeMB = filesize/1048576.
        if self.estimatedSize:
            self.progressBar.maximum = 100
            self.progressBar.setValue(100*sizeMB/self.estimatedSize)
        else:
            self.progressBar.setValue(0)
        self.fileSizeLabel.setText("file size %.2f MB out of estimated %.2f"%(sizeMB, self.estimatedSize))


class PmvViewer(Viewer, QtGui.QWidget):

    levelLabels = [ 'M:', 'C:', 'R:', 'A:']

    def __init__(self, pmv, master=None, nogui=0, 
                 guiMaster=None, showViewerGUI=True,
                 autoRedraw=True, cnf={}, **kw):

        print 'PMVVIEWER MASTER', id(master)
        QtGui.QWidget.__init__(self, master)
        Viewer.__init__(self, nogui=nogui, guiMaster=guiMaster,
                        showViewerGUI=showViewerGUI,
                        autoRedraw=autoRedraw, cnf=cnf, **kw)
        self.pmv = pmv
        self.statusBar = None      
        # self.processPicking will turn the pick into atoms and call self.setDragSelectCommand
        self.setDragSelectCommand(self.pmv.select)

        # this allows to specify what should happen when we pick on nothing, Set to None
        # if no action is wanted
        self.setEmptyPickCommand(None)

        # create a crossSet geom for showing selections
        from DejaVu2.Points import CrossSet
        self.selectionCrosses = g = CrossSet(
            'selectionCrosses', shape=(0,3), materials=((1.0, 1.0, 0.),), lineWidth=2,
            inheritMaterial=0, protected=True, disableStencil=True, visible=False,
            transparent=True, animatable=False)
        g.pickable = 0
        self.AddObject(g)
        self.selectionCrosses = g

        self.registerListerners()
 
        # create icons
        self.icons = {}
        self.icons['atom'] = QtGui.QIcon(os.path.join(PMVICONPATH, "atom.png"))
        self.icons['residue'] = QtGui.QIcon(os.path.join(PMVICONPATH, "sidechain.png"))
        self.icons['chain'] = QtGui.QIcon(os.path.join(PMVICONPATH, "chain.png"))
        self.icons['molecule'] = QtGui.QIcon(os.path.join(PMVICONPATH, "molecule.png"))
        self.icons['selection'] = QtGui.QIcon(os.path.join(PMVICONPATH, "selection.png"))
        self.icons['group'] = QtGui.QIcon(os.path.join(PMVICONPATH, "group.png"))

        from PySide import QtOpenGL

        fmt = QtOpenGL.QGLFormat()
        fmt.setDoubleBuffer(True)
        fmt.setOverlay(True)
        fmt.setDepth(True)
        fmt.setStencil(True)
        fmt.setAccum(True)
        fmt.setStereo(False)
        fmt.setSampleBuffers(True)
        self.hiddenGLWidget = QtOpenGL.QGLWidget(fmt, parent=None)

        ## overriding AddCamera causes segfault :(
        self.AddCamera(master=master, verbose=True)
        self.cameras[-1]._mouseReleaseNoMotionActions[
           int(QtCore.Qt.RightButton)][
           int(QtCore.Qt.NoModifier)] = self.pickMenu

        self.cameras[-1]._mouseReleaseNoMotionActions[
           int(QtCore.Qt.RightButton)][
           int(QtCore.Qt.ControlModifier)] = self.pickCenter

        self.customizeReprDialog = None

    ## def AddCamera(self, master=None, screenName=None, classCamera=None,
    ##               stereo='none', num=None, verbose=True, cnf={}, **kw):

    ##     Viewer.AddCamera(self, master=master, screenName=screenName,
    ##                      classCamera=classCamera, stereo=stereo, num=num,
    ##                      verbose=verbose,cnf=cnf, **kw)
        
    ##     self.cameras[-1]._mouseReleaseNoMotionActions[
    ##         int(QtCore.Qt.RightButton)][
    ##         int(QtCore.Qt.NoModifier)] = self.pickMenu


    def getFilename(self, filenames=None):
        filenames, filter = QtGui.QFileDialog.getOpenFileNames(
            self, self.tr("Read molecule or session"), '',
            self.tr("PDB Files (*.pdb);; PDBQT Files (*.pdbqt);; MOL2 Files (*.mol2);; SDF Files (*.sdf);; MMTF Files (*.mmtf);; Session Files (*.psf);; Docking Results Object Files(*.dro);; All files (*) "))
        if isinstance(filenames, str):
            filenames = [filenames]
        return filenames
    
    def openFile(self, filenames=None):
        if filenames is None:
            filenames = self.getFilename()
        if filenames:
            molfiles = []
            mols = []
            for i, filename in enumerate(filenames):
                #fn = filename.encode('ascii', 'ignore')
                #FIXME .. how to handle unicode names ?
                filenames[i] = filename.encode('ascii', 'replace')
                ext = os.path.splitext(filename)[1].lower()
                if ext == ".psf":
                    # Load Pmv session file
                    self.pmv.readFullSession(filename)
                elif ext == ".dro":
                    from ADFR.PmvInterface.fileCmds import DockingResultReader
                    self.pmv.addCommand(DockingResultReader(), 'loadDRO')
                    mols = self.pmv.loadDRO([str(filename)])
                else:
                    molfiles.append(filenames[i])
            if len(molfiles):
                mols = self.pmv.readMolecules.guiCall(molfiles, header=True)
            return mols

    def saveSession(self, filename=None):
        ## filename, filtr = QtGui.QFileDialog.getSaveFileName(
        ##     self, self.tr("Save session"), '',
        ##     self.tr("Session Files (*.psf)"))
        filename, filtr = QtGui.QFileDialog.getSaveFileName(
            self, "Save session", '', "Session Files (*.psf)")
        if not filename: return
        fi = QtCore.QFileInfo(filename)
        if len(fi.suffix())==0: filename += ".psf"

        # check for multimolecules and compute size
        sizeMB = 0
        nbmm = 0
        for mol in self.pmv.Mols:
            if mol._multi=='molecules':
                filesize = os.path.getsize(mol.filename)
                # compute size in MB
                sizeMB += int(round(filesize/1048576.))
                nbmm += 1
                
        inlcudeMultiMol = True
        if sizeMB > self.pmv.userpref['Save session multiMolSize']['value']:
            msg = "Your session contains a multi-molecule file with %d molecule(s) amounting to %dMB.\nIncluding the molecule in the saved session will allow restoring the multi-molecule from the session file.\nIf you choose to not include it, only the currently select molecule from that file will be saved.\n\nDo you want to include the compelte file in the session file?"%(nbmm, sizeMB)
            reply = QtGui.QMessageBox.question(
                self.pmv.gui(), "save session", msg,
                QtGui.QMessageBox.Yes | QtGui.QMessageBox.No| QtGui.QMessageBox.Cancel)
            if reply==QtGui.QMessageBox.Cancel:
                return
            inlcudeMultiMol = reply==QtGui.QMessageBox.Yes

        dialog =  SaveSessionDialog(self.pmv.saveFullSession,
                                    (filename, inlcudeMultiMol), estimatedSize=sizeMB)
        dialog.setWindowTitle("saving session %s"%filename)
        result = dialog.exec_()
        #print 'GOTBACK', result==dialog.Rejected
        #if result==dialog.Rejected:
        #    try:
        #        print 'removing'
        #        os.remove(filename)
        #    except:
        #        pass
        #self.pmv.runInThread(self.pmv.saveFullSession,
        #                     "saving session %s"%filename,
        #                     (filename, inlcudeMultiMol))

    def runInThread(self, func, description, args):
        # put up a modal form showing progress
        Dialog = QtGui.QDialog()
        ui = Ui_Dialog()
        ui.setupUi(Dialog)
        Dialog.show()
        
    def saveAs_cb(self, obj):
        mol = None
        if isinstance(obj, SelectionSet):
            mol = obj[0].getAtomGroup().getMolecule()
        elif isinstance(obj, Molecule):
            mol = obj
        if mol:
            ext = os.path.splitext(mol.filename)[1]
            extensions = [".pdb", ".mol2"]
            if ext in extensions:
                extensions.pop(extensions.index(ext))
                extensions.insert(0, ext)
            filename, filtr = QtGui.QFileDialog.getSaveFileName(
                self, "Save As", '', "Molecule (*%s);;Molecule (*%s)"%tuple(extensions))
            fi = QtCore.QFileInfo(filename)
            print filtr, QtCore.QFileInfo(filtr)
            if len(fi.suffix())==0:
                for e in extensions:
                    if filtr.find(e) > 0:
                        filename += e
                        break
            print filename
            self.pmv.saveSelectionAsMolecule.guiCall(obj, filename)

    def openMapZipFile(self, filename=None):
        if not filename:
            filename, filter = QtGui.QFileDialog.getOpenFileName(
                self, self.tr("Read Maps"), '',
                self.tr("zip Files (*.zip);; All files (*)"))
        if filename:
            self.pmv.readGrid.guiCall(filename)
    
    def setStatusBar(self, bar):
        self.statusBar = bar
        
    def setEmptyPickCommand(self, cmd):
        # force loading the selection command
        if cmd is not None:
            try:
                cmd = cmd.loadCommand()
            except AttributeError:
                pass
        # set the selection command to be the command
        self.emptyPickFunc = cmd

    def setDragSelectCommand(self, cmd):
        # force loading the selection command
        try:
            cmd = cmd.loadCommand()
        except AttributeError:
            pass
        # set the selection command to be the command
        self.dragSelectFunc = cmd

    def processPicking(self, pick):
        cam = self.currentCamera
        objects = []
        if len(pick.hits):
            selections = self.findPickedAtoms(pick)
            if selections.nbAtoms():
                if pick.operation=='add': negate = False
                else: negate = True
                self.dragSelectFunc(selections, negate=negate)
        else:
            if self.emptyPickFunc:
                self.emptyPickFunc()
        
    def findPickedAtoms(self, pick):
        """
given a PickObject this function finds all corresponding atoms.
Each atom in the returned set has its attribute pickedInstances set to a list
of 2-tuples [(geom, instance),...].
"""
        selections = SelectionSet()

        # loop over object, i.e. geometry objects
        for obj, values in pick.hits.items():
            # build a list of vertices and list of instances
            vertInds, instances = zip(*values)

            # only geometry bound to molecules is pickable in PMV
            if not hasattr(obj, 'mol') or len(vertInds)<1:
                continue

            mol = obj.mol()
            # only vertices of geometry have a mapping to atoms
            # for other types we return an empty atom set
            if pick.type!='vertices':
                return mol.emptySelection()

            gc = mol.geomContainer

            # convert vertex indices into atoms
            if gc.geomPickToAtoms.has_key(obj.name):
                # the geometry obj has a function to convert to atoms
                # specified it he geomContainer[obj], e.g. MSMS surface
                func = gc.geomPickToAtoms[obj.name]
                if func:
                    atInds = func(obj, vertInds)
                    selections.append( Selection(mol._ag, atInds, 'picked',
                                             mol._ag.getACSIndex()) )
            else:
                # we assume a 1 to 1 mapping of vertices with atoms
                # e.g. the lines geometry
                selections.append( Selection(mol._ag, vertInds, 'picked',
                                             mol._ag.getACSIndex()) )

        return selections
 
    def findPickedBonds(self, pick):
        """do a pick operation and return a 2-tuple holding (the picked bond,
        the picked geometry)"""

        allbonds = BondSet( [] )
        for o, val in pick.hits.items(): #loop over geometries
            # loop over list of (vertices, instance) (e.g. (45, [0,0,2,0]))
            for instvert in val:
                primInd = instvert[0]
                if not hasattr(o, 'mol'): continue
                g = o.mol.geomContainer
                if g.geomPickToBonds.has_key(o.name):
                    func = g.geomPickToBonds[o.name]
                    if func: allbonds = allbonds + func(o, primInd)
                else:
                    l = []
                    bonds = g.atoms[o.name].bonds[0]
                    for i in range(len(primInd)):
                        l.append(bonds[int(primInd[i])])
                    allbonds = allbonds + BondSet(l)

        return allbonds

    def buildMenuForCamera(self, camera, parent=None):

        def askCamColor():
            def setCamCol(color):
                camera.Set(color=color.getRgb())
                camera.viewer.OneRedraw()
            oldCol = camera.backgroundColor
            w = QtGui.QColorDialog(parent)
            w.setCurrentColor(QtGui.QColor(*oldCol))
            w.currentColorChanged.connect(setCamCol)
            value = w.exec_()
            if value==0:
                camera.Set(color=oldCol)
                camera.viewer.OneRedraw()
        
        menu = QtGui.QMenu(self.tr('Camera'), parent)
        action = menu.addAction("background Color ...")
        self.connect(action, QtCore.SIGNAL('triggered()'), askCamColor)
        menu.popup(QtGui.QCursor.pos())
        
    def transformedCoordinatesWithInstances(self, hits):
        """ hist is pick.hits = {geom: [(vertexInd, intance),...]}
This function will use the instance information to return a list of transformed
coordinates
"""
        vt = []
        for geom, values in hits.items():
            coords = geom.vertexSet.vertices.array
            for vert, instance in values:
                M = geom.GetMatrix(geom.LastParentBeforeRoot(), instance[1:])
                pt = coords[vert]
                ptx = M[0][0]*pt[0]+M[0][1]*pt[1]+M[0][2]*pt[2]+M[0][3]
                pty = M[1][0]*pt[0]+M[1][1]*pt[1]+M[1][2]*pt[2]+M[1][3]
                ptz = M[2][0]*pt[0]+M[2][1]*pt[1]+M[2][2]*pt[2]+M[2][3]
                vt.append( (ptx, pty, ptz) )
        return vt

    def pickCenter(self, x, y, dx, dy, e):
        # set center of rotation to the picked atom
        cam = self.currentCamera
        pick = cam.DoPick(x, y, event=e)
        if len(pick.hits)==0: return
        vt = self.transformedCoordinatesWithInstances(pick.hits)
        g = [0.,0.,0.]
        i = 0
        for v in vt:
            g[0] += v[0]
            g[1] += v[1]
            g[2] += v[2]
            i+=1
        g[0] = g[0]/i
        g[1] = g[1]/i
        g[2] = g[2]/i

        vi = self
        if self.redirectTransformToRoot or self.currentObject==self.rootObject:
            self.currentObject.SetPivot(g)
        else:
            m = self.currentObject.GetMatrixInverse(root=self.currentObject)
            g.append(1.0)
            g = Numeric.dot(m, g)[:3]
            self.currentObject.SetPivot(g)
        
    def pickMenu(self, x, y, dx, dy, e):
        cam = self.currentCamera
        pick = cam.DoPick(x, y, event=e)
        objects = []
        if len(pick.hits):
            atoms = self.findPickedAtoms(pick)
            if atoms.nbAtoms():
                obj = atoms[0] # obj is a MolKit Selection
                atIndex = obj.getIndices()[0]
                mol = obj.getAtomGroup().getMolecule()
                objItems = []
                levels = []
                # now check for groups
                #ob = objItems[-1][0]
                #while ob._group:
                #    objItems.append( (ob._group, ob._group._treeItems.keys()[0]) )
                #    ob = ob._group

                # check for selections
                sele = self.pmv.curSelection.selectionFor(mol)
                if sele:
                    if atIndex in sele.atomsDict:
                        objItems.append( (self.pmv.curSelection, None) )
                for selSet in self.pmv.namedSelections.values():
                    sele = selSet.selectionFor(mol)
                    if atIndex in sele.atomsDict:
                        objItems.append( (selSet, None) )

                #molecule
                sele = SelectionSet( [mol.select()] )
                sele.name = mol.name
                objItems.append( (sele, None) )
                levels.append('molecule')

                # chains
                ag = mol._ag
                chid = mol._ag.getChids()[atIndex]
                chainSel = mol.select('chid "%s"'%chid)
                sele = SelectionSet( [chainSel] )
                sele.name = chid
                objItems.append( (sele, None) )
                levels.append('chain')

                # residue
                resname = ag.getResnames()[atIndex]
                resnum = ag.getResnums()[atIndex]
                resSel = mol.select("resname %s resnum %d"%(resname, resnum))
                sele = SelectionSet( [resSel] )
                sele.name = '%s%d'%(resname, resnum)
                objItems.append( (sele, None) )
                levels.append('residue')

                # atom
                atom = Atom(obj.iterAtoms().next())
                sele = SelectionSet([atom.select()])
                sele.name = atom.name
                objItems.append( (sele, None) )
                levels.append('atom')

                self.buildMenuForObjects(objItems, levels, parent=self)
                #print 'pick object', atom, pick.hits.keys()[0].name
        else:
            self.buildMenuForCamera(self.currentCamera, parent=self)

    ## def getVisibleGeoms(self, obj, visibleOnly=True):
    ##     geoms = {}
    ##     moleculesAndSelections = [[x.getAtomGroup().getMolecule(),x] for x in obj]
    ##     for mol, sele in moleculesAndSelections:
    ##         gc = mol.geomContainer
    ##         for name, ats in gc.atoms.items():
    ##             if gc.geoms.has_key(name):
    ##                 if visibleOnly and not gc.geoms[name].visible:
    ##                     continue
    ##                 if gc.displayedAs([name], sele, 'fast'):
    ##                     if name in ['balls', 'sticks']:
    ##                         geoms['sticksAndBalls'] = gc.geoms[name]
    ##                     elif name in ['noBond', 'singleBonds', 'doubleBonds', 'tripleBonds', 'aromaticBonds']:
    ##                         geoms['lines'] = gc.geoms[name]
    ##                     elif name.startswith("chain_") and gc.geoms[name].parent.name == "cartoon":
    ##                         geoms["cartoon %s"%mol._basename] = gc.geoms[name]
    ##                     else:
    ##                         geoms[name] = gc.geoms[name]
    ##     return geoms

    def deleteSelection_cb(self, objItems):
        # items is a list of items from the dashboard
        # the items all correspond to object of the same type
        rstr = ""
        for ob, it in objItems:
            rstr += ob.name+', '
        rstr = rstr[:-2]
        reply = QtGui.QMessageBox.question(self, "Delete Selection",
                                           "do you want want to delete the selection(s)\n%s?"%rstr,
                                        QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        if reply == QtGui.QMessageBox.Yes:
            for ob, it in objItems:
                self.pmv.deleteNamedSelection(ob)

    def deleteObjects_cb(self, objItems):
        # items is a list of items from the dashboard
        # the items all correspond to object of the same type
        rstr = ""
        obj0, it0 = objItems[-1]
        #if isinstance(obj0, (Protein, PmvSelection, MoleculeGroup)):
        if isinstance(obj0, (Molecule, SelectionSet, Group)):
            for obj, it in objItems:
                rstr += obj.name+', '
        else:
            rstr = ''
            for obj, it in objItems:
                rstr += obj.name + ' '

        reply = QtGui.QMessageBox.question(self, "Delete",
                                           "do you want to delete %s ?"%rstr,
                                        QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        if reply == QtGui.QMessageBox.Yes:
            self.delete(objItems)

    def delete(self, objItems):
        app = self.pmv
        allObjects = [x[0] for x in objItems]
        #allObjects = allObjects[0].setClass(allObjects) # ???

        #if isinstance(allObjects[0], Protein):
        if isinstance(allObjects[0], Molecule):
            self.pmv.deleteMolecule.guiCall(MoleculeSet("molset", allObjects))
        elif isinstance(allObjects[0], SelectionSet):
            ## atoms = allObjects[0].atoms
            ## for obj in allObjects[1:]:
            ##     atoms += obj.atoms
            ## self.pmv.deleteAtoms.guiCall(atoms)
            self.pmv.deleteAtoms.guiCall(allObjects[0])
        elif isinstance(allObjects[0], Group):
            self.pmv.deleteGroups(allObjects)
        else:
            sel = allObjects[0].selSet
            for obj in allObjects[1:]:
                sel |= obj.selSet
            self.pmv.deleteAtoms.guiCall(sel)

    def rename_cb(self, objItems):
        # items is a list of items from the dashboard
        # the items all correspond to object of the same type

        w = GetNewName(objItems, self)
        ok = w.exec_()
	if ok==QtGui.QDialog.Accepted:
	    name = w.textValue()
	    name = name.encode('ascii', 'ignore')
            names = name.split(',')
            if len(names)==1:
                names = names*len(objItems)            
            # remove leading and ending spaces
            names = [n.strip() for n in names]
	else:
            names = None
	# this version does not allow me to overwrite done()
        #name, ok = GetSelectionName.getText(self, 'Name', 'Selection name:')
        if ok and names:
            self.rename(names, objItems)

    def rename(self, names, objItems):
        app = self.pmv
        for name, objItem in zip(names, objItems):
            obj, item = objItem
            self.pmv.rename(obj, name, item)

    def buildMenuForObjects(self, objItems, levels=None, parent=None):
        # this function is call from processPick with a list of 4 (obj,item)
        # for Molecule, Chain, Residue, Atom that was picked
        # We build a menu that will cascade the dashboard menu for each of these objects
        # assuming the the item picked was the Molecule object in the dashboard
        # We also build menu entries for the current selection if the picked atom is selected
        # and antries for each names selection the picked atom belings to
        menu = QtGui.QMenu(self.tr('Levels'), parent)

        for objIt in objItems[:-4]: # last 4 are Mol, Chain, res, atom
            # ones ebfore can be seelections or group
            # Current Selection Entry
            ob, it = objIt
            if ob == self.pmv.curSelection: #activeSelection.isAnySelected(ob):
                submenu = menu.addMenu(self.icons['selection'], 'Current Selection')
                item = None #self.pmv.activeSelection._treeItems.keys()[0]
                ob = self.pmv.activeSelection
                self.buildMenuForObject( [(ob, item)], level=4, parent=parent, parentMenu=submenu)
            elif isinstance(ob, SelectionSet):
                submenu = menu.addMenu(self.icons['selection'], ob.name)
                item = None #ob._treeItems.keys()[0]
                self.buildMenuForObject( [(ob, item)], level=4, parent=parent, parentMenu=submenu)
            elif isinstance(ob, Group):
                raise
            
        level = None
        if levels: level = levels[-1]

        # build dashboard menu for Molecule, Chain, Residue and Atom
        iconNames = ['molecule', 'chain', 'residue', 'atom']
        for num, objIt in enumerate(objItems[-4:]):
            if levels: level = levels[num]
            obj, item = objIt
            if num < 4:
                label = self.levelLabels[num]
            else:
                label = 'G:'

            submenu = menu.addMenu(self.icons[iconNames[num]], '%s %s'%(label, getattr(obj, 'name')))
            self.buildMenuForObject( [(obj, item)], level=level, parent=parent, parentMenu=submenu )
            
        menu.popup(QtGui.QCursor.pos())
        return menu
            
    def getMoleculesInGroups(self, items):
        # gets all the molecule in the subtrees rooted at items
        # items are supposed to be for groups / molecules
        molecules = MoleculeSet([])
        for item in items:
            molecules.extend( self.getMoleculesInGroup(item) )
        return molecules
            
    def getMoleculesInGroup(self, item):
        molecules = MoleculeSet([])
        obj = item._pmvObject
        if isinstance(obj, Protein):
            molecules.append(obj)
        elif isinstance(obj, Group):
            # for groups return all molecules in subtree
            for n in range(item.childCount()):
                child = item.child(n)
                molecules.extend( self.getMoleculesInGroup(child) )
        return molecules

    def repostMenus(self, menu):
        ## to repost a menu we need to repost all parent menus at
        ## their respective locations and register a callback which
        ## will remove all reposted menus once a menu entry that does
        ## not trigger a repost is triggered

        ## create a list of all menus
        parent = menu
        parents = []
        while isinstance(parent, QtGui.QMenu):
            parents.append(parent)
            parent = parent.parent()

        # reverse the list to post them in the order the user posted
        # the menus (as they overlap slightly
        parents.reverse()

        # post the menus
        for parent in parents:
            parent.popup(parent.pos())

        # define a function that will unpost them all
        def hideAll(menu):
            parent = menu
            while isinstance(parent, QtGui.QMenu):
                parent.setVisible(False)
                parent = parent.parent()

        # assign the function to hide all menus
        cb = CallbackFunction(hideAll, menu)
        menu.aboutToHide.connect(cb)

        ## menu.popup(menu.pos()) # this reposts the menu at its position but not parent menus
    def showHide(self, items):
        for item in items:
            self.pmv.gui().objTree.showHide(item, 0, expandCollapse=False)
            
    def getAddHOpt_cb(self, molSel):
        from .addHGUI import AddHydrogenGUI
        molname = molSel[0].getAtomGroup().getMolecule().name
        opt, ok = AddHydrogenGUI.getAddHPrams(molname)
        if ok:
            if opt[0]=='reduce':
                self.pmv.protonateWithReduce.guiCall(molSel, **opt[1])
            else:
                self.pmv.protonateWithOpenBabel.guiCall(molSel, **opt[1])

    def getGeomParams(self, selection, objItems, gname=None):
        from .advancedParamGui import AdvancedParamGUI
        obj0, item0 = objItems[0]
        rootItem = item0
        if rootItem:
            while rootItem.parent() is not None:
                rootItem = rootItem.parent()
            # for Selections, there can only be one selection we are operating on
            pmvObject = None
            if isinstance(rootItem._pmvObject, PmvSelection):
                pmvObject = rootItem._pmvObject
        else:
            pmvObject = obj0
        opt, ok = AdvancedParamGUI.getParams(selection, func=self.setGeomProperties, gname=gname, pmvObject=pmvObject)
        #print "getGeomParams", opt
        if ok:
            self.setGeomProperties(selection, pmvObject, opt, gname)
        #print opt, ok

    def setGeomProperties(self, selection, pmvObject, opt, gname=None):
        if not len(opt): return
        if opt[0] == "Properties":
            cmdDict = {"Lines": self.pmv.displayLines.guiCall,
                       "CPK" : self.pmv.displayCPK.guiCall,
                       "Sticks And Balls": self.pmv.displaySB.guiCall,
                       }
            geom = opt[1] #geometry category from AdvancedParamGUI "Geom representation" list box; corresponds to one of the kyes from cmdDict.
            #print "setGeomProperties:", geom, gname
            kw = opt[2] # geometry parameters from AdvancedParamGUI.
            if cmdDict.has_key(geom):
                cmdDict[geom](selection, **kw)
            elif geom == "Surface":
                recompute=False
                for sele in selection:
                    mol = sele.getAtomGroup().getMolecule()
                    if gname is not None:
                        surfName = gname
                    else:
                        if pmvObject and pmvObject is not self.pmv.curSelection:
                            if pmvObject._multi=='conformations':
                                surfName = 'surface-%s'%pmvObject._basename
                            else:
                                surfName = 'surface-%s'%pmvObject.name
                        else:
                            if mol._multi=='conformations':
                                surfName = 'surface-%s'%mol._basename
                            else:
                                surfName = 'surface-%s'%mol.name
                    # find last parameters used to compute this surface
                    
                    if not hasattr(mol.geomContainer, "msms"):
                        recompute = True
                    else:
                        try:
                            molInd = mol._ag.getACSIndex()
                            srf = mol.geomContainer.msms[surfName][molInd][0]
                            lastValues = {
                                'pRadius': srf.probeRadius,
                                'density': srf.density,
                                #'perMol': srf.perMol,
                                #'surfName': srf.surfName,
                                #'noHetatm': srf.noHetatm,
                                'hdset': srf.hdset,
                                'hdensity': srf.hdensity,
                            }
                        except KeyError:
                            lastValues = {}
                            recompute = False
                        if len(lastValues)==0: # surface was not computed yet
                            recompute=True
                        else: # compare last parameters with ones in form
                            for k,v in lastValues.items():
                                nv = kw.get(k, None)
                                if nv!=v:
                                    #print 'MSMS recompute: new param', surfNme, k, nv, v
                                    recompute=True
                    if recompute:
                        perMol = True
                        if pmvObject and pmvObject is not self.pmv.curSelection:
                            # molecule or named selection
                            perMol=False
                            #surfName = 'surface-%s'%pmvObject.name
                        kw['perMol'] = perMol
                        kw['surfName']=surfName
                        kw['display']=False
                        self.pmv.computeMSMS.guiCall( *(sele,), **kw)
                        mol.geomContainer.geoms[surfName]._msmsType = True
                    self.pmv.displayMSMS.guiCall(sele, surfNames=[surfName])
            elif geom == "CoarseMS":
                recompute=False
                for sele in selection:
                    mol = sele.getAtomGroup().getMolecule()
                    if gname is not None:
                        surfName = gname
                    else:
                        if pmvObject and pmvObject is not self.pmv.curSelection:
                            surfName = 'coarseMolSurf-%s'%pmvObject.name
                        else:
                            surfName = 'coarseMolSurf-%s'%mol.name
                    # find last parameters used to compute this surface
                    if not hasattr(mol.geomContainer, "boundGeom"):
                        recompute = True
                    elif not mol.geomContainer.boundGeom.has_key(surfName):
                        recompute = True
                    else:
                        try:
                            srfParams = mol._coarseMolSurfParams[surfName]
                            lastValues = { 'gridSize':srfParams['gridSize'],
                                           'padding':srfParams['padding'],
                                           'resolution':srfParams['resolution'],
                                           'isovalue':srfParams['isovalue']}

                        except:
                            lastValues = {}
                            recompute = True
                        if len(lastValues): # compare last parameters with ones in form
                            for k,v in lastValues.items():
                                nv = kw.get(k, None)
                                if nv!=v:
                                    #print 'coarseMS recompute: new param', surfNme, k, nv, v
                                    recompute=True
                                    break
                    if recompute:
                        kw['surfName']=surfName
                        # if pmvObject is a molecule or currentSelection, compute surface for all molecule's atoms.
                        # Display surface for the selected atoms.
                        atoms = sele
                        if pmvObject and isinstance(pmvObject, PmvSelection):
                            if pmvObject is self.pmv.curSelection:
                                atoms = mol.select()
                        else:
                            atoms = mol.select()
                        self.pmv.computeCoarseMolecularSurface.guiCall(*(atoms,), **kw)
                    self.pmv.displayBoundGeom.guiCall(sele, surfName)
                    
            elif geom in ["atomLabels", "residueLabels"]:
                for sele in selection:
                    mol = sele.getAtomGroup().getMolecule()
                    if mol.geomContainer.geoms.get(geom, None):
                        labgeom = mol.geomContainer.geoms[geom]
                        #print "label geom:", labgeom, kw
                        if kw.has_key('labelPlacement'):
                            labPlacement = kw.pop('labelPlacement')
                            recompute = True
                            if labgeom.params.has_key('mode'):
                                if labgeom.params['mode'] != labPlacement:
                                    recompute = True
                                else:
                                    recompute = False
                            if recompute:
                                self.pmv.labelResidues.guiCall(sele, font=kw['font'], mode=labPlacement)
                        labgeom.Set(**kw)

    def buildMenuForObject(self, objItems, level=None, parent=None, parentMenu=None,
                           dynamicContextMenu=True):
        # objItems is a list of (obj, item) pairs. item can be None when this is called
        # obj is a SelectionSet or a MolKit2.molecule.Molecule Chain Residue or Atom that has a .selSet
        
        # the 3D Viewer pick as the item might not exist for a residue for instance when
        # the molecule was not yet expanded.
        # The list can be of length one (e.g. right click on single entry)
        # or have N items if right click on entry that is part of dashboard selection
        # parent menu is used to decide if we post the menu or it will be a cascade of the
        # the parent menu (used by context menu in 3D viewer)
        t0 = time()
        if parentMenu is None:
            menu = QtGui.QMenu(parent)
        else:
            menu = parentMenu

        curSel = self.pmv.curSelection

        object0, item0 = objItems[0]

        # add showHide entry
        #import pdb; pdb.set_trace()
        if item0: # could be None which righ click in 3d View
            root = self.pmv.gui().objTree.getRootItem(item0)
            if root is item0: # we are at the top
                action = menu.addAction("show/Hide")
                cb = CallbackFunction(self.showHide, [it[1] for it in objItems])
                action.triggered.connect(cb)
                if not isinstance(object0, Group):
                    actionSubmenu = menu.addMenu('actions')
                    selSet = SelectionSet()
                    for obj, it in objItems:
                        if isinstance(obj, SelectionSet):
                            selSet = selSet | obj
                        elif isinstance(obj, Molecule):
                            selSet.append(obj.select())
                        
                    ## action = actionSubmenu.addAction('add hydrogen')
                    ## cb = CallbackFunction(self.getAddHOpt_cb, selSet)
                    ## action.triggered.connect(cb)

                    action = actionSubmenu.addAction("assign atom and boond types")
                    cb = CallbackFunction(self.pmv.typeAtomsAndBonds.guiCall, selSet)
                    action.triggered.connect(cb)

                    action = actionSubmenu.addAction("biological Unit")
                    cb = CallbackFunction(self.pmv.biologicalUnit, selSet)
                    action.triggered.connect(cb)

                    ## action = actionSubmenu.addAction("minimize")
                    ## cb = CallbackFunction(self.pmv.minimize.guiCall, selSet)
                    ## action.triggered.connect(cb)

        if item0 is None: # use the molecule's root item
            try:
                item0 = object0[0].getAtomGroup().getMolecule()._treeItems.keys()[0]
                # find the level of objects we are building a menu for
                if level is None:
                    level = self.objTree.getLevel(item0)
            except AttributeError: # happend in PmvViewer without dashboard
                level = 3 # molecule
        else:
            root = self.pmv.gui().objTree.getRootItem(item0)

        isSelection = False

        # build the PmvSelection for all lines highlighted in the dashboard
        # all highlighted lines are of the same type (eg. all molecules or PmvSelection)
        selected = PmvSelection(self.pmv)
        if isinstance(object0, (PmvSelection, SelectionSet)):
            isSelection = True
            for ob, it in objItems:
                selected = selected | ob
                
        elif isinstance(object0, Group):
            #raise
            # at this point only show/Hide menu has created
            if parentMenu is None:
                menu.popup(QtGui.QCursor.pos())
            return menu
            #items = [ob[0]._treeItems.keys()[0] for ob in objItems]
            #molecules = self.getMoleculesInGroups(items)
            #atoms = molecules.allAtoms
        else:
            for ob, it in objItems:
                selected = selected | SelectionSet([ob.select()])
 
        # if the item is in a selection, restrict atoms to selection
        #allAtoms = SelectionSet()
        #for sel in selected:
        #    allAtoms.append(sel)

        #if isinstance(rootItem._pmvObject, PmvSelection):
        #    selected = selected & rootItem._pmvObject

        ##
        ## Select/Deselect menu entry
        if selected.nbAtoms():
            AddDeselected = self.pmv.activeSelection.allSelected(selected)
            if AddDeselected is True:
                action = menu.addAction("Deselect")
                cb = CallbackFunction(self.pmv.select.guiCall, selected, negate=1)
                self.connect(action, QtCore.SIGNAL('triggered()'), cb)
            else:
                if self.pmv.activeSelection.isAnySelected(selected): # at least one selected
                    # objects in selection add to the selection.e.g. complete a residue
                    action = menu.addAction("Complete Partial Selection")
                    cb = CallbackFunction(self.pmv.select.guiCall, selected, negate=0, only=False)
                    self.connect(action, QtCore.SIGNAL('triggered()'), cb)

                # set as selection
                action = menu.addAction("Set as Selection")
                cb = CallbackFunction(self.pmv.select.guiCall, selected, negate=0, only=True)
                self.connect(action, QtCore.SIGNAL('triggered()'), cb)

                ##
                ## selection and deselection subsmenus
                selsubmenu = menu.addMenu('select')
                self.connect(selsubmenu, QtCore.SIGNAL("hovered(QAction *)"),
                             self._showmenuStatusTip)
                if self.pmv.activeSelection.nbAtoms() > 0:
                    deselsubmenu = menu.addMenu('deselect')
                else:
                    deselsubmenu = None

                # handle 'all'
                action = selsubmenu.addAction('all')
                cb = CallbackFunction(self.pmv.select.guiCall, selected, negate=0)
                self.connect(action, QtCore.SIGNAL('triggered()'), cb)
                selsubmenu.addSeparator()

                if deselsubmenu:
                    action = deselsubmenu.addAction('all')
                    cb = CallbackFunction(self.pmv.select.guiCall, selected, negate=1)
                    self.connect(action, QtCore.SIGNAL('triggered()'), cb)
                    deselsubmenu.addSeparator()

                if dynamicContextMenu:
                    for menuName, whatList in selectionTypes.items():
                        targets = []
                        nb = 0
                        for wlist in whatList:
                            if wlist is not 'Separator':
                                what, name, minLevel, doc = wlist
                                if level >= minLevel:
                                    atoms = selected.select(what)
                                    if atoms is not None:
                                        nb += atoms.nbAtoms()
                                    targets.append(atoms)
                                else:
                                   targets.append(None) 
                            else:
                                targets.append(None)
                        if nb:
                            selsubsubmenu = selsubmenu.addMenu(menuName)
                            if deselsubmenu:
                                deselsubsubmenu = deselsubmenu.addMenu(menuName)
                            for n, wlist in enumerate(whatList):
                                if targets[n]==None: continue
                                if wlist is 'Separator':
                                    selsubsubmenu.addSeparator()
                                    if deselsubmenu:
                                        deselsubsubmenu.addSeparator()
                                else:
                                    what, name, minLevel, doc = wlist
                                    action = selsubsubmenu.addAction(what)
                                    if doc:
                                        action.setStatusTip('select '+doc)
                                    cb = CallbackFunction(self.pmv.select.guiCall, targets[n], negate=0)
                                    self.connect(action, QtCore.SIGNAL('triggered()'), cb)

                                    if deselsubmenu:
                                        atoms = self.pmv.activeSelection.select(what)
                                        if atoms and len(atoms):
                                            action = deselsubsubmenu.addAction(what)
                                            if doc:
                                                action.setStatusTip('deselect '+doc)
                                            cb = CallbackFunction(self.pmv.select.guiCall, targets[n], negate=1)
                                            self.connect(action, QtCore.SIGNAL('triggered()'), cb)
                else:
                    for menuName, whatList in selectionTypes.items():
                        nb = 0
                        for wlist in whatList:
                            if wlist is not 'Separator':
                                what, name, minLevel, doc = wlist
                                if level >= minLevel:
                                    nb += 1
                        if nb:
                            selsubsubmenu = selsubmenu.addMenu(menuName)
                            if deselsubmenu:
                                deselsubsubmenu = deselsubmenu.addMenu(menuName)
                            for n, wlist in enumerate(whatList):
                                if wlist is 'Separator':
                                    selsubsubmenu.addSeparator()
                                    if deselsubmenu:
                                        deselsubsubmenu.addSeparator()
                                else:
                                    what, name, minLevel, doc = wlist
                                    if level < minLevel: continue
                                    action = selsubsubmenu.addAction(what)
                                    if doc:
                                        action.setStatusTip('select '+doc)
                                    cb = CallbackFunction(self.select, selected, what, negate=0)
                                    self.connect(action, QtCore.SIGNAL('triggered()'), cb)

                                    if deselsubmenu:
                                        atoms = self.pmv.activeSelection.select(what)
                                        if atoms and len(atoms):
                                            action = deselsubsubmenu.addAction(what)
                                            if doc:
                                                action.setStatusTip('deselect '+doc)
                                            cb = CallbackFunction(self.select, selected, what, negate=1)
                                            self.connect(action, QtCore.SIGNAL('triggered()'), cb)
                        
        ##
        ## add to selections menu entry
        seleCount = len(self.pmv.namedSelections)
        curSelInTree = hasattr(self.pmv.curSelection, '_treeItems')
        if selected.nbAtoms() and (seleCount or curSelInTree):
            addSubmenu = menu.addMenu('Add to Selection')
            invertAction = menu.addAction('Invert Selection')
            if curSelInTree:
                pmvsel = SelectionSet(curSel)
                pmvsel.extend(selected)
                action1 = addSubmenu.addAction(self.tr('Current Selection'))
                cb1 = CallbackFunction(self.pmv.select.guiCall, pmvsel, negate=0, only=False)
                self.connect(action1, QtCore.SIGNAL('triggered()'), cb1)
                if selected != self.pmv.curSelection:
                    # SelectAround and ExpandSelection commands expand current selection only to atoms that belong to the "selected". So it makes no sense to add these actions to the menu if "selected" is the Current Selection.
                    expandSubmenu = menu.addMenu('Expand Selection')
                    selaroundSubmenu = menu.addMenu('Select Around')
                    action2 = expandSubmenu.addAction(self.tr('Current Selection'))
                    action3 = selaroundSubmenu.addAction(self.tr('Current Selection'))
                    cb2 = CallbackFunction(self.setCutOff, self.pmv.expandSelection, SelectionSet(curSel), selected)

                    self.connect(action2, QtCore.SIGNAL('triggered()'), cb2)
                    #cb3 = CallbackFunction(self.setCutOff, self.pmv.selectAround.guiCall, SelectionSet(curSel), selected)
                    cb3 = CallbackFunction(self.setCutOff, self.pmv.selectAround, SelectionSet(curSel), selected)
                    self.connect(action3, QtCore.SIGNAL('triggered()'), cb3)

                cb4 = CallbackFunction(self.pmv.invertSelection.guiCall, SelectionSet(curSel))
                self.connect(invertAction, QtCore.SIGNAL('triggered()'), cb4)
            for name, sel in self.pmv.namedSelections.items():
                action = addSubmenu.addAction(self.tr(name))
                pmvsel = SelectionSet(sel)
                pmvsel.extend(selected)
                cb = CallbackFunction(self.pmv.select.guiCall, pmvsel, negate=0, only=False)
                self.connect(action, QtCore.SIGNAL('triggered()'), cb)
                
        menu.addSeparator()

        ##
        ## Name menu entry
        action = menu.addAction("Name ...")
        cb = CallbackFunction(self.rename_cb, objItems)
        self.connect(action, QtCore.SIGNAL('triggered()'), cb)

        ##
        ## Delete menu entry
        action = menu.addAction("Delete ...")
        cb = CallbackFunction(self.deleteObjects_cb, objItems)
        self.connect(action, QtCore.SIGNAL('triggered()'), cb)

        ##
        ## Clone menu entry
        action = menu.addAction("Clone ...")
        cb = CallbackFunction(self.pmv.clone, objItems[0][0])
        self.connect(action, QtCore.SIGNAL('triggered()'), cb)

        ##
        ## Save As
        action = menu.addAction("Save As ...")
        cb = CallbackFunction(self.saveAs_cb, objItems[0][0])
        self.connect(action, QtCore.SIGNAL('triggered()'), cb)

        menu.addSeparator()
        action = menu.addAction("Focus ")
        cb = CallbackFunction(self.pmv.focusScene, selected, self)
        self.connect(action, QtCore.SIGNAL('triggered()'), cb)
        
        ## Delete Selection menu entry
        if isSelection:
            namedSelections = [x for x in objItems if x[0]!=curSel]
            if len(namedSelections):
                action = menu.addAction("Delete Selection")
                cb = CallbackFunction(self.deleteSelection_cb, namedSelections)
                self.connect(action, QtCore.SIGNAL('triggered()'), cb)
            
        menu.addSeparator()
        ## add all visible geometries
        ##
        geomsDict = self.pmv.getGeoms(selected)
        geomNames = geomsDict.keys()
        geomNames.sort()

        ## for gname in geomNames:
        ##     geom = geomsDict[gname]
        ##     mol = geom.mol()
        ##     submenu = menu.addMenu(gname)
        ##     action = submenu.addAction(self.tr('Remove'))
        ##     if hasattr(geom, '_msmsType'):
        ##         cb = CallbackFunction(self.pmv.displayMSMS.guiCall, selected, surfName=gname, negate=True)
        ##     elif gname == 'atomLabels':
        ##         cb = CallbackFunction(self.pmv.labelAtoms.guiCall, selected, negate=True)
        ##     elif gname == 'residuelabels':
        ##         cb = CallbackFunction(self.pmv.labelResidues.guiCall, selected, negate=True)
        ##     elif hasattr(geom, '_boundGeomType'):
        ##         cb = CallbackFunction(self.pmv.displayBoundGeom.guiCall, selected, gname, negate=True)
        ##     elif geom.parent.name == "cartoon": # (ribbon geom)
        ##         for sel in selected:
        ##             if sel.getAtomGroup().getMolecule() == mol:
        ##                 cb = CallbackFunction(self.pmv.displayCartoon.guiCall, sel, negate=True)
        ##     else:
        ##         cb = CallbackFunction(self.displayfor, selected, gname, negate=True)
        ##     self.connect(action, QtCore.SIGNAL('triggered()'), cb)

        ##     # add color submenu
        ##     colorSubMenu = submenu.addMenu(self.tr('Color by'))

        ##     # Ca only for coloring action
        ##     COnly = self.COnly = colorSubMenu.addAction(self.tr('Carbon only'))
        ##     COnly.setStatusTip('Apply coloring to Carbon atoms only')
        ##     COnly.setCheckable(True)
        ##     COnly.setChecked(False)
        ##     colorSubMenu.addSeparator()

        ##     cb = CallbackFunction(self.repostMenus, colorSubMenu)
        ##     self.connect(COnly, QtCore.SIGNAL('triggered()'), cb)
            
        ##     if gname != 'residuelabels':
        ##         action = colorSubMenu.addAction(self.tr('Atom Type'))
        ##         cb = CallbackFunction(self.color, selected, [gname], 'atomType')
        ##         self.connect(action, QtCore.SIGNAL('triggered()'), cb)

        ##     ## action = colorSubMenu.addAction(self.tr('Rasmol Residue Colors'))
        ##     ## cb = CallbackFunction(self.color, selected, [gname], 'rasmol')
        ##     ## self.connect(action, QtCore.SIGNAL('triggered()'), cb)

        ##     ## action = colorSubMenu.addAction(self.tr('Shapely Residues Colors'))
        ##     ## cb = CallbackFunction(self.color, selected, [gname], 'shapely')
        ##     ## self.connect(action, QtCore.SIGNAL('triggered()'), cb)

        ##     ## action = colorSubMenu.addAction(self.tr('Polarity'))
        ##     ## cb = CallbackFunction(self.color, selected, [gname], 'polarity')
        ##     ## self.connect(action, QtCore.SIGNAL('triggered()'), cb)

        ##     action = colorSubMenu.addAction(self.tr('Chain'))
        ##     cb = CallbackFunction(self.color, selected, [gname], 'chain')
        ##     self.connect(action, QtCore.SIGNAL('triggered()'), cb)

        ##     action = colorSubMenu.addAction(self.tr('Molecule'))
        ##     cb = CallbackFunction(self.color, selected, [gname], 'molecule')
        ##     self.connect(action, QtCore.SIGNAL('triggered()'), cb)

        ##     if gname != 'residuelabels':
        ##         action = colorSubMenu.addAction(self.tr('N to C Rainbow'))
        ##         cb = CallbackFunction(self.color, selected, [gname], 'rainbow')
        ##         self.connect(action, QtCore.SIGNAL('triggered()'), cb)
        ##         action = colorSubMenu.addAction(self.tr('Sec Structure'))
        ##         cb = CallbackFunction(self.color, selected, [gname], 'secstructure')
        ##         self.connect(action, QtCore.SIGNAL('triggered()'), cb)

        ##     action = colorSubMenu.addAction(self.tr('Custom color'))
        ##     cb = CallbackFunction(self.color, selected, [gname], 'customColor')
        ##     self.connect(action, QtCore.SIGNAL('triggered()'), cb)

        ##     action = colorSubMenu.addAction(self.tr('Molecules Rainbow'))
        ##     cb = CallbackFunction(self.color, selected, [gname], 'moleculesRainbow')
        ##     self.connect(action, QtCore.SIGNAL('triggered()'), cb)

        ##     action = colorSubMenu.addAction(self.tr('Rainbow By Chain'))
        ##     cb = CallbackFunction(self.color, selected, [gname], 'rainbowChain')
        ##     self.connect(action, QtCore.SIGNAL('triggered()'), cb)

        ##     action = submenu.addAction(self.tr('Customize'))
        ##     cb = CallbackFunction(self.getGeomParams, selected, objItems, gname=gname)
        ##     self.connect(action, QtCore.SIGNAL('triggered()'), cb)

        # add submenu to add representation
        submenu = menu.addMenu(self.tr('Add Representation'))

        action = submenu.addAction(self.tr('lines'))
        cb = CallbackFunction(self.displayfor, selected, 'lines', negate=0)
        self.connect(action, QtCore.SIGNAL('triggered()'), cb)

        action = submenu.addAction(self.tr('Sticks And Balls'))
        cb = CallbackFunction(self.displayfor, selected, 'sb', negate=0)
        self.connect(action, QtCore.SIGNAL('triggered()'), cb)

        action = submenu.addAction(self.tr('CPK'))
        cb = CallbackFunction(self.displayfor, selected, 'cpk', negate=0)
        self.connect(action, QtCore.SIGNAL('triggered()'), cb)
        #  Hyper balls (do not add for now)
        ## subsubmenu = submenu.addMenu('HyperBalls')
        ## for name, kw in {"S&B":{'shrink':0.01, 'scaleFactor':0.0 ,  'bScale':0.5, 'cpkRad':0.5 },
        ##                  "CPK":{'shrink':0.01, 'scaleFactor':1   ,  'bScale':0.01, 'cpkRad':0.0},
        ##                  "LIC":{'shrink':0.01, 'scaleFactor':0.0,  'bScale':1.0, 'cpkRad':0.3},
        ##                  "HBL":{'shrink':0.3,  'scaleFactor':0.0 ,  'bScale':1.0, 'cpkRad':0.6 }}.items():
        ##     action = subsubmenu.addAction(self.tr(name))
        ##     cb = CallbackFunction(self.displayfor, selected, 'hpballs', negate=0, **kw)
        ##     self.connect(action, QtCore.SIGNAL('triggered()'), cb)

        # FIXME only add of atoms is suitable for ribbon (see old dashboard)
        action = submenu.addAction(self.tr('cartoon'))
        cb = CallbackFunction(self.displayfor, selected, 'cartoon', negate=0)
        self.connect(action, QtCore.SIGNAL('triggered()'), cb)

        action = submenu.addAction(self.tr('surface'))
        cb = CallbackFunction(self.displayMSMSfor, objItems, selected, negate=0)
        self.connect(action, QtCore.SIGNAL('triggered()'), cb)
        
        #action = submenu.addAction(self.tr('coarse surface'))
        #cb = CallbackFunction(self.displayCoarseMSfor, objItems, selected, negate=0)
        #self.connect(action, QtCore.SIGNAL('triggered()'), cb)
        
        #submenu.addSeparator()
        #action = submenu.addAction(self.tr('Advanced...'))
        #cb = CallbackFunction(self.getGeomParams, selected, objItems)
        #self.connect(action, QtCore.SIGNAL('triggered()'), cb)
        
        # add submenu to remove representation
        submenu = menu.addMenu(self.tr('Remove Representation'))

        #action = submenu.addAction(self.tr('None'))
        for gname in geomNames:
            geom = geomsDict[gname]
            mol = geom.mol()
            action = submenu.addAction(self.tr(gname))
            #print "REMOVE MENU:" , gname
            if gname == 'atomLabels':
                cb = CallbackFunction(self.pmv.unlabelAtoms.guiCall, selected)
            elif gname == 'residueLabels':
                cb = CallbackFunction(self.pmv.unlabelResidues.guiCall, selected)
            elif hasattr(geom, '_boundGeomType'):
                cb = CallbackFunction(self.pmv.displayBoundGeom.guiCall, selected, gname, negate=True)
            elif gname.startswith("msms_"):
                cb = CallbackFunction(self.pmv.undisplayMSMS.guiCall, selected, surfNames=[gname.split("msms_")[1]])

            #elif geom.parent.name == "cartoon": # (ribbon geom)
            #    for sel in selected:
            #        if sel.getAtomGroup().getMolecule() == mol:
            #            cb = CallbackFunction(self.pmv.displayCartoon.guiCall, sel, negate=True)
            else:
                cb = CallbackFunction(self.displayfor, selected, gname, negate=True)
            self.connect(action, QtCore.SIGNAL('triggered()'), cb)

        # add submenu to Customize representation
        action = menu.addAction("Customize Representations ...")
        cb = CallbackFunction(self.customizeRepr, objItems, selected)
        self.connect(action, QtCore.SIGNAL('triggered()'), cb)

        # add submenu for labeling
        submenu = menu.addMenu(self.tr('Label'))
        ## if level >= 3:
        ##     action = submenu.addAction(self.tr('molecule'))
        ##     cb = CallbackFunction(self.pmv.labelMolecules.guiCall, selected)
        ##     self.connect(action, QtCore.SIGNAL('triggered()'), cb)
        ## if level >= 2:
        ##     action = submenu.addAction(self.tr('chains'))
        ##     cb = CallbackFunction(self.pmv.labelChains.guiCall, selected)
        ##     self.connect(action, QtCore.SIGNAL('triggered()'), cb)
        if level >= 1:
            action = submenu.addAction(self.tr('residues'))
            cb = CallbackFunction(self.pmv.labelResidues.guiCall, selected)
            self.connect(action, QtCore.SIGNAL('triggered()'), cb)

        # label atoms
        subsubmenu = submenu.addMenu(self.tr('Atoms'))
        for prop in ['name', 'element', 'atomicNumber', 'charge', 'atomType', 'hbType', 'serial', 'index']:
            action = subsubmenu.addAction(self.tr(prop))
            cb = CallbackFunction(self.pmv.labelAtoms.guiCall, selected, propName=prop)
            self.connect(action, QtCore.SIGNAL('triggered()'), cb)

        # add AD_element if the root of picked tree entry is a molecule
        if item0 and isinstance(item0._pmvObject, Molecule):
            if item0._pmvObject._ag._data.get('AD_element', None) is not None:
                action = subsubmenu.addAction(self.tr('AD_element'))
                cb = CallbackFunction(self.pmv.labelAtoms.guiCall, selected, propName='AD_element')
                self.connect(action, QtCore.SIGNAL('triggered()'), cb)
                
        if parentMenu is None:
            menu.popup(QtGui.QCursor.pos())
        #print 'TIME TO BUILD MENU', time()-t0
        return menu

    def customizeRepr(self, objItems, selSet):
        from .advancedParamGui import CustomizeRepresentationsDialog
        #import pdb; pdb.set_trace()
        molNames = numpy.unique([x.getAtomGroup().getMolecule().name for x in selSet])
        selNames = [x[0].name for x in objItems]
        if len(selNames) > 1:
            selName = "%s..%s"%(selNames[0], selNames[-1])
        else:
            selName = selNames[0]
        title = "Customize representations for " + "%s "*len(molNames)%tuple(molNames)
        ##     rootNames = numpy.unique([self.pmv.gui().objTree.getRootItem(x[1]).text(0) for x in objItems])
        ##     selName = ' '.join(rootNames)
        ##     title = "Customize representations"
        if not self.customizeReprDialog:
            self.customizeReprDialog = dialog = CustomizeRepresentationsDialog(selSet, selName, self.pmv)#, parent=self)
            dialog.setWindowModality(QtCore.Qt.NonModal)
        else:
            self.customizeReprDialog.updateGUI(selSet, selName)
        self.customizeReprDialog.setWindowTitle(title)
        self.customizeReprDialog.show()
        
    def _showmenuStatusTip(self, action):
        tip = action.statusTip()
        if self.statusBar:
            self.statusBar.showMessage(tip)#self.tr(tip))

    def select(self, selection, what, negate):
        # perform a selection.select(waht) and pmv.select the result
        sel = selection.select(what)
        if sel:
            self.pmv.select.guiCall(sel, negate=negate)
       
        
    def color(self, selections, gnames, mode):
            
        #for selection in selections:
        selection = selections
        if mode=='atomType':
            self.pmv.colorByAtomType.guiCall(selection, geomsToColor=gnames, carbonsOnly=self.COnly.isChecked())
        elif mode=='rasmol':
            self.pmv.colorByResidueType.guiCall(selection, geomsToColor=gnames,
                                        carbonsOnly=self.COnly.isChecked())
        elif mode=='shapely':
            self.pmv.colorResiduesUsingShapely.guiCall(selection, geomsToColor=gnames,
                                               carbonsOnly=self.COnly.isChecked())
        elif mode=='polarity':
            self.pmv.colorAtomsUsingDG.guiCall(selection, geomsToColor=gnames)
        elif mode=='molecule':
            self.pmv.colorByMolecules.guiCall(selection, geomsToColor=gnames,
                                      carbonsOnly=self.COnly.isChecked())
        elif mode=='chain':
            self.pmv.colorByChains.guiCall(selection, geomsToColor=gnames,
                                   carbonsOnly=self.COnly.isChecked())
        elif mode=='rainbow':
            self.pmv.colorRainbow.guiCall(selection, geomsToColor=gnames,
                                  carbonsOnly=self.COnly.isChecked())
        elif mode=='rainbowChain':
            self.pmv.colorRainbowChain.guiCall(selection, geomsToColor=gnames,
                                         carbonsOnly=self.COnly.isChecked())
        elif mode=='secstructure':
            self.pmv.colorBySS.guiCall(selection, geomsToColor=gnames)

        elif mode=='customColor':
            self.undoColorCmd = None # will hold the color command and its arguments needed
                                       # to be executed to go back to the current geometry color
                                       # if "Cancel" button of the color dialog is pressed.
            colorDialog = QtGui.QColorDialog(QtCore.Qt.green, self)
            cb = CallbackFunction(self.setCustomColor_cb, selection, gnames=gnames, carbonsOnly=self.COnly.isChecked() )
            colorDialog.currentColorChanged.connect(cb)
            colorDialog.finished.connect(self.customColor_cb)
            cb = CallbackFunction(self.setCustomColorOK_cb, selection, gnames=gnames, carbonsOnly=self.COnly.isChecked())
            colorDialog.colorSelected.connect(cb)
            colorDialog.open()

        elif mode == 'moleculesRainbow':
            self.pmv.moleculesRainbow.guiCall(selection, geomsToColor=gnames,
                                  carbonsOnly=self.COnly.isChecked())

    def setCustomColor_cb(self, selection, color, gnames=[], carbonsOnly=False):
        """Callback of the 'Custom color' dialog. Is invoked when a new color is selected."""
        if color.isValid():
            custColor = numpy.array([color.getRgbF()[:3]])
            #print custColor
            self.pmv.customColor(selection, custColor, gnames, carbonsOnly=carbonsOnly)
            if self.undoColorCmd is None:
                if len(self.pmv._topUndoCmds):
                    self.undoColorCmd = self.pmv._topUndoCmds[:]

    def customColor_cb(self, res):
        """Callback of the 'Custom color' dialog. Is called when 'Ok' or 'Cancel' button is pressed."""
        if res == 0:  # Cancel button
            if self.undoColorCmd:
                # undo the Color Dialog changes, restore the original color of the geometry.
                for item in self.undoColorCmd:
                    cmd, args, kw = item
                    if isinstance(cmd, self.pmv.color.__class__):
                        kw['setupUndo'] = False
                        cmd(*args, **kw)
                self.undoColorCmd = None

    def setCustomColorOK_cb(self, selection, color, gnames=[], carbonsOnly=False, ):
        """Callback of the 'Custom color' dialog. Is called when 'OK' button is pressed."""
        if self.undoColorCmd:
            self.pmv.undo.addUndoCall( *(self.undoColorCmd, self.pmv.customColor.name) )
            self.undoColorCmd = None
        else:
            if color.isValid():
                custColor = numpy.array([color.getRgbF()[:3]])
                self.pmv.customColor(selection, custColor, gnames, carbonsOnly=carbonsOnly)

    def displayfor(self, obj, geomName, negate=0, **kw):
        if isinstance(obj, Group):
            objects = obj.findType(Protein)
        else:
            objects = [obj]

        for obj in objects:
            if geomName=='lines':
                self.pmv.undisplayLines.guiCall(obj) if negate else self.pmv.displayLines.guiCall(obj)
            elif geomName=='sb':
                self.pmv.undisplaySB.guiCall(obj) if negate else self.pmv.displaySB.guiCall(obj)
            elif geomName=='cpk':
                self.pmv.undisplayCPK.guiCall(obj) if negate else self.pmv.displayCPK.guiCall(obj)
            elif geomName=='hpballs':
                self.pmv.displayHyperBalls.guiCall(obj, negate=negate, **kw)
                self.Redraw()
            elif geomName.startswith("cartoon"):
                for sel in obj:
                    mol = sel.getAtomGroup().getMolecule()
                    g = mol.geomContainer.geoms.get('cartoon', None)
                    if g is None: # not cartoon geometry means no SS was computed so far
                        self.pmv.computeCartoon.guiCall(obj)
                    self.pmv.undisplayCartoon.guiCall(sel) if negate else self.pmv.displayCartoon.guiCall(sel)

    def _displayMSMS(self, atoms, mol, surfName, negate=None, closedSurface=True):
        # here all atoms are in mol
        app = self.pmv
        gc = mol.geomContainer
        if app.computeMSMS.isLoader():
            app.computeMSMS.loadCommand()
        cidx = mol._ag.getACSIndex()
        if not hasattr(mol, '_msmsData') or mol._msmsData['msms'].get(surfName, None) is None:
            app.computeMSMS.guiCall(atoms, surfName=surfName, closedSurface=closedSurface)
            app.displayMSMS.guiCall(atoms, surfNames=[surfName])
        elif mol._msmsData['msms'][surfName][cidx] is None:
            app.computeMSMS.guiCall(atoms, surfName=surfName, closedSurface=closedSurface)
            app.displayMSMS.guiCall(atoms, surfNames=[surfName])
        else:
            srf = mol._msmsData['msms'][surfName][cidx]
            srfAtoms = mol._msmsData['atoms'][surfName]
            if negate is None:
                negate = gc.displayedAs([surfName], atoms, 'fast')
            if negate:
                app.undisplayMSMS.guiCall(atoms, surfNames=[surfName])
            else:
                app.displayMSMS.guiCall(atoms, surfNames=[surfName])
        
    def displayMSMSfor(self, objItems, atoms, negate=None):
        """Display (and compute if needed) molecular surfaces.
        """
        app = self.pmv
        obj0, item0 = objItems[0]
        rootItem = item0
        if rootItem: # rootItem can be None when we right click on atom in camera
            while rootItem.parent() is not None:
                rootItem = rootItem.parent()
            molFrag = rootItem._pmvObject
        else:
            molFrag = atoms
        #mols, atms = self.pmv.getNodesByMolecule(atoms)

        # for Selections, there can only be one selection we are operating on
        #import pdb; pdb.set_trace()
        if isinstance(molFrag, PmvSelection):
            # if it is the current selection we operate on 'surface-%s'%mol.name
            if molFrag is self.pmv.curSelection:
                for ats in atoms:
                    mol = ats.getAtomGroup().getMolecule()
                    surfName = '%s_surface'%mol._basename
                    if not hasattr(mol, '_msmsData') or mol._msmsData['msms'].get(surfName, None) is None:
                        self.pmv.computeMSMS(mol.select(), surfName=surfName)
                    if negate:
                        self.pmv.undisplayMSMS(ats, surfNames=[surfName])
                    else:
                        self.pmv.displayMSMS(ats, surfNames=[surfName]) 
            else:
                surfName = '%s_surface'%molFrag.name
                for ats in atoms:
                    mol = ats.getAtomGroup().getMolecule()
                    if not hasattr(mol, '_msmsData') or mol._msmsData['msms'].get(surfName, None) is None:
                        self.pmv.computeMSMS(ats, surfName=surfName)
                    if negate:
                        self.pmv.undisplayMSMS(ats, surfNames=[surfName])
                    else:
                        self.pmv.displayMSMS(ats, surfNames=[surfName]) 
        else:
            # if it is the current selection we operate on '%s_surface'%mol.name
            #for mol, ats in zip(mols, atms):
            for ats in atoms:
                mol = ats.getAtomGroup().getMolecule()
                surfName = '%s_surface'%mol._basename
                if not hasattr(mol, '_msmsData') or mol._msmsData['msms'].get(surfName, None) is None:
                    self.pmv.computeMSMS(mol, surfName=surfName)
                if negate:
                    self.pmv.undisplayMSMS(ats, surfNames=[surfName])
                else:
                    self.pmv.displayMSMS(ats, surfNames=[surfName]) 

        ## if isinstance(molFrag, PmvSelection):
        ##     # if it is the current selection we operate on 'surface-%s'%mol.name
        ##     if molFrag is self.pmv.curSelection:
        ##         #for mol, ats in zip(mols, atms):
        ##         for ats in atoms:
        ##             mol = ats.getAtomGroup().getMolecule()
        ##             self._displayMSMS(mol.select(), mol, '%s_surface'%mol._basename, negate=negate, closedSurface=False)
        ##     else:
        ##         #for mol, ats in zip(mols, atms):
        ##         #    self._displayMSMS(ats, mol, 'surface-%s'%sele.name, negate=negate, perMol=False)
        ##         for ats in atoms:
        ##             mol = ats.getAtomGroup().getMolecule()
        ##             self._displayMSMS(ats, mol, 'surface-%s'%molFrag.name, negate=negate, closedSurface=True)
        ## else:
        ##     # if it is the current selection we operate on '%s_surface'%mol.name
        ##     #for mol, ats in zip(mols, atms):
        ##     for ats in atoms:
        ##         mol = ats.getAtomGroup().getMolecule()
        ##         surfName = '%s_surface'%mol._basename

        ##         self._displayMSMS(ats, mol, surfName, negate=negate, closedSurface=False)
        
        

    def displayCoarseMSfor(self, objItems, atoms, negate=None):
        """Display coarse molecular surfaces.
        """
        app = self.pmv
        obj0, item0 = objItems[0]
        rootItem = item0
        if rootItem: # rootItem can be None when we right click on atom in camera
            while rootItem.parent() is not None:
                rootItem = rootItem.parent()
            molFrag = rootItem._pmvObject
        else:
            molFrag = atoms
        # for Selections, there can only be one selection we are operating on

        if isinstance(molFrag, PmvSelection):
            # if it is the current selection we operate on 'surface-%s'%mol.name
            for ats in atoms:
                recompute = False
                mol = ats.getAtomGroup().getMolecule()
                if molFrag is self.pmv.curSelection:
                    surfName = "coarseMolSurf-%s"%mol.name
                    selection = mol
                else:
                    surfName = "coarseMolSurf-%s"%molFrag.name
                    selection = ats
                gc = mol.geomContainer
                if not gc.geoms.has_key(surfName) or len(ats)!=len(gc.atoms[surfName]):
                    recompute = True
                if recompute:
                    app.computeCoarseMolecularSurface.guiCall(selection, surfName=surfName, bind_surface_to_molecule=True)
                    negate = False
                else:
                    if negate is None:
                        negate = gc.displayedAs([surfName], ats, 'fast')
                app.displayBoundGeom.guiCall(ats, surfName, negate=negate)
        else: # molecule
            for ats in atoms:
                mol = ats.getAtomGroup().getMolecule()
                gc = mol.geomContainer
                surfName = "coarseMolSurf-%s"%mol.name
                if not gc.geoms.has_key(surfName):
                   self.pmv.computeCoarseMolecularSurface.guiCall(mol, surfName=surfName, bind_surface_to_molecule=True)
                   negate = False
                if negate is None:
                    negate = gc.displayedAs([surfName], ats, 'fast')   
                app.displayBoundGeom.guiCall(ats, surfName, negate=negate)

    def setCutOff(self, cmd, selection, obj):
        """This sets cutoff for Expand Selection and Select Around commands.
        Expand selection cmd will select atoms located within a user-specified distance
        from the currently selected atoms (selection). Select around will do the same as
        expand selection and will also de-select the current selection. The selection expands
        only to atoms belonging to the molecular fragment (obj) on which the selection menu was
        activated.
        As the cutoff is modified, cyan colored spheres show what would be selected for the current
        cutoff value. """
        app = self.pmv
        if hasattr(cmd , "loadCommand"):
            cmd = cmd.loadCommand()
        centers = selection[0].getCoords()
        if len(selection) > 1:
            for sele in selection[1:]:
                centers = numpy.concatenate((centers, sele.getCoords() ))
        # add geometry to show what would be selected with this cut off
        from DejaVu2.Spheres import Spheres
        self.showCutoffSph = Spheres(
            'cutOffFeedBack', inheritMaterial=0, radii=(0.3,),
            materials=((0,1,1, 0.5),), transparent=1)

        #self.pmv.gui().viewer.AddObject(self.showCutoffSph)
        self.AddObject(self.showCutoffSph)
        
        def showCutOffSelection(dist):
            # show what this cutoff will select
            if dist <= 0: return 
            atsets = cmd.getAtoms(centers, dist, obj)
            if len(atsets):
                coords = atsets[0].getCoords()
                if len(atsets) > 1:
                    for ats in atsets[1:]:
                        coords = numpy.concatenate((coords, ats.getCoords() ))
                self.showCutoffSph.Set(visible=True, vertices=coords)
            else:
                self.showCutoffSph.Set(visible=False)
                
        from .setDistanceDialog import SetDistanceDialog
        if hasattr(obj, "name") and obj.name != "NoName":
            title="%s for %s"% (cmd.name, obj.name)
        elif hasattr(cmd, "name"):
            title = cmd.name
        else:
            title = 'select around'
        opt, ok = SetDistanceDialog.getDialogParams(title=title, callBack=showCutOffSelection)
        self.pmv.gui().viewer.RemoveObject(self.showCutoffSph)
        if ok:    
            dist = opt['distance']
            if dist > 0:
                cmd(selection, dist, obj)

    def registerListerners(self):
        ##
        evh = self.pmv.eventHandler

        from AppFramework.App import AddGeometryEvent, RemoveGeometryEvent, \
             RedrawEvent
        evh.registerListener(AddGeometryEvent, self.addGeometryEventHandler)
        evh.registerListener(RemoveGeometryEvent,
                             self.removeGeometryEventHandler)
        evh.registerListener(RedrawEvent, self.redrawEventHandler)

        ##
        ## Selection Events        
        from PmvApp.selectionCmds import SelectionEvent, RefreshSelectionEvent
        evh.registerListener(RefreshSelectionEvent, self.selectionEventHandler)
        evh.registerListener(SelectionEvent, self.selectionEventHandler)
        from PmvApp.Pmv import EditGeomsEvent
        evh.registerListener(EditGeomsEvent, self.highlightSelection)

    def addGeometryEventHandler(self, event):
        #obj, parent=None, redo=False):
        #print 'Handling add geom event for', event.object.name
        self.AddObject(event.object, parent=event.parent, redo=event.redo)
        
    def removeGeometryEventHandler(self, event):
        #obj, parent=None, redo=False):
        #print 'Handling add geom event for', event.object.name
        self.RemoveObject(event.object)

    def redrawEventHandler(self, event):
        #print 'Handling redraw event'
        self.redraw()

    def selectionEventHandler(self, event=None):
        self.updateSelectionIcons(event)
        self.highlightSelection(event)

    def updateSelectionIcons(self, event=None):
        """update selection icons"""
        #if SelectionSpheres is currently turned on: 
        app = self.pmv
        l = app.activeSelection.nbAtoms()
        if app.gui is not None and hasattr(app.gui(), "activeSelection") and app.gui().activeSelection is None or l==0:
                self.selectionCrosses.Set(visible=0, tagModified=False)
        else:
            coords = numpy.zeros( (l, 3), 'f')
            off = 0
            for sel in app.activeSelection:
                l1 = len(sel)
                coords[off:off+l1] = sel.getCoords()
                off += l1
            self.selectionCrosses.Set(vertices=coords, visible=1, tagModified=False)
        self.Redraw()

    def highlightSelection(self, event=None):
        #print "highlight selection"
        app = self.pmv

        ## if self.activeSelection is None or len(self.pmv.activeSelection)==0:
        ##     selMols = selAtms = []
        ## else:
        ##     selMols, selAtms = app.getNodesByMolecule(self.pmv.activeSelection.get())
        ## allMols = set( app.Mols[:] )
        ## unselectedMols = allMols.difference(selMols)
        #for sel in event.setOff:
        #    mol = sel.getAtomGroup().getMolecule()
        ## reset highlight 
        for mol in app.Mols:
            geomC = mol.geomContainer
            for geomName, lGeom in geomC.geoms.items():
                if isinstance(lGeom, Spheres) \
                  or isinstance(lGeom, Cylinders)\
                  or (isinstance(lGeom, IndexedPolygons) \
                      and hasattr(lGeom,'_isMSMS') ) \
                   or (isinstance(lGeom, IndexedPolygons) \
                       and hasattr(geomC, "boundGeom") and \
                       geomC.boundGeom.has_key(geomName)):
                    lGeom.Set(highlight=[])
                elif isinstance(lGeom, IndexedPolygons) \
                  and lGeom.parent.name == 'secondarystructure':
                    lGeom.Set(highlight=[])

        ## if self.activeSelection is None or len(self.pmv.activeSelection)==0:
        ##     selMols2 = selResidueSets = []
        ## else:
        ##     selMols2, selResidueSets = app.getNodesByMolecule(self.pmv.activeSelection.get(Residue))
        ## molSelectedResiduesDict = dict( zip( selMols2, selResidueSets) )
        ## for mol, selectedResidueSet in molSelectedResiduesDict.items():
        ##     for lGeom in mol.geomContainer.geoms.values():
        ##         if isinstance(lGeom, IndexedPolygons) and lGeom.parent.name == 'secondarystructure':
        ##             highlight = [0] * len(lGeom.vertexSet.vertices.array)
        ##             for selectedResidue in selectedResidueSet:
        ##                 if hasattr(lGeom, 'resfacesDict') and lGeom.resfacesDict.has_key(selectedResidue):
        ##                     for lFace in lGeom.resfacesDict[selectedResidue]:
        ##                         for lVertexIndex in lFace:
        ##                             highlight[int(lVertexIndex)] = 1
        ##             lGeom.Set(highlight=highlight)

        ## set highlight 
        for sel in app.activeSelection:
            mol = sel.getAtomGroup().getMolecule()
            geomC = mol.geomContainer
            for geomName, lGeom in geomC.geoms.items():
                if isinstance(lGeom, Spheres) \
                  or isinstance(lGeom, Cylinders):
                    highlight = numpy.zeros( (len(lGeom.vertexSet),), 'i')
                    if len(highlight):
                        indices = sel.getIndices()
                        bonds = None
                        if geomName == "doubleBondsSticks":
                            bonds = sel.getBonds()
                            doubleBonds = mol._bondOrderData['doubleBonds']
                            for i1,i2 in bonds[2]:
                                k = "%d %d" %(i1, i2)
                                if doubleBonds.has_key(k):
                                    n = doubleBonds[k] * 4
                                    highlight[[n,n+1,n+2,n+3]] = 1
                        elif geomName == 'tripleBondsSticks':
                            if not bonds: bonds = sel.getBonds()
                            tripleBonds = mol._bondOrderData['tripleBonds']
                            for i1,i2 in bonds[3]:
                                k = "%d %d" %(i1, i2)
                                if tripleBonds.has_key(k):
                                    n = tripleBonds[k] * 4
                                    highlight[[n, n+1, n+2, n+3]] = 1 
                        elif geomName == 'aromaticSpheres':
                            aromaticArcs = mol._bondOrderData['aromaticArcs']
                            if aromaticArcs is not None and len(aromaticArcs):
                                aindices = aromaticArcs.keys()
                                arinds = numpy.intersect1d(indices, aindices)
                                for index in arinds:
                                    if aromaticArcs.has_key(index):
                                        for ii in aromaticArcs[index]:
                                            highlight[ii[0]] = 1
                        else:
                            highlight[indices] = 1
                        lGeom.Set(highlight=highlight.tolist())

                elif isinstance(lGeom, IndexedPolygons):
                    if hasattr(lGeom,'_isMSMS'):
                        lAtomSet = geomC.atoms[geomName]
                        if len(lAtomSet) > 0:
                            highlight = numpy.zeros(len(lGeom.vertexSet.vertices), 'i')
                            if len(highlight):
                                name = geomName[5:] # skip "msms_"
                                surfinds = mol._msmsData['indexMaps'][name]
                                al = sel.select('not deleted')
                                indices = al.getIndices()
                                atomindices = surfinds[indices]
                                srf = mol._msmsData['msms'][name][mol._ag.getACSIndex()]
                                globProp = mol._renderingProp['msms'][name]
                                nbvInComp = lGeom._nbvInComp

                                highlight = numpy.zeros( (len(lGeom.vertexSet),), 'i')
                                fs = []
                                if globProp['components']=='all':
                                    compIndices = range(srf.rsr.nb)
                                else:
                                    compIndices = globProp['components']
                                for cn, compNum in enumerate(compIndices):
                                    dum, dum, fsc = srf.getTriangles(atomindices, component=compNum,
                                                                     selnum=3, keepOriginalIndices=1)
                                    if len(fsc):
                                        if cn>0:
                                            fsc = fsc+nbvInComp[cn-1]
                                        fs.extend(fsc.tolist())
                                        fs = numpy.array(fs)
                                        highlight[numpy.unique(fs.flatten())] = 1

                                #print 'MINI3', min(atomindices)
                                ## lvf, lvint, lTri = srf.getTriangles(
                                ##     atomindices,selnum=1, keepOriginalIndices=1)
                                ## for lThreeIndices in lTri:
                                ##     highlight[int(lThreeIndices[0])] = 1
                                ##     highlight[int(lThreeIndices[1])] = 1
                                ##     highlight[int(lThreeIndices[2])] = 1
                                        lGeom.Set(highlight=highlight)
                    
                    elif hasattr(geomC, "boundGeom") and geomC.boundGeom.has_key(geomName):
                        lSelectedAtoms = sel.getIndices()
                        if len(lSelectedAtoms):
                            lAtomsSet = geomC.boundGeom[geomName]['atoms']
                            cl_atoms = geomC.boundGeom[geomName]['cl_atoms']
                            highlight = numpy.zeros(len(lGeom.vertexSet.vertices), 'i')
                            mol_lookup = geomC.boundGeom[geomName]['mol_lookup']
                            lAtomVerticesDict = {}
                            for i, at  in enumerate(cl_atoms):
                                if not lAtomVerticesDict.has_key(at):
                                    lAtomVerticesDict[at]=[]
                                lAtomVerticesDict[at].append(i)
                            for at in lSelectedAtoms:
                                _at = mol_lookup.get(at, None)
                                if _at is not None:
                                    vertInds = lAtomVerticesDict.get(_at, [])
                                    for i in vertInds:
                                        highlight[i] = 1
                            lGeom.Set(highlight=highlight)
            
        ## for mol, atoms in map(None, selMols, selAtms):
        ##     for geomName, lGeom in mol.geomContainer.geoms.items():
        ##         if   isinstance(lGeom, Spheres) \
        ##           or isinstance(lGeom, Cylinders):
        ##             lAtomSet = mol.geomContainer.atoms[geomName]
        ##             if len(lAtomSet) > 0:
        ##                 lAtomSetDict = dict(zip(lAtomSet, range(len(lAtomSet))))
        ##                 highlight = [0] * len(lAtomSet)
        ##                 for i in range(len(atoms)):
        ##                     lIndex = lAtomSetDict.get(atoms[i], None)
        ##                     if lIndex is not None:
        ##                         highlight[lIndex] = 1
        ##                 lGeom.Set(highlight=highlight)
        ##         elif isinstance(lGeom, IndexedPolygons):
        ##           if hasattr(mol.geomContainer,'msmsAtoms') and mol.geomContainer.msmsAtoms.has_key(geomName):
        ##             lAtomSet = mol.geomContainer.msmsAtoms[geomName]
        ##             if len(lAtomSet) > 0:
        ##                 lAtomSetDict = dict(zip(lAtomSet, range(len(lAtomSet))))
        ##                 lAtomIndices = []
        ##                 for i in range(len(atoms)):
        ##                     lIndex = lAtomSetDict.get(atoms[i], None)
        ##                     if lIndex is not None:
        ##                         lAtomIndices.append(lIndex)
        ##                 lSrfMsms = mol.geomContainer.msms[geomName][0]
        ##                 lvf, lvint, lTri = lSrfMsms.getTriangles(
        ##                     lAtomIndices, 
        ##                     selnum=numOfSelectedVerticesToSelectTriangle,
        ##                     keepOriginalIndices=1)
        ##                 highlight = [0] * len(lGeom.vertexSet.vertices)
        ##                 for lThreeIndices in lTri:
        ##                     highlight[int(lThreeIndices[0])] = 1
        ##                     highlight[int(lThreeIndices[1])] = 1
        ##                     highlight[int(lThreeIndices[2])] = 1
        ##                 lGeom.Set(highlight=highlight)
        ##           elif app.bindGeomToMolecularFragment.data.has_key(lGeom.fullName) \
        ##             and app.bindGeomToMolecularFragment.data[lGeom.fullName].has_key('atomVertices'):
        ##               bindcmd = app.bindGeomToMolecularFragment
        ##               lSelectedAtoms = atoms
        ##               if len(lSelectedAtoms) > 0:
        ##                 lAtomVerticesDict = bindcmd.data[lGeom.fullName]['atomVertices']
        ##                 highlight = [0] * len(lGeom.vertexSet.vertices)
        ##                 for lSelectedAtom in lSelectedAtoms:
        ##                     lVertexIndices = lAtomVerticesDict.get(lSelectedAtom, [])
        ##                     for lVertexIndex in lVertexIndices:
        ##                         highlight[lVertexIndex] = 1
        ##                 lGeom.Set(highlight=highlight)

class PmvGUI(QtGui.QMainWindow):
    
    _NoERROR = 0
    _WARNING = 1
    _ERROR = 2
    _EXCEPTION = 3

    windowShown = QtCore.Signal()

    def closeEvent(self, event):
        alive = self.pmv.checkForThreads()
        if len(alive):
            msg = "the following threads are still running:\n"
            for descr in alive: msg += "    %s\n"%descr
            msg +="\n\nPmv will exit after these threads finish or you hit Ctrl-C int eh shell"
            QtGui.QMessageBox.information(
                self.pmv.gui(), "Pmv has running threads", msg,
                QtGui.QMessageBox.Ok)
        ## #import pdb; pdb.set_trace()
        ## if reply == QtGui.QMessageBox.No:
        ##     event.ignore()
        ##     return
        self.pmv.cleanup()
        self.settings = QtCore.QSettings("TSRI", "PmvApp")
        self.settings.beginGroup("MainWindow")
        self.settings.setValue("Pmv/geometry", self.saveGeometry())
        self.settings.setValue("Pmv/windowState", self.saveState())
        self.settings.endGroup()
        QtGui.QMainWindow.closeEvent(self, event)
        
    def readSettings(self):
        self.settings = QtCore.QSettings("TSRI", "PmvApp")
        self.settings.beginGroup("MainWindow")
        self.restoreGeometry(self.settings.value("Pmv/geometry"))
        self.restoreState(self.settings.value("Pmv/windowState"))
        self.settings.endGroup()

    def event(self, event):
        ret_val = QtGui.QMainWindow.event(self, event)
        if not self._functionAfterShownCalled and event.type() == QtCore.QEvent.Paint:
            self.windowShown.emit()
            self._functionAfterShownCalled = True
        return ret_val
    
    def __init__(self, pmv, parent=None, classCamera=None, autoRedraw=True):
 	"""PmvGui constructor
"""
        self.pmv = pmv
        pmv.gui = weakref.ref(self)
        self.settings = None
        self._functionAfterShownCalled = False
        
        # if None no selection is active in the GUI
        # else is points to the active Pmv Selection object
        self.activeSelection = None
        
        QtGui.QMainWindow.__init__(self, parent)
        sb = self.statusBar()
        #b = self.execStatusWidget = QtGui.QToolButton()
        b = self.execStatusWidget = QtGui.QPushButton()
        b.setStyleSheet("QPushButton { border: none; background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1, stop: 0 #a6a6a6, stop: 0.08 #7f7f7f, stop: 0.39999 #717171, stop: 0.4 #626262, stop: 0.9 #4c4c4c, stop: 1 #333333);}")
        #b.setContentsMargins(0,0,0,0)
        #b.setFlat(True)
        #b.setIcon(QtGui.QIcon(os.path.join(PMVICONPATH, 'logWindow.png')))
        #b.setIconSize(QtCore.QSize(16,16))
        #b.resize(16, 16)
        sb.addPermanentWidget ( self.execStatusWidget, stretch = 0 )
        self.setFont(contentFont())
        self.groups = {}
        import mglutil
        from Support import version
        self.setWindowTitle('PMV2_Qt %s (built: %s)'%(
            version.__version__, mglutil.__revision__))
        self.createDockedtWidgets()
        self.viewer.currentCamera.keyPressSignal.connect(self.keyPressed_cb)
        self.createActions()
        b.released.connect(self.displayReports)
        self.createMenus()
        self.createToolBar()
        self.createStatusBar()
        self.registerListerners()
        # dont show here because on mac event() will be called and emt windowShow too early
        #self.show()

        self.unseenReports = [] # list of execution reports not yet seen by the user
        self._worstError = self._NoERROR
        self.readSettings()
        self.trgMapGui = None
        
    def createStatusBar(self):
        self.statusBar().showMessage(self.tr("Ready"))

    def createActions(self):
        # open file action
        open = self.openAct = QtGui.QAction(
            QtGui.QIcon(os.path.join(PMVICONPATH, 'open.png')),
            'Open', self)
        open.setShortcut('Ctrl+O')
        open.setStatusTip('Open File')
        self.connect(open, QtCore.SIGNAL('triggered()'), self.viewer.openFile)

        # save session
        saveSession = self.saveSessionAct =  QtGui.QAction(
            QtGui.QIcon(os.path.join(PMVICONPATH, 'save.png')),
            'Save session', self)
        saveSession.setShortcut('Ctrl+S')
        saveSession.setStatusTip('Save session')
        self.connect(saveSession,  QtCore.SIGNAL('triggered()'), self.viewer.saveSession)


        # undo Action
        undo = self.undoAct = QtGui.QAction(
            QtGui.QIcon(os.path.join(PMVICONPATH, 'undo.png')),
            'Undo', self)
        undo.setShortcut('Ctrl+Z')
        undo.setStatusTip('Undo last command')
        self.connect(undo, QtCore.SIGNAL('triggered()'), self.pmv.undo)

        # redo Action
        redo = self.redoAct = QtGui.QAction(
            QtGui.QIcon(os.path.join(PMVICONPATH, 'redo.png')),
           'Redo', self)
        redo.setShortcut('Ctrl+Shift+Z')
        redo.setStatusTip('Redo last command')
        self.connect(redo, QtCore.SIGNAL('triggered()'), self.pmv.redo)

        # clearSelection Action
        clearSel = self.clearSelAct = QtGui.QAction(
            QtGui.QIcon(os.path.join(PMVICONPATH, 'eraser.png')),
            'Clear selection', self)
        clearSel.setStatusTip('Clear Selection')
        self.connect(clearSel, QtCore.SIGNAL('triggered()'), self.pmv.clearSelection)

        # focus 3D Scene
        focusScene = self.focusSceneAct = QtGui.QAction(
            QtGui.QIcon(os.path.join(PMVICONPATH, 'crosshair.png')),
            'Normalize and center 3D Scene', self)
        focusScene.setStatusTip('Normalize and center 3D Scene using current selection')
        self.connect(focusScene, QtCore.SIGNAL('triggered()'), self.pmv.focusScene)

        # report action
        report = self.reportAct = QtGui.QAction(
            QtGui.QIcon(os.path.join(PMVICONPATH, 'report.png')),
            'Execution Report', self)
        report.setStatusTip('Display Execution Report')
        self.connect(report, QtCore.SIGNAL('triggered()'), self.displayReports)

        # Python shell
        pyShell = self.togglePyShellAct =  QtGui.QAction(
            QtGui.QIcon(os.path.join(PMVICONPATH, 'PyShell.png')),
            'Python shell', self)
        pyShell.setStatusTip('Show/hide Python shell')
        self.connect(pyShell, QtCore.SIGNAL('triggered()'), self.showPythonShell)

        #snapshots= self.snapshots =  QtGui.QAction(
        #    QtGui.QIcon(os.path.join(PMVICONPATH, 'snapshot.png')),
        #    'Snapshots', self)
        #self.connect(snapshots, QtCore.SIGNAL('triggered()'), self.showSnapshotsGUI)      
        
        # exit action
        exit = self.exitAct = QtGui.QAction(
            QtGui.QIcon(os.path.join(PMVICONPATH, 'quit.png')),
            'Exit', self)
        exit.setShortcut('Ctrl+Q')
        exit.setStatusTip('Exit application')
        self.connect(exit, QtCore.SIGNAL('triggered()'), QtCore.SLOT('close()'))

        ## # select action
        ## select = self.selectAct = QtGui.QAction(QtGui.QIcon('icons/config.png'),
        ##                                       'Select ...', self)
        ## select.setShortcut('Ctrl+S')
        ## select.setStatusTip('Select')
        ## self.connect(select, QtCore.SIGNAL('triggered()'), self.selectionPanel)

    def selectionPanel(self):
        print 'DISPLAY SELECTION WINDOW'
        # windows visiblity actions
        #self.toggleLogViewAct = QtGui.QAction('showHide_Log', self)
        #connect(openAct, QtCore.SIGNAL("triggered()"),
        #        self.logWidget, SLOT("setVisible(bool)"))

    def createMenus(self):
        self.fileMenu = self.menuBar().addMenu(self.tr("&File"))
        self.fileMenu.addAction(self.openAct)
        self.fileMenu.addAction(self.saveSessionAct)
        self.fileMenu.addSeparator()

        self.fetchMenu = self.fileMenu.addAction("Fetch File")
        self.connect(self.fetchMenu, QtCore.SIGNAL('triggered()'), self.fetchFile_cb)
        
        submenu = self.fileMenu.addMenu(self.tr('Recent Files'))
        cb = CallbackFunction(self.getRecentFilesMenu, submenu)
        self.connect(submenu, QtCore.SIGNAL('aboutToShow()'), cb)
        self.fileMenu.addSeparator()

        self.showMapGui = self.fileMenu.addAction('Trg Map GUI')
        self.connect(self.showMapGui, QtCore.SIGNAL('triggered()'), self.showTrgMapGui_cb)
        self.fileMenu.addSeparator()
        
        ## mapMenu = self.mapMenu= self.fileMenu.addAction("Read maps")
        ## self.connect(mapMenu, QtCore.SIGNAL('triggered()'), self.viewer.openMapZipFile)
        ## self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.exitAct)

        self.editMenu = self.menuBar().addMenu(self.tr("&Edit"))
        self.editMenu.addAction(self.focusSceneAct)
        self.editMenu.addAction(self.undoAct)
        self.editMenu.addAction(self.redoAct)
        self.editMenu.addAction(self.clearSelAct)
        
        self.newGrpAct = self.editMenu.addAction('Create New Group')
        self.connect(self.newGrpAct, QtCore.SIGNAL('triggered()'),
                     self.createNewGroup)

        self.windowsMenu = self.menuBar().addMenu(self.tr("&Windows"))
        self.windowsMenu.addAction(self.toggleShellViewAct)
        self.windowsMenu.addAction(self.toggleLogViewAct)
#        self.windowsMenu.addAction(self.toggleCmdlineViewAct)
        self.windowsMenu.addAction(self.toggleDashboardViewAct)
        self.windowsMenu.addAction(self.toggle3DViewViewAct)
        self.windowsMenu.addAction(self.reportAct)
        
        self.helpMenu = self.menuBar().addMenu(self.tr("&Help"))
        self.aboutAct = self.helpMenu.addAction(self.tr("&About"))
        self.aboutAct.triggered.connect(self.about)

        msmsEnabled = False
        for p in sys.path:
            if os.path.exists(os.path.join(p, 'mslib', 'msms.py')):
                msmsEnabled = True
                break

        if not msmsEnabled:
            self.enableMSMSAct = self.helpMenu.addAction(
                self.tr("&Enable MSMS-based molecular surfaces"))
            self.enableMSMSAct.triggered.connect(self.enableMSMS)

    def about(self):
        #print 'ABOUT'
        msg = QtGui.QMessageBox()
        #msg.setIcon(QtGui.QMessageBox.Information)
        
        msg.setWindowTitle("About PMV2")
        from Support.version import __version__
        from mglutil import __revision__
        vtxt = "Python Molecule Viewer, version: %s(%s)." %(__version__, __revision__)
        import platform
        uname = platform.uname()
        nametxt = "%s "*len(uname) % uname
        txt = """%s \nRunning on %s.""" %(vtxt, nametxt)
        msg.setText(txt)
        #msg.setDetailedText(txt)
        msg.setStandardButtons(QtGui.QMessageBox.Ok)
        #msg.buttonClicked.connect(msgbtn)
        retval = msg.exec_()
        
    def enableMSMS(self):
        reply = QtGui.QInputDialog.getText(
            None, "Enable MSMS/MSLIB",
            "Please enter our license key.\nContact Dr. Sanner (sanner@scripps.edu) for a license key")
        if reply[1]:
            # user clicked OK
            from mglkey import MGL_check_key
            if MGL_check_key(reply[0].encode('ascii', 'replace')):
                mslibACA = None
                for p in sys.path:
                    if os.path.exists(os.path.join(p, 'mslibACA')):
                        mslibACA = os.path.join(p, 'mslibACA')
                        break
                mslib = None
                for p in sys.path:
                    if os.path.exists(os.path.join(p, 'mslib')):
                        mslib = os.path.join(p, 'mslib')
                        break
                if mslibACA is None or mslib is None:
                    reply = QtGui.QMessageBox.information(
                        None, 'ERROR locating mslib',
                        "The mslib package could not be located.\nPlease contact Dr Sanner at sanner@scripps.edu") 
                else:
                    if platform.system()=='Darwin':
                       sudoPassword = QtGui.QInputDialog().getText(
                           None, "Password", "Enter your password",
                           QtGui.QLineEdit.EchoMode.Password)
                       sudoPassword = sudoPassword[0].encode('ascii', 'replace')
                       command = 'mv %s %s'%(mslib, mslib+'COM')
                       p = os.system('echo %s|sudo -S %s' % (sudoPassword,
                                                             command))
                       command = 'mv %s %s'%(mslibACA, mslibACA[:-3])
                       p = os.system('echo %s|sudo -S %s' % (sudoPassword,
                                                             command))
                    else: # for windows and Linux we can simply rename
                        os.rename(mslib, mslib+'COM')
                        os.rename(mslibACA, mslibACA[:-3])

                    try:
                        import mslib
                        if not self.pmv.computeMSMS.isLoader():
                            self.pmv.computeMSMS.MSMSdisabled = False
                        self.enableMSMSAct.setDisabled(True)
                    except RuntimeError:
                        reply = QtGui.QMessageBox.information(
                            None, 'ERROR enabling MSMS',
                            "Failed to enable mslib.\nPlease contact Dr Sanner at sanner@scripps.edu") 
            else:
                reply = QtGui.QMessageBox.question(
                    None, "", "The key you entered is not valid.\n try again ?",
                    QtGui.QMessageBox.Yes | QtGui.QMessageBox.No,
                    QtGui.QMessageBox.No)
                if reply == QtGui.QMessageBox.Yes:
                    self.enableMSMS()

    def createNewGroup(self, name=None, parent=None):
        if name is None:
            name, ok = QtGui.QInputDialog.getText(
                self, 'Group Name', 'Group name:')
            name = name.encode('ascii', 'ignore')
        else:
            ok = True
            
        if ok and name:
            self.pmv.createGroup(name)

    def fetchFile_cb(self, name=None):
        if name is None:
            res = False
            while not res:
                name, ans = FetchGUI.getName(parent=self)
                res = True
                names = [n for n in name[0].split(" ") if len(n)]
                ext = name[1]
                #print "fetchFile:", names, ext
                if ans:
                    mols = self.pmv.fetch(names, ext)
                    #print "fetchFile_cb", mols
                    #import pdb; pdb.set_trace()
                    errmsg = ""
                    if self.pmv.trapExceptions:
                        #find out if the command raised an exception  
                        exreport =  self.pmv._executionReport
                        if exreport.cmd.name=='fetch' and exreport.numberOf['errors']>0:
                            err = exreport.getErrors()[0]
                            if err.obj in names:
                                from urllib2 import URLError, HTTPError
                                if isinstance(err.exception, HTTPError):
                                    errmsg =  'ERROR: "%s" is not a valid pdb id.\nError code %s(%s)'%(err.obj, err.exception.code, err.exception.reason)
                                elif isinstance(err.exception, URLError):
                                    errmsg = "ERROR: the server could not be reached.\nPlease check your internet connection.\n"
                                    if hasattr(err.exception, 'reason'):
                                        msg = msg+"%s" % (err.exception.reason)
                                else:
                                    errmsg = "ERROR: failed to fetch %s" % (err.obj)
                                ## exc = err.getException("")
                                ## if not len(exc):
                                ##     fexc = err.formatedException
                                ##     if len(fexc):
                                ##        for i in range(len(fexc)-1, -1, -1): 
                                ##            if len(fexc[i])>0:
                                ##                exc = [fexc[i]]
                                ##                break
                                ##     else:
                                ##        exc = [""] 
                                ## errmsg = "Fetch failed: %s %s" %(err.obj, exc[-1])
                    else:
                        #if trapExceptions is set to False, the command will not raise an exception if fetching the molecule fails. It will return an error message along with an empty list for molecules.
                        for mol, msg in mols:
                            if len(msg):
                                errmsg = errmsg+"\n" + msg
                    if len(errmsg):        
                        # pop up a warning dialog, ask the user to retry
                        msgBox = QtGui.QMessageBox(QtGui.QMessageBox.Warning,
                                                   "Fetch Failed", errmsg,
                                                   QtGui.QMessageBox.NoButton, self)
                        msgBox.addButton("Try Again", QtGui.QMessageBox.AcceptRole)
                        msgBox.addButton("Close", QtGui.QMessageBox.RejectRole)
                        if msgBox.exec_() == QtGui.QMessageBox.AcceptRole:
                            res = False
                            #continue

                    
    def getRecentFilesMenu(self, submenu):
        actions = [act.text() for act in submenu.actions()]
        for category, files in  self.pmv.recentFiles.categories.items():
            for fname, cmdname in files:
                if fname not in actions:
                    self.recentFilesAct = submenu.addAction(self.tr(fname))
                    cmd = eval("self.pmv.%s"%cmdname)
                    cb = CallbackFunction(cmd, fname)
                    self.connect(self.recentFilesAct, QtCore.SIGNAL('triggered()'), cb)

    def showTrgMapGui_cb(self):
        if not self.trgMapGui:
            from ADFR.GUI.trgMapsGUI import TrgMapsGui
            self.trgMapGui = TrgMapsGui(parent=self, app=self)
            self.trgMapGui.viewer = self.viewer
            self.trgMapGui.show()
            #self.trgMapGui.addButton.setDisabled(True)
            #self.trgMapGui.removeButton.setDisabled(True)
        else:
            # the GUI has been created. Update its tree widget with the
            if not self.trgMapGui.isVisible():
                self.trgMapGui.show()
                
    def createToolBar(self):
        self.toolBar = self.addToolBar(self.tr("PmvToolBar"))
        self.toolBar.setObjectName("ToolBar")
        self.toolBar.addAction(self.openAct)
        self.toolBar.addAction(self.saveSessionAct)
        self.toolBar.addAction(self.undoAct)
        self.toolBar.addAction(self.redoAct)
        self.toolBar.addAction(self.clearSelAct)
        self.toolBar.addAction(self.focusSceneAct)        
        self.toolBar.addAction(self.reportAct)
        self.toolBar.addAction(self.togglePyShellAct)
        #self.toolBar.addAction(self.snapshots)
        
        # put the Exit icon to the right side of the toolbar by inserting
        # this "spacer" widget before adding the exitAct.
        spacer = QtGui.QWidget()
        spacer.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        self.toolBar.addWidget(spacer)
        self.toolBar.addAction(self.exitAct)

        self.toolBar.setStyleSheet("QToolBar { background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1, stop: 0 #a6a6a6, stop: 0.08 #7f7f7f, stop: 0.39999 #717171, stop: 0.4 #626262, stop: 0.9 #4c4c4c, stop: 1 #333333);}")

    def selectString_cb(self, text):
        self.pmv.clearSelection()
        if text[-1] in ['#', '/', '+', '-', '&'] or text.endswith('to'):
            return
        sel = self.pmv.Mols.selectMolKit(text)
        if len(sel):
            self.pmv.select(sel)
        #except:
        #    pass
        
    def createDockedtWidgets(self):
        # create a viewer and make it the central widget
        dock = QtGui.QDockWidget('3DViewer', self)
        dock.setMinimumSize(300, 200)
        from DejaVu2.Qt.Viewer import Viewer
        #self.viewer = vi = Viewer(master=dock)
        self.viewer = vi = PmvViewer(self.pmv)#, master=dock)
        vi.setStatusBar(self.statusBar())
        dock.setWidget(vi.cameras[0])
        self.addDockWidget(QtCore.Qt.RightDockWidgetArea, dock)
        #self.setCentralWidget(vi.cameras[0])
        self.setCentralWidget(dock)
        self.toggle3DViewViewAct = dock.toggleViewAction()

        glayout = QtGui.QGridLayout()
        dock.setLayout(glayout)
        glayout.addWidget(vi, 0, 0, 1, 1)
        w = self.stringSelectorWidget = QtGui.QLineEdit(dock)
        w.textChanged.connect(self.selectString_cb)
        #glayout.addWidget(w, 1, 1, 1, 1)
        #w.show()
        # override picking

        # the Viewer.processPicking is called with the pick object
        #self.viewer.processPicking = self.processPicking

        # self.processPicking will turn the pick into atoms and call self.setDragSelectCommand
        #self.setDragSelectCommand(self.pmv.select)

        # this allows to specify what should happen when we pick on nothing, Set to None
        # if no action is wanted
        #self.setEmptyPickCommand(None)
        
        # override picking
        ## vi.cameras[0]._mouseReleaseNoMotionActions[
        ##     int(QtCore.Qt.LeftButton)][
        ##     int(QtCore.Qt.NoModifier)] = self.processPicking

        ## vi.cameras[0]._mouseReleaseWithMotionActions[
        ##     int(QtCore.Qt.LeftButton)][
        ##     int(QtCore.Qt.ShiftModifier)] = self.pick

        # create the dashboard widget
        dock = QtGui.QDockWidget('Dashboard', self)
        dock.setObjectName('dock')
        from dashboard import Dashboard
        self.objTree = Dashboard(self, parent=dock)
        self.objTree.setObjectName('Dashboard')
        dock.setWidget(self.objTree)
        dock.setMinimumSize(50, 50)
        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, dock)
        self.toggleDashboardViewAct = dock.toggleViewAction()
        
        dock = self.pyShellDockWidget = QtGui.QDockWidget(self.tr("Python Shell"), self)
        dock.setObjectName("PythonShell")

        if not use_ipython_shell:        
            # create the PyShell widget
            from pyshell import PyShell
            self.pyShellWidget = PyShell(dock)
        else :
            #create ipython embeded shell
            kernel_manager = QtInProcessKernelManager()
            kernel_manager.start_kernel()
            kernel = kernel_manager.kernel
            kernel.gui = 'qt4'
            #other variable to pass to the context ?
            kernel.shell.push({'pmv':self.pmv,'pmvgui': self})
            #do you want pylab ?
            #kernel.shell.run_cell("import matplotlib")
            #kernel.shell.run_cell("matplotlib.rcParams['backend.qt4']='PySide'")
            #kernel.shell.run_cell("%pylab")
#            kernel.shell
            #matplotlib.rcParams['backend.qt4']='PySide'
            kernel_client = kernel_manager.client()
            kernel_client.start_channels()
            app = guisupport.get_app_qt4()
            def stop():
                kernel_client.stop_channels()
                kernel_manager.shutdown_kernel()
                app.exit()
    
            self.pyShellWidget = RichIPythonWidget()
            self.pyShellWidget.kernel_manager = kernel_manager
            self.pyShellWidget.kernel_client = kernel_client
            self.pyShellWidget.exit_requested.connect(stop)
        
        dock.setWidget(self.pyShellWidget)
        self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, dock)
        self.toggleShellViewAct = dock.toggleViewAction()
        dock.hide()

        # create the Log widget
        dock = QtGui.QDockWidget(self.tr("Log"), self)
        dock.setObjectName("Log")
        #dock.setAllowedAreas(QtCore.Qt.LeftDockWidgetArea | QtCore.Qt.RightDockWidgetArea)
        # prevent widget from being closeable
        #dock.setFeatures(QtGui.QDockWidget.DockWidgetMovable |
        #                 QtGui.QDockWidget.DockWidgetFloatable)
        self.logWidget = QtGui.QListWidget(dock)
        self.logWidget.addItems(["## Welcome to PMV",
                                 "## commands issues in PMV will log themselves here"])
        dock.setWidget(self.logWidget)
        self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, dock)
        self.toggleLogViewAct = dock.toggleViewAction()
        dock.hide()
        
        # create cmd line widget
        ## dock = QtGui.QDockWidget(self.tr("cmdline"), self)
        ## self.cmdLineWidget = QtGui.QLineEdit(dock)
        ## self.cmdLineWidget.setFocus()
        ## dock.setWidget(self.cmdLineWidget)
        ## self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, dock)
        ## self.toggleCmdlineViewAct = dock.toggleViewAction()

        ## dock = QtGui.QDockWidget(self.tr("Python Shell"), self)
        ## self.PyShellWidget =PyShell(dock)
        ## self.PyShellWidget.setFocus()
        ## dock.setWidget(self.PyShellWidget)
        ## self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, dock)
        ## self.togglePyShellViewAct = dock.toggleViewAction()

    def keyPressed_cb(self, key):        
        """This method is connected to Qt Signal defined in self.viewer.camera.
        When 'f' key is pressed the emit('f') of the Signal is issued.
        """
        if key=='f' or key=="F":
            self.pmv.focusScene()
        

    ## def remove_cb(self, obj, root):
    ##     if isinstance(obj, (Molecule, Atom)):
    ##         if root:
    ##             if isinstance(root._pmvObject, Molecule):
    ##                 print 'Delete Molecular fragment', obj.full_name()
    ##             elif isinstance(root._pmvObject, Selection):
    ##                 print 'Deselect %s under %s'%(obj.full_name(), root.text(0))
    ##             elif isinstance(root._pmvObject, Group):
    ##                 item = obj._treeItems[root]
    ##                 print 'Remove %s from group %s'%(obj.full_name(), item.parent().text(0))
    ##         else:
    ##             print 'Delete Molecular fragment', obj.full_name()
                
    ##     elif isinstance(obj, PmvSelection):
    ##         print 'Delete selectino', obj.name
    ##     elif isinstance(obj, Group):
    ##         if obj == root._pmvObject:
    ##             print 'Delete group %s'%(obj.name)
    ##         else:
    ##             print 'Delete group %s under %s'%(obj.name, root.text(0))
        #self.objTree.removeItem(items)

    ## def openFile(self):
    ##     filename, filter = QtGui.QFileDialog.getOpenFileName(
    ##         self, self.tr("Read molecule or session"),
    ##         '',
    ##         self.tr("Molecule Files (*.pdb *.mol2);; Session Files (*.psf)"))
    ##     if filename:
    ##         # FIXME .. how to handle unicode names ?
    ##         filename = filename.encode('ascii', 'ignore')
    ##         #filename = filename.encode('ascii', 'replace')
    ##         mols = self.pmv.readMolecules.guiCall([filename], header=True)#)
        
    ##         ###for mol in mols:
    ##         ###    self.pmv.buildBondsByDistance(mol)
    ##         # FIXME we hardwire the the commands used upon loading
    ##         #if mols:
    ##             #self.pmv.displayLines.onAddObjectToViewer(mols[0])
    ##             #mols[0].geomContainer.geoms['bonded'].Set(
    ##             #    faces=mols[0].patoms._bonds, visible=1)
    ##             #self.pmv.displayLines.guiCall(mols[0].makeSelection())
    ##             #self.pmv.displayLines.guiCall(mol.name)
    ##             ###self.pmv.colorByAtomType.guiCall(mols, geomsToColor=['lines'])
   
    def registerListerners(self):
        ##
        evh = self.pmv.eventHandler

        ## from AppFramework.App import AddGeometryEvent, RemoveGeometryEvent, \
        ##      RedrawEvent
        ## evh.registerListener(AddGeometryEvent, self.addGeometryEventHandler)
        ## evh.registerListener(RemoveGeometryEvent,
        ##                      self.removeGeometryEventHandler)
        ## evh.registerListener(RedrawEvent, self.redrawEventHandler)

        from PmvApp.Pmv import AfterAddMoleculeEvent
        evh.registerListener(AfterAddMoleculeEvent,
                             self.onNewMolDisplayed)

        from AppFramework.notOptionalCommands import NewUndoEvent, \
             AfterUndoEvent
        evh.registerListener(NewUndoEvent, self.undoEventHandler)
        evh.registerListener(AfterUndoEvent, self.undoEventHandler)

        ## reports
        from AppFramework.AppCommands import ExecutionReportEvent
        evh.registerListener(ExecutionReportEvent,
                             self.executionReportEventHandler)
        #from execReport import RemoveReportsEvent
        #evh.registerListener(RemoveReportsEvent,
        #                     self.removeReportsEventHandler)
        
        ##
        ## Selection Events        
        ##from PmvApp.selectionCmds import SelectionEvent

        from PmvApp.Pmv import \
             ActiveSelectionChangedEvent, DeleteNamedSelectionEvent, \
             RenameSelectionEvent, RenameGroupEvent, RenameTreeNodeEvent

        evh.registerListener(ActiveSelectionChangedEvent,
                             self.activeSelectionChangedEventHandler)
        evh.registerListener(DeleteNamedSelectionEvent,
                             self.deleteNamedSelectionEventHandler)
        evh.registerListener(RenameSelectionEvent,
                             self.renameSelectionEventHandler)
        evh.registerListener(RenameTreeNodeEvent,
                             self.renameTreeNodeEventHandler)

        ##
        ## groups
        from PmvApp.Pmv import AddGroupEvent, DeleteGroupsEvent, ReparentGroupObject
        evh.registerListener(AddGroupEvent, self.addGroupEventHandler)
        evh.registerListener(DeleteGroupsEvent, self.deleteGroupsEventHandler)
        evh.registerListener(ReparentGroupObject, self.reparentGroupObjectHandler)
        evh.registerListener(RenameGroupEvent,
                             self.renameGroupEventHandler)

        ##
        ## Delete Events
        ## from PmvApp.deleteCmds import BeforeDeleteMoleculeEvent
        ## from PmvApp.Pmv import AfterDeleteAtomsEvent

        ## evh.registerListener(BeforeDeleteMoleculeEvent,
        ##                      self.beforeDeleteMoleculeEventHandler)

        ## evh.registerListener(AfterDeleteAtomsEvent,
        ##                      self.afterDeleteAtomsEventHandler)

    def updateMultiMoleculeCount(self, mol):
        name = '%s %s %d/%d'%(mol._basename, mol.name,
                              mol.currentMoleculeIndex+1, len(mol.index))
        item = mol._treeItems.keys()[0]
        item.setText(0, name)
        if mol._hasIndex:
            self._timer.stop()

    def onNewMolDisplayed(self, event):
        mol = event.molecule
        if hasattr(event, 'atIndex'):
            atIndex = event.atIndex
        else:
            atIndex = None # add at the end
        if hasattr(event, 'select'):
            wasSelected = event.select
        else:
            wasSelected = False # not selected
        item = self.objTree.addObject( mol, None, mol._basename,
                                       atIndex=atIndex)
        item.setSelected(wasSelected)

        if mol._group is not None:
            self.reparentGroupObject(mol, mol._group)
            
        if mol._multi=='molecules':
            if mol._hasIndex is False:
                self._timer = QtCore.QTimer(self)
                cb = CallbackFunction(self.updateMultiMoleculeCount, mol)
                self._timer.timeout.connect(cb)
                self._timer.start(1000)
            else:
                item = mol._treeItems.keys()[0]
                item.setText(0, "%s %s %d/%d"%(
                    mol._basename, mol.name, mol.currentMoleculeIndex+1, len(mol.index)))
        elif mol._multi=='conformations':
            mol.name = '%s %d/%d'%(mol._basename, mol._ag.getACSIndex()+1,
                                   mol._ag._coords.shape[0])
            item = mol._treeItems.keys()[0]
            item.setText(0, mol.name)
            
        app = self.pmv

        if not hasattr(event, 'ignoreCenterPref') or event.ignoreCenterPref==False:
            center = False
            if app.userpref['Center Scene']['value']=='always':
                center = True
            elif app.userpref['Center Scene']['value']=='firstMoleculeOnly' and \
                 len(app.Mols)==1:
                center = True
            if center:
                vi = self.viewer
                vi.Reset_cb()
                vi.Normalize_cb()
                vi.Center_cb()
        
        
    ## def addMoleculeEventHandler(self, event):
    ##     import pdb
    ##     pdb.set_trace()
    ##     mol = event.kw['object']
    ##     #name = event.kw['name']
    ##     self.objTree.addObject( mol, None, mol.name)
        
    def undoEventHandler(self, event):
        # set the button tool tip and the status bar message to the
        # the name of the next action undo(redo) could trigger
        cmd = event.command
        if len(event.objects):
            if cmd.name == "undo":
            #if event.objects == self.pmv.undo.cmdStack:
                self.undoAct.setStatusTip('Undo '+event.objects[-1][1])
                self.undoAct.setToolTip('Undo '+event.objects[-1][1])
            #elif event.objects == self.pmv.redo.cmdStack:
            elif cmd.name == "redo":
                self.redoAct.setStatusTip('Redo '+event.objects[-1][1])
                self.redoAct.setToolTip('Redo '+event.objects[-1][1])
            else:
                raise # FIXME
        else:
            if cmd.name == "undo":
                self.undoAct.setStatusTip('Undo (empty stack)')
                self.undoAct.setToolTip('Undo (empty stack)')
            elif cmd.name == "redo":
                self.redoAct.setStatusTip('Redo (empty stack)')
                self.redoAct.setToolTip('Redo (empty stack)')
            #print 'Undo event on empty stack'

    def executionReportEventHandler(self, event):
        report = event.report
        self.updateReportIcon(report)
        
    def updateReportIcon(self, report):
        if report.numberOf['exceptions'] or report._requestUserConfirmation:
            self.unseenReports.append(report)
            self.reportAct.setIcon(QtGui.QIcon(os.path.join(
                PMVICONPATH, 'reportException.png')))
            self._worstError = self._EXCEPTION

        elif report.numberOf['errors']:
            if self._worstError < self._ERROR:
                self._worstError = self._ERROR
                self.reportAct.setIcon(QtGui.QIcon(os.path.join(
                    PMVICONPATH, 'reportError.png')))
            self.unseenReports.append(report)
            
        elif report.numberOf['warnings'] > 1:
            if self._worstError < self._WARNING:
                self._worstError = self._WARNING
                self.reportAct.setIcon(QtGui.QIcon(os.path.join(
                    PMVICONPATH, 'reportWarning.png')))
            self.unseenReports.append(report)
        elif report.numberOf['successes']:
            self.unseenReports.append(report)
            

    def displayReports(self):
        # reset icons and state
        self._worstError = self._NoERROR
        self.reportAct.setIcon(QtGui.QIcon(os.path.join(PMVICONPATH,
                                                        'report.png')))
        
        if len(self.unseenReports)==0: return
        from execReport import ExecutionReports
        w = ExecutionReports(self.pmv, self.unseenReports)
        ## w.setModal(False)
        ## w.show()
        ## w.raise_()
        ## w.activateWindow()
        ok = w.exec_()
        #report.printReport()

    def showPythonShell(self):
        """toggles show/hide of the Python Shell. PyShell is
        displayed in an undocked window."""
        isFloating = self.pyShellDockWidget.isFloating()
        isVisible =  self.pyShellDockWidget.isVisible()
        #show = not (isFloating and isVisible)
        show = not isVisible
        #print "pyshell floating", isFloating, "visible:", isVisible, "show:", show
        #self.pyShellDockWidget.setFloating(show)
        self.pyShellDockWidget.setVisible(show)

    def showSnapshotsGUI(self):
        from DejaVu2.Qt.snapshots import SnapshotDialog
        ok = SnapshotDialog.addSnapshots(self)

    ## def removeReportsEventHandler(self, event):
    ##     self._worstError = self._NoERROR
    ##     self.reportAct.setIcon(QtGui.QIcon(os.path.join(PMVICONPATH,
    ##                                                     'report.png')))
    ##     for report in self.unseenReports:
    ##         self.updateReportIcon(report)
            
    def addGroupEventHandler(self, event):
        parent = event.parent
        if parent:
            parent = parent._treeItems.keys()[0]
        self.objTree.addObject(event.group, parent, event.group.name)

    def deleteGroupsEventHandler(self, event):
        for group in event.groups:
            groupItem = group._treeItems.keys()[0]
            if groupItem.parent() is None:
                self.objTree.invisibleRootItem().removeChild(groupItem)
            else:
                groupItem.parent().removeChild(groupItem)
            
    def reparentGroupObjectHandler(self, event):
        self.reparentGroupObject(event.object, event.newGroup)

    def reparentGroupObject(self, obj, group):
        # remove object's item from old group
        if hasattr(obj, '_treeItems'):
            objItem = obj._treeItems.keys()[0]
            parent = objItem.parent()
            if parent is None:
                self.objTree.invisibleRootItem().removeChild(objItem)
            else:
                parent.removeChild(objItem)

        # add object's item to new group
        if group is None:
            self.objTree.invisibleRootItem().addChild(objItem)
            self.objTree.setItemsExpandable(True)
        elif hasattr(group, '_treeItems'):
            groupItem = group._treeItems.keys()[0]
            collapse = False
            if not groupItem.isExpanded():
                collapse = True
            groupItem.addChild(objItem)
            if collapse:
                groupItem.setExpanded(False)
                
    def renameGroupEventHandler(self, event):
        obj = event.object
        item = event.item
        if item is None:
            item = obj._treeItems.keys()[0]
        item.setText(0, obj.name)
            
    def afterDeleteAtomsEventHandler(self, event):
        roots = {}
        for atom in event.objects:
            if hasattr(atom, '_treeItems'):
                for rootItem, item in [atom._treeItems.keys()[0]]:
                    parentItem = item.parent()
                    parentItem.removeChild(item)
                    while parentItem is not None:
                        if parentItem.childCount()==0:
                            newparentItem = parentItem.parent()
                            newparentItem.removeChild(parentItem)
                            parentItem = newparentItem
            else:
                # this happens if we delete all the atoms in a bunch of residues
                # for instance. The atoms are not in the dahsboard, but the residues might be
                # and need to be removed
                
                # go up the tree until we find an parent in the Dashboard
                parent = atom.parent
                while not hasattr(parent, '_treeItems'):
                    parent = parent.parent
                roots[parent] = True
                
        # now traverse all molecules and delete all entries with no children
        nodes = roots.keys()
        parentItems = {}
        for node in roots.keys():
            if len(node.children)==0:
                for it in node._treeItems.values():
                    parent = it.parent()
                    if parent:
                        parent.removeChild(it)
                        parentItems[parent] = True

        # node clean up all parent in the tree
        items = parentItems.keys()
        while len(items):
            parentItems = {}
            for item in items:
                if item.childCount()==0:
                    parent = item.parent()
                    if parent:
                        parent.removeChild(item)
                        parentItems[parent] = True
                    elif isinstance(item._pmvObject, (Protein, Group)):
                        self.objTree.invisibleRoot().removeChild(item)
                    elif isinstance(item._pmvObject, (Selection)):
                        if item._pmvObject==self.pmv.curSelection:
                            self.objTree.invisibleRoot().removeChild(item)
                            self.activeSelection = None
            items = parentItems.keys()

            
    def beforeDeleteMoleculeEventHandler(self, event):
        mol = event.object
        # remove the molecule form the tree
        item = mol._treeItems.keys()[0]
        if item.parent() is None:
            self.objTree.invisibleRootItem().removeChild(item)
        else:
            item.parent().removeChild(item)
        # remove item's ._pmvObject attribute to break cyclic reference
        def delPmvObject(item):
            del item._pmvObject
            for n in range(item.childCount()):
                child = item.child(n)
                if hasattr(child, '_pmvObject'): # dummy child of unexpanded node does not have _pmvObject
                    delPmvObject(child)
        delPmvObject(item)

    def renameTreeNodeEventHandler(self, event):
        obj = event.object
        if obj.alias:
            newName = "%s (%s)"%(obj.alias, obj.name)
        else:
            newName = obj.name
        if hasattr(obj, '_treeItems'):
            obj._treeItems.keys()[0].setText(0, newName)

    def renameSelectionEventHandler(self, event):
        item = event.item
        obj =  event.selection
        if item is None:
            if hasattr(obj, '_treeItems'):
                item = obj._treeItems.keys()[0]
            else:
                return
        item.setText(0, obj.name)
        if hasattr(event, 'setCurrent'): # we renamed the curSelection
            if event.setCurrent:
                self.pmv.setActiveSelection(obj)

    def deleteNamedSelectionEventHandler(self, event):
        sele = event.selection
        if self.activeSelection == sele:
            self.activeSelection = None
        if hasattr(sele, '_treeItems'):
            self.objTree.invisibleRootItem().removeChild(sele._treeItems.keys()[0])
            
    def activeSelectionChangedEventHandler(self, event):
        old = event.old
        new = event.new
        self.activeSelection = new
        if hasattr(new, '_treeItems'):
            self.objTree.setCurrentSelection(new._treeItems.keys()[0])
