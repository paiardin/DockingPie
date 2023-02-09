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
# Copyright: M. Sanner TSRI 2013
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/PmvApp/GUI/Qt/execReport.py,v 1.7.4.1 2017/07/13 20:46:53 annao Exp $
#
# $Id: execReport.py,v 1.7.4.1 2017/07/13 20:46:53 annao Exp $
#
import sys, os
from mglutil.util.packageFilePath import findFilePath
PMVICONPATH = findFilePath('Icons', 'PmvApp.GUI')

def contentFont():
        font = QtGui.QFont()
        font.setStyleStrategy(QtGui.QFont.PreferAntialias)
 
        #if sys.platform == 'darwin':
        #    font.setPixelSize(9)
        #    font.setFamily('Arial')
        #else:
        #    font.setPixelSize(11)
        #    font.setFamily('Arial')
        font.setPixelSize(12)
        font.setFamily('Arial')
 
        return font

from PySide import QtGui, QtCore

class SummaryTab(QtGui.QWidget):
    def __init__(self, report, parent=None):
        super(SummaryTab, self).__init__(parent)

        title = QtGui.QLabel("Command execution report for %s:"%report.cmd.name)
        
        succLab = QtGui.QLabel("    %d success(es)"%report.numberOf['successes'])
        warnLab = QtGui.QLabel("    %d warning(s)"%report.numberOf['warnings'])
        errorLab = QtGui.QLabel("    %d error(s)"%report.numberOf['errors'])
        exceptLab = QtGui.QLabel("    %d exception(s)"%report.numberOf['exceptions'])
        timeLab = QtGui.QLabel("    %.2f total time (s)"%(report.endTime-report.startTime))

        mainLayout = QtGui.QVBoxLayout()
        mainLayout.addWidget(title)
        mainLayout.addWidget(succLab)
        mainLayout.addWidget(warnLab)
        mainLayout.addWidget(errorLab)
        mainLayout.addWidget(exceptLab)
        mainLayout.addWidget(timeLab)
        #mainLayout.addStretch(1)
        self.setLayout(mainLayout)

class SummaryTab(QtGui.QWidget):
    def __init__(self, reports, parent=None):
        super(SummaryTab, self).__init__(parent)

        title = QtGui.QLabel("%d command execution reports"%len(reports))


        mainLayout = QtGui.QVBoxLayout()
        mainLayout.addWidget(table)
        self.setLayout(mainLayout)

class WizardTab(QtGui.QWizard):
    def __init__(self, objects, app, name, parent=None):
        super(WizardTab, self).__init__(parent)
        for num, obj in enumerate(objects):
            page = QtGui.QWizardPage()
            page.setTitle("%s %d"%(name, num+1))
            text = QtGui.QTextEdit()
            for line in obj.getAll(app, "", details=False):
                text.append(line)
            text.setReadOnly(True)
            
            layout = QtGui.QVBoxLayout()
            layout.addWidget(text)
            page.setLayout(layout)
            self.addPage(page)
        self.setWindowTitle("Trivial Wizard")
        self.show()
           
class EventsTab(QtGui.QWidget):
    def __init__(self, objects, details=True, parent=None):
        super(EventsTab, self).__init__(parent)
        table = QtGui.QTableWidget()
	header = table.horizontalHeader()
        if details:
            table.setColumnCount(2)
            #table.resize(800, 250)
            table.setSizePolicy(QtGui.QSizePolicy.Expanding,  QtGui.QSizePolicy.Expanding);
        else:
            table.setColumnCount(1)
            header.setResizeMode(0, QtGui.QHeaderView.Stretch)
        table.setRowCount(len(objects))
        table.setHorizontalHeaderLabels(['message', 'exception'])
        if details:
            header.setResizeMode(1, QtGui.QHeaderView.Stretch)
        else:
            header.setResizeMode(0, QtGui.QHeaderView.Stretch)
            
        row = 0
        for obj in objects:
            item =  QtGui.QTableWidgetItem(obj.msg)
            table.setItem(row, 0, item)
            if obj.exception:
                msg = obj.exception.message
                if msg:
                    item =  QtGui.QTableWidgetItem(msg)
                    item._obj = obj
                table.setItem(row, 1, item)
            row += 1

        #table.resizeColumnsToContents()

        if details:
            self.details = QtGui.QTextEdit()
            self.details.setReadOnly(True)
            table.currentItemChanged.connect(self.onSetCurrentItem)
        
        layout = QtGui.QVBoxLayout()
        layout.addWidget(table)
        if details:
            layout.addWidget(self.details)
        self.setLayout(layout)

    def onSetCurrentItem(self, item):
        self.details.clear()
        if hasattr(item, '_obj'):
            for line in item._obj.getFormatedException():
                self.details.append(line)
            for line in item._obj.getContext():
                self.details.append(line)
        
class TimingTab(QtGui.QWidget):
    def __init__(self, report, parent=None):
        super(TimingTab, self).__init__(parent)
        t0 = report.startTime
        topLabel = QtGui.QLabel("Timing:")

        timesTable = QtGui.QTableWidget()
        timesTable.setColumnCount(2)
        timesTable.setRowCount(len(report.timeStamps))
        timesTable.setHorizontalHeaderLabels(['Time', 'time (s)'])

        header = timesTable.horizontalHeader()
        header.setResizeMode(1, QtGui.QHeaderView.Stretch)
        
        row = 0
        for t, msg in report.timeStamps:
            item =  QtGui.QTableWidgetItem('%.2f'%(t-t0))
            timesTable.setItem(row, 0, item)
            t0 = t
            item =  QtGui.QTableWidgetItem(msg)
            timesTable.setItem(row, 1, item)
            row += 1

        timesTable.resizeColumnsToContents()
        
        layout = QtGui.QVBoxLayout()
        layout.addWidget(topLabel)
        layout.addWidget(timesTable)
        self.setLayout(layout)
    

class ExecutionReport(QtGui.QWidget):

    def __init__(self, report, title, parent=None):
        super(ExecutionReport, self).__init__(parent)

        title = QtGui.QLabel(title, parent)
        
        tabWidget = QtGui.QTabWidget(parent)
        successes = report.getSuccesses()
        if len(successes):
            tabWidget.addTab(EventsTab(successes, details=False), "Success (%d)"%len(successes))

        warnings = report.getWarnings()
        if len(warnings):
            tabWidget.addTab(EventsTab(warnings), "Warnings (%d)"%len(warnings))

        errors = report.getErrors()
        if len(errors):
            tabWidget.addTab(EventsTab(errors), "Errors (%d)"%len(errors))

        exceptions = report.getExceptions()
        if len(exceptions):
            tabWidget.addTab(EventsTab(exceptions), "Exceptions (%d)"%len(exceptions))

        tabWidget.addTab(TimingTab(report), 'Timing')
        
        mainLayout = QtGui.QVBoxLayout()
        mainLayout.addWidget(title)
        mainLayout.addWidget(tabWidget)
        self.setLayout(mainLayout)

        #self.adjustSize()

class ExecutionReports(QtGui.QDialog):

    def __init__(self, app, reports, parent=None):
        super(ExecutionReports, self).__init__(parent)
        self.setFont(contentFont())
        self.reports = reports
        self.app = app
        # Widget that contains the table and the stack of reports
        stackAndTable = QtGui.QWidget()

        # widget taht contains table and button to delete reports in table
        tableAndControls = QtGui.QWidget(stackAndTable)

        title = QtGui.QLabel('Reports', parent)

        # create a table of reports
        table = QtGui.QTableWidget(tableAndControls)
        table.setRowCount(len(reports))
        table.setColumnCount(5)
        table.setHorizontalHeaderLabels(['command', 'Ex', 'Er', 'W', 'S'])

        header = table.horizontalHeader()
        header.setResizeMode(0, QtGui.QHeaderView.Stretch)
        row = 0
        for report in reports:
            item0 =  QtGui.QTableWidgetItem(report.cmd.name)
            table.setItem(row, 0, item0)
            item1 =  QtGui.QTableWidgetItem(str(report.numberOf['exceptions']))
            table.setItem(row, 1, item1)
            item2 =  QtGui.QTableWidgetItem(str(report.numberOf['errors']))
            table.setItem(row, 2, item2)
            item3 =  QtGui.QTableWidgetItem(str(report.numberOf['warnings']))
            table.setItem(row, 3, item3)
            item4 =  QtGui.QTableWidgetItem(str(report.numberOf['successes']))
            table.setItem(row, 4, item4)
            report.items = [item0, item1, item2, item3, item4]
            row += 1

        table.setRangeSelected(QtGui.QTableWidgetSelectionRange(0, 0, 0, 4),
                               True)
        table.cellClicked.connect(self.onCellClick)
        self.table = table

        # create controls for the table
        buttonBox = QtGui.QWidget(tableAndControls)
        ## deleteButton = QtGui.QPushButton("Current", buttonBox)
        ## deleteButton.setIcon(QtGui.QIcon(os.path.join(PMVICONPATH, 'removex.png')))
        ## deleteButton.setToolTip('Delete curently selected report')
        ## deleteButton.released.connect(self.deleteCurrent)
            
        #deleteButtonAll = QtGui.QPushButton("Clear All", buttonBox)
        deleteButtonAll = QtGui.QPushButton(buttonBox)
        deleteButtonAll.setIcon(QtGui.QIcon(os.path.join(PMVICONPATH, 'trash_empty.png')))
        deleteButtonAll.setToolTip('Delete all reports from the table')
        deleteButtonAll.released.connect(self.deleteAll)
        deleteButtonAll.setIconSize(QtCore.QSize(24,24))
        
        #saveButton = QtGui.QPushButton("Save...", buttonBox)
        saveButton = QtGui.QPushButton(buttonBox)
        saveButton.setIcon(QtGui.QIcon(os.path.join(PMVICONPATH, 'open.png')))
        saveButton.setToolTip('Save all Reports')
        saveButton.released.connect(self.saveReports)
        saveButton.setIconSize(QtCore.QSize(24,24))

        #submitButton = QtGui.QPushButton("submit ...", buttonBox)
        submitButton = QtGui.QPushButton(buttonBox)
        submitButton.setIcon(QtGui.QIcon(os.path.join(PMVICONPATH, 'mail.png')))
        submitButton.setToolTip('Send these reports to the developers')
        submitButton.released.connect(self.submitReports)
        submitButton.setIconSize(QtCore.QSize(24,24))

        layout = QtGui.QHBoxLayout()
        #layout.addWidget(deleteButton)
        layout.addWidget(deleteButtonAll)
        layout.addWidget(saveButton)
        layout.addWidget(submitButton)
        buttonBox.setLayout(layout)

        layout = QtGui.QVBoxLayout()
        layout.addWidget(title)
        layout.addWidget(table)
        layout.addWidget(buttonBox)
        tableAndControls.setLayout(layout)
        
        ## I can not find how to make the table fit the data nicely :(
        ## best hack I could come up with so far
        
        #table.resize(table.maximumViewportSize())
        # compute width
        width = 0
        vheader = table.verticalHeader()
        width = 50 #vheader.sizeHintForColumn(0)
        for i in range(5):
            header.setResizeMode(i, QtGui.QHeaderView.ResizeToContents)
            width += table.sizeHintForColumn(i)
        #print 'WIDTH', width
        table.setMaximumWidth(width)
            
        sizePol = QtGui.QSizePolicy()
        sizePol.setHorizontalPolicy(QtGui.QSizePolicy.Minimum)
        sizePol.setVerticalPolicy(QtGui.QSizePolicy.Expanding)
        table.setSizePolicy(sizePol)

        ## create the stacked widget
        self.stackedWidget = stack = QtGui.QStackedWidget(stackAndTable)
        for i, report in enumerate(reports):
            title = 'Details for execution report %d for command %s'%(
                i, report.cmd.name)
            w = ExecutionReport(report, title, parent=stack)
            stack.addWidget( w )
        buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok)

        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)

        layout = QtGui.QHBoxLayout()
        layout.addWidget(tableAndControls)
        layout.addWidget(stack)
        stackAndTable.setLayout(layout)

        mainLayout = QtGui.QVBoxLayout()
        mainLayout.addWidget(stackAndTable)
        mainLayout.addWidget(buttonBox)
        self.setLayout(mainLayout)

        if sys.platform=='linux2':
            self.resize(600, 400)

    def saveReports(self):
        print 'NOT IMPLEMENTED YET'

    def submitReports(self):
        print 'NOT IMPLEMENTED YET'
        
    ## def deleteCurrent(self):
    ##     #reportIndex = self.table.selectedRanges()[0].topRow()
    ##     reportIndex = self.table.currentRow()
    ##     report = self.reports.pop(reportIndex)
    ##     self.table.removeRow(reportIndex)
    ##     if reportIndex > 0:
    ##         newRow = reportIndex-1
    ##     elif self.table.rowCount() >=  reportIndex:
    ##         newRow = reportIndex
    ##     else:
    ##         newRow = None

    ##     if newRow is not None:
    ##         self.stackedWidget.setCurrentIndex(newRow)
    ##         w = self.stackedWidget.widget(reportIndex)
    ##         self.stackedWidget.removeWidget(w)
    ##         self.table.setRangeSelected(QtGui.QTableWidgetSelectionRange(
    ##             newRow, 0, newRow, 4), True)
    ##     event = RemoveReportsEvent(reports=[report])
    ##     self.app.eventHandler.dispatchEvent(event)
        
    def deleteAll(self):
        for i in range(len(self.reports)):
            w = self.stackedWidget.widget(0)
            self.stackedWidget.removeWidget(w)
            self.table.removeRow(0)
        reports = self.reports[:]
        self.app.gui().unseenReports = []
        #event = RemoveReportsEvent(reports=reports)
        #self.app.eventHandler.dispatchEvent(event)
        
    def onCellClick(self, row, column):
        for range in self.table.selectedRanges():
            self.table.setRangeSelected(range, False)
        self.table.setRangeSelected(QtGui.QTableWidgetSelectionRange(
            row, 0, row, 4), True)
        self.stackedWidget.setCurrentIndex(row)
       
