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
# Author: Michel Sanner
#
# $Header: /mnt/raid/services/cvs/ADFR/bin/ADFRgui.py,v 1.1 2016/12/08 18:40:23 sanner Exp $
#
# $Id: ADFRgui.py,v 1.1 2016/12/08 18:40:23 sanner Exp $
#
import sys
from optparse import OptionParser
from PySide import QtCore, QtGui

from ADFR.utils.optParser import ArgParser
parser = ArgParser('ADFR', check=False)
kw = parser.parse_args()

app = QtGui.QApplication(sys.argv)

from PmvApp import mkPmvApp
pmv = mkPmvApp()
    
from PmvApp.GUI.Qt.PmvGUI import PmvViewer
viewer = PmvViewer(pmv)
viewer.cameras[0].setMinimumWidth(250)

from ADFR.GUI.ADFRgui import SingleDockingWidget
widget = SingleDockingWidget(pmvViewer=viewer)    
widget.inputWidget.pmvViewer=viewer

ligFilename = kw['ligandFile']
if ligFilename:
    widget.inputWidget.ligandEntryWidget.setText(ligFilename)

targetFilename = kw['receptorMapsFile']
if targetFilename:
    widget.inputWidget.mapsEntryWidget.setText(targetFilename)

widget.resize(150,300)
widget.show()
sys.exit(app.exec_())
