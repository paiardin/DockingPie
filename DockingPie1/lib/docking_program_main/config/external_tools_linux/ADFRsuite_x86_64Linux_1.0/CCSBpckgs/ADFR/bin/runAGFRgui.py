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
# $Header: /mnt/raid/services/cvs/ADFR/bin/runAGFRgui.py,v 1.3 2017/05/16 17:50:01 sanner Exp $
#
# $Id: runAGFRgui.py,v 1.3 2017/05/16 17:50:01 sanner Exp $
#

import os, sys
from optparse import OptionParser
from PySide import QtCore, QtGui

from ADFR.GUI.gridGUI import GridGUI


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

# use 4.2 by default
#from ADFRcc import setForceFieldVersion
#parameters = setForceFieldVersion('4.2')

(options, args) = parser.parse_args()
app = QtGui.QApplication(sys.argv)
qtpath = os.path.join(os.getenv('MGL_ROOT'), 'plugins')
app.addLibraryPath(qtpath)
widget = GridGUI()

if len(sys.argv)==2:
    widget.pmv.readMolecules(sys.argv[1])
#widget.resize(1000,400)
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
