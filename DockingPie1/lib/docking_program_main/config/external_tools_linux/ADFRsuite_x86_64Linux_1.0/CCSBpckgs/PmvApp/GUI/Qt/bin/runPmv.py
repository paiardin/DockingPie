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

# $Header: /mnt/raid/services/cvs/PmvApp/GUI/Qt/bin/runPmv.py,v 1.10.4.1 2017/07/13 20:48:47 annao Exp $
# $Id: runPmv.py,v 1.10.4.1 2017/07/13 20:48:47 annao Exp $

import os, sys, atexit
#sys.path.insert(0,'.')
args = sys.argv
# remove -i from the argv list:
ind = None
debug = False
for i, arg in enumerate(args):
    if arg == "-i":
        ind = i
    elif arg == "-d":
        debug = True

if ind is not None:
    args.pop(ind)
    
from PySide import QtGui, QtCore
app = QtGui.QApplication(args)

from PmvApp import mkPmvApp, loadDefaultCommands
pmv = mkPmvApp()
loadDefaultCommands(pmv)

pmv.trapExceptions = False

if debug:
    pmv._stopOnError = True

#pmv.lazyLoad('bondsCmds', package='PmvApp')
#pmv.lazyLoad('fileCmds', package='PmvApp')
#pmv.lazyLoad('displayCmds', package='PmvApp')
#pmv.lazyLoad('editCmds', package='PmvApp')
#pmv.typeAtomsAndBonds.loadCommand()
#pmv.displayLines.loadCommand()
#pmv.displaySB.loadCommand()
#pmv.lazyLoad("colorCmds", package="PmvApp")
#pmv.lazyLoad("selectionCmds", package="PmvApp")
#pmv.lazyLoad('deleteCmds', package='PmvApp')
#pmv.lazyLoad('secondaryStructureCmds', package='PmvApp')
#pmv.lazyLoad('labelCmds', package='PmvApp')

## from PmvApp.msmsCmds import ComputeMSMS, DisplayMSMS, UndisplayMSMS
## pmv.addCommand(ComputeMSMS(), 'computeMSMS')
## pmv.computeMSMS.loadCommand() # load the command
## #pmv.userpref.set('Compute cavities by default', 'yes')
## pmv.addCommand(DisplayMSMS(), 'displayMSMS')
## pmv.addCommand(UndisplayMSMS(), 'undisplayMSMS')
    
## pmv.lazyLoad('displayHyperBallsCmds', package='PmvApp')
## pmv.lazyLoad('cartoonCmds', package='PmvApp')
## pmv.lazyLoad('interactionsCmds', package='PmvApp')
## pmv.lazyLoad('coarseMolecularSurfaceCmds', package='PmvApp')
## pmv.setOnAddObjectCmd('Molecule', [pmv.displayLines,
##                                    pmv.colorByAtomType,
##                                    pmv.colorByMolecules],
##                       kwList=[{}, {}, {'carbonsOnly':True}])

from PmvApp.GUI.Qt.PmvGUI import PmvGUI,use_ipython_shell
pmvgui = PmvGUI(pmv)

molNames = []
dros = []
scripts = []
session = []
for arg in sys.argv[1:]:
    if os.path.splitext(arg)[1] in ['.pdb','.mol2','.sdf','.pdbqt','.pqr','.pdbq', '.mmtf']:
	molNames.append(arg)
    elif os.path.splitext(arg)[1] in ['.dro']:
	dros.append(arg)
    elif os.path.splitext(arg)[1] in ['.py']:
        scripts.append(arg)
    elif os.path.splitext(arg)[1] in ['.psf']:
        session.append(arg)

if len(molNames):
    mols = pmv.readMolecules(molNames)

if len(dros):
    from ADFR.PmvInterface.fileCmds import DockingResultReader
    pmv.addCommand(DockingResultReader(), 'loadDRO')
    dros = pmv.loadDRO(dros)
    
def runScripts():
    for arg in session:
        pmv.readFullSession(arg)

    for arg in scripts:
        glob = {'pmv':pmv}
        execfile(arg, glob)
        del glob['pmv']
        sys.modules['__main__'].__dict__.update(glob)

# run scripts and session after the gui is displayed
#import pdb; pdb.set_trace()
pmvgui.windowShown.connect(runScripts)
pmvgui.show()
# make sure we delete opened session files if Ctrl-D is used in the Python shell
atexit.register(pmv.cleanup)

#pmv.displayHyperBalls(pmv.Mols[0],shrink=0.01, scaleFactor = 0.26, bScale = 0.26,cpkRad=0.0)
#shrink=0.01, scaleFactor = 0.0, bScale = 1.0,cpkRad=0.3) LIC ?
if use_ipython_shell :
    from IPython.lib import guisupport
    guisupport.start_event_loop_qt4(app)
    #pylab?
else :
    sys.exit(app.exec_())
