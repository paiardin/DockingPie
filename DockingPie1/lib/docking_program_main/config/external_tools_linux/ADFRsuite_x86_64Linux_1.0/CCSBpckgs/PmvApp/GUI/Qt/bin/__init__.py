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
# Copyright: M. Sanner TSRI 2014
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/PmvApp/GUI/Qt/bin/__init__.py,v 1.2.4.1 2017/07/13 20:48:47 annao Exp $
#
# $Id: __init__.py,v 1.2.4.1 2017/07/13 20:48:47 annao Exp $
#
import sys

def runPmv2():
    from PySide import QtGui, QtCore
    app = QtGui.QApplication(sys.argv)

    from PmvApp.GUI.Qt.PmvGUI import PmvGUI
    from PmvApp.Pmv import MolApp

    pmv = MolApp()
    pmv.lazyLoad('bondsCmds', commands=[
        'buildBondsByDistance',],
                 package='PmvApp')
    pmv.lazyLoad('fileCmds', commands=[
        'readMolecules', 'readPmvSession', 'fetch', 'readAny', 'writePDB'],
                 package='PmvApp')
    pmv.lazyLoad('displayCmds', commands=[
        'displayLines', 'undisplayLines', 'displayCPK', 'undisplayCPK',
        'displaySticksAndBalls', 'undisplaySticksAndBalls',
        'displayBackboneTrace', 'undisplayBackboneTrace',
        'displayBoundGeom', 'undisplayBoundGeom','bindGeomToMolecularFragment',
        'showMolecules'],
                 package='PmvApp')

    cmds = ['color', 'colorByAtomType', 'colorByResidueType',
            'colorAtomsUsingDG', 'colorResiduesUsingShapely',
            'colorByChains', 'colorByMolecules', 'colorByInstance',
            'colorByProperty', 'colorRainbow', 'colorRainbowByChain',
            'colorByExpression', 'colorByLinesColor']
    pmv.lazyLoad("colorCmds", commands=cmds, package="PmvApp")
    pmv.lazyLoad("selectionCmds", package="PmvApp")


    pmv.lazyLoad('msmsCmds', package='PmvApp')

    pmvgui = PmvGUI(pmv)
    pmvgui.resize(800,600)

    ## cam = pmvgui.viewer.cameras[0]
    ## cam.Set(color=(.7,.7,.7))

    ## cursor_px = QtGui.QPixmap('/mgl/ms1/people/sanner/python/lazy/PmvApp/GUI/Icons/selectionAdd.png')
    ## #cursor_px = QtGui.QPixmap('/mgl/ms1/people/sanner/python/lazy/PmvApp/GUI/Icons/colorChooser32.png')
    ## cursor_px.setMask(cursor_px.mask())
    ## cursor = QtGui.QCursor(cursor_px)
    ## cam.setCursor(cursor)



    for arg in sys.argv[1:]:
        mols = pmv.readMolecule(arg)
        for mol in mols:
            pmv.buildBondsByDistance(mol)
        pmv.displayLines(mols)

    ## print 'DDDD', len(pmv.Mols[0].allAtoms[:-18])

    ## from time import time
    ## t0 = time()
    ## pmv.selection.set(pmv.Mols[0].allAtoms[:-18])
    ## print 'FASFD', time()-t0

    #pmvgui.createNewGroup('Group 1')
    #pmvgui.createNewGroup('Group 2')
    #pmvgui.createNewGroup('Group 3')
    #pmvgui.createNewGroup('Group 4')
    #pmvgui.createNewGroup('Group 5')
    sys.exit(app.exec_())
