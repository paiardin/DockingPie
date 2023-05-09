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
# $Header: /mnt/raid/services/cvs/PmvApp/__init__.py,v 1.2.4.1 2017/07/13 20:55:28 annao Exp $
#
# $Id: __init__.py,v 1.2.4.1 2017/07/13 20:55:28 annao Exp $
#
def mkPmvApp(eventHandler=None):
    # create PmvApp
    from PmvApp.Pmv import MolApp
    pmv = MolApp()
    pmv.trapExceptions = False
    return pmv

def loadDefaultCommands(pmv):
    from PmvApp.msmsCmds import ComputeMSMS, DisplayMSMS, UndisplayMSMS
    pmv.addCommand(ComputeMSMS(), 'computeMSMS')
    pmv.computeMSMS.loadCommand() # load the command
    #pmv.userpref.set('Compute cavities by default', 'yes')
    pmv.addCommand(DisplayMSMS(), 'displayMSMS')
    pmv.addCommand(UndisplayMSMS(), 'undisplayMSMS')

    #pmv.lazyLoad('displayHyperBallsCmds', package='PmvApp')
    pmv.lazyLoad('cartoonCmds', package='PmvApp')
    pmv.lazyLoad('interactionsCmds', package='PmvApp')
    #pmv.lazyLoad('coarseMolecularSurfaceCmds', package='PmvApp')
    pmv.setOnAddObjectCmd('Molecule', [pmv.displayLines,
                                       pmv.colorByAtomType,
                                       pmv.colorByMolecules],
                          kwList=[{}, {}, {'carbonsOnly':True}])
    
    ## pmv.lazyLoad('bondsCmds', package='PmvApp')
    ## pmv.lazyLoad('fileCmds', package='PmvApp')
    ## pmv.lazyLoad('displayCmds', package='PmvApp')
    ## pmv.lazyLoad('editCmds', package='PmvApp')
    ## pmv.displayLines.loadCommand()
    ## pmv.lazyLoad("colorCmds", package="PmvApp")
    ## pmv.color.loadCommand()
    ## pmv.lazyLoad("selectionCmds", package="PmvApp")
    ## pmv.lazyLoad('deleteCmds', package='PmvApp')
    ## pmv.lazyLoad('labelCmds', package='PmvApp')
    
    ## pmv.lazyLoad('msmsCmds', package='PmvApp')
    ## pmv.lazyLoad('displayHyperBallsCmds', package='PmvApp')
    ## pmv.lazyLoad('interactionsCmds', package='PmvApp')
    ## pmv.lazyLoad('coarseMolecularSurfaceCmds', package='PmvApp')

    ## pmv.setOnAddObjectCmd('Molecule', [pmv.displayLines, pmv.colorByAtomType])
