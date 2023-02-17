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
# Author: Michel F. SANNER, Anna Omelchenko
#
# Copyright: M. Sanner TSRI 2013
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/PmvApp/displayCmdsProxy.py,v 1.1.4.1 2017/07/13 20:55:28 annao Exp $
#
# $Id: displayCmdsProxy.py,v 1.1.4.1 2017/07/13 20:55:28 annao Exp $
#

def getGUI(GUITK):
    if GUITK=='Tk':
        from guiTK.displayCmds import DisplayLinesGUI, DisplayCPKGUI,\
             DisplaySticksAndBallsGUI, DisplayBackboneTraceGUI, DisplayBoundGeomGUI,\
             BindGeomToMolecularFragmentGUI, ShowMoleculesGUI
        return {
            'displayLines':[(DisplayLinesGUI, (), {})],
            'undisplayLines':[(None, (), {})],
            'displayCPK': [(DisplayCPKGUI, (), {})],
            'undisplayCPK':[(None, (), {})],
            'displaySticksAndBalls':[(DisplaySticksAndBallsGUI, (), {})],
            'undisplaySticksAndBalls':[(None, (), {})],
            'displayBoundGeom':[(DisplayBoundGeomGUI, (), {})],
            'undisplayBoundGeom':[(None, (), {})],
            'showMolecules':[(ShowMoleculesGUI, (), {})]
            }
    elif GUITK=='Qt':
        return {}
    elif GUITK=='Wx':
        return {}
    else:
        return {'displayLines':[(None, (), {})],
                'undisplayLines':[(None, (), {})],
                'displayCPK': [(None, (), {})],
                'undisplayCPK':[(None, (), {})],
                'displaySB':[(None, (), {})],
                'undisplaySB':[(None, (), {})],
                'displayBoundGeom':[(None, (), {})],
                'undisplayBoundGeom':[(None, (), {})],
                'showMolecules':[(None, (), {})]
                }
    
commandsInfo = {
    'icoms' : {
        }
    }
