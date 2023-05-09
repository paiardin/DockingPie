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
# $Header: /mnt/raid/services/cvs/PmvApp/selectionCmdsProxy.py,v 1.1.4.1 2017/07/13 20:55:28 annao Exp $
#
# $Id: selectionCmdsProxy.py,v 1.1.4.1 2017/07/13 20:55:28 annao Exp $
#

def getGUI(GUITK):
    if GUITK=='Tk':
        from guiTK.selectionCmds import MVClearSelectionGUI, MVSelectSetCommandGUI,\
       MVSelectFromStringCommandGUI,  MVDirectSelectCommandGUI,\
       MVSelectSphericalRegionGUI, MVInvertSelectionGUI,\
       MVSelectNoWaterHeteroAtomsCommandGUI, MVSaveSetCommandGUI
        return {
    'select' : [(None, (), {})],
    'deselect' : [(None, (), {})],
    'clearSelection' : [(MVClearSelectionGUI, (), {})],
    'expandSelection' : [(None, (), {})],
    'selectAround' : [(None, (), {})],
    'saveSet' : [(MVSaveSetCommandGUI, (), {}) ],
    'createSetIfNeeded' : [(None, (), {})],
    'invertSelection' : [(MVInvertSelectionGUI, (), {})],
    'selectSet' : [(MVSelectSetCommandGUI, (), {})],
    'selectFromString' : [(MVSelectFromStringCommandGUI, (), {})],
    'directSelect' : [(MVDirectSelectCommandGUI, (), {})],
    'selectInSphere' : [(MVSelectSphericalRegionGUI, (), {})],
    'selectHeteroAtoms' : [(MVSelectNoWaterHeteroAtomsCommandGUI, (), {})],

}
    elif GUITK=='Qt':
        return {}
    elif GUITK=='Wx':
        return {}
    else:
        return {
            'select' : [(None, (), {})],
            'deselect' : [(None, (), {})],
            'clearSelection' : [(None, (), {})],
            'expandSelection' : [(None, (), {})],
            'selectAround' : [(None, (), {})],
            'saveSet' : [(None, (), {})],
            'createSetIfNeeded' : [(None, (), {})],
            'invertSelection' : [(None, (), {})],
            'selectSet' : [(None, (), {})],
            'selectFromString' : [(None, (), {})],
            'directSelect' : [(None, (), {})],
            'selectInSphere' : [(None, (), {})],
            'selectHeteroAtoms' : [(None, (), {})]
            }

commandsInfo = {
    'icoms' : {
        }
    }
