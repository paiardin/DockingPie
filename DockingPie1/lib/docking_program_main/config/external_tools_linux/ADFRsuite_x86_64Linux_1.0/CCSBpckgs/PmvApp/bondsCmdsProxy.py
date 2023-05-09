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
# Author: Michel F. SANNER,  Anna Omelchenko
#
# Copyright: M. Sanner TSRI 2014
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/PmvApp/bondsCmdsProxy.py,v 1.1.4.1 2017/07/13 20:55:28 annao Exp $
#
# $Id: bondsCmdsProxy.py,v 1.1.4.1 2017/07/13 20:55:28 annao Exp $
#
#


def getGUI(GUITK):
    if GUITK=='Tk':
        from Pmv.guiTK.bondsCmds import BuildBondsByDistanceGUI,\
             AddBondsGUICommandGUI, RemoveBondsGUICommandGUI
        return {#'buildBondsByDistance' : [(BuildBondsByDistanceGUI, (), {})],
                'buildBondsByDistance' : [(None , (), {})],
                'addBonds' : [(None , (), {})],
                'addBondsGC' : [(AddBondsGUICommandGUI, (), {})],
                'removeBonds' : [(None , (), {})],
                'removeBondsGC' : [(RemoveBondsGUICommandGUI, (), {})]
            }
    elif GUITK=='Qt':
        return {}
    elif GUITK=='Wx':
        return {}
    else:
        return {'buildBondsByDistance' : [(None , (), {})],
                'addBonds' : [(None , (), {})],
                'addBondsGC' :  [(None , (), {})],
                'removeBonds' : [(None , (), {})],
                'removeBondsGC' : [(None , (), {})]}
    
commandsInfo = {
    'icoms' : {}
    }
