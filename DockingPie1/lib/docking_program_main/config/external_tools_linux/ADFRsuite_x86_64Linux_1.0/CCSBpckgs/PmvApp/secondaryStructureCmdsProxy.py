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
# $Header: /mnt/raid/services/cvs/PmvApp/secondaryStructureCmdsProxy.py,v 1.1.4.1 2017/07/13 20:55:28 annao Exp $
#
# $Id: secondaryStructureCmdsProxy.py,v 1.1.4.1 2017/07/13 20:55:28 annao Exp $
#

def getGUI(GUITK):
    if GUITK=='Tk':
        from guiTK.secondaryStructureCmds import ComputeSecondaryStructureCommandGUI, \
             ExtrudeSecondaryStructureCommandGUI, \
             ExtrudeSecondaryStructureCommandUnicGUI, RibbonCommandGUI,\
             ColorBySSElementTypeGUI,  DisplayExtrudedSSCommandGUI
        return {
    'computeSecondaryStructure' : [(ComputeSecondaryStructureCommandGUI, (), {})],
    'extrudeSecondaryStructure' : [(ExtrudeSecondaryStructureCommandGUI, (), {})],
    'extrudeSecondaryStructureUnic' : [(ExtrudeSecondaryStructureCommandUnicGUI, (), {})],
    'displayExtrudedSS' : [(DisplayExtrudedSSCommandGUI, (), {})],
    'colorBySecondaryStructure' : [(ColorBySSElementTypeGUI, (), {})],
    'undisplayExtrudedSS' : [(None, (), {})],
    'ribbon' : [(RibbonCommandGUI, (), {})],
            }
    elif GUITK=='Qt':
        return {}
    elif GUITK=='Wx':
        return {}
    else:
        return {'computeSecondaryStructure' : [(None, (), {})],
    'extrudeSecondaryStructure' : [(None, (), {})],
    'extrudeSecondaryStructureUnic' : [(None, (), {})],
    'displayExtrudedSS' : [(None, (), {})],
    'colorBySecondaryStructure' : [(None, (), {})],
    'undisplayExtrudedSS' : [(None, (), {})],
    'ribbon' : [(None, (), {})]
                }

commandsInfo = {
    'icoms' : {
        }
    }
