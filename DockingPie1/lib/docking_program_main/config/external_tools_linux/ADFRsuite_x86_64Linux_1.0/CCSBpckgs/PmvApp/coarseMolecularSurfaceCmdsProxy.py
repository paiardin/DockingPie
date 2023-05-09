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
# Author: Michel F. SANNER, Ruth Huey
#
# Copyright: M. Sanner TSRI 2015
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/PmvApp/coarseMolecularSurfaceCmdsProxy.py,v 1.1.4.1 2017/07/13 20:52:52 annao Exp $
#
# $Id: coarseMolecularSurfaceCmdsProxy.py,v 1.1.4.1 2017/07/13 20:52:52 annao Exp $
#
def getGUI(GUITK):
    if GUITK=='Tk':
        from guiTK.computeCoarseMolSurfaceCmds import coarseMolSurfaceGUI
        return {
            'computeCoarseMolecularSurface': [(coarseMolSurfaceGUI, (), {})],
            }
    elif GUITK=='Qt':
        return {}
    elif GUITK=='Wx':
        return {}
    else:
        #from .GUI.QT.coarseMolSurfGUI import coarseMolSurfaceGUI
        return {
            #'computeCoarseMolecularSurface', [(coarseMolSurfaceGUI, (), {})],
            'computeCoarseMolecularSurface':[(None, (), {})],
            }

commandsInfo = {
    'icoms' : {
        }
    }
    
        
