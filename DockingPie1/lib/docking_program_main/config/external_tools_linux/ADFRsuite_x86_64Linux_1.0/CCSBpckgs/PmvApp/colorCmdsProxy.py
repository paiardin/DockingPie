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
# $Header: /mnt/raid/services/cvs/PmvApp/colorCmdsProxy.py,v 1.1.4.1 2017/07/13 20:55:28 annao Exp $
#
# $Id: colorCmdsProxy.py,v 1.1.4.1 2017/07/13 20:55:28 annao Exp $
#

def getGUI(GUITK):
    print "GUITK:", GUITK
    if GUITK=='Tk':
        from guiTK.colorCmds import ColorCommandGUI, ColorByAtomTypeGUI,\
             ColorByDGGUI, ColorByResidueTypeGUI, ColorShapelyGUI, \
             ColorByChainGUI, ColorByMoleculeGUI, ColorByInstanceGUI,\
             ColorByPropertiesGUI, RainbowColorGUI, RainbowColorByChainGUI,\
             ColorByExpressionGUI, ColorByLineGeometryColorGUI
        return {
            'color':[(ColorCommandGUI, (), {})],
            'colorByAtomType' :[(ColorByAtomTypeGUI, (), {})],
            'colorByResidueType' :[(ColorByResidueTypeGUI, (), {})],
            'colorAtomsUsingDG':[(ColorByDGGUI, (), {})],
            'colorResiduesUsingShapely' :[(ColorShapelyGUI, (), {})] ,
            'colorByChains' :[(ColorByChainGUI, (), {})],
            'colorByMolecules':[(ColorByMoleculeGUI, (), {})],
            'colorByInstance' :[(ColorByInstanceGUI, (), {})],
            'colorByProperty': [(ColorByPropertiesGUI, (), {})],
            'colorRainbow' : [(RainbowColorGUI, (), {})],
            'colorRainbowByChain' : [(RainbowColorByChainGUI,(), {})],
            'colorByExpression' : [(ColorByExpressionGUI, (), {})],
            'colorByLinesColor': [(ColorByLineGeometryColorGUI, (), {})]
            }
    elif GUITK=='Qt':
        return {}
    elif GUITK=='Wx':
        return {}
    else:
        return {'color':[(None, (), {})],
                'colorByAtomType' :[(None, (), {})],
                'colorByMolecules':[(None, (), {})],   
                'colorByChains' :[(None, (), {})],
                'colorRainbow' : [(None, (), {})],
                'customColor': [(None, (), {})],
                'setOpacity':  [(None, (), {})],
                'moleculesRainbow': [(None, (), {})],
                'colorBySS': [(None, (), {})],
                'colorRainbowChain' : [(None, (), {})],
                #'colorByResidueType' :[(None, (), {})],
                #'colorAtomsUsingDG':[(None, (), {})],
                #'colorResiduesUsingShapely' :[(None, (), {})],
                #'colorByInstance' :[(None, (), {})],
                #'colorByProperty': [(None, (), {})],
                #'colorByExpression' : [(None, (), {})],
                #'colorByLinesColor': [(None, (), {})],
                            }

commandsInfo = {
    'icoms' : {
        }
    }
