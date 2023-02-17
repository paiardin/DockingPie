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
## Copyright (c) MGL TSRI 2016
##
################################################################################

#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner and TSRI 2016
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/ADFRccDIST/ADFRcc/__init__.py,v 1.5 2016/12/07 01:08:39 sanner Exp $
#
# $Id: __init__.py,v 1.5 2016/12/07 01:08:39 sanner Exp $
#
import ADFRcc, os
from ADFRcc.adfr import Parameters

#set log level to ERROR to avoid too many messages printed
from ADFRcc import _adfrcc
_adfrcc.Logger_setLogLevelError()

# create a parameter object to be used
_parameters = Parameters.getParameters()

def setForceFieldVersion(version):
    """
    Set the forcefield parameters using a version number (string) and return
    a handle to the object holding them

    current version are:
    '4.2': emulating AutoDock4.2
    'default': future AutoDock5
    
    Parameters <- setForceFieldVersion(version)
    """
    _parameters.clearAll()
    if version=='default':
        _parameters.loadDefaults()
    elif version=='4.2':
        _parameters.loadFromDatFile(os.path.join(ADFRcc.__path__[0], 'Data',
                                                'AD42PARAM.DAT'))
    else:
        raise ValueError("ERROR: unknown forcefield version %s"%version)
    return _parameters

def setForceFieldFromFile(filename):
    """
    Set the forcefield parameters from a file and return a handle to the object
    holding them
    
    Parameters <- setForceFieldFromFile(filename)
    """
    _parameters.clearAll()
    _parameters.loadFromDatFile(version)
    return _parameters

def getFFParameters():
    """
    get a handle to the forcefield parameter object
    
    Parameters <- getParameters()
    """
    return _parameters
