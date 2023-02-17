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
# Copyright: M. Sanner and TSRI 2015
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/ADFR/bin/runAGFR.py,v 1.5 2017/05/05 00:21:02 annao Exp $
#
# $Id: runAGFR.py,v 1.5 2017/05/05 00:21:02 annao Exp $
#
#
# Usage: pythonsh runADFR.py lig_random.pdbqt mapFolder -o file.log -ref lig_xtal.pdbqt
#
from ADFR.utils.runAGFR import runAGFR
from ADFR.utils.optParser import AGFRParser

parser = AGFRParser()
cnfFile = None
import sys
# check if argument list contains configuration file
args = sys.argv[1:]
if len(args):
    for i, arg in enumerate(args):
        if arg in ['--config', '-F']:
            cnfFile = args[i+1]
            break
if cnfFile:
    # read arguments from the file, pass them to the parser
    from ADFR.utils.optParser import readConfigFile
    cnfArgs = readConfigFile(cnfFile)
    kw =  vars(parser.parse_args(cnfArgs))
else:
    # process command line arguments
    kw =  vars(parser.parse_args())

runner = runAGFR()
runner(**kw)

