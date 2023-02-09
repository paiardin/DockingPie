#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner and TSRI 2016
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/mslibCOM/__init__.py,v 1.3 2016/12/09 20:11:24 sanner Exp $
#
# $Id: __init__.py,v 1.3 2016/12/09 20:11:24 sanner Exp $
#
# while MGLTools2 is under the GNU Lesser Public License, 
# MSMS is not open source. it is freely avaialble for academic research but
# for commercial research a license is needed.
# Please contact Dr. Michel F. Sanner (sanner@scripps.edu )for a commercial
# license.
# This package is installed for MGLTools2 installation in commercial settings
# without an MSMS license.
#
msg = """MSMS/MSLIB is not enabled in your installation of MGLTools2.\nHence MSMS-based molecular surfaces cannot be calculated.\nMSMS/MSLIB is free for academic but requires a license for commercial use.\nPlease contact Dr. Michel F. Sanner (sanner@scripps.edu) for a license dor commercial use of MSMS/MSLIB.\nOnce you obtain a license key, you can enable the molecular surface features through the help menu in the Pmv application"""

print msg

import sys, os
# if there is a display we display a graphical version of the message
# not sure how to dected the display on windos so we assume that on windows
# we always have a display

if os.name=='nt' or "DISPLAY" in os.environ:
    from PySide import QtGui
    try:
        app = QtGui.QApplication(sys.argv)
        makeApp = True
        
    except RuntimeError:
        makeApp = False
    msgBox = QtGui.QMessageBox()
    msgBox.setText(msg)
    msgBox.exec_()

raise RuntimeError(msg)
