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

########################################################################
#
# Date: 2014 Authors: Michel Sanner
#
#    sanner@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Michel Sanner and TSRI
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/DejaVu2/Tk/__init__.py,v 1.1.1.1.4.1 2017/07/13 22:22:07 annao Exp $
#
# $Id: __init__.py,v 1.1.1.1.4.1 2017/07/13 22:22:07 annao Exp $
#

from Tkinter import _default_root, Tk

def loadTogl(master):
    # simulate the setting of TCLLIPATH

    import sys, os
    from os import path
    # Togl is expected to be 

    # build path to directory containing Togl
    from opengltk.OpenGL import Tk
    ToglPath = path.dirname(path.abspath(Tk.__file__))
    # get TCL interpreter auto_path variable
    tclpath = master.tk.globalgetvar('auto_path')

    # ToglPath not already in there, add it
    from string import split
    if ToglPath not in tclpath:
        tclpath = (ToglPath,) + tclpath
        master.tk.globalsetvar('auto_path', tclpath )
    # load Togl extension into TCL interpreter

    #if os.name == 'nt':
    #    toglVersion = master.tk.call('package', 'require', 'Togl','1.7')  
    #else:
    #    toglVersion = master.tk.call('package', 'require', 'Togl','2.1')
    toglVersion = master.tk.call('package', 'require', 'Togl','2.1')

    return toglVersion


