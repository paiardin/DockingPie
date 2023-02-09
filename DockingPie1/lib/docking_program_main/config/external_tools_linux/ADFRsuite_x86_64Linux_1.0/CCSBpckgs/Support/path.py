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

# This modules sets MGLTools path
# Author: Sargis Dallakyan (sargis at scripps.edu)

import os, sys
from user import home
from Support.version import __version__
rc = home + os.sep + ".mgltools" + os.sep + __version__

if rc:
    updates_rc_dir = rc + os.sep + 'update'
else:
    updates_rc_dir = 'update'
tested_rc = updates_rc_dir + os.sep + 'tested'
nightly_rc = updates_rc_dir + os.sep + 'nightly' #may or may not be tested
release_path = ''
path_text = ''

def setSysPath(path):
    """Sets sys.path for MGLTools"""
    global release_path, path_text
    sys.path.insert(1,path + os.sep + 'PIL')
    release_path = path

    if os.path.exists(tested_rc):
        p = open(tested_rc).read().split('\t')
        if os.path.exists(p[0]):
            if os.listdir(p[0]):
                sys.path.insert(0, p[0])
                path_text = '  Tested '+ p[1]+ '  : ' + p[0] +'\n' + path_text
        else:
            os.remove(tested_rc)
            
    if os.path.exists(nightly_rc):
        p = open(nightly_rc).read().split('\t')
        if os.path.exists(p[0]):
            if os.listdir(p[0]):
                sys.path.insert(0, p[0])
                path_text = ' Nightly '+ p[1]+ '  : ' + p[0] +'\n' + path_text
