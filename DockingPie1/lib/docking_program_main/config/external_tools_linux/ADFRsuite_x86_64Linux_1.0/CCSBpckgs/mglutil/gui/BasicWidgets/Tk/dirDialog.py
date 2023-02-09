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

#########################################################################
#
# Date: May 2003 Authors: Sophie COON
# 

#########################################################################
# This is a fix to use the askdirectory widget developped by written by
# Fredrik Lundh but only available for Python2.2 and higher.
#########################################################################

from tkCommonDialog import Dialog

class Directory(Dialog):
    "Ask for a directory"

    command = "tk_chooseDirectory"

    def _fixresult(self, widget, result):
        self.widget = widget
        if result:
            # keep directory until next time
            self.options["initialdir"] = result
        self.directory = result # compatibility
        return result

class CreateDirectory(Dialog):
    command = "tk_getSaveFile"
    def _fixresult(self, widget, result):
        if result:
            # keep directory until next time
            self.options["initialdir"] = result
        self.directory = result # compatibility
        return result


def askdirectory (**options):
    "Ask for a directory, and return the file name"
    return apply(Directory, (), options).show()

def createdirectory(**options):
    "Ask for a directory, and return the file name"
    return apply(CreateDirectory, (), options).show()
    
