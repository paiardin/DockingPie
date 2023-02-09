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

from Tkinter import _default_root, Tk

if _default_root is None:
	_default_root = Tk()
        _default_root.withdraw()
        
#toglInstallDir = '/tsri/python/sun4SunOS5/lib/ '
#path = _default_root.tk.globalgetvar('auto_path')
#_default_root.tk.globalsetvar('auto_path', (toglInstallDir,) + path )

#_default_root.tk.call('lappend', 'auto_path',
#                      '/mgl/ms1/python/dev/opengltk/OpenGL/Tk/')
#
#_default_root.tk.call('package', 'require', 'Togl')
#
#                      os.path.join( \
#				  os.path.dirname(__file__), \
#				  sys.platform + "-tk" + _default_root.getvar("tk_version")))
