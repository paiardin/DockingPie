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

# $Header: /mnt/raid/services/cvs/python/packages/share1.5/mglutil/splashregister/license.py,v 1.15.12.1 2017/07/26 22:35:43 annao Exp $
# $Id: license.py,v 1.15.12.1 2017/07/26 22:35:43 annao Exp $
#

import Tkinter

tk_root = Tkinter.Tk()
tk_root.title("Commercial Usage")
txt = """
 The software component for computing molecular surfaces (MSMS) 
 is not free for commercial usage. If you plan to use MSMS for commercial 
 research please contact sanner@scripps.edu

 Some software components such at the volume rendering and 
 isocontouring were developed at UT Austin.

 If you publish scientific results generated using this software 
 please cite the appropriate software components.
 A list of papers is provided under Help -> Citation Information 
 menu in PMV and ADT. 
"""
Tkinter.Label(tk_root, text=txt, justify=Tkinter.LEFT).pack()
Tkinter.Button(tk_root, text="OK", command=tk_root.quit).pack()
tk_root.mainloop()
