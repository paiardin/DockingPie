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

import user
import tkFileDialog, os


def fileOpenAsk(master, idir=None, ifile=None, types=None,
                title='Open'):
    if types==None: types = [ ('All files', '*') ]
    file = tkFileDialog.askopenfilename( filetypes=types,
                                         initialdir=idir,
                                         initialfile=ifile,
                                         title=title)
    if file=='': file = None
    return file
2
def fileSaveAsk(master, idir=None, ifile=None, types = None,
                title='Save'):
    if types==None: types = [ ('All files', '*') ]
    file = tkFileDialog.asksaveasfilename( filetypes=types,
                                           initialdir=idir,
                                           initialfile=ifile,
                                           title=title)
    if file=='': file = None
    return file

from mglutil.gui.BasicWidgets.Tk.dirDialog import askdirectory

class DirOpenBrowser:
    
    def __init__(self, lastDir=None, title=None, parent=None):
        self.lastDir = lastDir
        if lastDir is None:
            self.lastDir = user.home
        self.title = title
        if title is None:
            self.title = 'Choose Directory'
        self.parent = parent
            

    def get(self):
        folder = tkFileDialog.askdirectory(parent = self.parent,
                                           initialdir = self.lastDir,
                                           title=self.title)
        if folder:
            self.lastDir = os.path.split(folder)[0]
            return folder
        else:
            return None


class FileOpenBrowser:
    
    def __init__(self, lastDir=None, title=None, filetypes=None, parent=None):
        self.lastDir = lastDir
        if lastDir is None:
            self.lastDir = user.home
        self.title = title
        if title is None:
            self.title = 'Choose File'
        self.filetypes = filetypes
        self.parent = parent
        if filetypes is None:
            self.filetypes = [('all', '*')]
            

    def get(self):
        file = tkFileDialog.askopenfilename(parent = self.parent,
            initialdir = self.lastDir, filetypes=self.filetypes,
            title=self.title)

        if file:
            self.lastDir = os.path.split(file)[0]
            return file
        else:
            return None
        
class FileSaveBrowser:
    
    def __init__(self, lastDir=None, title=None, filetypes=None):
        self.lastDir = lastDir
        if lastDir is None:
            self.lastDir = user.home
        self.title = title
        if title is None:
            self.title = 'Choose File'
        self.filetypes = filetypes
        if filetypes is None:
            self.filetypes = [('all', '*')]
            

    def get(self):
        file = tkFileDialog.asksaveasfilename(
            initialdir = self.lastDir, filetypes=self.filetypes,
            title=self.title)

        if file:
            self.lastDir = os.path.split(file)[0]
            return file
        else:
            return None
        
