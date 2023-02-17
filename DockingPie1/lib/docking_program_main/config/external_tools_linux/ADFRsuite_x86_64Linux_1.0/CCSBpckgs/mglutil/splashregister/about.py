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

# $Header: /mnt/raid/services/cvs/python/packages/share1.5/mglutil/splashregister/about.py,v 1.10.10.1 2017/07/26 22:35:43 annao Exp $
# $Id: about.py,v 1.10.10.1 2017/07/26 22:35:43 annao Exp $
import Tkinter, os
from mglutil.util.misc import ensureFontCase

try:
    from PIL import Image, ImageTk
except:
    pass

        
class About:
    """
    package    : mglutil
    module     : splashregister.about
    class      : About
    description:
        Displays information needed for About widget
    """
    def __init__(self, title=None, image_dir='.', version='', revision=None, authors=None,
                 icon=None, copyright=None, third_party='', path_data=''):
        self.title = title
        self.image_dir = image_dir
        self.version = version
        self.revision = revision
        self.authors = authors
        self.icon = icon
        self.copyright = copyright    
        self.third_party = third_party
        self.path_data = path_data
    
    def gui(self, master):
        
        Tkinter.Label(master, text=self.title, font =(ensureFontCase('helvetica'), 16, 'bold') ).\
                      pack(side='top')
        text = 'Version ' + self.version
        if self.revision is not None:
            text += '  ' + self.revision
        night = self.path_data.find(" Nightly ")
        if night != -1:
            tmpTxt = self.path_data[night:].split()
            text += " - Update Nightly Build " + tmpTxt[1]
        else:
            tested = self.path_data.find(" Tested ")
            if tested != -1:
                tmpTxt = self.path_data[tested:].split()
                text += " - Update Tested Build " + tmpTxt[1]
            
        Tkinter.Label(master, text=text).pack(side='top')

        files = os.listdir(self.image_dir)
        import fnmatch
        files = fnmatch.filter(files,'*.jpg') + fnmatch.filter(files,'*.png')
        import random
        rand = random.randint(0,len(files)-1)
        image_file = os.path.join(os.path.join(self.image_dir ,files[rand]))
        image = Image.open(image_file)
        self.image1 = ImageTk.PhotoImage(image, master=master)
        self.imageTk = Tkinter.Label(master,image=self.image1 )
        self.imageTk.pack()
        Tkinter.Label(master, text=self.copyright, relief='sunken' ).pack()
        logoFrame = Tkinter.Frame(master, bg='white')
        logoFrame.pack(fill='x',expand=True)
        basepath = os.path.split(__file__)[0]
        
        NBCR = Image.open(os.path.join(basepath,'NBCR.jpg'))
        self.NBCR1 = ImageTk.PhotoImage(NBCR, master=master)
        self.NBCRTk = Tkinter.Label(logoFrame,image=self.NBCR1, bd=0 )
        self.NBCRTk.pack(side='left', padx=40, expand=True)
        NIH = Image.open(os.path.join(basepath,'NIH.gif'))
        self.NIH1 = ImageTk.PhotoImage(NIH, master=master)
        self.NIHTk = Tkinter.Label(logoFrame,image=self.NIH1, bd=0)
        self.NIHTk.pack(side='left', padx=40,expand=True)
        NSF = Image.open(os.path.join(basepath,'NSF.gif'))
        self.NSF1 = ImageTk.PhotoImage(NSF, master=master)
        self.NSFTk = Tkinter.Label(logoFrame,image=self.NSF1, bd=0)
        self.NSFTk.pack(side='left', padx=40, expand=True)

if __name__ == '__main__':
    root = Tkinter.Tk()
    about = About(image_dir='../../Pmv/Icons/Images')
    about.gui(root)
    root.mainloop()