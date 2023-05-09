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
# Copyright: M. Sanner and TSRI 2017
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/ADFR/GUI/statusBar.py,v 1.4 2017/05/19 01:34:00 sanner Exp $
#
# $Id: statusBar.py,v 1.4 2017/05/19 01:34:00 sanner Exp $
#
import sys
from PySide import QtCore, QtGui
from mglutil.util.callback import CallbackFunction

class MyStatusBar(QtGui.QStatusBar):

    _levelIcons = {}
    
    def __init__(self, parent=None):
        super(MyStatusBar, self).__init__(parent=parent)
        # make label appear flat by default
        self.setStyleSheet("QStatusBar::item {border: 0;}")

        self._okPixmap = QtGui.QPixmap("ADFR/GUI/Icons/ok24.png")
        self._notokPixmap = QtGui.QPixmap("ADFR/GUI/Icons/stop24.png")
        self._warningPixmap = QtGui.QPixmap("ADFR/GUI/Icons/warning24.png")
        self._nothingPixmap = QtGui.QPixmap("ADFR/GUI/Icons/nothing24.png")

        self._busy = QtGui.QMovie("ADFR/GUI/Icons/animatedGears2_24.gif")
        self._busy.setCacheMode(QtGui.QMovie.CacheAll) 
        self._busy.setSpeed(100) 
        #print 'supported formats', self._busy.supportedFormats()

        self._levelIcons['ok'] = self._okPixmap
        self._levelIcons['error'] = self._notokPixmap
        self._levelIcons['warning'] = self._warningPixmap
        self._levelIcons['nothing'] = self._nothingPixmap
        self._levelIcons['busy'] = self._busy
        self.timer = QtCore.QTimer(self)
        self._currentLevel = None
        self._animated = False
        self._animatedIcons = ['busy']
        
        # add icon on the left for icon
        self._statusIcon = QtGui.QLabel("statusIcon", self)
        self._statusIcon.setPixmap(self._nothingPixmap)
        self._statusIcon.setFrameStyle(QtGui.QFrame.NoFrame)
        self._statusIcon.setStyleSheet("QLabel {border: 0;}")
        self.addWidget(self._statusIcon)
        #p = self._statusIcon.parent()
        
        # add stretchable text box
        self._message = QtGui.QLabel("Welcome")
        self._message.setFrameStyle(QtGui.QFrame.NoFrame)
        self.addPermanentWidget(self._message, 1)

    def showMessage(self, msg, level='nothing', timeout=0):
        # update message level icon
        #import traceback; traceback.print_stack()
        #print 'LOG', msg, level, timeout, self._message.text(), self._currentLevel
        if self._currentLevel in self._animatedIcons:
            self._levelIcons[self._currentLevel].stop()
        if level in self._animatedIcons:
            movie = self._levelIcons[level]
            self._statusIcon.setMovie(movie)
            movie.start()
        else:
            pixmap = self._levelIcons.get(level, self._levelIcons['nothing'])
            self._statusIcon.setPixmap(pixmap)

        # register timer if needed
        if not self.timer.isActive() and timeout>0:
            # save previous message
            curMsg = self._message.text()
            #print 'CB', curMsg, self._currentLevel, timeout
            cb = CallbackFunction( self.showMessage, curMsg, self._currentLevel)
            self.timer.singleShot(timeout, cb)
            self._currentLevel = level

        # set message
        self._message.setText(msg)

if __name__=='__main__':       
    import os
    app = QtGui.QApplication(sys.argv)
    qtpath = os.path.join(os.getenv('MGL_ROOT'), 'plugins')
    app.addLibraryPath(qtpath)#'/mgl/ms1/people/sanner/python/mglTools2/MGLTools2-1.1/plugins')
    bar = MyStatusBar()

    label = QtGui.QLabel("Box")
    label.setFrameStyle(QtGui.QFrame.StyledPanel | QtGui.QFrame.Sunken)
    sample_palette = QtGui.QPalette() 
    sample_palette.setColor(QtGui.QPalette.Window, QtCore.Qt.red);
    label.setAutoFillBackground(True)
    label.setPalette(sample_palette)
    bar.addPermanentWidget(label)

    label1 = QtGui.QLabel("TP")
    label1.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Plain)
    sample_palette = QtGui.QPalette() 
    sample_palette.setColor(QtGui.QPalette.Window, QtCore.Qt.gray);
    label1.setAutoFillBackground(True)
    label1.setPalette(sample_palette)
    label1.setDisabled(True)
    bar.addPermanentWidget (label1)

    label2 = QtGui.QLabel("FR")
    label2.setFrameStyle(QtGui.QFrame.StyledPanel | QtGui.QFrame.Sunken)
    sample_palette = QtGui.QPalette() 
    sample_palette.setColor(QtGui.QPalette.Window, QtCore.Qt.green);
    label2.setAutoFillBackground(True)
    label2.setPalette(sample_palette)
    bar.addPermanentWidget (label2)

    #bar.showMessage('Ready', level='warning')
    bar.showMessage('Ready')
    #bar.showMessage('there was an error', 'error', 1000)
    bar.showMessage('I am busy', 'busy', 0)
    #bar.showMessage('there was NO error', 'ok', 1000)

    bar.show()
    sys.exit(app.exec_())
