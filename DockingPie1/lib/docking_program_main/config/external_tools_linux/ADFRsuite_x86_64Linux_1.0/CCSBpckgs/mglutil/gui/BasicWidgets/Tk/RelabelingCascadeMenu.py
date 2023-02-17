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
# Date: Novembre 2005 Authors: Guillaume Vareille, Michel Sanner
#
#    vareille@scripps.edu
#    sanner@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Guillaume Vareille, Michel Sanner and TSRI
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/python/packages/share1.5/mglutil/gui/BasicWidgets/Tk/RelabelingCascadeMenu.py,v 1.4.22.1 2017/07/26 22:35:41 annao Exp $
#
# $Id: RelabelingCascadeMenu.py,v 1.4.22.1 2017/07/26 22:35:41 annao Exp $
#

import Tkinter
import weakref

class RelabelingCascadeMenu(Tkinter.Menu):
    """
"""
    def __init__(self, label, variable, master=None, cnf={}, **kw):
        #print "RelabelingCascadeMenu.__init__", cnf, kw
        Tkinter.Menu.__init__(self,  *(master, cnf), **kw)
        self.baseLabel = label
        self.cascadeVariable = variable
        self._upperMenu = weakref.ref(master)
        self._cascadeMenuIndex = None
        self._externalCallbackFunction = None
        self._valuesLabels = {}

        if hasattr(self._upperMenu(), 'relabelingCascadeMenus') is False:
            self._upperMenu().relabelingCascadeMenus = {}
        self._upperMenu().relabelingCascadeMenus[self.baseLabel] = self


    def add_radiobutton(self, cnf={}, **kw):
        #print "add_radiobutton", cnf, kw
        if self._cascadeMenuIndex is None:
            self._cascadeMenuIndex = self._upperMenu().index(self.baseLabel)

        self._valuesLabels[kw['value']] = kw['label']

        if kw.has_key('command') and kw['command'] is not None:
            self._externalCallbackFunction = kw['command']
        kw['command'] = self._envelopeCallbackFunction
        Tkinter.Menu.add_radiobutton( self, *(cnf,), **kw)


    def setWithoutCallbackFunction(self, value):
        #print "setWithoutCallbackFunction", self.baseLabel, value
        self.cascadeVariable.set(value)
        self._relabelCascade()


    def setWithCallbackFunction(self, value=None):
        #print "setWithCallbackFunction", value
        if value is None:
            value = self.cascadeVariable.get()
        lLabel = self._valuesLabels[value]
        if lLabel == 'none': # because self.index can't deal with label 'none'
            lIndex = 0 # be sure to put label 'none' in index 0 !!!!
        else: # self.index also have problems with numeric label such as '1'
              # use ' 1' instead
            lIndex = self.index(lLabel)
        self.invoke(lIndex)


    def _relabelCascade(self):
        #print "_relabelCascade"
        self._upperMenu().entryconfigure(
          self._cascadeMenuIndex,
          label=self.baseLabel \
                +' [ ' \
                +str(self._valuesLabels[self.cascadeVariable.get()])+' ]')


    def _envelopeCallbackFunction(self):
        #print "_envelopeCallbackFunction"
        self._relabelCascade()
        if self._externalCallbackFunction is not None:
            self._externalCallbackFunction()




