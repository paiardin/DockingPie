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
# $Header: /mnt/raid/services/cvs/DejaVu2/Qt/EventHandler.py,v 1.1.1.1.4.1 2017/07/13 22:25:45 annao Exp $
#
# $Id: EventHandler.py,v 1.1.1.1.4.1 2017/07/13 22:25:45 annao Exp $
#

from DejaVu2.EventHandler import EventManagerBase

class EventManager(EventManagerBase):
    """Object used to manage callback functions for the events of a Widget

    Public Methods:
    ValideEvent(eventType)
    AddCallback(eventType, callback)
    SetCallback(func)
    RemoveCallback(eventType, func)
    ListBindings(event=None)
    """

# NOT USED, imstead I simply try to bind to a dummy widget to check if
# the given event type is valid
#
#    eventTypes = ('Key', 'KeyPress', 'KeyPress', 
#		  'Button', 'ButtonPress', 'ButtonRelease',
#		  'Enter', 'Leave', 'Motion')
#		  
#    eventModifiers = ('Control' 'Shift', 'Lock', 
#		      'Button1', 'B1', 'Button2', 'B2','Button3', 'B3',
#		      'Button4', 'B4', 'Button5', 'B5',
#		      'Any', 'Double', 'Triple',
#		      'Mod1', 'M1', 'Meta', 'M',
#		      'Mod2', 'M2', 'Alt',
#		      'Mod3', 'M3', 'Mod4', 'M4', 'Mod5', 'M5' )
#    buttonDetails = ( '1', '2', '3' )
#    keyDetails = any keysym

    def __init__(self, widget):
        EventManagerBase.__init__(self, widget)
        

    def ValidEvent(self, eventType):
	"""Check if an event is valid"""

	#try: self.dummyFrame.bind(eventType, self.DummyCallback) 
	#except: return 0
	return 1


    def AddCallback(self, eventType, callback):
	"""Add a callback fuction"""

	EventManagerBase.AddCallback(self, eventType, callback)
	self.widget.bind(eventType, callback, '+')


    def BindFuncList(self,eventType, funcList):
	"""Bind a list of functions to an event"""

	self.widget.bind(eventType, funcList[0])
	for f in funcList[1:]:
	    self.widget.bind(eventType, f, '+')
	    

    def SetCallback(self, eventType, callback):
	"""Set func as the callback or list of callback functions"""

	EventManagerBase.SetCallback(self, eventType, callback)

	if callable(callback):
	    self.widget.bind(eventType, callback)
	elif len(callback)>0:
	    self.BindFuncList(eventType, callback)
	else:
	    raise ValueError('First argument has to be a function or a list of\
functions')
	    
	return funcList


    def RemoveCallback(self, eventType, func):
	"""Delete function func from the list of callbacks for eventType"""

	EventManagerBase.RemoveCallback(self, eventType, func)

	if len(self.eventHandlers[eventType])==0:
	    self.widget.bind(eventType, self.DummyCallback)
	else:
	    self.BindFuncList(eventType, self.eventHandlers[eventType])
	return func


