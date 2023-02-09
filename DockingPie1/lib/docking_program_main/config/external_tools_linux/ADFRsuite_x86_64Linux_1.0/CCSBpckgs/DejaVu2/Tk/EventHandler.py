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
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/DejaVu2/Tk/EventHandler.py,v 1.1.1.1.4.1 2017/07/13 22:22:07 annao Exp $
#
# $Id: EventHandler.py,v 1.1.1.1.4.1 2017/07/13 22:22:07 annao Exp $
#

from Tkinter import Frame

class EventManager:
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

       # keys are Tk events, values are lists of callback functions
	self.eventHandlers = {}

	self.widget = widget

        # create a dummy frame to try to bind event to check for event validity
        self.dummyFrame = Frame(widget,width=1, height=1)


    def DummyCallback(self, event):
	"""dummy function used to check event validity"""
	pass


    def ValideEvent(self, eventType):
	"""Check if an event is valide"""

	try: self.dummyFrame.bind(eventType, self.DummyCallback) 
	except: return 0
	return 1


    def AddCallback(self, eventType, callback):
	"""Add a callback fuction"""

	assert type(eventType) == type('string')
	assert callable(callback)
	# Here we should also check that callback has exactly 1 argument

	if not self.ValideEvent(eventType):
	    raise ValueError('%s is not a valide event type' % eventType)

	if eventType in self.eventHandlers.keys():
	    self.eventHandlers[eventType].append(callback)
	else:
	    self.eventHandlers[eventType] = [callback,]

	self.widget.bind(eventType, callback, '+')


    def BindFuncList(self,eventType, funcList):
	"""Bind a list of functions to an event"""

	self.widget.bind(eventType, funcList[0])
	for f in funcList[1:]:
	    self.widget.bind(eventType, f, '+')


    def HasCallback(self, eventType, callback):
	"""Check whether a function is registered as a callback for an event
	"""

	assert callable(callback)
	if self.eventHandlers.has_key(eventType):
	    for f in self.eventHandlers[eventType]:
	    	if f==callback: return 1
	return 0
	    

    def SetCallback(self, eventType, callback):
	"""Set func as the callback or list of callback functions"""

	assert type(eventType) == type('string')

	if self.eventHandlers.has_key(eventType):
	    funcList = self.eventHandlers[eventType]
	else: funcList = None

	if callable(callback):
	    self.eventHandlers[eventType] = [callback,]
	    self.widget.bind(eventType, callback)
	elif len(callback)>0:
	    self.eventHandlers[eventType] = callback
	    self.BindFuncList(eventType, callback)
	else:
	    raise ValueError('First argument has to be a function or a list of\
functions')
	    
	return funcList


    def FindFunctionByName(self, funcName, funcList):
	"""find a function with a given name in a list of functions"""

	for f in funcList:
	    if f.__name__==funcName: return f
	return None
    
    
    def RemoveCallback(self, eventType, func):
	"""Delete function func from the list of callbacks for eventType"""

	if not self.eventHandlers.has_key(eventType):
            return None
#	    raise ValueError('Widget %s has no event %s registered' % \
#			     self.widget, eventType)

	if type(func)==type('string'):
	    func = self.FindFunctionByName(func, self.eventHandlers[eventType])
	    if not func: return

        try:
            self.eventHandlers[eventType].remove(func)
        except:
            pass
	if len(self.eventHandlers[eventType])==0:
	    del self.eventHandlers[eventType]
	    self.widget.bind(eventType, self.DummyCallback)
	else:
	    self.BindFuncList(eventType, self.eventHandlers[eventType])
	return func


    def ListOneBinding(self, key):
	"""List all bindings for one events"""

	print 'Event', key
	if self.eventHandlers.has_key(key):
	    for f in self.eventHandlers[key]:
		print '\t %s' % f.__name__
	else:
	    print '\tNone'


    def ListBindings(self, event=None):
	"""List all bindings for one or all events"""

	if event is None:
	    for k in self.eventHandlers.keys():
		self.ListOneBinding(k)
	else:
	    self.ListOneBinding(event)



class CallbackFunctions:
    """Base class for objects with callback functions"""

    def __init__(self):

	self.callbacks = []


    def Callbacks(self):
	"""Call all callback functions, to be implemented by sub class"""
	pass
#	for f in self.callbacks: f(self.val)


    def SetCallback(self, func):
	"""Delete all and set a callback fuction"""

	assert callable(func)
	self.callbacks = [func, ]


    def AddCallback(self, func):
	"""Add a callback fuction"""

	assert callable(func)
	self.callbacks.append(func)


    def RemoveCallback(self, func):
	"""Delete a callback fuction"""

	self.callbacks.remove(func)


    def MouseUp(self, event):
	"""Call callbak function for non immediate sliders"""
	self.Callbacks()

