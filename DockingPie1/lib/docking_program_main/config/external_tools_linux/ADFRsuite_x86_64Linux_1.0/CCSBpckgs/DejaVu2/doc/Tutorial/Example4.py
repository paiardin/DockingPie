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

import sys; sys.path.insert(0, '.')

#Example Four: Callbacks and the event Manager
#This example illustrates the event Manager property of the camera.

#First get a handle to the event manager:

from DejaVu import Viewer
MyViewer=Viewer()

ehm=MyViewer.cameras[0].eventManager

#The eventManager is a property of an individual camera. You can
#get a complete listing off all the callback functions bound to keys strokes
#existing in the event Manager at any time:

ehm.ListBindings()

#or the list of all callbacks bound to any specific key stroke:

ehm.ListBindings('<Button-1>')

#Predefined callback functions can be added to specific key strokes in the camera:

def mycallback1(event):
	print 'mycallback1 Event at %d %d' % (event.x, event.y)

def mycallback2(event):
	print 'mycallback2 Event at %d %d' % (event.x, event.y)

ehm.AddCallback("<F1>", mycallback1)

#you can check:
ehm.ListBindings("<F1>") 


#Note that the callback function must have 'event' as its parameter and
#the key stroke must occur with the cursor over the camera. AddCallback
#adds this callback to whatever else may have been previously bound to this
#keystroke. Another function, SetCallback, replaces the previously bound
#funtions with the new one. It returns a list of the previous callbacks
#which can be stored and restored.

funcList = ehm.SetCallback("<F1>", mycallback2)

#Now, funcList is a list:
funcList

#and mycallback2 is bound to F1. mycallback1 could be restored
#as follows:

ehm.SetCallback("<F1>", funcList)

#Callback functions can be removed if they are currently listed:
#if mycallback2 is one of the current callbacks for F1, it can be removed

if ehm.HasCallback("<F1>", mycallback2):
	ehm.RemoveCallback("<F1>", mycallback2)

print "end of Example4"
