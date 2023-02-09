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

#
# Author Michel F. Sanner  (may 2001)          Copyright M. Sanner, TSRI
#
# $Id: callback.py,v 1.10.18.1 2017/07/26 22:35:43 annao Exp $
#
# $Author: annao $
#
import traceback
import types

class CallbackManager:
    """Class to manage a list of callback functions"""

    def __init__(self):
	self.callbacks = []


    def FindFunctionByName(self, funcName):
	"""find a function with a given name in a list of functions"""

        for f in self.callbacks:
	    if f.__name__==funcName: return f
	return None


    def SetCallback(self, func):
	"""Delete all and set a callback fuction"""

	assert func is None or callable(func)
        if func is None:
            self.callbacks = []
        else:
            self.callbacks = [func, ]


    def AddCallback(self, func):
	"""Add a callback fuction"""

	assert callable(func)
	self.callbacks.append(func)


    def CallCallbacks(self, *args, **kw):
	"""call all callback fuctions"""

        results = []
	for func in self.callbacks:
            try:
                results.append( apply(func, args, kw) )
                
            except:
                print 'ERROR ************************************************'
                traceback.print_exc()
                return 'ERROR'

        return results

            
    def ListCallbacks(self):
        for func in self.callbacks:
            print func.__name__,func


    def RemoveCallback(self, func):
	"""Delete a callback fuction"""
        if type(func)==types.StringType:
            func = self.FindFunctionByName(func)
            if func is None: return "function %s not found"%func
        if func in self.callbacks:
            self.callbacks.remove(func)
        else:
            return "function %s not found"%func.__name__


class CallbackFunction:
    """Class to allow to specify arguments to a callback function"""

    def __init__(self, function, *args, **kw):
        self.function = function
        self.args = args
        self.kw = kw

    def __call__(self, *args, **kw):
        args = self.args + args
        kw.update(self.kw)
        return apply(self.function, args, kw)

CallBackFunction = CallbackFunction
