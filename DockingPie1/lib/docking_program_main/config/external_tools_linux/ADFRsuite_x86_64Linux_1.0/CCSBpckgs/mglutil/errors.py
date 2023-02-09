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
# Authors: Michel F. SANNER, 
#
# Copyright: M. Sanner TSRI 2013
#
#############################################################################
# $Header: /mnt/raid/services/cvs/python/packages/share1.5/mglutil/errors.py,v 1.5.6.1 2017/07/26 22:35:40 annao Exp $
#
# $Id: errors.py,v 1.5.6.1 2017/07/26 22:35:40 annao Exp $
#

from mglutil.events import Event
import traceback

_debug = False # set to True to avoid displayErrors from being blocking
_hasDisplayed = False # gets set to True when errors display while _debug is True


def hasDisplayed():
    return _hasDisplayed

def setHasDisplayed(val):
    global _hasDisplayed
    assert val in (True, False, 1, 0)
    _hasDisplayed = val

def resetHasDisplayed():
    global _hasDisplayed
    _hasDisplayed = False

def setDebug(val):
    global _debug
    assert val in (True, False, 1, 0)
    _debug = val
    
class MGLException(Exception):
    pass

class ErrorEvent(Event):
    pass


"""
Error management component
"""

def getErrors(errorManager=None):
    if errorManager is None:
        errorManager = _ErrorManager
    return errorManager.errorStack


def hasErrors(errorManager=None):
    if errorManager is None:
        errorManager = _ErrorManager
    return len(errorManager.errorStack)>0

    
def clearErrors(errorManager=None):
    if errorManager is None:
        errorManager = _ErrorManager
    errorManager.clear()

    
def displayErrors(errorManager=None):
    if errorManager is None:
        errorManager = _ErrorManager
    errorManager.displayErrors()

    
def MGLError(exception, msg, errorManager=None, raiseException=False, **kw):
    if errorManager is None:
        errorManager = _ErrorManager
    #stack = traceback.format_stack()[:-2]
    formatedException = traceback.format_exc().split('\n')
    error = Error(exception, msg, formatedException, errorManager.appName, **kw)
    errorManager.errorStack.append(error)
    if raiseException:
        raise MGLException
    return error
    
class Error:

    def __init__(self, exception, msg, formatedException, appName, **kw):
        self.exception = exception
        self.msg = msg
        self.kw = kw
        self.appName = appName
        self.formatedException = formatedException

    def getMsg(self, indent):
        lines = []
        lines.append(indent+'Message: ')
        nbc = 0
        line = indent+"  "
        for oneLine in self.msg.split('\n'):
            for word in oneLine.split():
                l = len(word)
                if nbc+l >= 75:
                    lines.append(line)
                    line = indent+"  "
                    nbc = 0
                line += word+' '
                nbc += l+1
            if line != indent+"  ":
                lines.append(line)
                nbc = 0
                line = indent+"  "
        if line != indent+"  ":
            lines.append(line)

        return lines

    def printMsg(self, indent):
        for l in self.getMsg(indent): print l


    def getException(self, indent):
        lines = []
        lines.append(indent+'Exception: ')
        msg = self.exception.message
        if not msg:
            if hasattr(self.exception, 'msg'): # assertion errors have no .msg
                msg = self.exception.msg
            else: return []
        lines.append(indent+'  '+msg)
        return lines

    def printException(self, indent):
        for l in self.getException(indent): print l

    def getFormatedException(self, indent):
        lines = []
        lines.append(indent+'Exception: ')
        return lines + [indent+l for l in self.formatedException]

    def printFormatedException(self, indent):
        for l in self.getFormatedException(indent): print l

    def getContext(self, indent):
        lines = []
        if len(self.kw):
            lines.append("  -----------------------------------------------------------------------------")
            lines.append(indent+'  Context:')
            nbc = max([len(k) for k in self.kw.keys()])
            fmt = "  %"+"%-ds"%(nbc+1)
            for k,v in self.kw.items():
                lines.append(indent+fmt%k+":"+str(v))
        return lines

    def printContext(self, indent):
        for l in self.getContext(indent): print l

    def getAll(self, indent, details=False):
        lines = []
        lines.append( "*******************************************************************************")
        lines.append( 'ERROR in %s'%self.appName)
        lines.append( "  -----------------------------------------------------------------------------")
        lines.extend(self.getMsg(indent="  "))
        lines.append( "  -----------------------------------------------------------------------------")
        lines.extend(self.getException(indent="  "))
        lines.extend(self.getContext(indent="  "))
        if details:
            lines.append( "  -----------------------------------------------------------------------------")
            lines.extend(self.getFormatedException(indent="  "))
        lines.append( "*******************************************************************************")
        return lines
                      
    def printAll(self, indent=''):
        for l in self.getAll(indent): print l



class ErrorManager:

    """
    Defines an object for handling errors
    """

    def __init__(self, appName):
        """
        Create an error manager for an application
        """
        self.appName = appName # name of the application for which the object
                               # is created
        self.Status = 'OKAY'
        self.errorStack = []
        self.handleError = None # function to be called when an error occurs
        
        self.showDetails = False
        
    def setAppName(self, name):
        self.appName = name

        
    def error(self, exception, msg, **kw):
        """
        Handle an error
        """
        error = Error(exception, msg, self.appName, **kw)
        self.errorStack.append(error)
        if self.handleError:
            self.handleError(error)
        else:
            error.printAll()

    def clear(self):
        self.errorStack = []


    def displayErrors(self, master=None):

        if len(self.errorStack)==0: return
        # fixme: we should show all errors or an interface to go betwen erros
        errMsg = self.errorStack[0].getAll(indent='', details=self.showDetails)
        self.currentlyDisplayedError = 0
        print errMsg

    def execute(self, result):
        if result == 'OK':
            self.dialog.deactivate(result)
            setHasDisplayed(False)
        elif result=='Next Error':
            if self.currentlyDisplayedError<len(self.errorStack)-1:
                self.currentlyDisplayedError += 1
                errMsg = self.errorStack[self.currentlyDisplayedError].getAll(indent='', details=self.showDetails)
                self.st.delete('0.0', 'end')
                for line in errMsg:
                    self.st.insert('end', line+'\n')
        elif result=='Previous Error':
            if self.currentlyDisplayedError>0:
                self.currentlyDisplayedError -= 1
                errMsg = self.errorStack[self.currentlyDisplayedError].getAll(indent='', details=self.showDetails)
                self.st.delete('0.0', 'end')
                for line in errMsg:
                    self.st.insert('end', line+'\n')
        elif result=='Details >>':
            self.showDetails = not self.showDetails
            errMsg = self.errorStack[self.currentlyDisplayedError].getAll(indent='', details=self.showDetails)
            self.st.delete('0.0', 'end')
            for line in errMsg:
                self.st.insert('end', line+'\n')
            
try:
    _ErrorManager
except:
    _ErrorManager = ErrorManager('MGLTools')


def MGLMessage(type, text, messageManager=None):
    if messageManager is None:
        messageManager = _MessageManager
    messageManager.message(type, text)

    
def addMessageCallback(func, messageManager=None):
    if messageManager is None:
        messageManager = _MessageManager
    messageManager.addCallback(func)

def deleteMessageCallback(func, messageManager=None):
    if messageManager is None:
        messageManager = _MessageManager
    messageManager.deleteCallback(func)


class MessageManager:

    """
    Defines an object for handling messages
    """

    def __init__(self, appName):
        """
        Create an error manager for an application
        """
        self.appName = appName # name of the application for which the object
                               # is created
        self.callbacks = []

    def message(self, type, text):
        #print 'OK'
        for func in self.callbacks:
            #print func, type, text
            func(type, text)

    def addCallback(self, func):
        assert callable(func)
        self.callbacks.append(func)

    def deleteCallback(self, func):
        if func in self.callbacks:
            self.callbacks.remove(func)

try:
    _MessageManager
except:
    _MessageManager = MessageManager('MGLTools')



if __name__=='__main__':
    #pm = ErrorManager('test')
    kw = {'value':35, 'line': 'Hello World 35', 'file':'myTestFile'}
    msg = "Unexpected value while parsing file"
    try:
        1/0.0
    except Exception, e:
        #print dir(e)
        #print e.args
        #print e.message
        #stack = traceback.format_stack()
        #pm.error(e, msg, stack, **kw)
        MGLError(e, msg, **kw)
        
    print
    print
    msg = "Let's make this error message very long to test if handling messages spanning multiple lines works effectively and the words get broken in the right place and it all look good.\nAnd lets throw in a new line to see what happens"
    l = range(24)
    try:
        l[100] = 3
    except Exception, e:
        #stack = traceback.format_stack()
        #pm.error(e, msg, stack)
        MGLError(e, msg)

    _ErrorManager.displayErrors()
