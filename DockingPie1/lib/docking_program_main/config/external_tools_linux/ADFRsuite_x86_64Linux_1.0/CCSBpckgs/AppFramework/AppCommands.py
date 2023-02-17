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
# Copyright: M. Sanner TSRI 2013
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/AppFramework/AppCommands.py,v 1.11.4.1 2017/07/28 00:58:00 annao Exp $
#
# $Id: AppCommands.py,v 1.11.4.1 2017/07/28 00:58:00 annao Exp $
#

import weakref, sys, numpy, traceback
from time import time

###########################################################################
##
##   EVENTS
##
###########################################################################

from mglutil.events import Event

###########################################################################
##
##   ERROR ANSD MESSAGES
##
###########################################################################
class AppSuccess:
    """Object used to describe a successful operation
    """
    def __init__(self, time, msg, obj, **kw):
        self.time = time
        self.msg = msg
        self.obj = obj
        self.kw = kw
	self.exception = None
	self.formatedException = None
	
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
	if self.exception is None: return lines
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

    def getFormatedException(self, indent=""):
        lines = []
	if self.formatedException is None: return lines
        lines.append(indent+'Exception: ')
        return lines + [indent+l for l in self.formatedException]

    def printFormatedException(self, indent):
        for l in self.getFormatedException(indent): print l

    def getContext(self, indent=''):
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

    def getAll(self, app, indent, details=False):
        lines = []
        lines.append( "*******************************************************************************")
        lines.append( 'ERROR in %s'%app.name)
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
    
class AppWarning(AppSuccess):
    """Object used to describe a warning encountered during the execution of a command

    Warning should be issued if something happens but the command still is able to complete
    """
    def __init__(self, time, msg, obj, exception=None, formatedException=None, **kw):
        AppSuccess.__init__(self, time, msg, obj, **kw)
        self.exception = exception
        self.formatedException = formatedException
    
class AppError(AppWarning):
    """Object used to describe an error encountered during the execution of a command

    Errors should be issued if something happens and the command not able to complete
    """
    def __init__(self, time, msg, obj, exception=None, formatedException=None, **kw):
        AppWarning.__init__(self, time, msg, obj, exception, formatedException, **kw)
    
class AppException(AppWarning):
    """Object used to describe an unexpected exception encountered during the execution of a command
    """
    def __init__(self, time, msg, obj, exception, formatedException, **kw):
        AppWarning.__init__(self, time, msg, obj, exception, formatedException, **kw)
    
class AppExecutionReport:
    """Object used to store information about a command execution
"""
    def __init__(self, cmd):
        self.cmd = cmd
        self.app = cmd.app()
        self.initialize()
        
    def initialize(self):
        self.startTime = time()
        self.endTime = None
        self.timeStamps = [] # list of time:msg
        self.reports = [] # this list will hold a instances of MGLsuccess, MGLwarning, MGLErrors, and
                          # MGLexecption objects
        self._requestUserConfirmation = False
        self.numberOf = {
            'successes' : 0,
            'warnings' : 0,
            'errors' : 0,
            'exceptions': 0
            }

    def getSuccesses(self):
	return [x for x in self.reports if x.__class__ == AppSuccess]

    def getWarnings(self):
	return [x for x in self.reports if x.__class__ == AppWarning]

    def getErrors(self):
	return [x for x in self.reports if x.__class__ == AppError]

    def getExceptions(self):
	return [x for x in self.reports if x.__class__ == AppException]

    def setRequestUserConfirmation(self, boolean):
        self._requestUserConfirmation = boolean

    def addSuccess(self, msg, obj=None):
        self.reports.append(AppSuccess(time(), msg, obj))
        self.numberOf['successes'] += 1

    def throw(self, error):
        if self.app._stopOnError:
            raise error
        
    def addException(self, msg, error, obj=None):
        self.reports.append(AppException(time(), msg, obj, error, traceback.format_exc().split('\n')))
        self.setRequestUserConfirmation(True)
        self.numberOf['exceptions'] += 1
        self.throw(error)
        
    def addError(self, msg, error, obj=None):
        if error:
            formatedError = traceback.format_exc().split('\n')
        else:
            formatedError = None
        self.reports.append(AppError(time(), msg, obj, error, formatedError))
        self.numberOf['errors'] += 1
        if error: self.throw(error)
        
    def addWarnings(self, msg, error=None, obj=None):
        if error:
            formatedError = traceback.format_exc().split('\n')
        else:
            formatedError = None
        self.reports.append(AppWarning(time(), msg, obj, error, formatedError))
        self.numberOf['warnings'] += 1
        if error: self.throw(error)

    def timeStamp(self, message):
        self.timeStamps.append( (time(), message) )

    def finalize(self):
        self.endTime = time()
        event = ExecutionReportEvent(report = self)
        self.app.eventHandler.dispatchEvent(event)

    def printReport(self):
        print 'Execution Report for %s (%.2fs)'%(self.cmd.name, self.endTime - self.startTime), self.numberOf


class ExecutionReportEvent(Event):
    # event is created after a command is executed
    # event.report is an instance of ExecutionReportEvent
    pass

class BeforeDoitEvent(Event):
    # event is created before a command is executed
    # event.app is an istance of the application
    pass

class AfterDoitEvent(Event):
    # event is created after a command is executed
    # event.app is an instance of the application
    pass
    
class CmdLoader:
    """
    This object is used to create a gui and provide call backs
    for a cmd to load the actual command the first time the command is called
    """
    def isLoader(self):
        return True

    def __init__(self, app, package, module, cmdName, gui, _callable=True):
        """
        create an object to load a command when hte command is called

        app - instance of an app (i.e. sub class of ViewerFramework
        package - (string) nme of Python package (i.e. 'pmv')
        module - (string) name of a python file in package (i.e. 'selectionCommands.py')
        cmdName - (string) name of the command (i.e. 'select')
             NOTE: cmdName has to be a key in package.module.commandClassFromName{}
        gui - instance of CommandGUI     
        """
        self.package = package
        self.module = module
        self.cmdName = cmdName
        self.name = cmdName # used by color undo for instance
        self.app = weakref.ref(app)
        self.gui = gui
        self.cmdIsLoaded = False
        self.callable = _callable
        self.cmd = None
        
    def loadCommand(self):
        # actually load the command and add it to the application
        if self.cmdIsLoaded: return self.cmd
        importName = self.package + '.' + self.module
        # import the module
        mod = __import__(importName, globals(), locals(), [self.module])
        cmd = mod.commandClassFromName[self.cmdName][0]()
        cmd.GUI = self.gui
        
        # add the command to the application
        del self.app().commands[self.cmdName]
        self.app().addCommand( cmd, self.cmdName)
        self.cmdIsLoaded = True
        self.cmd = cmd
        self.cmd.proxy = weakref.ref(self)
        return cmd

    def safeLoadCommand(self):
        # load a command with error handling
        # Returns a cmd (instance if Command) or an Error instance
        # if loading failed
        if self.app().trapExceptions:
            if not self.cmdIsLoaded:
                errorMsg = 'Failed to load command %s'%self.name
                return self.app().GUI.safeCall(self.loadCommand, errorMsg)
        else:
            self.loadCommand()
        
    ## most Commands have no arg for guiCallback
    ## but icoms for instance have one args
    def guiCallback(self, *args, **kw):
        # returns an Error instance if somethign went wrong or the result
        # returned by calling guiCallback

        # safely load the cmd if needed
        if self.app().trapExceptions:
            cmd = self.safeLoadCommand()
            if isinstance(cmd, Error): return cmd

            # safely call cmd.guiCallback()
            errorMsg = 'Error while running command %s'%self.name
            self.app().GUI.safeCall(self.cmd.guiCallback, errorMsg, *args, **kw) 
        else:
            cmd = self.loadCommand()
            self.cmd.guiCallback(*args, **kw)
            
    def guiCall(self, *args, **kw):
        # this method is to be used to run this command from the GUI
        if not self.cmdIsLoaded: self.loadCommand()
        result = self.cmd.guiCall( *args, **kw)
        return result

    def __call__(self, *args, **kw):
        # This method is NOT SAFE

        if not self.cmdIsLoaded: self.loadCommand()

        #print '  CALLING __call__ for', self.cmd, args, kw
        result = self.cmd( *args, **kw)
        return result

    ## def getLastUsedValues(self, formName='default', **kw):
    ##     # safely load the cmd if needed
    ##     cmd = self.safeLoadCommand()
    ##     if isinstance(cmd, Error): return cmd

    ##     # safely callgetLastUsedValues
    ##     # print 'CALLING getLastUsedValues for', self.cmd
    ##     return self.app().GUI.safeCall(self.cmd.getLastUsedValues,
    ##                                 *(formName,), **kw)

    def onAddObjectToViewer(self, obj):
        # safely load the cmd is needed
        cmd = self.safeLoadCommand()
        if isinstance(cmd, Error): return cmd

        # safely cmd.onAddObjectToViewer
        return self.app().GUI.safeCall( self.cmd.onAddObjectToViewer, *(obj,) )


class AppCommand:
    """
    Base class for commands to be added to an AppFramework-based application

    commands are stored in app.commands[cmdName]
    """
    def isLoader(self):
        return False
    
    def __init__(self):
        self.app = None # will be a weakref to the application
        self.busyCursor = 'watch' # FIXME TK dependency

        self.lastUsedValues = {} # key is formName, values is dict of defaults
        self.lastUsedValues['default'] = self.getValNamedArgs()

        self.managedGeometries = [] # list of geometries 'owned' by this
                                    # command, used in updateGeoms
        self.initializedFor = {}  # this is used to check if command's initializeFor(object) has been called.
        self.createEvents = True # set to False to avoid command from issuing events
        self.iterateOverFirstArgument = False

    def expandArg0(self, obj):
        return obj

    def getValNamedArgs(self):
        """
         """
        from inspect import getargspec
        #args = getargspec(self.__call__.im_func)
        args = getargspec(self.checkArguments.im_func)
        allNames = args[0][1:]
        defaultValues = args[3]
        if defaultValues is None:
            defaultValues = []
        nbNamesArgs = len(defaultValues)
        posArgsNames = args[0][1:-nbNamesArgs]
        d = {}
        for name, val in zip(args[0][-nbNamesArgs:], defaultValues):
            d[name] = val
        return d

        
    def strArg(self, arg):
        before = ""
        if isinstance(arg, (list, tuple)):
            seenBefore = {} # used to save only once each before line needed
            if isinstance(arg, list):
                argstr = "["
                endstr = "], "
            else:
                argstr = "("
                endstr = ",), "
            for a in arg:
                astr, before = self._strArg(a)
                argstr += astr
                if before is not None:
                    seenBefore[before] = True
            # if condition is needed to fix bug #734 "incomplete log string"
            if len(argstr) > 1:
                argstr = argstr[:-2]+endstr
            else:
                argstr = argstr+endstr
            for s in seenBefore.keys():
                before += s+'\n'
            return argstr, before
        elif isinstance(arg, dict):
            seenBefore = {} # used to save only once each before line needed
            # d = {'k1':5, 'k2':6, 'k3':7, 'k8':14}
            argstr = "{"
            endstr = "}, "
            if len(arg)==0:
                #special handling for empty dictionary
                return "{}, ", before
            for key, value in arg.items():
                astr, before = self._strArg(key)
                if before is not None:
                    seenBefore[before] = True
                argstr += astr[:-2] + ':'
                astr, before = self._strArg(value)
                if before is not None:
                    seenBefore[before] = True
                argstr += astr[:-2] + ','
            argstr = argstr[:-1]+endstr
            return argstr, before
        else: # not a sequence
            return self._strArg(arg)
        
        
    def _strArg(self, arg):
        """
        Method to turn a command argument into a string,
FIXME describe what types of arguments are handled
        """
        from mglutil.util.misc import isInstance
        import inspect
        if inspect.isclass(arg):
            before = 'from %s import %s'%(arg.__module__, arg.__name__)
            return arg.__name__+', ', before

        elif isInstance(arg) is True:
            if isinstance(arg, Command):
                return 'self.'+arg.name+', ', None
            elif isinstance(arg, Geom):
                return "'"+arg.fullName+"', ", None
            elif isinstance(arg, ColorMap):
                return "'"+arg.name+"', ", None
            elif hasattr(arg, 'returnStringRepr'):
                # the returnStringRepr method has to be implemented by
                # the instance class and needs to return
                # the before string which can be None but usually is the
                # from module import class string
                # and the argst which is also a string allowing the
                # instanciation of the object.
                before, argst = arg.returnStringRepr()
                return argst+', ', before

            try:
                import cPickle
                pickle = cPickle
                before = 'import cPickle; pickle = cPickle'
            except:
                import pickle
                before = 'import pickle'
                self.app().log( before )
            sp = pickle.dumps(arg)
            # Add a \ so when the string is written in a file the \ or \n
            # are not interpreted.
            pl1 =  sp.replace('\\', '\\\\')
            picklelog = pl1.replace('\n', '\\n')
            return 'pickle.loads("' + picklelog + '"), ', before

        elif isinstance(arg, str):
            arg1 = arg.replace('\n', '\\n')
            if arg.find("'") != -1:
                return '"'+ arg1 + '",', None
            else:
                return "'" + arg1 + "', ", None

        elif isinstance(arg, numpy.ndarray):
            before = 'from numpy import array\n'
            #arg = arg.tolist()
            return repr(arg)+ ', ', before

        else:
            return str(arg) + ', ', None


    def buildLogArgList(self, args, kw):
        """build and return the log string representing the arguments
        a list of python statements called before is also built. This list
        has to be exec'ed to make sure the log can be played back"""
        if self.app() is None: return
        argString = ''
##         args, kw = self.getLogArgs( args, kw )
        before = []
        for arg in args:
            s, bef = self.strArg(arg)
            argString = argString + s
            if bef is not None: before.append(bef)
        for name,value in kw.items():
            s, bef = self.strArg(value)
            argString = argString + '%s=%s'%(name, s)
            if bef is not None: before.append(bef)
        return '('+argString[:-2]+')', before # remove last ", "
    

    def logString(self, *args, **kw):
        """build and return the log string"""

        argString, before = self.buildLogArgList(args, kw)
        log = ''
        for l in before: log = log + l + '\n'
        log = log + 'self.' + self.name + argString
        return log

    def getName(self, args, kw):
        return self.name

    def getUndoName(self, cmds):
        return self.getName(self._args, self._kw)
    
    def guiCall(self, *args, **kw):
        # this method is to be used to run this command from the GUI

        # returns an Error instance if something went wrong or the result
        # returned by calling guiCallback

        if self.app().gui is not None:
            # make cursor busy if there is a gui
            from PySide import QtGui, QtCore
            QtGui.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
            if self.app().trapExceptions:
                try:
                # safely load the cmd if needed
                    cmd = self.safeLoadCommand()
                    if isinstance(cmd, Error): return cmd

                    # safely call cmd.guiCallback()
                    errorMsg = 'Error while running command %s'%self.name
                    self.app().GUI.safeCall(self.cmd.guiCallback, errorMsg, *args, **kw) 
                except Exception as e:
                    raise e
                    print("Error {}".format(e.args[0]))
                finally:
                    QtGui.QApplication.restoreOverrideCursor()
            else:
                self(*args, **kw)
                QtGui.QApplication.restoreOverrideCursor()
        else:
            self(*args, **kw)

    def __call__(self, *args, **kw):
        t0 = time()
        try:
            app = self.app()
            app._cmdLevel += 1

            #e = kw.pop('topCommand', True)
            #assert e in [True, False, 0, 1], \
            #"topCommand can only be True or False, got %g"%e
            if app._cmdLevel == 1:
                self._topCommand = True
                app._topCommand = self # # remember who the top command is used for undo redo
                execRep = app._executionReport = AppExecutionReport(self)
                #execRep.initialize() 
            else:
                self._topCommand = False
                execRep = app._executionReport

            #print 'CALLING', self.name, self._topCommand, app._cmdLevel
            e = self._createEvents = kw.pop('createEvents', True)
            e = self._failForTesting = kw.pop('failForTesting', False)
            e = self._setupUndo = kw.pop('setupUndo', True)
            e = self._callDoitWrappers = kw.pop('callDoitWrappers', True)
            e = self._redraw = kw.pop('redraw', False)
            #import pdb; pdb.set_trace()
            if len(args):
                arg0 = self.expandArg0(args[0])
                args = (arg0,)+args[1:]

            args, kw = self.checkArguments(*args, **kw)
            self._args = args
            self._kw = kw
            execRep.timeStamp('%s.check arguments'%self.name) 

            if self._callDoitWrappers:
                self.beforeDoit()
                if self._topCommand:
                    execRep.timeStamp('%s.before execution callback'%self.name) 

            if self._topCommand:
                self.app().clearTopUndoCmds()

            if self.iterateOverFirstArgument:
                results = [] # list of results for each object in the list args[0]
                args0 = self._args[0][:]
                args1 = self._args[1:]
                for obj in args0:
                    if isinstance(obj, str):
                        objname = obj
                    elif hasattr(obj, "__repr__"):
                        objname = obj.__repr__()
                    else:
                        objname = str(obj)
                    try:
                        nbUndoCmds = len(app._topUndoCmds)
                        result = self.doit(*((obj,)+args1), **self._kw)
                        if result is not None:
                            results.append(result)
                        app._executionReport.addSuccess(
                            'cmd %s successful for %s'%(
                                self.getName(args, kw), objname),obj=obj)
                    except:
                        if len(app._topUndoCmds)>nbUndoCmds:
                            app._topUndoCmds = app._topUndoCmds[:nbUndoCmds]
                        if app.trapExceptions is False:
                            exc_info = sys.exc_info()
                            raise exc_info[1], None, exc_info[2]
                        else:
                            msg = 'Error while executing command %s for object %s'%(
                                self.getName(args, kw), objname)
                            app.errorMsg(sys.exc_info(), msg, obj=obj)
            else:
                try:
                    nbUndoCmds = len(app._topUndoCmds)
                    results = self.doit(*self._args, **self._kw)
                    app._executionReport.addSuccess(
                        'cmd %s successful for %s'%(
                            self.getName(args, kw), ''),obj=None)
                except:
                    if len(app._topUndoCmds)>nbUndoCmds:
                        app._topUndoCmds = app._topUndoCmds[:nbUndoCmds]
                    if app.trapExceptions is False:
                        exc_info = sys.exc_info()
                        raise exc_info[1], None, exc_info[2]
                    else:
                        msg = 'Error while executing command %s'%(
                            self.getName(args, kw))
                        app.errorMsg(sys.exc_info(), msg, obj=args)
                
            if self._topCommand:
                execRep.timeStamp('%s.execute command'%self.name) 

            if self._callDoitWrappers:
                self.afterDoit(results)
                if self._topCommand:
                    execRep.timeStamp('%s.after execution callback'%self.name) 

            if self._topCommand:
                execRep.finalize() 

            #print 'DONE', self.name, self._topCommand, app._cmdLevel
            #try:
            #    if app._cmdLevel==1:
            #        del app._cmdLevel
            #    else:
            #        app._cmdLevel -= 1
            #except:
            #    import pdb
            #    pdb.set_trace()
                
            #print 'Time for command %s: %f.4(s)'%(self.getName(args, kw), time()-t0)
            return results

        except:
            exc_info = sys.exc_info()
            if not app.trapExceptions:
                raise exc_info[1], None, exc_info[2]
            else:
                # FIXME define message
                app._executionReport.addException(exc_info[1].message, exc_info[1])
                if self._topCommand:
                    app._executionReport.finalize() 
        finally:
            app._cmdLevel -= 1

    def nameCmd(self, *args, **kw):
        return self.name
            
    def beforeDoit(self):
        """called before specialized doit method is called"""
        topCommand = self._topCommand

        event = BeforeDoitEvent(object=self.app())
        self.app().eventHandler.dispatchEvent(event)
        
        # Update self.lastUsedValues
        defaultValues = self.lastUsedValues['default']
        for key, value in self._kw.items():
            if defaultValues.has_key(key):
                defaultValues[key] = value

    
    def afterDoit(self, result):
        """called after specialized doit method is called"""
        # calls the cleanup methods after the doit has been called.

        topCommand = self._topCommand

        execRep = self.app()._executionReport
        if execRep.numberOf['errors'] + execRep.numberOf['exceptions'] > 0:
            self.cleanupUndoArgs()

        userpref = self.app().userpref
        if topCommand and len(self.app()._topUndoCmds) and \
               userpref.has_key('Number of Undo') and \
               userpref['Number of Undo']['value']>0 and self._setupUndo:

            # the last issued comamnd should be the first undo
            #self.app()._topUndoCmds.reverse()

            topCom = self.app()._topCommand
            
            if topCom == self.app().undo: # performing Undo. Negation command go into redo
                name = self.getUndoName(self.app()._topUndoCmds)
                self.app().redo.addUndoCall( *(self.app()._topUndoCmds, name) )
                
            else:  # performing cmd or Redo. Negation command go into Undo
                name = self.getUndoName(self.app()._topUndoCmds)
                self.app().undo.addUndoCall( *(self.app()._topUndoCmds, name) )

        if topCommand:
            kw = self._kw
            kw['createEvents'] = self._createEvents
            kw['setupUndo'] = self._setupUndo
            kw['callDoitWrappers'] = self._callDoitWrappers
            kw['redraw'] = self._redraw
            self.app().addCmdToHistory( self, self._args, kw)

        del self._args
        del self._kw
        
        ### The following code should go to an AppGUI method. The method
        ### will be registered with AfterDoitEvent()
        event = AfterDoitEvent(object=self.app())
        self.app().eventHandler.dispatchEvent(event)

    ##
    ## to be subclassed
    ##
    def checkArguments(*args, **kw):
        return args, kw

    ## def undoCmdBefore(self, result=None, *args, **kw):
    ##     """A chance to modify undo commands after the command was carried
    ##     out"""
    ##     return None

    ## def undoCmdAfter(self, *args, **kw):
    ##     """A chance to modify undo commands after the command was carried
    ##     out"""
    ##     return None

    def onAddCmdToApp(self, *args, **kw):
        return
    
    def doit(self):
        pass
    
    def cleanupUndoArgs(self):
        # call after doit was called if there was errors or exceptions
        # should clean the list of arguments of self.undoCmds
        pass
    
