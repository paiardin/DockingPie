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
#########################################################################
#
# $Header: /mnt/raid/services/cvs/AppFramework/App.py,v 1.21.4.1 2017/07/28 00:58:00 annao Exp $
#
# $Id: App.py,v 1.21.4.1 2017/07/28 00:58:00 annao Exp $
#

"""
define a base class AppFramework that supports lazy loading of commands
"""
import sys, os, thread, weakref, traceback
from time import time

from mglutil.events import Event, EventHandler
from mglutil.preferences import UserPreference
#from mglutil.errors import MGLException, MGLError, Error, ErrorManager, MessageManager
from mglutil.errors import MGLException, MessageManager

from AppFramework.AppCommands import CmdLoader, AppCommand

###########################################################################
##
##   EVENTS
##
###########################################################################
class StartAddObjectEvent(Event):
    # event is created at the begging of the addition of an object
    #
    # event.name is the name of the object added
    # event.object is the object being added
    pass


class AddObjectCmdEvent(Event):
    # event is created each time we call a command to be carried out on
    #  an object that is added
    #
    # event.number number of command currently called
    # event.total total number of commands to be called
    # event.cmdName name of the command currently called
    pass


class EndAddObjectEvent(Event):
    # event is created after the object has been added
    # 
    # event.object is the object that has been added
    pass


class BeforeDeleteObjectEvent(Event):
    # event is created before the object is deleted
    # 
    # event.object is the object that is being delete
    # event.objectType is the type of the object that is being delete
    pass


class AfterDeleteObjectEvent(Event):
    # event is created after the object is deleted
    # 
    # event.object is the object that is being deleted
    # event.objectType is the type of the object that is being delete
    pass

class DeleteObjectEvent(Event):
    # event.object is the object that is being deleted
    # event.objectType is the type of the object that is being deleted
    pass

class AddCmdEvent(Event):
    # event is created after command has been loaded
    # 
    # event.command is the command
    pass



class AddGeometryEvent(Event):
    # event.objects is a list of geometries to be added.
    # event.parent is the geometries parent.
    # event should contain just one func in its eventListeners
    # list of functions. This function should add the geometry to the
    # Viewer
    def __init__(self, obj, parent=None, redo=False):
        Event.__init__(self)
        self.object = obj
        self.parent = parent
        self.redo = redo

class RemoveGeometryEvent(Event):
    def __init__(self, obj):
        Event.__init__(self)
        self.object = obj

class RedrawEvent(Event):
    # event to redraw the Viewer
    # event.object is the app (from app we can find out the Viewer)
    def __init__(self, obj):
        Event.__init__(self)
        self.object = obj


class AppFramework:
    """
    Application that support lazy loading of commands

    
    """
    def __del__(self):
        print 'destructor called'

    def exit(self):
        # the MessageCallback creates a reference to the app that needs
        # to be removed to ensure garbage collection of the app
        from mglutil.errors import deleteMessageCallback
        deleteMessageCallback(self.handleMessage)
        del self.objectTypes # release validators
    
    def __init__(self, name='MyApp', eventHandler=None, errorManager=None,
                 messageManager=None):
        """
        Construct an instance of a AppFramework object
        """
        self._stopOnError = False
        
        assert isinstance(name, str)
        if eventHandler:
            assert isinstance(eventHandler, EventHandler)
            self.eventHandler = eventHandler
        else:
            self.eventHandler = EventHandler()

        self.name = name
        self.embeded = False # will be set to True if the App is embeded in another program

        # the history of command [(cmd, args, kw) or ("MSG", type, string)]
	self.cmdHistory = [] # history of command [(cmd, args, kw)]
        self._executionReport = None
        
        self.trapExceptions = True # set to false to not catch unexpected exceptions
        
        # MS I think we use undo.cmdStack and redo.cmdStack now
        #self.undoCmdStack = [] # list of strings used to undo

        self.commands = {}    # dictionnary of command added to a Viewer
        self._cmdLevel = 0
        
        self.objectTypes = {} # typeName: {
                              #   'validator':validator,
                              #   'events': {eventName:eventClass}
                              #   }
        # dict of objType: List of command that can be carried out each time
        # an object is added to the application
        # every entry is a tuple (function, args_tuple, kw_dict)
        self.onAddObjectCmds = {} 
        self.objects = {}     # typeName: list of objects loaded in the App

        self.colorMaps = {}   # available colorMaps
        self.colorMapCt = 0   # used to make sure names are unique

        # lock needs to be acquired before object can be added
        self.objectsLock = thread.allocate_lock()

        # lock needs to be acquired before topcommands can be run
        self.commandsLock = thread.allocate_lock()

        # nexted commands counter
        self.commandNestingLevel = 0

        self.createUserPreferences()
       
        ## # register a callback for handling messages
        from mglutil.errors import addMessageCallback
        addMessageCallback(self.handleMessage)

        ## load not optional commands
        from notOptionalCommands import UndoCommand, RedoCommand
        self.addCommand(UndoCommand(), 'undo')
        self.addCommand(RedoCommand(), 'redo')
        
        ## private attributes
        self._colorMapCt = 0   # used to make sure names are unique
        self._topUndoCmds = [] # used to accumulate undo commands
                               # for sub commands of a top command
        self._currentlySourcedFiles = [] # used to avoid recusive sourcing

        ## end private attributes
        #if errorManager:
        #    assert isinstance(errorManager, ErrorManager)
        #    self.errorManager = errorManager
        #else:
        #    self.errorManager = ErrorManager(name)

        if messageManager:
            self.messageManager = messageManager
        else:
            self.messageManager = MessageManager(name)
        from mglutil.util.recentFiles import RecentFiles
        self.recentFiles = RecentFiles(
                    self, None)

    def pushUndoCmd(self, cmd, args, kw):
        assert isinstance(cmd, (AppCommand, CmdLoader))
        assert isinstance(args, tuple)
        assert isinstance(kw, dict)
        self._topUndoCmds.append( (cmd, args, kw) )

    def clearTopUndoCmds(self):
        self._topUndoCmds = []
    
    def lazyLoad(self, module, commands=None, package="PmvApp"):
        """
    \nSynopsis:\n
        None <-- lazyLoad(module, commands=None, package='PmvApp')
    \nRequired Arguments:\n
        module --- name of the module (eg:basicCommands)
    \nOptional Arguments:\n
        commands --- one command name or a list of command names to load,
                     None means load all commands
        package --- name of the package to which module belongs(eg:PmvApp,Vision)

    The purpose of lazy loading commands is to speed up the initialization
    of the application and save memory by avoiding importing and instanciating
    commands before they are used.

    This method will install a command loader for the specified commands.
    A command loader is an instance of a lightweight class that will load the
    actual command once the command is invoked.
    """

        ## first try to create a gui without importing command
        ## only possible if package.moduleProxy is found
        try:
            importName = package + '.' + module+'Proxy'
            mod = __import__(importName, globals(), locals(),
                             [module])
        except ImportError:
            raise RuntimeError, 'Proxy file not found for %s.%s'%(
                package, module)

        guiClasses = mod.getGUI(None)
        if commands is None: # no particular commands, we get all
            commands = guiClasses.keys()
        elif isinstance(commands, str):
            commands = [commands]

        # at this point commands is a list of command names
        # we loop over the commands
        for cmdName in commands:
            # if the commands is already loaded do nothing
            if self.commands.has_key(cmdName): # already loaded
                continue
            # REPORT LOADING
            #print 'lazy', package, cmdName
            #guiClassList = mod.commandGUIs[cmdName]
            #if not isinstance(guiClassList, list): continue
            
            guiClassList = guiClasses.get(cmdName, None)

            #create the command loader
            if mod.commandsInfo['icoms'].has_key(cmdName):
                lClass = ICmdLoader
            else:
                lClass =  CmdLoader
            cmdLoader = lClass(self, package, module, cmdName, None)

            self.commands[cmdName] = cmdLoader
            setattr(self, cmdName, cmdLoader)


    def addCommand(self, command, name):
        """
        Add a command to a app.
        arguments:
           command: Command instance
           name: string

        name is used to create an alias for the command in the viewer
        """
        #print "addCommand", name
        from AppCommands import AppCommand
        assert isinstance(command, AppCommand)

        # if the commands is laread loaded return
        # happens because of dependencies
        if name in self.commands.keys():
            return self.commands[name]

        # REPORT LOADING
        #print "addCommand", name

        ## set the command's .app attribute
        command.app = weakref.ref(self)

        ## name sure the command name has no spaces
        name = name.strip()
        name = name.replace(' ', '_')

        # add the command to the app's commands dict
        self.commands[name] = command
        command.name = name
        command.undoMenuString = name   # string used for menu entry for Undo
        command.undoMenuStringPrefix='' # prefix used to change menu entry for Undo
        # create an alias in the App
        # FIXME this should be replaced using __getattr__
        setattr(self, name, command)

        # call the command's callback when adding to an app
        command.onAddCmdToApp()

        event = AddCmdEvent(command=command)
        self.eventHandler.dispatchEvent(event)


    def addObjectType(self, typeName, typeValidator, events):
        # add a type of object for the application and define
        # a validator for this type and event
        assert callable(typeValidator)
        assert isinstance(typeName, str)
        assert not self.objectTypes.has_key(typeName)
        assert isinstance(events, dict)
        for key in ['beforeAddObject', 'afterAddObject',
                    'beforeDeleteObject', 'afterDeleteObject']:
            assert events.has_key(key)
            assert events[key] is None or issubclass(events[key], Event)

        self.objectTypes[typeName] = {
            'validator':typeValidator,
            'events':events
            }
        self.objects[typeName] = []
        self.onAddObjectCmds[typeName] = []
        
        
    def setOnAddObjectCmd(self, objType, cmds, argsList=None, kwList=None):
        """
        define the list of commands that have to be carried out when an
        object is added to the App.

        oldCmds <- setOnAddObjectCmd(objType, cmds, argsList=None, kwList=None)

        objType --- the type of objects for which we define the list
        cmds ---    a list of commands
        args ---    has to be a list of tuples of the same length as cmds
        kw ---      has to be a list of dict of the same length as cmds

        raises assertErrors if somethign goes wrong

        returns oldList of commands
        """

        assert self.objectTypes.has_key(objType), "bad objecttype %s"%objType

        if argsList is None:
            argsList = [()]*len(cmds)
        else:
            assert len(argsList)==len(cmds)
            for value in argsList: assert isinstance(value, tuple)
            
        if kwList is None:
            kwList = [{}]*len(cmds)
        else:
            assert len(kwList)==len(cmds)
            for value in kwList: assert isinstance(value, dict)
            
        oldList = self.onAddObjectCmds[objType][:]
        #self.onAddObjectCmds[objType] = []

        for cmd, args, kw, in zip(cmds, argsList, kwList):
            if isinstance(cmd, CmdLoader):
                cmd = cmd.loadCommand()
            
            assert callable(cmd), "%s is not callable"%str(cmd)
            #assert cmd.flag & Command.objArgOnly

            # set parameter to call this command when the molecule is loaded
            #kw['undo'] = 0
            self.onAddObjectCmds[objType].append( (cmd, args, kw) )

        return oldList

    ##
    ## method for adding an object to the apllication
    ##
    def addObject(self, name, obj, objType, geomContainer=None):
        """
        Add an object to a Viewer
        None <-  name, obj, objType, geomContainer=None)

        name is the name of he object
        obj is the actual object
        objType is the type of the object. It has to be in self.objectTypes
            and the validator for this object has to return True hen called
            with obj
        """

        # make sure we add an object of the proper type
        validator = self.objectTypes[objType]['validator']
        assert validator(obj), "%s is not of type %s"%(name, objType)

        events = self.objectTypes[objType]['events']
        # create StartAddObjectEvent
        event = events['beforeAddObject'](name=name, object=obj)
        self.eventHandler.dispatchEvent(event)
        
        #print 'acquiring addObject lock'
        self.objectsLock.acquire()
        self.objects[objType].append(obj)
        self.objectsLock.release()
        #print 'releasing addObject lock'

        if geomContainer is None:
            geomContainer = GeomContainer( self )
        else:
            assert isinstance(geomContainer, GeomContainer)

        obj.geomContainer = geomContainer

        # create add object event
        if events['afterAddObject']:
            event = events['afterAddObject'](object=obj)
            self.eventHandler.dispatchEvent(event)

        ## # call commands registered to be called when obj is added
        ## #t0 = time()
        #l = len(self.onAddObjectCmds[objType])
        #eventClass = events['addObjectCmd']
        #import pdb
        #pdb.set_trace()

    def applyDefaultCommands(self, obj, objType):
        # apply the commands set using setOnAddObjectCmd() for this type
        # of object. Object types are registers using addObjectType()
        if self.gui:
            old = self.gui().viewer.suspendRedraw
            self.gui().viewer.suspendRedraw = True
        for com in self.onAddObjectCmds[objType]:
            com[2]['redraw']=0
            com[0]( *((obj,)+com[1]), **com[2] )            
        if self.gui:
            self.gui().viewer.suspendRedraw = old

    def removeObject(self, obj, objType):
        """Remove an object from a Viewer"""
        #1 Delete the obj from the list of objects.
        validator = self.objectTypes[objType]['validator']
        events = self.objectTypes[objType]['events']
        
        event = events['beforeDeleteObject'](object=obj)
        self.eventHandler.dispatchEvent(event)
        del(self.objects[objType][self.objects[objType].index(obj)])

        # clean up the managedGeometries list
        if obj.geomContainer:
            for cmd in self.commands.values():
                if isinstance(cmd, CmdLoader): continue
                if len(cmd.managedGeometries)==0: continue
                geomList = []
                for g in cmd.managedGeometries:
                    if hasattr(g, 'mol') and g.mol==obj:
                        continue
                    geomList.append(g)
                cmd.managedGeometries = geomList

            # remove everything created in the geomContainer associated to the
            # mol we want to destroy,
            obj.geomContainer.delete()

        # create remove object event
        event = events['afterDeleteObject'](object=obj)
        self.eventHandler.dispatchEvent(event)


    def handleMessage(self, type, message):
        """
        print the message
        """
        msgstr = 'Message: (%s) %s'%(type, message)
        self.cmdHistory.append( ('MSG', type, message) )
        print msgstr


    def getHistory(self, commands=True, messages=True):
        """
        generate strings for all commands and/or messages in history list
        """
        lines = []
        i = 0
        for cmd, args, kw in self.cmdHistory:
            if cmd=='MSG' and Messages is True:
                lines.append('Message: (%s) %s\n'%(args, kw))
            elif commands:
                try:
                    log = cmd.logString( *args, **kw)+'\n'
                except Exceptions, e:
                    lines.append('Failed to create log for (%d) %s %s %s\n'%(
                        i, str(cmd), str(args), str(kw)))
                    lines.expand(traceback.format_exc().split('\n'))
            i += 1       
        return lines


    def addCmdToHistory(self, cmd, args, kw):
        """
        append a command to the history of commands
        """
        #print "ADDING Command to history", cmd.name
        self.cmdHistory.append( (cmd, args, kw))
        maxLen = self.userpref['Command History Depth']['value']
        lenCmds = len(self.cmdHistory)
        if maxLen>0 and lenCmds > maxLen:
            #print "maxLen", maxLen, lenCmds
            self.cmdHistory = self.cmdHistory[-maxLen:]


    def source(self, filename, globalNames=1, globDict={}):
        # if globalNames==1 the objects defined in the file will be
        # visible in the main interpreter

        if not os.path.exists(filename):
            raise ValueError, "file %s does not exist"%filename

        if filename in self._currentlySourcedFiles:
            from mglutil.errors import MGLMessage
            message = 'WARNING: %s is already being sourced\n'%filename
            message += '         skipping to avoid endless loop'
            MGLMessage('usererror', message)
            return
        
        self._currentlySourcedFiles.append(filename)
        # make self know while executing filename
        globDict['self'] = self
        globDict['pmv'] = self
       
        try:
            # if we pass loaclDict Vision macros defined explicitely in 
            # alog file do not build .. strange scope problem
            execfile( filename, globDict)#, localDict)
        finally:
            self._currentlySourcedFiles.pop()# = self._currentlySourcedFiles[:-1]
        
        if globalNames:
            # all global objects created in filename are now in glob
            # we remove self and update the interpreter's main dict
            del globDict['self']
            sys.modules['__main__'].__dict__.update(globDict)


    def createUserPreferences(self):
        self.userpref = UserPreference()

        ## start up folder. By default it is current folder
        ##
        def validateStartupDir(value):
            # return True if the folder exists
            value = os.path.expanduser(value)
            return os.path.exists(value) and os.path.isdir(value)

        def setStartupDir(folder, old, new):
            folder = os.path.expanduser(new)
            if not os.path.isdir(folder):
                msg = "try to set startup folder: %s is not a folder"%new
                raise ValueError, msg
            else:
                os.chdir(new)
            
        self.userpref.add(
            'Startup Folder', os.getcwd(), validateFunc=validateStartupDir,
            callbackFunc = [setStartupDir],
            doc="Startup Folder. Defines the working directory for the app")

        ## cmd history depth
        ##
        def validateCmdHistDepth(value):
            # needs to be positive integer
            return isinstance(value, int) and value >-1

        self.userpref.add(
            'Command History Depth', 20, validateFunc=validateCmdHistDepth,
            doc="Set Command History Depth, i.e. number of commands kept in the command history list")

        ## cmd history depth
        ##
        self.userpref.add(
            'Sharp Color Boundaries for polygonal surfaces', 'sharp',
            ('sharp', 'blur'),
            doc="""Color boundaries for polygonal surface [sharp or blur]
(will not modify already displayed surfaces when set,
but will be used a new surfaces are created)""", category="Graphics")
        
    ## def handleError(self, exc_info, message, obj=None, errorHandling='raise'):
    ##     # gets call by commands to decide what to do when an exception
    ##     # was raised during the execution of the command
    
    ##     # handle errors
    ##     #  'raise': raise the exception
    ##     #  'reportAndRaise': add exception to report and raise
    ##     #  'report': add exception to report but do not raise exception

    ##     if errorHandling=='raise':
    ##         # re-throw exception so that traceback and pdb.pm takes you back
    ##         # to the original exception
    ##         raise exc_info[1], None, exc_info[2]
    ##     elif errorHandling=='reportAndRaise':
    ##         self._executionReport.addError(message, exc_info[1], obj=obj)
    ##         raise exc_info[1], None, exc_info[2]
    ##     elif errorHandling=='report':
    ##         self._executionReport.addError(message, exc_info[1], obj=obj)

    def errorMsg(self, exc_info, message, obj=None):
        """None <- warningMsg(msg)"""
        self._executionReport.addError(message, exc_info[1], obj=obj)

    def warningMsg(self, message):
        """None <- warningMsg(msg)"""
        self._executionReport.addWarnings(message, None)

    def successMsg(self, message):
        """None <- warningMsg(msg)"""
        self._executionReport.addError(message, None)

    #def error(self, exception, msg, errorManager=None, raiseException=False,
    #          **kw):
    #    return MGLError(exception, msg, errorManager=self.errorManager,
    #             raiseException=raiseException, **kw)

    #def getErrors(self):
    #    return self.errorManager.errorStack

    #def hasErrors(self):
    #    return len(self.errorManager.errorStack)>0

    #def clearErrors(self):
    #    self.errorManager.clear()

    #def displayErrors(self):
    #    self.errorManager.displayErrors()

    def message(self, type, text):
        self.messageManager.message(type, text)
    
    def addMessageCallback(self, func):
        self.messageManager.addCallback(func)

    def deleteMessageCallback(self, func):
        self.messageManager.deleteCallback(func)

    def addColorMap(self, colorMap):
        from DejaVu.colorMap import ColorMap
        assert isinstance(colorMap, ColorMap)
        if self.colorMaps.has_key('colorMap.name'):
            self.warningMsg('invalid attemp to replace an existing colormap')
        else:
            self.colorMaps[colorMap.name] = colorMap
            
    def loadColormap(self, filename=None):
        # load the rgb256_map  by default
        from DejaVu.colorMap import ColorMap
        if not filename:
            from mglutil.util.packageFilePath import findFilePath
            cmapPath = findFilePath('ColorMaps', 'AppFramework')
            filename = os.path.join(cmapPath, "rgb256_map.py")
        #colormap can be built from a filename
        name = os.path.splitext(os.path.basename(filename))[0]
        l = {}
        g = {}
        execfile(filename, g, l)
        newColorMap = None
        for name, obj in l.items():
            if isinstance(obj, ColorMap):
                newColorMap = obj
                break
            
        #newColorMap = ColorMap(name, filename=filename)
        if newColorMap:
            #name should be uniq already
            for k in self.colorMaps.keys():
                if k==name:
                    newColorMap.name = name + '_' + str(self.colorMapCt)
                    self.colorMapCt = self.colorMapCt + 1
            self.addColorMap(newColorMap)
            
        return newColorMap



class GeomContainer:
    """Class to hold geometries.
    This class provides 'geoms' dictionary --
    {name: DejaVu geometry object}. 
    Geometries can be added using addGeometry() method.
    """

    def __init__(self, app):
        """constructor of the geometry container"""

        ## Dictionary of geometries(sticks, balls, CPK spheres, etc)
        ## used to display atoms of a molecule
        self.app = app
        self.geoms = {}
        self.masterGeom = None

        ## Dictionary linking geom names to cmds which update texture coords
        ## for the current set of coordinates
        self.texCoordsLookup = {}
        self.updateTexCoords = {}


    def delete(self):
        """Method to remove self.geoms['master'] and
        self.geoms['selectionSpheres'] from the viewer when deleted"""

        # switch the object and descendant to protected=False
        for c in self.geoms['master'].AllObjects():
            c.protected = False

        event = RemoveGeometryEvent(self.geoms['master'])
        self.app.eventHandler.dispatchEvent(event)


    def addGeom(self, geom, parent=None, redo=False):
        """
        This method should be called to add a molecule-specific geometry.
        geom     -- DejaVu Geom instance
        parent   -- parent geometry, if not specified we use self.masterGeom
        """
        if parent is None:
            parent = self.masterGeom

        # we need to make sure the geometry name is unique in self.geoms
        # and in parent.children
        nameUsed=False
        geomName = geom.name
        for object in parent.children:
            if object.name==geomName:
                nameUsed=True
                break

        if nameUsed or self.geoms.has_key(geomName):
            newName = geomName+str(len(self.geoms))
            geom.name = newName
            import warnings
            warnings.warn("renaming geometry %s to %s"%(geomName, newName))#, stacklevel=14)
        self.geoms[geomName]=geom

        # add the geometry to the viewer.  At this point the name should be
        # unique in both the parent geoemtry and the geomContainer.geoms dict
        
        event = AddGeometryEvent(geom, parent=parent, redo=redo)
        self.app.eventHandler.dispatchEvent(event)

        if geom not in parent.children:
            parent.children.append(geom)
            geom.parent = parent
            geom.fullName = parent.fullName+'|'+geom.name
