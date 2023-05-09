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

############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner TSRI 201
#
#############################################################################
"""
Module implementing the commands that are present when instanciating
an AppFramework class or AppFramework derived class.
   - loadModuleCommand
   - UndoCommand
   - RedoCommand
   - BrowseCommandsCommand
"""

# $Header: /mnt/raid/services/cvs/AppFramework/notOptionalCommands.py,v 1.7.4.1 2017/07/28 00:58:01 annao Exp $
#
# $Id: notOptionalCommands.py,v 1.7.4.1 2017/07/28 00:58:01 annao Exp $
#

## FIXME these should become part of the AppFramework rather than commands
##

import os, sys
from string import join

from mglutil.util.packageFilePath import findFilePath, findModulesInPackage
from AppFramework.AppCommands import AppCommand
from mglutil.events import Event

class NewUndoEvent(Event):
    pass

class AfterUndoEvent(Event):
    pass

class AfterRedoEvent(Event):
    pass


commandslist=[]
cmd_docslist={}
def findAllAppPackages():
    """Returns a list of package names found in sys.path"""
    packages = {}
    for p in ['.']+sys.path:
        flagline = []
        if not os.path.exists(p) or not os.path.isdir(p):
            continue
        files = os.listdir(p)
        for f in files:
            pdir = os.path.join(p, f)
            if not os.path.isdir(pdir):
                continue
            if os.path.exists( os.path.join( pdir, '__init__.py')) :
            
                fptr =open("%s/__init__.py" %pdir)
                Lines = fptr.readlines()
                flagline =filter(lambda x:x.startswith("packageContainsVFCommands"),Lines)
                if not flagline ==[]:
                    if not packages.has_key(f):
                        packages[f] = pdir
    return packages




class UndoCommand(AppCommand):
    """pops undo string from the stack and executes it in the AppFrameworks
    scope
    \nPackage : AppFramework
    \nModule : notOptionalCommands.py
    \nClass : UndoCommand
    \nCommand : Undo
    \nSynopsis:\n
        None <- Undo()
    """

    def validateUserPref(self, value):
        try:
            val = int(value)
            if val >-1:
                return 1
            else:
                return 0
        except:
            return 0


    def onAddCmdToApp(self):
        doc = """Number of commands that can be undone"""
        self.app().userpref.add( 'Number of Undo', 100,
                              validateFunc=self.validateUserPref,
                              doc=doc)


    def addUndoCall(self, cmdList, name):
        #print self.name, "addUndoCall for:", name
        # FIXME handle user pref
        self.cmdStack.append( (cmdList, name) )
        maxLen = self.app().userpref['Number of Undo']['value']
        if maxLen>0 and len(self.cmdStack)>maxLen:
            forget = self.cmdStack[:-maxLen]
            self.cmdStack = self.cmdStack[-maxLen:]

            for cmdList, name in forget:
                for cmd, args, kw in cmdList:
                    if hasattr(cmd, "handleForgetUndo"):
                        cmd.handleForgetUndo(*args, **kw)
        #the gui part of the application should register the following
        # event listener that will update the label if necessary
        event = NewUndoEvent(objects=self.cmdStack, command=self)
        self.app().eventHandler.dispatchEvent(event)

    def getUndoName(self, cmds):
        name = ""
        for cmd in cmds:
            name = name + cmd[0].getName(cmd[1], cmd[2]) + '\n'
        return name
    
    def doit(self, **kw):
        """
        pop cmdList from stack and execute each cmd in cmdlList    
        """

        stack = self.cmdStack
        if stack:
            cmdList, name = stack.pop()
            ncmds = len(cmdList)
            self._cmdList = ([], name) # this list will gather undoCommands generated during the undo
            for i, (cmd, args, kw) in enumerate(cmdList):
                self.inUndo = ncmds-i-1
                if hasattr(cmd, 'name'):
                    name = cmd.name # this is a command
                else:
                    #a method or a function
                    if hasattr(cmd, "im_class"):
                        name = "%s.%s" % (cmd.im_class, cmd.__name__)
                    else:
                        name = cmd.__name__
                #msg = "Failed to run %s from %s"%(name, self.name)
                cmd( *args, **kw)
                #self.app().GUI.safeCall( cmd, msg, *args, **kw)
            self._cmdList = () # this list will gather undoCommands generated during the undo
            #self.inUndo = True
            #for cmd, args, kw in cmdList:
            #    cmd( *args, **kw)
            #self.inUndo = False
            self.inUndo = -1
        else:
            self.app().warningMsg('ERROR: Undo called for %s when undo stack is empty'%\
                               self.name)
        event = AfterUndoEvent(objects=self.cmdStack, command=self)
        self.app().eventHandler.dispatchEvent(event)

    def __init__(self):
        AppCommand.__init__(self)
        # cmdStack is a list of tuples providing 1-a list of commands to execute and 2 a name for this operation
        # the list of commands is in the following format [ (cmd, *args, **kw) ]
        self.cmdStack = []
        self.inUndo = -1 # will be 0 or a positive integer  while we are executing command(s) to undo last operation.
        self._cmdList = () # this tuple will contain a list that will collect negation of commands during a loop over commands
                           # corresponding to an Undo (or Redo in subclassed command)

    def checkArguments(self, **kw):
        """None<---NEWundo()
        """
        kw['topCommand'] = 0
        return (), kw


    def resetCmdStack(self):
        #remove all items from self.cmdStack
        if len(self.cmdStack):
            del(self.cmdStack)
            self.cmdStack = []
            event = AfterUndoEvent(objects=self.cmdStack, command=self)
            self.app().eventHandler.dispatchEvent(event)


class  RedoCommand(UndoCommand):
    """pops redo cmdList from the stack and executes it in the AppFrameworks
    scope
    \nPackage : AppFramework
    \nModule : notOptionalCommands.py
    \nClass : RedoCommand
    \nCommand : Undo
    \nSynopsis:\n
        None <- Undo()
    """
    pass


class BrowseCommandsCommand(AppCommand):
    """Command to load dynamically either modules or individual commands
    in the Application.
    \nPackage : AppFramework
    \nModule : notOptionalCommands.py
    \nClass : BrowseCommandsCommand
    \nCommand : browseCommands
    \nSynopsis:\n
        None <-- browseCommands(module, commands=None, package=None, **kw)
    \nRequired Arguements:\n
        module --- name of the module(eg:colorCommands)
    \nOptional Arguements:\n
        commnads --- one list of commands to load
        \npackage --- name of the package to which module belongs(eg:Pmv,Vision)
    """
    def __init__(self):
        AppCommand.__init__(self)
        self.allPack = {}
        self.packMod = {}
        self.allPackFlag = False
        self.txtGUI = ""


    def doit(self, module, commands=None, package=None, removable=False, gui=False):
#        if removable:
#            self.app().removableCommands.settings[module] = [commands, package]
#            self.app().removableCommands.saveAllSettings()
        # If the package is not specified the default is the first library

        #global commandslist,cmd_docslist
        #import pdb
        #pdb.set_trace()
        if package is None: package = self.app().libraries[0]

        importName = package + '.' + module
        try:
            # try to execute import Pmv.colorCommands
            mod = __import__(importName, globals(), locals(),
                            [module])
        except:
            if self.cmdForms.has_key('loadCmds') and \
                self.cmdForms['loadCmds'].f.winfo_toplevel().wm_state() == \
                                                                       'normal':
                   self.app().errorMsg(sys.exc_info(),
                                       "ERROR: Could not load module %s"%module,
                                       obj=module )
            elif self.app().loadModule.cmdForms.has_key('loadModule') and \
                self.app().loadModule.cmdForms['loadModule'].f.winfo_toplevel().wm_state() == \
                                                                       'normal':
                   self.app().errorMsg(sys.exc_info(),
                                       "ERROR: Could not load module %s"%module,
                                       obj=module)
            else:
                self.app().errorMsg(sys.exc_info(),
                                    "ERROR: Could not load module %s"%module,
                                    obj=module)
            #import traceback
            #traceback.print_exc()

        if commands is None:
            # no particular commmand is asked for, so we try
            # to run the initModule
            if hasattr(mod,"initModule"):
                mod.initModule(self.app(), gui=gui)
            elif hasattr(mod, 'commandList'):
                for d in mod.commandList:
                    cmd = d['cmd'].__class__()
                    self.app().addCommand( cmd, d['name'], None)
            elif hasattr(mod, 'commandClassFromName'):
                for name, values in mod.commandClassFromName.items():
                    cmd = values[0]()
                    self.app().addCommand( cmd, name, None)
            else :
                raise RuntimeError, "cannot load module %s, missing init"%importName

        else:  # a single com,mand or a list of commands was given
            if isinstance(commands, str):
                commands = [commands,]

            elif hasattr(mod, 'commandList'):
                for cmdName in commands:
                    found = False
                    for d in mod.commandList:
                        if d['name']==cmdName:
                            cmd = d['cmd'].__class__()
                            self.app().addCommand( cmd, d['name'], d['gui'])
                            found = True
                            break
                    if not Found:
                        raise RuntimeError, 'ERROR: cmd %s not found in %s'%(cmdName, importName)
                        
            elif hasattr(mod, 'commandClassFromName'):
                for cmdName in commands:
                    values = mod.commandClassFromName.get(cmdName, None)
                    if values:
                        cmd = values[0]()
                        # FIXME gui are instances, that measn that 2 PMV would share
                        # these instances :(. Lazy loading fixes this since the GUI is
                        # created independently
                        gui = values[1]
                        self.app().addCommand( cmd, cmdName, gui)
                    else:
                        raise RuntimeError, 'ERROR: cmd %s not found in %s'%(cmdName, importName)
            else :
                raise RuntimeError, "cannot load module %s, missing init"%importName


    def checkArguments(self, module, commands=None, package=None, **kw):
        """None<---browseCommands(module, commands=None, package=None, **kw)
        \nmodule --- name of the module(eg:colorCommands)
        \ncommnads --- one list of commands to load
        \npackage --- name of the package to which module belongs(eg:Pmv,Vision)
        """
        kw['commands'] = commands
        kw['package'] = package
        return (module,), kw 
    
 
class loadModuleCommand(AppCommand):
    """Command to load dynamically modules to the App. import the file called name.py and execute the function initModule defined in that file Raises a ValueError exception if initModule is not defined
    \nPackage : AppFramework
    \nModule : notOptionalCommands.py
    \nClass : loadModuleCommand
    \nCommand : loadModule
    \nSynopsis:\n
        None<--loadModule(filename, package=None, **kw)
    \nRequired Arguements:\n
        filename --- name of the module
    \nOptional Arguments:\n   
        package --- name of the package to which filename belongs
    """

    active = 0
    
    def doit(self, filename, package):
        # This is NOT called because we call browseCommand()"
        if package is None:
            _package = filename
        else:
            _package = "%s.%s"%(package, filename)
        try:
            mod = __import__( _package, globals(), locals(), ['initModule'])
            if hasattr(mod, 'initModule') or not callable(mod.initModule):
                mod.initModule(self.app())
            else:
                self.app().errorMsg(sys.exc_info(), '%s:Module %s does not have initModule function'%(self.name, filename))    
        except:
            self.app().errorMsg(sys.exc_info(), '%s:Module %s could not be imported'%(self.name, _package))
        
            

    def checkArguments(self, filename, package=None, **kw):
        """None<---loadModule(filename, package=None, **kw)
        \nRequired Arguements:\n
            filename --- name of the module
        \nOptional Arguements:\n   
            package --- name of the package to which filename belongs
        """
        if package==None:
            package=self.app().libraries[0]
        if not kw.has_key('redraw'):
            kw['redraw'] = 0
        kw['package'] = package
        return (filename,), kw




    def loadModules(self, package, library=None):
        modNames = []
        doc = []
        self.filenames={}
        self.allPack={}
        self.allPack=findAllVFPackages()
        if package is None: return [], []
        if not self.filenames.has_key(package):
            pack=self.allPack[package]
            #finding modules in a package
            self.filenames[pack] =findModulesInPackage(pack,"^def initModule",fileNameFilters=['Command'])
        # dictionary of files keys=widget, values = filename
        for key, value in self.filenames[pack].items():
            pathPack = key.split(os.path.sep)
            if pathPack[-1] == package:
                newModName = map(lambda x: x[:-3], value)
                #for mname in newModName:
                   #if not modulename has Command in it delete from the
                   #modules list  
                   #if "Command" not in mname :
                       #ind = newModName.index(mname)
                       #del  newModName[ind]
                   #if "Command"  in mname :
                if hasattr(newModName,"__doc__"):
                    doc.append(newModName.__doc__)
                else:
                    doc.append(None)
                modNames = modNames + newModName
                
            else:
                pIndex = pathPack.index(package)
                prefix = join(pathPack[pIndex+1:], '.')
                newModName = map(lambda x: "%s.%s"%(prefix, x[:-3]), value)
                #for mname in newModName:
                   #if not modulename has Command in it delete from the
                   #modules list
                   #if "Command" not in mname :
                       #ind = newModName.index(mname)
                       #del  newModName[ind]
                if hasattr(newModName,"__doc__"):
                    doc.append(newModName.__doc__)
                else:
                    doc.append(None)      
                modNames = modNames + newModName
            modNames.sort()
            return modNames, doc     
 
