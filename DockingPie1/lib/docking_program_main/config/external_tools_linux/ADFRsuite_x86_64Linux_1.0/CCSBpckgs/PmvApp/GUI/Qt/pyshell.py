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

#------------------------------------------------------------------------------ 
# Copyright (c) 2007, Riverbank Computing Limited 
# All rights reserved. 
#  
# This software is provided without warranty under the terms of the GPL v2 
# license. 
#  
# Author: Riverbank Computing Limited 
# Description: <Enthought pyface package component> 
#------------------------------------------------------------------------------ 
 
 
# Standard library imports. 
import sys 
import code 
from ununicode import toascii

# Major package imports. 
from PySide import QtCore, QtGui 
 
# Enthought library imports. 
#from enthought.traits.api import Event, implements 
 
# Private Enthought library imports. 
#from enthought.util.clean_strings import python_name 
#from enthought.traits.ui.qt4.clipboard import PyMimeData 
 
# Local imports. 
#from enthought.pyface.i_python_shell import IPythonShell, MPythonShell 
#from enthought.pyface.key_pressed_event import KeyPressedEvent 
#from widget import Widget 
 
 
## class PythonShell(MPythonShell, Widget): 
##     """ The toolkit specific implementation of a PythonShell.  See the 
##     IPythonShell interface for the API documentation. 
##     """ 
 
## #    implements(IPythonShell) 
 
##     #### 'IPythonShell' interface ############################################# 
 
##     command_executed = Event 
 
##     key_pressed = Event(KeyPressedEvent) 
 
##     ########################################################################### 
##     # 'object' interface. 
##     ########################################################################### 
 
##     # FIXME v3: Either make this API consistent with other Widget sub-classes 
##     # or make it a sub-class of HasTraits. 
##     def __init__(self, parent, **traits): 
##         """ Creates a new pager. """ 
 
##         # Base class constructor. 
##         super(PythonShell, self).__init__(**traits) 
 
##         # Create the toolkit-specific control that represents the widget. 
##         self.control = self._create_control(parent) 
 
##         # Set up to be notified whenever a Python statement is executed: 
##         self.control.exec_callback = self._on_command_executed 
 
##     ########################################################################### 
##     # 'IPythonShell' interface. 
##     ########################################################################### 
 
##     def interpreter(self): 
##         return self.control.interpreter 
 
##     def execute_command(self, command, hidden=True): 
##         self.control.run(command, hidden) 
 
##     ########################################################################### 
##     # Protected 'IWidget' interface. 
##     ########################################################################### 
 
##     def _create_control(self, parent): 
##         # FIXME v3: Note that we don't (yet) support the zoom(?) of the wx 
##         # version. 
##         return PyShell(parent) 
 
 
class PyShell(QtGui.QTextEdit): 
    """ A simple GUI Python shell until we do something more sophisticated. """ 
 
    def __init__(self, parent=None): 
        """ Initialise the instance. """ 
 
        QtGui.QTextEdit.__init__(self, parent) 
 
        self.setAcceptDrops(True) 
        self.setAcceptRichText(False) 
        self.setWordWrapMode(QtGui.QTextOption.WrapAnywhere) 

        # make '__main__' variables known in '__console__'
        mod = __import__('__main__')
        local = {}
        for k,v in mod.__dict__.items():
            if k[:2]!='__':
                local[k] = v
        self.interpreter = code.InteractiveInterpreter(local)
 
        self.exec_callback = None 
 
        self._line = ""
        self._lines = [] 
        self._more = False 
        self._history = [] 
        self._pointer = 0 
        self._reading = False 
        self._point = 0 
 
        # Interpreter prompts. 
        try: 
            sys.ps1 
        except AttributeError: 
            sys.ps1 = ">>> " 
 
        try: 
            sys.ps2 
        except AttributeError: 
            sys.ps2 = "... " 
 
        # Interpreter banner. 
        self.write('Python %s on %s.\n' % (sys.version, sys.platform)) 
        self.write('Type "copyright", "credits" or "license" for more information.\n') 
        self.write(sys.ps1) 
 
    def flush(self): 
        """ Emulate a file object. """ 
 
        pass 
 
    def isatty(self): 
        """ Emulate a file object. """ 
 
        return 1 
 
    def readline(self): 
        """ Emulate a file object. """ 
 
        self._reading = True 
        self._clear_line() 
        self.moveCursor(QtGui.QTextCursor.EndOfLine) 
 
        while self._reading: 
            QtCore.QCoreApplication.processEvents() 
 
        if len(self._line):
            return str(self._line)  
 
        return '\n' 
     
    def write(self, text): 
        """ Emulate a file object. """ 
 
        self.insertPlainText(text) 
        self.ensureCursorVisible() 
 
    def writelines(self, text): 
        """ Emulate a file object. """ 
 
        map(self.write, text) 
 
    def run(self, command, hidden=False): 
        """ Run a (possibly partial) command without displaying anything. """ 
 
        self._lines.append(command) 
        source = '\n'.join(self._lines) 
        if not hidden: 
            self.write(source + '\n') 
 
        # Save the current std* and point them here. 
        old_stdin, old_stdout, old_stderr = sys.stdin, sys.stdout, sys.stderr 
        sys.stdin = sys.stdout = sys.stderr = self 
 
        self._more = self.interpreter.runsource(source) 
 
        # Restore std* unless the executed changed them. 
        if sys.stdin is old_stdin: 
            sys.stdin = old_stdin 
 
        if sys.stdout is old_stdout: 
            sys.stdout = old_stdout 
 
        if sys.stderr is old_stderr: 
            sys.stderr = old_stderr 
 
        if not self._more: 
            self._lines = [] 
 
            if self.exec_callback: 
                self.exec_callback() 
 
        if not hidden: 
            self._pointer = 0 
            self._history.append(command) 
 
            self._clear_line() 
 
            if self._more: 
                self.write(sys.ps2) 
            else: 
                self.write(sys.ps1) 
 
    def dragEnterEvent(self, e): 
        """ Handle a drag entering the widget. """ 
 
        if self._dragged_object(e) is not None: 
            # Make sure the users knows we will only do a copy. 
            e.setDropAction(QtCore.Qt.CopyAction) 
            e.accept() 
 
    def dragMoveEvent(self, e): 
        """ Handle a drag moving across the widget. """ 
 
        if self._dragged_object(e) is not None: 
            # Make sure the users knows we will only do a copy. 
            e.setDropAction(QtCore.Qt.CopyAction) 
            e.accept() 
 
    def dropEvent(self, e): 
        """ Handle a drop on the widget. """ 
 
        obj = self._dragged_object(e) 
        if obj is None: 
            return 
 
        # If we can't create a valid Python identifier for the name of an 
        # object we use this instead. 
        name = 'dragged' 
 
        if hasattr(obj, 'name') \
           and isinstance(obj.name, basestring) and len(obj.name) > 0: 
            py_name = python_name(obj.name) 
 
            # Make sure that the name is actually a valid Python identifier. 
            try: 
                if eval(py_name, {py_name : True}): 
                    name = py_name 
            except: 
                pass 
 
        self.interpreter.locals[name] = obj 
        self.run(name) 
        self.setFocus() 
 
        e.setDropAction(QtCore.Qt.CopyAction) 
        e.accept() 
 
    @staticmethod 
    def _dragged_object(e): 
        """Return the Python object being dragged or None if there isn't one. 
        """ 
 
        md = e.mimeData() 
 
        if isinstance(md, PyMimeData): 
            obj = md.instance() 
        else: 
            obj = None 
 
        return obj 
 
    def keyPressEvent(self, e): 
        """ Handle user input a key at a time. """ 
 
        text = e.text()
 
        if len(text) and 32 <= ord(toascii(text)[0]) < 127: 
            self._insert_text(text) 
            return 
 
        key = e.key() 
 
        if e.matches(QtGui.QKeySequence.Copy): 
            text = self.textCursor().selectedText()
            if len(text):
                QtGui.QApplication.clipboard().setText(text) 
        elif e.matches(QtGui.QKeySequence.Paste): 
            self._insert_text(QtGui.QApplication.clipboard().text())
            self.run(QtGui.QApplication.clipboard().text())
        elif key == QtCore.Qt.Key_Backspace: 
            if self._point: 
                cursor = self.textCursor() 
                cursor.deletePreviousChar() 
                self.setTextCursor(cursor) 
 
                self._point -= 1
                self._line = self._line [:self._point] + \
                             self._line[self._point+1:]
        elif key == QtCore.Qt.Key_Delete: 
            cursor = self.textCursor() 
            cursor.deleteChar() 
            self.setTextCursor(cursor) 
 
            self._line.remove(self._point, 1) 
        elif key == QtCore.Qt.Key_Return or key == QtCore.Qt.Key_Enter: 
            self.write('\n') 
            if self._reading: 
                self._reading = False 
            else: 
                self.run(str(self._line)) 
        elif key == QtCore.Qt.Key_Tab: 
            self._insert_text(text) 
        elif e.matches(QtGui.QKeySequence.MoveToPreviousChar): 
            if self._point: 
                self.moveCursor(QtGui.QTextCursor.Left) 
                self._point -= 1 
        elif e.matches(QtGui.QKeySequence.MoveToNextChar): 
            if self._point < len(self._line):
                self.moveCursor(QtGui.QTextCursor.Right) 
                self._point += 1 
        elif e.matches(QtGui.QKeySequence.MoveToStartOfLine): 
            while self._point: 
                self.moveCursor(QtGui.QTextCursor.Left) 
                self._point -= 1 
        elif e.matches(QtGui.QKeySequence.MoveToEndOfLine): 
            self.moveCursor(QtGui.QTextCursor.EndOfBlock) 
            self._point = len(self._line)
        elif e.matches(QtGui.QKeySequence.MoveToPreviousLine): 
            if len(self._history): 
                if self._pointer == 0: 
                    self._pointer = len(self._history) 
                self._pointer -= 1 
                self._recall() 
        elif e.matches(QtGui.QKeySequence.MoveToNextLine): 
            if len(self._history): 
                self._pointer += 1 
                if self._pointer == len(self._history): 
                    self._pointer = 0 
                self._recall() 
        else: 
            e.ignore() 
 
    def focusNextPrevChild(self, next): 
        """ Suppress tabbing to the next window in multi-line commands. """ 
 
        if next and self._more: 
            return False 
 
        return QtGui.QTextEdit.focusNextPrevChild(self, next) 
 
    def mousePressEvent(self, e): 
        """ 
        Keep the cursor after the last prompt. 
        """
        QtGui.QTextEdit.mousePressEvent(self, e)
        if e.button() == QtCore.Qt.LeftButton: 
            self.moveCursor(QtGui.QTextCursor.EndOfLine) 

    def mouseReleaseEvent(self, e):
        if e.button() == QtCore.Qt.MiddleButton:
            #self._lines.append(toascii(QtGui.QApplication.clipboard().text()))
            self.run(QtGui.QApplication.clipboard().text())
            #self._insert_text(QtGui.QApplication.clipboard().text())
 
    def contentsContextMenuEvent(self, ev): 
        """ Suppress the right button context menu. """ 
 
        pass 
 
    def _clear_line(self): 
        """ Clear the input line buffer. """ 
 
        self._line = ""
        self._point = 0 
         
    def _insert_text(self, text): 
        """ Insert text at the current cursor position. """ 
 
        self.insertPlainText(text) 
        self._line = self._line[:self._point]+ text + self._line[self._point:]
        self._point += len(text)
 
    def _recall(self): 
        """ Display the current item from the command history. """ 
 
        self.moveCursor(QtGui.QTextCursor.EndOfBlock) 
 
        while self._point: 
            self.moveCursor(QtGui.QTextCursor.Left, QtGui.QTextCursor.KeepAnchor) 
            self._point -= 1 
 
        self.textCursor().removeSelectedText() 
 
        self._clear_line() 
        self._insert_text(self._history[self._pointer]) 
 
if __name__=='__main__':
    app = QtGui.QApplication(sys.argv)
    def foo():
        print 'selectionChanged'
        
    QtGui.QApplication.clipboard().selectionChanged.connect(foo)
    sh = PyShell()
    sh.show()
    app.exec_()
    
