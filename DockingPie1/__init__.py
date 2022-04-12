###########################################################################
# Copyright (C) 2022 Serena Rosignoli, Alessandro Paiardini

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
###########################################################################


# Module first imported by the PyMOL plugin system. Contains:
#     - code to check for the presence of Python libraries needed for DockingPie to work.
#     - code to initialize DockingPie as a PyMOL plugin.


#--------------------------
# Check Python libraries. -
#--------------------------

import os
import sys
from importlib import util
import shutil

from pymol import cmd


# Gets the Python version.
python_version = sys.version_info.major
python_minor_version = "%s.%s" % (sys.version_info.major, sys.version_info.minor)
python_micro_version = "%s.%s.%s" % (sys.version_info.major, sys.version_info.minor, sys.version_info.micro)


# Checks for dependencies.
try:
    # Checks for some Qt bindings in PyMOL.
    from pymol.Qt import QtWidgets

    def showerror(title, message):
        QtWidgets.QMessageBox.critical(None, title, message)

    has_gui = "qt"
    pyqt_found = True

except ImportError:

    # Checks for Tkinter.
    try:
        if python_version == 3: # Python 3.
            from tkinter.messagebox import showerror
        else: # Python 2.
            from tkMessageBox import showerror
        has_gui = "tkinter"
    except ImportError: # On some open source builds, tkinter is missing.
        has_gui = None

    pyqt_found = False

try:
    import numpy
    numpy_found = True
except ImportError:
    numpy_found = False

try:
    import Bio
    biopython_found = True
except ImportError:
    biopython_found = False


# Sets the version of the plugin.
__docking_program_version__ = "1.0"
__revision__ = "1"
__version__ = float(__docking_program_version__ + __revision__.replace(".", ""))
docking_program_plugin_name = "DockingPie " + __docking_program_version__


#----------------------------------------
# Initialize the plugin in PyMOL. -
#----------------------------------------

def __init_plugin__(app):
    """
    Initializes the plugin in the plugin menu of PyMOL.
    """
    from pymol.plugins import addmenuitemqt
    addmenuitemqt(docking_program_plugin_name, startup_docking_program)


def startup_docking_program(app=None):

    """
    Executed when clicking on the plugin item in PyMOL's plugin menu.
    """
    QtWidgets.QApplication.setStyle('Fusion')

    if has_gui is None:
        print("\n# No GUI library (either Tkinter or Qt bindings) was found. DockingPie Plugin"
              " can not be launched.")
        return None

    # Check if a RxDock main window is already open.
    if pyqt_found:
        try:
            for widget in QtWidgets.QApplication.instance().topLevelWidgets():
                if hasattr(widget, "is_docking_program_main_window") and widget.isVisible():
                    title = "DockingPie Plugin Error"
                    message = ("DockingPie Plugin is already running. Please close its main"
                               " window or restart PyMOL in order to launch it again.")
                    showerror(title, message)
                    return None
        except Exception as e:
            pass

    # Checks if Python 3 is available.
    if python_version != 3:
        title = "Python Version Error"
        message = "DockingPie Plugin %s requires Python 3. Your current Python version is %s." % (__docking_program_version__, python_micro_version)
        showerror(title, message)
        return None

    # Checks the PyMOL version.
    pymol_version = float(".".join(cmd.get_version()[0].split(".")[0:2]))
    if pymol_version < 2.3:
        title = "PyMOL Version Error"
        message = "DockingPie Plugin %s requires a PyMOL version of 2.3 or higher. Your current PyMOL version is %s." % (__docking_program_version__, pymol_version)
        showerror(title, message)
        return None

    # Checks for PyQt.
    if not pyqt_found:
        title = "Import Error"
        message = "PyQt5 is not installed on your system. Please install it in order to use DockingPie Plugin."
        showerror(title, message)
        return None

    # Checks if NumPy and Biopython are present
    if not numpy_found:
        title = "Import Error"
        message = "NumPy is not installed on your system. Please install it in order to use DockingPie Plugin."
        showerror(title, message)
        return None

    if not biopython_found:
        title = "Import Error"
        message = "Biopython is not installed on your system. Please install it in order to use DockingPie Plugin."
        showerror(title, message)
        return None

    # Adds to the sys.path the directory where the plugin module is located.
    docking_program_plugin_dirpath = os.path.dirname(__file__)
    if os.path.isdir(docking_program_plugin_dirpath):
        sys.path.append(docking_program_plugin_dirpath)

    from lib import docking_program_main

    docking_program_main.docking_program_launcher(app=app,
                              docking_program_plugin_name=docking_program_plugin_name,
                              docking_program_version=__docking_program_version__,
                              docking_program_revision=__revision__)
