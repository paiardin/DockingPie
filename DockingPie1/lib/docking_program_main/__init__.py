# Copyright 2022 by Serena Rosignoli. All rights reserved.
# This code is part of the DockingPie package and governed by its license. Please
# see the LICENSE file.


import os
import sys
import shutil
import re
import json
import datetime

# PyMOL.
from pymol import cmd

from pymol.Qt import QtWidgets, QtCore, QtGui

import os
import shutil
import warnings
import math

from .docking_program_gui.main_window import DockingProgram_main_window_qt

def docking_program_launcher(app, docking_program_plugin_name, docking_program_version, docking_program_revision):
    docking_program = DockingProgram(app, docking_program_plugin_name, docking_program_version, docking_program_revision)

class DockingProgram():

    def __init__(self, app, docking_program_plugin_name, docking_program_version, docking_program_revision):

        self.docking_program_plugin_name = docking_program_plugin_name
        self.docking_program_version = docking_program_version
        self.docking_program_revision = docking_program_revision

        # Set to 'True' when developing, useful for debugging.
        self.DEVELOP = False
        # Set to 'True' to perform some tests on sequences/structures from the GUI.
        self.TEST = False

        self.app = app

         # If set to 'True' the most time consuming protocols will be run in a thread so
        # that the GUI is not freezed. When developing the code, it is better to set it to 'False',
        # in order to better track exceptions.
        if self.DEVELOP:
            self.use_protocol_threads = False
        else:
            self.use_protocol_threads = True


      #### MAIN WINDOW OF THE PLUGIN ####

        self.main_window = DockingProgram_main_window_qt(self)
