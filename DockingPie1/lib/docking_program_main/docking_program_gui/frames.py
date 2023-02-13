# Copyright 2022 by Serena Rosignoli. All rights reserved.
# This code is part of the DockingPie package and governed by its license. Please
# see the LICENSE file.


# Utilities
from pathlib import Path
import os
import sys
import shutil
import re
import json
import datetime
import fileinput
import warnings
import subprocess
import itertools
import time
import platform
import struct
import signal
import urllib.request

# PyMOL.
import pymol
from pymol import cmd
from pymol.cgo import *
from pymol.vfont import plain
from pymol.Qt import QtWidgets, QtCore, QtGui

from pymol import Qt
from pymol import stored
from pymol import viewing
from PyQt5.QtCore import Qt
from PyQt5.QtCore import QProcess

# Statistic
import math
import statistics

# Functionality modules
from lib.docking_program_main.Functions.pymol_interactions import Import_from_Pymol, PyMOLInteractions, ObjectParser
from lib.docking_program_main.Functions.pymol_interactions import *
from lib.docking_program_main.Functions.rxdock_functions import RxDock_Functions, RxDock_docking, RxDock_Cavity, RxDock_parse_results
from lib.docking_program_main.Functions.smina_functions import Smina_docking, Smina_parse_results
from lib.docking_program_main.Functions.adfr_functions import ADFR_docking, ADFR_parse_results
from lib.docking_program_main.Functions.vina_functions import Vina_docking, Vina_Parse_Results
from lib.docking_program_main.Functions.handle_widgets import HandleWidgets
from lib.docking_program_main.Functions.consensus_protocol import *
from lib.docking_program_main.docking_program_gui.new_windows import NewWindow, Import_from_pymol_window_qt
from lib.docking_program_main.docking_program_gui.dialogs import *
from lib.docking_program_main.Functions.threads import Protocol_exec_dialog
from lib.docking_program_main.Functions.general_docking_func import Generate_Object, Calculate_RMSD
from lib.docking_program_main.Functions.general_functions import OpenFromFile, Check_current_tab, Save_to_Csv, SelectAll
from lib.docking_program_main.Functions.installer import Installation, External_tools_installation_thread, External_tools_download_thread, External_components_dialog
from lib.docking_program_main.Functions.installer import Installation, External_tools_installation_thread, External_tools_download_thread, External_components_dialog
from lib.docking_program_main.docking_program_gui.frames import *
# Plotting Tools
from lib.docking_program_main.plots.plots_window import Plot_window_qt
from lib.docking_program_main.plots import pyqtgraph
import lib.docking_program_main.plots.plots_window as plot

# Table Tools
from lib.docking_program_main.tables.tables import *

# BioPython modules
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.PDBIO import Select
import Bio

# Numpy
import numpy as np

# RMSD calulating module
try:
    from spyrmsd import io, rmsd
except:
    pass

# csv module
import csv


class ConsensusFrame(QtWidgets.QFrame, PyMOLInteractions):

    def __init__(self, parent, main_window,
                 program,
                 *args, **configs):

        super(ConsensusFrame, self).__init__(main_window, *args, **configs)

        self.main_window = main_window

        self.program = program

        # Sets the layout of the frame.
        self.structure_frame_layout = QtWidgets.QGridLayout()
        self.setLayout(self.structure_frame_layout)
        self.tot_rows = 0

        # Builds a frame for each template structure and all its options.
        self.build_use_structure_frame()

        # Set the style
        self.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.structure_frame_layout.setAlignment(QtCore.Qt.AlignLeft)


    def build_use_structure_frame(self):

        self.name = QtWidgets.QLabel(str(self.program))

        self.label_weigth = QtWidgets.QLabel("Weigth: ")
        self.weight_slider = QtWidgets.QDoubleSpinBox()
        self.weight_slider.setRange(1, 5)
        self.weight_slider.setSingleStep(0.5)
        self.weight_slider.setValue(1)

        self.label_poses = QtWidgets.QLabel("Poses threshold: ")
        self.label_box = QtWidgets.QComboBox()
        # self.label_box.setRange(1, 100)
        # self.label_box.setSingleStep(1)
        # self.label_box.setValue(1)

        self.label_box.addItems(["All", "1", "3", "10"])

        self.structure_frame_layout.addWidget(self.name, 0, 0)

        self.structure_frame_layout.addWidget(self.label_poses, 1, 0)
        self.structure_frame_layout.addWidget(self.label_box, 1, 1)

        # self.structure_frame_layout.addWidget(self.label_weigth, 2, 0)
        # self.structure_frame_layout.addWidget(self.weight_slider, 2, 1)



class DockingsFrame(QtWidgets.QFrame, PyMOLInteractions):

    def __init__(self, parent, main_window,
                 results_name,
                 *args, **configs):

        super(DockingsFrame, self).__init__(main_window, *args, **configs)

        self.tab = parent

        self.results_name = results_name

        # Sets the layout of the frame.
        self.structure_frame_layout = QtWidgets.QGridLayout()
        self.setLayout(self.structure_frame_layout)
        self.tot_rows = 0

        # Builds a frame for each template structure and all its options.
        self.build_use_structure_frame()

        # Set the style
        self.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.structure_frame_layout.setAlignment(QtCore.Qt.AlignLeft)


    def build_use_structure_frame(self):

        # Checkboxes
        self.docking_checkbox = QtWidgets.QCheckBox(self.results_name)
        self.structure_frame_layout.addWidget(self.docking_checkbox, 0, 0)

        # Show details btn
        self.show_details_btn = QtWidgets.QRadioButton("Show Details")
        self.show_details_btn.clicked.connect(self.show_details_func)
        self.structure_frame_layout.addWidget(self.show_details_btn, 0, 1)

        # self.docking_run = self.main_window.docking_programs_child_tabs.docking_programs.all_runs[self.results_name]["docking_run"]

        # Docking run details
        ligand = self.tab.last_docking.ligand_to_dock
        receptor = self.tab.last_docking.receptor_to_dock
        cavity = self.tab.last_docking.cavity_name
        poses = str(cmd.count_states(self.results_name))
        self.details_label = QtWidgets.QLabel(str("Ligand: " + ligand +"\nReceptor: " + receptor + "\nCavity: " + cavity + "\nPoses: " + poses))
        self.structure_frame_layout.addWidget(self.details_label, 1, 1, 1, 2)
        self.details_label.hide()
        self.details_label.setStyleSheet("QLabel"
        "{"
        "font-size: 6px;"
        "}")


    def show_details_func(self):

        if self.show_details_btn.isChecked():
            self.details_label.show()
        else:
            self.details_label.hide()



class InstallationFrames(QtWidgets.QFrame, PyMOLInteractions):

    def __init__(self, parent, main_window,
                 program_name, run_configuration = False,
                 *args, **configs):

        super(InstallationFrames, self).__init__(main_window, *args, **configs)

        self.main_window = main_window

        self.program_name = program_name

        # Sets the layout of the frame.
        self.structure_frame_layout = QtWidgets.QGridLayout()
        self.setLayout(self.structure_frame_layout)
        self.tot_rows = 0

        # Builds a frame for each template structure and all its options.
        self.build_use_structure_frame()

        # Set the style
        self.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.structure_frame_layout.setAlignment(QtCore.Qt.AlignLeft)


    def build_use_structure_frame(self):

        self.program_name_label = QtWidgets.QLabel(self.program_name)

        self.installation_btn = QtWidgets.QPushButton("Install")

        self.installation_info_text_area = QtWidgets.QPlainTextEdit()
        self.installation_info_text_area.setReadOnly(True)

        config_path = self.main_window.docking_programs.config_path

        if self.program_name == "Smina":
            self.installation_info_text_area.setPlainText(open(str(os.path.join(config_path, "info_smina.txt")), encoding = "utf8").read())
        elif self.program_name == "Vina":
            self.installation_info_text_area.setPlainText(open(str(os.path.join(config_path, "info_vina.txt")), encoding = "utf8").read())
        elif self.program_name == "RxDock":
            self.installation_info_text_area.setPlainText(open(str(os.path.join(config_path, "info_rxdock.txt")), encoding = "utf8").read())
        elif self.program_name == "Openbabel":
            self.installation_info_text_area.setPlainText(open(str(os.path.join(config_path, "info_openbabel.txt")), encoding = "utf8").read())
        elif self.program_name == "sPyRMSD":
            self.installation_info_text_area.setPlainText(open(str(os.path.join(config_path, "info_spyrmsd.txt")), encoding = "utf8").read())
        elif self.program_name == "sdsorter":
            self.installation_info_text_area.setPlainText(open(str(os.path.join(config_path, "info_sdsorter.txt")), encoding = "utf8").read())
        elif self.program_name == "ADFR":
            self.installation_info_text_area.setPlainText(open(str(os.path.join(config_path, "info_adfr.txt")), encoding = "utf8").read())

        if self.program_name == "Config":
            self.configure_text = QtWidgets.QLabel("External Tools: ")
            self.configure_line_edit = QtWidgets.QLineEdit("Not Found")

            self.configure_btn = QtWidgets.QPushButton("Configure")
            self.show_info_btn = QtWidgets.QPushButton("Show Info")
            self.check_for_updates = QtWidgets.QPushButton("Check for Updates")

            self.show_info_btn.clicked.connect(self.show_detailed_info_window)
            self.configure_btn.clicked.connect(self.main_window.configure_external_tools_directory)
            self.check_for_updates.clicked.connect(self.check_for_updates_func)

            self.configure_line_edit.setEnabled(False)

            self.structure_frame_layout.addWidget(self.configure_text, 0, 0)
            self.structure_frame_layout.addWidget(self.configure_line_edit, 0, 2)
            self.structure_frame_layout.addWidget(self.configure_btn, 1, 0)
            self.structure_frame_layout.addWidget(self.show_info_btn, 1, 1)
            self.structure_frame_layout.addWidget(self.check_for_updates, 1, 2)

            self.check_external_tools()

        else:

            self.structure_frame_layout.addWidget(self.installation_info_text_area, 0, 2, 3, 1)
            self.structure_frame_layout.addWidget(self.program_name_label, 0, 0)
            self.structure_frame_layout.addWidget(self.installation_btn, 1, 0)


    def _get_path_string(self, path):
        _path = path
        if os.path.isdir(_path):
            return _path
        else:
            return _path + " (not found)"


    def get_python_architecture(self): # "32" for x86
        """
        Gets the architecture of the PyMOL built in which PyMod is running. "64" for x86-64.
        """
        return str(8*struct.calcsize("P"))


    def get_os_architecture(self):
        if sys.platform == "win32":
            if 'PROGRAMFILES(X86)' in os.environ:
                return "64"
            else:
                return "32"
        else:
            return platform.architecture()[0][0:2]


    def show_detailed_info_window(self):

        plugin_path = (os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

        try:
            import PyQt5
            pyqt5_version = PyQt5.QtCore.PYQT_VERSION_STR
        except:
            pyqt5_version = "-"

        try:
            from pymol import Qt
            pymol_pyqt_name = Qt.PYQT_NAME
        except:
            pymol_pyqt_name = "-"

        try:
            import Bio
            biopython_version = Bio.__version__
        except:
            biopython_version = "-"

        try:
            import numpy
            numpy_version = numpy.__version__
        except:
            numpy_version = "-"

        try:
            import conda
            import conda.cli.python_api as conda_api
            conda_version = conda.__version__
        except:
            conda_version = "-"
            conda_info_dict = {}
            conda_info_text = ""

        has_pymol_conda = str(hasattr(pymol, "externing") and hasattr(pymol.externing, "conda"))

        config_path = self.main_window.docking_programs.config_path
        with open(str(os.path.join(config_path, "version.txt"))) as f:
            read_version = f.readline().rstrip()

        version = (read_version.split("_")[1]).replace("v", "")


        additional_text = ("INFO and CONTACTS\n"
                           "Copyright (C): 2022 Serena Rosignoli, Alessandro Paiardini\n"
                           "Contacts: serena.rosignoli@uniroma1.it, alessandro.paiardini@uniroma1.it\n"
                           "For information about DockingPie visit:\n"
                           "https://github.com/paiardin/DockingPie\n\n"

                           "# DockingPie\n"
                           "- Version: " + version + "\n"
                           "- Plugin path: " + self._get_path_string(plugin_path) + " \n"
                           "- Config directory: " + self._get_path_string(self.main_window.docking_programs.config_path) + "\n\n"

                           "# PyMOL\n"
                           "- Version: " + str(cmd.get_version()[0]) + "\n"
                           "- Path: " + sys.executable + "\n"
                           "- Qt: " + str(pymol_pyqt_name) + "\n"
                           "- Has Conda: " + has_pymol_conda + "\n\n"

                           "# Python\n"
                           "- Version: " + str(sys.version) + "\n"
                           "- Arch: " + self.get_python_architecture() + "\n"
                           "- Path: " + sys.executable + "\n\n"

                           "# Operating system\n"
                           "- Platform: " + sys.platform + "\n"
                           "- Arch: " + self.get_os_architecture() + "\n\n"

                           "# Python libs\n"
                           "- PyQt5: " + pyqt5_version + "\n"
                           "- Conda version: " + conda_version + "\n"
                           "- Numpy version: " + numpy_version + "\n"
                           "- Biopython version: " + biopython_version + "\n"
                          )

        self.about_text_area = QtWidgets.QPlainTextEdit()
        self.about_text_area.setReadOnly(True)

        self.about_text_area.setPlainText(additional_text)

        self.detailed_info_window = NewWindow(parent = self.main_window,
        title = "DockingPie 1", upper_frame_title = "About",
        submit_command = None, submit_button_text= None,
        with_scroll = True)

        self.detailed_info_window.middle_layout_type.addWidget(self.about_text_area, 0, 0)

        self.detailed_info_window.show()

    #
    def check_external_tools(self):

        if sys.platform == "linux":
            dir_name = "external_tools_linux"
        if sys.platform == "darwin":
            dir_name = "external_tools_macOS"
        if sys.platform == "win32":
            dir_name = "external_tools_windows"

        ext_tools_path = os.path.join(self.main_window.docking_programs.dockingpie_extdir, dir_name)

        if os.path.isdir(ext_tools_path):
            self.configure_line_edit.setText(ext_tools_path)
            self.configure_btn.setEnabled(False)
            self.main_window.continue_check_installation = True
            self.check_for_updates.setEnabled(True)

        else:
            self.configure_line_edit.setText("Not Found")
            self.configure_btn.setEnabled(True)
            self.main_window.continue_check_installation = False
            self.check_for_updates.setEnabled(False)


    # def configure_external_tools_directory(self):
    #
    #     config_path = self.main_window.docking_programs.config_path
    #
    #     with open(str(os.path.join(config_path, "version.txt"))) as f:
    #         self.files_version = f.readline().rstrip()
    #
    #     try:
    #         urllib.request.urlopen("https://github.com/paiardin/DockingPie")
    #         github_accessible = True
    #     except:
    #         github_accessible = False
    #
    #     if github_accessible:
    #
    #         if sys.platform == "linux":
    #             dir_link = "https://github.com/paiardin/DockingPie/releases/download/" + str(self.files_version) + "/external_tools_linux.zip"
    #             dir_name = "external_tools_linux"
    #         if sys.platform == "darwin":
    #             dir_link = "https://github.com/paiardin/DockingPie/releases/download/" + str(self.files_version) + "/external_tools_macOS.zip"
    #             dir_name = "external_tools_macOS"
    #         if sys.platform == "win32":
    #             dir_link = "https://github.com/paiardin/DockingPie/releases/download/" + str(self.files_version) + "/external_tools_windows.zip"
    #             dir_name = "external_tools_windows"
    #
    #     else:
    #
    #         if sys.platform == "linux":
    #             dir_link = "http://schubert.bio.uniroma1.it/temp/DockingPie_config_file/" + str(self.files_version) + "/external_tools_linux.zip"
    #             dir_name = "external_tools_linux"
    #         if sys.platform == "darwin":
    #             dir_link = "http://schubert.bio.uniroma1.it/temp/DockingPie_config_file/" + str(self.files_version) +"/external_tools_macOS.zip"
    #             dir_name = "external_tools_macOS"
    #         if sys.platform == "win32":
    #             dir_link = "http://schubert.bio.uniroma1.it/temp/DockingPie_config_file/" + str(self.files_version) +"/external_tools_windows.zip"
    #             dir_name = "external_tools_windows"
    #
    #
    #     # Show the external components download dialog.
    #     install_dialog = External_components_dialog(self,
    #                                                 url=dir_link,
    #                                                 os_arch="64",
    #                                                 dir_name = dir_name,
    #                                                 local_mode=False)
    #     install_dialog.setModal(True)
    #     install_dialog.exec_()
    #
    #     # Finishes the installation.
    #     if install_dialog.complete_status:
    #         self.main_window.check_installation()
    #         self.check_external_tools()
    #

    def check_for_updates_func(self):


        ## This function check for external tools updates ##

        # Read current version of the plugin #
        config_path = self.main_window.docking_programs.config_path
        with open(str(os.path.join(config_path, "version.txt"))) as f:
            self.files_version = f.readline().rstrip()

        available_updates = False

        try:
            urllib.request.urlopen("https://github.com/paiardin/DockingPie")
            github_accessible = True
        except:
            github_accessible = False

        if github_accessible:

            # Download a file from GitHub were the version of the external tools is stored
            response = urllib.request.urlopen("https://github.com/paiardin/DockingPie/releases/download/versioning/last_version.txt")
            data = response.read()
            filename = os.path.join(self.main_window.docking_programs.config_path, "last_version.txt")
            file_ = open(filename, 'wb')
            file_.write(data)
            file_.close()

        else:

            # Download a file from GitHub were the version of the external tools is stored
            response = urllib.request.urlopen("http://schubert.bio.uniroma1.it/temp/DockingPie_config_file/" + str(self.files_version) + "/last_version.txt")
            data = response.read()
            filename = os.path.join(self.main_window.docking_programs.config_path, "last_version.txt")
            file_ = open(filename, 'wb')
            file_.write(data)
            file_.close()

        # Check last version available
        last_version_file = str(os.path.join(self.main_window.docking_programs.config_path, "last_version.txt"))
        with open(last_version_file) as f:
            last_version = f.readline().rstrip()

        current_version_file = str(os.path.join(self.main_window.docking_programs.config_path, "version.txt"))
        # Check current version
        with open(current_version_file) as f:
            current_version = f.readline().rstrip()

        # If the version is the same, delete the "last_version.txt" file
        if last_version == current_version:
            available_updates = False
            os.remove(last_version_file)

            QtWidgets.QMessageBox.about(self.main_window, "DockingPie", "External tools are up to date")

        # If the version is not the same
        else:

            # Ask if to update the external tools
            qm = QtWidgets.QMessageBox.question(self,'', "An new version of the external tools is available\n Do you want to update?", QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

            # If the user wants to update the external tools
            if qm == QtWidgets.QMessageBox.Yes:

                if sys.platform == "linux":
                    shutil.rmtree(str(os.path.join(self.main_window.docking_programs.dockingpie_extdir, "external_tools_linux")))
                if sys.platform == "win32":
                    shutil.rmtree(str(os.path.join(self.main_window.docking_programs.dockingpie_extdir, "external_tools_windows")))
                if sys.platform == "darwin":
                    shutil.rmtree(str(os.path.join(self.main_window.docking_programs.dockingpie_extdir, "external_tools_macOS")))

                # Update
                shutil.rmtree(str(os.path.join(self.main_window.docking_programs.dockingpie_extdir, "AutoDockTools")))
                os.remove(str(os.path.join(self.main_window.docking_programs.dockingpie_extdir, "prepare_flexreceptor4.py")))
                os.remove(str(os.path.join(self.main_window.docking_programs.dockingpie_extdir, "prepare_ligand4.py")))
                os.remove(str(os.path.join(self.main_window.docking_programs.dockingpie_extdir, "prepare_receptor4.py")))

                # Delete "version.txt" file
                os.remove(current_version_file)

                # Rename "last_version.txt" to "version.txt"
                os.rename(last_version_file, os.path.join(self.main_window.docking_programs.dockingpie_extdir, "version.txt"))

                self.configure_external_tools_directory()



            else: # If the user does not want to update the external tools

                # Delete the "last_version.txt" file
                os.remove(last_version_file)


class StructuresFrame(QtWidgets.QFrame, PyMOLInteractions):

    def __init__(self, parent, main_window,
                 strc_path, dict,
                 *args, **configs):

        super(StructuresFrame, self).__init__(main_window, *args, **configs)

        self.frame = main_window

        # name of the PyMOL obj, which is the same of the pdb file in tmp dir
        self.strc_path = strc_path

        # Dict of Loaded Objects
        self.dict = dict

        ### Dict keys are PyMOl objs names. parsed_object inside Dict contains Biopython parsed info about the PyMOl obj ###

        # Number of states detected in PyMOL
        self.num_states = str(self.dict[self.strc_path]["parsed_object"].num_states)

        # List of heteroresidues
        self.heteroresidues = self.dict[self.strc_path]["parsed_object"].heteroresidues_list

        # Only for receptor objects, it is detected if PROTEIN, RNA or DNA.
        self.receptor_type = self.dict[self.strc_path]["parsed_object"].receptor_type

        # Sets the layout of the frame.
        self.structure_frame_layout = QtWidgets.QGridLayout()
        self.setLayout(self.structure_frame_layout)
        self.tot_rows = 2

        Check_current_tab.check_docking_program_current_tab(self.frame)

        # Builds a frame for each template structure and all its options.
        self.build_use_structure_frame()

        # Set the style
        self.setFrameShape(QtWidgets.QFrame.StyledPanel)


    def build_use_structure_frame(self):

        #chbox_text = '\n'.join([self.strc_path, self.receptor_type]=
        # Create checkboxes corresponding to the structure
        self.strc_checkbox = QtWidgets.QCheckBox(self.strc_path)
        self.strc_checkbox.setStyleSheet('QCheckBox { font-size: 13pt}')
        self.strc_checkbox.toggled.connect(self.click_on_structure_checkbutton) # When checkbox is clicked it shows other options
        self.structure_frame_layout.addWidget(self.strc_checkbox, 0, 0)

        self.rec_type_label = QtWidgets.QLabel(self.receptor_type)
        self.structure_frame_layout.addWidget(self.rec_type_label, 1, 0, 4, 1)

        self.show_from_frame = QtWidgets.QPushButton("Zoom in PyMOL")
        self.structure_frame_layout.addWidget(self.show_from_frame, 0, 1)
        self.show_from_frame.clicked.connect(self.show_from_frame_func)
        self.show_from_frame.setEnabled(False)

        self.delete = QtWidgets.QPushButton("Remove")
        self.structure_frame_layout.addWidget(self.delete, 0, 2)
        self.delete.clicked.connect(self.delete_from_frame_func)
        self.delete.setEnabled(False)
        self.delete.setStyleSheet("QPushButton"
                             "{"
                             "background-color : rgb(255,51,51);"
                             "}"
                             "QPushButton::pressed"
                             "{"
                             "background-color : rgb(255,51,51);"
                             "}"
                             "QPushButton:disabled"
                             "{"
                             "background-color : rgb(160,160,160);"
                             "}"
                             )

        #self.delete.hide()

    def show_het_func(self):

        for i in self.het_text_list:
            if self.show_het_btn.isChecked() and self.show_het_btn.isEnabled():
                i.show()
            else:
                i.hide()

        if self.show_het_btn.isChecked() and self.show_het_btn.isEnabled():
            self.het_id_text.show()
        else:
            self.het_id_text.hide()


    def delete_from_frame_func(self):

        Check_current_tab.check_docking_program_current_tab(self.frame)

        if self.frame.is_vina_tab:
            HandleWidgets.remove_frame(self = self,
        dict = self.frame.docking_programs_child_tabs.docking_programs.vina_receptors_dict)

        if self.frame.is_smina_tab:
            HandleWidgets.remove_frame(self = self,
        dict = self.frame.docking_programs_child_tabs.docking_programs.smina_receptors_dict)

        if self.frame.is_rxdock_tab:
            HandleWidgets.remove_frame(self = self,
        dict = self.frame.docking_programs_child_tabs.docking_programs.rxdock_receptors_dict)


    def click_on_structure_checkbutton(self):

        Check_current_tab.check_docking_program_current_tab(self.frame)

        if self.get_use_as_template_var():

            self.show_from_frame.setEnabled(True)
            self.delete.setEnabled(True)
            self.show_het_btn.setEnabled(True)
            self.show_het_btn.setChecked(False)

            # if self.frame.is_vina_tab or self.frame.is_smina_tab:
            #
            #     self.add_h.setEnabled(True)
            #     self.remove_nonstd.setEnabled(True)
            #     self.remove_water.setEnabled(True)
            #
            # if self.frame.is_adfr_tab:
            #     self.add_h.setEnabled(False)
            #     self.add_h.setChecked(True)
            #     self.remove_nonstd.setEnabled(True)
            #     self.remove_water.setEnabled(True)

        else:
            self.show_from_frame.setEnabled(False)
            self.delete.setEnabled(False)
            self.show_het_btn.setEnabled(False)
            self.show_het_btn.setChecked(False)

            for i in self.het_text_list:
                i.hide()

            self.het_id_text.hide()

            # if self.frame.is_vina_tab or self.frame.is_smina_tab or self.frame.is_adfr_tab:
            #     self.add_h.setEnabled(False)
            #     self.remove_nonstd.setEnabled(False)
            #     self.remove_water.setEnabled(False)


    def get_use_as_template_var(self):
        return self.strc_checkbox.isChecked()


    def show_het_lig_list(self):
        pass


    def show_from_frame_func(self):
        text = self.strc_checkbox.text()
        new_text = text.split()[0]

        cmd.zoom(new_text)


class LigandsFrame(QtWidgets.QFrame, PyMOLInteractions):

    def __init__(self, parent, main_window,
                 strc_path, dict,
                 *args, **configs):

        super(LigandsFrame, self).__init__(main_window, *args, **configs)

        self.frame = main_window

        # name of the PyMOL obj, which is the same of the pdb file in tmp dir
        self.strc_path = strc_path

        # Dict of Loaded Objects
        self.dict = dict

        ### Dict keys are PyMOl objs names. parsed_object inside Dict contains Biopython parsed info about the PyMOl obj ###

        # Number of states detected in PyMOL
        self.num_states = str(self.dict[self.strc_path]["parsed_object"].num_states)

        # List of heteroresidues
        self.heteroresidues = self.dict[self.strc_path]["parsed_object"].heteroresidues_list

        # Only for receptor objects, it is detected if PROTEIN, RNA or DNA.
        self.receptor_type = self.dict[self.strc_path]["parsed_object"].receptor_type

        # Sets the layout of the frame.
        self.structure_frame_layout = QtWidgets.QGridLayout()
        self.setLayout(self.structure_frame_layout)
        self.tot_rows = 2

        Check_current_tab.check_docking_program_current_tab(self.frame)

        # Builds a frame for each template structure and all its options.
        self.build_use_structure_frame()

        # Set the style
        self.setFrameShape(QtWidgets.QFrame.StyledPanel)


    def build_use_structure_frame(self):

        self.checkbox_text = str(self.strc_path + " (" + self.num_states + ") " + self.receptor_type)

        self.strc_checkbox = QtWidgets.QCheckBox(self.checkbox_text) # Create checkboxes corresponding to the structure
        self.strc_checkbox.toggled.connect(self.click_on_structure_checkbutton) # When checkbox is clicked it shows other options
        self.structure_frame_layout.addWidget(self.strc_checkbox, 0, 0)

        self.show_from_frame = QtWidgets.QPushButton("Zoom in PyMOL")
        self.structure_frame_layout.addWidget(self.show_from_frame, 0, 1)
        self.show_from_frame.clicked.connect(self.show_from_frame_func)
        self.show_from_frame.setEnabled(False)

        self.delete = QtWidgets.QPushButton("Remove")
        self.structure_frame_layout.addWidget(self.delete, 0, 2)
        self.delete.clicked.connect(self.delete_from_frame_func)
        self.delete.setEnabled(False)
        self.delete.setStyleSheet("QPushButton"
                             "{"
                             "background-color : rgb(255,51,51);"
                             "}"
                             "QPushButton::pressed"
                             "{"
                             "background-color : rgb(255,51,51);"
                             "}"
                             "QPushButton:disabled"
                             "{"
                             "background-color : rgb(160,160,160);"
                             "}"
                             )
        #self.delete.hide()

    def delete_from_frame_func(self):

        Check_current_tab.check_docking_program_current_tab(self.frame)

        if self.frame.is_vina_tab:
            HandleWidgets.remove_frame(self = self,
        dict = self.frame.docking_programs_child_tabs.docking_programs.vina_ligands_dict)

        if self.frame.is_smina_tab:
            HandleWidgets.remove_frame(self = self,
        dict = self.frame.docking_programs_child_tabs.docking_programs.smina_ligands_dict)

        if self.frame.is_rxdock_tab:
            HandleWidgets.remove_frame(self = self,
        dict = self.frame.docking_programs_child_tabs.docking_programs.rxdock_ligands_dict)

        if self.frame.is_adfr_tab:
            HandleWidgets.remove_frame(self = self,
        dict = self.frame.docking_programs_child_tabs.docking_programs.adfr_ligands_dict)


    def click_on_structure_checkbutton(self):

        if self.get_use_as_template_var():
            self.show_from_frame.setEnabled(True)
            self.delete.setEnabled(True)

            # if self.frame.is_vina_tab or self.frame.is_smina_tab:
            #     self.add_h.setEnabled(True)
            #     self.add_h.setChecked(True)
            #     self.active_torsions_group.setEnabled(True)
            #     self.all_but_torsions.setChecked(True)
            #     self.all_but_torsions.setEnabled(True)
            #     self.all_torsions.setEnabled(True)
            #     self.none_torsions.setEnabled(True)
            #
            # if self.frame.is_adfr_tab:
            #     self.add_h.setEnabled(False)
            #     self.add_h.setChecked(True)
            #     self.active_torsions_group.setEnabled(True)
            #     self.all_but_torsions.setChecked(True)
            #     self.all_but_torsions.setEnabled(True)
            #     self.all_torsions.setEnabled(True)
            #     self.none_torsions.setEnabled(True)

        else:
            self.show_from_frame.setEnabled(False)
            self.delete.setEnabled(False)

            # if self.frame.is_vina_tab or self.frame.is_smina_tab or self.frame.is_adfr_tab:
            #     self.add_h.setEnabled(False)
            #     self.active_torsions_group.setEnabled(False)
            #     self.all_but_torsions.setEnabled(False)
            #     self.all_torsions.setEnabled(False)
            #     self.none_torsions.setEnabled(False)


    def get_use_as_template_var(self):
        return self.strc_checkbox.isChecked()


    def show_het_lig_list(self):
        pass


    def show_from_frame_func(self):
        text = self.strc_checkbox.text()
        new_text = text.split()[0]

        cmd.zoom(new_text)



class OptionsFrame(QtWidgets.QFrame, PyMOLInteractions):


    """
    A class to represent a Frame layout, in which docking options are specified.
    """


    def __init__(self, parent, main_window,
                 *args, **configs):

        super(OptionsFrame, self).__init__(main_window, *args, **configs)

        self.main_window = main_window

        # Sets the layout of the frame.
        self.options_frame_layout = QtWidgets.QGridLayout()
        self.setLayout(self.options_frame_layout)
        self.tot_rows = 0

        self.options_frame_layout.setAlignment(QtCore.Qt.AlignTop)
        self.setFrameShape(QtWidgets.QFrame.StyledPanel)



class OptionsFrameCanonical(QtWidgets.QFrame, PyMOLInteractions):


    """
    A class to represent a Frame layout, in which canonical docking options are specified.
    """


    def __init__(self, parent, main_window,
                 *args, **configs):

        super(OptionsFrameCanonical, self).__init__(main_window, *args, **configs)

        self.main_window = main_window

        # Sets the layout of the frame.
        self.options_frame_layout = QtWidgets.QGridLayout()
        self.setLayout(self.options_frame_layout)
        self.tot_rows = 0

        self.options_frame_layout.setAlignment(QtCore.Qt.AlignCenter)
        self.setFrameShape(QtWidgets.QFrame.StyledPanel)



class ResultsFrame(QtWidgets.QFrame, PyMOLInteractions):

    def __init__(self, parent, main_window, results_file, results_data,
                 results_file_name, results_obj_name, program, columns_names, log_file,
                 *args, **configs):

        super().__init__(parent, *args, **configs)

        self.docking_programs_child_tabs = main_window
        self.results_file = results_file # ResultFile object
        self.results_file_name = results_file_name # name of the results file with extension
        self.results_obj_name = results_obj_name # name of the results file w/o extension
        self.results_data = results_data
        self.log_file = log_file
        self.program = program
        self.summary_file = str(self.results_obj_name + "_summary.txt")

        self.has_RMSD = False

        # Sets the layout of the frame.
        self.results_frame_layout = QtWidgets.QGridLayout()
        self.setLayout(self.results_frame_layout)

        if self.results_file and self.results_data is not None:
            # for plotting
            if self.program == "Vina":
                self.dict = self.docking_programs_child_tabs.docking_programs.vina_runs_dict
                self.score_header = "Affinity"

            elif self.program == "Smina":
                self.dict = self.docking_programs_child_tabs.docking_programs.smina_runs_dict
                self.score_header = "Affinity"

            elif self.program == "RxDock":
                self.dict = self.docking_programs_child_tabs.docking_programs.rxdock_runs_dict
                self.score_header = "SCORE"

            elif self.program == "ADFR":
                self.dict = self.docking_programs_child_tabs.docking_programs.adfr_runs_dict
                self.score_header = "Affinity"

            self.view_pose_btn_list = []
            # Create a table with the score results
            for data in self.results_data:

                self.results_table = TableView(parent=self.docking_programs_child_tabs, data=self.results_data,
                                           row_labels=None, row_labels_height=25,
                                           column_labels=columns_names,
                                           sortable=False)

                self.results_table.itemDoubleClicked.connect(self.on_click_results_table)

            # Create a label for the docking process
            label = QtWidgets.QLabel(self.results_obj_name)

            # Add the Table to the results frame
            self.results_frame_layout.addWidget(label, 0, 0)
            self.results_frame_layout.addWidget(self.results_table, 1, 0)

            self.options_frame_all = OptionsFrame(parent=None,
            main_window=self.docking_programs_child_tabs)
            self.results_frame_layout.addWidget(self.options_frame_all, 1, 1)

            # RMSD pushbutton
            self.rmsd_checkbox = QtWidgets.QPushButton("Calculate RMSD")
            self.rmsd_checkbox.clicked.connect(self.show_results_advanced_options_func)
            self.options_frame_all.options_frame_layout.addWidget(self.rmsd_checkbox, 2, 0)

            # Save to file
            self.save_to_file = QtWidgets.QPushButton("Save to File")
            self.save_to_file.clicked.connect(self.save_to_file_choise_func)
            self.options_frame_all.options_frame_layout.addWidget(self.save_to_file, 3, 0)

            # Plot
            self.single_plot_btn = QtWidgets.QPushButton("Plot")
            self.single_plot_btn.clicked.connect(self.open_plot_window_choice)
            self.options_frame_all.options_frame_layout.addWidget(self.single_plot_btn, 4, 0)

            # Show LOG button
            self.log_btn = QtWidgets.QPushButton("LOG files")
            self.log_btn.clicked.connect(self.show_log_file_window)
            self.options_frame_all.options_frame_layout.addWidget(self.log_btn, 6, 0)

            # # Plot all (out of OptionsFrame)
            # self.all_plot_btn = QtWidgets.QPushButton("Plot All")
            # self.all_plot_btn.clicked.connect(self.plot_run_score_scatter)
            # #results_frame_layout
            # self.options_frame_all.options_frame_layout.addWidget(self.all_plot_btn, 5, 0)

        else:
            self.results_table = None

            self.options_frame_all = OptionsFrame(parent=None,
            main_window=self.docking_programs_child_tabs)
            self.results_frame_layout.addWidget(self.options_frame_all, 1, 1)

            if self.log_file is not None:

                if Path(self.log_file).is_file():

                    self.log_label = QtWidgets.QLabel("Log file")
                    self.options_frame_all.options_frame_layout.addWidget(self.log_label, 0, 0)
                    self.log_text_area = QtWidgets.QPlainTextEdit()
                    self.options_frame_all.options_frame_layout.addWidget(self.log_text_area, 1, 0)
                    self.log_text_area.setPlainText(open(str(self.log_file)).read())
                    self.log_text_area.setReadOnly(True)

                if Path(self.summary_file).is_file():

                    self.summary_label = QtWidgets.QLabel("Summary")
                    self.options_frame_all.options_frame_layout.addWidget(self.summary_label, 0, 1)
                    self.summary_text_area = QtWidgets.QPlainTextEdit()
                    self.options_frame_all.options_frame_layout.addWidget(self.summary_text_area, 1, 1)
                    self.summary_text_area.setPlainText(open(str(self.summary_file)).read())
                    self.summary_text_area.setReadOnly(True)

            else:

                if Path(self.summary_file).is_file():

                    self.summary_label = QtWidgets.QLabel("Summary")
                    self.options_frame_all.options_frame_layout.addWidget(self.summary_label, 0, 1)
                    self.summary_text_area = QtWidgets.QPlainTextEdit()
                    self.options_frame_all.options_frame_layout.addWidget(self.summary_text_area, 1, 1)
                    self.summary_text_area.setPlainText(open(str(self.summary_file)).read())
                    self.summary_text_area.setReadOnly(True)



    def show_log_file_window(self):

        self.log_window = NewWindow(parent = self.docking_programs_child_tabs.docking_programs,
        title = "LOG", upper_frame_title = "",
        submit_command = self.close_window, submit_button_text= "Close",
        with_scroll = True)

        if self.program == "Vina":
            dir = self.docking_programs_child_tabs.docking_programs.vina_tmp_dir

        elif self.program == "Smina":
            dir = self.docking_programs_child_tabs.docking_programs.smina_tmp_dir

        elif self.program == "RxDock":
            dir = self.docking_programs_child_tabs.docking_programs.rxdock_tmp_dir

        elif self.program == "ADFR":
            dir = self.docking_programs_child_tabs.docking_programs.adfr_tmp_dir

        os.chdir(dir)

        if self.log_file is not None:

            self.log_label = QtWidgets.QLabel("Log file")
            self.log_window.middle_layout_type.addWidget(self.log_label, 0, 0)
            self.log_text_area = QtWidgets.QPlainTextEdit()
            self.log_window.middle_layout_type.addWidget(self.log_text_area, 1, 0)
            self.log_text_area.setPlainText(open(str(self.log_file)).read())
            self.log_text_area.setReadOnly(True)

        self.summary_label = QtWidgets.QLabel("Summary")
        self.log_window.middle_layout_type.addWidget(self.summary_label, 0, 1)
        self.summary_text_area = QtWidgets.QPlainTextEdit()
        self.log_window.middle_layout_type.addWidget(self.summary_text_area, 1, 1)
        self.summary_text_area.setPlainText(open(str(self.summary_file)).read())
        self.summary_text_area.setReadOnly(True)

        self.log_window.show()


    def close_window(self):
        self.log_window.close()


    def on_click_results_table(self):

        row_index = self.results_table.currentIndex().row()
        #row_index = index.row()

        column_index = (TableView.get_column_index_from_header(self, table = self.results_table, header = "NAME"))[0]

        model = self.results_table.model()
        index = model.index(row_index, column_index)
        name = model.data(index)
        HideEverythingPyMOL(self, to_show = [self.results_obj_name, self.docking_programs_child_tabs.docking_tab_ui.last_docking.receptor_to_dock])
        cmd.frame(row_index+1)
        cmd.orient(self.results_obj_name)


    def plot_run_score_scatter(self):

        # Initialize plotting area
        cp = Plot_window_qt(self)
        cp.initialize(parent=self, title="Score vs Runs")
        cp.build_plotting_area(use_controls=True,
                           x_label_text="Runs", y_label_text= self.score_header,
                           hide_y_ticks=False,
                           use_all_controls_buttons=True,
                           use_save_to_csv=False,
                           discrete_x_ticks=True,
                           discrete_x_ticks_label = list(self.dict.keys()))

        for i, key in enumerate(self.dict):
            # Get table
            table = self.dict[key]["results_table"]
            # Get index of "Score/Affinity" column
            y_index = TableView.get_column_index_from_header(self, table = table, header = self.score_header)
            # Get data from index
            y = TableView.get_column_data_from_index(self, table = table, column_header_index = y_index[0])
            # Create repeats for each run. Workaround for discrete variables
            num = list(itertools.repeat(i, len(y)))
            # Add a ScatterPlot for each Docking Run
            cp.add_scatter_plot(num, y, label = key, corr=False)

        cp.show()


    def open_plot_window_choice(self):

        self.plot_rmsd_score_scatter()


    def plot_rmsd_score_scatter(self):

        if self.has_RMSD == False:
            QtWidgets.QMessageBox.warning(self.docking_programs_child_tabs, "RMSD not found", "To plot the results, please calculate the RMSD for the Run of interest")

        else:
            # Get index of "Score/Affinity" column
            x_index = TableView.get_column_index_from_header(self, table = self.results_table, header = self.score_header)
            # Get data from index
            x = TableView.get_column_data_from_index(self, table = self.results_table, column_header_index = int(x_index[0]))
            # Get index of "RMSD" column
            y_index = TableView.get_column_index_from_header(self, table = self.results_table, header = "RMSD")


            # Initialize plotting area
            cp = Plot_window_qt(self)
            cp.initialize(parent=self, title=str("RMSD vs " + self.score_header))
            cp.build_plotting_area(use_controls=True,
                               x_label_text= self.score_header, y_label_text="RMSD",
                               hide_y_ticks=False,
                               use_all_controls_buttons=True,
                               use_save_to_csv=False,
                               update_messagebar = True)


            for y_ind in y_index:
                # Get data from index
                y = TableView.get_column_data_from_index(self, table = self.results_table, column_header_index = int(y_ind))
                # Add a ScatterPlot for each RMSD column found
                cp.add_scatter_plot(x, y,
                label = self.results_table.column_labels[y_ind],
                corr = True,
                additional_data = x)

            cp.show()


    def save_to_file_choise_func(self):

        self.save_to_file_choise_window = NewWindow(parent = self.docking_programs_child_tabs,
        title = "Save to File", upper_frame_title = "Select the saving option",
        submit_command = self.get_saving_choice, submit_button_text= "Start",
        with_scroll = True)

        self.save_table_to_csv = QtWidgets.QRadioButton("Save table to \".csv\" file")
        self.save_docking_files = QtWidgets.QRadioButton("Save Docking results file")
        self.save_log = QtWidgets.QRadioButton("Save Log File")

        # Add Options to the New Window
        self.save_to_file_choise_window.middle_layout_type.addWidget(self.save_table_to_csv, 0, 0)
        self.save_to_file_choise_window.middle_layout_type.addWidget(self.save_docking_files, 1, 0)
        self.save_to_file_choise_window.middle_layout_type.addWidget(self.save_log, 2, 0)

        self.save_to_file_choise_window.show()


    def get_saving_choice(self):

        if self.save_docking_files.isChecked():

            if self.program == "Vina":
                # Save the grd file
                filepath = asksaveasfile_qt("Save Docking Results", name_filter="*.pdbqt")

                if not filepath:
                    return None

                else:
                    os.chdir(self.docking_programs_child_tabs.docking_programs.vina_tmp_dir)
                    shutil.copyfile(self.results_file_name, filepath)

            elif self.program == "RxDock":

                filepath = asksaveasfile_qt("Save Docking Results", name_filter="*.sd")

                if not filepath:
                    return None

                else:
                    os.chdir(self.docking_programs_child_tabs.docking_programs.rxdock_tmp_dir)
                    shutil.copyfile(self.results_file_name, filepath)


            elif self.program == "Smina":

                filepath = asksaveasfile_qt("Save Docking Results", name_filter="*.pdbqt")

                if not filepath:
                    return None

                else:
                    os.chdir(self.docking_programs_child_tabs.docking_programs.smina_tmp_dir)
                    shutil.copyfile(self.results_file_name, filepath)

            elif self.program == "ADFR":

                filepath = asksaveasfile_qt("Save Docking Results", name_filter="*.pdbqt")

                if not filepath:
                    return None

                else:
                    os.chdir(self.docking_programs_child_tabs.docking_programs.adfr_tmp_dir)
                    shutil.copyfile(self.results_file_name, filepath)


        if self.save_table_to_csv.isChecked():

            Save_to_Csv.save_to_csv_event(self, table = self.results_table)


        if self.save_log.isChecked():

            if self.program == "Vina":
                # Save the grd file
                filepath = asksaveasfile_qt("Save Log File", name_filter="*.txt")

                if not filepath:
                    return None

                else:
                    os.chdir(self.docking_programs_child_tabs.docking_programs.vina_tmp_dir)
                    shutil.copyfile(self.log_file, filepath)

            elif self.program == "RxDock":

                 print("RxDock LOG File is not avilable")

            elif self.program == "Smina":

                filepath = asksaveasfile_qt("Save Log File", name_filter="*.txt")

                if not filepath:
                    return None

                else:
                    os.chdir(self.docking_programs_child_tabs.docking_programs.smina_tmp_dir)
                    shutil.copyfile(self.log_file, filepath)

            elif self.program == "ADFR":

                filepath = asksaveasfile_qt("Save Log File", name_filter="*.txt")

                if not filepath:
                    return None

                else:
                    os.chdir(self.docking_programs_child_tabs.docking_programs.adfr_tmp_dir)
                    shutil.copyfile(self.log_file, filepath)


        self.save_to_file_choise_window.close()


    def show_results_advanced_options_func(self):

        if self.program == "Smina":
            path = self.docking_programs_child_tabs.docking_programs.smina_tmp_dir
            # The name of the Docking Result obj in PyMOL (ex. Run_0)
            self.reference_obj_name = self.docking_programs_child_tabs.docking_programs.SMINA.docking_tab_ui.results_obj_name
            self.PROGRAM = self.docking_programs_child_tabs.docking_programs.SMINA

        if self.program == "RxDock":
            path = self.docking_programs_child_tabs.docking_programs.rxdock_tmp_dir
            # The name of the Docking Result obj in PyMOL (ex. Run_0)
            self.reference_obj_name = self.docking_programs_child_tabs.docking_programs.RXDOCK.docking_tab_ui.results_obj_name
            self.PROGRAM = self.docking_programs_child_tabs.docking_programs.RXDOCK

        if self.program == "Vina":
            path = self.docking_programs_child_tabs.docking_programs.vina_tmp_dir
            # The name of the Docking Result obj in PyMOL (ex. Run_0)
            self.reference_obj_name = self.docking_programs_child_tabs.docking_programs.VINA.docking_tab_ui.results_obj_name
            self.PROGRAM = self.docking_programs_child_tabs.docking_programs.VINA

        if self.program == "ADFR":
            path = self.docking_programs_child_tabs.docking_programs.adfr_tmp_dir
            # The name of the Docking Result obj in PyMOL (ex. Run_0)
            self.reference_obj_name = self.docking_programs_child_tabs.docking_programs.ADFR.docking_tab_ui.results_obj_name
            self.PROGRAM = self.docking_programs_child_tabs.docking_programs.ADFR

        os.chdir(path)

        self.results_advanced_options_window = NewWindow(parent = self.docking_programs_child_tabs,
        title = "Advanced Options", upper_frame_title = "Set Options to calculate RMSD",
        submit_command = self.calculate_rmsd_func, submit_button_text= "Calculate RMSD",
        with_scroll = True)

        self.docked_objs_label = QtWidgets.QLabel("Choose a docked ligand")
        self.results_advanced_options_window.middle_layout_type.addWidget(self.docked_objs_label)

        self.docked_objs_combobox = QtWidgets.QComboBox()
        self.results_advanced_options_window.middle_layout_type.addWidget(self.docked_objs_combobox)
        self.fill_docked_objs_combobox()

        self.rmsd_label = QtWidgets.QLabel("Choose an object to use as a reference")
        self.results_advanced_options_window.middle_layout_type.addWidget(self.rmsd_label)

        self.rmsd_pushbtn = QtWidgets.QPushButton("Import from PyMOL")
        self.results_advanced_options_window.middle_layout_type.addWidget(self.rmsd_pushbtn)
        self.rmsd_pushbtn.clicked.connect(self.import_from_pymol)

        self.rmsd_combobox = QtWidgets.QComboBox()
        self.results_advanced_options_window.middle_layout_type.addWidget(self.rmsd_combobox)

        self.results_advanced_options_window.show()


    def fill_docked_objs_combobox(self):

        self.docked_objs_combobox.addItem(self.results_obj_name)

    def calculate_rmsd_func(self):

        # The state of the Docking Result obj in PyMOL we are interested in. (ex. obj01_pose_1 --> reference_obj_state = 1)
        # self.reference_obj_state = self.docked_objs_combobox.currentText().split("_")[2]

        # The name of the reference file saved in tmp dir
        self.reference_obj_name_state = str(self.reference_obj_name + "_reference.pdb")

        # The name of the docked file saved in tmp dir
        self.docked_obj_name = str(self.reference_obj_name + "_docked.sdf")

        HandleWidgets.combobox_check_if_empty(self=self,
        widgets_list=[self.rmsd_combobox])

        if self.is_empty:
            QtWidgets.QMessageBox.warning(self.docking_programs_child_tabs, "Warning", "There aren't enough parameters")

        else:
            self.save_and_calculate_RMSD()


    def save_and_calculate_RMSD(self):

        self.rmsd_array = [[]]

        cmd.save(self.reference_obj_name_state, self.rmsd_combobox.currentText(), format = "pdb")
        cmd.save(self.docked_obj_name, self.docked_objs_combobox.currentText(), state = 0, format = "sdf")

        self.rmsd = Calculate_RMSD(self, self.reference_obj_name_state, self.docked_obj_name)

        if self.rmsd.rmsd_computed:

            # list_of_list = []
            #
            # for i in range(len(self.rmsd.rmsd_list)):
            #     new_list = [self.rmsd_combobox.currentText(), str(self.rmsd.rmsd_list[i])]
            #     list_of_list.append(new_list)
            #
            # self.rmsd_array = list_of_list


            listof = []
            for i in range(len(self.rmsd.rmsd_list)):
                new_list = [str(self.rmsd.rmsd_list[i])]
                listof.append(new_list)

            column_position = self.results_table.columnCount()
            self.results_table.insertColumn(column_position)

            TableView.add_data(self, existing_table = self.results_table,
            new_data = listof,
            row_index = 0,
            column_index = column_position,
            new_row_label = None,
            new_column_label = [str("RMSD\nReference: " + self.rmsd_combobox.currentText())],
            by_row = False)

            # Update has_RMSD
            self.has_RMSD = True

            # Update dict with info about Docking runs
            self.docking_programs_child_tabs.docking_programs.all_runs[self.results_obj_name]["rmsd_list"] = self.rmsd.rmsd_list

        else:
            pass

        # close window
        self.results_advanced_options_window.close()

    def rmsd_boxplot_func(self):
        pass


    def import_from_pymol(self):

        # Create a list of all the importable Objects and selections from PyMOL
        self.selections = [str(obj) for obj in cmd.get_names("objects") + cmd.get_names("selections")]

        # Check for importable objects
        if self.selections == []:
            QtWidgets.QMessageBox.warning(self.docking_programs_child_tabs, "PyMOL is empty", str("There isn't any object to import"))

        # If present, update the Combo Box
        else:
            self.rmsd_combobox.clear()
            for i in self.selections:
                type = cmd.get_type(i)
                if type == str("object:molecule") and cmd.count_atoms(i) < 99:
                    self.rmsd_combobox.addItem(i)
                elif type == str("selection") and cmd.count_atoms(i) < 99:
                    self.rmsd_combobox.addItem(i)
                else:
                    pass


    def view_pose_func(self):
        # index = self.results_frame_layout.indexOf(self.view_pose_btn)
        # row, column, cols, rows = self.results_frame_layout.getItemPosition(index)
        # self.results_frame_layout.itemSelected.emit(self.view_pose_btn, pos(row, column))
        sender = self.sender()
        self.index = self.view_pose_btn_list.index(sender)
        cmd.zoom(self.results_obj_name, 0.0, int(self.index), 0)
