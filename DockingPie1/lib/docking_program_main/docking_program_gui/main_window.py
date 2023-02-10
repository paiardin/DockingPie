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
import pymol
from pymol import cmd

from pymol.Qt import QtWidgets, QtCore, QtGui

import os
import shutil
import warnings
import math

from .tabs import GridTab_RxDock
from .tabs import GridTab
from .tabs import ReceptorTab
from .tabs import LigandTab
from .tabs import DockingTab
from .tabs import ResultsTab
from .tabs import ConfigurationTab, DataAnalysisTab

def welcomewindow(title, message, parent=None, buttons_text=None):
    """
    Wrapper to a Yes/no dialog in PyQt. If 'buttons_text' is 'None', the default
    "Yes" and "No" buttons will be used. If 'buttons_text' is a list with two
    strings, the first string will be the text of the "Yes" button and the second
    one will be the text of the "No" button.
    """

    dialog = QtWidgets.QMessageBox(parent)
    dialog.setWindowTitle(title)
    dialog.setText(message)

    count = dialog.layout().count()

    # for i in range(count):
    #     item = dialog.layout().itemAt(i)
    #     print(item)

    if len(buttons_text) == 1:
        start_button = dialog.addButton(buttons_text[0], QtWidgets.QMessageBox.AcceptRole)
        dialog.setStyleSheet('QLabel { font-size: 13pt; padding: 15px}')
        answer = dialog.exec_()

    else:

        yesbutton = dialog.addButton(buttons_text[0], QtWidgets.QMessageBox.YesRole)
        nobutton = dialog.addButton(buttons_text[1], QtWidgets.QMessageBox.NoRole)
        dialog.setStyleSheet('QLabel { font-size: 13pt; padding: 15px}')
        yesbutton.setStyleSheet('QPushButton {font-weight: bold; border-color: green}')

        answer = dialog.exec_()
        return dialog.clickedButton() is yesbutton






class DockingProgram_main_window_main_menu:

    is_rxdock_main_window = True

    def make_main_menu(self):

        """
        A method to create the Main Window Main Menu.

        note: it is currently not used.
        """

        self.menubar = self.menuBar()
        self.menubar.setNativeMenuBar(False)

        #---------------
        # "File" menu. -
        #---------------

        self.file_menu = self.menubar.addMenu('File')

        # Workspace submenu.
        self.sessions_submenu = QtWidgets.QMenu('Sessions', self)
        self.file_menu.addMenu(self.sessions_submenu)

        self.file_menu.addSeparator()
        self.exit_submenu = QtWidgets.QMenu('Exit', self)
        self.file_menu.addMenu(self.exit_submenu)



class DockingProgram_main_window_qt(QtWidgets.QMainWindow, DockingProgram_main_window_main_menu):

    is_docking_program_main_window = True

    def __init__(self, docking_program, parent=None):
        super(DockingProgram_main_window_qt, self).__init__(parent)

        """

        MAIN WINDOW

        """

        # Initial settings.
        self.docking_program = docking_program
        self.title = self.docking_program.docking_program_plugin_name + "." + self.docking_program.docking_program_revision
        self.statusBar_message = "Welcome to DockingPie"
        self.left = 550
        self.top = 50
        self.width = 550
        self.height = 400

        self.vertical_spacing = 1

        # Creating a menu bar.
        #self.make_main_menu()

        # Create the Tabwidget for all the Docking Programs and Configuration Tab. Each Docking Program has its own Tab.
        self.main_docking_programs_tabs = Docking_Programs(self)
        main_tab_bar = self.main_docking_programs_tabs.docking_programs_tabs.tabBar()

        # Creating central widget and setting it as central.
        self.central_widget = Centralwid(self)
        self.setCentralWidget(self.central_widget)

        # Add the Docking Program Tab Widget and the Bottom Widget to the Central Widget
        self.central_widget.central_layout.addWidget(self.main_docking_programs_tabs.docking_programs_tabs)
        self.central_widget.central_layout.addWidget(self.central_widget.bottom_widg)

        # Set the layout.
        self.central_widget.central_layout.setFormAlignment(QtCore.Qt.AlignLeft)
        self.central_widget.central_layout.setVerticalSpacing(self.vertical_spacing)
        self.central_widget.setLayout(self.central_widget.central_layout)

        # Creating status bar.
        self.statusBar().showMessage(self.statusBar_message, 3000)
        self.statusBar().setSizeGripEnabled(1)

        if os.path.isdir(self.main_docking_programs_tabs.tmp_dir_path):
            pass
        else:
            os.mkdir(self.main_docking_programs_tabs.tmp_dir_path)
            os.mkdir(self.main_docking_programs_tabs.general_tmp_dir)
            os.mkdir(self.main_docking_programs_tabs.rxdock_tmp_dir)
            os.mkdir(self.main_docking_programs_tabs.vina_tmp_dir)
            os.mkdir(self.main_docking_programs_tabs.smina_tmp_dir)
            os.mkdir(self.main_docking_programs_tabs.adfr_tmp_dir)
            os.mkdir(self.main_docking_programs_tabs.consensus_tmp_dir)

        # Clean and make a tmp directory
        shutil.rmtree(self.main_docking_programs_tabs.general_tmp_dir)
        os.mkdir(self.main_docking_programs_tabs.general_tmp_dir)

        shutil.rmtree(self.main_docking_programs_tabs.rxdock_tmp_dir)
        os.mkdir(self.main_docking_programs_tabs.rxdock_tmp_dir)

        shutil.rmtree(self.main_docking_programs_tabs.vina_tmp_dir)
        os.mkdir(self.main_docking_programs_tabs.vina_tmp_dir)

        shutil.rmtree(self.main_docking_programs_tabs.smina_tmp_dir)
        os.mkdir(self.main_docking_programs_tabs.smina_tmp_dir)

        shutil.rmtree(self.main_docking_programs_tabs.adfr_tmp_dir)
        os.mkdir(self.main_docking_programs_tabs.adfr_tmp_dir)

        shutil.rmtree(self.main_docking_programs_tabs.consensus_tmp_dir)
        os.mkdir(self.main_docking_programs_tabs.consensus_tmp_dir)

        # Initialize User Interface.
        self.initUI()

        ## Check if the configuration is completed
        configuration = self.check_configuration(self.main_docking_programs_tabs)
        configuration_completed = configuration[1]
        self.run_configuration = False
# """ We are thrilled to have you as a user. Our goal is to provide you with a seamless and efficient experience with Molecular Docking. If you need any assistance, please contact us at: serena.rosignoli@uniroma1.it. Follow our GitHub to be always up to date with the new releases and bug fixes: https://github.com/paiardin/DockingPie Thank you for choosing DockingPie!"""

        if not configuration_completed:
            ## Check if all the dependencies are installed - only if the configuration is completed
            # TODO
            print(configuration_completed)

            # Add Welcome Message to CONFIGURATION TAB
            title = "Welcome to DockingPie"
            message = """We are thrilled to have you as a user! Our goal is to provide you with a seamless and efficient experience with Molecular Docking.

It seems that DockingPie need to be configured,
do you want to proceed? """
            choice = welcomewindow(title,
                                 message,
                                 parent=self,
                                 buttons_text = ["Yes, CONFIGURE IT! (Best choice)", "No, I will do it later (You will forget ...)"])

            if choice:
                self.run_configuration = True
                # print(dir(self.main_docking_programs_tabs.docking_programs_tabs.tabText(0)))
                # print(self.main_docking_programs_tabs.docking_programs_tabs.tabBar().tabText(0))
                # #self.main_docking_programs_tabs.CONFIGURATION.installation_frame.configure_external_tools_directory()

        else:
            title = "Welcome to DockingPie"
            message = """We are thrilled to have you as a user! Our goal is to provide you with a seamless and efficient experience with Molecular Docking.

DockingPie is properly configured!
Please check that any dependency is installed before continuing.
 """

            choice = welcomewindow(title,
                                 message,
                                 parent=self,
                                 buttons_text = ["Let's Start"])

        if self.run_configuration:
            self.main_docking_programs_tabs.CONFIGURATION_TAB.configure_external_tools_directory()

    def check_configuration(self, main):

        configuration = {"et": False,
                        "adt": False,
                        "vina_tree": False,
                        "adfr_tree": False,
                        "agfr_tree": False,
                        "smina_tree": False,
                        "sdsorter_tree": False,
                        "pfr": False,
                        "pl": False,
                        "pr": False}

        config_path = main.config_path

        path_to_adt = os.path.join(config_path, "AutoDockTools")
        path_to_et = os.path.join(config_path, main.et)
        path_to_pfr = os.path.join(config_path, "prepare_flexreceptor4.py")
        path_to_pl = os.path.join(config_path, "prepare_ligand4.py")
        path_to_pr = os.path.join(config_path, "prepare_receptor4.py")

        if os.path.isdir(path_to_adt):
            configuration["adt"] = True
        if os.path.isdir(path_to_et):
            configuration["et"] = True

        if os.path.isfile(path_to_pfr):
            configuration["pfr"] = True
        if os.path.isfile(path_to_pl):
            configuration["pl"] = True
        if os.path.isfile(path_to_pr):
            configuration["pr"] = True

        vina_tree_score = 1
        for i in range(len(main.vina_tree)):
            if i == 0:
                new_path = os.path.join(path_to_et, main.vina_tree[i])
            else:
                new_path = os.path.join(new_path, main.vina_tree[i])

            if os.path.isdir(new_path):
                vina_tree_score += 1

        if vina_tree_score == len(main.vina_tree):
            configuration["vina_tree"] = True

        adfr_tree_score = 1
        for i in range(len(main.adfr_tree)):
            if i == 0:
                new_path = os.path.join(path_to_et, main.adfr_tree[i])
            else:
                new_path = os.path.join(new_path, main.adfr_tree[i])

            if os.path.isdir(new_path):
                adfr_tree_score += 1

        if adfr_tree_score == len(main.adfr_tree):
            configuration["adfr_tree"] = True

        if main.agfr_tree is not None:
            agfr_tree_score = 1
            for i in range(len(main.agfr_tree)):
                if i == 0:
                    new_path = os.path.join(path_to_et, main.agfr_tree[i])
                else:
                    new_path = os.path.join(new_path, main.agfr_tree[i])

                if os.path.isdir(new_path):
                    agfr_tree_score += 1

            if agfr_tree_score == len(main.agfr_tree):
                configuration["agfr_tree"] = True

        if main.smina_tree is not None:
            smina_tree_score = 1
            for i in range(len(main.smina_tree)):
                if i == 0:
                    new_path = os.path.join(path_to_et, main.smina_tree[i])
                else:
                    new_path = os.path.join(new_path, main.smina_tree[i])

                if os.path.isdir(new_path):
                    smina_tree_score += 1

            if smina_tree_score == len(main.smina_tree):
                configuration["smina_tree"] = True

        if main.sdsorter_tree is not None:
            sdsorter_tree_score = 1
            for i in range(len(main.sdsorter_tree)):
                if i == 0:
                    new_path = os.path.join(path_to_et, main.sdsorter_tree[i])
                else:
                    new_path = os.path.join(new_path, main.sdsorter_tree[i])

                if os.path.isdir(new_path):
                    sdsorter_tree_score += 1

            if sdsorter_tree_score == len(main.sdsorter_tree):
                configuration["sdsorter_tree"] = True

        configuration_completed = True

        for key in configuration:
            if configuration[key] == False:
                configuration_completed = False

        if not configuration_completed:
            if os.path.isdir(path_to_adt):
                shutil.rmtree(path_to_adt)
            if os.path.isdir(path_to_et):
                shutil.rmtree(path_to_et)

            if os.path.isfile(path_to_pfr):
                os.remove(path_to_pfr)
            if os.path.isfile(path_to_pl):
                os.remove(path_to_pl)
            if os.path.isfile(path_to_pr):
                os.remove(path_to_pr)

        return configuration, configuration_completed


    def set_pymol_visualization_options(self):

        pymol_version = float(".".join(cmd.get_version()[0].split(".")[0:2]))

        if pymol_version >= 2.5:
            try:
                cmd.undo_disable()
            except:
                pass

        cmd.set("cartoon_fancy_helices", 1)
        cmd.set("cartoon_highlight_color", "grey50")
        cmd.set("sphere_scale", 0.15)
        cmd.set("orthoscopic", "on")
        cmd.set("cartoon_dumbbell_length", 1.5)
        cmd.set("cartoon_dumbbell_width", 0.4)
        cmd.set("cartoon_dumbbell_radius", 0.3)
        cmd.set("cartoon_ring_mode", 1)
        cmd.set("cartoon_ring_transparency", 0.5)


    def initUI(self):

        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        # Loads and sets the Qt stylesheet.
        module_path = sys.modules[__name__].__file__
        self.qss = QSSHelper.open_qss(os.path.join(os.path.dirname(module_path), 'aqua', 'aqua.qss'))
        self.setStyleSheet(self.qss)
        self.show()


class QSSHelper:
    def __init__(self):
        pass

    @staticmethod
    def open_qss(path):
        """
        opens a Qt stylesheet with a path relative to the project

        Note: it changes the urls in the Qt stylesheet (in memory), and makes these urls relative to the project
        Warning: the urls in the Qt stylesheet should have the forward slash ('/') as the pathname separator
        """
        with open(path) as f:
            qss = f.read()
            pattern = r'url\((.*?)\);'
            for url in sorted(set(re.findall(pattern, qss)), key=len, reverse=True):
                directory, basename = os.path.split(path)
                new_url = os.path.join(directory, *url.split('/'))
                new_url = os.path.normpath(new_url)
                new_url = new_url.replace(os.path.sep, '/')
                qss = qss.replace(url, new_url)
            return qss



class Centralwid(QtWidgets.QWidget):

    """
    A class to reppresent the Central Widget of the Main Window - Work in Progress
    """

    def __init__(self, main_window):
        super(Centralwid, self).__init__(main_window)
        self.style = "background-color: rgb(0, 0, 0); color: rgb(255, 255, 255); font-weight: bold"
        self.main_window = main_window
        self.initUI()

    def initUI(self):

        # Work in Progress

        self.central_layout = QtWidgets.QFormLayout()

        self.reinitialize_everything = QtWidgets.QPushButton("Reinitialize Everything")
        self.reinitialize_everything.hide()
        self.prev_btn = QtWidgets.QPushButton("Vina")
        self.next_btn = QtWidgets.QPushButton("RxDock")

        self.bottom_widg = QtWidgets.QWidget()
        self.bottom_widg_layout = QtWidgets.QHBoxLayout()
        self.bottom_widg.setLayout(self.bottom_widg_layout)

        # self.bottom_widg_layout.addWidget(self.prev_btn)
        self.bottom_widg_layout.addStretch()
        self.bottom_widg_layout.addWidget(self.reinitialize_everything)
        # self.bottom_widg_layout.addStretch()
        # self.bottom_widg_layout.addWidget(self.next_btn)
        # self.next_btn.clicked.connect(self.next_tab_func)
        # self.prev_btn.clicked.connect(self.prev_tab_func)


    def next_tab_func(self):
        cur_index = self.main_docking_programs_tabs.currentIndex()
        if cur_index < len(self.main_docking_programs_tabs)-1:
            self.main_docking_programs_tabs.setCurrentIndex(cur_index+1)

    def prev_tab_func(self):
        cur_index = self.main_docking_programs_tabs.currentIndex()
        if cur_index > 0:
            self.main_docking_programs_tabs.setCurrentIndex(cur_index-1)



class Docking_Programs(QtWidgets.QWidget):

    """
    A class to reppresent the Tab Widget of the Docking Programs.
    """

    def __init__(self, main_window):
        super().__init__(main_window)
        self.main_window = main_window

        self.create_child_tabs_for_docking_programs()

    def create_child_tabs_for_docking_programs(self):

        # Customize PyMOl visualization.
        self.main_window.set_pymol_visualization_options()

                ### Inizialize paths to directories ###

        # Path where the plugin is located
        self.current_path = (os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

        # Path to tmp
        self.tmp_dir_path = (os.path.join(self.current_path, "tmp"))

        # Path to general_tmp
        self.general_tmp_dir = os.path.join(self.tmp_dir_path, "General_tmp")

        # Path to consensus_tmp
        self.consensus_tmp_dir = os.path.join(self.tmp_dir_path, "consensus_tmp")

        # Path to RxDock_tmp
        self.rxdock_tmp_dir = os.path.join(self.tmp_dir_path, "RxDock_tmp")

        # Path to Vina_tmp
        self.vina_tmp_dir = os.path.join(self.tmp_dir_path, "Vina_tmp")

        # Path to Smina_tmp
        self.smina_tmp_dir = os.path.join(self.tmp_dir_path, "Smina_tmp")

        # Path to ADFR_tmp
        self.adfr_tmp_dir = os.path.join(self.tmp_dir_path, "ADFR_tmp")

        # Path to Docking_results
        self.docking_results = (os.path.join(self.current_path, "Docking_results"))

        # Path to config
        self.config_path = (os.path.join(self.current_path, "config"))

        # Paths
        if sys.platform == "win32":
            self.et = "external_tools_windows"
            self.vina_tree = ["vina_win32", "bin", "vina.exe"]
            self.adfr_tree = ["adfr_win32"]
            self.agfr_tree = None
            self.smina_tree = None
            self.sdsorter_tree = None

            self.path_sep = "\\"
            self.path_to_vina = (os.path.join(self.config_path, "external_tools_windows", "vina_win32", "bin", "vina.exe"))
            self.path_to_ADFR = (os.path.join(self.config_path, "external_tools_windows", "adfr_win32"))
            # self.path_to_ADFR = (os.path.join(self.config_path, "external_tools_windows", "ADFRsuite_win", "bin", "adfr.bat"))
            # self.path_to_agfr = (os.path.join(self.config_path, "external_tools_windows", "ADFRsuite_win", "bin", "agfr.bat"))

        elif sys.platform == "linux":
            self.et = "external_tools_linux"
            self.vina_tree = ["vina_linux", "bin", "vina"]
            self.adfr_tree = ["ADFRsuite_x86_64Linux_1.0", "bin", "adfr"]
            self.agfr_tree = ["ADFRsuite_x86_64Linux_1.0", "bin", "agfr"]
            self.smina_tree = ["smina_linux", "bin", "smina.static"]
            self.sdsorter_tree = ["sdsorter_linux", "bin", "sdsorter.static"]

            self.path_sep = "/"
            self.path_to_vina = (os.path.join(self.config_path, "external_tools_linux", "vina_linux", "bin", "vina"))
            self.path_to_ADFR = (os.path.join(self.config_path, "external_tools_linux", "ADFRsuite_x86_64Linux_1.0", "bin", "adfr"))
            self.path_to_agfr = (os.path.join(self.config_path, "external_tools_linux", "ADFRsuite_x86_64Linux_1.0", "bin", "agfr"))
            self.path_to_smina = (os.path.join(self.config_path, "external_tools_linux", "smina_linux", "bin", "smina.static"))
            self.path_to_sdsorter = (os.path.join(self.config_path, "external_tools_linux", "sdsorter_linux", "bin", "sdsorter.static"))

        elif sys.platform == "darwin":
            self.et = "external_tools_macOS"
            self.vina_tree = ["vina_darwin", "bin", "vina"]
            self.adfr_tree = ["ADFRsuite_x86_64Darwin_1.0", "bin", "adfr"]
            self.agfr_tree = ["ADFRsuite_x86_64Darwin_1.0", "bin", "agfr"]
            self.smina_tree = ["smina_darwin", "bin", "smina.osx"]
            self.sdsorter_tree = ["sdsorter_darwin", "bin", "sdsorter.osx"]

            self.path_sep = "/"
            self.path_to_vina = (os.path.join(self.config_path, "external_tools_macOS", "vina_darwin", "bin", "vina"))
            self.path_to_ADFR = (os.path.join(self.config_path, "external_tools_macOS", "ADFRsuite_x86_64Darwin_1.0", "bin", "adfr"))
            self.path_to_agfr = (os.path.join(self.config_path, "external_tools_macOS", "ADFRsuite_x86_64Darwin_1.0", "bin", "agfr"))
            self.path_to_smina = (os.path.join(self.config_path, "external_tools_macOS", "smina_darwin", "bin", "smina.osx"))
            self.path_to_sdsorter = (os.path.join(self.config_path, "external_tools_macOS", "sdsorter_darwin", "bin", "sdsorter.osx"))

        # Path to RxDock config files
        self.path_to_cavity = (os.path.join(self.config_path, "RxDock", "Cavity"))
        self.path_to_pharma = (os.path.join(self.config_path, "RxDock", "pharma.const"))
        self.path_to_prm_file = (os.path.join(self.config_path, "RxDock", "standard_prm_file"))

        # Dictionaries to store informations about the ligands and structures that have been loaded.
        self.vina_receptors_dict = {}
        self.vina_ligands_dict = {}

        self.rxdock_receptors_dict = {}
        self.rxdock_ligands_dict = {}

        self.smina_receptors_dict = {}
        self.smina_ligands_dict = {}

        self.adfr_receptors_dict = {}
        self.adfr_ligands_dict = {}

        # Dictionary to keep track of the created grids {name: (x_center, y_center, x_center, spacing, x_spacing, y_spacing, z_spacing))}
        self.grid_center = {}
        self.ready_grid_centers = {}

        # Lists to keep track of the names and amount of loaded and prepared structures and ligands
        self.vina_receptors_prepared = []
        self.vina_ligands_prepared = []

        self.rxdock_receptors_prepared = []
        self.rxdock_ligands_prepared = []

        self.smina_receptors_prepared = []
        self.smina_ligands_prepared = []

        self.adfr_receptors_prepared = []
        self.adfr_ligands_prepared = []

        # Lists to keep track of the names and amount of combined objects created in PyMOl
        self.combined_objs_list = []

        # List of generated cavities with rxdock
        self.generated_cavity = []

        # Dict to keep track of the docking processes that have been run
        self.rxdock_runs_dict = {}
        self.smina_runs_dict = {}
        self.vina_runs_dict = {}
        self.adfr_runs_dict = {}

        self.rxdock_runs = 1
        self.smina_runs = 1
        self.vina_runs = 1
        self.adfr_runs = 1

        self.all_programs_dict = [self.rxdock_runs_dict, self.smina_runs_dict, self.vina_runs_dict, self.adfr_runs_dict]

        self.all_runs = {}

        # Dict to keep track of the results - Work in Progress
        self.results_dict = {}

        # Create the Docking Programs' Tab Widget and set the Layout
        self.docking_prog_tabs_widg = QtWidgets.QWidget()
        self.docking_programs_tabs = QtWidgets.QTabWidget()

        self.docking_programs_layout = QtWidgets.QFormLayout()

        # Create configuration tab
        self.CONFIGURATION = QtWidgets.QWidget()
        self.docking_programs_tabs.addTab(self.CONFIGURATION, "CONFIGURATION")
        self.CONFIGURATION_TAB = ConfigurationTab(self)
        self.CONFIGURATION.setLayout(self.CONFIGURATION_TAB.layout_config_tab)


                ### Create Child Tabs for each Docking Program ###

                # (*1) note: RxDock or Vina specification in Child_Tabs is needed to set some peculiarity of a single Docking Program in DockingTab layout

        self.RXDOCK = Child_Tabs(self, tab = "RxDock")
        self.docking_programs_tabs.addTab(self.RXDOCK.child_tabs_widget, "RxDock")

        self.VINA = Child_Tabs(self, tab = "Vina")
        self.docking_programs_tabs.addTab(self.VINA.child_tabs_widget, "Vina")

        self.SMINA = Child_Tabs(self, tab = "Smina")
        self.docking_programs_tabs.addTab(self.SMINA.child_tabs_widget, "Smina")

        self.ADFR = Child_Tabs(self, tab = "ADFR")
        self.docking_programs_tabs.addTab(self.ADFR.child_tabs_widget, "ADFR")

        if sys.platform == "win32":
            self.docking_programs_tabs.setTabEnabled(1, False)
            self.docking_programs_tabs.setTabEnabled(3, False)


                ### Create Data Analyses Tab ###

        self.DATA_ANALYSIS = QtWidgets.QWidget()
        self.docking_programs_tabs.addTab(self.DATA_ANALYSIS, "DATA ANALYSIS")
        self.data_analysis_layout = DataAnalysisTab(self)
        self.DATA_ANALYSIS.setLayout(self.data_analysis_layout.layout_data_analysis_tab)


                ### Set the Layouts for GridTab ###

        # Different Docking Programs use the same Layouts for the ReceptorTab, the LigandsTab, the DockingTab(*1) and the ResultsTab
        # The GridTab is peculiar for each Docking Program

        self.RXDOCK.grid_settings.setLayout(self.RXDOCK.grid_tab_rxdock_ui.layout_grid_tab_rxdock)
        self.SMINA.grid_settings.setLayout(self.SMINA.grid_tab_smina_ui.layout_grid_tab)
        self.VINA.grid_settings.setLayout(self.VINA.grid_tab_vina_ui.layout_grid_tab)
        self.ADFR.grid_settings.setLayout(self.ADFR.grid_tab_adfr_ui.layout_grid_tab)

        # Rxdock additional Tabs
        # self.htvs = QtWidgets.QWidget()
        # self.RXDOCK.child_tabs.addTab(self.htvs, "HTVS")
        # self.RXDOCK.child_tabs_layout.addWidget(self.RXDOCK.child_tabs)

        # Gnina and Smina additional tabs
        # self.scoring = QtWidgets.QWidget()
        # self.SMINA.child_tabs.addTab(self.scoring, "Scoring")
        # self.SMINA.child_tabs_layout.addWidget(self.SMINA.child_tabs)
        #
        # self.scoring = QtWidgets.QWidget()

        # Add the Nested Tab Widget to the Layout
        self.docking_programs_layout.addWidget(self.docking_programs_tabs)

        return self.docking_programs_tabs


class Child_Tabs(QtWidgets.QWidget):

    """
    A class to reppresent the Child Tabs of each Docking Program.
    """

    def __init__(self, docking_programs, tab):
        super().__init__(docking_programs)
        self.docking_programs = docking_programs

        self.create_child_tabs(tab)

    def create_child_tabs(self, tab, *additional_tabs):

        # Create the Tab Widget and set the Layout
        self.child_tabs_widget = QtWidgets.QWidget()
        self.child_tabs = QtWidgets.QTabWidget()

        self.child_tabs_layout = QtWidgets.QHBoxLayout()
        self.child_tabs_widget.setLayout(self.child_tabs_layout)

        # ReceptorTab
        self.receptor = QtWidgets.QWidget()
        self.child_tabs.addTab(self.receptor, "Receptors")
        self.receptor.setLayout(ReceptorTab(self, current_tab = tab).layout_receptor_tab)

        # LigandTab
        self.ligands = QtWidgets.QWidget()
        self.child_tabs.addTab(self.ligands, "Ligands")
        self.ligands.setLayout(LigandTab(self, current_tab = tab).layout_ligand_tab)

        # GridTab
        self.grid_settings = QtWidgets.QWidget()
        self.child_tabs.addTab(self.grid_settings, "Grid settings")

        # GridTab layouts
        self.grid_tab_rxdock_ui = GridTab_RxDock(self)
        self.grid_tab_vina_ui = GridTab(self, current_tab = "Vina")
        self.grid_tab_smina_ui = GridTab(self, current_tab = "Smina")
        self.grid_tab_adfr_ui = GridTab(self, current_tab = "ADFR")

        # DockingTab
        self.docking = QtWidgets.QWidget()
        self.child_tabs.addTab(self.docking, "Docking")
        self.docking_tab_ui = DockingTab(self, current_tab = tab)
        self.docking.setLayout(self.docking_tab_ui.layout_docking_tab)

        # ResultsTab
        self.results = QtWidgets.QWidget()
        self.child_tabs.addTab(self.results, "Results")
        self.results_tab_ui = ResultsTab(self)
        self.results.setLayout(self.results_tab_ui.layout_tab6)

        self.child_tabs_layout.addWidget(self.child_tabs)

        if additional_tabs is not None:

            for tab in additional_tabs:
                self.tab_name = tab
                self.tab_name = QtWidgets.QWidget()
                self.child_tabs.addTab(self.tab_name, str(tab))
                self.child_tabs_layout.addWidget(self.child_tabs)

        return self.child_tabs_widget
