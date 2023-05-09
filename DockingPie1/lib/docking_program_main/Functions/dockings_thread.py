# Copyright 2022 by Serena Rosignoli. All rights reserved.
# This code is part of the DockingPie package and governed by its license.
# Please see the LICENSE file.


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


##
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
from lib.docking_program_main.docking_program_gui.frames import *



def catch_errors_installer_threads(function):
    """
    Catches errors in the installer threads and closes the dialogs.
    """

    def wrapper(self, *args, **kwargs):
        try:
            return function(self, *args, **kwargs)
        except Exception as e:
            self.emit_exception.emit(e)

    return wrapper



class _dialog_mixin:
    """
    Mixin class to be incorporated in all the installation process dialogs.
    """

    def keyPressEvent(self, event):
        """
        By overriding this method, the dialog will not close when pressing the "esc" key.
        """
        if event.key() == QtCore.Qt.Key_Escape:
            pass
        else:
            QtWidgets.QDialog.keyPressEvent(self, event)


class DialogThread(_dialog_mixin, QtWidgets.QDialog):

    """
    A dialog from within to initialize Threads
    """

    def __init__(self, tab, main):

        self.tab = tab
        self.main = main

        QtWidgets.QDialog.__init__(self, parent=self.tab)


    def setup_thread(self, starting_function):

        self.jobthread = JobThread(self)

        # Setup event-signal actions
        self.jobthread.set_params(self, self.tab, self.main, starting_function)


    def on_button_click(self):

        text = self.start_button.text()

        if text == "Close":
            self.on_closing_thread()
        else:
            self.jobthread.start()
            self.setWindowTitle('Docking in process. Please wait ...')
            #self.docking_progressbar.setFormat("Running...")
            #self.initial_label.setText("")
            self.start_button.setEnabled(False)
            self.close_button.hide()


    def close_dialog_prior_process(self):

        self._terminate_threads()
        self.close()

    def _terminate_threads(self):
        if self.jobthread.isRunning():
            self.jobthread.terminate()

    def closeEvent(self, evnt):
        self._terminate_threads()
        self.close()


class JobThread(QtCore.QThread):

    """
    A Thread for Docking Processes
    """

    def __init__(self, tab):
        super().__init__(tab)

    # Signals.
    emit_exception = QtCore.pyqtSignal(Exception)
    update_single_docking_progressbar = QtCore.pyqtSignal(int)
    update_consensus_docking_progressbar = QtCore.pyqtSignal(int)
    process_completed = QtCore.pyqtSignal(object)
    update_progress_text = QtCore.pyqtSignal(str)
    user_asks_to_interrupt = QtCore.pyqtSignal()
    interrupt = QtCore.pyqtSignal(object)


    def set_params(self, dialog, tab, main, starting_function):

        self.dialog = dialog
        self.tab = tab
        self.main = main

        # The function to run when starting thread
        self.starting_function = starting_function


    @catch_errors_installer_threads
    def run(self):

        self.starting_function()


class Dockings_dialog(_dialog_mixin, QtWidgets.QDialog):

    """
    A dialog from within to initialize Thread for Docking Processes
    """

    def __init__(self, tab, number_of_dockings_to_do, ligands_to_dock, receptors_to_dock, dockings_to_do = None):

        # To assign to DockingTab object
        self.tab = tab

        QtWidgets.QDialog.__init__(self, parent=self.tab)

        self.complete_status = False
        self.maximum_prog_bar = number_of_dockings_to_do

        ####
        # Initialize UI of the Dockings Dialog (Window, progressbar etc...)
        ####
        self.initUI()

        self.docking_thread = Dockings_thread(self)

        self.docking_thread.update_progressbar.connect(self.on_update_progressbar)
        self.docking_thread.update_single_docking_progressbar.connect(self.on_update_single_docking_progressbar)
        self.docking_thread.set_params(self.tab, self.tab.docking_programs_child_tabs.docking_programs,
        number_of_dockings_to_do = number_of_dockings_to_do,
        ligands_to_dock = ligands_to_dock,
        receptors_to_dock = receptors_to_dock,
        dockings_to_do = dockings_to_do)
        self.docking_thread.emit_exception.connect(self.on_emit_exception)
        self.docking_thread.update_results_tab.connect(self.on_update_results_tab)
        self.docking_thread.check_docking_completed.connect(self.on_checking_docking_completed)
        self.docking_thread.update_progress_text.connect(self.on_update_progress_text)
        self.docking_thread.process_completed.connect(self.on_process_completed)
        self.docking_thread.user_asks_to_interrupt.connect(self.on_interrupting_process)
        self.docking_thread.update_data_analysis_tab.connect(self.on_update_data_analysis_tab)

        ### When multiple dockings must be run:
        # this variable is used to stop prior the starting of the next process if the user clicks on 'Cancel' button
        self.tab.user_asks_to_interrupt = False


    def initUI(self):

        Check_current_tab.check_docking_program_current_tab(self.tab)

        self.setWindowTitle('Initializing Docking Process ...')
        self.setWindowFlags(self.windowFlags() | QtCore.Qt.CustomizeWindowHint)
        self.setWindowFlags(self.windowFlags() & ~QtCore.Qt.WindowMaximizeButtonHint)
        self.setWindowFlags(self.windowFlags() & ~QtCore.Qt.WindowCloseButtonHint)

        vertical_layout = QtWidgets.QVBoxLayout()
        self.docking_progressbar = QtWidgets.QProgressBar(self)
        # self.progress.setGeometry(0, 0, 340, 25)
        self.docking_progressbar.setMaximum(self.maximum_prog_bar)

        progressbar_text = "Click the button to start running the Docking Processes"

        self.docking_progressbar.setValue(0)
        vertical_layout.addWidget(self.docking_progressbar)
        self.progressbar_label = QtWidgets.QLabel(progressbar_text)
        vertical_layout.addWidget(self.progressbar_label)

        ## VINA DOCKING SETUP
        if self.tab.is_rxdock_tab:
            pass
        else:
            self.single_docking_progressbar = QtWidgets.QProgressBar(self)
            self.single_docking_progressbar.setMaximum(51)
            self.single_docking_progressbar.setValue(0)
            vertical_layout.addWidget(self.single_docking_progressbar)

        # Button for starting the installation.
        horizontal_layout = QtWidgets.QHBoxLayout()
        start_button_text = 'Start'

        self.start_button = QtWidgets.QPushButton(start_button_text, self)
        # self.start_button.setStyleSheet(label_style_2)
        self.start_button.clicked.connect(self.on_button_click)
        horizontal_layout.addWidget(self.start_button)

        horizontal_layout.addStretch(1)

        # Button for closing dialog before the process starts.
        self.close_button = QtWidgets.QPushButton('Close', self)
        # self.cancel_button.setStyleSheet(label_style_2)
        self.close_button.clicked.connect(self.close_dialog_prior_process)
        horizontal_layout.addWidget(self.close_button)

        horizontal_layout.addStretch(1)

        # Button for canceling the installation.
        self.cancel_button = QtWidgets.QPushButton('Interrupt', self)
        # self.cancel_button.setStyleSheet(label_style_2)
        self.cancel_button.clicked.connect(self.on_cancel_button_click)
        self.cancel_button.setEnabled(False)
        horizontal_layout.addWidget(self.cancel_button)

        vertical_layout.addLayout(horizontal_layout)

        self.setLayout(vertical_layout)


    def close_dialog_prior_process(self):

        self._terminate_threads()
        self.close()

    def on_process_completed(self):
        self.setWindowTitle('Done')
        self.start_button.setText("Close")
        self.start_button.setEnabled(True)
        self.complete_status = True

    def on_update_progress_text(self, text):

        self.progressbar_label.setText(text)


    def on_checking_docking_completed(self, runs, tmp_dir):

        """
        To check if the docking process was completed or interrupted
        """

        ## To think about: the warning message was implemented when Docking Thread was not. Is the warning message still necessary?

        # Path to results file
        # file_path = os.path.join(tmp_dir, self.tab.last_docking.results_file_name_ext)
        #
        # self.tab.docking_completed = True
        #
        # # Check first if the file exists
        # if Path(file_path).is_file():
        #
        #     # If the file exists, check if it contains data.
        #     if os.path.getsize(file_path):
        #         self.tab.docking_completed = True
        #
        #     else:
        #         # If the file exists but it is empty, the docking was not completed
        #         self.tab.docking_completed = False
        #
        #         # If the file exists but it is empty, if the docking runs are more than one, don't show the Warning Message
        #         if len(self.docking_thread.ligands_to_dock) > 1 or len(self.docking_thread.receptors_to_dock) > 1:
        #             os.remove(file_path)
        #         elif len(self.docking_thread.dockings_to_do) > 1:
        #             os.remove(file_path)
        #         else:
        #             os.remove(file_path)
        #             QtWidgets.QMessageBox.warning(self.tab.docking_programs_child_tabs, "", str("Something went wrong during Docking. \nPlease check LOG files."))
        #
        # else:
        #     self.tab.docking_completed = False
        #     QtWidgets.QMessageBox.warning(self.tab.docking_programs_child_tabs, "", str("Something went wrong during Docking. \nPlease check LOG files."))

        file_path = os.path.join(tmp_dir, self.tab.last_docking.results_file_name_ext)

        self.tab.docking_completed = True

        # Check first if the file exists
        if Path(file_path).is_file():

            # If the file exists, check if it contains data.
            if os.path.getsize(file_path):
                self.tab.docking_completed = True

            else:
                # If the file exists but it is empty, the docking was not completed
                self.tab.docking_completed = False

                os.remove(file_path)

        else:
            self.tab.docking_completed = False


    def on_update_data_analysis_tab(self, results_file_name, all_runs_dict, last_docking, results_file):

        # Create checkbox in DATA ANALYSIS tab
        self.tab.dockings_frame = DockingsFrame(parent=self.tab, main_window=self.tab.docking_programs_child_tabs, results_name=results_file_name)
        self.tab.docking_programs_child_tabs.docking_programs.data_analysis_layout.docking_runs_scroll_layout.addRow(self.tab.dockings_frame)

        # Update dict with info about Docking runs
        all_runs_dict[self.tab.last_docking.results_file_name] = {}
        all_runs_dict[self.tab.last_docking.results_file_name]["docking_run"] = last_docking
        all_runs_dict[self.tab.last_docking.results_file_name]["docking_run_results_file"] = results_file
        all_runs_dict[self.tab.last_docking.results_file_name]["dockings_frame"] = self.tab.dockings_frame

        self.cancel_button.setEnabled(False)


    def on_update_results_tab(self, results_tab, docking_process, program, columns_names):

        """
        When the Docking process is completed, a signal is emitted to update the UI of the 'Results' Tab and 'DATA ANALYSIS' Tab
        n.b.: signal is needed when updating a UI from a thread
        """

        # The ResultsTab is updated both in case of successful and unsuccessful Docking
        docking_programs = self.tab.docking_programs_child_tabs.docking_programs

        if program == "Vina":
            name = "Vina"
            main_tab = docking_programs.VINA
            dict = docking_programs.vina_runs_dict

        if program == "RxDock":
            name = "RxDock"
            main_tab = docking_programs.RXDOCK
            dict = docking_programs.rxdock_runs_dict

        if program == "Smina":
            name = "Smina"
            main_tab = docking_programs.SMINA
            dict = docking_programs.smina_runs_dict

        if program == "ADFR":
            name = "ADFR "
            main_tab = docking_programs.ADFR
            dict = docking_programs.adfr_runs_dict

        if self.tab.docking_completed:

            results_tab.tab_widget = QtWidgets.QWidget()

            results_tab.tab_layout = QtWidgets.QGridLayout()
            results_tab.tab_widget.setLayout(results_tab.tab_layout)
            results_tab.result_tabs.addTab(results_tab.tab_widget, docking_process.results_file_name)
            results_tab.results_group_layout.addWidget(results_tab.result_tabs, 0, 0)

            self.tab.results_obj_name = docking_process.results_file_name
            self.tab.results_file_name = docking_process.results_file_name_ext
            self.tab.log_file_name = docking_process.log_file_name

            if type(self.tab.results_file_name) is list:
                for i in self.tab.results_file_name:
                    self.tab.docking_programs_child_tabs.docking_programs.results_dict[i] = {}

            else:
                self.tab.docking_programs_child_tabs.docking_programs.results_dict[self.tab.results_file_name] = {}

            self.tab.results_frame = ResultsFrame(parent=self.tab,
            main_window=self.tab.docking_programs_child_tabs,
            results_file = self.tab.results_file,
            results_file_name=self.tab.results_file_name,
            results_obj_name=self.tab.results_obj_name,
            results_data = self.tab.results_file.results_data,
            log_file = self.tab.log_file_name,
            program = program,
            columns_names = columns_names)

            results_tab.tab_layout.addWidget(self.tab.results_frame)

            #Load the results in PyMOL
            if type(self.tab.results_file_name) is list:
                for i in self.tab.results_file_name:
                    cmd.load(i, self.tab.results_obj_name)

            else:
                cmd.load(self.tab.results_file_name, self.tab.results_obj_name)

            # Group the results in PyMOL
            cmd.group(name, members=self.tab.results_obj_name, action='auto')

        else:

            results_tab.tab_widget = QtWidgets.QWidget()

            results_tab.tab_layout = QtWidgets.QGridLayout()
            results_tab.tab_widget.setLayout(results_tab.tab_layout)
            results_tab.result_tabs.addTab(results_tab.tab_widget, docking_process.results_file_name)
            results_tab.results_group_layout.addWidget(results_tab.result_tabs, 0, 0)
            results_tab.save_all_files.setEnabled(True)
            results_tab.save_all_files.show()

            self.tab.results_obj_name = docking_process.results_file_name
            self.tab.results_file_name = docking_process.results_file_name_ext
            self.tab.log_file_name = docking_process.log_file_name

            self.tab.docking_programs_child_tabs.docking_programs.results_dict[self.tab.results_file_name] = {}

            self.tab.results_frame = ResultsFrame(parent=self.tab,
            main_window=self.tab.docking_programs_child_tabs,
            results_file = None,
            results_file_name=self.tab.results_file_name,
            results_obj_name=self.tab.results_obj_name,
            log_file = self.tab.log_file_name,
            results_data = None,
            program = program,
            columns_names = columns_names)

            results_tab.tab_layout.addWidget(self.tab.results_frame)

        main_tab.results_tab_ui.save_all_files.setEnabled(True)
        main_tab.results_tab_ui.all_plot_btn.setEnabled(True)
        main_tab.results_tab_ui.save_all_files.show()
        main_tab.results_tab_ui.all_plot_btn.show()

        # Update dict with info about Docking runs
        dict[self.tab.last_docking.results_file_name] = {}
        dict[self.tab.last_docking.results_file_name]["docking_run"] = self.tab.last_docking
        dict[self.tab.last_docking.results_file_name]["program"] = name
        dict[self.tab.last_docking.results_file_name]["results_table"] = self.tab.results_frame.results_table

        return self.tab.results_frame


    # Interactions with the buttons.
    def on_button_click(self):

        text = self.start_button.text()

        if text == "Close":
            self.on_closing_thread()
        else:
            self.docking_thread.start()
            self.setWindowTitle('Running Job. Please wait ...')
            #self.docking_progressbar.setFormat("Running...")
            self.progressbar_label.setText("")
            self.start_button.setEnabled(False)
            self.close_button.hide()

    def on_closing_thread(self):

        # When the 'Close' button is enabled
        #self._terminate_threads()
        self.close()

    def on_interrupting_process(self):
        self.setWindowTitle('Interrupted')
        self.start_button.setText("Close")
        self.start_button.setEnabled(True)
        self.complete_status = True

    def pid_exists(self, pid):
        if pid < 0: return False #NOTE: pid == 0 returns True

        if sys.platform == "win32":
            try:
                os.kill(pid, signal.CTRL_C_EVENT)
            except ProcessLookupError: # errno.ESRCH
                return False # No such process
            except PermissionError: # errno.EPERM
                return True # Operation not permitted (i.e., process exists)
            else:
                return True # no error, we can send a signal to the process

        else:
            try:
                os.kill(pid, 0)
            except ProcessLookupError: # errno.ESRCH
                return False # No such process
            except PermissionError: # errno.EPERM
                return True # Operation not permitted (i.e., process exists)
            else:
                return True # no error, we can send a signal to the process


    def on_cancel_button_click(self):

        if sys.platform == "win32":
            qm = QtWidgets.QMessageBox.question(self.tab.docking_programs_child_tabs,'Interrupting Docking Process', str('Are your really sure you want to interrupt the Docking Process?'), QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

            if qm == QtWidgets.QMessageBox.Yes:

                os.system("TASKKILL /F /PID " + str(self.docking_thread.docking_subprocess.pid))
                self.tab.user_asks_to_interrupt = True
                self.cancel_button.setEnabled(False)

            elif qm == QtWidgets.QMessageBox.No:
                pass

        else:
            # process_exists = self.pid_exists(self.docking_thread.docking_subprocess.pid)
            #
            # if process_exists:

            qm = QtWidgets.QMessageBox.question(self.tab.docking_programs_child_tabs,'Interrupting Docking Process', str('Are your really sure you want to interrupt the Docking Process?'), QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

            if qm == QtWidgets.QMessageBox.Yes:

                try:
                    self.tab.user_asks_to_interrupt = True
                    self.cancel_button.setEnabled(False)
                    os.killpg(os.getpgid(self.docking_thread.docking_subprocess.pid), signal.SIGHUP)
                except:
                    pass

            elif qm == QtWidgets.QMessageBox.No:
                pass

            # else:
            #     self.on_closing_thread()


    def closeEvent(self, evnt):
        self._terminate_threads()


    def _terminate_threads(self):
        if self.docking_thread.isRunning():
            self.docking_thread.terminate()


    def on_update_single_docking_progressbar(self, value):
        self.single_docking_progressbar.setValue(value)


    def on_update_progressbar(self, value):

        """
        Updates the percentage in the progressbar.
        """
        self.docking_progressbar.setValue(value)


    def on_emit_exception(self, e):
        """
        Quit the threads and close the installation dialog.
        """
        self._terminate_threads()

        message = "There was an error: %s" % str(e)
        if hasattr(e, "url"):
            message += " (%s)." % e.url
        else:
            message += "."
        message += " Quitting the Docking process."
        print(message)
        ##self.tab.show_error_message("Installation Error", message)

        self.close()



class Dockings_thread(QtCore.QThread):

    """
    A Thread for Docking Processes
    """

    def __init__(self, tab):
        super().__init__(tab)

    # Signals.
    update_progressbar = QtCore.pyqtSignal(int)
    emit_exception = QtCore.pyqtSignal(Exception)
    update_single_docking_progressbar = QtCore.pyqtSignal(int)

    update_progress_text = QtCore.pyqtSignal(str)
    update_results_tab = QtCore.pyqtSignal(object, object, str, list)
    check_docking_completed = QtCore.pyqtSignal(int, str)
    docking_process = QtCore.pyqtSignal(object)
    process_completed = QtCore.pyqtSignal()
    user_asks_to_interrupt = QtCore.pyqtSignal()
    update_data_analysis_tab = QtCore.pyqtSignal(str, dict, object, object)


    def set_params(self, tab, main, number_of_dockings_to_do, ligands_to_dock, receptors_to_dock, dockings_to_do = []):

        self.tab = tab
        self.main = main

        self.number_of_dockings_to_do = number_of_dockings_to_do

        self.ligands_to_dock = ligands_to_dock
        self.receptors_to_dock = receptors_to_dock

        self.dockings_to_do = dockings_to_do


    @catch_errors_installer_threads
    def run(self):

        if self.tab.protocol_box.currentText() == "Custom":

            i = 1
            self.update_progress_text.emit("Running docking " + str(int(i)) + " of " + str(len(self.dockings_to_do)))

            for docks in self.dockings_to_do:
                if not self.tab.user_asks_to_interrupt:

                    receptor = docks["receptor"]
                    cavity = docks["cavity"]
                    ligand = docks["ligand"]

                    self.starting_docking_process(cavity = cavity, receptor = receptor, ligand = ligand, num = i)
                    self.update_progressbar.emit(int(i))

                    i += 1
                    if i > len(self.dockings_to_do):
                        self.update_progress_text.emit("Process Completed")
                    elif i <= len(self.dockings_to_do):
                        self.update_progress_text.emit("Running docking  " + str(int(i)) + " of " + str(len(self.dockings_to_do)))

                else:
                    self.update_progress_text.emit("Process Interrupted")

            # Completes and updates the GUI.
            if self.tab.user_asks_to_interrupt == False:
                self.process_completed.emit()
                self.update_progressbar.emit(len(self.dockings_to_do))
            else:
                self.user_asks_to_interrupt.emit()

            # Wait a little bit of time.
            time.sleep(0.5)

        if self.tab.protocol_box.currentText() == "ALLvsALL":

            i = 1
            self.update_progress_text.emit("Running docking  " + str(int(i)) + " of " + str(self.number_of_dockings_to_do))

            for receptor in self.receptors_to_dock:
                for ligand in self.ligands_to_dock:
                    if not self.tab.user_asks_to_interrupt:

                        self.starting_docking_process(receptor = receptor, ligand = ligand, num = i)
                        self.update_progressbar.emit(int(i))

                        i += 1

                        if i > int(self.number_of_dockings_to_do):
                            self.update_progress_text.emit("Process Completed")
                            self.update_progressbar.emit(int(self.number_of_dockings_to_do))
                        elif i <= int(self.number_of_dockings_to_do):
                            self.update_progress_text.emit("Running docking  " + str(int(i)) + " of " + str(int(self.number_of_dockings_to_do)))

                    else:
                        self.update_progress_text.emit("Process Interrupted")

            # Completes and updates the GUI.
            if self.tab.user_asks_to_interrupt == False:
                self.process_completed.emit()
            else:
                self.user_asks_to_interrupt.emit()
                self.update_progress_text.emit("Process Interrupted")

            # Wait a little bit of time.
            time.sleep(0.5)


    def starting_docking_process(self, receptor, ligand, cavity = None, num = 0):

        Check_current_tab.check_docking_program_current_tab(self.tab)

        docking_programs = self.main

        self.tab.docking_completed = False

        ## VINA DOCKING SETUP
        if self.tab.is_vina_tab:

            if cavity is not None:

                self.tab.last_docking = Vina_docking(self.tab,
                ligand = ligand,
                receptor = receptor,
                cavity = cavity,
                tmp_dir = self.main.vina_tmp_dir,
                cavity_list = self.main.ready_grid_centers[cavity],
                main = self.main,
                ligands_to_dock = self.ligands_to_dock,
                receptors_to_dock = self.receptors_to_dock,
                use_flex_vina = self.tab.use_flex_vina_cb.isChecked(),
                flex_residues = self.tab.flex_arg_edit.text(),
                poses = str(self.tab.poses_box.value()),
                exhaustiveness = self.tab.exhaustiveness_box.value(),
                energy = self.tab.energy_box.value(),
                scoring_function = self.tab.scoring_box.currentText())

            else:

                self.tab.last_docking = Vina_docking(self.tab,
                ligand = ligand,
                receptor = receptor,
                tmp_dir = self.main.vina_tmp_dir,
                cavity = self.tab.loaded_cavities.currentText(),
                cavity_list = self.main.ready_grid_centers[self.tab.loaded_cavities.currentText()],
                main = self.main,
                ligands_to_dock = self.ligands_to_dock,
                receptors_to_dock = self.receptors_to_dock,
                use_flex_vina = self.tab.use_flex_vina_cb.isChecked(),
                flex_residues = self.tab.flex_arg_edit.text(),
                poses = str(self.tab.poses_box.value()),
                exhaustiveness = self.tab.exhaustiveness_box.value(),
                energy = self.tab.energy_box.value(),
                scoring_function = self.tab.scoring_box.currentText())

            ferr = open('stdout.txt','w+')
            # Run Docking Process in a different environment, to facilitate the interruption of the protocol
            if sys.platform == "win32":

                CREATE_NO_WINDOW = 0x08000000

                self.docking_subprocess = subprocess.Popen(self.tab.last_docking.run_docking_vina_settings,
                creationflags=subprocess.CREATE_NEW_PROCESS_GROUP | CREATE_NO_WINDOW,
                stdout = ferr)
            else:
                self.docking_subprocess = subprocess.Popen(self.tab.last_docking.run_docking_vina_settings,
                preexec_fn=os.setsid,
                stdout=ferr,
                )

            # 'Cancel' button to interrupt the process is enabled after a bit, to wait for effective begin of the process
            time.sleep(2)
            self.tab.docking_dialog.cancel_button.setEnabled(True)

            # To interrupt the process.
            # process.Popen.communicate() would wait until the end of the process before to continue.
            # Here process.Popen.communicate() is not used to ensure for the possibility to interrupt the protocol
            # In any case, to wait until the end of the process before to continue, process.poll() is used

            # In while loop, the file where the stdout is stored, is read to update progressbar
            #
            # For future developers: Use 'time.sleep(5)' if nothing is done in while loop
            # 'time.sleep(5)' instead of 'pass' is used to reduce the number of time the process is checked.

            # for line in ferr:
            #     print(line)

            ferr.flush()

            bar_value = 0
            current_size = os.path.getsize('stdout.txt')
            while self.docking_subprocess.poll() is None:
                time.sleep(0.7)
                writing = open('stdout.txt','r')
                for line in writing:
                    if line.startswith("*"):
                        bar_value = line.count("*")
                        self.update_single_docking_progressbar.emit(int(bar_value))
                writing.close()


            self.update_single_docking_progressbar.emit(int(51))

            ferr.close()

            if os.path.isfile(self.tab.last_docking.log_file_name):
                os.remove('stdout.txt')
            else:
                os.rename('stdout.txt', self.tab.last_docking.log_file_name)


            if self.tab.last_docking.interrupt == False:

                self.check_docking_completed.emit(docking_programs.vina_runs,
                docking_programs.vina_tmp_dir)

                time.sleep(1)

                self.tab.summary_file = self.write_summary_file(receptor = receptor, ligand = ligand)

                if self.tab.docking_completed:

                    self.tab.results_file = Vina_Parse_Results(self.tab,
                    main = self.main,
                    last_docking = self.tab.last_docking,
                    results_file_name = self.tab.last_docking.results_file_name,
                    results_dict = docking_programs.results_dict,
                    poses = self.tab.last_docking.poses,
                    ligand = ligand)

                self.update_results_tab.emit(docking_programs.VINA.results_tab_ui,
                self.tab.last_docking,
                "Vina",
                ["NAME", "POSE", "Affinity\n(kcal/mol)", "dist from best mode\nRMSD l.b.", "dist from best mode\nRMSD u.b."])


        ## RXDOCK DOCKING SETUP
        if self.tab.is_rxdock_tab:

            self.get_rxdock_docking_protocol()

            # Run RxDock Docking

            if cavity is not None:

                self.tab.last_docking = RxDock_docking(self.tab,
                ligand = ligand,
                receptor = receptor,
                cavity = cavity,
                tmp_dir = self.main.rxdock_tmp_dir,
                trans_mode = self.tab.trans_mode_combo.currentText(),
                rot_mode = self.tab.rot_mode_combo.currentText(),
                die_mode = self.tab.die_mode_combo.currentText(),
                pharma_restrains = self.tab.pharma_restrains,
                tethered_docking = self.tab.tethered_docking,
                poses_box = str(self.tab.poses_box.value()),
                cavity_to_dock = self.tab.loaded_cavities.currentText(),
                cavity_name = self.tab.loaded_cavities.currentText(),
                protein_segments_to_exclude = self.tab.protein_segments_to_exclude,
                use_water = self.tab.receptor_water_checkbtn.isChecked(),
                main = self.main)

            else:

                self.tab.last_docking = RxDock_docking(self.tab,
                ligand = ligand,
                receptor = receptor,
                tmp_dir = self.main.rxdock_tmp_dir,
                trans_mode = self.tab.trans_mode_combo.currentText(),
                rot_mode = self.tab.rot_mode_combo.currentText(),
                die_mode = self.tab.die_mode_combo.currentText(),
                pharma_restrains = self.tab.pharma_restrains,
                tethered_docking = self.tab.tethered_docking,
                poses_box = str(self.tab.poses_box.value()),
                cavity_to_dock = self.tab.loaded_cavities.currentText(),
                cavity_name = self.tab.loaded_cavities.currentText(),
                protein_segments_to_exclude = self.tab.protein_segments_to_exclude,
                use_water = self.tab.receptor_water_checkbtn.isChecked(),
                main = self.main)



            # Run Docking Process in a different environment, to facilitate the interruption of the protocol
            ferr = open('stdout.txt','w+')
            self.docking_subprocess = subprocess.Popen(self.tab.last_docking.run_docking_rxdock_settings, preexec_fn=os.setsid, stdout = ferr)
            ferr.close()

            if os.path.isfile(self.tab.last_docking.log_file_name):
                os.remove('stdout.txt')
            else:
                os.rename('stdout.txt', self.tab.last_docking.log_file_name)

            # 'Cancel' button to interrupt the process, is enabled after a bit, to wait for effective begin of the process
            self.tab.docking_dialog.cancel_button.setEnabled(True)

            # To interrupt the process.
            # process.Popen.communicate() would wait until the end of the process before to continue.
            # Here process.Popen.communicate() is not used to ensure for the possibility to interrupt the protocol
            # In any case, to wait until the end of the process before to continue, process.poll() is used
            while self.docking_subprocess.poll() is None:
                # 'time.sleep' instead of 'pass' is used to reduce the number of time the process is checked
                time.sleep(5)


            if Path("results_tmp_name.sd").is_file():

                if sys.platform == "darwin":

                    path_to_sdsorter = os.path.join(self.tab.docking_programs_child_tabs.docking_programs.path_to_sdsorter)
                    subprocess.run([path_to_sdsorter, "-sort", "'>  <SCORE>'", "results_tmp_name.sd", str(self.tab.last_docking.results_file_name + ".sd")])

                else:
                # subprocess.run(["obabel", "-i", "sd", "results_tmp_name.sd", "-o", "sdf", "-O", "results_tmp_name.sdf"])
                # subprocess.run(["obabel", "results_tmp_name.sdf", "-O", str(self.results_file_name + ".sdf"), "--sort", "SCORE"])

                    try:
                        f = open(str(self.tab.last_docking.results_file_name + ".sd"), "w")
                        subprocess.run(["sdsort", "-n", "-fSCORE", "results_tmp_name.sd"], stdout = f, check = True)
                        f.close()
                    except:
                        f.close()
                        os.remove(str(self.tab.last_docking.results_file_name + ".sd"))
                        os.rename("results_tmp_name.sd", str(self.tab.last_docking.results_file_name + ".sd"))

            if self.tab.last_docking.interrupt == False:

                self.check_docking_completed.emit(docking_programs.rxdock_runs,
                docking_programs.rxdock_tmp_dir)

                time.sleep(1)

                self.tab.summary_file = self.write_summary_file(receptor = receptor, ligand = ligand)

                if self.tab.docking_completed:

                    self.tab.results_file = RxDock_parse_results(self.tab,
                    main = self.main,
                    results_file_name = self.tab.last_docking.results_file_name,
                    results_dict = docking_programs.results_dict,
                    poses = self.tab.last_docking.poses,
                    ligand = ligand)

                self.update_results_tab.emit(docking_programs.RXDOCK.results_tab_ui,
                self.tab.last_docking,
                "RxDock",
                ["NAME", "POSE", "SCORE", "SCORE-inter", "SCORE-intra"])


        if self.tab.is_smina_tab:

            if cavity is not None:

                self.tab.last_docking = Smina_docking(self.tab,
                ligand = ligand,
                receptor = receptor,
                cavity = cavity,
                tmp_dir = self.main.smina_tmp_dir,
                cavity_list = self.main.ready_grid_centers[cavity],
                poses = str(self.tab.poses_box.value()),
                use_flex_smina = self.tab.use_flex_vina_cb.isChecked(),
                flex_residues = self.tab.flex_arg_edit.text(),
                exhaustiveness = self.tab.exhaustiveness_box.value(),
                buffer = self.tab.buffer_box.value(),
                energy = self.tab.energy_box.value(),
                rmsd = self.tab.rmsd_box.value(),
                scoring = self.tab.scoring_box.currentText(),
                main = self.main)

            else:
                self.tab.last_docking = Smina_docking(self.tab,
                ligand = ligand,
                receptor = receptor,
                tmp_dir = self.main.smina_tmp_dir,
                cavity = self.tab.loaded_cavities.currentText(),
                cavity_list = self.main.ready_grid_centers[self.tab.loaded_cavities.currentText()],
                poses = str(self.tab.poses_box.value()),
                use_flex_smina = self.tab.use_flex_vina_cb.isChecked(),
                flex_residues = self.tab.flex_arg_edit.text(),
                exhaustiveness = self.tab.exhaustiveness_box.value(),
                buffer = self.tab.buffer_box.value(),
                energy = self.tab.energy_box.value(),
                rmsd = self.tab.rmsd_box.value(),
                scoring = self.tab.scoring_box.currentText(),
                main = self.main)

            ferr = open('stdout.txt','w')

            # Run Docking Process in a different environment, to facilitate the interruption of the protocol
            self.docking_subprocess = subprocess.Popen(self.tab.last_docking.run_docking_smina_settings,
            preexec_fn=os.setsid,
            stdout=ferr)

            # 'Cancel' button to interrupt the process is enabled after a bit, to wait for effective begin of the process
            time.sleep(2)
            self.tab.docking_dialog.cancel_button.setEnabled(True)

            # To interrupt the process.
            # process.Popen.communicate() would wait until the end of the process before to continue.
            # Here process.Popen.communicate() is not used to ensure for the possibility to interrupt the protocol
            # In any case, to wait until the end of the process before to continue, process.poll() is used

            # In while loop, the file where the stdout is stored, is read to update progressbar
            #
            # For future developers: Use 'time.sleep(5)' if nothing is done in while loop
            # 'time.sleep(5)' instead of 'pass' is used to reduce the number of time the process is checked.

            ferr.flush()

            bar_value = 0
            current_size = os.path.getsize('stdout.txt')
            while self.docking_subprocess.poll() is None:
                time.sleep(0.7)
                writing = open('stdout.txt','r')
                for line in writing:
                    if line.startswith("*"):
                        bar_value = line.count("*")
                        self.update_single_docking_progressbar.emit(int(bar_value))
                writing.close()

            ferr.close()
            os.remove('stdout.txt')

            # 'Cancel' button to interrupt the process, is enabled after a bit, to wait for effective begin of the process
            time.sleep(2)
            self.tab.docking_dialog.cancel_button.setEnabled(True)

            if self.tab.last_docking.interrupt == False:

                self.check_docking_completed.emit(docking_programs.smina_runs,
                docking_programs.smina_tmp_dir)

                time.sleep(1)

                self.tab.summary_file = self.write_summary_file(receptor = receptor, ligand = ligand)

                if self.tab.docking_completed:

                    self.tab.results_file = Smina_parse_results(self.tab,
                    main = self.main,
                    results_file_name = self.tab.last_docking.results_file_name,
                    results_dict = docking_programs.results_dict,
                    poses = self.tab.last_docking.poses,
                    ligand = ligand)

                self.update_results_tab.emit(docking_programs.SMINA.results_tab_ui,
                self.tab.last_docking,
                "Smina",
                ["NAME", "POSE", "Affinity (kcal/mol)"])

        if self.tab.is_adfr_tab:

            self.update_single_docking_progressbar.emit(int(0))

            if cavity is not None:

                self.tab.last_docking = ADFR_docking(self.tab,
                ligand = ligand,
                receptor = receptor,
                cavity = cavity,
                cavity_list = self.main.ready_grid_centers[cavity],
                tmp_dir = self.main.adfr_tmp_dir,
                ga_evol = str(self.tab.ga_evol.value()),
                ga_threshold = str(self.tab.ga_threshold.value()),
                max_gen = str(self.tab.max_gen.value()),
                buffer = str(self.tab.buffer_box.value()),
                use_flex = self.tab.use_flex_vina_cb.isChecked(),
                flex_residues = self.tab.flex_arg_edit.text(),
                main = self.main)

            else:
                self.tab.last_docking = ADFR_docking(self.tab,
                ligand = ligand,
                receptor = receptor,
                tmp_dir = self.main.adfr_tmp_dir,
                cavity = self.tab.loaded_cavities.currentText(),
                cavity_list = self.main.ready_grid_centers[self.tab.loaded_cavities.currentText()],
                ga_evol = str(self.tab.ga_evol.value()),
                ga_threshold = str(self.tab.ga_threshold.value()),
                max_gen = str(self.tab.max_gen.value()),
                buffer = str(self.tab.buffer_box.value()),
                use_flex = self.tab.use_flex_vina_cb.isChecked(),
                flex_residues = self.tab.flex_arg_edit.text(),
                main = self.main)

            self.update_progress_text.emit("Generating Grid Map. Docking " + str(int(num)) + " of " + str(self.number_of_dockings_to_do))

            ferr = open('stdout.txt','w')

            # Run Grid creation Process in a different environment, to facilitate the interruption of the protocol
            if sys.platform == "win32":
                self.docking_subprocess = subprocess.Popen(self.tab.last_docking.generate_grid_adfr_settings,
                creationflags=subprocess.CREATE_NEW_PROCESS_GROUP,
                stdout=ferr,
                shell = True)

            else:
                self.docking_subprocess = subprocess.Popen(self.tab.last_docking.generate_grid_adfr_settings,
                preexec_fn=os.setsid,
                stdout=ferr)

                self.tab.docking_dialog.cancel_button.setEnabled(True)

            # 'Cancel' button to interrupt the process is enabled after a bit, to wait for effective begin of the process
            time.sleep(2)

            while self.docking_subprocess.poll() is None:
                time.sleep(2)

            ferr.close()
            os.remove('stdout.txt')


            self.update_progress_text.emit("Running docking " + str(int(num)) + " of " + str(self.number_of_dockings_to_do))

            ferr = open('stdout.txt','w')
            # Run Docking Process in a different environment, to facilitate the interruption of the protocol
            if sys.platform == "win32":
                self.docking_subprocess = subprocess.Popen(self.tab.last_docking.run_docking_adfr_settings,
                creationflags=subprocess.CREATE_NEW_PROCESS_GROUP,
                stdout=ferr,
                shell = True)
            else:
                self.docking_subprocess = subprocess.Popen(self.tab.last_docking.run_docking_adfr_settings,
                preexec_fn=os.setsid,
                stdout=ferr,
                )

                self.tab.docking_dialog.cancel_button.setEnabled(True)

            # 'Cancel' button to interrupt the process is enabled after a bit, to wait for effective begin of the process
            time.sleep(2)

            # To interrupt the process.
            # process.Popen.communicate() would wait until the end of the process before to continue.
            # Here process.Popen.communicate() is not used to ensure for the possibility to interrupt the protocol
            # In any case, to wait until the end of the process before to continue, process.poll() is used

            # In while loop, the file where the stdout is stored, is read to update progressbar
            #
            # For future developers: Use 'time.sleep(5)' if nothing is done in while loop
            # 'time.sleep(5)' instead of 'pass' is used to reduce the number of time the process is checked.

            ferr.flush()

            bar_value = 0
            current_size = os.path.getsize('stdout.txt')
            while self.docking_subprocess.poll() is None:
                time.sleep(0.7)
                writing = open('stdout.txt','r')
                for line in writing:
                    if line.startswith("*"):
                        bar_value = line.count("*")
                        self.update_single_docking_progressbar.emit(int(bar_value))
                writing.close()

            self.update_single_docking_progressbar.emit(int(51))

            ferr.close()
            os.remove('stdout.txt')

            # 'Cancel' button to interrupt the process, is enabled after a bit, to wait for effective begin of the process
            time.sleep(2)
            self.change_results_name()

            if self.tab.last_docking.interrupt == False:

                self.check_docking_completed.emit(docking_programs.adfr_runs,
                docking_programs.adfr_tmp_dir)

                time.sleep(1)

                self.tab.summary_file = self.write_summary_file(receptor = receptor, ligand = ligand)

                if self.tab.docking_completed:

                    self.tab.results_file = ADFR_parse_results(self.tab,
                    main = self.main,
                    results_file_name = self.tab.last_docking.results_file_name,
                    ligand = ligand)

                self.update_results_tab.emit(docking_programs.ADFR.results_tab_ui,
                self.tab.last_docking,
                "ADFR",
                ["NAME", "POSE", "Affinity (kcal/mol)"])


        if self.tab.docking_completed:

            self.update_data_analysis_tab.emit(self.tab.last_docking.results_file_name,
            self.tab.docking_programs_child_tabs.docking_programs.all_runs,
            self.tab.last_docking,
            self.tab.results_file
            )

        time.sleep(1)


    def change_results_name(self):

        self.file_path = os.path.join(self.tab.docking_programs_child_tabs.docking_programs.adfr_tmp_dir, self.tab.last_docking.results_file_name_ext)

        for file in os.listdir(self.tab.docking_programs_child_tabs.docking_programs.adfr_tmp_dir):
            temp = re.search("_out", file)
            temp2 = re.search("_grid.log", file)
            if temp:
                os.rename(file, self.file_path)
            if temp2:
                os.rename(file, str(file.replace("grid.log", "")) + "log.txt")


    def check_if_docking_completed(self, runs, tmp_dir):

        file_path = os.path.join(tmp_dir, self.tab.last_docking.results_file_name_ext)

        if Path(file_path).is_file():

            if os.path.getsize(file_path):
                self.docking_completed = True
                runs += 1

            else:
                self.docking_completed = False

                if len(self.ligands_to_dock) > 1 or len(self.receptors_to_dock) > 1:

                    os.remove(file_path)
                    runs += 1

                else:
                    os.remove(file_path)
                    runs += 1
                    QtWidgets.QMessageBox.warning(self.tab.docking_programs_child_tabs, "", str("Something went wrong during Docking. \nPlease check LOG files."))

        else:
            self.docking_completed = False
            QtWidgets.QMessageBox.warning(self.tab.docking_programs_child_tabs, "", str("Something went wrong during Docking. \nPlease check LOG files."))
            runs += 1


    def write_summary_file(self, receptor, ligand):

        Check_current_tab.check_docking_program_current_tab(self.tab)

        summary_list = []

        #states = cmd.count_states(self.last_docking.results_file_name)

        ligand = str("Ligand: " + ligand)
        receptor = str("Receptor: " + receptor)
        cavity = str("Cavity: " + self.tab.last_docking.cavity_name)

        if self.tab.is_smina_tab and self.tab.last_docking.grid is None:
            cavity = str("Cavity: " + self.tab.last_docking.cavity_name + " with buffer " + str(self.tab.last_docking.buffer))
        #poses = str("Generated Poses: " + str(states))

        if self.tab.is_adfr_tab:
            ga_evol = str("Number of GA evolutions: " + self.tab.last_docking.ga_evol)
            ga_threshold = str("GA evaluation stoped after " + self.tab.last_docking.ga_threshold + " \nnumber of generations with no improvement in best energy in all clusters")
            max_gen = str("Maximum number of generations: " + self.tab.last_docking.max_gen)


        if self.tab.is_vina_tab or self.tab.is_smina_tab:
            exhaustiveness = str("Exhaustiveness: " + str(self.tab.last_docking.exhaustiveness))
            energy = str("Energy Range: " + str(self.tab.last_docking.energy))

            summary_list.extend([exhaustiveness, energy])

            if self.tab.last_docking.use_flex_protocol:
                flex = str("Flexible side chains: " + str(self.tab.last_docking.flex_residues))

                summary_list.extend([flex])

        summary_list.extend([ligand, receptor, cavity])

        textfilename = str(self.tab.last_docking.results_file_name + "_summary.txt")
        textfile = open(textfilename, "w")

        for element in summary_list:

            textfile.write(element + "\n")

        textfile.close()

        return textfilename


    def create_results_frame(self, results_tab, docking_process, program, columns_names):

        if self.docking_completed:

            results_tab.tab_widget = QtWidgets.QWidget()

            results_tab.tab_layout = QtWidgets.QGridLayout()
            results_tab.tab_widget.setLayout(results_tab.tab_layout)
            results_tab.result_tabs.addTab(results_tab.tab_widget, docking_process.results_file_name)
            results_tab.results_group_layout.addWidget(results_tab.result_tabs, 0, 0)

            self.results_obj_name = docking_process.results_file_name
            self.results_file_name = docking_process.results_file_name_ext
            self.log_file_name = docking_process.log_file_name

            if type(self.results_file_name) is list:
                for i in self.results_file_name:
                    self.tab.docking_programs_child_tabs.docking_programs.results_dict[i] = {}

            else:
                self.tab.docking_programs_child_tabs.docking_programs.results_dict[self.results_file_name] = {}

            self.results_frame = ResultsFrame(parent=None,
            main_window=self.tab.docking_programs_child_tabs,
            results_file = self.results_file,
            results_file_name=self.results_file_name,
            results_obj_name=self.results_obj_name,
            results_data = self.results_file.results_data,
            log_file = self.log_file_name,
            program = program,
            columns_names = columns_names)

            results_tab.tab_layout.addWidget(self.results_frame)

            ## LOAD THE RESULTS IN PYMOL
            if type(self.results_file_name) is list:
                for i in self.results_file_name:
                    cmd.load(i, self.results_obj_name)

            else:
                cmd.load(self.results_file_name, self.results_obj_name)

            if program == "Vina":
                cmd.group("Vina", members=self.results_obj_name, action='auto')

            if program == "RxDock":
                cmd.group("RxDock", members=self.results_obj_name, action='auto')

            if program == "Smina":
                cmd.group("Smina", members=self.results_obj_name, action='auto')

            if program == "ADFR":
                cmd.group("ADFR", members=self.results_obj_name, action='auto')

        else:

            results_tab.tab_widget = QtWidgets.QWidget()

            results_tab.tab_layout = QtWidgets.QGridLayout()
            results_tab.tab_widget.setLayout(results_tab.tab_layout)
            results_tab.result_tabs.addTab(results_tab.tab_widget, docking_process.results_file_name)
            results_tab.results_group_layout.addWidget(results_tab.result_tabs, 0, 0)
            results_tab.save_all_files.setEnabled(True)
            results_tab.save_all_files.show()

            self.results_obj_name = docking_process.results_file_name
            self.results_file_name = docking_process.results_file_name_ext
            self.log_file_name = docking_process.log_file_name

            self.tab.docking_programs_child_tabs.docking_programs.results_dict[self.results_file_name] = {}

            self.results_frame = ResultsFrame(parent=None,
            main_window=self.tab.docking_programs_child_tabs,
            results_file = None,
            results_file_name=self.results_file_name,
            results_obj_name=self.results_obj_name,
            log_file = self.log_file_name,
            results_data = None,
            program = program,
            columns_names = columns_names)

            results_tab.tab_layout.addWidget(self.results_frame)

        return self.results_frame

    def check_flex_vina_input_func(self):
            return self.is_valid_input

    def get_rxdock_docking_protocol(self):

        # Set RxDock options
        self.tab.simple_docking = False
        self.tab.pharma_restrains = False
        self.tab.tethered_docking = False

        if self.tab.use_pharma_restrain_cb.isChecked():
            self.tab.pharma_restrains = True
