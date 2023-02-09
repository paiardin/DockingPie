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

from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout
import random
import numpy as np

import statistics as stat


class Import_from_pymol_window_qt(QtWidgets.QMainWindow):

    middle_layout_type = "qform"

    def __init__(self, parent, selections_list,
                 title="New Window",
                 upper_frame_title="New Window Sub-title",
                 submit_command=None, submit_button_text="Import",
                 with_scroll=True,
                 # geometry=None
                 ):

        super().__init__(parent)

        #------------------------
        # Configure the window. -
        #------------------------

        # Command executed when pressing on the main button of the window.
        self.submit_command = submit_command
        self.selections_list = selections_list

        # Configure the window.
        self.setWindowTitle(title)
        # if geometry is not None:
        #     self.setGeometry(*geometry)

        # Sets the central widget.
        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)

        # The window has a main vbox layout.
        self.main_vbox = QtWidgets.QVBoxLayout()


        #---------------
        # Upper frame. -
        #---------------

        self.upper_frame_title = QtWidgets.QLabel(upper_frame_title)
        self.main_vbox.addWidget(self.upper_frame_title)


        #----------------
        # Middle frame. -
        #----------------

        # Widget that contains the collection of Vertical Box.
        self.middle_widget = QtWidgets.QWidget()
        # The Vertical Box that contains other widgets to be displayed in the window.
        self.middle_vbox = QtWidgets.QVBoxLayout()

        self.middle_scroll = QtWidgets.QScrollArea()
        self.middle_scroll.setWidgetResizable(True)
        self.middle_scroll.setWidget(self.middle_widget)

        self.main_vbox.addWidget(self.middle_scroll)

        self.middle_layout_type = QtWidgets.QFormLayout()
        self.middle_widget.setLayout(self.middle_layout_type)

        self.sele_checkbox_list = []
        for sele in self.selections_list:
            checkbox = QtWidgets.QCheckBox(sele)
            self.sele_checkbox_list.append(checkbox)
            self.middle_layout_type.addRow(checkbox)

        #----------------
        # Bottom frame. -
        #----------------

        self.submit_command = submit_command
        if self.submit_command is not None:
            self.main_button = QtWidgets.QPushButton(submit_button_text)
            self.main_button.clicked.connect(lambda a=None: self.submit_command())
            self.main_vbox.addWidget(self.main_button)
            self.main_button.setFixedWidth(self.main_button.sizeHint().width())

        self.select_all_btn = QtWidgets.QCheckBox("All")
        self.middle_layout_type.addWidget(self.select_all_btn)
        self.select_all_btn.clicked.connect(self.get_state)


        # Sets the main vertical layout.
        self.central_widget.setLayout(self.main_vbox)
        self.main_vbox.setAlignment(self.main_button, QtCore.Qt.AlignCenter)



    def get_objects_to_import(self):
        sele_list = []
        for sele, checkbox in zip(self.selections_list, self.sele_checkbox_list):
            if checkbox.isChecked():
                sele_list.append(sele)
        return sele_list


    def get_state(self):

        if self.select_all_btn.isChecked():
            self.all_check = True
        else:
            self.all_check = False

        self.all_func()


    def all_func(self):

        for sele, checkbox in zip(self.selections_list, self.sele_checkbox_list):
            if self.all_check:
                if checkbox.isChecked():
                    pass
                else:
                    checkbox.setChecked(True)
            else:
                if checkbox.isChecked():
                    checkbox.setChecked(False)
                else:
                    pass




class NewWindow(QtWidgets.QMainWindow):


    middle_layout_type = "qform"
    is_rxdock_window = True


    def __init__(self, parent,
                 title="New Window",
                 upper_frame_title="New Window Sub-title",
                 submit_command=None, submit_button_text="Submit",
                 with_scroll=True,
                 # geometry=None
                 ):

        super().__init__(parent)

        #------------------------
        # Configure the window. -
        #------------------------

        # Command executed when pressing on the main button of the window.
        self.submit_command = submit_command

        # Configure the window.
        self.setWindowTitle(title)
        # if geometry is not None:
        #     self.setGeometry(*geometry)

        # Sets the central widget.
        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)

        # The window has a main vbox layout.
        self.main_vbox = QtWidgets.QVBoxLayout()


        #---------------
        # Upper frame. -
        #---------------

        self.upper_frame_title = QtWidgets.QLabel(upper_frame_title)
        self.main_vbox.addWidget(self.upper_frame_title)

        #----------------
        # Middle frame. -
        #----------------

        # Widget that contains the collection of Vertical Box.
        self.middle_widget = QtWidgets.QWidget()
        # The Vertical Box that contains other widgets to be displayed in the window.
        self.middle_vbox = QtWidgets.QVBoxLayout()

        self.middle_scroll = QtWidgets.QScrollArea()
        self.middle_scroll.setWidgetResizable(True)
        self.middle_scroll.setWidget(self.middle_widget)

        self.main_vbox.addWidget(self.middle_scroll)

        self.middle_layout_type = QtWidgets.QGridLayout()
        self.middle_widget.setLayout(self.middle_layout_type)

        #----------------
        # Bottom frame. -
        #----------------

        self.submit_command = submit_command
        if self.submit_command is not None:
            self.main_button = QtWidgets.QPushButton(submit_button_text)
            self.main_button.clicked.connect(lambda a=None: self.submit_command())
            self.main_vbox.addWidget(self.main_button)
            self.main_button.setFixedWidth(self.main_button.sizeHint().width())
            self.main_vbox.setAlignment(self.main_button, QtCore.Qt.AlignCenter)


        # Sets the main vertical layout.
        self.central_widget.setLayout(self.main_vbox)


class InfoWindow(QtWidgets.QMainWindow):


    middle_layout_type = "qform"
    is_rxdock_window = True


    def __init__(self, parent,
                 title="New Window",
                 upper_frame_title="New Window Sub-title",
                 submit_command=None, submit_button_text="Submit",
                 with_scroll=True,
                 # geometry=None
                 ):

        super().__init__(parent)

        #------------------------
        # Configure the window. -
        #------------------------

        # Command executed when pressing on the main button of the window.
        self.submit_command = submit_command

        # Configure the window.
        self.setWindowTitle(title)
        # if geometry is not None:
        #     self.setGeometry(*geometry)

        # Sets the central widget.
        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)

        # The window has a main vbox layout.
        self.main_vbox = QtWidgets.QVBoxLayout()


        #---------------
        # Upper frame. -
        #---------------

        self.upper_frame_title = QtWidgets.QLabel(upper_frame_title)
        self.main_vbox.addWidget(self.upper_frame_title)

        #----------------
        # Middle frame. -
        #----------------

        # Widget that contains the collection of Vertical Box.
        self.middle_widget = QtWidgets.QWidget()
        # The Vertical Box that contains other widgets to be displayed in the window.
        self.middle_vbox = QtWidgets.QVBoxLayout()

        self.middle_scroll = QtWidgets.QScrollArea()
        self.middle_scroll.setWidgetResizable(True)
        self.middle_scroll.setWidget(self.middle_widget)

        self.main_vbox.addWidget(self.middle_scroll)

        self.middle_layout_type = QtWidgets.QGridLayout()
        self.middle_widget.setLayout(self.middle_layout_type)

        #----------------
        # Bottom frame. -
        #----------------

        self.submit_command = submit_command
        if self.submit_command is not None:
            self.main_button = QtWidgets.QPushButton(submit_button_text)
            self.main_button.clicked.connect(lambda a=None: self.submit_command())
            self.main_vbox.addWidget(self.main_button)
            self.main_button.setFixedWidth(self.main_button.sizeHint().width())
            self.main_vbox.setAlignment(self.main_button, QtCore.Qt.AlignCenter)

        # Sets the main vertical layout.
        self.central_widget.setLayout(self.main_vbox)


class WelcomeWindow(QtWidgets.QMainWindow):

    middle_layout_type = "qform"

    def __init__(self, parent,
                 title="New Window",
                 upper_frame_title="New Window Sub-title",
                 submit_command=None, submit_button_text="Submit",
                 with_scroll=True,
                 # geometry=None
                 ):

        super().__init__(parent)

        #------------------------
        # Configure the window. -
        #------------------------

        # Command executed when pressing on the main button of the window.
        self.submit_command = submit_command

        # Configure the window.
        self.setWindowTitle(title)
        # if geometry is not None:
        #     self.setGeometry(*geometry)

        # Sets the central widget.
        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)

        # The window has a main vbox layout.
        self.main_vbox = QtWidgets.QVBoxLayout()


        #---------------
        # Upper frame. -
        #---------------

        self.upper_frame_title = QtWidgets.QLabel(upper_frame_title)
        self.main_vbox.addWidget(self.upper_frame_title)

        #----------------
        # Middle frame. -
        #----------------

        # Widget that contains the collection of Vertical Box.
        self.middle_widget = QtWidgets.QWidget()
        # The Vertical Box that contains other widgets to be displayed in the window.
        self.middle_vbox = QtWidgets.QVBoxLayout()

        self.middle_scroll = QtWidgets.QScrollArea()
        self.middle_scroll.setWidgetResizable(True)
        self.middle_scroll.setWidget(self.middle_widget)

        self.main_vbox.addWidget(self.middle_scroll)

        self.middle_layout_type = QtWidgets.QGridLayout()
        self.middle_widget.setLayout(self.middle_layout_type)

        #----------------
        # Bottom frame. -
        #----------------

        self.submit_command = submit_command
        if self.submit_command is not None:
            self.main_button = QtWidgets.QPushButton(submit_button_text)
            self.main_button.clicked.connect(lambda a=None: self.submit_command())
            self.main_vbox.addWidget(self.main_button)
            self.main_button.setFixedWidth(self.main_button.sizeHint().width())


        # Sets the main vertical layout.
        self.central_widget.setLayout(self.main_vbox)
        self.main_vbox.setAlignment(self.main_button, QtCore.Qt.AlignCenter)



class WelcomeWindow2(QtWidgets.QDialog):

    """
    A dialog from within to initialize Thread for Docking Processes
    """

    def __init__(self, parent):

        QtWidgets.QDialog.__init__(self, parent=parent)

        self.initUI()

        # self.docking_thread = Dockings_thread(self)
        #
        # self.docking_thread.update_progressbar.connect(self.on_update_progressbar)
        # self.docking_thread.update_single_docking_progressbar.connect(self.on_update_single_docking_progressbar)
        # self.docking_thread.set_params(self.tab,
        # number_of_dockings_to_do = number_of_dockings_to_do,
        # ligands_to_dock = ligands_to_dock,
        # receptors_to_dock = receptors_to_dock,
        # dockings_to_do = dockings_to_do)
        # self.docking_thread.emit_exception.connect(self.on_emit_exception)
        # self.docking_thread.update_results_tab.connect(self.on_update_results_tab)
        # self.docking_thread.check_docking_completed.connect(self.on_checking_docking_completed)
        # self.docking_thread.update_progress_text.connect(self.on_update_progress_text)
        # self.docking_thread.process_completed.connect(self.on_process_completed)
        # self.docking_thread.user_asks_to_interrupt.connect(self.on_interrupting_process)
        # self.docking_thread.update_data_analysis_tab.connect(self.on_update_data_analysis_tab)
        #
        # ### When multiple dockings must be run:
        # # this variable is used to stop prior the starting of the next process if the user clicks on 'Cancel' button
        # self.tab.user_asks_to_interrupt = False

    def initUI(self):

        # Check_current_tab.check_docking_program_current_tab(self.tab)
        #
        # self.setWindowTitle('Initializing Docking Process ...')
        # self.setWindowFlags(self.windowFlags() | QtCore.Qt.CustomizeWindowHint)
        # self.setWindowFlags(self.windowFlags() & ~QtCore.Qt.WindowMaximizeButtonHint)
        # self.setWindowFlags(self.windowFlags() & ~QtCore.Qt.WindowCloseButtonHint)
        #
        vertical_layout = QtWidgets.QVBoxLayout()
        # self.docking_progressbar = QtWidgets.QProgressBar(self)
        # # self.progress.setGeometry(0, 0, 340, 25)
        # self.docking_progressbar.setMaximum(self.maximum_prog_bar)
        #
        # progressbar_text = "Click the button to start running the Docking Processes"
        #
        # self.docking_progressbar.setValue(0)
        # vertical_layout.addWidget(self.docking_progressbar)
        # self.progressbar_label = QtWidgets.QLabel(progressbar_text)
        # vertical_layout.addWidget(self.progressbar_label)
        #
        # ## VINA DOCKING SETUP
        # if self.tab.is_rxdock_tab:
        #     pass
        # else:
        #     self.single_docking_progressbar = QtWidgets.QProgressBar(self)
        #     self.single_docking_progressbar.setMaximum(51)
        #     self.single_docking_progressbar.setValue(0)
        #     vertical_layout.addWidget(self.single_docking_progressbar)
        #
        # # Button for starting the installation.
        horizontal_layout = QtWidgets.QHBoxLayout()
        # start_button_text = 'Start'

        self.start_button = QtWidgets.QPushButton("prova", self)
        # self.start_button.setStyleSheet(label_style_2)
        # self.start_button.clicked.connect(self.on_button_click)
        horizontal_layout.addWidget(self.start_button)
        #
        # horizontal_layout.addStretch(1)

        # Button for closing dialog before the process starts.
        # self.close_button = QtWidgets.QPushButton('Close', self)
        # # self.cancel_button.setStyleSheet(label_style_2)
        # self.close_button.clicked.connect(self.close_dialog_prior_process)
        # horizontal_layout.addWidget(self.close_button)
        #
        # horizontal_layout.addStretch(1)
        #
        # # Button for canceling the installation.
        # self.cancel_button = QtWidgets.QPushButton('Interrupt', self)
        # # self.cancel_button.setStyleSheet(label_style_2)
        # self.cancel_button.clicked.connect(self.on_cancel_button_click)
        # self.cancel_button.setEnabled(False)
        # horizontal_layout.addWidget(self.cancel_button)

        vertical_layout.addLayout(horizontal_layout)

        self.setLayout(vertical_layout)
