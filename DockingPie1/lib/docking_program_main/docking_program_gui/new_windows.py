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


        # Sets the main vertical layout.
        self.central_widget.setLayout(self.main_vbox)
        self.main_vbox.setAlignment(self.main_button, QtCore.Qt.AlignCenter)
