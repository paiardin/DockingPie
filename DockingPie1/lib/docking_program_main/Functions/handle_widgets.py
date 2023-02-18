## Copyright 2022 by Serena Rosignoli. All rights reserved.
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
from pymol import cmd
from pymol.cgo import *
from pymol.vfont import plain
from pymol import stored
from pymol import viewing

from pymol.Qt import QtWidgets, QtCore, QtGui

import os
import shutil
import warnings
import math
import subprocess
import statistics
import re

from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.PDBIO import Select

from lib.docking_program_main.Functions.pymol_interactions import *


class HandleWidgets:

    """
    A class to represent handling of PyQt5 Widgets - Work in progress
    """

    def __init__(self,
    tab,
    dict = {},
    widget = None,
    item_name = None,
    widgets_list = [],
    input_widget = None,
    output_widget = None,
    prepared_objects_list = [],
    format = None,
    tmp_dir = None
    ):

        self.tab = tab


    def remove_frame(self, dict):

        qm = QtWidgets.QMessageBox.question(self,'', "Are you sure you want to remove this item?", QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

        if qm == QtWidgets.QMessageBox.Yes:
            text = self.strc_checkbox.text()
            new_text = text.split()[0]
            frame = dict[new_text]["frame"]
            self.layout = frame.structure_frame_layout

            for i in reversed(range(self.layout.count())):
                self.layout.itemAt(i).widget().deleteLater()

            dict.pop(new_text, None)
            frame.deleteLater()

        else:
            pass


    def remove_general_frame(self, layout, frame):

        for i in reversed(range(layout.count())):
            if str(type(layout.itemAt(i))) == "<class 'PyQt5.QtWidgets.QSpacerItem'>":
                pass
            else:
                layout.itemAt(i).widget().deleteLater()

        frame.deleteLater()


    def combobox_orient_current_item_pymol(self, widget):

        ### To zoom/orient QComboBox CurrentItem in PyMOL ###

        HandleWidgets.combobox_check_if_empty(self = self,
        widgets_list = [widget])

        if self.is_empty:
            pass

        else:
            item = widget.currentText()
            grid = re.search("Grid Center_", item)

            if grid:
                coord_dict = self.docking_programs_child_tabs.docking_programs.ready_grid_centers[item]
                self.show_crisscross_changed(1, 1, 1)
                self.calculate_box_changed(coord_dict[0], coord_dict[1], coord_dict[2], coord_dict[6], int(float(coord_dict[3])), int(float(coord_dict[4])), int(float(coord_dict[5])))

            # else:
            #     PyMOL_Zoom_Orient_Show_Hide(self, item, show_only = False)


    def combobox_check_existing_item(self, widget, item_name):

        ### To check if an item is already present in a QComboBox ###

        index = widget.findText(item_name, QtCore.Qt.MatchFixedString)

        if index >= 0:
            widget.removeItem(index)
            widget.addItem(item_name)

        else:
            widget.addItem(item_name)

        widget.setCurrentText(item_name)


    def combobox_check_if_empty(self, widgets_list):

        ### To check if QComboBox is empty ###

        self.is_empty = False

        for widgets in widgets_list:
            if str(type(widgets)) == "<class 'PyQt5.QtWidgets.QScrollArea'>":
                list = [str(widgets.widget().layout().itemAt(i).widget().text()) for i in range(widgets.widget().layout().count())]
            else:
                list = [widgets.itemText(i) for i in range(widgets.count())]

            if not list:
                self.is_empty = True

        return self.is_empty


    def add_to_other_tab(self, input_widget, output_widget):

        ### To add objects from an input_widget to an output_widget in another tab ###


        typeof = str(type(input_widget))
        output_typeof = (str(type(output_widget)))

        if typeof == "<class 'PyQt5.QtWidgets.QListWidget'>" and output_typeof == "<class 'PyQt5.QtWidgets.QListWidget'>":
            for items in input_widget.selectedItems():
                output_widget.addItem(items.text())

        elif typeof == "<class 'PyQt5.QtWidgets.QListWidget'>" and output_typeof == "<class 'PyQt5.QtWidgets.QComboBox'>":

            for items in input_widget.selectedItems():
                index = output_widget.findText(items.text(), QtCore.Qt.MatchFixedString)

                if index >= 0:
                    output_widget.removeItem(index)
                    output_widget.addItem(items.text())
                else:
                    output_widget.addItem(items.text())


        elif typeof == "<class 'PyQt5.QtWidgets.QComboBox'>" and output_typeof == "<class 'PyQt5.QtWidgets.QComboBox'>":

            index = output_widget.findText(input_widget.currentText(), QtCore.Qt.MatchFixedString)

            if index >= 0:
                output_widget.removeItem(index)
                output_widget.addItem(input_widget.currentText())

            else:
                output_widget.addItem(input_widget.currentText())


        elif typeof == "<class 'PyQt5.QtWidgets.QListWidget'>" and output_typeof == "<class 'PyQt5.QtWidgets.QScrollArea'>":

            layout = output_widget.widget().layout()

            for items in input_widget.selectedItems():
                if layout.count() == 0:
                    layout.addWidget(QtWidgets.QCheckBox(items.text()))

                else:
                    tmp_list = []
                    for i in reversed(range(layout.count())):
                        text = layout.itemAt(i).widget().text()
                        tmp_list.append(text)

                    exist_count = tmp_list.count(items.text())

                    if exist_count > 0:
                        pass
                    else:
                        layout.addWidget(QtWidgets.QCheckBox(items.text()))


    def remove_dialog(self, prepared_objects_list, widgets_list, format, tmp_dir):

        ### To remove items from listwidgets of Receptor and Ligand tabs ###
        # When an item is removed from the list widget it is removed also from:
        # - the list of prepared object,
        # - the docking tab,
        # - tmp directory
        # - PyMOL

        self.check_list_widget_empty(self.listwidget)

        if self.is_empty:
            pass

        else:

            qm = QtWidgets.QMessageBox.question(self,'', "Are you sure to remove?", QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

            if qm == QtWidgets.QMessageBox.Yes:
                for items in self.listwidget.selectedItems():

                    # remove from listwidget
                    self.listwidget.takeItem(self.listwidget.indexFromItem(items).row())

                    # remove from the list of prepared structures
                    prepared_objects_list.remove(items.text())

                    # remove from tmp_dir
                    name = str(items.text() + format)
                    path = os.path.join(tmp_dir, name)
                    if os.path.isfile(path):
                        os.remove(path)

                    # remove from docking tab
                    for widgets in widgets_list:
                        if str(type(widgets)) == "<class 'PyQt5.QtWidgets.QComboBox'>":
                            index = widgets.findText(items.text())
                            widgets.removeItem(index)
                        else:
                            for i in reversed(range(widgets.widget().layout().count())):

                                widget = widgets.widget().layout().itemAt(i).widget()

                                if str(type(widget)) == "<class 'PyQt5.QtWidgets.QWidget'>":
                                    pass
                                else:
                                    if str(widget.text()) == items.text():
                                        widgets.widget().layout().itemAt(i).widget().deleteLater()

                    # remove prepared receptor from PyMOL
                    try:
                        cmd.delete(items.text())
                    except pymol.CmdException as error:
                        print("File not found")
            else:
                pass


    def check_list_widget_empty(self, widget):

        ### To check if a QListWidget is empty

        self.is_empty = False

        n = widget.count()
        if n == 0:
            self.is_empty = True
