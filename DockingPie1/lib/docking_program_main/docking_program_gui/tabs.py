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
from lib.docking_program_main.docking_program_gui.new_windows import NewWindow, Import_from_pymol_window_qt, InfoWindow, WelcomeWindow2
from lib.docking_program_main.docking_program_gui.dialogs import *
from lib.docking_program_main.docking_program_gui.widgets_utilities import *
from lib.docking_program_main.Functions.threads import Protocol_exec_dialog
from lib.docking_program_main.Functions.general_docking_func import Generate_Object, Calculate_RMSD, update_widget_with_pymol_object, get_current_sele, PDBQT_OptionsWindows
from lib.docking_program_main.Functions.general_functions import OpenFromFile, Check_current_tab, Save_to_Csv, SelectAll, check_configuration
from lib.docking_program_main.Functions.installer import Installation, External_tools_installation_thread, External_tools_download_thread, External_components_dialog
from lib.docking_program_main.Functions.dockings_thread import _dialog_mixin, Dockings_dialog, Dockings_thread
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


class ConsensusScoringTab(QtWidgets.QWidget, PyMOLInteractions):

    def __init__(self, main_window):
        super().__init__(main_window)
        self.docking_programs = main_window

        self.initUI()

    def initUI(self):
        # Tab Layout
        self.layout_data_analysis_tab = QtWidgets.QGridLayout()

        # Receptor QGroupBox
        self.receptor_group_box = QtWidgets.QGroupBox("Receptor")
        self.layout_data_analysis_tab.addWidget(self.receptor_group_box, 0, 0)
        self.receptor_group_box.setLayout(QtWidgets.QGridLayout())

        self.receptor_cb = QtWidgets.QComboBox()
        self.rec_update = QtWidgets.QPushButton("Update from PyMOL")
        self.rec_ao = QtWidgets.QPushButton("Advanced Options")
        self.rec_clear = QtWidgets.QPushButton("Clear")
        self.rec_clear.clicked.connect(self.clear_rec_widg)
        self.rec_update.clicked.connect(self.rec_update_func)
        self.receptor_group_box.layout().addWidget(self.receptor_cb, 0, 0, 1, 2)
        self.receptor_group_box.layout().addWidget(self.rec_update, 1, 0)
        self.receptor_group_box.layout().addWidget(self.rec_clear, 1, 1)
        self.receptor_group_box.layout().addWidget(self.rec_ao, 2, 0, 1, 2)

         # Ligands QGroupBox
        self.ligand_group_box = QtWidgets.QGroupBox("Ligand(s)")
        self.layout_data_analysis_tab.addWidget(self.ligand_group_box, 1, 0, 2, 1)
        self.ligand_group_box.setLayout(QtWidgets.QGridLayout())

        self.ligand_cb = QtWidgets.QListWidget()
        self.lig_update = QtWidgets.QPushButton("Update from PyMOL")
        self.lig_ao = QtWidgets.QPushButton("Advanced Options")
        self.lig_clear = QtWidgets.QPushButton("Clear")
        self.lig_clear.clicked.connect(self.clear_lig_widg)
        self.lig_update.clicked.connect(self.lig_update_func)
        self.ligand_group_box.layout().addWidget(self.ligand_cb, 0, 0, 1, 2)
        self.ligand_group_box.layout().addWidget(self.lig_update, 1, 0)
        self.ligand_group_box.layout().addWidget(self.lig_clear, 1, 1)
        self.ligand_group_box.layout().addWidget(self.lig_ao, 2, 0, 1, 2)

         # Grid QGroupBox
        self.grid_group_box = QtWidgets.QGroupBox("Grid Settings")
        self.layout_data_analysis_tab.addWidget(self.grid_group_box, 3, 0, 1, 3)
        self.grid_layout = QtWidgets.QGridLayout()
        self.grid_layout.setAlignment(QtCore.Qt.AlignCenter)
        self.grid_group_box.setLayout(self.grid_layout)

        self.spacing_label, self.spacing_widg, self.x_dim_label, self.x_dim_widg, self.y_dim_label, self.y_dim_widg, self.z_dim_label, self.z_dim_widg = grid_dimensions(self.scroll_changed)
        self.x_label, self.x_widg, self.y_label, self.y_widg, self.z_label, self.z_widg = grid_position(self.scroll_changed)

        self.update_grid = QtWidgets.QPushButton("Get coords from PyMOL object")
        self.grid_selections = QtWidgets.QComboBox()
        self.update_grid.clicked.connect(self.get_pymol_obj_for_grid)
        self.grid_layout.addWidget(self.update_grid, 3, 0)
        self.grid_layout.addWidget(self.grid_selections, 3, 1)
        self.grid_selections.currentTextChanged.connect(self.get_current_sele)
        self.grid_selections.view().pressed.connect(self.get_current_sele)

        # Grid Position widgets
        self.grid_layout.addWidget(self.x_label, 1, 1)
        self.grid_layout.addWidget(self.x_widg, 1, 2)

        self.grid_layout.addWidget(self.y_label, 1, 3)
        self.grid_layout.addWidget(self.y_widg, 1, 4)

        self.grid_layout.addWidget(self.z_label, 1, 5)
        self.grid_layout.addWidget(self.z_widg, 1, 6)

        # Grid Dimension widgets
        self.grid_layout.addWidget(self.spacing_label, 0, 0)
        self.grid_layout.addWidget(self.spacing_widg, 0, 1)

        self.grid_layout.addWidget(self.x_dim_label, 0, 2)
        self.grid_layout.addWidget(self.x_dim_widg, 0, 3)

        self.grid_layout.addWidget(self.y_dim_label, 0, 4)
        self.grid_layout.addWidget(self.y_dim_widg, 0, 5)

        self.grid_layout.addWidget(self.z_dim_label, 0, 6)
        self.grid_layout.addWidget(self.z_dim_widg, 0, 7)

        ####

         # Settings
        self.settings_group_box = QtWidgets.QGroupBox("Settings")
        self.layout_data_analysis_tab.addWidget(self.settings_group_box, 0, 1, 3, 2)
        self.settings_group_box.setLayout(QtWidgets.QGridLayout())

        self.consensus_score_cb = QtWidgets.QComboBox()
        self.settings_group_box.layout().addWidget(self.consensus_score_cb, 0, 0)

        self.docking_program_check_list = []
        idx = 0
        for dp in ["RxDock", "Vina", "Smina", "ADFR"]:
            self.cb = QtWidgets.QCheckBox(dp)
            self.docking_program_check_list.append(self.cb)
            idx += 1
            self.settings_group_box.layout().addWidget(self.cb, idx, 0)

        # Run Consensus Docking
        self.run_consensus_docking_button = QtWidgets.QPushButton("Run Consensus Button")
        self.show_results_summary = QtWidgets.QPushButton("Show Results Summary")
        self.layout_data_analysis_tab.addWidget(self.run_consensus_docking_button, 4, 0)
        self.layout_data_analysis_tab.addWidget(self.show_results_summary, 4, 1)
        self.run_consensus_docking_button.clicked.connect(self.run_consensus_job)


    def run_consensus_job(self):

        ### Get general parameters
        self.get_general_parameters()


    def get_general_parameters(self):

        self.receptor = [self.receptor_cb.currentText()] ## do as also the receptor is a list, so if in the future we would like to add the possibility to run with multiple receptors it is easier
        self.ligands = [i.text() for i in self.ligand_cb.selectedItems()]

        self.selected_docking_programs = [i.text() for i in self.docking_program_check_list if i.isChecked()]

        ### Make consensus_job_<index> directory
        self.cs_job_index = len(self.docking_programs.consensus_job_dict)
        cs_job_dir_name = "consensus_job_" + str(self.cs_job_index)
        self.consensus_job_dir = os.path.join(self.docking_programs.consensus_tmp_dir, cs_job_dir_name)
        if os.path.isdir(self.consensus_job_dir):
            shutil.rmtree(self.consensus_job_dir)
        os.mkdir(self.consensus_job_dir)

        ### Update consensus_job_dict
        self.docking_programs.consensus_job_dict[str(self.cs_job_index)] = {}
        self.docking_programs.consensus_job_dict[str(self.cs_job_index)]["directory"] = self.consensus_job_dir
        self.docking_programs.consensus_job_dict[str(self.cs_job_index)]["docking_programs"] = self.selected_docking_programs
        self.docking_programs.consensus_job_dict[str(self.cs_job_index)]["receptors"] = {}
        self.docking_programs.consensus_job_dict[str(self.cs_job_index)]["ligands"] = {}

        for rec in self.receptor:
            self.docking_programs.consensus_job_dict[str(self.cs_job_index)]["receptors"][rec] = []
        for lig in self.ligands:
            self.docking_programs.consensus_job_dict[str(self.cs_job_index)]["ligands"][lig] = []

        ### Prepare specifics for each docking program
        self.dp_specifics_dict = {}
        for dp in self.selected_docking_programs:
            self.dp_specifics_dict[dp] = {}
            self.dp_specifics_dict[dp]["directory"] = os.path.join(self.consensus_job_dir, dp)
            if dp == "Vina" or dp == "Smina" or dp == "ADFR":
                self.dp_specifics_dict[dp]["format_ligand"] = "pdb"
                self.dp_specifics_dict[dp]["format_receptor"] = "pdb"
                self.dp_specifics_dict[dp]["generate_pdbqt"] = True
            if dp == "RxDock":
                self.dp_specifics_dict[dp]["format_ligand"] = "sdf"
                self.dp_specifics_dict[dp]["format_receptor"] = "mol2"
                self.dp_specifics_dict[dp]["generate_pdbqt"] = False

            if os.path.isdir(os.path.join(self.consensus_job_dir, dp)):
                shutil.rmtree(os.path.join(self.consensus_job_dir, dp))
            os.mkdir(os.path.join(self.consensus_job_dir, dp))

        ### Run consensus
        for dp in self.selected_docking_programs:
            self.run_single_consensus(dp)


    def run_single_consensus(self, dp):
        #### TO DO - IN CONSENUS_PROTOCOL.PY

        for rec in self.receptor:

            self.prepare_receptors(rec = rec,
                                   docking_program = dp,
                                   directory = self.dp_specifics_dict[dp]["directory"],
                                   format = self.dp_specifics_dict[dp]["format_receptor"],
                                   generate_pdbqt = self.dp_specifics_dict[dp]["generate_pdbqt"])

        for lig in self.ligands:

            self.prepare_ligands(lig = lig,
                                 docking_program = dp,
                                 directory = self.dp_specifics_dict[dp]["directory"],
                                 format = self.dp_specifics_dict[dp]["format_ligand"],
                                 generate_pdbqt = self.dp_specifics_dict[dp]["generate_pdbqt"])


    def prepare_ligands(self, lig, docking_program, directory, format, generate_pdbqt):

        self.prepared_ligands = self.docking_programs.consensus_job_dict[str(self.cs_job_index)]["ligands"][lig]

        # Keep track of the imported receptors that have multiple states, just to give a warning to the user
        states = cmd.count_states(lig)
        if states > 1:
            warning = "Multiple states have been detected in PyMOL for the selected receptor. Only the current state is going to be taken into consideration"

        tmp_path_name = os.path.join(directory, str(lig + ".pdb"))

        # This save-delete-load is done to avoid errors while parsing  the files. In this way, despite the type of file that was previously loaded in PyMOL, now it is loaded as a PDB file.
        cmd.save(tmp_path_name, lig, format = 'pdb')
        cmd.delete(lig)
        cmd.load(tmp_path_name, lig)

        self.pdbqt_options_dict_lig = {}
        self.pdbqt_options_dict_lig["add_h"] = True
        self.pdbqt_options_dict_lig["none_torsions"] = False
        self.pdbqt_options_dict_lig["all_torsions"] = False
        self.pdbqt_options_dict_lig["all_but_ga"] = True

        self.generated_receptor = Generate_Object(self, main = self.docking_programs,
        prepared_objects_list = self.prepared_ligands,
        tmp_path = directory,
        format = format,
        docking_program_name = docking_program,
        generate_pdbqt = generate_pdbqt,
        object_names = [lig],
        pdbqt_dict = self.pdbqt_options_dict_lig,
        is_receptor = False,
        from_gui = False)

        # Update dict
        self.docking_programs.consensus_job_dict[str(self.cs_job_index)]["ligands"][lig].append(self.generated_receptor.new_strc_name)


    def prepare_receptors(self, rec, docking_program, directory, format, generate_pdbqt):

        self.prepared_receptors = self.docking_programs.consensus_job_dict[str(self.cs_job_index)]["receptors"][rec]

        # Keep track of the imported receptors that have multiple states, just to give a warning to the user
        states = cmd.count_states(rec)
        if states > 1:
            warning = "Multiple states have been detected in PyMOL for the selected receptor. Only the current state is going to be taken into consideration"

        tmp_path_name = os.path.join(directory, str(rec + ".pdb"))

        # This save-delete-load is done to avoid errors while parsing  the files. In this way, despite the type of file that was previously loaded in PyMOL, now it is loaded as a PDB file.
        cmd.save(tmp_path_name, rec, format = 'pdb')
        cmd.delete(rec)
        cmd.load(tmp_path_name, rec)

        # Extract some information from the input
        self.object = ObjectParser(file_name = rec,
        file_path = tmp_path_name,
        is_receptor = True)

        ###
        self.pdbqt_options_dict = {}
        self.pdbqt_options_dict["add_h"] = True
        self.pdbqt_options_dict["bonds"] = False
        self.pdbqt_options_dict["add_gast"] = False
        self.pdbqt_options_dict["remove_nonstd"] = False
        self.pdbqt_options_dict["remove_water"] = True
        self.pdbqt_options_dict["remove_lone_pairs"] = False
        self.pdbqt_options_dict["remove_non_polar_H"] = False
        self.pdbqt_options_dict["remove_non_protein"] = False

        self.generated_receptor = Generate_Object(self, main = self.docking_programs,
        prepared_objects_list = self.prepared_receptors,
        tmp_path = directory,
        format = format,
        docking_program_name = docking_program,
        generate_pdbqt = generate_pdbqt,
        object_names = [rec],
        pdbqt_dict = self.pdbqt_options_dict,
        is_receptor = True,
        from_gui = False)

        # Update dict
        self.docking_programs.consensus_job_dict[str(self.cs_job_index)]["receptors"][rec].append(self.generated_receptor.new_strc_name)


    def get_current_sele(self):

        get_current_sele(self, self.docking_programs, self.grid_selections,
                         self.x_widg,
                         self.y_widg,
                         self.z_widg,
                         self.spacing_widg,
                         self.x_dim_widg,
                         self.y_dim_widg,
                         self.z_dim_widg
                         )

    def clear_rec_widg(self):

        self.receptor_cb.clear()

    def clear_lig_widg(self):

        self.ligand_cb.clear()

    def get_pymol_obj_for_grid(self):

        update_widget_with_pymol_object(self.docking_programs, self.grid_selections, selections = True)

    def rec_update_func(self):

        update_widget_with_pymol_object(self.docking_programs, self.receptor_cb, small_molecule = False)

    def lig_update_func(self):

        update_widget_with_pymol_object(self.docking_programs, self.ligand_cb, polymer = False)


    def scroll_changed(self):

        # Docking Box is updated dinamically
        x = self.x_widg.value()
        y = self.y_widg.value()
        z = self.z_widg.value()
        spacing = self.spacing_widg.value()
        x_vis = self.x_dim_widg.value()
        y_vis = self.y_dim_widg.value()
        z_vis = self.z_dim_widg.value()

        spinbox = self.sender()

        if spinbox is self.x_widg:
            x = self.x_widg.value()
        elif spinbox is self.y_widg:
            y = self.y_widg.value()
        elif spinbox is self.z_widg:
            z = self.z_widg.value()
        elif spinbox is self.spacing_widg:
            spacing = self.spacing_widg.value()
        elif spinbox is self.x_dim_widg:
            x_vis = self.x_dim_widg.value()
        elif spinbox is self.y_dim_widg:
            y_vis = self.y_dim_widg.value()
        elif spinbox is self.z_dim_widg:
            z_vis = self.z_dim_widg.value()

        # Show the center and the box each time the values are updated
        self.show_crisscross_changed(x, y, z)
        self.calculate_box_changed(x, y, z, spacing, x_vis, y_vis, z_vis)



class Consensus_layout(QtWidgets.QWidget):

    def __init__(self, main_window):
        super().__init__(main_window)
        self.consensus_window = main_window

        self.initUI()

    def initUI(self):

        self.widget = QtWidgets.QWidget()
        self.layout = QtWidgets.QGridLayout()
        self.widget.setLayout(self.layout)

        self.list_widget = QtWidgets.QListWidget()
        self.list_widget.setSelectionMode(3)
        self.layout.addWidget(self.list_widget, 0, 0, 1, 2)

        self.remove_btn = QtWidgets.QPushButton("Remove")
        self.layout.addWidget(self.remove_btn, 1, 0, 1, 2)
        self.remove_btn.clicked.connect(lambda: self.remove(self.list_widget))

        #self.remove_btn.hide()

        self.label_poses = QtWidgets.QLabel("Poses threshold: ")
        self.label_box = QtWidgets.QComboBox()
        # self.label_box.setRange(1, 100)
        # self.label_box.setSingleStep(1)
        # self.label_box.setValue(1)

        self.label_box.addItems(["Top-ranked only", "3", "10", "All"])

        self.layout.addWidget(self.label_poses, 2, 0)
        self.layout.addWidget(self.label_box, 2, 1)


    def remove(self, list):

        for items in list.selectedItems():
            self.consensus_window.central_list_widget.addItem(items.text())
            list.takeItem(list.indexFromItem(items).row())



class DataAnalysisTab(QtWidgets.QWidget):

    """
    Layout for DataAnalysisTab.
    """

    def __init__(self, main_window):
        super().__init__(main_window)
        self.docking_programs = main_window

        self.initUI()


    def initUI(self):

        # Tab Layout
        self.layout_data_analysis_tab = QtWidgets.QGridLayout()

        # Docking Runs QGroupBox
        self.docking_runs_group = QtWidgets.QGroupBox("All Docking Processes")
        self.layout_data_analysis_tab.addWidget(self.docking_runs_group, 0, 0, 5, 1)
        self.docking_runs_group_layout = QtWidgets.QVBoxLayout()
        self.docking_runs_group.setLayout(self.docking_runs_group_layout)

        # Widgets for Docking Runs group
        self.docking_runs_widget = QtWidgets.QWidget()
        self.docking_runs_scroll = QtWidgets.QScrollArea()
        self.docking_runs_scroll.setWidgetResizable(True)
        self.docking_runs_scroll.setWidget(self.docking_runs_widget)
        # Set the layout of the Scroll Area for the Docking runs
        self.docking_runs_scroll_layout = QtWidgets.QFormLayout()
        self.docking_runs_widget.setLayout(self.docking_runs_scroll_layout)
        self.docking_runs_group_layout.addWidget(self.docking_runs_scroll)

        self.all_cb = QtWidgets.QCheckBox("All")
        self.docking_runs_scroll_layout.addRow(self.all_cb)
        self.all_cb.clicked.connect(self.get_list_of_current_cb)

        # Right part widget
        self.consensus_btn = QtWidgets.QPushButton("Consensus Scoring")
        self.layout_data_analysis_tab.addWidget(self.consensus_btn, 0, 1, 1, 3)
        self.consensus_btn.clicked.connect(self.open_consensus_scoring_window)

        ### Plot RMSD button currently not used

        # self.plot_rmsd_btn = QtWidgets.QPushButton("Plot RMSD")
        # self.layout_data_analysis_tab.addWidget(self.plot_rmsd_btn, 1, 1, 1, 3)
        # self.plot_rmsd_btn.clicked.connect(self.plot_rmsd_boxplot_func)

        # Create the Scroll Area for the Table
        self.table_widget = QtWidgets.QWidget()
        self.table_scroll = QtWidgets.QScrollArea()
        self.table_scroll.setWidgetResizable(True)
        self.table_scroll.setWidget(self.table_widget)
        self.layout_data_analysis_tab.addWidget(self.table_scroll, 2, 1, 1, 3)
        # Set the layout of the Scroll Area for the Table
        self.table_scroll_layout = QtWidgets.QGridLayout()
        self.table_widget.setLayout(self.table_scroll_layout)

        # Export btn
        self.to_csv_btn = QtWidgets.QPushButton("Export Table to \".csv\"")
        self.to_csv_btn.clicked.connect(self.save_to_csv)
        self.to_csv_btn.setEnabled(False)
        self.layout_data_analysis_tab.addWidget(self.to_csv_btn, 3, 1)

        # Save to File btn
        self.save_cons_btn = QtWidgets.QPushButton("Save Consensus Results")
        self.save_cons_btn.clicked.connect(self.save_consensus_results)
        self.save_cons_btn.setEnabled(False)
        self.layout_data_analysis_tab.addWidget(self.save_cons_btn, 4, 1, 1, 2)

        # Show matrix
        self.show_consensus_matrix_btn = QtWidgets.QPushButton("Show Consensus Matrix")
        self.show_consensus_matrix_btn.clicked.connect(self.show_consensus_matrix)
        self.show_consensus_matrix_btn.setEnabled(False)
        self.layout_data_analysis_tab.addWidget(self.show_consensus_matrix_btn, 3, 2)


    def save_consensus_results(self):

        """
        Save results of consensus.
        """

        self.table = self.consensus_run.consensus_new_table

        tmp_dir = self.docking_programs.consensus_tmp_dir
        tmp_path = os.path.join(tmp_dir, "temporary")
        os.mkdir(tmp_path)

        first_index = TableView.get_column_index_from_header(self, table = self.table, header = "Program 1")
        second_index = TableView.get_column_index_from_header(self, table = self.table, header = "Program 2")
        # Get data from index
        first_list = TableView.get_column_data_from_index(self, table = self.table, column_header_index = first_index[0])
        second_list = TableView.get_column_data_from_index(self, table = self.table, column_header_index = second_index[0])

        joined_list = first_list + second_list

        for file in joined_list:
            shutil.copy(os.path.join(tmp_dir, file), os.path.join(tmp_path, file))

        filepath = asksaveasfile_qt("Save Docking Results", name_filter="*.zip")

        if not filepath:
            shutil.rmtree(tmp_path)
            return None

        else:
            os.chdir(tmp_path)
            split_filepath = os.path.splitext(filepath)[0]
            shutil.make_archive(split_filepath, 'zip')

        os.chdir(self.docking_programs.tmp_dir_path)
        shutil.rmtree(tmp_path)


    def get_list_of_current_cb(self):

        self.list_of_cb = []

        for i in reversed(range(self.docking_runs_scroll_layout.count())):
            obj = self.docking_runs_scroll_layout.itemAt(i).widget()
            if str(type(obj)) == "<class 'lib.docking_program_main.docking_program_gui.frames.DockingsFrame'>":
                self.list_of_cb.append(obj.docking_checkbox)

        SelectAll(self, all_cb = self.all_cb,
        list_of_cb = self.list_of_cb)

    def show_consensus_matrix(self):

        self.consensus_run.consensus_matrix_window.show()


    def open_consensus_scoring_window(self):

        self.consensus_scoring_window = NewWindow(parent = self.docking_programs,
        title = "Consensus Scoring", upper_frame_title = "",
        submit_command = self.consensus_func_dialog, submit_button_text= "Start",
        with_scroll = True)

        #self.consensus_scoring_window.setFixedSize(450, 480)

        # Parameters Setting GroupBox
        self.internal_consensus_box = QtWidgets.QGroupBox("Parameters Setting")
        self.consensus_scoring_window.middle_layout_type.addWidget(self.internal_consensus_box, 1, 1)
        self.internal_consensus_layout = QtWidgets.QGridLayout()
        self.internal_consensus_box.setLayout(self.internal_consensus_layout)

        # self.consensus_info = QtWidgets.QPushButton("Info")
        # self.consensus_info.clicked.connect(self.show_consensus_info)
        # self.consensus_scoring_window.middle_layout_type.addWidget(self.consensus_info, 0, 0, 1, 2)

        self.label_rmsd = QtWidgets.QLabel("RMSD threshold: ")
        self.label_rmsd_box = QtWidgets.QSpinBox()
        self.label_rmsd_box.setRange(1, 100)
        self.label_rmsd_box.setSingleStep(1)
        self.label_rmsd_box.setValue(1)
        self.label_rmsd.setEnabled(False)
        self.label_rmsd_box.setEnabled(False)

        self.internal_consensus_layout.addWidget(self.label_rmsd, 4, 0)
        self.internal_consensus_layout.addWidget(self.label_rmsd_box, 4, 1)

        self.label_consensus_type = QtWidgets.QLabel("Consensus Score: ")
        self.box_consensus_type = QtWidgets.QComboBox()
        self.box_consensus_type.addItems(["Rank by Rank", "Average of Auto-Scaled Scores", "Z-scores"])

        self.filter_by_rmsd = QtWidgets.QCheckBox("Filter by RMSD")
        self.filter_by_rmsd.clicked.connect(self.show_rmsd_protocol_options)
        self.filter_by_rmsd.toggled.connect(self.show_rmsd_threshold_box)

        self.paired_rmsd_protocol = QtWidgets.QRadioButton("Paired RMSD protocol")
        self.clustered_rmsd_protocol = QtWidgets.QRadioButton("Clustered RMSD protocol")
        self.paired_rmsd_protocol.setEnabled(False)
        self.clustered_rmsd_protocol.setEnabled(False)
        self.paired_rmsd_protocol.setChecked(True)

        self.internal_consensus_layout.addWidget(self.label_consensus_type, 3, 0)
        self.internal_consensus_layout.addWidget(self.box_consensus_type, 3, 1)
        self.internal_consensus_layout.addWidget(self.filter_by_rmsd, 2, 0)
        # self.internal_consensus_layout.addWidget(self.paired_rmsd_protocol, 5, 0)
        # self.internal_consensus_layout.addWidget(self.clustered_rmsd_protocol, 6, 0)

        self.rxdock_frame = ConsensusFrame(parent=None, main_window=self, program = "RxDock")
        self.smina_frame = ConsensusFrame(parent=None, main_window=self, program = "Smina")
        self.vina_frame = ConsensusFrame(parent=None, main_window=self, program = "Vina")
        self.adfr_frame = ConsensusFrame(parent=None, main_window=self, program = "ADFR")

        self.internal_consensus_layout.addWidget(self.rxdock_frame, 7, 0, 1, 2)
        self.internal_consensus_layout.addWidget(self.smina_frame, 8, 0, 1, 2)
        self.internal_consensus_layout.addWidget(self.vina_frame, 9, 0, 1, 2)
        self.internal_consensus_layout.addWidget(self.adfr_frame, 10, 0, 1, 2)

        self.consensus_scoring_window.show()


    def show_rmsd_threshold_box(self):
        if self.filter_by_rmsd.isChecked():
            self.label_rmsd.setEnabled(True)
            self.label_rmsd_box.setEnabled(True)
        else:
            self.label_rmsd.setEnabled(False)
            self.label_rmsd_box.setEnabled(False)

    def show_rmsd_protocol_options(self):

        if self.filter_by_rmsd.isChecked():
            self.paired_rmsd_protocol.setEnabled(True)
            self.clustered_rmsd_protocol.setEnabled(True)
            self.paired_rmsd_protocol.setChecked(True)

        else:
            self.paired_rmsd_protocol.setEnabled(False)
            self.clustered_rmsd_protocol.setEnabled(False)
            self.paired_rmsd_protocol.setChecked(False)
            self.clustered_rmsd_protocol.setChecked(False)

    def show_consensus_info(self):

        self.consensus_info_window = NewWindow(parent = self.docking_programs,
        title = "Information", upper_frame_title = "What is a Consensus Scoring analysis?",
        submit_command = self.close_func, submit_button_text= "Close",
        with_scroll = True)

        path_to_consensus_file = os.path.join(self.docking_programs.config_path, "consensus_info.txt")

        self.consensus_info_text_area = QtWidgets.QPlainTextEdit()

        self.consensus_info_text_area.setReadOnly(True)

        self.consensus_info_text_area.setPlainText(open(str(path_to_consensus_file)).read())
        self.consensus_info_window.middle_layout_type.addWidget(self.consensus_info_text_area)

        self.consensus_info_window.show()


    def close_func(self):
        self.consensus_info_window.hide()


    def consensus_func_dialog(self):
        # Get Consensus Protocol
        if self.filter_by_rmsd.isChecked():

            if self.clustered_rmsd_protocol.isChecked():
                self.rmsd_protocol = "clustered"
            else:
                self.rmsd_protocol = "paired"

        else:
            self.rmsd_protocol = "normsd"

        # Among the runs checked by the user, divide_by_program
        self.smina_runs_list = []
        self.vina_runs_list = []
        self.rxdock_runs_list = []
        self.adfr_runs_list = []

        for key in self.docking_programs.all_runs:
            cb = self.docking_programs.all_runs[key]["dockings_frame"].docking_checkbox
            if cb.isChecked():
                self.divide_by_program(cb)

        self.consensus_scoring_window.hide()

        # Check number of runs
        s = 0
        if self.smina_runs_list:
            s += 1
        if self.vina_runs_list:
            s +=1
        if self.rxdock_runs_list:
            s += 1
        if self.adfr_runs_list:
            s+=1

        if s < 2:
            QtWidgets.QMessageBox.warning(self.docking_programs, "Warning", "At least 2 runs with different docking programs are needed to perform 'Consensus Scoring' analysis")
        else:

            self.consensus_func()


    def consensus_func(self):

        # All lists are provided to the ConsensusProtocol, even if empty
        self.consensus_run = ConsensusProtocol(self,
        smina_poset = self.smina_frame.label_box.currentText(),
        rxdock_poset = self.rxdock_frame.label_box.currentText(),
        vina_poset = self.vina_frame.label_box.currentText(),
        adfr_poset = self.adfr_frame.label_box.currentText(),
        smina_runs_list = self.smina_runs_list,
        vina_runs_list = self.vina_runs_list,
        rxdock_runs_list = self.rxdock_runs_list,
        adfr_runs_list = self.adfr_runs_list,
        rmsd_protocol = self.rmsd_protocol)


    def divide_by_program(self, cb):

        text = cb.text()
        if re.search("Smina", text):
            self.smina_runs_list.append(text)
        if re.search("RxDock", text):
            self.rxdock_runs_list.append(text)
        if re.search("Vina", text):
            self.vina_runs_list.append(text)
        if re.search("ADFR", text):
            self.adfr_runs_list.append(text)


    def save_to_csv(self):

        """
        Saves the table data to a .csv file.
        """

        self.table = self.consensus_run.consensus_new_table

        # Let the user select the filepath.
        filepath = asksaveasfile_qt("Save CSV file", name_filter="*.csv")

        if not filepath:
            return None

        try:
            # Writes a .csv file on that path.
            with open(filepath, 'w') as csv_fh:

                writer = csv.writer(csv_fh, delimiter=',', quoting=csv.QUOTE_MINIMAL)

                if self.table.row_labels is not None:
                    writer.writerow([" "] + self.table.column_labels)
                else:
                    writer.writerow(self.table.column_labels)
                for row_idx, row in enumerate(self.table.data):
                    if self.table.row_labels is not None:
                        writer.writerow([self.table.row_labels[row_idx]] + [str(v) for v in row])
                    else:
                        writer.writerow([str(v) for v in row])

        except Exception as e:
            print("- WARNING: could not write a csv file: %s" % str(e))



    def plot_rmsd_boxplot_func(self):
        pass


    def consensus_scoring_func(self, first_list_consensus, second_list_consensus):

        ## NOT USED

        os.chdir(self.docking_programs.tmp_dir_path)

        # Save the Checked Docking Runs in tmp dir
        for runs in first_list_consensus:
            cmd.save(str(os.path.join(self.docking_programs.tmp_dir_path, runs) + ".sdf"), runs, format = 'sdf', state = 0)

        for runs in second_list_consensus:
            cmd.save(str(os.path.join(self.docking_programs.tmp_dir_path, runs) + ".sdf"), runs, format = 'sdf', state = 0)

        list_of_list = []
        for first_runs in first_list_consensus:

            ref = io.loadallmols(str(first_runs + ".sdf"))						# Read crystal pose
            #The spyrmsd class needs atomic coordinates, atomic number and the molecular adjacency
            # matrix to compute the standard RMSD with spyrmsd.rmsd.symmrmsd.
            for refs in ref:
                refs.strip()

            for references in ref:
                coords_ref = references.coordinates
                anum_ref = references.atomicnums
                adj_ref = references.adjacency_matrix

                for second_runs in second_list_consensus:

                    first_ligand = self.docking_programs.all_runs[first_runs]["docking_run"].ligand_to_dock.split("_")[1]
                    second_ligand = self.docking_programs.all_runs[second_runs]["docking_run"].ligand_to_dock.split("_")[1]

                    if first_ligand == second_ligand:

                        mols = io.loadallmols(str(second_runs + ".sdf"))					# Read poses
                        for mol in mols:									# Strip Hydrogens of pose
                            mol.strip()

                        # Create Adj matrix for rxdock otuput
                        coords = [mol.coordinates for mol in mols]
                        anum = mols[0].atomicnums
                        adj = mols[0].adjacency_matrix

                    # With this information we can easily
                    # compute the RMSD between the reference molecule and all other molecules
                        RMSD = rmsd.symmrmsd(coords_ref, coords, anum_ref, anum, adj_ref, adj)
                        i = 0
                        for RMSD_Value in RMSD:

                            if RMSD_Value <= 1:
                                rmsd_value = round(RMSD_Value, 3)
                                tmp_list = []

                                consensus_score = ((ref.index(references)+1) + (i+1))/2
                                tmp_list.extend([first_runs, second_runs, first_ligand, str(rmsd_value), ref.index(references)+1, i+1, consensus_score])
                                i += 1

                                list_of_list.append(tmp_list)

                            else:
                                i += 1

        if list_of_list:
            for data in list_of_list:

                new_window = TableView(parent=self.docking_programs, data=list_of_list,
                                           row_labels=None, row_labels_height=25,
                                           column_labels=["RUN 1", "RUN2", "NAME", "RMSD", "RANKING 1", "RANKING 2", "SCORE"],
                                           sortable=False)

                new_window.itemDoubleClicked.connect(self.on_click_show_consensus)


            self.table_scroll_layout.addWidget(new_window, 2, 0)
        else:
            pass


    def on_click_show_consensus(self):
        print("ciao")



class ConfigurationTab(QtWidgets.QWidget):

    """
    Layout for ConfigurationTab.
    """

    def __init__(self, main_window):
        super().__init__(main_window)
        self.docking_programs = main_window

        self.initUI()


    def initUI(self):

        if sys.platform == "linux":
            self.dir_name = "external_tools_linux"
        if sys.platform == "darwin":
            self.dir_name = "external_tools_macOS"
        if sys.platform == "win32":
            self.dir_name = "external_tools_windows"

        self.external_tools_error = False

        # Grid Tab Layout
        self.layout_config_tab = QtWidgets.QGridLayout()

        # Widget that contains the scroll Area
        self.config_widget = QtWidgets.QWidget()

        # Create the Scroll Area and add it to the widget
        self.config_scroll = QtWidgets.QScrollArea()
        self.config_scroll.setWidgetResizable(True)
        self.config_scroll.setWidget(self.config_widget)
        self.layout_config_tab.addWidget(self.config_scroll)

        # Set the layout of the Scroll Area
        self.config_scroll_layout = QtWidgets.QGridLayout()
        self.config_widget.setLayout(self.config_scroll_layout)

        # Create a Frame for each program to install
        self.programs_to_install = {}

        self.continue_check_installation = False

        for programs in ["Config", "RxDock", "Vina", "Smina", "ADFR", "Openbabel", "sPyRMSD", "sdsorter"]:

            self.installation_frame = InstallationFrames(parent=None, main_window=self, program_name=programs)

            # Add frames to the Scroll Area
            self.config_scroll_layout.addWidget(self.installation_frame)
            self.programs_to_install[programs] = self.installation_frame

        # Functions for installation
        self.programs_to_install["Openbabel"].installation_btn.clicked.connect(lambda: self.installations_func("Openbabel"))
        self.programs_to_install["RxDock"].installation_btn.clicked.connect(lambda: self.installations_func("RxDock"))
        self.programs_to_install["Vina"].installation_btn.clicked.connect(lambda: self.installations_func("Vina"))
        self.programs_to_install["Smina"].installation_btn.clicked.connect(lambda: self.installations_func("Smina"))
        self.programs_to_install["sPyRMSD"].installation_btn.clicked.connect(lambda: self.installations_func("sPyRMSD"))
        self.programs_to_install["ADFR"].installation_btn.clicked.connect(lambda: self.installations_func("ADFR"))
        self.programs_to_install["sdsorter"].installation_btn.clicked.connect(lambda: self.installations_func("sdsorter"))

        if self.continue_check_installation:

            # Check if the programs are already installed
            try:
                self.check_installation()
            except:
                # An exception may occur if for some reasons the Config directory was found but is not complete
                self.external_tools_error = True
                QtWidgets.QMessageBox.warning(self.docking_programs, "External Tools Error", str("This error is due to a wrong Configuration of the external tools.\nA possible cause could be an interruption during download.\nPlease Configure again External Tools.\nPress \"Ok\" to continue"))
                path_to_config = os.path.join(self.docking_programs.config_path)
                adt_dir = os.path.join(path_to_config, "AutoDockTools")
                et = os.path.join(path_to_config, self.dir_name)
                flex = os.path.join(path_to_config, "prepare_flexreceptor4.py")
                ligand = os.path.join(path_to_config, "prepare_ligand4.py")
                receptor = os.path.join(path_to_config, "prepare_receptor4.py")

                if os.path.exists(adt_dir):
                    shutil.rmtree(adt_dir)
                if os.path.exists(et):
                    shutil.rmtree(et)
                if os.path.exists(flex):
                    os.remove(flex)
                if os.path.exists(ligand):
                    os.remove(ligand)
                if os.path.exists(receptor):
                    os.remove(receptor)

                self.programs_to_install["Config"].configure_btn.setEnabled(True)


        if self.continue_check_installation == False or self.external_tools_error == True:

            self.programs_to_install["RxDock"].installation_btn.setEnabled(False)
            self.programs_to_install["RxDock"].installation_btn.setText("Configure External Tools to proceed with \nRxDock installation")

            self.programs_to_install["Openbabel"].installation_btn.setEnabled(False)
            self.programs_to_install["Openbabel"].installation_btn.setText("Configure External Tools to proceed with \nOpenbabel installation")

            self.programs_to_install["Vina"].installation_btn.setEnabled(False)
            self.programs_to_install["Vina"].installation_btn.setText("Configure External Tools to proceed with \nVina installation")

            self.programs_to_install["Smina"].installation_btn.setEnabled(False)
            self.programs_to_install["Smina"].installation_btn.setText("Configure External Tools to proceed with \nSmina installation")

            self.programs_to_install["ADFR"].installation_btn.setEnabled(False)
            self.programs_to_install["ADFR"].installation_btn.setText("Configure External Tools to proceed with \nADFR installation")

            self.programs_to_install["sPyRMSD"].installation_btn.setEnabled(False)
            self.programs_to_install["sPyRMSD"].installation_btn.setText("Configure External Tools to proceed with \nsPyRMSD installation")

            self.programs_to_install["sdsorter"].installation_btn.setEnabled(False)
            self.programs_to_install["sdsorter"].installation_btn.setText("Configure External Tools to proceed with \nsdsorter installation")


    def configure_external_tools_directory(self):

        config_path = self.docking_programs.config_path

        with open(str(os.path.join(config_path, "version.txt"))) as f:
            self.files_version = f.readline().rstrip()

        try:
            urllib.request.urlopen("https://github.com/paiardin/DockingPie")
            github_accessible = True
        except:
            github_accessible = False

        if github_accessible:

            if sys.platform == "linux":
                dir_link = "https://github.com/paiardin/DockingPie/releases/download/" + str(self.files_version) + "/external_tools_linux.zip"
                dir_name = "external_tools_linux"
            if sys.platform == "darwin":
                dir_link = "https://github.com/paiardin/DockingPie/releases/download/" + str(self.files_version) + "/external_tools_macOS.zip"
                dir_name = "external_tools_macOS"
            if sys.platform == "win32":
                dir_link = "https://github.com/paiardin/DockingPie/releases/download/" + str(self.files_version) + "/external_tools_windows.zip"
                dir_name = "external_tools_windows"

        else:

            if sys.platform == "linux":
                dir_link = "http://schubert.bio.uniroma1.it/temp/DockingPie_config_file/" + str(self.files_version) + "/external_tools_linux.zip"
                dir_name = "external_tools_linux"
            if sys.platform == "darwin":
                dir_link = "http://schubert.bio.uniroma1.it/temp/DockingPie_config_file/" + str(self.files_version) +"/external_tools_macOS.zip"
                dir_name = "external_tools_macOS"
            if sys.platform == "win32":
                dir_link = "http://schubert.bio.uniroma1.it/temp/DockingPie_config_file/" + str(self.files_version) +"/external_tools_windows.zip"
                dir_name = "external_tools_windows"

        # Show the external components download dialog.
        install_dialog = External_components_dialog(self,
                                                    url=dir_link,
                                                    os_arch="64",
                                                    dir_name = dir_name,
                                                    local_mode=False)
        install_dialog.setModal(True)
        install_dialog.exec_()

        # Finishes the installation.
        if install_dialog.complete_status:
            self.check_installation()
            self.check_external_tools()

    def check_external_tools(self):

        if sys.platform == "linux":
            dir_name = "external_tools_linux"
        if sys.platform == "darwin":
            dir_name = "external_tools_macOS"
        if sys.platform == "win32":
            dir_name = "external_tools_windows"

        ext_tools_path = os.path.join(self.docking_programs.config_path, dir_name)

        if os.path.isdir(ext_tools_path):
            config_frame = self.programs_to_install['Config']
            config_frame.configure_line_edit.setText(ext_tools_path)
            config_frame.configure_btn.setEnabled(False)
            self.continue_check_installation = True
            config_frame.check_for_updates.setEnabled(True)

        else:
            config_frame.configure_line_edit.setText("Not Found")
            config_frame.configure_btn.setEnabled(True)
            self.continue_check_installation = False
            config_frame.check_for_updates.setEnabled(False)


    def close_welcome_window(self):
        self.welcome_message_window.close()

    def installations_func(self, program_to_install):

        installation = Installation(self, program_to_install)


    def check_installation(self, check_only = False):

        self.openbabel_installed = True
        self.rxdock_installed = True
        self.vina_installed = True
        self.smina_installed = True
        self.adfr_installed = True
        self.spyrmsd_installed = True
        self.sdsorter_installed = True

        ### Check Openbabel ###
        try:
            from openbabel import pybel
        except:
            self.openbabel_installed = False


        ### Check RxDock ###
        if sys.platform == "win32":
            self.rxdock_installed = False

            if not check_only:
                self.programs_to_install["RxDock"].installation_btn.setEnabled(False)
                self.programs_to_install["RxDock"].installation_btn.setText("RxDock is not available for this OS")

        else:
            try:
                result = subprocess.run(["rbdock"], stdout=subprocess.PIPE, stderr = subprocess.PIPE)
            except FileNotFoundError:
                self.rxdock_installed = False


        ### Check Smina ###
        if sys.platform == "win32":
            self.smina_installed = False

            if not check_only:
                self.programs_to_install["Smina"].installation_btn.setEnabled(False)
                self.programs_to_install["Smina"].installation_btn.setText("Smina is not available for this OS")

        else:
            path_to_smina = os.path.join(self.docking_programs.path_to_smina)
            try:
                result = subprocess.run([path_to_smina], stdout=subprocess.PIPE, stderr = subprocess.PIPE)
            except PermissionError:
                self.smina_installed = False

        ### Check Vina ###
        path_to_vina = os.path.join(self.docking_programs.path_to_vina)
        try:
            result = subprocess.run([path_to_vina], stdout=subprocess.PIPE, stderr = subprocess.PIPE)
        except PermissionError:
            self.vina_installed = False

        ### Check sdsorter ###

        if sys.platform == "win32":
            self.sdsorter_installed = False

            if not check_only:
                self.programs_to_install["sdsorter"].installation_btn.setEnabled(False)
                self.programs_to_install["sdsorter"].installation_btn.setText("sdsorter is not available for this OS")

        else:
            path_to_sdsorter = os.path.join(self.docking_programs.path_to_sdsorter)
            try:
                result = subprocess.run([path_to_sdsorter], stdout=subprocess.PIPE, stderr = subprocess.PIPE)
            except PermissionError:
                self.sdsorter_installed = False

        ### Check ADFR ###
        if sys.platform == "win32":

            a = shutil.which('adfr')
            if a is None:
                self.adfr_installed = False

        else:
            path_to_adfr = os.path.join(self.docking_programs.path_to_ADFR)

            try:
                result = subprocess.run([path_to_adfr], stdout=subprocess.PIPE, stderr = subprocess.PIPE)
            except PermissionError:
                self.adfr_installed = False


        ### Check spyrmsd ###
        try:
            from spyrmsd import io, rmsd
        except:
            self.spyrmsd_installed = False


        self.check_installation_dict = {}
        self.check_installation_dict["Openbabel"] = self.openbabel_installed
        self.check_installation_dict["RxDock"] = self.rxdock_installed
        self.check_installation_dict["Vina"] = self.vina_installed
        self.check_installation_dict["Smina"] = self.smina_installed
        self.check_installation_dict["ADFR"] = self.adfr_installed
        self.check_installation_dict["sPyRMSD"] = self.spyrmsd_installed
        self.check_installation_dict["sdsorter"] = self.sdsorter_installed


        if not check_only:
            # If the programs are installed the checkbuttons are disabled
            if self.rxdock_installed:
                self.programs_to_install["RxDock"].installation_btn.setEnabled(False)
                self.programs_to_install["RxDock"].installation_btn.setText("RxDock \nalready installed")
            else:
                if not sys.platform == "win32":
                    self.programs_to_install["RxDock"].installation_btn.setEnabled(True)
                    self.programs_to_install["RxDock"].installation_btn.setText("Install")

            if self.openbabel_installed:
                self.programs_to_install["Openbabel"].installation_btn.setEnabled(False)
                self.programs_to_install["Openbabel"].installation_btn.setText("Openbabel \nalready installed")
            else:
                self.programs_to_install["Openbabel"].installation_btn.setEnabled(True)
                self.programs_to_install["Openbabel"].installation_btn.setText("Install")

            if self.vina_installed:
                self.programs_to_install["Vina"].installation_btn.setEnabled(False)
                self.programs_to_install["Vina"].installation_btn.setText("Vina \nalready installed")
            else:
                self.programs_to_install["Vina"].installation_btn.setEnabled(True)
                self.programs_to_install["Vina"].installation_btn.setText("Install")

            if self.smina_installed:
                self.programs_to_install["Smina"].installation_btn.setEnabled(False)
                self.programs_to_install["Smina"].installation_btn.setText("Smina \nalready installed")
            else:
                if not sys.platform == "win32":
                    self.programs_to_install["Smina"].installation_btn.setEnabled(True)
                    self.programs_to_install["Smina"].installation_btn.setText("Install")

            if self.spyrmsd_installed:
                self.programs_to_install["sPyRMSD"].installation_btn.setEnabled(False)
                self.programs_to_install["sPyRMSD"].installation_btn.setText("sPyRMSD \nalready installed")
            else:
                self.programs_to_install["sPyRMSD"].installation_btn.setEnabled(True)
                self.programs_to_install["sPyRMSD"].installation_btn.setText("Install")

            if self.adfr_installed:
                self.programs_to_install["ADFR"].installation_btn.setEnabled(False)
                self.programs_to_install["ADFR"].installation_btn.setText("ADFR \nalready installed")
            else:
                self.programs_to_install["ADFR"].installation_btn.setEnabled(True)
                self.programs_to_install["ADFR"].installation_btn.setText("Install")

            if self.sdsorter_installed:
                self.programs_to_install["sdsorter"].installation_btn.setEnabled(False)
                self.programs_to_install["sdsorter"].installation_btn.setText("sdsorter \nalready installed")
            else:
                if not sys.platform == "win32":
                    self.programs_to_install["sdsorter"].installation_btn.setEnabled(True)
                    self.programs_to_install["sdsorter"].installation_btn.setText("Install")


        return self.check_installation_dict


class GridTab(QtWidgets.QWidget, PyMOLInteractions, HandleWidgets):

    """
    Layout for Grid Tab.
    It can be used for all those Docking Programs that set the Grid with the usage of
    center coordinates, spacing and x,y,z translations in space (e.g. Vina, Smina, Gnina)
    """

    def __init__(self, main_window, current_tab):
        super().__init__(main_window)
        self.docking_programs_child_tabs = main_window

        self.initUI(current_tab)


    def initUI(self, current_tab):

        # Grid Tab Layout
        self.layout_grid_tab = QtWidgets.QGridLayout()

        # Create Group Boxes
        self.grid_from_selection_group = QtWidgets.QGroupBox("Calculate Grid by selection")
        self.group_grid_center = QtWidgets.QGroupBox("Grid Center Coordinates")
        self.group_config_file = QtWidgets.QGroupBox("Open Config File")
        self.visualize_grid_group = QtWidgets.QGroupBox("Grid Dimension")

        # Add the group boxes to the Grid Tab Layout
        self.layout_grid_tab.addWidget(self.visualize_grid_group, 1, 0, 1, 2)
        self.layout_grid_tab.addWidget(self.grid_from_selection_group, 0, 0)
        self.layout_grid_tab.addWidget(self.group_grid_center, 2, 0, 1, 2)
        self.layout_grid_tab.addWidget(self.group_config_file, 0, 1)

        # GRID FROM SELECTION

        self.grid_selection_layout = QtWidgets.QVBoxLayout()
        self.grid_selection_layout.setAlignment(QtCore.Qt.AlignCenter)
        self.grid_from_selection_group.setLayout(self.grid_selection_layout)

        self.sele_label = QtWidgets.QLabel("Selection:")
        self.imported_sele = QtWidgets.QComboBox()
        self.imported_sele.view().pressed.connect(self.get_current_sele)
        self.imported_sele.currentTextChanged.connect(self.get_current_sele)

        self.import_sele = QtWidgets.QPushButton("Import Objects from PyMOL")
        self.import_sele.clicked.connect(self.import_objs_sele_func)

        self.set_as_reference_btn = QtWidgets.QPushButton("Use as a Reference")
        self.set_as_reference_btn.setToolTip('By clicking this Button\n    the Selection is ready to be used as a reference object to automatically generate the grid')
        self.set_as_reference_btn.clicked.connect(self.set_as_reference_func)

        self.grid_selection_layout.addStretch()
        self.grid_selection_layout.addStretch()
        self.grid_selection_layout.addWidget(self.import_sele)
        self.grid_selection_layout.addStretch()
        self.grid_selection_layout.addWidget(self.sele_label)
        self.grid_selection_layout.addWidget(self.imported_sele)
        self.grid_selection_layout.addStretch()
        self.grid_selection_layout.addWidget(self.set_as_reference_btn)
        self.grid_selection_layout.addStretch()

        # VISUALIZE GRID GROUP

        self.spacing, self.spacing_scroll_vis, self.x, self.x_scroll_vis, self.y, self.y_scroll_vis, self.z, self.z_scroll_vis = grid_dimensions(self.scroll_changed)

        self.visualize_layout = QtWidgets.QHBoxLayout()
        self.visualize_layout.setAlignment(QtCore.Qt.AlignCenter)
        self.visualize_grid_group.setLayout(self.visualize_layout)

        self.visualize_layout.addWidget(self.spacing)
        self.visualize_layout.addWidget(self.spacing_scroll_vis)
        self.visualize_layout.addStretch()

        self.visualize_layout.addWidget(self.x)
        self.visualize_layout.addWidget(self.x_scroll_vis)
        self.visualize_layout.addStretch()
        #
        self.visualize_layout.addWidget(self.y)
        self.visualize_layout.addWidget(self.y_scroll_vis)
        self.visualize_layout.addStretch()

        self.visualize_layout.addWidget(self.z)
        self.visualize_layout.addWidget(self.z_scroll_vis)
        self.visualize_layout.addStretch()

        #  GRID FROM CENTER COORDINATES

        self.grid_from_coord_layout = QtWidgets.QHBoxLayout()
        self.grid_from_coord_layout.setAlignment(QtCore.Qt.AlignCenter)
        self.group_grid_center.setLayout(self.grid_from_coord_layout)

        self.x, self.x_scroll, self.y, self.y_scroll, self.z, self.z_scroll = grid_position(self.scroll_changed)

        self.set_coords_btn = QtWidgets.QPushButton("Set in Docking Tab")
        self.set_coords_btn.clicked.connect(self.set_coords_func)

        self.grid_from_coord_layout.addStretch()
        self.grid_from_coord_layout.addWidget(self.x)
        self.grid_from_coord_layout.addWidget(self.x_scroll)
        self.grid_from_coord_layout.addStretch()
        #
        self.grid_from_coord_layout.addWidget(self.y)
        self.grid_from_coord_layout.addWidget(self.y_scroll)
        self.grid_from_coord_layout.addStretch()

        self.grid_from_coord_layout.addWidget(self.z)
        self.grid_from_coord_layout.addWidget(self.z_scroll)
        self.grid_from_coord_layout.addStretch()
        self.grid_from_coord_layout.addWidget(self.set_coords_btn)
        self.grid_from_coord_layout.addStretch()

        # GRID VINA CONFIG txt

        self.vina_layout = QtWidgets.QGridLayout()
        self.group_config_file.setLayout(self.vina_layout)
        self.open_vina_config = QtWidgets.QPushButton("Open .txt file")
        self.text_config_area = QtWidgets.QPlainTextEdit()
        self.opened_config_file_text = QtWidgets.QLabel("Opened files")
        self.opened_config_files = QtWidgets.QComboBox()

        text= " ---- Open Config File ----"
        self.text_config_area.setPlainText(text)
        self.text_config_area.setReadOnly(True)

        self.vina_layout.addWidget(self.open_vina_config, 0, 0, 1, 2)
        self.vina_layout.addWidget(self.text_config_area, 1, 0, 1, 2)
        self.vina_layout.addWidget(self.opened_config_files, 2, 0, 1, 2)

        self.open_vina_config.clicked.connect(lambda: self.open_vina_config_func(".txt"))
        self.opened_config_files.currentTextChanged.connect(self.on_config_changed)
        self.opened_config_files.view().pressed.connect(self.on_config_changed)

        if current_tab == "Vina":
            self.set_as_reference_btn.setEnabled(False)


    def get_current_sele(self):

        get_current_sele(self, self.docking_programs_child_tabs.docking_programs, self.imported_sele,
                         self.x_scroll,
                         self.y_scroll,
                         self.z_scroll,
                         self.x_scroll_vis,
                         self.y_scroll_vis,
                         self.z_scroll_vis,
                         self.spacing_scroll_vis
                         )

    def set_as_reference_func(self):

        Check_current_tab.check_docking_program_current_tab(self)

        reference = self.imported_sele.currentText()

        if reference:

            if self.is_smina_tab:
                self.docking_programs_child_tabs.docking_programs.SMINA.docking_tab_ui.loaded_cavities.addItem(reference)

            if self.is_adfr_tab:
                self.docking_programs_child_tabs.docking_programs.ADFR.docking_tab_ui.loaded_cavities.addItem(reference)

        self.docking_programs_child_tabs.docking_programs.main_window.statusBar().showMessage("Reference object added to Docking tab", 3000)


    def scroll_changed(self):

        # Docking Box is updated dinamically
        x = self.x_scroll.value()
        y = self.y_scroll.value()
        z = self.z_scroll.value()
        spacing = self.spacing_scroll_vis.value()
        x_vis = self.x_scroll_vis.value()
        y_vis = self.y_scroll_vis.value()
        z_vis = self.z_scroll_vis.value()

        spinbox = self.sender()

        if spinbox is self.x_scroll:
            x = self.x_scroll.value()
        elif spinbox is self.y_scroll:
            y = self.y_scroll.value()
        elif spinbox is self.z_scroll:
            z = self.z_scroll.value()
        elif spinbox is self.spacing_scroll_vis:
            spacing = self.spacing_scroll_vis.value()
        elif spinbox is self.x_scroll_vis:
            x_vis = self.x_scroll_vis.value()
        elif spinbox is self.y_scroll_vis:
            y_vis = self.y_scroll_vis.value()
        elif spinbox is self.z_scroll_vis:
            z_vis = self.z_scroll_vis.value()

        # Show the center and the box each time the values are updated
        self.show_crisscross_changed(x, y, z)
        self.calculate_box_changed(x, y, z, spacing, x_vis, y_vis, z_vis)


    def on_config_changed(self):

        # Get current Vina Config File
        self.current_file = self.opened_config_files.currentText()

        HandleWidgets.combobox_check_if_empty(self = self,
        widgets_list = [self.opened_config_files])

        if self.is_empty:
            pass

        else: # display current Config File parameters
            self.text_config_area.setPlainText(open(str(self.current_file + ".txt")).read())
            self.x_scroll.setValue(self.docking_programs_child_tabs.docking_programs.grid_center[self.current_file][0])
            self.y_scroll.setValue(self.docking_programs_child_tabs.docking_programs.grid_center[self.current_file][1])
            self.z_scroll.setValue(self.docking_programs_child_tabs.docking_programs.grid_center[self.current_file][2])
            self.x_scroll_vis.setValue(self.docking_programs_child_tabs.docking_programs.grid_center[self.current_file][3])
            self.y_scroll_vis.setValue(self.docking_programs_child_tabs.docking_programs.grid_center[self.current_file][4])
            self.z_scroll_vis.setValue(self.docking_programs_child_tabs.docking_programs.grid_center[self.current_file][5])

            #self.show_crisscross(self.current_file)


    def set_coords_func(self):

        Check_current_tab.check_docking_program_current_tab(self)

        # Take current info
        self.tmp_coord_list = [1, 1, 1, 1, 1, 1, 1]
        self.tmp_coord_list[0] = str(self.x_scroll.value())
        self.tmp_coord_list[1] = str(self.y_scroll.value())
        self.tmp_coord_list[2] = str(self.z_scroll.value())
        self.tmp_coord_list[3] = str(self.x_scroll_vis.value())
        self.tmp_coord_list[4] = str(self.y_scroll_vis.value())
        self.tmp_coord_list[5] = str(self.z_scroll_vis.value())
        self.tmp_coord_list[6] = str(self.spacing_scroll_vis.value())

        # Give a name to the cavity
        num = str(len(self.docking_programs_child_tabs.docking_programs.ready_grid_centers) + 1)
        name_cav = str("Grid Center_" + num)

        # Set the name as a key for the 'ready_grid_centers' DICTIONARY, where coordinates are stored
        self.docking_programs_child_tabs.docking_programs.ready_grid_centers[name_cav] = self.tmp_coord_list

        if self.is_vina_tab:
            current_program = self.docking_programs_child_tabs.docking_programs.VINA
        if self.is_adfr_tab:
            current_program = self.docking_programs_child_tabs.docking_programs.ADFR
        if self.is_smina_tab:
            current_program = self.docking_programs_child_tabs.docking_programs.SMINA

        # Add these info to the Docking Tab
        current_program.docking_tab_ui.loaded_cavities.addItem(name_cav)
        self.docking_programs_child_tabs.docking_programs.main_window.statusBar().showMessage("Grid parameters added in other tabs", 3000)
        self.docking_programs_child_tabs.docking_programs.main_window.statusBar().setSizeGripEnabled(1)


    def open_vina_config_func(self, extension):

        Check_current_tab.check_docking_program_current_tab(self)

        if self.is_vina_tab:
            current_program = self.docking_programs_child_tabs.docking_programs.VINA
            tmp_dir = self.docking_programs_child_tabs.docking_programs.vina_tmp_dir
        if self.is_adfr_tab:
            current_program = self.docking_programs_child_tabs.docking_programs.ADFR
            tmp_dir = self.docking_programs_child_tabs.docking_programs.adfr_tmp_dir
        if self.is_smina_tab:
            current_program = self.docking_programs_child_tabs.docking_programs.SMINA
            tmp_dir = self.docking_programs_child_tabs.docking_programs.smina_tmp_dir

        valid = True

        # Open a txt File
        opened_file = OpenFromFile(self,
        file_type_string = "Vina Config File (*.txt)")

        # Update the layout if the file is valid and if it is not already loaded
        if opened_file.valid_config and opened_file.is_already_loaded == False:

            # Change dir --> Vina tmp_dir
            os.chdir(tmp_dir)

            # Read the info in Config File
            self.tmp_coord_dict = {"center_x": 1, "center_y": 1, "center_z":1, "size_x":1, "size_y":1, "size_z":1, "exhaustiveness":8}

            self.tmp_coord_list = [1, 1, 1, 1, 1, 1, 1]

            input = open(str(opened_file.file_name + ".txt"), "rt")

            for line in input:
                for parameter in opened_file.config_parameters_list:
                    if line.startswith(parameter):
                        value = re.findall('[-+]?([0-9]*\.[0-9]+|[0-9]+)', line)
                        self.tmp_coord_dict[parameter] = float(value[0])

            input.close()

            self.tmp_coord_list[6] = float(self.spacing_scroll_vis.value())

            values = self.tmp_coord_dict.values()
            self.tmp_coord_list = list(values)

            # Display text
            self.text_config_area.setPlainText(open(str(opened_file.file_name + ".txt")).read())

            # Update the grid_center_dictionary
            self.docking_programs_child_tabs.docking_programs.grid_center[str(opened_file.file_name)] = self.tmp_coord_list

            # Set the info
            self.x_scroll.setValue(self.tmp_coord_dict["center_x"])
            self.y_scroll.setValue(self.tmp_coord_dict["center_y"])
            self.z_scroll.setValue(self.tmp_coord_dict["center_z"])
            self.x_scroll_vis.setValue(self.tmp_coord_dict["size_x"])
            self.y_scroll_vis.setValue(self.tmp_coord_dict["size_y"])
            self.z_scroll_vis.setValue(self.tmp_coord_dict["size_z"])

            current_program.docking_tab_ui.exhaustiveness_box.setValue(self.tmp_coord_dict["exhaustiveness"])

            # Update ComboBox of opened files
            HandleWidgets.combobox_check_existing_item(self = self, widget = self.opened_config_files, item_name = opened_file.file_name)

            # Show grid_center, box and wirebox
            #self.show_crisscross(self.opened_config_files.currentText())

            self.docking_programs_child_tabs.docking_programs.main_window.statusBar().showMessage(str("Loaded file '" + opened_file.file_name + "'"), 3000)


    def import_objs_sele_func(self):

        update_widget_with_pymol_object(self.docking_programs_child_tabs, self.imported_sele, selections = True)

        # # Create a list of all the importable Objects and selections from PyMOL
        # self.selections = [str(obj) for obj in cmd.get_names("objects") + cmd.get_names("selections")]
        #
        # # Check for importable objects
        # if self.selections == []:
        #     QtWidgets.QMessageBox.warning(self.docking_programs_child_tabs, "PyMOL is empty", str("There isn't any object to import"))
        #
        # # If importable objects are present, update the ComboBox
        # else:
        #     self.imported_sele.clear()
        #     for i in self.selections:
        #         type = cmd.get_type(i)
        #         if type == str("object:molecule"):
        #             if re.search("Run_", i):
        #                 pass
        #             else:
        #                 self.imported_sele.addItem(i)
        #         elif type == str("selection"):
        #             self.imported_sele.addItem(i)
        #         else:
        #             pass


class GridTab_RxDock(QtWidgets.QWidget, PyMOLInteractions):

    """
    Layout for Grid Tab.
    It is specific for RxDock
    """

    def __init__(self, main_window):
        super().__init__(main_window)
        self.style = "background-color: rgb(0, 0, 0); color: rgb(255, 255, 255); font-weight: bold"
        self.docking_programs_child_tabs = main_window
        self.initUI()

    def initUI(self):

        # Grid Tab Layout
        self.layout_grid_tab_rxdock = QtWidgets.QGridLayout()

        # Create Group Boxes
        self.grid_from_file_group = QtWidgets.QGroupBox("Open Grid from file")
        self.grid_reference_group = QtWidgets.QGroupBox("Calculate Grid")
        self.visualize_grid_group = QtWidgets.QGroupBox("Visualize")

        # Add the group boxes to the Grid Tab Layout
        self.layout_grid_tab_rxdock.addWidget(self.grid_from_file_group, 0, 0)
        self.layout_grid_tab_rxdock.addWidget(self.grid_reference_group, 1, 0)
        self.layout_grid_tab_rxdock.addWidget(self.visualize_grid_group, 2, 0)

        # GRID FROM FILE

        self.grid_from_file_layout = QtWidgets.QGridLayout()
        self.grid_from_file_layout.setAlignment(QtCore.Qt.AlignLeft)
        self.grid_file_btn = QtWidgets.QPushButton("Open from file")
        self.grid_from_file_group.setLayout(self.grid_from_file_layout)
        self.grid_from_file_layout.addWidget(self.grid_file_btn, 0, 0)
        self.grid_file_btn.clicked.connect(lambda: self.open_grid_file(".grd"))

        # GRID WITH A REFERENCE OBJECT

        self.grid_reference_layout = QtWidgets.QGridLayout()
        self.grid_reference_group.setLayout(self.grid_reference_layout)

        self.two_spheres_radiobtn = QtWidgets.QRadioButton("Two Spheres Method")
        self.reference_radiobtn = QtWidgets.QRadioButton("Reference ligand method")
        self.reference_radiobtn.setChecked(True)
        self.two_spheres_radiobtn.clicked.connect(self.change_text)
        self.reference_radiobtn.clicked.connect(self.change_text)

        self.grid_reference_layout.addWidget(self.two_spheres_radiobtn, 0, 0)
        self.grid_reference_layout.addWidget(self.reference_radiobtn, 0, 1)

        self.ready_receptors_list = QtWidgets.QComboBox()
        self.ready_receptors_list_text = QtWidgets.QLabel("Receptor")
        self.ready_receptors_list_text.setBuddy(self.ready_receptors_list)
        self.grid_reference_layout.addWidget(self.ready_receptors_list_text, 1, 0)
        self.grid_reference_layout.addWidget(self.ready_receptors_list, 2, 0)
        self.ready_receptors_list.currentTextChanged.connect(lambda: HandleWidgets.combobox_orient_current_item_pymol(self = self, widget = self.ready_receptors_list))
        self.ready_receptors_list.view().pressed.connect(lambda: HandleWidgets.combobox_orient_current_item_pymol(self = self, widget = self.ready_receptors_list))

        self.ready_ligands_list = QtWidgets.QComboBox()
        self.ready_ligands_list_text = QtWidgets.QLabel("Reference Ligand")
        self.grid_reference_layout.addWidget(self.ready_ligands_list_text, 1, 1)
        self.grid_reference_layout.addWidget(self.ready_ligands_list, 2, 1)
        self.ready_ligands_list.currentTextChanged.connect(lambda: HandleWidgets.combobox_orient_current_item_pymol(self = self, widget = self.ready_ligands_list))
        self.ready_ligands_list.view().pressed.connect(lambda: HandleWidgets.combobox_orient_current_item_pymol(self = self, widget = self.ready_ligands_list))

        self.ready_ligands_list.currentTextChanged.connect(self.get_reference_coord)
        self.ready_ligands_list.view().pressed.connect(self.get_reference_coord)


        self.generate_cavity_btn = QtWidgets.QPushButton("Generate Cavity")
        self.generate_cavity_btn.setToolTip("Generate the grid that will be used during Docking \nusing the selected Receptor and Reference Ligand")
        self.generate_cavity_btn.clicked.connect(self.generate_cavity_func)
        self.grid_reference_layout.addWidget(self.generate_cavity_btn, 3, 0, 1, 2)

        self.radius_label = QtWidgets.QLabel("        Radius:")
        self.radius_label.setEnabled(False)
        self.radius_spinbox = QtWidgets.QDoubleSpinBox()
        self.radius_spinbox.setEnabled(False)
        self.radius_spinbox.setRange(10.00, 20.00)
        self.radius_spinbox.setSingleStep(1)
        self.radius_spinbox.setValue(10.00)

        self.small_sphere_label = QtWidgets.QLabel("        Small Sphere:")
        self.small_sphere_label.setEnabled(False)
        self.small_sphere_spinbox = QtWidgets.QDoubleSpinBox()
        self.small_sphere_spinbox.setEnabled(False)
        self.small_sphere_spinbox.setRange(1.00, 2.00)
        self.small_sphere_spinbox.setSingleStep(00.01)
        self.small_sphere_spinbox.setValue(1.50)

        self.large_sphere_label = QtWidgets.QLabel("        Large Sphere:")
        self.large_sphere_label.setEnabled(False)
        self.large_sphere_spinbox = QtWidgets.QDoubleSpinBox()
        self.large_sphere_spinbox.setEnabled(False)
        self.large_sphere_spinbox.setRange(3.50, 6.00)
        self.large_sphere_spinbox.setSingleStep(00.01)
        self.large_sphere_spinbox.setValue(4.00)

        self.grid_reference_layout.addWidget(self.radius_label, 0, 2)
        self.grid_reference_layout.addWidget(self.radius_spinbox, 0, 3)
        self.grid_reference_layout.addWidget(self.small_sphere_label, 1, 2)
        self.grid_reference_layout.addWidget(self.small_sphere_spinbox, 1, 3)
        self.grid_reference_layout.addWidget(self.large_sphere_label, 2, 2)
        self.grid_reference_layout.addWidget(self.large_sphere_spinbox, 2, 3)

        self.coord_frame = OptionsFrame(parent=None,
        main_window=self.docking_programs_child_tabs)
        self.coord_frame.options_frame_layout.setAlignment(QtCore.Qt.AlignCenter)

        self.grid_reference_layout.addWidget(self.coord_frame, 4, 0, 1, 4)

        self.x = QtWidgets.QLabel("X:")
        self.x_scroll = QtWidgets.QDoubleSpinBox()
        self.x.setEnabled(False)
        self.x_scroll.setEnabled(False)
        self.x_scroll.setRange(-1000, 1000)
        self.x_scroll.setSingleStep(00.10)
        self.x_scroll.valueChanged.connect(self.scroll_changed)

        self.y = QtWidgets.QLabel("             Y:")
        self.y_scroll = QtWidgets.QDoubleSpinBox()
        self.y.setEnabled(False)
        self.y_scroll.setEnabled(False)
        self.y_scroll.setRange(-1000, 1000)
        self.y_scroll.setSingleStep(00.10)
        self.y_scroll.valueChanged.connect(self.scroll_changed)

        self.z = QtWidgets.QLabel("             Z:")
        self.z_scroll = QtWidgets.QDoubleSpinBox()
        self.z.setEnabled(False)
        self.z_scroll.setEnabled(False)
        self.z.setBuddy(self.z_scroll)
        self.z_scroll.setRange(-1000, 1000)
        self.z_scroll.setSingleStep(00.10)
        self.z_scroll.valueChanged.connect(self.scroll_changed)

        self.coord_frame.options_frame_layout.addWidget(self.x, 4, 0)
        self.coord_frame.options_frame_layout.addWidget(self.x_scroll, 4, 1)

        self.coord_frame.options_frame_layout.addWidget(self.y, 4, 2)
        self.coord_frame.options_frame_layout.addWidget(self.y_scroll, 4, 3)

        self.coord_frame.options_frame_layout.addWidget(self.z, 4, 4)
        self.coord_frame.options_frame_layout.addWidget(self.z_scroll, 4, 5)

        # VISUALIZE and save
        self.visualize_layout = QtWidgets.QGridLayout()
        self.visualize_grid_group.setLayout(self.visualize_layout)

        # self.grid_file_btn = QtWidgets.QPushButton("Open from file")
        # self.grid_file_btn.clicked.connect(lambda: self.open_grid_file(".grd"))
        # self.visualize_layout.addWidget(self.grid_file_btn, 0, 0)

        self.cavity_listwidget = QtWidgets.QComboBox()
        self.visualize_layout.addWidget(self.cavity_listwidget, 0, 0)
        self.cavity_listwidget.currentTextChanged.connect(lambda: HandleWidgets.combobox_orient_current_item_pymol(self = self, widget = self.cavity_listwidget))
        self.cavity_listwidget.view().pressed.connect(lambda: HandleWidgets.combobox_orient_current_item_pymol(self = self, widget = self.cavity_listwidget))

        self.set_cav_btn = QtWidgets.QPushButton("Set")
        self.set_cav_btn.clicked.connect(self.add_to_other_tabs)
        self.visualize_layout.addWidget(self.set_cav_btn, 1, 0)

        self.export_btn = QtWidgets.QPushButton("Save to file")
        self.visualize_layout.addWidget(self.export_btn, 1, 1)
        self.export_btn.clicked.connect(self.save_grd_file_func)


    def get_reference_coord(self):

        try:

            stored.xyz = []
            cmd.iterate_state(1, self.ready_ligands_list.currentText(),"stored.xyz.append([x,y,z])")
            xx = statistics.mean(map(lambda a: a[0], stored.xyz))
            yy = statistics.mean(map(lambda a: a[1], stored.xyz))
            zz = statistics.mean(map(lambda a: a[2], stored.xyz))
            x = str(round(xx,2))
            y = str(round(yy,2))
            z = str(round(zz,2))
            coords = str("(" + x + "," + y + "," + z + ")")

            self.x_scroll.setValue(float(x))
            self.y_scroll.setValue(float(y))
            self.z_scroll.setValue(float(z))

        except:
            pass


    def scroll_changed(self):

        # Docking Box is updated dinamically
        x = self.x_scroll.value()
        y = self.y_scroll.value()
        z = self.z_scroll.value()

        spinbox = self.sender()

        if spinbox is self.x_scroll:
            x = self.x_scroll.value()
        elif spinbox is self.y_scroll:
            y = self.y_scroll.value()
        elif spinbox is self.z_scroll:
            z = self.z_scroll.value()

        # Show the center and the box each time the values are updated
        self.show_crisscross_changed(x, y, z)



    def save_grd_file_func(self):


        # Save the grd file
        if self.cavity_listwidget.currentText():

            filepath = asksaveasfile_qt("Save Grid file", name_filter="*.grd")

            if not filepath:
                return None

            else:
                os.chdir(self.docking_programs_child_tabs.docking_programs.rxdock_tmp_dir)
                shutil.copyfile(str(self.cavity_listwidget.currentText() + ".grd"), filepath)

            # Save the as file
            filepath = asksaveasfile_qt("Save as file", name_filter="*.as")

            if not filepath:
                return None

            else:
                os.chdir(self.docking_programs_child_tabs.docking_programs.rxdock_tmp_dir)
                name = self.cavity_listwidget.currentText().replace("_cav1", "")
                shutil.copyfile(str(name + ".as"), filepath)



    def change_text(self):

        prepared_ligands = self.docking_programs_child_tabs.docking_programs.rxdock_ligands_prepared

        if self.two_spheres_radiobtn.isChecked():
            self.update_with_pymol_objects(self.ready_ligands_list)
            self.ready_ligands_list_text.setText("Reference to extract coordinates")
            self.radius_label.setEnabled(True)
            self.radius_spinbox.setEnabled(True)
            self.small_sphere_label.setEnabled(True)
            self.small_sphere_spinbox.setEnabled(True)
            self.large_sphere_label.setEnabled(True)
            self.large_sphere_spinbox.setEnabled(True)
            self.x.setEnabled(True)
            self.x_scroll.setEnabled(True)
            self.y.setEnabled(True)
            self.y_scroll.setEnabled(True)
            self.z.setEnabled(True)
            self.z_scroll.setEnabled(True)


        if self.reference_radiobtn.isChecked():
            self.ready_ligands_list_text.setText("Reference Ligand")
            self.ready_ligands_list.clear()
            self.ready_ligands_list.addItems(prepared_ligands)
            self.radius_label.setEnabled(False)
            self.radius_spinbox.setEnabled(False)
            self.small_sphere_label.setEnabled(False)
            self.small_sphere_spinbox.setEnabled(False)
            self.large_sphere_label.setEnabled(False)
            self.large_sphere_spinbox.setEnabled(False)
            self.x.setEnabled(False)
            self.x_scroll.setEnabled(False)
            self.y.setEnabled(False)
            self.y_scroll.setEnabled(False)
            self.z.setEnabled(False)
            self.z_scroll.setEnabled(False)



    def update_with_pymol_objects(self, widget):

        # Create a list of all the importable Objects and selections from PyMOL
        self.selections = [str(obj) for obj in cmd.get_names("objects") + cmd.get_names("selections")]

        # Check for importable objects
        if self.selections == []:
            pass
        # If importable objects are present, update the ComboBox
        else:
            widget.clear()
            for i in self.selections:
                type = cmd.get_type(i)
                if type == str("object:molecule"):
                    if re.search("Run_", i):
                        pass
                    else:
                        widget.addItem(i)
                elif type == str("selection"):
                    widget.addItem(i)
                else:
                    pass


    def add_to_other_tabs(self):

        HandleWidgets.combobox_check_if_empty(self = self,
        widgets_list = [self.cavity_listwidget])

        if self.is_empty:
            pass
        else:
            HandleWidgets.add_to_other_tab(self = self, input_widget = self.cavity_listwidget, output_widget = self.docking_programs_child_tabs.docking_programs.RXDOCK.docking_tab_ui.loaded_cavities)


    def generate_cavity_func(self):

        self.two_spheres_method = False
        self.reference_ligand = False

        widget = self.cavity_listwidget

        HandleWidgets.combobox_check_if_empty(self=self,
        widgets_list=[self.ready_ligands_list, self.ready_receptors_list])

        if self.is_empty:
            QtWidgets.QMessageBox.warning(self.docking_programs_child_tabs.docking_programs.main_window, "Warning", "There aren't enough parameters")

        else:
            if self.two_spheres_radiobtn.isChecked():
                self.two_spheres_method = True
            else:
                self.reference_ligand = True

            cavity = RxDock_Cavity(self,
                                   main = self.docking_programs_child_tabs.docking_programs,
                                   reference_receptor = self.ready_receptors_list.currentText(),
                                   reference_ligand = self.ready_ligands_list.currentText(),
                                   two_spheres_method = self.two_spheres_method,
                                   reference_ligand_method = self.reference_ligand,
                                   radius_spinbox = self.radius_spinbox.value(),
                                   small_sphere_value = self.small_sphere_spinbox.value(),
                                   large_sphere_value = self.large_sphere_spinbox.value()
                                   )

            # The cavity generated by RxDock takes the name by the prm_file name + "_cav1" in ".grd" format
            self.file_to_load = str(cavity.prm_file.prm_file_name + "_cav1.grd")
            self.object_name = str(cavity.prm_file.prm_file_name + "_cav1")

            # To check if PyMOL can find the file
            try:
                RxDock_Functions.import_grid_in_pymol(self, self.file_to_load, self.object_name)

                # update the list of generated cavities
                self.docking_programs_child_tabs.docking_programs.generated_cavity.append(cavity.prm_file.prm_file_name)

                # check for file with the same name already loaded and add it to the combobox
                HandleWidgets.combobox_check_existing_item(self = self, widget = self.cavity_listwidget, item_name = self.object_name)
                self.cavity_listwidget.setCurrentText(self.object_name)

            except:
                QtWidgets.QMessageBox.warning(self.docking_programs_child_tabs, "Grid computing failed", str("RxDock was not able to compute grid with the input parameters.\n\nHints:\n- check if the specified coordinates are suited for RxDock to find a cavity in that position (the 'grid_center' object in PyMOL could assist you)"))

        try:
            self.docking_programs_child_tabs.docking_programs.main_window.statusBar().showMessage(str("Generated Cavity '" + self.object_name + "'"), 3000)
        except:
            pass


    def show_grid_in_pymol_func(self, widget):

        self.grid_name = widget.currentText()
        self.grid_name_ext = str(self.grid_name + ".grd")

        os.chdir(self.docking_programs_child_tabs.docking_programs.rxdock_tmp_dir)

        RxDock_Functions.import_grid_in_pymol(self = self, file_to_load = self.grid_name_ext, pymol_object_name = self.grid_name)

        cmd.zoom(self.grid_name)


    def open_grid_file(self, extension):

        # warn user that should be loaded two different types of files
        QtWidgets.QMessageBox.warning(self.docking_programs_child_tabs, "", str("You will be asked to load two files: \".as\" and \".grd\" files. \nThese two formats are needed by RxDock to carry out the docking process."))

        # open ".grd" file
        opened_file = OpenFromFile(self,
        file_type_string = "Grid File (*.grd)")

        if opened_file.is_valid and opened_file.is_already_loaded == False:

            RxDock_Functions.import_grid_in_pymol(self = self, file_to_load = opened_file.new_file_path, pymol_object_name = opened_file.file_name)

            HandleWidgets.combobox_check_existing_item(self = self, widget = self.cavity_listwidget, item_name = opened_file.file_name)

            self.docking_programs_child_tabs.docking_programs.main_window.statusBar().showMessage(str("Loaded file '" + opened_file.file_name + "'"), 3000)

        # open ".as" file
        as_opened_file = OpenFromFile(self,
        file_type_string = "Grid File (*.as)")

        if as_opened_file.is_valid and as_opened_file.is_already_loaded == False:

            self.docking_programs_child_tabs.docking_programs.main_window.statusBar().showMessage(str("Loaded file '" + as_opened_file.file_name + "'"), 3000)



class ReceptorTab(QtWidgets.QWidget, PyMOLInteractions, HandleWidgets, Import_from_Pymol):

    def __init__(self, main_window, current_tab):
        super().__init__(main_window)
        self.docking_programs_child_tabs = main_window

        self.initUI(current_tab)

    def initUI(self, current_tab):

        # Receptor Tab Layout
        self.layout_receptor_tab = QtWidgets.QVBoxLayout()

        # Create Group Boxes
        self.load_receptors_group = QtWidgets.QGroupBox("Receptors")
        self.prepare_receptors_group = QtWidgets.QGroupBox("Receptors settings")

        # Add the Group Boxes to the Receptor Tab Layout
        self.layout_receptor_tab.addWidget(self.load_receptors_group)
        self.layout_receptor_tab.addWidget(self.prepare_receptors_group)

        # Create and set Group Boxes' layouts
        self.load_receptors_layout = QtWidgets.QHBoxLayout()
        self.prepare_receptors_layout = QtWidgets.QHBoxLayout()

        self.load_receptors_group.setLayout(self.load_receptors_layout)
        self.prepare_receptors_group.setLayout(self.prepare_receptors_layout)

        # Create widgets for "Load receptors group Box"
        self.import_btn = QtWidgets.QPushButton("Import from PyMOL")

        self.load_receptors_layout.addWidget(self.import_btn)
        self.load_receptors_layout.addStretch()

        # Create widgets for "Prepare Receptors Group"
        self.main_page_scroll = QtWidgets.QScrollArea()
        self.main_page_interior = QtWidgets.QWidget()
        self.main_page_scroll.setWidgetResizable(True)
        self.main_page_scroll.setWidget(self.main_page_interior)

        self.middle_left_layout = QtWidgets.QGridLayout()
        self.prepare_receptors_layout.addLayout(self.middle_left_layout)
        self.middle_left_layout.addWidget(self.main_page_scroll, 1, 0)

        self.prepare_structures_frame_layout = QtWidgets.QFormLayout()
        self.main_page_interior.setLayout(self.prepare_structures_frame_layout)

        # Add pdbqt options
        if current_tab != "RxDock":
            self.pdbqt_options_dict = {}
            self.pdbqt_options_dict["add_h"] = True
            self.pdbqt_options_dict["bonds"] = False
            self.pdbqt_options_dict["add_gast"] = False
            self.pdbqt_options_dict["remove_nonstd"] = False
            self.pdbqt_options_dict["remove_water"] = True
            self.pdbqt_options_dict["remove_lone_pairs"] = False
            self.pdbqt_options_dict["remove_non_polar_H"] = False
            self.pdbqt_options_dict["remove_non_protein"] = False

            self.pdbqt_options_pushbutton = QtWidgets.QPushButton("PDBQT Advanced Options")
            self.prepare_structures_frame_layout.addRow(self.pdbqt_options_pushbutton)
            self.pdbqt_options_pushbutton.clicked.connect(self.show_pdbqt_options_window)

        # Add Select All Option
        self.select_all_btn = QtWidgets.QCheckBox("All")
        self.prepare_structures_frame_layout.addRow(self.select_all_btn)
        self.select_all_btn.clicked.connect(self.get_list_of_current_cb)

        self.middle_btn_layout = QtWidgets.QVBoxLayout()
        self.prepare_receptors_layout.addLayout(self.middle_btn_layout)
        self.rec_text_set2 = QtWidgets.QPushButton("Generate Receptor")

        self.middle_btn_layout.addStretch()
        self.middle_btn_layout.addWidget(self.rec_text_set2)
        self.middle_btn_layout.addStretch()

        self.middle_right_layout = QtWidgets.QVBoxLayout()
        self.prepare_receptors_layout.addLayout(self.middle_right_layout)

        self.listwidget = QtWidgets.QListWidget()
        self.listwidget.setFixedWidth(70)
        self.listwidget.setSelectionMode(3)
        self.middle_right_layout.addWidget(self.listwidget)

        self.set = QtWidgets.QPushButton("Set")
        self.remove = QtWidgets.QPushButton("Remove")

        self.middle_right_layout.addWidget(self.set)
        self.middle_right_layout.addWidget(self.remove)

        #self.remove.hide()

        # Functions
        self.rec_text_set2.clicked.connect(self.generate_receptor_func)
        self.import_btn.clicked.connect(self.import_rec_from_pymol_current_tab)
        self.set.clicked.connect(self.add_to_other_tabs)
        self.remove.clicked.connect(self.remove_rec_from_list)


    def get_list_of_current_cb(self):

        self.list_of_cb = []

        for i in reversed(range(self.prepare_structures_frame_layout.count())):
            obj = self.prepare_structures_frame_layout.itemAt(i).widget()
            if str(type(obj)) == "<class 'lib.docking_program_main.docking_program_gui.frames.StructuresFrame'>":
                self.list_of_cb.append(obj.strc_checkbox)

        SelectAll(self, all_cb = self.select_all_btn,
        list_of_cb = self.list_of_cb)


    def import_rec_from_pymol_current_tab(self):

        a = self.parent().docking_programs.docking_programs_tabs.currentIndex()
        current_tab = self.parent().docking_programs.docking_programs_tabs.tabText(a)

        Import_from_Pymol(self, current_tab, is_receptor = True)


    def add_to_other_tabs(self):

        Check_current_tab.check_docking_program_current_tab(self)

        if self.is_vina_tab:
            HandleWidgets.add_to_other_tab(self = self,
            input_widget = self.listwidget,
            output_widget = self.docking_programs_child_tabs.docking_programs.VINA.docking_tab_ui.ready_receptors_list)

        if self.is_rxdock_tab:
            HandleWidgets.add_to_other_tab(self = self,
            input_widget = self.listwidget,
            output_widget = self.docking_programs_child_tabs.docking_programs.RXDOCK.docking_tab_ui.ready_receptors_list)

            HandleWidgets.add_to_other_tab(self = self,
            input_widget = self.listwidget,
            output_widget = self.docking_programs_child_tabs.docking_programs.RXDOCK.grid_tab_rxdock_ui.ready_receptors_list)

        if self.is_smina_tab:
            HandleWidgets.add_to_other_tab(self = self,
            input_widget = self.listwidget,
            output_widget = self.docking_programs_child_tabs.docking_programs.SMINA.docking_tab_ui.ready_receptors_list)

        if self.is_adfr_tab:
            HandleWidgets.add_to_other_tab(self = self,
            input_widget = self.listwidget,
            output_widget = self.docking_programs_child_tabs.docking_programs.ADFR.docking_tab_ui.ready_receptors_list)

        if self.listwidget.count() == 0:
            pass
        else:
           self.docking_programs_child_tabs.docking_programs.main_window.statusBar().showMessage("Selected Items added to other tabs", 3000)


    def show_pdbqt_options_window(self):

        pdbqt_options_window = PDBQT_OptionsWindows(tab = self, main = self.docking_programs_child_tabs,
                             obj_type = "receptor", options_dict = self.pdbqt_options_dict)

        # self.pdbqt_options_window = NewWindow(parent = self.docking_programs_child_tabs,
        # title = "PDBQT options window", upper_frame_title = "Select Options",
        # submit_command = self.apply_pdbqt_options, submit_button_text= "Set",
        # with_scroll = True)
        #
        # self.add_h = QtWidgets.QCheckBox("Add Hydrogens")
        # self.bonds = QtWidgets.QCheckBox("Repair Bonds")
        # self.add_gast = QtWidgets.QCheckBox("Add Gasteiger Charges")
        # self.remove_nonstd = QtWidgets.QCheckBox("Remove ALL non-standard residues")
        # self.remove_water = QtWidgets.QCheckBox("Remove Water")
        # self.remove_lone_pairs = QtWidgets.QCheckBox("Remove Lone Pairs")
        # self.remove_non_polar_H = QtWidgets.QCheckBox("Remove Non-Polar Hydrogens")
        # self.remove_non_protein = QtWidgets.QCheckBox("Remove All Non-Protein chains")
        #
        # self.pdbqt_options_window.middle_layout_type.addWidget(self.add_h)
        # self.pdbqt_options_window.middle_layout_type.addWidget(self.bonds)
        # self.pdbqt_options_window.middle_layout_type.addWidget(self.add_gast)
        # self.pdbqt_options_window.middle_layout_type.addWidget(self.remove_nonstd)
        # self.pdbqt_options_window.middle_layout_type.addWidget(self.remove_water)
        # self.pdbqt_options_window.middle_layout_type.addWidget(self.remove_lone_pairs)
        # self.pdbqt_options_window.middle_layout_type.addWidget(self.remove_non_polar_H)
        # self.pdbqt_options_window.middle_layout_type.addWidget(self.remove_non_protein)
        #
        # self.add_h.setChecked(self.pdbqt_options_dict["add_h"])
        # self.bonds.setChecked(self.pdbqt_options_dict["bonds"])
        # self.bonds.setChecked(self.pdbqt_options_dict["add_gast"])
        # self.remove_nonstd.setChecked(self.pdbqt_options_dict["remove_nonstd"])
        # self.remove_water.setChecked(self.pdbqt_options_dict["remove_water"])
        # self.remove_lone_pairs.setChecked(self.pdbqt_options_dict["remove_lone_pairs"])
        # self.remove_non_polar_H.setChecked(self.pdbqt_options_dict["remove_non_polar_H"])
        # self.remove_non_polar_H.setChecked(self.pdbqt_options_dict["remove_non_polar_H"])
        # self.remove_non_protein.setChecked(self.pdbqt_options_dict["remove_non_protein"])
        #
        # self.pdbqt_options_window.show()


    def apply_pdbqt_options(self):

        self.pdbqt_options_dict["add_h"] = self.add_h.isChecked()
        self.pdbqt_options_dict["bonds"] = self.bonds.isChecked()
        self.pdbqt_options_dict["remove_nonstd"] = self.remove_nonstd.isChecked()
        self.pdbqt_options_dict["remove_water"] = self.remove_water.isChecked()
        self.pdbqt_options_dict["remove_lone_pairs"] = self.remove_lone_pairs.isChecked()
        self.pdbqt_options_dict["remove_non_polar_H"] = self.remove_non_polar_H.isChecked()
        self.pdbqt_options_dict["remove_non_protein"] = self.remove_non_protein.isChecked()

        self.pdbqt_options_window.close()


    def generate_receptor_func(self):

        if not check_configuration(self, self.docking_programs_child_tabs.docking_programs):
            QtWidgets.QMessageBox.warning(self.docking_programs_child_tabs, "Warning", "Please check DockingPie configuration!")

        else:
            Check_current_tab.check_docking_program_current_tab(self)

            if self.is_vina_tab:
                self.generated_receptor = Generate_Object(self, main = self.docking_programs_child_tabs.docking_programs,
                dict = self.docking_programs_child_tabs.docking_programs.vina_receptors_dict,
                prepared_objects_list = self.docking_programs_child_tabs.docking_programs.vina_receptors_prepared,
                tmp_path = self.docking_programs_child_tabs.docking_programs.vina_tmp_dir,
                format = "pdb",
                docking_program_name = "Vina",
                generate_pdbqt = True,
                pdbqt_dict = self.pdbqt_options_dict,
                is_receptor = True)

                self.generated_object_update_gui()
                cmd.group("Vina", members=self.generated_receptor.new_strc_name, action='auto', quiet=1)

            if self.is_rxdock_tab:
                self.generated_receptor = Generate_Object(self, main = self.docking_programs_child_tabs.docking_programs,
                dict = self.docking_programs_child_tabs.docking_programs.rxdock_receptors_dict,
                prepared_objects_list = self.docking_programs_child_tabs.docking_programs.rxdock_receptors_prepared,
                tmp_path = self.docking_programs_child_tabs.docking_programs.rxdock_tmp_dir,
                format = "mol2",
                docking_program_name = "RxDock",
                generate_pdbqt = False,
                is_receptor = True)

                self.generated_object_update_gui()
                cmd.group("RxDock", members=self.generated_receptor.new_strc_name, action='auto', quiet=1)

            if self.is_smina_tab:
                self.generated_receptor = Generate_Object(self, main = self.docking_programs_child_tabs.docking_programs,
                dict = self.docking_programs_child_tabs.docking_programs.smina_receptors_dict,
                prepared_objects_list = self.docking_programs_child_tabs.docking_programs.smina_receptors_prepared,
                tmp_path = self.docking_programs_child_tabs.docking_programs.smina_tmp_dir,
                format = "pdb",
                docking_program_name = "Smina",
                generate_pdbqt = True,
                pdbqt_dict = self.pdbqt_options_dict,
                is_receptor = True)

                self.generated_object_update_gui()
                cmd.group("Smina", members=self.generated_receptor.new_strc_name, action='auto', quiet=1)

            if self.is_adfr_tab:
                self.generated_receptor = Generate_Object(self, main = self.docking_programs_child_tabs.docking_programs,
                dict = self.docking_programs_child_tabs.docking_programs.adfr_receptors_dict,
                prepared_objects_list = self.docking_programs_child_tabs.docking_programs.adfr_receptors_prepared,
                tmp_path = self.docking_programs_child_tabs.docking_programs.adfr_tmp_dir,
                format = "pdb",
                docking_program_name = "ADFR",
                generate_pdbqt = True,
                pdbqt_dict = self.pdbqt_options_dict,
                is_receptor = True)

                self.generated_object_update_gui()
                cmd.group("ADFR", members=self.generated_receptor.new_strc_name, action='auto', quiet=1)


    def generated_object_update_gui(self):

        if os.path.isfile(self.generated_receptor.new_strc_name_format):

            # Add to the listwidget
            self.listbtn = self.listwidget.addItem(self.generated_receptor.new_strc_name)

            # Add to the list of prepared files
            self.generated_receptor.prepared_objects_list.append(self.generated_receptor.new_strc_name)

            # Load in PyMOL
            cmd.load(self.generated_receptor.new_receptor_file_path, self.generated_receptor.new_strc_name)

            self.docking_programs_child_tabs.docking_programs.main_window.statusBar().showMessage(str("Generated Receptor '" + self.generated_receptor.new_strc_name + "'"), 3000)

        else:
            if self.generated_receptor.generate_pdbqt:

                self.about_text_area = QtWidgets.QPlainTextEdit()
                self.about_text_area.setReadOnly(True)

                self.about_text_area.setPlainText(open(self.pdbqt_log_file_name).read())

                self.pdbqt_log = NewWindow(parent = self.tab,
                title = "Warning", upper_frame_title = "PDBQT generation Error",
                submit_command = None, with_scroll = True)

                self.pdbqt_log.middle_layout_type.addWidget(self.about_text_area, 0, 0)

                self.pdbqt_log.show()


    def remove_rec_from_list(self):

        Check_current_tab.check_docking_program_current_tab(self)

        if self.is_vina_tab:
            HandleWidgets.remove_dialog(self = self,
        prepared_objects_list = self.docking_programs_child_tabs.docking_programs.vina_receptors_prepared,
        widgets_list = [self.docking_programs_child_tabs.docking_programs.VINA.docking_tab_ui.ready_receptors_list],
        format = ".pdbqt",
        tmp_dir = self.docking_programs_child_tabs.docking_programs.vina_tmp_dir)

        if self.is_smina_tab:
            HandleWidgets.remove_dialog(self = self,
        prepared_objects_list = self.docking_programs_child_tabs.docking_programs.smina_receptors_prepared,
        widgets_list = [self.docking_programs_child_tabs.docking_programs.SMINA.docking_tab_ui.ready_receptors_list],
        format = ".pdb",
        tmp_dir = self.docking_programs_child_tabs.docking_programs.smina_tmp_dir)

        if self.is_adfr_tab:
            HandleWidgets.remove_dialog(self = self,
        prepared_objects_list = self.docking_programs_child_tabs.docking_programs.adfr_receptors_prepared,
        widgets_list = [self.docking_programs_child_tabs.docking_programs.ADFR.docking_tab_ui.ready_receptors_list],
        format = ".pdb",
        tmp_dir = self.docking_programs_child_tabs.docking_programs.adfr_tmp_dir)

        if self.is_rxdock_tab:
            HandleWidgets.remove_dialog(self = self,
        prepared_objects_list = self.docking_programs_child_tabs.docking_programs.rxdock_receptors_prepared,
        widgets_list = [self.docking_programs_child_tabs.docking_programs.RXDOCK.docking_tab_ui.ready_receptors_list,
        self.docking_programs_child_tabs.docking_programs.RXDOCK.grid_tab_rxdock_ui.ready_receptors_list,
        self.docking_programs_child_tabs.docking_programs.RXDOCK.grid_tab_rxdock_ui.ready_ligands_list],
        format = ".mol2",
        tmp_dir = self.docking_programs_child_tabs.docking_programs.rxdock_tmp_dir)


    def build_structures_frame(self, strc, dict):

        Check_current_tab.check_docking_program_current_tab(self)

        # Create a Frame for each loaded Receptor
        self.structure_frame = StructuresFrame(parent=None, main_window=self, strc_path=strc, dict = dict)

        # Update info in Receptors dictionary
        dict[strc]["frame"] = self.structure_frame

        # if self.is_vina_tab or self.is_smina_tab or self.is_adfr_tab:
        #
        #     # Add some options only to Vina Receptors' frames
        #     self.structure_frame.add_h = QtWidgets.QCheckBox("Add Hydrogens")
        #     self.structure_frame.structure_frame_layout.addWidget(self.structure_frame.add_h, 1, 1)
        #     self.structure_frame.add_h.setEnabled(False)
        #
        #     self.structure_frame.remove_nonstd = QtWidgets.QCheckBox("Remove ALL non-standard residues")
        #     self.structure_frame.structure_frame_layout.addWidget(self.structure_frame.remove_nonstd, 2, 1)
        #     self.structure_frame.remove_nonstd.setEnabled(False)
        #
        #     self.structure_frame.remove_water = QtWidgets.QCheckBox("Remove Water")
        #     self.structure_frame.structure_frame_layout.addWidget(self.structure_frame.remove_water, 3, 1)
        #     self.structure_frame.remove_water.setEnabled(False)

        # Add some options only to Receptors' frames
        self.tot_rows = 5

        self.structure_frame.het_id_text = QtWidgets.QLabel(" IDs: ")
        self.structure_frame.het_id_text.setStyleSheet("QLabel"
            "{"
            "font-size: 8px;"
            "}")
        self.structure_frame.structure_frame_layout.addWidget(self.structure_frame.het_id_text, 4, 2)
        self.structure_frame.het_id_text.hide()

        self.structure_frame.show_het_btn = QtWidgets.QRadioButton("Show Heteroatoms")
        self.structure_frame.show_het_btn.setEnabled(False)
        self.structure_frame.show_het_btn.clicked.connect(self.structure_frame.show_het_func)
        self.structure_frame.structure_frame_layout.addWidget(self.structure_frame.show_het_btn, 4, 1)

        self.structure_frame.het_text_list= []

        for het in self.structure_frame.heteroresidues:
            self.structure_frame.het_text = QtWidgets.QLabel(het.split("_")[1])
            self.structure_frame.structure_frame_layout.addWidget(self.structure_frame.het_text, self.tot_rows, 2)
            self.structure_frame.het_text.hide()
            self.structure_frame.het_text_list.append(self.structure_frame.het_text)
            self.structure_frame.het_text.setStyleSheet("QLabel"
            "{"
            "font-size: 6px;"
            "}")
            self.tot_rows += 1


        # Add the frame to the layout
        self.prepare_structures_frame_layout.addRow(dict[strc]["frame"])



class LigandTab(QtWidgets.QWidget, PyMOLInteractions, HandleWidgets):

    def __init__(self, main_window, current_tab):
        super().__init__(main_window)
        self.docking_programs_child_tabs = main_window
        self.initUI(current_tab)


    def initUI(self, current_tab):

        # Ligand Tab Layout
        self.layout_ligand_tab = QtWidgets.QVBoxLayout()

        # Create Group Boxes
        self.load_ligands_group = QtWidgets.QGroupBox("Ligands")
        self.prepare_ligands_group = QtWidgets.QGroupBox("Ligands settings")

        # Add the Group Boxes to the Receptor Tab Layout
        self.layout_ligand_tab.addWidget(self.load_ligands_group)
        self.layout_ligand_tab.addWidget(self.prepare_ligands_group)

        # Create and set Group Boxes' layouts
        self.load_ligands_layout = QtWidgets.QHBoxLayout()
        self.prepare_ligands_layout = QtWidgets.QHBoxLayout()

        self.load_ligands_group.setLayout(self.load_ligands_layout)
        self.prepare_ligands_group.setLayout(self.prepare_ligands_layout)

        # Create widgets for "Load Ligands group Box"
        self.import_btn = QtWidgets.QPushButton("Import from PyMOL")

        self.load_ligands_layout.addWidget(self.import_btn)
        self.load_ligands_layout.addStretch()

        # Create widgets for "Prepare Ligands group Box"
        self.main_page_scroll = QtWidgets.QScrollArea()
        self.main_page_interior = QtWidgets.QWidget()
        self.main_page_scroll.setWidgetResizable(True)

        self.main_page_scroll.setWidget(self.main_page_interior)
        self.prepare_ligands_layout.addWidget(self.main_page_scroll)
        # self.combine_button = QtWidgets.QPushButton("Combine Ligands")
        # self.combine_button.clicked.connect(self.combine_ligands_func)

        self.prepare_ligands_frame_layout = QtWidgets.QFormLayout()
        self.main_page_interior.setLayout(self.prepare_ligands_frame_layout)
        # self.prepare_ligands_frame_layout.addWidget(self.combine_button)

        # Add Select All Option
        self.select_all_btn = QtWidgets.QCheckBox("All")
        self.prepare_ligands_frame_layout.addWidget(self.select_all_btn)
        self.select_all_btn.clicked.connect(self.get_list_of_current_cb)

        # Add pdbqt options
        if current_tab != "RxDock":
            self.pdbqt_options_dict_lig = {}
            self.pdbqt_options_dict_lig["add_h"] = True
            self.pdbqt_options_dict_lig["none_torsions"] = False
            self.pdbqt_options_dict_lig["all_torsions"] = False
            self.pdbqt_options_dict_lig["all_but_ga"] = True

            self.pdbqt_options_pushbutton = QtWidgets.QPushButton("PDBQT Advanced Options")
            self.prepare_ligands_frame_layout.addRow(self.pdbqt_options_pushbutton)
            self.pdbqt_options_pushbutton.clicked.connect(self.show_pdbqt_options_window)

        self.middle_btn_layout = QtWidgets.QVBoxLayout()
        self.prepare_ligands_layout.addLayout(self.middle_btn_layout)
        self.gen_lig_btn = QtWidgets.QPushButton("Generate Ligand")
        self.gen_lig_btn.clicked.connect(self.generate_ligand_func)

        self.middle_btn_layout.addStretch()
        self.middle_btn_layout.addWidget(self.gen_lig_btn)
        self.middle_btn_layout.addStretch()

        self.middle_right_layout = QtWidgets.QVBoxLayout()
        self.prepare_ligands_layout.addLayout(self.middle_right_layout)

        self.listwidget = QtWidgets.QListWidget()
        self.listwidget.setFixedWidth(70)
        self.listwidget.setSelectionMode(3)
        self.middle_right_layout.addWidget(self.listwidget)

        self.set = QtWidgets.QPushButton("Set")
        self.remove = QtWidgets.QPushButton("Remove")
        self.middle_right_layout.addWidget(self.set)
        self.middle_right_layout.addWidget(self.remove)

        #self.remove.hide()

        # Functions
        self.set.clicked.connect(self.add_to_other_tabs)
        self.remove.clicked.connect(self.remove_lig_from_list)
        self.import_btn.clicked.connect(self.import_lig_from_pymol_current_tab)


    def show_pdbqt_options_window(self):

        pdbqt_options_window = PDBQT_OptionsWindows(tab = self, main = self.docking_programs_child_tabs,
                             obj_type = "ligand", options_dict = self.pdbqt_options_dict_lig)

        self.pdbqt_options_dict_lig = pdbqt_options_window.options_dict

        # self.pdbqt_options_window = NewWindow(parent = self.docking_programs_child_tabs,
        # title = "PDBQT options window", upper_frame_title = "Select Options",
        # submit_command = self.apply_pdbqt_options, submit_button_text= "Set",
        # with_scroll = True)
        #
        # self.add_h = QtWidgets.QCheckBox("Add Hydrogens")
        # self.pdbqt_options_window.middle_layout_type.addWidget(self.add_h)
        #
        # self.active_torsions_group = QtWidgets.QGroupBox("Active Torsions")
        # self.pdbqt_options_window.middle_layout_type.addWidget(self.active_torsions_group)
        # self.active_torsions_group_layout = QtWidgets.QVBoxLayout()
        # self.active_torsions_group.setLayout(self.active_torsions_group_layout)
        #
        # self.none_torsions = QtWidgets.QRadioButton("None")
        # self.all_torsions = QtWidgets.QRadioButton("All")
        # self.all_but_ga = QtWidgets.QRadioButton("All But Guanidinium and Amide")
        #
        # self.active_torsions_group_layout.addWidget(self.none_torsions)
        # self.active_torsions_group_layout.addWidget(self.all_torsions)
        # self.active_torsions_group_layout.addWidget(self.all_but_ga)
        #
        # self.add_h.setChecked(self.pdbqt_options_dict_lig["add_h"])
        # self.none_torsions.setChecked(self.pdbqt_options_dict_lig["none_torsions"])
        # self.all_torsions.setChecked(self.pdbqt_options_dict_lig["all_torsions"])
        # self.all_but_ga.setChecked(self.pdbqt_options_dict_lig["all_but_ga"])
        #
        # self.pdbqt_options_window.show()
        #

    def apply_pdbqt_options(self):

        self.pdbqt_options_dict_lig["add_h"] = self.add_h.isChecked()
        self.pdbqt_options_dict_lig["none_torsions"] = self.none_torsions.isChecked()
        self.pdbqt_options_dict_lig["all_torsions"] = self.all_torsions.isChecked()
        self.pdbqt_options_dict_lig["all_torsions"] = self.all_but_ga.isChecked()

        self.pdbqt_options_window.close()


    def get_list_of_current_cb(self):

        self.list_of_cb = []

        for i in reversed(range(self.prepare_ligands_frame_layout.count())):
            obj = self.prepare_ligands_frame_layout.itemAt(i).widget()
            if str(type(obj)) == "<class 'lib.docking_program_main.docking_program_gui.frames.LigandsFrame'>":
                self.list_of_cb.append(obj.strc_checkbox)

        SelectAll(self, all_cb = self.select_all_btn,
        list_of_cb = self.list_of_cb)


    def combine_ligands_func(self):

        Check_current_tab.check_docking_program_current_tab(self)

        if self.is_vina_tab:
            dict = self.docking_programs_child_tabs.docking_programs.vina_ligands_dict

        elif self.is_rxdock_tab:
            dict = self.docking_programs_child_tabs.docking_programs.rxdock_ligands_dict

        self.ligands_to_combine = []

        # Get checked ligands
        for strc in dict:
            # For each object that is checked
            if dict[strc]["frame"].strc_checkbox.isChecked():
                self.ligands_to_combine.append(dict[strc]["frame"].strc_checkbox.text())

        # Define the combined ligands object name
        self.combined_ligands_name = str("Object_combined" + str(len(self.docking_programs_child_tabs.docking_programs.combined_objs_list)))

        # If the list of combined ligands is not empty, combine the objects into a single one
        if self.ligands_to_combine:
            qm = QtWidgets.QMessageBox.question(self.docking_programs_child_tabs,'Combine Ligands', str("Do you want to combine these ligands into a single object?\n\n" +
            ', '.join(self.ligands_to_combine) + "will become " + self.combined_ligands_name), QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

            if qm == QtWidgets.QMessageBox.Yes:

                for ligands in self.ligands_to_combine:
                    pymol.creating.join_states(self.combined_ligands_name, ligands.split()[0])

                self.docking_programs_child_tabs.docking_programs.combined_objs_list.append(self.combined_ligands_name)

            elif qm == QtWidgets.QMessageBox.No:
                self.is_already_loaded = True

        else:
            QtWidgets.QMessageBox.warning(self.docking_programs_child_tabs, "Warning", "Check the ligands you wish to combine")
            combine = False


    def remove_lig_from_list(self):

        Check_current_tab.check_docking_program_current_tab(self)

        if self.is_vina_tab:
            os.chdir(self.docking_programs_child_tabs.docking_programs.vina_tmp_dir)
            HandleWidgets.remove_dialog(self = self,
        prepared_objects_list = self.docking_programs_child_tabs.docking_programs.vina_ligands_prepared,
        widgets_list = [self.docking_programs_child_tabs.docking_programs.VINA.docking_tab_ui.ready_ligands_list],
        format = ".pdbqt",
        tmp_dir = self.docking_programs_child_tabs.docking_programs.vina_tmp_dir)

        if self.is_rxdock_tab:
            os.chdir(self.docking_programs_child_tabs.docking_programs.rxdock_tmp_dir)
            HandleWidgets.remove_dialog(self = self,
        prepared_objects_list = self.docking_programs_child_tabs.docking_programs.rxdock_ligands_prepared,
        widgets_list = [self.docking_programs_child_tabs.docking_programs.RXDOCK.docking_tab_ui.ready_ligands_list,
        self.docking_programs_child_tabs.docking_programs.RXDOCK.grid_tab_rxdock_ui.ready_receptors_list,
        self.docking_programs_child_tabs.docking_programs.RXDOCK.grid_tab_rxdock_ui.ready_ligands_list],
        format = ".sdf",
        tmp_dir = self.docking_programs_child_tabs.docking_programs.rxdock_tmp_dir)

        if self.is_smina_tab:
            os.chdir(self.docking_programs_child_tabs.docking_programs.smina_tmp_dir)
            HandleWidgets.remove_dialog(self = self,
        prepared_objects_list = self.docking_programs_child_tabs.docking_programs.smina_ligands_prepared,
        widgets_list = [self.docking_programs_child_tabs.docking_programs.SMINA.docking_tab_ui.ready_receptors_list],
        format = ".pdbqt",
        tmp_dir = self.docking_programs_child_tabs.docking_programs.smina_tmp_dir)

        if self.is_adfr_tab:
            os.chdir(self.docking_programs_child_tabs.docking_programs.adfr_tmp_dir)
            HandleWidgets.remove_dialog(self = self,
        prepared_objects_list = self.docking_programs_child_tabs.docking_programs.adfr_ligands_prepared,
        widgets_list = [self.docking_programs_child_tabs.docking_programs.ADFR.docking_tab_ui.ready_receptors_list],
        format = ".pdbqt",
        tmp_dir = self.docking_programs_child_tabs.docking_programs.adfr_tmp_dir)


    def import_lig_from_pymol_current_tab(self):

        a = self.parent().docking_programs.docking_programs_tabs.currentIndex()
        current_tab = self.parent().docking_programs.docking_programs_tabs.tabText(a)

        imported_ligand = Import_from_Pymol(self, current_tab, is_receptor = False)


    def build_structures_frame(self, strc, dict):

        Check_current_tab.check_docking_program_current_tab(self)

        # Build a Frame for each imported Ligand
        self.structure_frame = LigandsFrame(parent=None, main_window=self, strc_path=strc, dict = dict)

        # Update info in ligands dictionary
        dict[strc]["frame"] = self.structure_frame

        # if self.is_vina_tab or self.is_smina_tab or self.is_adfr_tab:
        #
        #     # Add some options only to Vina Ligands' frames
        #     self.structure_frame.add_h = QtWidgets.QCheckBox("Add Hydrogens")
        #     self.structure_frame.structure_frame_layout.addWidget(self.structure_frame.add_h, 2, 1)
        #     self.structure_frame.add_h.setEnabled(False)
        #
        #     self.structure_frame.active_torsions_group = QtWidgets.QGroupBox("Active Torsions")
        #     self.structure_frame.structure_frame_layout.addWidget(self.structure_frame.active_torsions_group, 1, 1)
        #     self.structure_frame.active_torsions_layout = QtWidgets.QVBoxLayout()
        #     self.structure_frame.active_torsions_group.setLayout(self.structure_frame.active_torsions_layout)
        #
        #     self.structure_frame.active_torsions_group.setEnabled(False)
        #
        #     self.structure_frame.all_but_torsions = QtWidgets.QRadioButton("All but Guanidinium and Amide")
        #     self.structure_frame.all_torsions = QtWidgets.QRadioButton("All")
        #     self.structure_frame.none_torsions = QtWidgets.QRadioButton("None")
        #     self.structure_frame.active_torsions_layout.addWidget(self.structure_frame.all_but_torsions)
        #     self.structure_frame.active_torsions_layout.addWidget(self.structure_frame.all_torsions)
        #     self.structure_frame.active_torsions_layout.addWidget(self.structure_frame.none_torsions)
        #     self.structure_frame.all_but_torsions.setEnabled(False)
        #     self.structure_frame.all_torsions.setEnabled(False)
        #     self.structure_frame.none_torsions.setEnabled(False)

        # Add to the layout
        self.prepare_ligands_frame_layout.addRow(dict[strc]["frame"])


    def generate_ligand_func(self):

        if not check_configuration(self, self.docking_programs_child_tabs.docking_programs):
            QtWidgets.QMessageBox.warning(self.docking_programs_child_tabs, "Warning", "Please check DockingPie configuration!")
        else:

            Check_current_tab.check_docking_program_current_tab(self)

            if self.is_vina_tab:
                self.generated_ligand = Generate_Object(self, main = self.docking_programs_child_tabs.docking_programs,
                dict = self.docking_programs_child_tabs.docking_programs.vina_ligands_dict,
                prepared_objects_list = self.docking_programs_child_tabs.docking_programs.vina_ligands_prepared,
                tmp_path = self.docking_programs_child_tabs.docking_programs.vina_tmp_dir,
                format = "sdf",
                docking_program_name = "Vina",
                generate_pdbqt = True,
                pdbqt_dict = self.pdbqt_options_dict_lig,
                is_receptor = False)

                self.generated_object_update_gui()
                cmd.group("Vina", members=self.generated_ligand.new_strc_name, action='auto', quiet=1)

            if self.is_rxdock_tab:
                self.generated_ligand = Generate_Object(self, main = self.docking_programs_child_tabs.docking_programs,
                dict = self.docking_programs_child_tabs.docking_programs.rxdock_ligands_dict,
                prepared_objects_list = self.docking_programs_child_tabs.docking_programs.rxdock_ligands_prepared,
                tmp_path = self.docking_programs_child_tabs.docking_programs.rxdock_tmp_dir,
                format = "sdf",
                docking_program_name = "RxDock",
                generate_pdbqt = False,
                is_receptor = False)

                self.generated_object_update_gui()
                cmd.group("RxDock", members=self.generated_ligand.new_strc_name, action='auto', quiet=1)

            if self.is_smina_tab:
                self.generated_ligand = Generate_Object(self, main = self.docking_programs_child_tabs.docking_programs,
                dict = self.docking_programs_child_tabs.docking_programs.smina_ligands_dict,
                prepared_objects_list = self.docking_programs_child_tabs.docking_programs.smina_ligands_prepared,
                tmp_path = self.docking_programs_child_tabs.docking_programs.smina_tmp_dir,
                format = "sdf",
                docking_program_name = "Smina",
                generate_pdbqt = True,
                pdbqt_dict = self.pdbqt_options_dict_lig,
                is_receptor = False)

                self.generated_object_update_gui()
                cmd.group("Smina", members=self.generated_ligand.new_strc_name, action='auto', quiet=1)

            if self.is_adfr_tab:
                self.generated_ligand = Generate_Object(self, main = self.docking_programs_child_tabs.docking_programs,
                dict = self.docking_programs_child_tabs.docking_programs.adfr_ligands_dict,
                prepared_objects_list = self.docking_programs_child_tabs.docking_programs.adfr_ligands_prepared,
                tmp_path = self.docking_programs_child_tabs.docking_programs.adfr_tmp_dir,
                format = "sdf",
                docking_program_name = "ADFR",
                generate_pdbqt = True,
                pdbqt_dict = self.pdbqt_options_dict_lig,
                is_receptor = False)

                self.generated_object_update_gui()
                cmd.group("ADFR", members=self.generated_ligand.new_strc_name, action='auto', quiet=1)


    def generated_object_update_gui(self):

        if os.path.isfile(self.generated_ligand.new_strc_name_format):

            # Add to the listwidget
            self.listbtn = self.listwidget.addItem(self.generated_ligand.new_strc_name)

            # Add to the list of prepared files
            self.generated_ligand.prepared_objects_list.append(self.generated_ligand.new_strc_name)

            # Load in PyMOL
            cmd.load(self.generated_ligand.new_receptor_file_path, self.generated_ligand.new_strc_name)

            self.docking_programs_child_tabs.docking_programs.main_window.statusBar().showMessage(str("Generated Receptor '" + self.generated_ligand.new_strc_name + "'"), 3000)

        else:
            if self.generated_ligand.generate_pdbqt:

                self.about_text_area = QtWidgets.QPlainTextEdit()
                self.about_text_area.setReadOnly(True)

                self.about_text_area.setPlainText(open(self.pdbqt_log_file_name).read())

                self.pdbqt_log = NewWindow(parent = self.tab,
                title = "Warning", upper_frame_title = "PDBQT generation Error",
                submit_command = None, with_scroll = True)

                self.pdbqt_log.middle_layout_type.addWidget(self.about_text_area, 0, 0)

                self.pdbqt_log.show()

    def add_to_other_tabs(self):

        Check_current_tab.check_docking_program_current_tab(self)

        if self.is_vina_tab:
            HandleWidgets.add_to_other_tab(self = self, input_widget = self.listwidget, output_widget = self.docking_programs_child_tabs.docking_programs.VINA.docking_tab_ui.ready_ligands_list)

        if self.is_rxdock_tab:
            HandleWidgets.add_to_other_tab(self = self, input_widget = self.listwidget, output_widget = self.docking_programs_child_tabs.docking_programs.RXDOCK.grid_tab_rxdock_ui.ready_ligands_list)
            HandleWidgets.add_to_other_tab(self = self, input_widget = self.listwidget, output_widget = self.docking_programs_child_tabs.docking_programs.RXDOCK.docking_tab_ui.ready_ligands_list)

        if self.is_smina_tab:
            HandleWidgets.add_to_other_tab(self = self, input_widget = self.listwidget, output_widget = self.docking_programs_child_tabs.docking_programs.SMINA.docking_tab_ui.ready_ligands_list)

        if self.is_adfr_tab:
            HandleWidgets.add_to_other_tab(self = self, input_widget = self.listwidget, output_widget = self.docking_programs_child_tabs.docking_programs.ADFR.docking_tab_ui.ready_ligands_list)

        if self.listwidget.count() == 0:
            pass
        else:
           self.docking_programs_child_tabs.docking_programs.main_window.statusBar().showMessage("Selected Items added to other tabs", 3000)


class DockingTab(QtWidgets.QWidget, PyMOLInteractions, HandleWidgets):


    """
    A class to reppresent the Docking Tab Layout.
    It is the same for each Docking Program, except for the OptionsFrame which is specific.
    """


    def __init__(self, main_window, current_tab):
        super().__init__(main_window)
        self.docking_programs_child_tabs = main_window
        self.initUI(current_tab)


    def initUI(self, current_tab):

        # Docking Tab Layout
        self.layout_docking_tab = QtWidgets.QVBoxLayout()

        # Create Group box
        self.group_docking = QtWidgets.QGroupBox("Docking")

        # Add the Group Box to the Docking Tab Layout and set its layout
        self.group_docking_layout = QtWidgets.QGridLayout()
        self.group_docking.setLayout(self.group_docking_layout)
        self.layout_docking_tab.addWidget(self.group_docking)

        ## WIDGETS FOR DOCKING SETTINGS
        self.options_frame_canonical = OptionsFrameCanonical(parent=None,
    main_window=self.docking_programs_child_tabs)

        # Add to the Docking Tab
        self.group_docking_layout.addWidget(self.options_frame_canonical, 1, 0, 7, 1)

        # Add Combobox to choose the protocol
        self.protocol_box = QtWidgets.QComboBox()
        self.protocol_box.addItems(["ALLvsALL", "Custom"])
        self.group_docking_layout.addWidget(self.protocol_box, 0, 0)
        self.protocol_box.currentTextChanged.connect(self.get_general_docking_protocol)

        ### GENERAL DOCKING OPTIONS

        self.options_frame_canonical = OptionsFrameCanonical(parent=None,
    main_window=self.docking_programs_child_tabs)
        self.group_docking_layout.addWidget(self.options_frame_canonical, 1, 0, 7, 1)

        self.loaded_cavities = QtWidgets.QComboBox()
        self.loaded_cavities_text = QtWidgets.QLabel("Available cavities")
        self.loaded_cavities.currentTextChanged.connect(lambda: HandleWidgets.combobox_orient_current_item_pymol(self = self, widget = self.loaded_cavities))
        self.loaded_cavities.view().pressed.connect(lambda: HandleWidgets.combobox_orient_current_item_pymol(self = self, widget = self.loaded_cavities))

        self.options_frame_canonical.options_frame_layout.addWidget(self.loaded_cavities_text, 0, 0, 1, 2)
        self.options_frame_canonical.options_frame_layout.addWidget(self.loaded_cavities, 1, 0, 1, 2)

        self.ready_receptors_list_text = QtWidgets.QLabel("Receptor(s)")

        # Create the Scroll Area update_current_ready_objectsfor the Table
        widget = QtWidgets.QWidget()
        self.ready_receptors_list = QtWidgets.QScrollArea()
        self.ready_receptors_list.setWidgetResizable(True)
        self.ready_receptors_list.setWidget(widget)

        # Set the layout of the Scroll Area for the Table
        self.receptors_scroll_layout = QtWidgets.QGridLayout()
        widget.setLayout(self.receptors_scroll_layout)

        self.options_frame_canonical.options_frame_layout.addWidget(self.ready_receptors_list_text, 2, 0)
        self.options_frame_canonical.options_frame_layout.addWidget(self.ready_receptors_list, 3, 0)

        self.ready_ligands_list_text = QtWidgets.QLabel("Ligand(s)")

        # Create the Scroll Area for the Table
        widget = QtWidgets.QWidget()
        self.ready_ligands_list = QtWidgets.QScrollArea()
        self.ready_ligands_list.setWidgetResizable(True)
        self.ready_ligands_list.setWidget(widget)

        # Set the layout of the Scroll Area for the Table
        self.ligand_scroll_layout = QtWidgets.QGridLayout()
        widget.setLayout(self.ligand_scroll_layout)

        self.options_frame_canonical.options_frame_layout.addWidget(self.ready_ligands_list_text, 2, 1)
        self.options_frame_canonical.options_frame_layout.addWidget(self.ready_ligands_list, 3, 1)

        self.run_btn = QtWidgets.QPushButton("Run Docking")
        self.run_btn.clicked.connect(self.run_docking_func)
        self.group_docking_layout.addWidget(self.run_btn, 9, 0, 1, 2)

        self.create_options_frame(current_tab)


    def get_general_docking_protocol(self):

        protocol = self.protocol_box.currentText()

        self.get_current_ready_objects()

        HandleWidgets.remove_general_frame(self = self,
        layout = self.options_frame_canonical.options_frame_layout,
        frame = self.options_frame_canonical)

        if protocol == "Normal":

            self.options_frame_canonical = OptionsFrameCanonical(parent=None,
        main_window=self.docking_programs_child_tabs)

            # Add to the Docking Tab
            self.group_docking_layout.addWidget(self.options_frame_canonical, 1, 0, 7, 1)

            self.loaded_cavities = QtWidgets.QComboBox()
            self.loaded_cavities_text = QtWidgets.QLabel("Available cavities")
            self.loaded_cavities.currentTextChanged.connect(lambda: HandleWidgets.combobox_orient_current_item_pymol(self = self, widget = self.loaded_cavities))
            self.loaded_cavities.view().pressed.connect(lambda: HandleWidgets.combobox_orient_current_item_pymol(self = self, widget = self.loaded_cavities))

            self.options_frame_canonical.options_frame_layout.addWidget(self.loaded_cavities_text, 0, 0, 1, 2)
            self.options_frame_canonical.options_frame_layout.addWidget(self.loaded_cavities, 1, 0, 1, 2)

            self.ready_receptors_list = QtWidgets.QComboBox()
            self.ready_receptors_list_text = QtWidgets.QLabel("Receptor(s)")
            self.ready_receptors_list.currentTextChanged.connect(lambda: HandleWidgets.combobox_orient_current_item_pymol(self = self, widget = self.ready_receptors_list))
            self.ready_receptors_list.view().pressed.connect(lambda: HandleWidgets.combobox_orient_current_item_pymol(self = self, widget = self.ready_receptors_list))

            self.options_frame_canonical.options_frame_layout.addWidget(self.ready_receptors_list_text, 2, 0, 1, 2)
            self.options_frame_canonical.options_frame_layout.addWidget(self.ready_receptors_list, 3, 0, 1, 2)

            self.ready_ligands_list = QtWidgets.QComboBox()
            self.ready_ligands_list_text = QtWidgets.QLabel("Ligand(s)")
            self.ready_ligands_list.currentTextChanged.connect(lambda: HandleWidgets.combobox_orient_current_item_pymol(self = self, widget = self.ready_receptors_list))
            self.ready_ligands_list.view().pressed.connect(lambda: HandleWidgets.combobox_orient_current_item_pymol(self = self, widget = self.ready_receptors_list))

            self.options_frame_canonical.options_frame_layout.addWidget(self.ready_ligands_list_text, 4, 0, 1, 2)
            self.options_frame_canonical.options_frame_layout.addWidget(self.ready_ligands_list, 5, 0, 1, 2)

        if protocol == "Custom":

            self.dockings_to_do_dict = {}

            self.options_frame_canonical = OptionsFrameCanonical(parent=None,
        main_window=self.docking_programs_child_tabs)

            # Add to the Docking Tab
            self.group_docking_layout.addWidget(self.options_frame_canonical, 1, 0, 7, 1)

            self.create_new_group_button = QtWidgets.QPushButton("Create New Group")
            self.options_frame_canonical.options_frame_layout.addWidget(self.create_new_group_button, 0, 0, 1, 2)
            self.create_new_group_button.clicked.connect(self.open_new_group_window)

            # Create the Scroll Area for the Table
            widget = QtWidgets.QWidget()
            self.group_scroll = QtWidgets.QScrollArea()
            self.group_scroll.setWidgetResizable(True)
            self.group_scroll.setWidget(widget)
            # Set the layout of the Scroll Area for the Table
            self.group_scroll_layout = QtWidgets.QVBoxLayout()
            widget.setLayout(self.group_scroll_layout)

            self.loaded_cavities = QtWidgets.QComboBox()
            self.loaded_cavities_text = QtWidgets.QLabel("Available cavities")
            self.loaded_cavities.currentTextChanged.connect(lambda: HandleWidgets.combobox_orient_current_item_pymol(self = self, widget = self.loaded_cavities))
            self.loaded_cavities.view().pressed.connect(lambda: HandleWidgets.combobox_orient_current_item_pymol(self = self, widget = self.loaded_cavities))

            self.ready_receptors_list = QtWidgets.QComboBox()
            self.ready_receptors_list_text = QtWidgets.QLabel("Receptor")
            self.ready_receptors_list.currentTextChanged.connect(lambda: HandleWidgets.combobox_orient_current_item_pymol(self = self, widget = self.ready_receptors_list))
            self.ready_receptors_list.view().pressed.connect(lambda: HandleWidgets.combobox_orient_current_item_pymol(self = self, widget = self.ready_receptors_list))

            self.ready_ligands_list = QtWidgets.QComboBox()
            self.ready_ligands_list_text = QtWidgets.QLabel("Ligand(s)")
            self.ready_ligands_list.currentTextChanged.connect(lambda: HandleWidgets.combobox_orient_current_item_pymol(self = self, widget = self.ready_receptors_list))
            self.ready_ligands_list.view().pressed.connect(lambda: HandleWidgets.combobox_orient_current_item_pymol(self = self, widget = self.ready_receptors_list))

            self.options_frame_canonical.options_frame_layout.addWidget(self.group_scroll, 4, 0)

        if protocol == "ALLvsALL":

            self.options_frame_canonical = OptionsFrameCanonical(parent=None,
        main_window=self.docking_programs_child_tabs)
            self.group_docking_layout.addWidget(self.options_frame_canonical, 1, 0, 7, 1)

            self.loaded_cavities = QtWidgets.QComboBox()
            self.loaded_cavities_text = QtWidgets.QLabel("Available cavities")
            self.loaded_cavities.currentTextChanged.connect(lambda: HandleWidgets.combobox_orient_current_item_pymol(self = self, widget = self.loaded_cavities))
            self.loaded_cavities.view().pressed.connect(lambda: HandleWidgets.combobox_orient_current_item_pymol(self = self, widget = self.loaded_cavities))

            self.options_frame_canonical.options_frame_layout.addWidget(self.loaded_cavities_text, 0, 0, 1, 2)
            self.options_frame_canonical.options_frame_layout.addWidget(self.loaded_cavities, 1, 0, 1, 2)

            self.ready_receptors_list_text = QtWidgets.QLabel("Receptor(s)")
            # Create the Scroll Area update_current_ready_objectsfor the Table
            widget = QtWidgets.QWidget()
            self.ready_receptors_list = QtWidgets.QScrollArea()
            self.ready_receptors_list.setWidgetResizable(True)
            self.ready_receptors_list.setWidget(widget)
            # Set the layout of the Scroll Area for the Table
            self.receptors_scroll_layout = QtWidgets.QGridLayout()
            widget.setLayout(self.receptors_scroll_layout)

            self.options_frame_canonical.options_frame_layout.addWidget(self.ready_receptors_list_text, 2, 0)
            self.options_frame_canonical.options_frame_layout.addWidget(self.ready_receptors_list, 3, 0)

            self.ready_ligands_list_text = QtWidgets.QLabel("Ligand(s)")

            # Create the Scroll Area for the Table
            widget = QtWidgets.QWidget()
            self.ready_ligands_list = QtWidgets.QScrollArea()
            self.ready_ligands_list.setWidgetResizable(True)
            self.ready_ligands_list.setWidget(widget)
            # Set the layout of the Scroll Area for the Table
            self.ligand_scroll_layout = QtWidgets.QGridLayout()
            widget.setLayout(self.ligand_scroll_layout)

            self.options_frame_canonical.options_frame_layout.addWidget(self.ready_ligands_list_text, 2, 1)
            self.options_frame_canonical.options_frame_layout.addWidget(self.ready_ligands_list, 3, 1)

        if protocol == "Redocking":
            print("work in progess")

        self.update_current_ready_objects()


    def open_new_group_window(self):

        self.new_group_window = NewWindow(parent = self.docking_programs_child_tabs,
        title = "New Group", upper_frame_title = "Select Inputs",
        submit_command = self.create_new_group, submit_button_text= "Set",
        with_scroll = True)

        self.new_group_window.middle_layout_type.addWidget(self.loaded_cavities_text, 0, 0, 1, 2)
        self.new_group_window.middle_layout_type.addWidget(self.loaded_cavities, 1, 0, 1, 2)

        self.new_group_window.middle_layout_type.addWidget(self.ready_receptors_list_text, 2, 0, 1, 2)
        self.new_group_window.middle_layout_type.addWidget(self.ready_receptors_list, 3, 0, 1, 2)

        self.new_group_window.middle_layout_type.addWidget(self.ready_ligands_list_text, 4, 0, 1, 2)
        self.new_group_window.middle_layout_type.addWidget(self.ready_ligands_list, 5, 0, 1, 2)

        self.new_group_window.show()


    def create_new_group(self):

        HandleWidgets.combobox_check_if_empty(self = self,
        widgets_list = [self.loaded_cavities, self.ready_receptors_list, self.ready_ligands_list])

        list_of_groups = []

        for i in range(self.group_scroll_layout.count()):
            list_of_groups.append(self.group_scroll_layout.itemAt(i))

        if self.is_empty:
            QtWidgets.QMessageBox.warning(self.docking_programs_child_tabs.docking_programs.main_window, "Warning", "There aren't enough parameters")

        else:

            self.docking_group = QtWidgets.QGroupBox("Group " + str(len(list_of_groups)+1))
            self.group_scroll_layout.addWidget(self.docking_group)
            self.docking_group_layout = QtWidgets.QVBoxLayout()
            self.docking_group.setLayout(self.docking_group_layout)
            self.docking_group.setCheckable(True)

            cavity = self.loaded_cavities.currentText()
            ligand = self.ready_ligands_list.currentText()
            receptor = self.ready_receptors_list.currentText()

            cavity_label = QtWidgets.QLabel("Cavity: " + cavity)
            ligand_label = QtWidgets.QLabel("Ligand: " + ligand)
            receptor_label = QtWidgets.QLabel("Receptor: " + receptor)

            self.docking_group_layout.addWidget(cavity_label)
            self.docking_group_layout.addWidget(ligand_label)
            self.docking_group_layout.addWidget(receptor_label)

            self.dockings_to_do_dict[self.docking_group.title()] = {}
            self.dockings_to_do_dict[self.docking_group.title()]["receptor"] = receptor
            self.dockings_to_do_dict[self.docking_group.title()]["ligand"] = ligand
            self.dockings_to_do_dict[self.docking_group.title()]["cavity"] = cavity


    def get_current_ready_objects(self):

        self.tmp_list_of_cavities = []

        for cav in range(self.loaded_cavities.count()):
            self.tmp_list_of_cavities.extend([self.loaded_cavities.itemText(cav)])

        if str(type(self.ready_ligands_list)) == "<class 'PyQt5.QtWidgets.QComboBox'>":

            self.tmp_list_of_ligand = []

            for lig in range(self.ready_ligands_list.count()):

                self.tmp_list_of_ligand.extend([self.ready_ligands_list.itemText(lig)])

        if str(type(self.ready_ligands_list)) == "<class 'PyQt5.QtWidgets.QScrollArea'>":

            self.tmp_list_of_ligand = []

            for i in reversed(range(self.ligand_scroll_layout.count())):
                self.tmp_list_of_ligand.append(self.ligand_scroll_layout.itemAt(i).widget().text())

        if str(type(self.ready_receptors_list)) == "<class 'PyQt5.QtWidgets.QScrollArea'>":

            self.tmp_list_of_receptors = []

            for i in reversed(range(self.receptors_scroll_layout.count())):
                self.tmp_list_of_receptors.append(self.receptors_scroll_layout.itemAt(i).widget().text())

        else:

            self.tmp_list_of_receptors = []

            for rec in range(self.ready_receptors_list.count()):
                self.tmp_list_of_receptors.extend([self.ready_receptors_list.itemText(rec)])


    def update_current_ready_objects(self):

        self.loaded_cavities.addItems(self.tmp_list_of_cavities)

        if str(type(self.ready_ligands_list)) == "<class 'PyQt5.QtWidgets.QComboBox'>":

            self.ready_ligands_list.addItems(self.tmp_list_of_ligand)

        if str(type(self.ready_ligands_list)) == "<class 'PyQt5.QtWidgets.QScrollArea'>":

            for lig in self.tmp_list_of_ligand:
                cb = QtWidgets.QCheckBox(lig)
                self.ligand_scroll_layout.addWidget(cb)

        if str(type(self.ready_receptors_list)) == "<class 'PyQt5.QtWidgets.QComboBox'>":

            self.ready_receptors_list.addItems(self.tmp_list_of_receptors)

        if str(type(self.ready_receptors_list)) == "<class 'PyQt5.QtWidgets.QScrollArea'>":

            for rec in self.tmp_list_of_receptors:
                cb = QtWidgets.QCheckBox(rec)
                self.receptors_scroll_layout.addWidget(cb)


    def create_options_frame(self, current_tab):

        #---------------
        # It creates a specific OptionsFrame for each Docking Program
        #---------------


        if current_tab == "ADFR":

            self.options_frame_all = OptionsFrame(parent=None,
        main_window=self.docking_programs_child_tabs)

            # Add Both OptionsFrame to the Docking Tab
            self.group_docking_layout.addWidget(self.options_frame_all, 0, 1, 3, 1)

            self.buffer_label = QtWidgets.QLabel("Padding")
            self.buffer_box = QtWidgets.QSpinBox()
            self.buffer_box.setRange(1, 10)
            self.buffer_box.setSingleStep(1)
            self.buffer_box.setValue(4)

            self.options_frame_all.options_frame_layout.addWidget(self.buffer_label, 0, 0)
            self.options_frame_all.options_frame_layout.addWidget(self.buffer_box, 0, 1)

            # Num poses
            self.ga_evol_label = QtWidgets.QLabel("Number of GA evolutions")
            self.ga_evol = QtWidgets.QSpinBox()
            self.ga_evol.setRange(1, 60)
            self.ga_evol.setSingleStep(1)
            self.ga_evol.setValue(20)
            self.options_frame_all.options_frame_layout.addWidget(self.ga_evol_label, 1, 0)
            self.options_frame_all.options_frame_layout.addWidget(self.ga_evol, 1, 1)

            self.ga_threshold_label = QtWidgets.QLabel("GA no improvement stop threshold")
            self.ga_threshold = QtWidgets.QSpinBox()
            self.ga_threshold.setRange(1, 15)
            self.ga_threshold.setSingleStep(1)
            self.ga_threshold.setValue(5)

            self.options_frame_all.options_frame_layout.addWidget(self.ga_threshold_label, 2, 0)
            self.options_frame_all.options_frame_layout.addWidget(self.ga_threshold, 2, 1)

            # self.max_ga_eval_label = QtWidgets.QLabel("Max n° of evaluations per GA (mln)")
            # self.max_ga_eval = QtWidgets.QSpinBox()
            # self.max_ga_eval.setRange(1, 3)
            # self.max_ga_eval.setSingleStep(1)
            # self.max_ga_eval.setValue(2)

            # self.options_frame_all.options_frame_layout.addWidget(self.max_ga_eval_label, 2, 0)
            # self.options_frame_all.options_frame_layout.addWidget(self.max_ga_eval, 2, 1)

            self.max_gen_label = QtWidgets.QLabel("Max n° of generations (mln)")
            self.max_gen = QtWidgets.QSpinBox()
            self.max_gen.setRange(1, 15)
            self.max_gen.setSingleStep(1)
            self.max_gen.setValue(5)

            self.options_frame_all.options_frame_layout.addWidget(self.max_gen_label, 3, 0)
            self.options_frame_all.options_frame_layout.addWidget(self.max_gen, 3, 1)


        if current_tab == "Vina" or current_tab == "Smina":

            self.options_frame_all = OptionsFrame(parent=None,
        main_window=self.docking_programs_child_tabs)

            # Add Both OptionsFrame to the Docking Tab
            self.group_docking_layout.addWidget(self.options_frame_all, 0, 1, 3, 1)

            # Num poses
            self.poses_label = QtWidgets.QLabel("Poses")
            self.poses_box = QtWidgets.QSpinBox()
            self.poses_box.setRange(1, 100)
            self.poses_box.setSingleStep(1)
            self.poses_box.setValue(1)
            self.options_frame_all.options_frame_layout.addWidget(self.poses_label, 0, 0)
            self.options_frame_all.options_frame_layout.addWidget(self.poses_box, 0, 1)

            # self.group_docking_layout.addWidget(self.options_frame, 0, 1, 8, 1)

            self.exhaustiveness_label = QtWidgets.QLabel("Exhaustiveness")
            self.exhaustiveness_box = QtWidgets.QSpinBox()
            self.exhaustiveness_box.setRange(1, 40)
            self.exhaustiveness_box.setSingleStep(1)
            self.exhaustiveness_box.setValue(8)

            self.options_frame_all.options_frame_layout.addWidget(self.exhaustiveness_label, 1, 0)
            self.options_frame_all.options_frame_layout.addWidget(self.exhaustiveness_box, 1, 1)

            self.energy_label = QtWidgets.QLabel("Energy Range")
            self.energy_box = QtWidgets.QSpinBox()
            self.energy_box.setRange(1, 10)
            self.energy_box.setSingleStep(1)
            self.energy_box.setValue(3)
            self.options_frame_all.options_frame_layout.addWidget(self.energy_label, 2, 0)
            self.options_frame_all.options_frame_layout.addWidget(self.energy_box, 2, 1)

        if current_tab == "Vina" or current_tab == "Smina" or current_tab == "ADFR":

            self.options_frame = OptionsFrame(parent=None,
        main_window=self.docking_programs_child_tabs)

            self.group_docking_layout.addWidget(self.options_frame, 3, 1, 5, 1)

            self.use_flex_vina_cb = QtWidgets.QCheckBox("Use Flex")
            self.use_flex_vina_cb.clicked.connect(self.use_flex_vina_cb_func)

            self.options_frame.options_frame_layout.addWidget(self.use_flex_vina_cb, 2, 0)

            # Options to select manually
            self.flex_chain_text = QtWidgets.QLabel("Chain")
            self.flex_chain_box = QtWidgets.QComboBox()
            self.flex_arg_edit = QtWidgets.QLineEdit()

            self.options_frame.options_frame_layout.addWidget(self.flex_arg_edit, 4, 0, 1, 2)
            self.flex_arg_edit.hide()

            # Options to select from PyMOL
            self.flex_import_btn = QtWidgets.QPushButton("Import Selection")
            self.flex_imported_sele = QtWidgets.QComboBox()
            self.flex_expand_text = QtWidgets.QLabel("Expand by")
            self.flex_combobox = QtWidgets.QComboBox()

            self.options_frame.options_frame_layout.addWidget(self.flex_import_btn, 4, 0)
            self.options_frame.options_frame_layout.addWidget(self.flex_imported_sele, 4, 1)
            self.options_frame.options_frame_layout.addWidget(self.flex_expand_text, 5, 0)
            self.options_frame.options_frame_layout.addWidget(self.flex_combobox, 5, 1)
            self.flex_import_btn.hide()
            self.flex_imported_sele.hide()
            self.flex_expand_text.hide()
            self.flex_combobox.hide()

            self.flex_import_btn.clicked.connect(self.import_sele_func)

            for i in [4,5,6,8,12,20]:
                i = str(i)
                self.flex_combobox.addItem(str(i + " A° residues"))

            self.choose_flex_manually = QtWidgets.QRadioButton("Select Manually")
            self.choose_flex_selection = QtWidgets.QRadioButton("Select from PyMOL Object")
            self.choose_flex_manually.setEnabled(False)
            self.choose_flex_selection.hide()
            self.choose_flex_selection.setEnabled(False)
            self.choose_flex_manually.hide()

            self.options_frame.options_frame_layout.addWidget(self.choose_flex_manually, 3, 0)
            self.options_frame.options_frame_layout.addWidget(self.choose_flex_selection, 3, 1)

            self.choose_flex_manually.clicked.connect(self.show_choose_flex_manually_options)
            self.choose_flex_selection.clicked.connect(self.show_choose_flex_selection_options)


        if current_tab == "Vina":
            self.scoring_label = QtWidgets.QLabel("Scoring Function")
            self.scoring_box = QtWidgets.QComboBox()
            self.scoring_box.addItems(["Standard", "Vinardo"])

            self.options_frame_all.options_frame_layout.addWidget(self.scoring_label, 3, 0)
            self.options_frame_all.options_frame_layout.addWidget(self.scoring_box, 3, 1)


        if current_tab == "Smina":

            self.rmsd_label = QtWidgets.QLabel("RMSD filter")
            self.rmsd_box = QtWidgets.QSpinBox()
            self.rmsd_box.setRange(1, 10)
            self.rmsd_box.setSingleStep(1)
            self.rmsd_box.setValue(1)
            self.options_frame_all.options_frame_layout.addWidget(self.rmsd_label, 3, 0)
            self.options_frame_all.options_frame_layout.addWidget(self.rmsd_box, 3, 1)

            self.buffer_label = QtWidgets.QLabel("Amount of Buffer")
            self.buffer_box = QtWidgets.QSpinBox()
            self.buffer_box.setRange(1, 10)
            self.buffer_box.setSingleStep(1)
            self.buffer_box.setValue(4)

            self.options_frame_all.options_frame_layout.addWidget(self.buffer_label, 4, 0)
            self.options_frame_all.options_frame_layout.addWidget(self.buffer_box, 4, 1)

            self.scoring_label = QtWidgets.QLabel("Scoring Function")
            self.scoring_box = QtWidgets.QComboBox()
            self.scoring_box.addItems(["Standard", "Vinardo"])

            self.options_frame_all.options_frame_layout.addWidget(self.scoring_label, 5, 0)
            self.options_frame_all.options_frame_layout.addWidget(self.scoring_box, 5, 1)


        if current_tab == "RxDock":

            self.protein_segments_to_exclude = []

            # OptionsFrame where general docking options are specified
            self.options_frame_all = OptionsFrame(parent=None,
        main_window=self.docking_programs_child_tabs)

            # OptionsFrame where docking options specific for a certain RxDock Docking Protocol are specified
            self.options_frame = OptionsFrame(parent=None,
        main_window=self.docking_programs_child_tabs)

            # Add Both OptionsFrame to the Docking Tab
            self.group_docking_layout.addWidget(self.options_frame, 2, 1, 6, 1)
            self.group_docking_layout.addWidget(self.options_frame_all, 0, 1, 2, 1)

            ### GENERAL DOCKING OPTIONS

            # Num poses
            self.poses_label = QtWidgets.QLabel("Poses")
            self.poses_box = QtWidgets.QSpinBox()
            self.poses_box.setRange(1, 100)
            self.poses_box.setSingleStep(1)
            self.poses_box.setValue(1)
            self.options_frame_all.options_frame_layout.addWidget(self.poses_label, 0, 0)
            self.options_frame_all.options_frame_layout.addWidget(self.poses_box, 0, 1)

            # Structural water
            self.receptor_water_checkbtn = QtWidgets.QCheckBox("Use Structural Water")
            self.options_frame_all.options_frame_layout.addWidget(self.receptor_water_checkbtn, 1, 0)

            # protein segments
            self.protein_segments_btn = QtWidgets.QPushButton("Exclude Heteroatoms")
            self.protein_segments_btn.setToolTip("To select the heteroatoms to exclude during Docking process")
            #self.options_frame_all.options_frame_layout.addWidget(self.protein_segments_btn, 2, 0)
            self.protein_segments_btn.clicked.connect(self.show_protein_segments_window)

            # Docking Alternative Options RadioButtons
            self.use_pharma_restrain_cb = QtWidgets.QCheckBox("Pharmacophoric restrains")
            self.use_pharma_restrain_cb.toggled.connect(self.alternative_docking_func)
            # self.docking_alternatives_combobox = QtWidgets.QComboBox()
            # self.docking_alternatives_combobox.addItems(["Simple Docking", "Pharmacophoric restrains"])
            # self.docking_alternatives_combobox.currentTextChanged.connect(self.alternative_docking_func)
            self.options_frame.options_frame_layout.addWidget(self.use_pharma_restrain_cb, 0, 0, 1, 2)

            ### PHARMACOPHORIC RESTRAINS WIDGETS

            self.pharma_widgets_list = [] # A list used to show and hide easily the widgets

            self.open_pharma_file_btn = QtWidgets.QPushButton("Open From File")
            self.options_frame.options_frame_layout.addWidget(self.open_pharma_file_btn, 1, 0)
            self.pharma_widgets_list.append(self.open_pharma_file_btn)
            self.open_pharma_file_btn.clicked.connect(self.open_pharma_file_func)

            self.set_manually_btn = QtWidgets.QPushButton("Set Manually")
            self.options_frame.options_frame_layout.addWidget(self.set_manually_btn, 2, 0)
            self.pharma_widgets_list.append(self.set_manually_btn)
            self.set_manually_btn.clicked.connect(self.show_pharma_advanced_options)

            self.pharma_list = QtWidgets.QListWidget()
            self.options_frame.options_frame_layout.addWidget(self.pharma_list, 3, 0)
            self.pharma_widgets_list.append(self.pharma_list)

            self.remove_pharma_btn = QtWidgets.QPushButton("Remove")
            self.remove_pharma_btn.clicked.connect(self.remove_pharma_func)
            self.options_frame.options_frame_layout.addWidget(self.remove_pharma_btn, 4, 0)
            self.pharma_widgets_list.append(self.remove_pharma_btn)

            for widgets in self.pharma_widgets_list:
                widgets.hide()

            ### TETHERED DOCKING WIDGETS

            self.tethered_widgets_list = [] # A list used to show and hide easily the widgets

            self.import_tethered_ligand = QtWidgets.QPushButton("Import from PyMOL")
            self.options_frame.options_frame_layout.addWidget(self.import_tethered_ligand, 1, 0)
            self.tethered_widgets_list.append(self.import_tethered_ligand)
            self.import_tethered_ligand.clicked.connect(self.import_tethered_ligand_func)

            self.tethered_ligand_combo = QtWidgets.QComboBox()
            self.options_frame.options_frame_layout.addWidget(self.tethered_ligand_combo, 1, 1)
            self.tethered_widgets_list.append(self.tethered_ligand_combo)

            self.trans_mode = QtWidgets.QLabel("Position")
            self.trans_mode_combo = QtWidgets.QComboBox()
            self.tethered_widgets_list.append(self.trans_mode)
            self.tethered_widgets_list.append(self.trans_mode_combo)

            self.options_frame.options_frame_layout.addWidget(self.trans_mode, 2, 0)
            self.options_frame.options_frame_layout.addWidget(self.trans_mode_combo, 2, 1)

            self.rot_mode = QtWidgets.QLabel("Rotation")
            self.rot_mode_combo = QtWidgets.QComboBox()
            self.tethered_widgets_list.append(self.rot_mode)
            self.tethered_widgets_list.append(self.rot_mode_combo)

            self.options_frame.options_frame_layout.addWidget(self.rot_mode, 3, 0)
            self.options_frame.options_frame_layout.addWidget(self.rot_mode_combo, 3, 1)

            self.die_mode = QtWidgets.QLabel("Dihedral")
            self.die_mode_combo = QtWidgets.QComboBox()
            self.tethered_widgets_list.append(self.die_mode)
            self.tethered_widgets_list.append(self.die_mode_combo)

            self.options_frame.options_frame_layout.addWidget(self.die_mode, 4, 0)
            self.options_frame.options_frame_layout.addWidget(self.die_mode_combo, 4, 1)

            for widgets in [self.trans_mode_combo, self.rot_mode_combo, self.die_mode_combo]:
                widgets.addItems(["FREE", "TETHERED", "FIXED"])

            for widgets in self.tethered_widgets_list:
                widgets.hide()


    def import_tethered_ligand_func(self):

        # Create a list of all the importable Objects and selections from PyMOL
        self.selections = [str(obj) for obj in cmd.get_names("objects") + cmd.get_names("selections")]

        # Check for importable objects
        if self.selections == []:
            QtWidgets.QMessageBox.warning(self.docking_programs_child_tabs, "PyMOL is empty", str("There isn't any object to import"))

        # If importable objects are present, update the ComboBox
        else:
            self.tethered_ligand_combo.clear()
            for i in self.selections:
                type = cmd.get_type(i)
                if type == str("object:molecule"):
                    self.tethered_ligand_combo.addItem(i)
                elif type == str("selection"):
                    self.tethered_ligand_combo.addItem(i)
                else:
                    pass


    def show_protein_segments_window(self):

        HandleWidgets.combobox_check_if_empty(self = self,
        widgets_list = [self.ready_receptors_list])

        if self.is_empty:
            QtWidgets.QMessageBox.warning(self.docking_programs_child_tabs, "Warning", "There aren't enough parameters.\n Please generate a Receptor Object")

        else:
            self.protein_segments_window = NewWindow(parent = self.docking_programs_child_tabs,
            title = "Advanced Options", upper_frame_title = "Select Heteroatoms to Exclude",
            submit_command = self.exclude_protein_segments, submit_button_text= "Set",
            with_scroll = True)

            # Create options to select protein_segments
            self.label = QtWidgets.QLabel()
            self.label.setText(self.ready_receptors_list.currentText())

            name = self.ready_receptors_list.currentText().split("_")[1]
            het_list = self.docking_programs_child_tabs.docking_programs.rxdock_receptors_dict[name]["parsed_object"].heteroresidues_list

            self.checkbtn_list = []
            row = 1
            for i in het_list:
                self.checkbtn = QtWidgets.QCheckBox(str(i))
                self.protein_segments_window.middle_layout_type.addWidget(self.checkbtn, row, 0)
                self.checkbtn_list.append(self.checkbtn)
                row += 1

            self.protein_segments_window.show()


    def exclude_protein_segments(self):

        for btn in self.checkbtn_list:
            if btn.isChecked():
                self.protein_segments_to_exclude.append(btn.text())

        self.protein_segments_window.hide()


    def show_pharma_advanced_options(self):

        self.pharma_advanced_options_window = NewWindow(parent = self.docking_programs_child_tabs,
        title = "Advanced Options", upper_frame_title = "Set Options to add Pharmacophoric Restrains",
        submit_command = self.set_pharma_btn_func, submit_button_text= "Set",
        with_scroll = True)

        # Create Options to Set Pharmacophoric Restrains
        self.import_label = QtWidgets.QPushButton("Import Atomic Selection")
        self.import_label.clicked.connect(self.import_atomic_from_pymol)
        self.import_label_box = QtWidgets.QComboBox()
        self.import_label_box.currentTextChanged.connect(self.zoom_sele)
        self.import_label_box.view().pressed.connect(self.zoom_sele)

        self.pharma_widgets_list.append(self.import_label)
        self.pharma_widgets_list.append(self.import_label_box)

        self.tolerance_radius_label = QtWidgets.QLabel("Tolerance Radius")
        self.tolerance_radius_box = QtWidgets.QDoubleSpinBox() ## TODO --- QUALI VALORI POSSIBILI??
        self.pharma_widgets_list.append(self.tolerance_radius_label)
        self.pharma_widgets_list.append(self.tolerance_radius_box)

        self.restrain_type_label = QtWidgets.QLabel("Restrain Type")
        self.restrain_type_box = QtWidgets.QComboBox()
        self.pharma_widgets_list.append(self.restrain_type_label)
        self.pharma_widgets_list.append(self.restrain_type_box)

        restrain_types_list = ["Any", "Don", "Acc", "Aro", "Hyd", "Hal", "Har", "Ani", "Cat"]
        self.restrain_type_box.addItems(restrain_types_list)

        # Add Options to the New Window
        self.pharma_advanced_options_window.middle_layout_type.addWidget(self.import_label, 0, 0)
        self.pharma_advanced_options_window.middle_layout_type.addWidget(self.import_label_box, 0, 1)

        self.pharma_advanced_options_window.middle_layout_type.addWidget(self.tolerance_radius_label, 1, 0)
        self.pharma_advanced_options_window.middle_layout_type.addWidget(self.tolerance_radius_box, 1, 1)

        self.pharma_advanced_options_window.middle_layout_type.addWidget(self.restrain_type_label, 2, 0)
        self.pharma_advanced_options_window.middle_layout_type.addWidget(self.restrain_type_box, 2, 1)

        self.pharma_advanced_options_window.show()


    ## VINA SPECIFIC DOCKING FUNCTIONS

    def import_sele_func(self):

        # self.selections = [str(obj) for obj in cmd.get_names("objects") + cmd.get_names("selections")]
        self.selections = [str(obj) for obj in cmd.get_names("selections")]
        self.flex_imported_sele.clear()
        self.flex_imported_sele.addItems(self.selections)

        # import inspect
        # print(inspect.getargspec(cmd.select))
        # #cmd.select(self.selections[0], around = 4.5)

    def use_flex_func(self):

        if self.use_flex_vina_cb.isChecked():
            self.flex_arg_edit.show()
        else:
            self.flex_arg_edit.hide()


    def show_choose_flex_manually_options(self):

        if self.choose_flex_manually.isChecked():
            self.flex_arg_edit.show()
            self.flex_arg_edit.setPlaceholderText("CHAIN:RESNAMEPOSITION (example: A:ARG220)")

            self.flex_import_btn.hide()
            self.flex_imported_sele.hide()
            self.flex_expand_text.hide()
            self.flex_combobox.hide()


    def show_choose_flex_selection_options(self):

        if self.choose_flex_selection.isChecked():
            self.flex_arg_edit.hide()

            self.flex_import_btn.show()
            self.flex_imported_sele.show()
            self.flex_expand_text.show()
            self.flex_combobox.show()

    def use_flex_vina_cb_func(self):

        if self.use_flex_vina_cb.isChecked():
            self.choose_flex_manually.setChecked(True)
            self.choose_flex_manually.setEnabled(True)
            self.choose_flex_selection.setEnabled(True)
            self.show_choose_flex_manually_options()

        else:
            self.choose_flex_manually.setEnabled(False)
            self.choose_flex_selection.setEnabled(False)
            self.flex_arg_edit.hide()
            self.flex_import_btn.hide()
            self.flex_imported_sele.hide()
            self.flex_expand_text.hide()
            self.flex_combobox.hide()


    ## RXDOCK SPECIFIC DOCKING FUNCTIONS

    def open_pharma_file_func(self):

        # Open a .const File
        opened_file = OpenFromFile(self,
        file_type_string = "Pharma Restrains File (*.const)")

        if opened_file.is_valid and opened_file.is_already_loaded == False:

            # Change dir --> tmp_dir
            os.chdir(self.docking_programs_child_tabs.docking_programs.rxdock_tmp_dir)

            input = open(str(opened_file.file_name + ".const"), "rt")

            for line in input:
                if line == "\n":
                    pass
                else:
                    self.pharma_list.addItem(str(line))

            input.close()


            self.docking_programs_child_tabs.docking_programs.main_window.statusBar().showMessage(str("Loaded file '" + opened_file.file_name + "'"), 3000)


    def remove_pharma_func(self):

        for items in self.pharma_list.selectedItems():
            qm = QtWidgets.QMessageBox.question(self,'', str("Are you sure you want to remove " + items.text()), QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

            if qm == QtWidgets.QMessageBox.Yes:
                self.pharma_list.takeItem(self.pharma_list.indexFromItem(items).row())
            else:
                pass


    def alternative_docking_func(self):

        if self.use_pharma_restrain_cb.isChecked():
            self.pharma_restrains = True
            self.show_widgets([self.pharma_widgets_list])
        else:
            self.pharma_restrains = False
            self.hide_widgets([self.pharma_widgets_list])

        # if self.docking_alternatives_combobox.currentText() == "Simple Docking":
        #     self.simple_docking = True
        #     self.hide_widgets([self.pharma_widgets_list, self.tethered_widgets_list])
        #
        # elif self.docking_alternatives_combobox.currentText() == "Pharmacophoric restrains":
        #     self.pharma_restrains = True
        #     self.show_widgets([self.pharma_widgets_list])
        #     self.hide_widgets([self.tethered_widgets_list])
        #
        # elif self.docking_alternatives_combobox.currentText() == "Tethered Docking":
        #     self.tethered_docking = True
        #     self.hide_widgets([self.pharma_widgets_list])
        #     self.show_widgets([self.tethered_widgets_list])
        #

    def hide_widgets(self, list):

        for i in list:
            for items in i:
                items.hide()


    def show_widgets(self, list):

        for i in list:
            for items in i:
                items.show()


    def set_pharma_btn_func(self):

        coords = (self.import_label_box.currentText()).replace("sele ", "")
        tolerance_radius = str(self.tolerance_radius_box.value())
        restrain_type = self.restrain_type_box.currentText()

        HandleWidgets.combobox_check_if_empty(self = self,
        widgets_list = [self.import_label_box, self.restrain_type_box])

        if self.is_empty:
            QtWidgets.QMessageBox.warning(self.docking_programs_child_tabs, "There isn't any specified coordinate", "Please import an atomic \nselection from PyMOL")

        else:
            to_add = str(coords + " " + tolerance_radius + " " + restrain_type)
            self.pharma_list.addItem(to_add)


    def zoom_sele(self):

        current_text = self.import_label_box.currentText()
        current_text_pymol_name = current_text.split(" ")
        cmd.zoom(current_text_pymol_name[0])


    def import_atomic_from_pymol(self):

         #TODO - IN PYMOL INTERACTIONS
        tmp_list = [str(obj) for obj in cmd.get_names("selections")]

        for i in tmp_list:
            type = cmd.get_type(i)
            if type == str("object:molecule") or str("object:selection"):
                count = cmd.count_atoms(i)
                if count == 1:
                    coords = cmd.get_coords(str(i), 1)
                    x = str(coords[0][0])
                    y = str(coords[0][1])
                    z = str(coords[0][2])
                    self.import_label_box.clear()
                    self.import_label_box.addItem(str(i + " "+ x + " " + y + " " + z + " "))


    def loaded_cavities_show_func(self):

        self.cavity_to_show = self.loaded_cavities.currentText()

        try:
            self.show_crisscross(self.cavity_to_show)
            self.calculate_box(self.cavity_to_show)
        except:
            pass


    def initialize_docking_inputs(self):

        self.dockings_to_do = []

        self.ligands_to_dock = []
        self.receptors_to_dock = []

        if self.protocol_box.currentText() == "Custom":
            for i in reversed(range(self.group_scroll_layout.count())):
                cb = self.group_scroll_layout.itemAt(i).widget()

                if cb.isChecked():
                    tmp_dict = self.dockings_to_do_dict[cb.title()]
                    # tmp_list = []
                    # for j in reversed(range(self.docking_group_layout.count())):
                    #     print(self.docking_group_layout.itemAt(j).widget())
                    #     text = self.docking_group_layout.itemAt(j).widget().text()
                    #     print(text)
                    #
                    #     if re.search("Receptor", text):
                    #         receptor = text.split()[1]
                    #         tmp_list.extend([receptor])
                    #
                    #     if re.search("Ligand", text):
                    #         ligand = text.split()[1]
                    #         tmp_list.extend([ligand])
                    #
                    #     if re.search("Cavity:", text):
                    #         cavity = text.replace("Cavity: ", "")
                    #         tmp_list.extend([cavity])

                    self.dockings_to_do.append(tmp_dict)

            print(self.dockings_to_do)
            if not self.is_empty:
                self.ligands_to_dock.append(self.ready_ligands_list.currentText())
                self.receptors_to_dock.append(self.ready_receptors_list.currentText())

        if self.protocol_box.currentText() == "ALLvsALL":

            self.dockings_to_do = []

            self.cavity = self.loaded_cavities.currentText()

            for i in reversed(range(self.ligand_scroll_layout.count())):
                cb = self.ligand_scroll_layout.itemAt(i).widget()
                if cb.isChecked():
                    self.ligands_to_dock.append(cb.text())

            for i in reversed(range(self.receptors_scroll_layout.count())):
                cb = self.receptors_scroll_layout.itemAt(i).widget()
                if cb.isChecked():
                    self.receptors_to_dock.append(cb.text())


    def show_resume_window(self):

        if self.protocol_box.currentText() == "Custom":
            number_of_dockings_to_do = len(self.dockings_to_do)

        if self.protocol_box.currentText() == "ALLvsALL":
            number_of_dockings_to_do = int(len(self.receptors_to_dock) * len(self.ligands_to_dock))

        ####
        # Show Docking Dialog from within Thread is initialized
        ####

        self.docking_dialog = Dockings_dialog(self,
        number_of_dockings_to_do = number_of_dockings_to_do,
        receptors_to_dock = self.receptors_to_dock,
        ligands_to_dock = self.ligands_to_dock,
        dockings_to_do = self.dockings_to_do)

        self.docking_dialog.setModal(True)
        self.docking_dialog.exec_()

        # Finishes Docking Process.
        if self.docking_dialog.complete_status:
            pass

    def run_docking_func(self):

        # Check the Docking Program in use
        Check_current_tab.check_docking_program_current_tab(self)

        self.initialize_docking_inputs()

        dict = ConfigurationTab.check_installation(self = self.docking_programs_child_tabs, check_only = True)

        for key in dict:
            if dict[key] == False:
                if sys.platform == "win32":
                    if str(key) == "sdsorter":
                        pass
                    elif str(key) == "Smina":
                        pass
                    elif str(key) == "RxDock":
                        pass
                    else:
                        print(str(key) + " is not installed. Please check CONFIGURATION Tab")

                elif sys.platform == "linux":
                    if str(key) == "sdsorter":
                        pass
                    else:
                        print(str(key) + " is not installed. Please check CONFIGURATION Tab")

        # Check if some parameters are missing
        empty_dict = HandleWidgets.combobox_check_if_empty(self=self,
        widgets_list=[self.ready_ligands_list, self.ready_receptors_list, self.loaded_cavities])

        if self.receptors_to_dock == [] and self.ligands_to_dock == []:
            self.is_empty = True
            text = "Missing Receptor and Ligand"

        elif self.receptors_to_dock == []:
            self.is_empty = True
            text = "Missing Receptor"

        elif self.ligands_to_dock == []:
            self.is_empty = True
            text = "Missing Ligand"

        elif self.is_empty:
            text = "Missing Parameters"


        if self.is_empty:
            QtWidgets.QMessageBox.warning(self.docking_programs_child_tabs, "Warning", text)
        else:
            self.tab_docking_runs_scroll_layout = self.docking_programs_child_tabs.docking_programs.data_analysis_layout.docking_runs_scroll_layout

            self.show_resume_window()

    #
    # def starting_docking_process(self, receptor, ligand, cavity = None):
    #
    #     self.docking_completed = False
    #
    #     docking_programs = self.docking_programs_child_tabs.docking_programs
    #
    #     ## VINA DOCKING SETUP
    #     if self.is_vina_tab:
    #
    #         self.last_docking = Vina_docking(self,
    #         ligand = ligand,
    #         receptor = receptor)
    #
    #         if self.last_docking.interrupt == False:
    #
    #             self.check_if_docking_completed(runs = docking_programs.vina_runs,
    #             tmp_dir = docking_programs.vina_tmp_dir)
    #
    #             if self.docking_completed:
    #
    #                 self.results_file = Vina_Parse_Results(self, results_file_name = self.last_docking.results_file_name,
    #                 results_dict = self.docking_programs_child_tabs.docking_programs.results_dict,
    #                 poses = self.last_docking.poses)
    #
    #             self.summary_file = self.write_summary_file()
    #
    #             # Create a result tab for each Docking Run in RxDock Docking Tab
    #             results_frame = self.create_results_frame(results_tab = self.docking_programs_child_tabs.docking_programs.VINA.results_tab_ui,
    #             docking_process = self.last_docking,
    #             program = "Vina",
    #             columns_names = ["NAME", "POSE", "Affinity\n(kcal/mol)", "dist from best mode\nRMSD l.b.", "dist from best mode\nRMSD u.b."])
    #
    #             self.docking_programs_child_tabs.docking_programs.VINA.results_tab_ui.save_all_files.setEnabled(True)
    #             self.docking_programs_child_tabs.docking_programs.VINA.results_tab_ui.all_plot_btn.setEnabled(True)
    #
    #
    #             self.docking_programs_child_tabs.docking_programs.vina_runs_dict[self.last_docking.results_file_name] = {}
    #             self.docking_programs_child_tabs.docking_programs.vina_runs_dict[self.last_docking.results_file_name]["docking_run"] = self.last_docking
    #             self.docking_programs_child_tabs.docking_programs.vina_runs_dict[self.last_docking.results_file_name]["program"] = "Vina"
    #             self.docking_programs_child_tabs.docking_programs.vina_runs_dict[self.last_docking.results_file_name]["results_table"] = results_frame.results_table
    #
    #
    #     ## RXDOCK DOCKING SETUP
    #     if self.is_rxdock_tab:
    #
    #         self.get_rxdock_docking_protocol()
    #
    #         # Run RxDock Docking
    #         self.last_docking = RxDock_docking(self, ligand, receptor)
    #
    #         if self.last_docking.interrupt == False:
    #
    #             self.check_if_docking_completed(runs = docking_programs.rxdock_runs,
    #             tmp_dir = docking_programs.rxdock_tmp_dir)
    #
    #             if self.docking_completed:
    #
    #                 self.results_file = RxDock_parse_results(self, results_file_name = self.last_docking.results_file_name,
    #                 results_dict = self.docking_programs_child_tabs.docking_programs.results_dict,
    #                 poses = self.last_docking.poses,
    #                 ligand = ligand)
    #
    #             self.summary_file = self.write_summary_file()
    #
    #             # Create a result tab for each Docking Run in RxDock Docking Tab
    #             results_frame = self.create_results_frame(results_tab = self.docking_programs_child_tabs.docking_programs.RXDOCK.results_tab_ui,
    #             docking_process = self.last_docking,
    #             program = "RxDock",
    #             columns_names = ["NAME", "POSE", "SCORE", "SCORE-inter", "SCORE-intra"])
    #
    #             self.docking_programs_child_tabs.docking_programs.RXDOCK.results_tab_ui.save_all_files.setEnabled(True)
    #             self.docking_programs_child_tabs.docking_programs.RXDOCK.results_tab_ui.all_plot_btn.setEnabled(True)
    #
    #             self.docking_programs_child_tabs.docking_programs.rxdock_runs_dict[self.last_docking.results_file_name] = {}
    #             self.docking_programs_child_tabs.docking_programs.rxdock_runs_dict[self.last_docking.results_file_name]["docking_run"] = self.last_docking
    #             self.docking_programs_child_tabs.docking_programs.rxdock_runs_dict[self.last_docking.results_file_name]["program"] = "RxDock"
    #             self.docking_programs_child_tabs.docking_programs.rxdock_runs_dict[self.last_docking.results_file_name]["results_table"] = results_frame.results_table
    #
    #
    #     if self.is_smina_tab:
    #
    #         self.last_docking = Smina_docking(self, ligand, receptor)
    #
    #         if self.last_docking.interrupt == False:
    #
    #             self.check_if_docking_completed(runs = docking_programs.smina_runs,
    #             tmp_dir = docking_programs.smina_tmp_dir)
    #
    #             if self.docking_completed:
    #
    #                 self.results_file = Smina_parse_results(self, results_file_name = self.last_docking.results_file_name,
    #                 results_dict = self.docking_programs_child_tabs.docking_programs.results_dict,
    #                 poses = self.last_docking.poses,
    #                 ligand = ligand)
    #
    #             self.summary_file = self.write_summary_file()
    #
    #             # Create a result tab for each Docking Run in RxDock Docking Tab
    #             results_frame = self.create_results_frame(results_tab = self.docking_programs_child_tabs.docking_programs.SMINA.results_tab_ui,
    #             docking_process = self.last_docking, program = "Smina", columns_names = ["NAME", "POSE", "Affinity (kcal/mol)"])
    #
    #             self.docking_programs_child_tabs.docking_programs.SMINA.results_tab_ui.save_all_files.setEnabled(True)
    #             self.docking_programs_child_tabs.docking_programs.SMINA.results_tab_ui.all_plot_btn.setEnabled(True)
    #
    #             self.docking_programs_child_tabs.docking_programs.smina_runs_dict[self.last_docking.results_file_name] = {}
    #             self.docking_programs_child_tabs.docking_programs.smina_runs_dict[self.last_docking.results_file_name]["docking_run"] = self.last_docking
    #             self.docking_programs_child_tabs.docking_programs.smina_runs_dict[self.last_docking.results_file_name]["program"] = "Smina"
    #             self.docking_programs_child_tabs.docking_programs.smina_runs_dict[self.last_docking.results_file_name]["results_table"] = results_frame.results_table
    #
    #     if self.is_adfr_tab:
    #
    #         self.last_docking = ADFR_docking(self, ligand, receptor)
    #
    #         if self.last_docking.interrupt == False:
    #
    #             self.check_if_docking_completed(runs = docking_programs.adfr_runs,
    #             tmp_dir = docking_programs.adfr_tmp_dir)
    #
    #             if self.docking_completed:
    #
    #                 self.results_file = ADFR_parse_results(self, results_file_name = self.last_docking.results_file_name,
    #                 results_dict = self.docking_programs_child_tabs.docking_programs.results_dict,
    #                 ligand = ligand)
    #
    #             self.summary_file = self.write_summary_file()
    #
    #             # Create a result tab for each Docking Run in RxDock Docking Tab
    #             results_frame = self.create_results_frame(results_tab = self.docking_programs_child_tabs.docking_programs.ADFR.results_tab_ui,
    #             docking_process = self.last_docking, program = "ADFR", columns_names = ["NAME", "POSE", "Affinity (kcal/mol)"])
    #
    #             self.docking_programs_child_tabs.docking_programs.ADFR.results_tab_ui.save_all_files.setEnabled(True)
    #             self.docking_programs_child_tabs.docking_programs.ADFR.results_tab_ui.all_plot_btn.setEnabled(True)
    #
    #             self.docking_programs_child_tabs.docking_programs.adfr_runs_dict[self.last_docking.results_file_name] = {}
    #             self.docking_programs_child_tabs.docking_programs.adfr_runs_dict[self.last_docking.results_file_name]["docking_run"] = self.last_docking
    #             self.docking_programs_child_tabs.docking_programs.adfr_runs_dict[self.last_docking.results_file_name]["program"] = "ADFR"
    #             self.docking_programs_child_tabs.docking_programs.adfr_runs_dict[self.last_docking.results_file_name]["results_table"] = results_frame.results_table
    #
    #     if self.docking_completed:
    #
    #         # Create checkbox in DATA ANALYSIS tab
    #         self.dockings_frame = DockingsFrame(parent=None, main_window=self, results_name=self.last_docking.results_file_name)
    #         self.docking_programs_child_tabs.docking_programs.data_analysis_layout.docking_runs_scroll_layout.addRow(self.dockings_frame)
    #
    #         # Update dict with info about Docking runs
    #         self.docking_programs_child_tabs.docking_programs.all_runs[self.last_docking.results_file_name] = {}
    #         self.docking_programs_child_tabs.docking_programs.all_runs[self.last_docking.results_file_name]["docking_run"] = self.last_docking
    #         self.docking_programs_child_tabs.docking_programs.all_runs[self.last_docking.results_file_name]["docking_run_results_file"] = self.results_file
    #         self.docking_programs_child_tabs.docking_programs.all_runs[self.last_docking.results_file_name]["dockings_frame"] = self.dockings_frame


    def check_if_docking_completed(self, runs, tmp_dir):

        file_path = os.path.join(tmp_dir, self.last_docking.results_file_name_ext)

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
                    QtWidgets.QMessageBox.warning(self.docking_programs_child_tabs, "", str("Something went wrong during Docking. \nPlease check LOG files."))

        else:
            self.docking_completed = False
            QtWidgets.QMessageBox.warning(self.docking_programs_child_tabs, "", str("Something went wrong during Docking. \nPlease check LOG files."))
            runs += 1


    def write_summary_file(self):

        Check_current_tab.check_docking_program_current_tab(self)

        summary_list = []

        #states = cmd.count_states(self.last_docking.results_file_name)

        ligand = str("Ligand: " + self.last_docking.ligand_to_dock)
        receptor = str("Receptor: " + self.last_docking.receptor_to_dock)
        cavity = str("Cavity: " + self.last_docking.cavity_name)

        if self.is_smina_tab and self.last_docking.grid is None:
            cavity = str("Cavity: " + self.last_docking.cavity_name + " with buffer " + str(self.last_docking.buffer))
        #poses = str("Generated Poses: " + str(states))

        if self.is_adfr_tab:
            ga_evol = str("Number of GA evolutions: " + self.last_docking.ga_evol)
            ga_threshold = str("GA evaluation stoped after " + self.last_docking.ga_threshold + " \nnumber of generations with no improvement in best energy in all clusters")
            max_gen = str("Maximum number of generations: " + self.last_docking.max_gen)


        if self.is_vina_tab or self.is_smina_tab:
            exhaustiveness = str("Exhaustiveness: " + str(self.last_docking.exhaustiveness))
            energy = str("Energy Range: " + str(self.last_docking.energy))

            summary_list.extend([exhaustiveness, energy])

            if self.last_docking.use_flex_protocol:
                flex = str("Flexible side chains: " + str(self.last_docking.flex_residues))

                summary_list.extend([flex])

        summary_list.extend([ligand, receptor, cavity])

        textfilename = str(self.last_docking.results_file_name + "_summary.txt")
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
                    self.docking_programs_child_tabs.docking_programs.results_dict[i] = {}

            else:
                self.docking_programs_child_tabs.docking_programs.results_dict[self.results_file_name] = {}

            self.results_frame = ResultsFrame(parent=None,
            main_window=self.docking_programs_child_tabs,
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

            self.results_obj_name = docking_process.results_file_name
            self.results_file_name = docking_process.results_file_name_ext
            self.log_file_name = docking_process.log_file_name

            self.docking_programs_child_tabs.docking_programs.results_dict[self.results_file_name] = {}

            self.results_frame = ResultsFrame(parent=None,
            main_window=self.docking_programs_child_tabs,
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


    # def check_docking_program_current_tab(self):
    #
    #     a = self.parent().main_docking_programs_tabs.currentIndex()
    #     current_tab = self.parent().main_docking_programs_tabs.tabText(a)
    #
    #     self.is_vina_tab = False
    #     self.is_rxdock_tab = False
    #
    #     if current_tab == "Vina":
    #         self.is_vina_tab = True
    #
    #     elif current_tab == "RxDock":
    #         self.is_rxdock_tab = True


    def get_rxdock_docking_protocol(self):

        # Set RxDock options
        self.simple_docking = False
        self.pharma_restrains = False
        self.tethered_docking = False

        if self.use_pharma_restrain_cb.isChecked():
            self.pharma_restrains = True

        # if self.docking_alternatives_combobox.currentText() == "Simple Docking":
        #     self.simple_docking = True
        #
        # elif self.docking_alternatives_combobox.currentText() == "Pharmacophoric restrains":
        #     self.pharma_restrains = True
        #
        # elif self.docking_alternatives_combobox.currentText() == "Tethered Docking":
        #     self.tethered_docking = True



class ResultsAdvancedOptions(QtWidgets.QMainWindow):


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



class ResultsTab(QtWidgets.QWidget, PyMOLInteractions):

    def __init__(self, main_window):
        super().__init__(main_window)
        self.style = "background-color: rgb(0, 0, 0); color: rgb(255, 255, 255); font-weight: bold"
        self.docking_programs_child_tabs = main_window
        self.initUI()

    def initUI(self):

        self.layout_tab6 = QtWidgets.QGridLayout()

        # Group for results
        self.results_group_layout = QtWidgets.QGridLayout()
        self.results_group = QtWidgets.QGroupBox("Results")
        self.results_group.setLayout(self.results_group_layout)
        self.layout_tab6.addWidget(self.results_group)

        # Widgets to work with all runs
        self.save_all_files = QtWidgets.QPushButton("Save All Runs")
        self.results_group_layout.addWidget(self.save_all_files, 1, 0)
        self.save_all_files.clicked.connect(self.save_all_files_func)
        self.results_group_layout.setAlignment(self.save_all_files, QtCore.Qt.AlignBottom)
        self.save_all_files.setEnabled(False)
        self.save_all_files.hide()

        # Plot all (out of OptionsFrame)
        self.all_plot_btn = QtWidgets.QPushButton("Plot All")
        self.all_plot_btn.clicked.connect(self.plot_run_score_scatter)
        #results_frame_layout
        self.results_group_layout.addWidget(self.all_plot_btn, 2, 0)
        self.results_group_layout.setAlignment(self.all_plot_btn, QtCore.Qt.AlignBottom)
        self.all_plot_btn.setEnabled(False)
        self.all_plot_btn.hide()

        # Widgets to add in the Results Group Box
        # self.plot_distribution_btn = QtWidgets.QPushButton("Plot distribution")
        # # self.results_group_layout.addWidget(self.plot_distribution_btn)
        # self.plot_distribution_btn.clicked.connect(self.show_plot_distribution_window)
        # self.plot_distribution_btn.hide()

        # Tab

        self.result_tabs = QtWidgets.QTabWidget()


    def plot_run_score_scatter(self):

        Check_current_tab.check_docking_program_current_tab(self)

        if self.is_vina_tab:
            self.dict = self.docking_programs_child_tabs.docking_programs.vina_runs_dict
            self.score_header = "Affinity"

        if self.is_smina_tab:
            self.dict = self.docking_programs_child_tabs.docking_programs.smina_runs_dict
            self.score_header = "Affinity"

        if self.is_rxdock_tab:
            self.dict = self.docking_programs_child_tabs.docking_programs.rxdock_runs_dict
            self.score_header = "SCORE"

        if self.is_adfr_tab:
            self.dict = self.docking_programs_child_tabs.docking_programs.adfr_runs_dict
            self.score_header = "Affinity"

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
            if table is not None:
                # Get index of "Score/Affinity" column
                y_index = TableView.get_column_index_from_header(self, table = table, header = self.score_header)
                # Get data from index
                y = TableView.get_column_data_from_index(self, table = table, column_header_index = y_index[0])
                # Create repeats for each run. Workaround for discrete variables
                num = list(itertools.repeat(i, len(y)))
                # Add a ScatterPlot for each Docking Run
                cp.add_scatter_plot(num, y, label = key, corr=False)

        cp.show()


    def save_all_files_func(self):

        self.save_to_file_choise_window = NewWindow(parent = self.docking_programs_child_tabs,
        title = "Save to File", upper_frame_title = "Select the saving option",
        submit_command = self.get_saving_choice, submit_button_text= "Start",
        with_scroll = True)

        self.save_table_to_csv = QtWidgets.QRadioButton("Save table to \".csv\" file")
        self.save_docking_files = QtWidgets.QRadioButton("Save Docking results file")
        self.save_log = QtWidgets.QRadioButton("Save Log File")

        # Add Options to the New Window
        #self.save_to_file_choise_window.middle_layout_type.addWidget(self.save_table_to_csv, 0, 0)
        self.save_to_file_choise_window.middle_layout_type.addWidget(self.save_docking_files, 1, 0)
        self.save_to_file_choise_window.middle_layout_type.addWidget(self.save_log, 2, 0)

        self.save_to_file_choise_window.show()


    def get_saving_choice(self):

        Check_current_tab.check_docking_program_current_tab(self)

        self.docking_results_to_save_list = []
        self.log_files_to_save_list = []

        if self.is_vina_tab:
            tmp_dir = self.docking_programs_child_tabs.docking_programs.vina_tmp_dir

        if self.is_rxdock_tab:
            tmp_dir = self.docking_programs_child_tabs.docking_programs.rxdock_tmp_dir

        if self.is_smina_tab:
            tmp_dir = self.docking_programs_child_tabs.docking_programs.smina_tmp_dir

        if self.is_adfr_tab:
            tmp_dir = self.docking_programs_child_tabs.docking_programs.adfr_tmp_dir

        temp3 = None
        for element in os.listdir(tmp_dir):
            temp = re.search("_log", element)
            temp2 = re.search("Run_", element)
            temp3 = re.search(".txt", element)
            temp4 = re.search("docked", element)
            temp5 = re.search("reference", element)

            if temp:
                self.log_files_to_save_list.append(element)
            if temp2 and temp3 is None:
                if temp4 or temp5:
                    pass
                else:
                    self.docking_results_to_save_list.append(element)

        tmp_path = os.path.join(tmp_dir, "temporary")
        os.mkdir(tmp_path)

        if self.save_docking_files.isChecked():

            for file in self.docking_results_to_save_list:
                shutil.copy(os.path.join(tmp_dir, file), os.path.join(tmp_path, file))

            filepath = asksaveasfile_qt("Save Docking Results", name_filter="*.zip")

            if not filepath:
                shutil.rmtree(tmp_path)
                return None

            else:
                os.chdir(tmp_path)
                split_filepath = os.path.splitext(filepath)[0]
                shutil.make_archive(split_filepath, 'zip')


        # if self.save_table_to_csv.isChecked():
        #
        #     Save_to_Csv.save_to_csv_event(self, table = self.results_table)


        if self.save_log.isChecked():

            for file in self.log_files_to_save_list:
                shutil.copy(os.path.join(tmp_dir, file), os.path.join(tmp_path, file))

            filepath = asksaveasfile_qt("Save LOG Files", name_filter="*.zip")

            if not filepath:
                shutil.rmtree(tmp_path)
                return None

            else:
                os.chdir(tmp_path)
                split_filepath = os.path.splitext(filepath)[0]
                shutil.make_archive(split_filepath, 'zip')

        os.chdir(self.docking_programs_child_tabs.docking_programs.tmp_dir_path)

        shutil.rmtree(tmp_path)
        self.save_to_file_choise_window.close()


    def show_plot_distribution_window(self):
        pass

        self.plot_distribution_window = NewWindow(parent = self.main_window,
        title = "Analyses of the results", upper_frame_title = "",
        submit_command = self.set_pharma_btn_func, submit_button_text= "Set",
        with_scroll = True)

        # Create Options to Set Pharmacophoric Restrains
        self.import_label = QtWidgets.QPushButton("Import Atomic Selection")
        self.import_label.clicked.connect(self.import_atomic_from_pymol)
        self.import_label_box = QtWidgets.QComboBox()
        self.import_label_box.currentTextChanged.connect(self.zoom_sele)
        self.import_label_box.view().pressed.connect(self.zoom_sele)

        self.pharma_widgets_list.append(self.import_label)
        self.pharma_widgets_list.append(self.import_label_box)

        self.tolerance_radius_label = QtWidgets.QLabel("Tolerance Radius")
        self.tolerance_radius_box = QtWidgets.QDoubleSpinBox() ## TODO --- QUALI VALORI POSSIBILI??
        self.pharma_widgets_list.append(self.tolerance_radius_label)
        self.pharma_widgets_list.append(self.tolerance_radius_box)

        self.restrain_type_label = QtWidgets.QLabel("Restrain Type")
        self.restrain_type_box = QtWidgets.QComboBox()
        self.pharma_widgets_list.append(self.restrain_type_label)
        self.pharma_widgets_list.append(self.restrain_type_box)
