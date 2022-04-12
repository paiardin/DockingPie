# Copyright 2022 by Serena Rosignoli. All rights reserved.
# This code is part of the DockingPie package and governed by its license. Please
# see the LICENSE file.


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

# Statistic
import math
import statistics

# RMSD calulating module
try:
    from spyrmsd import io, rmsd
except:
    pass

#numpy
import numpy as np

#
from lib.docking_program_main.Functions.pymol_interactions import HideEverythingPyMOL
from lib.docking_program_main.Functions.general_docking_func import Calculate_RMSD
from lib.docking_program_main.tables.tables import *
from lib.docking_program_main.docking_program_gui.new_windows import NewWindow
from lib.docking_program_main.Functions.threads import Protocol_exec_dialog
# csv module
import csv


#
# class ClusteringConsensusProtocol():
#
#     def __init__(self, tab,
#     data,
#     ):
#
#         self.tab = tab
#
#         ##TODO Implement all types --> still depends on rmsd implementation
#         self.consensus_score_type = self.tab.box_consensus_type.currentText()
#
#         # Save the Checked Docking Runs in tmp dir
#         # In this way, all the states of the ligands are saved in the same file.
#         ## MAYBE: create a copy in consensus dir instead of saving directly from pymol. Ita voids problems if the user cancels or modify something, but at the same time an expert user may want to do so
#         for runs in first_list_consensus:
#             cmd.save(str(os.path.join(self.tab.docking_programs.consensus_tmp_dir, runs) + ".sdf"), runs, format = 'sdf', state = 0)
#
#         for runs in second_list_consensus:
#             cmd.save(str(os.path.join(self.tab.docking_programs.consensus_tmp_dir, runs) + ".sdf"), runs, format = 'sdf', state = 0)
#
#         for runs in third_list_consensus:
#             cmd.save(str(os.path.join(self.tab.docking_programs.consensus_tmp_dir, runs) + ".sdf"), runs, format = 'sdf', state = 0)
#
#         os.chdir(self.tab.docking_programs.consensus_tmp_dir)
#
#         self.cluster_one = []
#         self.cluster_two = []
#         self.cluster_three = []
#         self.cluster_four = []
#
#         # Iterate over the runs chosen in the first list
#         for first_runs in first_list_consensus:
#
#             ref = io.loadallmols(str(first_runs + ".sdf")) # Read all the poses of each run
#
#             # The spyrmsd class needs atomic coordinates, atomic number and the molecular adjacency
#             # matrix to compute the standard RMSD with spyrmsd.rmsd.symmrmsd.
#
#             # iterate over the poses of the first run
#             for refs in ref:
#                 refs.strip()
#
#             for references in ref:
#                 coords_ref = references.coordinates
#                 anum_ref = references.atomicnums
#                 adj_ref = references.adjacency_matrix
#
#                 # Iterate over the runs chosen in the second list
#                 for second_runs in second_list_consensus:
#
#                     # Check if the runs iterated were computed with the same ligand ### TODO: FIND ANOTHER METHOD
#                     first_ligand = self.tab.docking_programs.all_runs[first_runs]["docking_run"].ligand_to_dock.split("_")[1]
#                     second_ligand = self.tab.docking_programs.all_runs[second_runs]["docking_run"].ligand_to_dock.split("_")[1]
#
#                     # if the runs iterated were computed with the same ligand,
#                     # read all the poses of the second run
#                     if first_ligand == second_ligand:
#                         ligand_name = first_ligand
#
#                         mols = io.loadallmols(str(second_runs + ".sdf"))
#                         for mol in mols:
#                             mol.strip()
#                         # Create Adj matrix for rxdock otuput
#                         coords = [mol.coordinates for mol in mols]
#                         anum = mols[0].atomicnums
#                         adj = mols[0].adjacency_matrix
#
#                         # compute the RMSD between the iterated pose of the first run and all the poses of the second run
#                         RMSD = rmsd.symmrmsd(coords_ref, coords, anum_ref, anum, adj_ref, adj)
#
#                         # Check if the computer values of RMSD are less or equal than the threshold chosen by the user
#                         i = 0
#                         for RMSD_Value in RMSD:
#
#                             if RMSD_Value <= self.tab.label_rmsd_box.value():
#                                 rmsd_value = round(RMSD_Value, 3)
#                                 tmp_list = []
#
#                                 first_ligand_data = ConsensusData(self, element = first_runs, rank = ref.index(references)+1, list = first_list_consensus)
#                                 second_ligand_data = ConsensusData(self, element = second_runs, rank = i+1, list = second_list_consensus)
#
#                                 if self.consensus_score_type == "Rank by Rank":
#                                     self.score_name = "RbR"
#                                     consensus_score = ((ref.index(references)+1) + (i+1))/2
#
#                                 if self.consensus_score_type == "Rank by Vote":
#                                     self.score_name = "RbV"
#                                     consensus_score = ((ref.index(references)+1) + (i+1))/2
#
#                                 if self.consensus_score_type == "Average of Auto-Scaled Scores":
#                                     self.score_name = "AASS"
#
#                                     first_score = (first_ligand_data.POSE_SCORE - first_ligand_data.MIN)/(first_ligand_data.MAX - first_ligand_data.MIN)
#                                     second_score = (second_ligand_data.POSE_SCORE - second_ligand_data.MIN)/(second_ligand_data.MAX - second_ligand_data.MIN)
#
#                                     consensus_score = round((first_score + second_score)/2, 3)
#
#                                 if self.consensus_score_type == "Z-scores":
#                                     self.score_name = "Z-score"
#
#                                     first_score = (first_ligand_data.POSE_SCORE - first_ligand_data.MEAN)/(first_ligand_data.STD)
#                                     second_score = (second_ligand_data.POSE_SCORE - second_ligand_data.MEAN)/(second_ligand_data.STD)
#
#                                     consensus_score = round((first_score + second_score)/2, 3)
#
#                                 if self.consensus_score_type == "Exponential":
#                                     print("todo")
#
#                                 tmp_list.extend([first_runs, second_runs, ligand_name, str(rmsd_value), ref.index(references)+1, i+1, first_ligand_data.POSE_SCORE, second_ligand_data.POSE_SCORE, str(consensus_score)])
#                                 i += 1
#                                 list_of_list.append(tmp_list)
#
#                             else:
#                                 i += 1
#
#         if list_of_list:
#             for data in list_of_list:
#
#                 self.consensus_new_table = TableView(parent=self.tab.docking_programs, data=list_of_list,
#                                            row_labels=None, row_labels_height=25,
#                                            column_labels=["Program 1", "Program 2", "NAME", "RMSD", "RANKING 1", "RANKING 2", "SCORE 1", "SCORE 2", str("SCORE (" + self.score_name + ")")],
#                                            sortable=True)
#
#                 self.consensus_new_table.itemDoubleClicked.connect(self.on_click_show_consensus)
#
#
#             self.tab.table_scroll_layout.addWidget(self.consensus_new_table, 2, 0)
#         else:
#             pass
#


class ConsensusScoringFunctions():


    def _init__(self,
    tab,
    ):
        self.tab = tab


    def rank_by_rank(self, list_of_ranks):

        array_mean = np.array(list_of_ranks)

        self.consensus_score_type = "Rank by Rank"
        self.score_name = "RbR"
        self.consensus_score = round(np.mean(array_mean), 3)

        return self.consensus_score


    def rank_by_vote(self, list_of_ranks):

        array_mean = np.array(list_of_ranks)

        self.consensus_score_type = "Rank by Vote"
        self.score_name = "RbV"
        self.consensus_score = round(np.mean(array_mean), 3)

        return self.consensus_score


    def average_autoscaled_scores(self, x_data, y_data):

        self.consensus_score_type = "Average of Auto-Scaled Scores"
        self.score_name = "AASS"

        if x_data.MAX == x_data.MIN:
            first_score = np.NaN
            second_score = np.NaN
        elif y_data.MAX == y_data.MIN:
            first_score = np.NaN
            second_score = np.NaN
        else:
            first_score = (x_data.POSE_SCORE - x_data.MIN)/(x_data.MAX - x_data.MIN)
            second_score = (y_data.POSE_SCORE - y_data.MIN)/(y_data.MAX - y_data.MIN)

        array_mean = np.array([first_score, second_score])
        self.consensus_score = round(np.mean(array_mean), 3)

        return self.consensus_score


    def z_scores(self, x_data, y_data):

        self.consensus_score_type = "Z-scores"
        self.score_name = "Z-score"

        if x_data.MAX == x_data.MIN:
            first_score = np.NaN
            second_score = np.NaN
        elif y_data.MAX == y_data.MIN:
            first_score = np.NaN
            second_score = np.NaN
        else:
            first_score = (x_data.POSE_SCORE - x_data.MEAN)/(x_data.STD)
            second_score = (y_data.POSE_SCORE - y_data.MEAN)/(y_data.STD)

        array_mean = np.array([first_score, second_score])
        self.consensus_score = round(np.mean(array_mean), 3)

        return self.consensus_score


    def exponential(self):

        self.consensus_score_type = "Exponential"
        print("todo")



class ConsensusMatrixCell(ConsensusScoringFunctions):

    def __init__(self,
    tab,
    x_ind,
    y_ind,
    x,
    y,
    compute_rmsd = True
    ):

        self.tab = tab.tab

        self.x_index = x_ind
        self.y_index = y_ind

        self.x_data_name = x
        self.y_data_name = y

        self.x_rank = (self.x_data_name.split("_")[4]).replace(".sdf", "")
        self.y_rank = (self.y_data_name.split("_")[4]).replace(".sdf", "")

        self.x_run = self.x_data_name.replace("_pose_" + str(self.x_rank) + ".sdf", "")
        self.y_run = self.y_data_name.replace("_pose_" + str(self.y_rank) + ".sdf", "")

        self.x_data = ConsensusData(self, element = self.x_run, rank = self.x_rank)
        self.y_data = ConsensusData(self, element = self.y_run, rank = self.y_rank)

        os.chdir(self.tab.docking_programs.consensus_tmp_dir)

        self.calculate_rmsd()

        self.rank_by_rank = self.rank_by_rank(list_of_ranks = [float(self.x_rank), float(self.y_rank)])
        self.rank_by_vote = self.rank_by_vote(list_of_ranks = [float(self.x_rank), float(self.y_rank)])
        self.average_autoscaled_scores = self.average_autoscaled_scores(x_data = self.x_data, y_data = self.y_data)
        self.z_scores = self.z_scores(x_data = self.x_data, y_data = self.y_data)


    def calculate_rmsd(self):

        self.rmsd = Calculate_RMSD(self.tab,
        self.x_data_name,
        self.y_data_name,
        warning = False)

        self.rmsd_list = self.rmsd.rmsd_list

        if self.rmsd.rmsd_computed:

            return float(self.rmsd.rmsd_list[0])


class ConsensusProtocol(ConsensusMatrixCell):

    def __init__(self,
    tab,
    rmsd_protocol,
    smina_poset,
    rxdock_poset,
    vina_poset,
    adfr_poset,
    smina_runs_list = [],
    rxdock_runs_list = [],
    vina_runs_list = [],
    adfr_runs_list = []):


        self.tab = tab


        self.single_program_list = []
        self.complete_list = []


        self.rmsd_threshold_value = self.tab.label_rmsd_box.value()


        self.create_input_list_consensus_matrix(smina_runs_list, smina_poset, single_program = True)
        self.create_input_list_consensus_matrix(rxdock_runs_list, rxdock_poset)
        self.create_input_list_consensus_matrix(vina_runs_list, vina_poset)
        self.create_input_list_consensus_matrix(adfr_runs_list, adfr_poset)

        # self.docking_dialog = Dockings_dialog(self,
        # number_of_dockings_to_do = number_of_dockings_to_do,
        # receptors_to_dock = self.receptors_to_dock,
        # ligands_to_dock = self.ligands_to_dock,
        # dockings_to_do = self.dockings_to_do)
        #
        # self.docking_dialog.setModal(True)
        # self.docking_dialog.exec_()
        #

        if rmsd_protocol == "clustered":
            p_dialog = Protocol_exec_dialog(app=self.tab, docking_pie=self.tab,
                                            function=self.clustered_rmsd_protocol,
                                            args=(),
                                            wait_start=0.4, wait_end=0.4,
                                            lock=True,
                                            stdout_silence=False,
                                            title="Consensus Docking",
                                            label_text="Analysis of the results. \nPlease wait.")
            p_dialog.exec_()

            #self.clustered_rmsd_protocol()
        elif rmsd_protocol == "paired":
            p_dialog = Protocol_exec_dialog(app=self.tab, docking_pie=self.tab,
                                            function=self.paired_rmsd_protocol,
                                            args=(),
                                            wait_start=0.4, wait_end=0.4,
                                            lock=True,
                                            stdout_silence=False,
                                            title="Consensus Docking",
                                            label_text="Analysis of the results. \nPlease wait.")
            p_dialog.exec_()

            self.paired_protocol_update_table()


        elif rmsd_protocol == "normsd":
            p_dialog = Protocol_exec_dialog(app=self.tab, docking_pie=self.tab,
                                            function=self.no_rmsd_protocol,
                                            args=(),
                                            wait_start=0.4, wait_end=0.4,
                                            lock=True,
                                            stdout_silence=False,
                                            title="Consensus Docking",
                                            label_text="Analysis of the results. \nPlease wait.")
            p_dialog.exec_()

            self.no_rmsd_protocol_update_table()


    def create_input_list_consensus_matrix(self, list, pose_threshold, single_program = False):

        index = 1
        prev_run = "dummy"
        for run in list:
            for pose_value, pose in enumerate(range(cmd.count_states(run))):
                if run != prev_run:
                    index = 1
                else:
                    index = index
                name = run + "_pose_" + str(index) + ".sdf"
                consensus_tmp_dir_path = str(os.path.join(self.tab.docking_programs.consensus_tmp_dir, run) + "_pose_" + str(index) + ".sdf")
                cmd.save(consensus_tmp_dir_path, run, format = 'sdf', state = index)
                cmd.load(consensus_tmp_dir_path, name)
                cmd.group("Consensus", members=name, action='auto', quiet=1)

                if pose_threshold == "All":

                    self.complete_list.append(name)
                    if single_program:
                        self.single_program_list.append(name)
                    prev_run = run
                    index = index + 1

                else:

                    if pose_value+1 <= int(pose_threshold):

                        self.complete_list.append(name)
                        if single_program:
                            self.single_program_list.append(name)
                        prev_run = run
                        index = index + 1


    def paired_rmsd_protocol(self):

        # Create an empty array
        self.array = np.empty((len(self.complete_list), len(self.complete_list)))

        # Create an empty dictionary to store information of each cell
        self.consensus_matrix_dict = {}

        # Generate Matrix
        for x_ind, x in enumerate(self.complete_list):

            for y_ind, y in enumerate(self.complete_list):

                cell = ConsensusMatrixCell(self, x_ind, y_ind, x, y)

                if cell.rmsd.rmsd_computed == False:
                    rmsd_value = np.NaN
                else:
                    rmsd_value=float(cell.rmsd_list[0])

                self.consensus_matrix_dict[str(x_ind) + ":" + str(y_ind)] = {}
                self.consensus_matrix_dict[str(x_ind) + ":" + str(y_ind)]["cell"] = cell
                self.consensus_matrix_dict[str(x_ind) + ":" + str(y_ind)]["rmsd"] = rmsd_value

                self.array[x_ind, y_ind] = rmsd_value


    def paired_protocol_update_table(self):

        # Create Table with the generated matrix
        self.consensus_matrix_table = TableView(parent=self.tab.docking_programs, data=self.array,
                                   row_labels_height=25,
                                   column_labels= self.complete_list,
                                   row_labels = self.complete_list,
                                   sortable=True)

        self.consensus_matrix_table.itemDoubleClicked.connect(lambda: self.on_click_matrix_table(self.consensus_matrix_table))

        # Create window to put the generated matrix
        self.consensus_matrix_window = NewWindow(parent = self.tab,
        title = "Consensus Matrix", upper_frame_title = "",
        submit_command = self.update_consensus_table, submit_button_text= "Update Consensus Table",
        with_scroll = True)

        # Create widget for the matrix
        self.consensus_matrix_widget = QtWidgets.QWidget()
        self.consensus_matrix_scroll = QtWidgets.QScrollArea()
        self.consensus_matrix_scroll.setWidgetResizable(True)
        self.consensus_matrix_scroll.setWidget(self.consensus_matrix_widget)

        # Set the layout of the Scroll Area for the Table
        self.consensus_matrix_scroll_layout = QtWidgets.QGridLayout()
        self.consensus_matrix_widget.setLayout(self.consensus_matrix_scroll_layout)

        # Add matrix to window
        self.consensus_matrix_window.middle_layout_type.addWidget(self.consensus_matrix_table, 0, 0, 1, 3)

        score_label = QtWidgets.QLabel()
        self.consensus_matrix_window.middle_layout_type.addWidget(score_label, 1, 1, 1, 1)

        self.box_consensus_type = QtWidgets.QComboBox()
        self.box_consensus_type.addItems(["Rank by Rank", "Average of Auto-Scaled Scores", "Z-scores"])
        self.consensus_matrix_window.middle_layout_type.addWidget(self.box_consensus_type, 1, 2, 1, 1)

        self.rmsd_slider = QtWidgets.QSlider(Qt.Horizontal)
        self.rmsd_slider.setMinimum(1)
        self.rmsd_slider.setMaximum(10)
        self.rmsd_slider.setValue(1)
        self.rmsd_slider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        self.rmsd_slider.setTickInterval(0.5)
        self.rmsd_slider.valueChanged.connect(self.rmsd_slider_changed)

        new_rmsd_value = "1"

        self.rmsd_label = QtWidgets.QLabel("RMSD threshold: " + new_rmsd_value)

        self.consensus_matrix_window.middle_layout_type.addWidget(self.rmsd_slider, 2, 2, 1, 1)
        self.consensus_matrix_window.middle_layout_type.addWidget(self.rmsd_label, 2, 1)
        self.consensus_matrix_window.show()

        self.subset_matrix(self.rmsd_threshold_value, update_from_matrix = False)


    def subset_matrix(self, rmsd_threshold_value, update_from_matrix = False):

        list_of_new_entries = []
        list_of_list = []

        # # Subset array based on RMSD threshold value
        result = np.where(self.array <= rmsd_threshold_value)
        listOfCoordinates = list(zip(result[0], result[1]))

        # Exctract info from indexes
        list_of_list = []

        for cord in listOfCoordinates:
            tmp_list = []

            # Get index of cell
            cell_index = str(cord[0]) + ":" + str(cord[1])

            # Initialize cell
            cell = self.consensus_matrix_dict[cell_index]["cell"]

            # Get rmsd value from that cell
            rmsd = cell.rmsd.rmsd_list[0]

            # Get score value
            if update_from_matrix:
                score_box  = self.box_consensus_type.currentText()
            else:
                score_box = self.tab.box_consensus_type.currentText()

            score = self.get_consensus_score_type(cell, score_box)

            new_entry = (cell.x_data_name, cell.y_data_name)
            list_of_new_entries.append(new_entry)
            self.check_if_already_analyzed(new_entry, list_of_new_entries)

            if rmsd > 0.0 and self.already_analyzed == False:

                # Extend tmp_list to build table with info of that cell
                tmp_list.extend([cell.x_data_name,
                cell.y_data_name,
                rmsd,
                cell.x_rank,
                cell.y_rank,
                cell.x_data.POSE_SCORE,
                cell.y_data.POSE_SCORE,
                self.score_value])

                list_of_list.append(tmp_list)

        if list_of_list:

            self.consensus_new_table = TableView(parent=self.tab.docking_programs, data=list_of_list,
                                       row_labels=None, row_labels_height=25,
                                       column_labels=["Program 1", "Program 2", "RMSD", "RANKING 1", "RANKING 2", "SCORE 1", "SCORE 2", "CONSENSUS SCORE\n" + self.score_name],
                                       sortable=True)

            self.consensus_new_table.itemDoubleClicked.connect(lambda: self.on_click_consensus_table(self.consensus_new_table))

            self.tab.table_scroll_layout.addWidget(self.consensus_new_table, 2, 0)

            self.tab.show_consensus_matrix_btn.setEnabled(True)
            self.tab.to_csv_btn.setEnabled(True)
            self.tab.save_cons_btn.setEnabled(True)


    def check_if_already_analyzed(self, new_entry, list_of_new_entries):

        self.already_analyzed = False

        for x, y in list_of_new_entries:
            if y == new_entry[0] and x == new_entry[1]:
                self.already_analyzed = True


    def on_click_matrix_table(self, table):

        row_index = table.currentIndex().row()
        column_index_1 = table.currentIndex().column()

        name_1  = table.horizontalHeaderItem(row_index).text()
        name_2 = table.verticalHeaderItem(column_index_1).text()

        HideEverythingPyMOL(self, to_show = [name_1, name_2])
        cmd.orient(name_1)



    def on_click_consensus_table(self, table, cluster = False):

        if cluster:
            column_index_1 = (TableView.get_column_index_from_header(self, table = table, header = "Cluster\nReference"))[0]
            column_index_2 = (TableView.get_column_index_from_header(self, table = table, header = "NAME"))[0]

        else:
            column_index_1 = (TableView.get_column_index_from_header(self, table = table, header = "Program 1"))[0]
            column_index_2 = (TableView.get_column_index_from_header(self, table = table, header = "Program 2"))[0]


        row_index = table.currentIndex().row()

        model = table.model()

        index_1 = model.index(row_index, column_index_1)
        name_1 = model.data(index_1)

        index_2 = model.index(row_index, column_index_2)
        name_2 = model.data(index_2)

        HideEverythingPyMOL(self, to_show = [name_1, name_2])
        cmd.orient(name_1)


    def rmsd_slider_changed(self):
        new_rmsd_value = str(self.rmsd_slider.value())
        self.rmsd_label.setText("RMSD threshold: " + new_rmsd_value)


    def get_consensus_score_type(self, cell, score_box):

        # Get consensus score type
        self.consensus_score_type = score_box

        if self.consensus_score_type == "Rank by Rank":
            self.score_name = "RbR"
            self.score_value = cell.rank_by_rank

        if self.consensus_score_type == "Rank by Vote":
            self.score_name = "RbV"
            self.score_value = cell.rank_by_vote

        if self.consensus_score_type == "Average of Auto-Scaled Scores":
            self.score_name = "AASS"
            self.score_value = cell.average_autoscaled_scores

        if self.consensus_score_type == "Z-scores":
            self.score_name = "Z-score"
            self.score_value = cell.z_scores

        if self.consensus_score_type == "Exponential":
            print("todo")
            self.score_name = "Z-score"
            self.score_value = cell.z_scores

        return self.score_value


    def update_consensus_table(self):
        self.subset_matrix(self.rmsd_slider.value(), update_from_matrix = True)


    def no_rmsd_protocol(self):

        # Create an empty array
        self.array = np.empty((len(self.complete_list), len(self.complete_list)))

        # Create an empty dictionary to store information of each cell
        self.consensus_matrix_dict = {}

        # Generate Matrix
        for x_ind, x in enumerate(self.complete_list):

            for y_ind, y in enumerate(self.complete_list):

                cell = ConsensusMatrixCell(self, x_ind, y_ind, x, y, compute_rmsd = False)

                self.consensus_matrix_dict[str(x_ind) + ":" + str(y_ind)] = {}
                self.consensus_matrix_dict[str(x_ind) + ":" + str(y_ind)]["cell"] = cell

        score_box = self.tab.box_consensus_type.currentText()

        list_of_new_entries = []
        self.list_of_list = []

        for row_idx, row in enumerate(self.array):

            for col_idx, col in enumerate(self.array):

                tmp_list = []

                # Get index of cell
                cell_index = str(row_idx) + ":" + str(col_idx)

                # Initialize cell
                cell = self.consensus_matrix_dict[cell_index]["cell"]

                score = self.get_consensus_score_type(cell, score_box)

                new_entry = (cell.x_data_name, cell.y_data_name)
                list_of_new_entries.append(new_entry)

                self.check_if_already_analyzed(new_entry, list_of_new_entries)

                same_ligand = cell.rmsd.rmsd_computed

                if self.already_analyzed == False and same_ligand == True:

                    # Extend tmp_list to build table with info of that cell
                    tmp_list.extend([cell.x_data_name,
                    cell.y_data_name,
                    cell.x_rank,
                    cell.y_rank,
                    cell.x_data.POSE_SCORE,
                    cell.y_data.POSE_SCORE,
                    self.score_value])

                    self.list_of_list.append(tmp_list)


    def no_rmsd_protocol_update_table(self):

        self.consensus_new_table = TableView(parent=self.tab.docking_programs, data=self.list_of_list,
                                   row_labels=None, row_labels_height=25,
                                   column_labels=["Program 1", "Program 2", "NAME", "RANKING 1", "RANKING 2", "SCORE 1", "SCORE 2", "CONSENSUS SCORE\n" + self.score_name],
                                   sortable=True)

        self.consensus_new_table.itemDoubleClicked.connect(lambda: self.on_click_consensus_table(self.consensus_new_table))

        self.tab.table_scroll_layout.addWidget(self.consensus_new_table, 2, 0)

        if self.tab.show_consensus_matrix_btn.isEnabled():
            self.tab.show_consensus_matrix_btn.setEnabled(False)

        self.tab.to_csv_btn.setEnabled(True)
        self.tab.save_cons_btn.setEnabled(True)


    def clustered_rmsd_protocol(self):

        list_of_list = []
        tmp_list = []

        list_to_be_analized = self.complete_list

        threshold = 0
        cluster_dict = {}
        index = 0

        while list_to_be_analized:

            threshold += 1
            cluster_dict[threshold] = "cluster" + str(len(cluster_dict))

            if index+1 <= len(list_to_be_analized):
                x = list_to_be_analized[index]

            for y_ind, y in enumerate(list_to_be_analized):

                tmp_list = []
                cell = ConsensusMatrixCell(self, index, y_ind, x, y)

                if float(cell.rmsd_list[0]) < threshold:
                    tmp_list.extend([cluster_dict[threshold], cell.rmsd_list[0], cell.x_data_name, cell.y_data_name, cell.x_rank, cell.y_rank, cell.x_data.POSE_SCORE, cell.y_data.POSE_SCORE])
                    del list_to_be_analized[y_ind]

                    list_of_list.append(tmp_list)

            index += 1

        self.consensus_new_table = TableView(parent=self.tab.docking_programs, data=list_of_list,
                                   row_labels=None, row_labels_height=25,
                                   column_labels=["Cluster", "RMSD", "Cluster\nReference", "NAME", "RANKING 1", "RANKING 2", "SCORE 1", "SCORE 2"],
                                   sortable=True)

        self.tab.table_scroll_layout.addWidget(self.consensus_new_table, 2, 0)

        self.consensus_new_table.itemDoubleClicked.connect(lambda: self.on_click_consensus_table(self.consensus_new_table, cluster = True))

        if self.tab.show_consensus_matrix_btn.isEnabled():
            self.tab.show_consensus_matrix_btn.setEnabled(False)

        self.tab.to_csv_btn.setEnabled(True)
        self.tab.save_cons_btn.setEnabled(True)


        #
        # # Create an empty array
        # self.array = np.empty((len(self.single_program_list), len(self.complete_list)))
        #
        # # Create an empty dictionary to store information of each cell
        # self.consensus_matrix_dict = {}
        #
        # # Generate Matrix
        # for x_ind, x in enumerate(self.single_program_list):
        #
        #     for y_ind, y in enumerate(self.complete_list):
        #
        #         cell = ConsensusMatrixCell(self, x_ind, y_ind, x, y)
        #
        #         self.consensus_matrix_dict[str(x_ind) + ":" + str(y_ind)] = {}
        #         self.consensus_matrix_dict[str(x_ind) + ":" + str(y_ind)]["cell"] = cell
        #         self.consensus_matrix_dict[str(x_ind) + ":" + str(y_ind)]["rmsd"] = float(cell.rmsd_list[0])
        #
        #         self.array[x_ind, y_ind] = float(cell.rmsd_list[0])
        #
        # # Create Table with the generated matrix
        # self.consensus_matrix_table = TableView(parent=self.tab.docking_programs, data=self.array,
        #                            row_labels_height=25,
        #                            column_labels= self.complete_list,
        #                            row_labels = self.single_program_list,
        #                            sortable=True)
        #
        # # Create window to put the generated matrix
        # self.consensus_matrix_window = NewWindow(parent = self.tab,
        # title = "Consensus Matrix", upper_frame_title = "",
        # submit_command = self.update_consensus_table, submit_button_text= "Update Consensus Table",
        # with_scroll = True)
        #
        # # Create widget for the matrix
        # self.consensus_matrix_widget = QtWidgets.QWidget()
        # self.consensus_matrix_scroll = QtWidgets.QScrollArea()
        # self.consensus_matrix_scroll.setWidgetResizable(True)
        # self.consensus_matrix_scroll.setWidget(self.consensus_matrix_widget)
        #
        # # Set the layout of the Scroll Area for the Table
        # self.consensus_matrix_scroll_layout = QtWidgets.QGridLayout()
        # self.consensus_matrix_widget.setLayout(self.consensus_matrix_scroll_layout)
        #
        # # Add matrix to window
        # self.consensus_matrix_window.middle_layout_type.addWidget(self.consensus_matrix_table)
        #
        # self.box_consensus_type = QtWidgets.QComboBox()
        # self.box_consensus_type.addItems(["Rank by Rank", "Rank by Vote", "Average of Auto-Scaled Scores", "Z-scores"])
        # self.consensus_matrix_window.middle_layout_type.addWidget(self.box_consensus_type, 0, 1)
        #
        # self.rmsd_slider = QtWidgets.QSlider(Qt.Horizontal)
        # self.rmsd_slider.setMinimum(0)
        # self.rmsd_slider.setMaximum(10)
        # self.rmsd_slider.setValue(1)
        # self.rmsd_slider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        # self.rmsd_slider.setTickInterval(0.5)
        # self.rmsd_slider.valueChanged.connect(self.rmsd_slider_changed)
        #
        # self.consensus_matrix_window.middle_layout_type.addWidget(self.rmsd_slider, 1, 1)
        #
        # self.consensus_matrix_window.show()
        #
        # score_box = self.tab.box_consensus_type.currentText()
        #
        # list_of_list = []
        #
        # for row_idx, row in enumerate(self.array):
        #     single_row = np.array(row)
        #     result = np.where(single_row < self.rmsd_threshold_value)
        #     listOfCoordinates = list(result[0])
        #
        #     tmp_list = []
        #     mean_score = []
        #     mean_rmsd = []
        #
        #     if listOfCoordinates:
        #
        #         for cord in listOfCoordinates:
        #             cell_index = str(row_idx) + ":" + str(cord)
        #             cell = self.consensus_matrix_dict[cell_index]["cell"]
        #
        #             score = self.get_consensus_score_type(cell, score_box)
        #             mean_score.append(score)
        #
        #             rmsd = cell.rmsd.rmsd_list[0]
        #             mean_rmsd.append(rmsd)
        #
        #             tmp_list.extend([cell.x_data_name,
        #             cell.y_data_name,
        #             "nome_ligando"])
        #
        #         mean_score_array = np.array(mean_score).astype(np.float)
        #         mean_rmsd_array = np.array(mean_rmsd).astype(np.float)
        #
        #         tmp_list.extend([np.mean(mean_rmsd_array),
        #         np.mean(mean_score_array)])
        #
        #         list_of_list.append(tmp_list)
        #
        # list_of_new_entries = []
        # cluster_tree = {}
        # for row_idx, row in enumerate(self.array):
        #     for col_idx, col in enumerate(self.array):
        #
        #             cell_index = str(row_idx) + ":" + str(cord)
        #             cell = self.consensus_matrix_dict[cell_index]["cell"]
        #
        #             new_entry = (cell.x_data_name, cell.y_data_name)
        #             list_of_new_entries.append(new_entry)
        #
        #             self.check_if_already_analyzed(new_entry, list_of_new_entries)
        #
        #             score = self.get_consensus_score_type(cell, score_box)
        #
        #             rmsd = cell.rmsd.rmsd_list[0]
        #
        #             if rmsd < 1:
        #                 cluster_tree["cluster_1"] = cell
        #             elif rmsd > 1 and rmsd < 2:
        #                 cluster_tree["cluster_2"] = cell
        #             elif rmsd < 4:
        #                 cluster_tree["cluster_3"] = cell
        #
        #             tmp_list.extend([cell.x_data_name,
        #             cell.y_data_name,
        #             "nome_ligando"])
        #
        #         mean_score_array = np.array(mean_score).astype(np.float)
        #         mean_rmsd_array = np.array(mean_rmsd).astype(np.float)
        #
        #         tmp_list.extend([np.mean(mean_rmsd_array),
        #         np.mean(mean_score_array)])
        #
        #         list_of_list.append(tmp_list)



    def on_click_show_consensus(self):

        index = self.consensus_new_table.currentIndex()

        self.list_of_pymol_representations = ["lines", "spheres", "mesh", "ribbon", "cartoon", "sticks", "dots", "surface", "labels", "extent", "nonbonded", "nb_spheres", "slice", "extent", "slice", "dashes", "angles", "dihedrals", "cgo", "cell", "callback"]

        for representation in self.list_of_pymol_representations:
            cmd.hide(representation, "all")



class RMSD_Object():

    def __init__(self, protocol,
    file, format, all_poses = False
    ):

        self.protocol = protocol

        if all_poses:

            self.read_file = io.loadallmols(str(file + str("." + format)))

            for poses in self.read_file:
                poses.strip()

            for poses in self.read_file:

                # To consider all the poses of a single file together
                self.coords_list = [poses.coordinates for poses in self.read_file]
                self.anum_list = self.read_file[0].atomicnums
                self.adj_list = self.read_file[0].adjacency_matrix

        else:

            self.read_file = io.loadmol(str(file + str("." + format)))

            self.read_file.strip()

            # The spyrmsd class needs atomic coordinates, atomic number and the molecular adjacency
            # matrix to compute the standard RMSD with spyrmsd.rmsd.symmrmsd.

            # To treat each pose of a single file separatedly
            self.coords = self.read_file.coordinates
            self.anum = self.read_file.atomicnums
            self.adj = self.read_file.adjacency_matrix



class ConsensusData():

    def __init__(self, protocol,
    element, rank, list = ["dummy"]
    ):

        self.protocol = protocol

        if len(list) == 1:

            ### Get Docking Run info ###

            docking_run_results_file = self.protocol.tab.docking_programs.all_runs[element]["docking_run_results_file"]

                # all scores
            list_tot_scores = []
            list_tot_scores.append(docking_run_results_file.TOT_SCORE)
            array = np.array(list_tot_scores).astype(np.float)
            self.DATA = np.array(array)

                # score
            score = docking_run_results_file.TOT_SCORE[int(rank)-1]
            self.POSE_SCORE = np.float64(score)

                # mean, std, min and max of data

            self.MEAN = np.mean(self.DATA)
            self.STD = np.std(self.DATA)
            self.MIN = np.amin(self.DATA)
            self.MAX = np.amax(self.DATA)

                # number of poses generated

            self.NUM_POSES = len(docking_run_results_file.results_data)

        else:
            list_tot_scores = []

            for runs in list:
                list_tot_scores.extend(self.protocol.tab.docking_programs.all_runs[runs]["docking_run_results_file"].TOT_SCORE)

            self.POSE_SCORE = self.protocol.tab.docking_programs.all_runs[element]["docking_run_results_file"].TOT_SCORE[int(rank)-1]

            self.DATA = np.array(list_tot_scores)

            self.MEAN = np.mean(self.DATA)
            self.STD = np.std(self.DATA)
            self.MIN = np.amin(self.DATA)
            self.MAX = np.amax(self.DATA)




def cluster_from_tuple(mytuple):
    """
    Creates a cluster from a tupple formed by a first element being the prototype and
    a second element being the list of all elements.
    """
    prototype = mytuple[0]
    all_elements = list(mytuple[1])
    all_elements.extend([prototype])
    return Cluster(prototype, all_elements)


def get_cluster_sizes(clusters):
    """
    Calculates all the sizes of a clusters list and returns it in a tuple in which
    the first element is the total number of elements, and the second a list containing
    the size of each cluster (maintaining the ordering).
    """
    total_elements = 0
    cluster_sizes = []
    for c in clusters:
        size = c.get_size()
        total_elements = total_elements + size
        cluster_sizes.append(size)
    return total_elements,cluster_sizes


def gen_clusters_from_class_list(group_list,skip_list=[]):
    """
    Generates the clusters that describe a group list. A group list for a N elements clustering
    is defined as the list of N elements with the number of cluster to which each cluster
    belongs. Example: for 4 elements [1,2,3,4] a possible group list would be: [2,1,2,1] which
    means that element 0 and 2 belong to cluster 2 and the others to cluster 2. As it's not possible
    to define a centroid or medioid. ATENTION: the first element of the cluster will be defined as the
    centroid/medoid.
    """
    dic_clusters = {}
    for i in range(len(group_list)):
        if not group_list[i] in skip_list:
            if group_list[i] in dic_clusters.keys():
                dic_clusters[group_list[i]].append(i)
            else:
                dic_clusters[group_list[i]] = [i]
    clusters = []
    for k in dic_clusters.keys():
        clusters.append(Cluster(dic_clusters[k][0],dic_clusters[k]))
    return clusters

class Cluster(object):
    """
    A cluster object is defined a group of elements which have one or more characteristics in common
    and one element which is the most representative element of the cluster.
    """
    most_representative_element = None
    all_elements = []

    def __init__(self, prototype , elements):
        """
        Constructor, needs the prototype and the elements of the clustering.
        TODO: change it by (elements, [prototype]). Prototype must be calculated on demand
        and use bookkeeping
        """
        self.set_elements(elements)
        self.id = ""
        try:
            self.set_prototype(prototype)
        except TypeError:
            raise

    def set_prototype(self,this_one):
        """
        Adds a representative element which must already be inside the
        internal elements list.
        """
        if this_one == None:
            self.prototype = None
        else:
            if this_one in self.all_elements:
                self.prototype = this_one
            else:
                raise TypeError("[Error in Cluster::set_prototype] the prototype is not in the elements list.")

    def set_elements(self,elements):
        self.all_elements = elements

    def get_size(self):
        """
        Returns the size of the cluster (which is indeed the size of its elements list)
        """
        return len(self.all_elements)

    def __eq__(self, other):
        """
        Checks whether two clusters are equal or not. Returns True or False depending
        on it :P
        """
        if(self.get_size() != other.get_size()):
            return False
        else:
            elements = sorted(self.all_elements)
            other_elements = sorted(other.all_elements)
            for i in range(len(elements)):
                if elements[i] != other_elements[i]:
                    return False
            return True

    def __str__(self):
        return "["+str(self.prototype)+str(self.all_elements)+"]"

    def __getitem__(self, index):
        return self.all_elements[index]

    def calculate_biased_medoid(self, condensed_distance_matrix, elements_into_account):
        """
        Calculates the medoid (element with minimal distance to all other objects) of the
        elements of the cluster which are in elements_into_account.
        Note that, even if it is not intuitive, the medoid can be different from the most
        dense point of a cluster.
        """
        all_elems_set = set(self.all_elements)
        accountable_set = set(elements_into_account)

        # Check that elements_into_account is a subset of all_elements
        elem_inters = all_elems_set.intersection(accountable_set)
        if len(elem_inters) != len(elements_into_account):
            print("[ERROR Cluster::calculate_biased_medoid] 'elements_into_account' is not a subset of the elements of this cluster.")
            exit()

        if len(elements_into_account) == 0:
            print("[ERROR Cluster::calculate_biased_medoid] This cluster is empty.")
            return -1

        #average distance of medoid is maximal
        min_dist_pair = (sys.maxint, -1)
        for ei in elements_into_account:
            # Calculate distances for this vs all the others
            # Note that for comparing, the mean is not required,as
            # all have the same amount of elements
            summed_distance = 0
            for ej in elements_into_account:
                summed_distance = summed_distance + condensed_distance_matrix[ei,ej]
            min_dist_pair = min(min_dist_pair,(summed_distance,ei))

        medoid_element = min_dist_pair[1]

        return medoid_element

    def calculate_medoid(self, condensed_distance_matrix):
        """
        Calculates the medoid for all_elements of the cluster and updates the prototype.
        """
        return self.calculate_biased_medoid(condensed_distance_matrix, self.all_elements)

    def get_random_sample(self, n, rand_seed = None):
        """
        Returns a random sample of the elements.
        @param n: Number of random elements to get.
        @param rand_seed: Seed for the random package. Used for testing (repeteability)
        @return: A random sample of the cluster elements.
        """
        if not rand_seed is None:
            random.seed(rand_seed)
        temporary_list = list(self.all_elements)
        random.shuffle(temporary_list)
        return temporary_list[0:n]

    def to_dic(self):
        """
        Converts this cluster into a dictionary (to be used with json serializers).
        """
        json_dic = {}
        elements = sorted([int(self.all_elements[i]) for i in range(len(self.all_elements))])
        elements.append(-1) #capping
        str_elements = ""
        start = elements[0]
        for i in range(1,len(elements)):
            if elements[i-1] != elements[i]-1 :
                if elements[i-1] == start:
                    str_elements += str(start)+", "
                    start = elements[i]
                else:
                    str_elements += str(start)+":"+str(elements[i-1])+", "
                    start = elements[i]
        json_dic["elements"] =str_elements[:-2]

        if self.prototype is not None:
            json_dic["prototype"] = int(self.prototype)

        if self.id != "":
            json_dic["id"] = self.id

        return json_dic

    @classmethod
    def from_dic(cls, cluster_dic):
        """
        Creates a cluster from a cluster dictionary describing it (as reverse operation of
        'to_dic').
        @param cluster_dic: The cluster in dictionary form (output of 'to_dic')
        """
        if "prototype" in cluster_dic:
            proto = cluster_dic["prototype"]
        else:
            proto = None

        if "id" in cluster_dic:
            cid = cluster_dic["id"]
        else:
            cid = None

        values_string_parts = cluster_dic["elements"].split(",");
        elements = []
        for value_part in values_string_parts:
            if ":" in value_part:
                [ini,end] = value_part.split(":")
                elements.extend(range(int(ini),int(end)+1))
            else:
                elements.append(int(value_part))

        cluster = Cluster(proto, elements)
        cluster.id = cid

        return cluster
