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
from pymol import cmd
from pymol.cgo import *
from pymol.vfont import plain
from pymol import stored

from pymol.Qt import QtWidgets, QtCore, QtGui

import os
import shutil
import warnings
import math
import subprocess
import statistics
from pathlib import Path

from lib.docking_program_main.Functions.threads import Protocol_exec_dialog
from lib.docking_program_main.Functions.handle_widgets import HandleWidgets

from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.PDBIO import Select


class ADFR_docking():


    """
    A class to represent the ADFR docking process

    """


    def __init__(self, tab, main, ligand, receptor, cavity = None):

        self.tab = tab
        self.main = main
        self.thread = tab

        self.docking_completed = False
        self.interrupt = False
        self.use_flex_protocol = False

        # Initialize Standard Docking parameters
        self.receptor_to_dock = receptor
        self.ligand_to_dock = ligand

        self.ga_evol = str(self.tab.ga_evol.value())
        self.ga_threshold = str(self.tab.ga_threshold.value())
        #self.max_ga_eval = str(self.tab.max_ga_eval.value())
        self.max_gen = str(self.tab.max_gen.value())
        self.buffer = str(self.tab.buffer_box.value())

        if self.tab.use_flex_vina_cb.isChecked():
            self.use_flex_protocol = True

            self.flex_residues = self.tab.flex_arg_edit.text()
            self.check_valid_flex_input()

        if cavity is not None:

            # If cavity is specified in input -- Custom protocol

            name_cav = cavity

            # Check if Reference ligand or grid parameters

            self.grid = re.search("Grid Center_", name_cav)
            self.cavity_name = name_cav

        else:

            # If cavity is not specified in input -- AllvsAll protocol

            name_cav = self.tab.loaded_cavities.currentText()

            # Check if Reference ligand or grid parameters

            self.grid = re.search("Grid Center_", name_cav)
            self.cavity_name = name_cav


        if self.grid:
            self.cavity_to_dock = self.main.ready_grid_centers[name_cav]
            self.reference_cavity = False

        else:
            self.cavity_to_dock = name_cav
            cmd.save(str(self.cavity_to_dock + ".pdb"), self.cavity_to_dock, state = -1, format = 'pdb')
            self.reference_cavity = True

        # Initialize names and paths
        self.results_file_name = str("Run_" + str(self.main.adfr_runs) + "_ADFR")
        self.results_file_name_ext = str(self.results_file_name + ".pdbqt")
        self.log_file_name = str(self.results_file_name + "_log.txt")

        # Change directory --> RxDock tmp dir
        os.chdir(self.main.adfr_tmp_dir)

        self.show_resume_window()


    def check_valid_flex_input(self):
        pass


    def create_protein_segments_string(self):

        self.protein_segments_to_exclude_string = ""

        for items in self.tab.protein_segments_to_exclude:
            item_to_write = str(items.split("_")[1])
            self.protein_segments_to_exclude_string += item_to_write


    def show_resume_window(self):

        self.run_grid_creation_adfr()

        self.run_docking_adfr()


    def run_grid_creation_adfr(self):

        if sys.platform == "win32":
            path_to_agfr = "agfr"

        else:
            path_to_agfr = self.main.path_to_agfr

        os.chdir(self.main.adfr_tmp_dir)

        self.generate_grid_adfr_settings = [path_to_agfr,
        "-r",
        str(self.receptor_to_dock + ".pdbqt"),
        "-o",
        str(self.results_file_name + "_grid")]

        if self.reference_cavity:

            self.generate_grid_adfr_settings.extend([
            "-l",
            str(self.ligand_to_dock  + ".pdbqt"),
            "-P",
            str(self.buffer)])

        else:

            x_val = float(self.cavity_to_dock[3])*float(self.cavity_to_dock[6])
            y_val = float(self.cavity_to_dock[4])*float(self.cavity_to_dock[6])
            z_val = float(self.cavity_to_dock[5])*float(self.cavity_to_dock[6])

            # size_x = (x_val/float(self.cavity_to_dock[6]))
            # size_y = (y_val/float(self.cavity_to_dock[6]))
            # size_z = (z_val/float(self.cavity_to_dock[6]))

            # size_x = (x_val/0.375)
            # size_y = (y_val/0.375)
            # size_z = (z_val/0.375)
            #
            # size_x = (float(self.cavity_to_dock[3])/0.375)
            # size_y = (float(self.cavity_to_dock[4])/0.375)
            # size_z = (float(self.cavity_to_dock[5])/0.375)

            self.generate_grid_adfr_settings.extend(["-b",
            "user",
            self.cavity_to_dock[0],
            self.cavity_to_dock[1],
            self.cavity_to_dock[2],
            str(x_val),
            str(y_val),
            str(z_val)])

        if self.use_flex_protocol:

            # Extend options for docking
            self.generate_grid_adfr_settings.extend(["-f", self.flex_residues])


        # if sys.platform == "linux":
        #     try:
        #         output = subprocess.check_output(self.generate_grid_adfr_settings)
        #     except subprocess.CalledProcessError as e:
        #         output = e.output
        #         print(output)
        #         pass
        #     # output = subprocess.run(self.generate_grid_adfr_settings, check = True)
        #
        # elif sys.platform == "win32":
        #     try:
        #         output = subprocess.check_output(self.generate_grid_adfr_settings, shell = True)
        #     except subprocess.CalledProcessError as e:
        #         output = e.output
        #         print(output)
        #         pass
        #
        #     #subprocess.run(self.generate_grid_adfr_settings, shell = True)
        #
        # else:
        #     try:
        #         output = subprocess.check_output(self.generate_grid_adfr_settings)
        #     except subprocess.CalledProcessError as e:
        #         output = e.output
        #         print(output)
        #         pass

            #subprocess.run(self.generate_grid_adfr_settings)


    def run_docking_adfr(self):

                ### To run the ADFR Docking Process ###

        if sys.platform == "win32":
            path_to_ADFR = "adfr"

        else:
            path_to_ADFR = self.main.path_to_ADFR

        self.run_docking_adfr_settings = [path_to_ADFR,
        "-t", str(self.results_file_name + "_grid.trg"),
        "-l", str(self.ligand_to_dock  + ".pdbqt"),
        "-n", self.ga_evol,
        "-g", self.max_gen,
        "-s", self.ga_threshold,
        "-o", str(self.results_file_name + "_")]

        self.main.adfr_runs += 1


    #     if sys.platform == "linux":
    #         # try:
    #         #     output = subprocess.check_output(self.run_docking_adfr_settings)
    #         # except subprocess.CalledProcessError as e:
    #         #     output = e.output
    #         #     print(output)
    #         #     pass
    #
    #         try:
    #             subprocess.run(self.run_docking_adfr_settings, check = True)
    #         except Exception as e:
    #             print(e)
    #             pass
    #
    #     elif sys.platform == "win32":
    #         # try:
    #         #     output = subprocess.check_output(self.run_docking_adfr_settings, shell = True)
    #         # except subprocess.CalledProcessError as e:
    #         #     output = e.output
    #         #     print(output)
    #         #     pass
    #
    #         subprocess.run(self.run_docking_adfr_settings, shell = True)
    #
    #     else:
    #         # try:
    #         #     output = subprocess.check_output(self.run_docking_adfr_settings)
    #         # except subprocess.CalledProcessError as e:
    #         #     output = e.output
    #         #     print(output)
    #         #     pass
    #
    #         subprocess.run(self.run_docking_adfr_settings)
    #
    #     self.change_results_name()
    #     self.tab.docking_programs_child_tabs.docking_programs.adfr_runs += 1
    #
    #
    # def change_results_name(self):
    #
    #     self.file_path = os.path.join(self.tab.docking_programs_child_tabs.docking_programs.adfr_tmp_dir, self.results_file_name_ext)
    #
    #     for file in os.listdir(self.tab.docking_programs_child_tabs.docking_programs.adfr_tmp_dir):
    #         temp = re.search("_out", file)
    #         temp2 = re.search(".log", file)
    #         if temp:
    #             os.rename(file, self.file_path)
    #         if temp2:
    #             os.rename(file, str(file.replace("grid.log", "")) + "log.txt")


    def check_if_docking_completed(self):

        self.file_path = os.path.join(self.main.adfr_tmp_dir, self.results_file_name_ext)

        self.change_results_name()

        if Path(self.file_path).is_file():

            if os.path.getsize(self.file_path):
                self.docking_completed = True
                self.main.adfr_runs += 1

            else:
                self.docking_completed = False
                QtWidgets.QMessageBox.warning(self.tab, "", str("Something went wrong during Docking. \nPlease check LOG files."))
                os.remove(self.file_path)
                self.main.adfr_runs += 1

        else:
            self.docking_completed = False
            QtWidgets.QMessageBox.warning(self.tab, "", str("Something went wrong during Docking. \nPlease check LOG files."))
            self.main.adfr_runs += 1



class ADFR_parse_results:


    """
    A class to parse the ADFR results file
    """

    def __init__(self, tab, ligand, results_file_name, results_dict, results_data = [[]]):

        self.tab = tab.tab
        self.last_docking = self.tab.last_docking

        self.results_file_name = results_file_name

        self.results_dict = results_dict

        self.results_data = results_data

        self.SCORE = float(00.00)

        self.docked_ligands_list = []

        input = open(str(self.results_file_name + ".pdbqt"), "rt")

        list_of_list = []
        list_of_all_scores = []
        self.poses_list = []

        n = 0

        for line in input:

            if line.startswith("MODEL"):
                #name = str(line.strip())
                #name = self.last_docking.ligand_to_dock.split("_")[1] ## use this when the plugin will support multiligand files
                name = ligand

                self.poses_list.append(name)
                self.NAME = name
                self.POSE = str(len(self.poses_list))

                self.docked_ligands_list.append(str(name + "_all_poses"))
                self.docked_ligands_list.append(self.NAME)

            elif line.startswith("USER: SCORE"):
                value = str(line.strip())
                self.SCORE = value.split()[-1]

            elif line.startswith("ENDMDL"):
                new_list = [self.NAME, self.POSE, self.SCORE]
                list_of_list.append(new_list)
                list_of_all_scores.append(self.SCORE)


        self.results_data = list_of_list
        self.TOT_SCORE = list_of_all_scores

        input.close()
