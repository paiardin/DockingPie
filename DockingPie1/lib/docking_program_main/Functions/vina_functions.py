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


class Vina_docking():


    """
    A class to represent the Vina docking process

    """


    def __init__(self, tab, ligand, receptor, cavity = None):

        self.tab = tab.tab
        self.thread = tab

        self.docking_completed = False
        self.interrupt = False
        self.use_flex_protocol = False

        # Initialize Standard Docking parameters
        self.receptor_to_dock = receptor
        self.poses = str(self.tab.poses_box.value())
        # self.ligand_to_dock = self.tab.ready_ligands_list.currentText()
        self.ligand_to_dock = ligand

        self.exhaustiveness = self.tab.exhaustiveness_box.value()
        self.energy = self.tab.energy_box.value()

        if cavity is not None:

            name_cav = cavity
            self.cavity_to_dock = self.tab.docking_programs_child_tabs.docking_programs.ready_grid_centers[name_cav]
            self.cavity_name = name_cav

        else:

            name_cav = self.tab.loaded_cavities.currentText()
            self.cavity_to_dock = self.tab.docking_programs_child_tabs.docking_programs.ready_grid_centers[name_cav]
            self.cavity_name = name_cav

        # Initialize Flexible Docking
        if self.tab.use_flex_vina_cb.isChecked():
            self.use_flex_protocol = True

            self.flex_residues = self.tab.flex_arg_edit.text()
            self.check_valid_flex_input()

        # Initialize names and paths
        self.results_file_name = str("Run_" + str(self.tab.docking_programs_child_tabs.docking_programs.vina_runs) + "_Vina")
        self.results_file_name_ext = str(self.results_file_name + ".pdbqt")
        self.log_file_name = str(self.results_file_name + "_log.txt")

        # Change directory --> Vina tmp dir
        os.chdir(self.tab.docking_programs_child_tabs.docking_programs.vina_tmp_dir)

        self.show_resume_window()


    def check_valid_flex_input(self):
        pass


    def show_resume_window(self):

        self.run_docking_vina()


    def run_docking_vina(self):

        x_val = float(self.cavity_to_dock[3])*float(self.cavity_to_dock[6])
        y_val = float(self.cavity_to_dock[4])*float(self.cavity_to_dock[6])
        z_val = float(self.cavity_to_dock[5])*float(self.cavity_to_dock[6])

        self.path_to_vina = self.tab.docking_programs_child_tabs.docking_programs.path_to_vina

        self.run_docking_vina_settings = [self.path_to_vina,
        "--receptor", str(self.receptor_to_dock + ".pdbqt"),
        "--ligand", str(self.ligand_to_dock  + ".pdbqt"),
        "--center_x", self.cavity_to_dock[0],
        "--center_y", self.cavity_to_dock[1],
        "--center_z", self.cavity_to_dock[2],
        "--size_x", str(x_val),
        "--size_y", str(y_val),
        "--size_z", str(z_val),
        "--out", self.results_file_name_ext,
        "--exhaustiveness", str(self.exhaustiveness),
        "--num_modes", str(self.poses),
        "--energy_range", str(self.energy),
        "--log", self.log_file_name]

        if self.use_flex_protocol:

            # Generate pdbqt file with flexible side chains
            self.prepare_flex_receptor_path = os.path.join(self.tab.docking_programs_child_tabs.docking_programs.config_path, "prepare_flexreceptor4.py")
            self.preapre_flex_receptors_settings = ["python",
            self.prepare_flex_receptor_path,
            "-r", str(self.receptor_to_dock + ".pdbqt"),
            "-s", self.receptor_to_dock + ":" + self.flex_residues, "-v"]

            try:
                subprocess.run(self.preapre_flex_receptors_settings,
                check = True)

            except subprocess.CalledProcessError as error:
                print(error)
                self.docking_completed = False

            flexible_receptor_name = str(self.receptor_to_dock + "_flex.pdbqt")

            # Extend options for docking
            self.run_docking_vina_settings.extend(["--flex", flexible_receptor_name])

        self.tab.docking_programs_child_tabs.docking_programs.vina_runs += 1


    def check_if_docking_completed(self):

        file_path = os.path.join(self.tab.docking_programs_child_tabs.docking_programs.vina_tmp_dir, self.results_file_name_ext)

        if Path(file_path).is_file():

            if os.path.getsize(file_path):
                self.docking_completed = True
                self.tab.docking_programs_child_tabs.docking_programs.vina_runs += 1

            else:
                self.docking_completed = False

                if len(self.tab.ligands_to_dock) > 1 or len(self.tab.receptors_to_dock) > 1:

                    os.remove(file_path)
                    self.tab.docking_programs_child_tabs.docking_programs.vina_runs += 1

                else:
                    os.remove(file_path)
                    self.tab.docking_programs_child_tabs.docking_programs.vina_runs += 1
                    QtWidgets.QMessageBox.warning(self.tab, "", str("Something went wrong during Docking. \nPlease check LOG files."))

        else:
            self.docking_completed = False
            QtWidgets.QMessageBox.warning(self.tab, "", str("Something went wrong during Docking. \nPlease check LOG files."))
            self.tab.docking_programs_child_tabs.docking_programs.vina_runs += 1


class Vina_Parse_Results:


    """
    A class to parse the Smina results file
    """


    def __init__(self, tab, results_file_name, results_dict, poses, results_data = [[]]):


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

            if line.startswith("REMARK VINA RESULT: "):

                values_line = str(line.strip())

                self.SCORE = values_line.split()[3]
                self.RMSD_lb = values_line.split()[4]
                self.RMSD_ub = values_line.split()[5]

            elif line.startswith("MODEL "):

                # name_line = str(line.strip())
                # name = name_line.split("/")[-1]

                name = self.last_docking.ligand_to_dock
                self.poses_list.append(name)
                self.NAME = name
                self.POSE = str(len(self.poses_list))

            elif line.startswith("ENDMDL"):
                new_list = [self.NAME, self.POSE, self.SCORE, self.RMSD_lb, self.RMSD_ub]
                list_of_list.append(new_list)
                list_of_all_scores.append(self.SCORE)


        self.results_data = list_of_list
        self.TOT_SCORE = list_of_all_scores

        input.close()


    def check_name(self, name):

        if self.poses_list:

            for i in self.poses_list:
                if str(i) == name:
                    return True
                else:
                     return False
        else:
            return True
