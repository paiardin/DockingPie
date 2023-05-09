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


class Smina_docking():


    """
    A class to represent the Smina docking process

    """


    def __init__(self, tab, main,
                 ligand,
                 receptor,
                 exhaustiveness,
                 buffer,
                 energy,
                 rmsd,
                 tmp_dir,
                 cavity_list,
                 scoring = "",
                 use_flex_smina = False,
                 flex_residues = "",
                 poses = 1,
                 cavity = None):

        self.tab = tab
        self.main = main

        self.docking_completed = False
        self.interrupt = False
        self.use_flex_protocol = False
        self.grid = None

        # Initialize Standard Docking parameters
        self.receptor_to_dock = receptor
        self.poses = poses
        self.ligand_to_dock = ligand

        self.exhaustiveness = exhaustiveness
        self.buffer = buffer
        self.energy = energy
        self.rmsd_filter = rmsd
        self.scoring_function = scoring

        if use_flex_smina:
            self.use_flex_protocol = True

            self.flex_residues = ""

            flex_residues = flex_residues
            for res in flex_residues.split(","):
                residue = res.replace(" ", "")
                self.flex_residues += str(self.receptor_to_dock + ":" + residue + ",")

            self.check_valid_flex_input()

        self.grid = re.search("Grid Center_", cavity)
        self.cavity_name = cavity
        self.cavity_to_dock = cavity_list # [x_pos, y_pos, z_pos, x, y, z, spacing] e.g. ['191.48', '171.91', '15.94', '16', '1', '16', '16.0']

        # Check if Reference ligand or grid parameters
        if self.grid or cavity == "consensus_grid":
            self.cavity_to_dock = cavity_list
            self.reference_cavity = False

        else:
            self.cavity_to_dock = cavity
            cmd.save(str(self.cavity_to_dock + ".pdb"), self.cavity_to_dock, state = -1, format = 'pdb')
            self.reference_cavity = True

        # Initialize names and paths
        self.results_file_name = str("Run_" + str(self.main.smina_runs) + "_Smina")
        self.results_file_name_ext = str(self.results_file_name + ".pdb")
        self.log_file_name = str(self.results_file_name + "_log.txt")

        # Change directory --> RxDock tmp dir
        os.chdir(tmp_dir)

        self.show_resume_window()

    def check_valid_flex_input(self):
        pass

    # def create_protein_segments_string(self):
    #
    #     self.protein_segments_to_exclude_string = ""
    #
    #     for items in self.tab.protein_segments_to_exclude:
    #         item_to_write = str(items.split("_")[1])
    #         self.protein_segments_to_exclude_string += item_to_write


    def show_resume_window(self):

        self.run_docking_smina()


    def run_docking_smina(self):

                ### To run the Smina Docking Process ###

        path_to_smina = self.main.path_to_smina

        self.run_docking_smina_settings = [path_to_smina,
        "-r", str(self.receptor_to_dock + ".pdbqt"),
        "-l", str(self.ligand_to_dock  + ".pdbqt"),
        "-o", self.results_file_name_ext,
        "--exhaustiveness", str(self.exhaustiveness),
        "--num_modes", str(self.poses),
        "--energy_range", str(self.energy),
        "--min_rmsd_filter", str(self.rmsd_filter),
        "--log", self.log_file_name]

        if self.scoring_function == "Vinardo":
            self.run_docking_smina_settings.extend(["--scoring", "vinardo"])

        if self.reference_cavity:

            self.run_docking_smina_settings.extend(["--autobox_ligand",
            str(self.cavity_to_dock + ".pdb"),
            "--autobox_add",
            str(self.buffer)])

        else:

            x_val = float(self.cavity_to_dock[3])
            y_val = float(self.cavity_to_dock[4])
            z_val = float(self.cavity_to_dock[5])
            spacing = float(self.cavity_to_dock[6])
# *float(self.cavity_to_dock[6])
# *float(self.cavity_to_dock[6])
# *float(self.cavity_to_dock[6])

            self.run_docking_smina_settings.extend([
            "--center_x", self.cavity_to_dock[0],
            "--center_y", self.cavity_to_dock[1],
            "--center_z", self.cavity_to_dock[2],
            "--size_x", str(x_val),
            "--size_y", str(y_val),
            "--size_z", str(z_val)])


        if self.use_flex_protocol:
            # Generate pdbqt file with flexible side chains
            self.prepare_flex_receptor_path = os.path.join(self.main.config_path, "prepare_flexreceptor4.py")
            self.preapre_flex_receptors_settings = ["python",
            self.prepare_flex_receptor_path,
            "-r", str(self.receptor_to_dock + ".pdbqt"),
            "-s", self.flex_residues, "-v"]

            try:
                subprocess.run(self.preapre_flex_receptors_settings,
                check = True)

            except subprocess.CalledProcessError as error:
                print(error)
                self.docking_completed = False

            flexible_receptor_name = str(self.receptor_to_dock + "_flex.pdbqt")

            # Extend options for docking
            self.run_docking_smina_settings.extend(["--flex", flexible_receptor_name])

        self.main.smina_runs += 1


    def check_if_docking_completed(self):

        file_path = os.path.join(self.main.smina_tmp_dir, self.results_file_name_ext)

        if Path(file_path).is_file():

            if os.path.getsize(file_path):
                self.docking_completed = True
                self.main.smina_runs += 1

            else:
                self.docking_completed = False
                QtWidgets.QMessageBox.warning(self.tab, "", str("Something went wrong during Docking. \nPlease check LOG files."))
                os.remove(file_path)
                self.main.smina_runs += 1

        else:
            self.docking_completed = False
            QtWidgets.QMessageBox.warning(self.tab, "", str("Something went wrong during Docking. \nPlease check LOG files."))
            self.main.smina_runs += 1


class Smina_parse_results:


    """
    A class to parse the Smina results file
    """


    def __init__(self, tab, main, results_file_name, poses, ligand, results_dict = {}, results_data = [[]]):

        self.tab = tab
        self.main = main

        self.ligand_name = ligand

        self.results_file_name = results_file_name

        self.results_dict = results_dict

        self.results_data = results_data

        self.poses = poses

        self.SCORE = float(00.00)

        self.docked_ligands_list = []

        input = open(str(self.results_file_name + ".pdb"), "rt")

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

            elif line.startswith("REMARK minimizedAffinity"):
                value = str(line.strip())
                self.SCORE = float(value[24:])

            elif line.startswith("ENDMDL"):
                new_list = [self.NAME, self.POSE, self.SCORE]
                list_of_list.append(new_list)
                list_of_all_scores.append(self.SCORE)


        self.results_data = list_of_list
        self.TOT_SCORE = list_of_all_scores

        input.close()

#### CODE NEXT: USE WHEN SMINA WILL BE USED FOR THE DOCKING OF MULTIPLE LIGANDS

        # self.results_file_name = results_file_name
        #
        # self.results_dict = results_dict
        #
        # self.results_data = results_data
        #
        # self.poses = poses
        #
        # self.SCORE = float(00.00)
        #
        # self.docked_ligands_list = []
        #
        # input = open(str(self.results_file_name + ".pdb"), "rt")
        # print(input)
        #
        # list_of_list = []
        # list_of_all_scores = []
        # self.poses_list = []
        #
        # n = 0
        #
        # for line in input:
        #
        #     if line.startswith("MODEL"):
        #         name = str(line.strip())
        #         print(name)
        #
        #         if self.check_name(name): # if True, the ligand is the same, the state is different
        #             print("here")
        #             self.poses_list.append(name)
        #             self.NAME = str(name + "_pose_" + str(len(self.poses_list)))
        #
        #         else: # if False, the ligand is changed, thus the tmp list (poses_list) must be empty
        #             print("here2")
        #
        #             self.poses_list.append(name)
        #             self.NAME = str(name + "_pose_" + str(len(self.poses_list)))
        #
        #         self.docked_ligands_list.append(str(name + "_all_poses"))
        #         self.docked_ligands_list.append(self.NAME)
        #
        #     elif line.startswith("REMARK minimizedAffinity"):
        #         value = str(line.strip())
        #         print(value)
        #         print(value[24:])
        #         self.SCORE = float(value[24:])
        #         # self.MINIMIZED_AFFINITY = float(value[25:])
        #
        #     elif line.startswith("ENDMDL"):
        #         new_list = [self.NAME, self.SCORE]
        #         list_of_list.append(new_list)
        #         list_of_all_scores.append(self.SCORE)
        #
        #
        # self.results_data = list_of_list
        # self.TOT_SCORE = list_of_all_scores
        #
        # x = np.array([list_of_list])
        # print(x)
        #
        # y=np.array([np.array(xi) for xi in list_of_list])
        # print(y)
        #
        # print(self.docked_ligands_list)
        # # y = x.astype(np.str)
        #
        # np.savetxt("rxdock_results.csv", y, delimiter=",", fmt='%s')
        # #

    def check_name(self, name):

        if self.poses_list:

            for i in self.poses_list:
                if str(i) == name:
                    return True
                else:
                     return False
        else:
            return True
