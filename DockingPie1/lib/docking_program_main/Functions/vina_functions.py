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
from subprocess import check_output
import statistics
import re

from pathlib import Path

from lib.docking_program_main.Functions.threads import Protocol_exec_dialog
from lib.docking_program_main.Functions.handle_widgets import HandleWidgets

from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.PDBIO import Select


class Vina_docking():


    """
    A class to represent the Vina docking process

    """


    def __init__(self,
                 tab, main,
                 ligand,
                 receptor,
                 tmp_dir,
                 flex_residues = "",
                 ligands_to_dock = 1,
                 receptors_to_dock = 1,
                 use_flex_vina = False,
                 exhaustiveness = 8,
                 energy = 3,
                 scoring_function = "Standard",
                 poses = 1,
                 cavity = None,
                 cavity_list = [],
                 automatic = False):

        self.tab = tab
        self.main = main

        self.path_to_vina = self.main.path_to_vina

        ### CHECK VINA VERSION ###
        check_vina = check_output([self.path_to_vina, "--version"])
        self.vina_version = re.search("([\d.]+)", str(check_vina)).group(0)

        self.docking_completed = False
        self.interrupt = False
        self.use_flex_protocol = False

        # Initialize Standard Docking parameters
        self.receptor_to_dock = receptor
        self.poses = poses
        self.ligand_to_dock = ligand

        self.exhaustiveness = exhaustiveness
        self.energy = energy
        self.scoring_function = scoring_function

        self.cavity_to_dock = cavity_list # [x_pos, y_pos, z_pos, x, y, z, spacing] e.g. ['191.48', '171.91', '15.94', '16', '1', '16', '16.0']
        self.cavity_name = cavity

        # Initialize Flexible Docking
        if use_flex_vina:
            self.use_flex_protocol = True

            self.flex_residues = flex_residues
            self.check_valid_flex_input()

        # Initialize names and paths
        self.results_file_name = str("Run_" + str(self.main.vina_runs) + "_Vina")
        self.results_file_name_ext = str(self.results_file_name + ".pdbqt")
        self.results_file_name_sdf = str(self.results_file_name + ".sdf")

        # if float(self.vina_version[:3]) > 1.1:
        #     self.log_file_name = None
        # else:
        self.log_file_name = str(self.results_file_name + "_log.txt")

        self.ligands_to_dock = ligands_to_dock
        self.receptors_to_dock = receptors_to_dock
        # Change directory --> Vina tmp dir
        os.chdir(tmp_dir)

        self.show_resume_window()


    def check_valid_flex_input(self):
        pass


    def show_resume_window(self):

        self.run_docking_vina()


    def run_docking_vina(self):

        x_val = float(self.cavity_to_dock[3])
        y_val = float(self.cavity_to_dock[4])
        z_val = float(self.cavity_to_dock[5])
        spacing = float(self.cavity_to_dock[6])
# *float(self.cavity_to_dock[6])
# *float(self.cavity_to_dock[6])
# *float(self.cavity_to_dock[6])

        self.run_docking_vina_settings = [self.path_to_vina,
        "--receptor", str(self.receptor_to_dock + ".pdbqt"),
        "--ligand", str(self.ligand_to_dock  + ".pdbqt"),
        "--center_x", self.cavity_to_dock[0],
        "--center_y", self.cavity_to_dock[1],
        "--center_z", self.cavity_to_dock[2],
        "--size_x", str(x_val),
        "--size_y", str(y_val),
        "--size_z", str(z_val),
        "--spacing", str(spacing),
        "--out", self.results_file_name_ext,
        "--exhaustiveness", str(self.exhaustiveness),
        "--num_modes", str(self.poses),
        "--energy_range", str(self.energy),
        "--verbosity", str(2)]

        if self.scoring_function == "Vinardo":
            self.run_docking_vina_settings.extend(["--scoring", "vinardo"])

        if float(self.vina_version[:3]) > 1.1:
            pass
        else:
            self.run_docking_vina_settings.extend(["--log", self.log_file_name])

        if self.use_flex_protocol:

            # Generate pdbqt file with flexible side chains
            self.prepare_flex_receptor_path = os.path.join(self.main.config_path, "prepare_flexreceptor4.py")
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

        self.main.vina_runs += 1

        # If the


    def check_if_docking_completed(self):

        file_path = os.path.join(self.main.vina_tmp_dir, self.results_file_name_ext)

        if Path(file_path).is_file():

            if os.path.getsize(file_path):
                self.docking_completed = True
                self.main.vina_runs += 1

            else:
                self.docking_completed = False

                if len(self.ligands_to_dock) > 1 or len(self.receptors_to_dock) > 1:

                    os.remove(file_path)
                    self.main.vina_runs += 1

                else:
                    os.remove(file_path)
                    self.main.vina_runs += 1
                    QtWidgets.QMessageBox.warning(self.tab, "", str("Something went wrong during Docking. \nPlease check LOG files."))

        else:
            self.docking_completed = False
            QtWidgets.QMessageBox.warning(self.tab, "", str("Something went wrong during Docking. \nPlease check LOG files."))
            self.main.vina_runs += 1


class Vina_Parse_Results:


    """
    A class to parse the Vina results file
    """


    def __init__(self, tab, main, ligand, last_docking, results_file_name, poses, results_dict = {}, results_data = [[]]):


        self.tab = tab
        self.main = main
        self.last_docking = last_docking

        self.ligand_name = ligand

        self.results_file_name = results_file_name

        self.results_dict = results_dict

        self.results_data = results_data

        self.SCORE = float(00.00)

        self.docked_ligands_list = []

        self.check_output_file()

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


    def check_output_file(self):

        # Check whether some blank spaces in the file may lead to a misreading of the results

        src = (str(self.results_file_name + ".pdbqt"))
        dst = "temp_file.pdbqt"

        shutil.copyfile(src, dst)
        os.remove(src)

        out_file = open(str(self.results_file_name + ".pdbqt"), "wt")

        with open("temp_file.pdbqt",'r') as file:
            for line in file:
                if re.search("[a-zA-Z]", line):
                    out_file.write(line)
                else:
                    pass

        out_file.close()

        os.remove(dst)
