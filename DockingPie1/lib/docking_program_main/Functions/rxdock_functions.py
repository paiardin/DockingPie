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



class RxDock_Functions():


    """
    A class just to store functions specific for RxDock - WorkInProgress
    """


    def import_grid_in_pymol(self, file_to_load, pymol_object_name):

        cmd.load(file_to_load, pymol_object_name)
        cmd.group("RxDock", members=pymol_object_name, action='auto', quiet=1)

        isomesh_obj = str("Isomesh_" + pymol_object_name)
        cmd.isomesh(isomesh_obj, pymol_object_name, 0.99)
        cmd.group("RxDock", members=isomesh_obj, action='auto', quiet=1)

    def is_valid_grd_file(self, file_path):
        self.is_valid = True

    def is_valid_as_file(self, file_path):
        self.is_valid = True

    def is_valid_constrains_file(self, file_path):
        self.is_valid = True



class RxDock_Cavity():


    """
    A class to represent the process of generating a Cavity with RxDock

    - reference_receptor and reference_ligand get the current Receptor and Ligand to use as a reference to generate the Cavity

    - two_spheres_method and reference_ligand_method are boolean to specify the type of protocol to use to generate the Cavity

    - initial_prm_file_path specifies the path to the config directory where a standard prm_file is located
    """


    def __init__(self, tab, main,
                 reference_receptor,
                 two_spheres_method,
                 reference_ligand_method,
                 tmp_dir,
                 x,
                 y,
                 z,
                 progress_bar = True,
                 reference_ligand = "",
                 radius_spinbox = 10.00,
                 small_sphere_value = 1.50,
                 large_sphere_value = 4.00):

        self.tab = tab
        self.main = main

        # Initialize Parameters
        self.reference_receptor = reference_receptor
        self.reference_ligand = reference_ligand

        # Initialize Docking protocols to use
        self.two_spheres_method = two_spheres_method
        self.reference_ligand_method = reference_ligand_method

        # Get Parameters
        if self.two_spheres_method:
            self.radius_value = str(radius_spinbox)
            self.small_sphere_value = str(small_sphere_value)
            self.large_sphere_value = str(large_sphere_value)

        self.simple_docking = False
        self.pharma_restrains = False
        self.tethered_docking = False

        # Initialize names and paths
        self.initial_prm_file_path = self.main.path_to_cavity

        # Position
        self.x = x
        self.y = y
        self.z = z

        # Change directory --> RxDock tmp dir
        os.chdir(tmp_dir)

        if progress_bar:

            p_dialog = Protocol_exec_dialog(app=self.tab, docking_pie=self.tab,
                                            function=self.generate_cavity,
                                            args=(),
                                            wait_start=0.4, wait_end=0.4,
                                            lock=True,
                                            stdout_silence=True,
                                            title="Running",
                                            label_text="RxDock is generating the cavity. Please wait")
            p_dialog.exec_()

        else:
            self.generate_cavity()


    def generate_cavity(self):

        # Create the PrmFile
        self.prm_file = PrmFile(self,
        prm_file_name = "prm_file",
        counter = len(self.main.generated_cavity),
        x = self.x,
        y = self.y,
        z = self.z)

        f = open(self.prm_file.prm_file_name + "_LOG.txt", "w")
        # Run RxDock "rbcavity" func
        subprocess.run(["rbcavity", "-W", "-d", "-r", str(self.prm_file.prm_file_name)], stdout = f)
        f.close()



class RxDock_docking():


    """
    A class to represent the RxDock docking process

    """


    def __init__(self, tab, main,
                 receptor,
                 pharma_restrains,
                 tethered_docking,
                 poses_box,
                 cavity_to_dock,
                 cavity_name,
                 use_water,
                 ligand,
                 trans_mode,
                 rot_mode,
                 die_mode,
                 tmp_dir,
                 protein_segments_to_exclude = [],
                 pharma_list = [],
                 cavity = None):

        self.tab = tab
        self.main = main
        self.thread = tab

        self.rxdock_tmp_dir = tmp_dir

        self.docking_completed = False
        self.interrupt = False

        # Initialize Docking protocols to use
        self.two_spheres_method = False
        self.reference_ligand_method = False
        self.use_water = False

        self.pharma_restrains = pharma_restrains
        self.tethered_docking = tethered_docking

        self.set_docking_protocol()

        # Initialize Standard Docking parameters
        self.receptor_to_dock = receptor
        self.poses = poses_box
        self.cavity_to_dock = cavity_to_dock
        self.cavity_name = cavity_name
        self.ligand_to_dock = ligand

        # If Tethered Docking is performed the ligand for docking is modified, thus it is called as "name_of_the_ligand_tethered"
        if self.tethered_docking:
            self.ligand_to_dock_tethered = str(self.ligand_to_dock + "tethered")

        # If pharma_restrains protocol is used, create the needed file
        if self.pharma_restrains:
            self.pharma_list = pharma_list
            self.create_pharma_restrains_file()

        self.protein_segments_to_exclude = protein_segments_to_exclude
        if not self.protein_segments_to_exclude:
            self.exclude_segments = False
        else:
            self.exclude_segments = True
            self.create_protein_segments_string()

        if use_water:
            self.use_water = True
            self.create_structural_water_file()

        # Initialize names and paths
        r_name = str("Run_" + str(self.main.rxdock_runs) + "_RxDock")
        self.results_file_name = r_name

        self.log_file_name = str(self.results_file_name + "_log.txt")

        self.list_of_ligands_to_dock = []
        try:
            states = cmd.count_states(self.ligand_to_dock)
            self.list_of_ligands_to_dock = []

            if states > 1:
                self.split_input_ligand_file()
                self.results_file_name_ext = []
                for idx, i in enumerate(range(states)):
                    self.results_file_name_ext.extend([r_name + "_ML" + str(idx+1) + ".sd"])

            else:
                self.results_file_name_ext = str(self.results_file_name + ".sd")
        except:
            self.results_file_name_ext = str(self.results_file_name + ".sd")

        self.initial_prm_file_path = self.main.path_to_cavity

        # Initialize additional parameters
        self.trans_mode = trans_mode
        self.rot_mode = rot_mode
        self.die_mode = die_mode

        # Change directory --> RxDock tmp dir
        os.chdir(self.rxdock_tmp_dir)

        self.show_resume_window()


    def check_if_docking_completed(self):

        self.file_path = os.path.join(self.main.rxdock_tmp_dir, self.results_file_name_ext)

        if Path(self.file_path).is_file():

            if os.path.getsize(self.file_path):
                self.docking_completed = True
                self.main.rxdock_runs += 1

            else:
                self.docking_completed = False
                QtWidgets.QMessageBox.warning(self.tab, "", str("Something went wrong during Docking. \nPlease check LOG files."))
                os.remove(self.file_path)
                self.main.rxdock_runs += 1

        else:
            self.docking_completed = False
            QtWidgets.QMessageBox.warning(self.tab, "", str("Something went wrong during Docking. \nPlease check LOG files."))
            self.main.rxdock_runs += 1


    def split_input_ligand_file(self):

        rxdock_tmp_dir_path = self.main.rxdock_tmp_dir

        for idx, state in enumerate(range(cmd.count_states(self.ligand_to_dock))):

            save_to = str(os.path.join(rxdock_tmp_dir_path, self.ligand_to_dock) + "_ML" + str(idx+1))
            cmd.save(save_to, self.ligand_to_dock, format = 'sdf', state = idx+1)
            self.list_of_ligands_to_dock.append(save_to)


    def create_structural_water_file(self):

        list_of_waters = []

        input = open(self.receptor_to_dock.split("_")[1] + ".pdb")

        for line in input:
            if line.startswith("HETATM"):
                col = line.split()[3]
                if col == "HOH":
                    list_of_waters.append(line)

        input.close()

        output = open("inputreceptorwater.pdb", "w")

        for elements in list_of_waters:
            output.write(elements)

        output.close()

        if os.stat("inputreceptorwater.pdb").st_size == 0:
            self.use_water = False
            print("water not found")


    def create_protein_segments_string(self):

        list_of_ids = []

        # Open the Receptor file
        parsed_file_handle = open(self.receptor_to_dock.split("_")[1] + ".pdb", "r")
        # Creates a biopython 'Structure' object and starts to take informations from it.
        self.parsed_biopython_structure = PDBParser(PERMISSIVE=1, QUIET=True).get_structure(self.receptor_to_dock.split("_")[1] + ".pdb", parsed_file_handle)
        # Close the Receptor File
        parsed_file_handle.close()

        ### Store residues and heteroresidues information ###

        # Iterate over the biopython 'Structure' Object to get information
        for model in self.parsed_biopython_structure.get_list():
            for chain in model.get_list():
                for residue in chain:
                    # Gets the 3 letter name of the current residue.
                    resname = residue.get_resname()
                    # get_id() returns something like: ('H_SCN', 1101, ' '). The first item is
                    # the hetfield: 'H_SCN' for an HETRES, while ' ' for a normal residue. The
                    # second item is the id of the residue according to the PDB file.
                    hetfield, pdb_position = residue.get_id()[0:2]

                    list_of_ids.append(hetfield[0])

        #######

        self.protein_segments_to_exclude_string = ""

        for items in list_of_ids:
            for exclude_itm in self.tab.protein_segments_to_exclude:
                if items == exclude_itm:
                    pass
                else:
                    item_to_write = str(exclude_itm.split("_")[1])
                    self.protein_segments_to_exclude_string += item_to_write + ","



    def show_resume_window(self):

        self.run_docking_rxdock()


    def create_pharma_restrains_file(self):

        ### To create the file where the pharmacophoric restrains coordinates are stored ###
        tmp_list = []

        # Get pharmacophoric restrains info
        for i in range(self.pharma_list.count()):
            item_to_write = str((self.pharma_list.item(i)).text() + "\n")
            tmp_list.append(item_to_write)

        os.chdir(self.rxdock_tmp_dir)

        # Path to pharma file, it is created a new one each time, thus the name is always the same
        pharma_file_path = os.path.join(self.main.rxdock_tmp_dir, "pharma.const")

        # If the pharma file already exists, remove it.
        try:
            os.remove(pharma_file_path)
        except:
            pass

        # Create a new pharma file with pharmacophoric restrains' info
        file = open("pharma.const", "w")
        file.writelines(tmp_list)
        file.close()


    def set_docking_protocol(self):

                ### To see which protocol is specified by the user ###

        if self.pharma_restrains or self.tethered_docking:
            self.simple_docking = False
        else:
            self.simple_docking = True


    def run_docking_rxdock(self):

                ### To run the rxDock Docking Process ###

        # The name of the prm_file used for docking is taken from the name of the prm_file used to previously generate the cavity
        as_file_name = self.cavity_to_dock.replace("_cav1", "")

        # Create the prm_file
        self.prm_file = PrmFile(self,
        prm_file_name = as_file_name)

        os.chdir(self.rxdock_tmp_dir)

        if self.list_of_ligands_to_dock:

            for ligand in self.list_of_ligands_to_dock:

                subprocess.run(["rbdock",
                "-i", ligand,
                "-o",  "results_tmp_name",
                "-r", self.prm_file.prm_file_name,
                "-p",
                "dock.prm",
                "-n", self.poses],
                check=True)

                if Path("results_tmp_name.sd").is_file():

                    f = open(str(self.results_file_name + "_" + ligand.split("_")[-1] + ".sd"), "w")
                    subprocess.run(["sdsort", "-n", "-fSCORE", "results_tmp_name.sd"], stdout = f, check = True)
                    f.close()

        else:

            self.run_docking_rxdock_settings = ["rbdock",
            "-i",
            str(self.ligand_to_dock + ".sdf"),
            "-o",
            "results_tmp_name",
            "-r",
            self.prm_file.prm_file_name,
            "-p",
            "dock.prm",
            "-n",
            self.poses]

            # try:
            #     output = subprocess.check_output(self.run_docking_rxdock_settings)
            # except subprocess.CalledProcessError as e:
            #     output = e.output
            #     print(output)
            #     pass
            #
            # #subprocess.run(self.run_docking_rxdock_settings, check = True)



            # Run RxDock rbdock func
            # subprocess.run(["rbdock",
            # "-i", str(self.ligand_to_dock  + ".sdf"),
            # "-o",  "results_tmp_name",
            # "-r", self.prm_file.prm_file_name,
            # "-p",
            # "dock.prm",
            # "-n", self.poses],
            # check=True)

        #subprocess.run(["obabel", "results_tmp_name.sd", "-O", str(self.results_file_name + ".sd"), "--sort", "'>  <SCORE>'"])
        # # # # # sdsort -n -f'SCORE' <output>.sd > <sorted_output>.sd

        #
        # path_to_sdsorter = os.path.join(self.tab.docking_programs_child_tabs.docking_programs.config_path, "sdsorter.static")
        # subprocess.run([path_to_sdsorter, "-sort", "'>  <SCORE>'", "results_tmp_name.sd", str(self.results_file_name + ".sd")])

        # subprocess.run(["sdsort"])
            if Path("results_tmp_name.sd").is_file():

                if sys.platform == "darwin":

                    path_to_sdsorter = os.path.join(self.main.path_to_sdsorter)
                    subprocess.run([path_to_sdsorter, "-sort", "'>  <SCORE>'", "results_tmp_name.sd", str(self.results_file_name + ".sd")])

                else:
                # subprocess.run(["obabel", "-i", "sd", "results_tmp_name.sd", "-o", "sdf", "-O", "results_tmp_name.sdf"])
                # subprocess.run(["obabel", "results_tmp_name.sdf", "-O", str(self.results_file_name + ".sdf"), "--sort", "SCORE"])

                    f = open(str(self.results_file_name + ".sd"), "w")
                    subprocess.run(["sdsort", "-n", "-fSCORE", "results_tmp_name.sd"], stdout = f, check = True)
                    f.close()

        #self.docking_completed = True
        self.main.rxdock_runs += 1



class PrmFile():


    """
    A class to represent the prm file that is created by RxDock to set the desired parameters during the
    generation of cavities and the running of Docking processes.

    note: since the same file format is used for both type of processes it is necessary to set as 'True'
    the method that we wish to use.

    Cavity generation methods:
        - 'two-spheres method'
        - 'reference_ligand'

    Docking Protocols:
        - 'simple_docking'
        - 'tethered_docking'
        - 'pharma_restrains' ('pharma_list' is needed to generate the pharma.restr file where pharmacophoric restrains are defined)

    Additional options:

        - use_water -- To use structural water during the Docking process
        - input_file_path -- path to the standard_prm_file located in config directory
        - receptor_name -- Name of the Receptor. It is used both for the generation of the cavity with the Reference Ligand Method and for the Docking Process
        - reference_ligand_name -- Name of the ligand to use as a reference when creating the Cavity (Note: the input ligand for the Docking process is not specified in the prm file!)
    """


    def __init__(self, docking_process,
    x = "", y = "", z = "",
    prm_file_name = None,
    counter = 0):

        self.docking_process = docking_process

        # Initialize prm_file name
        self.prm_file_name = prm_file_name

        # To keep track of the number of cavities generated, number from which the prm_file name is set
        self.counter = str(counter)

        # Set the path to the constraints file
        self.pharma_file_path = self.docking_process

        self.generate_prm_file(prm_file_name, x, y, z)


    def generate_prm_file(self, prm_file_name, x, y, z):

        if self.docking_process.simple_docking:

            input = open(self.docking_process.initial_prm_file_path, "rt")
            output = open(self.prm_file_name, "wt")

            for line in input:
                if line.startswith("RECEPTOR_SEGMENT_NAME") and self.docking_process.exclude_segments:
                    output.write(line.replace("RECEPTOR_SEGMENT_NAME", str("RECEPTOR_SEGMENT_NAME " + self.docking_process.protein_segments_to_exclude_string)))
                elif line.startswith("RECEPTOR_FILE"):
                    output.write(line.replace("RECEPTOR_FILE inputreceptor", str('RECEPTOR_FILE ' + str(self.docking_process.receptor_to_dock))))
                elif line.startswith("    SMALL_SPHERE"):
                    output.write(line.replace("SMALL_SPHERE 1.0", "SMALL_SPHERE 1.0\n    LARGE_SPHERE 4.0"))
                elif line.startswith("<END FILE>") and self.docking_process.use_water:
                    output.write(line.replace("<END FILE>", "SECTION SOLVENT\n    FILE inputreceptorwater.pdb\n    TRANS_MODE TETHERED\n    ROT_MODE TETHERED\n    MAX_TRANS 1.0\n    MAX_ROT 30.0\n    OCCUPANCY 0.5\nEND_SECTION\n<END FILE>"))
                else:
                    output.write(line)


        if self.docking_process.pharma_restrains:

            input = open(self.docking_process.initial_prm_file_path, "rt")
            output = open(self.prm_file_name, "wt")

            for line in input:
                if line.startswith("RECEPTOR_FILE"):
                    output.write(line.replace("RECEPTOR_FILE inputreceptor", str('RECEPTOR_FILE ' + str(self.docking_process.receptor_to_dock))))
                elif line.startswith("    SMALL_SPHERE"):
                    output.write(line.replace("SMALL_SPHERE 1.0", "SMALL_SPHERE 1.0\n    LARGE_SPHERE 4.0"))
                elif line.startswith("<END FILE>"):
                    output.write(line.replace("<END FILE>", str("SECTION PHARMA\n    SCORING_FUNCTION RbtPharmaSF\n    WEIGHT 1.0\n    CONSTRAINTS_FILE pharma.const \nEND_SECTION\n<END FILE>")))
                else:
                    output.write(line)


        if self.docking_process.tethered_docking:

            self.prepare_ligand_for_tethered_docking()

            input = open(self.docking_process.initial_prm_file_path, "rt")
            output = open(self.prm_file_name, "wt")

            for line in input:
                if line.startswith("RECEPTOR_FILE"):
                    output.write(line.replace("RECEPTOR_FILE inputreceptor", str('RECEPTOR_FILE ' + str(self.docking_process.receptor_to_dock))))
                elif line.startswith("    SMALL_SPHERE"):
                    output.write(line.replace("SMALL_SPHERE 1.0", "SMALL_SPHERE 1.0\n    LARGE_SPHERE 4.0"))
                elif line.startswith("<END FILE>"):
                    output.write(line.replace("<END FILE>", str("SECTION LIGAND\n    TRANS_MODE " + self.docking_process.trans_mode + "\n    ROT_MODE " + self.docking_process.rot_mode + "\n    DIHEDRAL_MODE " + self.docking_process.die_mode + "\n    MAX_TRANS 1.0\n    MAX_ROT 30.0\n    MAX_DIHEDRAL 30.0\n    \nEND_SECTION\n<END FILE>")))
                else:
                    output.write(line)


        if self.docking_process.two_spheres_method:

            # Get the coordinates of the reference object, thus get the center of the cavity

            # x = str(self.docking_process.tab.x_scroll.value())
            # y = str(self.docking_process.tab.y_scroll.value())
            # z = str(self.docking_process.tab.z_scroll.value())

            coords = str("(" + str(x) + "," + str(y) + "," + str(z) + ")")

            # Name of the prm file
            self.prm_file_name = str(prm_file_name + self.counter)

            input = open(self.docking_process.initial_prm_file_path, "rt")
            output = open(self.prm_file_name, "wt")

            for line in input:

                if line.startswith("RECEPTOR_FILE"):
                    output.write(line.replace("RECEPTOR_FILE inputreceptor", str('RECEPTOR_FILE ' + str(self.docking_process.reference_receptor))))
                elif line.startswith("    SITE_MAPPER"):
                    output.write(line.replace("SITE_MAPPER algorithm", "SITE_MAPPER RbtSphereSiteMapper"))
                elif line.startswith("    REF_MOL"):
                    output.write(line.replace("REF_MOL referenceligand.sdf", str('CENTER ' + str(coords))))
                elif line.startswith("    RADIUS"):
                    output.write(line.replace("RADIUS 2.0", "RADIUS " + self.docking_process.radius_value))
                elif line.startswith("    SMALL_SPHERE"):
                    output.write(line.replace("SMALL_SPHERE 1.0", "SMALL_SPHERE " + self.docking_process.small_sphere_value + "\n    LARGE_SPHERE " + self.docking_process.large_sphere_value))
                elif line.startswith("    MIN_VOLUME 100"):
                    output.write(line.replace("MIN_VOLUME 100\n    MAX_CAVITIES 1\n    VOL_INCR 0.0\n    GRIDSTEP 0.5", ""))
                else:
                    output.write(line)

            input.close()
            output.close()


        if self.docking_process.reference_ligand_method:

            self.prm_file_name = str(prm_file_name + self.counter)

            input = open(self.docking_process.initial_prm_file_path, "rt")
            output = open(self.prm_file_name, "wt")

            for line in input:

                if line.startswith("RECEPTOR_FILE"):
                    output.write(line.replace("RECEPTOR_FILE inputreceptor", str('RECEPTOR_FILE ' + str(self.docking_process.reference_receptor))))
                elif line.startswith("    SITE_MAPPER"):
                    output.write(line.replace("SITE_MAPPER algorithm", "SITE_MAPPER RbtLigandSiteMapper"))
                elif line.startswith("    REF_MOL"):
                    output.write(line.replace("REF_MOL referenceligand", str('REF_MOL ' + str(self.docking_process.reference_ligand))))
                else:
                    output.write(line)

            input.close()
            output.close()


    def prepare_ligand_for_tethered_docking(self):

        self.reference_ligand_tether = self.docking_process.tab.tethered_ligand_combo.currentText()

        subprocess.run(["sdtether",
        self.reference_ligand_tether,
        self.docking_process.ligand_to_dock,
        self.docking_process.ligand_to_dock_tethered,
        str('SMARTS')
        ])



class RxDock_parse_results():


    """
    A class to parse the RxDock results file
    """


    def __init__(self, tab, main,
    results_file_name,
    results_dict,
    poses, ligand,
    results_data = [[]]):

        self.tab = tab
        self.main = main

        self.ligand_name = ligand

        self.results_file_name = results_file_name

        self.results_dict = results_dict

        self.results_data = results_data

        self.poses = poses

        self.SCORE = float(00.00),
        self.score_INTER = float(00.00),
        self.score_INTRA = float(00.00),

        self.docked_ligands_list = []

        self.list_of_list = []
        self.list_of_all_scores = []

        self.poses_list = []
        self.input = open(str(self.results_file_name + ".sd"), "rt")
        self.parse_file(ligand)

        # Input format for the table
        self.results_data = self.list_of_list
        # Save all the scores in a list
        self.TOT_SCORE = self.list_of_all_scores


    def parse_file(self, ligand):

        n = 0

        for line in self.input:

            if line.startswith(">  <Name>"):
                value = next(self.input)
                name = str(value.strip()) ## use this when the plugin will support multiligand files
                name = ligand

                if self.check_name(name): # if True, the ligand is the same, the state is different
                    self.poses_list.append(name)
                    self.NAME = name
                    self.POSE = str(len(self.poses_list))

                else: # if False, the ligand is changed, thus the tmp list (poses_list) must be empty
                    self.poses_list = []
                    self.poses_list.append(name)
                    self.NAME = name
                    self.POSE = str(len(self.poses_list))

                self.docked_ligands_list.append(self.NAME)
                self.docked_ligands_list.append(str(name + "_all_poses"))

            if line.startswith(">  <SCORE>"):
                value = next(self.input)
                self.SCORE = float(value.strip())

            elif line.startswith(">  <SCORE.INTER>"):
                value = next(self.input)
                self.score_INTER = float(value.strip())

            elif line.startswith(">  <SCORE.INTRA>"):
                value = next(self.input)
                self.score_INTRA = float(value.strip())

            elif line.startswith("$$$$"):
                new_list = [self.NAME, self.POSE, str(self.SCORE), str(self.score_INTER), str(self.score_INTRA)]
                self.list_of_list.append(new_list)

                self.list_of_all_scores.append(self.SCORE)


    def check_name(self, name):

        if self.poses_list:

            for i in self.poses_list:
                if str(i) == name:
                    return True
                else:
                     return False
        else:
            return True
