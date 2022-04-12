# Copyright 2022 by Serena Rosignoli. All rights reserved.
# This code is part of the DockingPie package and governed by its license. Please
# see the LICENSE file.


import os
import subprocess

# PyMOL.
import pymol
from pymol import cmd
from pymol.Qt import QtWidgets, QtCore, QtGui

# Biopython modules
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.PDBIO import Select

import math
#from openbabel import pybel
import sys

try:
    from spyrmsd import io, rmsd
except:
    pass

from lib.docking_program_main.Functions.threads import Protocol_exec_dialog



class ObjectParser():


    """
    When an Object is imported from PyMOL, some information are extracted and stored.
    """


    def __init__(self, file_name, file_path, is_receptor = False):

        # Define a list of aminoacids Biopython IDs
        self.aa_list = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]

        # Name of the PyMOL obj, which is the same of the pdb file in tmp dir
        self.file_name = file_name

        # Path where the file has been saved
        self.file_path = file_path

        # Inizialize the number of states found in PyMOL for the Object
        self.num_states = 1

        # It defines wether to parse as a receptor or as a ligand.
        self.is_receptor = is_receptor

        # Several Docking Programs are suited also for docking on RNA or DNA molecules (e.g. RxDock). The receptor_types defines if the Object is a PROTEIN or a POLYNUCLEOTIDE
        self.receptor_type = " "

        # Inizialize a list to store the residues of the Receptor
        self.residues_list = []

        # Inizialize a list to store the heteroresidues of the Receptor
        self.heteroresidues_list = []

        if self.is_receptor: # parse receptor
            self.parse_receptor_object()
        else: # parse ligand
            self.parse_ligand_object()


    def parse_ligand_object(self):
        self.num_states = cmd.count_states(self.file_name)


    def parse_receptor_object(self):

        # Counts the states of the PyMOL Object
        self.num_states = cmd.count_states(self.file_name)

        # Open the Receptor file
        parsed_file_handle = open(self.file_path, "r")
        # Creates a biopython 'Structure' object and starts to take informations from it.
        self.parsed_biopython_structure = PDBParser(PERMISSIVE=1, QUIET=True).get_structure(self.file_name, parsed_file_handle)
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

                    # For HETATM residues.
                    if hetfield[0] == "H":
                        self.heteroresidues_list.append(str(hetfield + str(pdb_position)))
                    # For water molecules.
                    elif hetfield == "W":
                        pass
                    else:
                        self.residues_list.append(resname)


        ### Get the receptor type ###

        found = False

        for aa in self.aa_list:
            if str(self.residues_list[0]) == str(aa):
                self.receptor_type = "PROTEIN"
                found = True
                if found:
                    break

        if not found:
            self.receptor_type = "POLYNUCLEOTIDE"


# class Calculate_RMSD_consensus:
#
#
#     def __init__(self, tab, pdb_reference_file, sdf_docked_file):
#
#         self.tab = tab
#
#         self.pdb_reference_file = pdb_reference_file
#         self.sdf_docked_file = sdf_docked_file
#
#         self.crystal = next(pybel.readfile("pdb", self.pdb_reference_file))
#
#         self.rmsd_list = []
#
#         self.find_automorphisms()
#
#
#     def find_automorphisms(self):
#
#         # Find automorphisms involving only non-H atoms
#         self.mappings = pybel.ob.vvpairUIntUInt()
#         bitvec = pybel.ob.OBBitVec()
#         self.lookup = []
#         for i, atom in enumerate(self.crystal):
#             if not atom.OBAtom.GetAtomicNum() == 1:
#                 bitvec.SetBitOn(i+1)
#                 self.lookup.append(i)
#         success = pybel.ob.FindAutomorphisms(self.crystal.OBMol, self.mappings, bitvec)
#
#         self.find_RMSD()
#
#
#     def find_RMSD(self):
#
#         # Find the RMSD between the crystal pose and each docked pose
#         xtalcoords = [atom.coords for atom in self.crystal if not atom.OBAtom.GetAtomicNum() == 1]
#         for i, dockedpose in enumerate(pybel.readfile("sdf", self.sdf_docked_file)):
#             #print("ligand %d" % (i+1))
#             posecoords = [atom.coords for atom in dockedpose  if not atom.OBAtom.GetAtomicNum() == 1]
#             minrmsd = 999999999999
#             for mapping in self.mappings:
#                 automorph_coords = [None] * len(xtalcoords)
#                 for x, y in mapping:
#                     automorph_coords[self.lookup.index(x)] = xtalcoords[self.lookup.index(y)]
#                 mapping_rmsd = self.rmsd(posecoords, automorph_coords)
#                 if mapping_rmsd < minrmsd:
#                     minrmsd = mapping_rmsd
#
#             print(minrmsd)
#             self.rmsd_list.append(minrmsd)
#             print("%d\t%.2f" % ((i+1), minrmsd))
#
#
#     def squared_distance(self, coordsA, coordsB):
#         """Find the squared distance between two 3-tuples"""
#         sqrdist = sum( (a-b)**2 for a, b in zip(coordsA, coordsB) )
#         return sqrdist
#
#
#     def rmsd(self, allcoordsA, allcoordsB):
#         """Find the RMSD between two lists of 3-tuples"""
#         deviation = sum(self.squared_distance(atomA, atomB) for
#                         (atomA, atomB) in zip(allcoordsA, allcoordsB))
#         return math.sqrt(deviation / float(len(allcoordsA)))
#

class Calculate_RMSD:


    def __init__(self, tab, pdb_reference_file, sdf_docked_file, warning = True):

        self.tab = tab

        self.pdb_reference_file = pdb_reference_file
        self.sdf_docked_file = sdf_docked_file

        self.rmsd_computed = False

        self.rmsd_list = []

        self.warning = warning

        try: 
            self.find_rmsd()
        except NameError:
            print("Openbabel is required for RMSD computing")


    def find_rmsd(self):


        ref = io.loadmol(self.pdb_reference_file)						# Read crystal pose
        mols = io.loadallmols(self.sdf_docked_file)					# Read poses
        ref.strip()											# Strip Hydrogens of X-ray
        for mol in mols:									# Strip Hydrogens of pose
            mol.strip()

        #The spyrmsd class needs atomic coordinates, atomic number and the molecular adjacency
        # matrix to compute the standard RMSD with spyrmsd.rmsd.symmrmsd.
        coords_ref = ref.coordinates
        anum_ref = ref.atomicnums
        adj_ref = ref.adjacency_matrix

        coords = [mol.coordinates for mol in mols]
        anum = mols[0].atomicnums
        adj = mols[0].adjacency_matrix

        # With this information we can easily
        # compute the RMSD between the reference molecule and all other molecules

        try:
            RMSD = rmsd.symmrmsd(coords_ref, coords, anum_ref, anum, adj_ref, adj)
            i = 0

            for RMSD_Value in RMSD:
                self.rmsd_list.append(round(RMSD_Value, 3))
                i += 1

            self.rmsd_computed = True

        except AssertionError:
            if self.warning:
                QtWidgets.QMessageBox.warning(self.tab.docking_programs_child_tabs, "Different Molecule", str("To compute RMSD value, the molecule must be the same"))
            self.rmsd_computed = False


class Generate_Object:


    def __init__(self, tab,
    dict,
    prepared_objects_list,
    format,
    docking_program_name,
    generate_pdbqt = False,
    is_receptor = False,
    ):

        # To check if the object has been generated
        self.generated_object = False

        self.tab = tab

        # To distinguish between Docking Programs
        self.dict = dict
        self.prepared_objects_list = prepared_objects_list

        # Inizialize the format to use
        self.format = format

        # Inizialize the name of the Docking program
        self.docking_program_name = docking_program_name

        # Set to True if the Docking program needs the generation of pdbqt files (es. Autodock4, Autodock-vina)
        self.generate_pdbqt = generate_pdbqt

        # Se to True if it is a receptor, else a ligand
        self.is_receptor = is_receptor

        self.generate_receptor()


    def generate_receptor(self):

        self.generated_receptor = False

        self.cwd = os.getcwd()

        for self.strc in self.dict:

            # For each object that is checked
            obj = self.dict[self.strc]["frame"].strc_checkbox

            if obj.isChecked():

                states = cmd.count_states(obj.text().split()[0])

                if states > 1 and self.generate_pdbqt:
                    self.generate_checked_object_multiple(states)

                else:
                    self.generate_checked_object()

                #self.generated_receptor = True
                self.generated_object = True

                if self.generated_object:

                    # Add to the listwidget
                    self.listbtn = self.tab.listwidget.addItem(self.new_strc_name)

                    # Add to the list of prepared files
                    self.prepared_objects_list.append(self.new_strc_name)

                    # Load in PyMOL
                    cmd.load(self.new_receptor_file_path, self.new_strc_name)

                    if self.docking_program_name == "Vina":
                        cmd.group("Vina", members=self.new_strc_name, action='auto', quiet=1)

                    if self.docking_program_name == "RxDock":
                        cmd.group("RxDock", members=self.new_strc_name, action='auto', quiet=1)

                    if self.docking_program_name == "Smina":
                        cmd.group("Smina", members=self.new_strc_name, action='auto', quiet=1)

                    if self.docking_program_name == "ADFR":
                        cmd.group("ADFR", members=self.new_strc_name, action='auto', quiet=1)


                    # self.tab.docking_programs_child_tabs.docking_programs.statusBar().showMessage(str("Generated Object: " + self.new_strc_name), 3000)


    def generate_checked_object_multiple(self, states):

        new_name = self.create_new_name(self.strc)

        index = str((len(self.prepared_objects_list)+1))
        self.new_strc_name = str("0" + index + "_" + new_name + "_" + self.docking_program_name) # Name of the prepared receptor

        for idx, i in enumerate(range(states)):
            save_to = str(os.path.join(self.cwd, self.new_strc_name) + "_ML" + str(idx+1) + str("." + self.format))
            self.new_receptor_file_path = str(os.path.join(self.cwd, self.new_strc_name) + "_ML" + str(idx+1) + str("." + self.format))

            cmd.save(save_to, self.strc, format = str(self.format), state = idx+1)

            self.prepare_receptor_path = os.path.join(self.tab.docking_programs_child_tabs.docking_programs.config_path, "prepare_receptor4.py")
            self.prepare_ligand_path = os.path.join(self.tab.docking_programs_child_tabs.docking_programs.config_path, "prepare_ligand4.py")
            save_to_pdb = str(os.path.join(self.cwd, self.new_strc_name) + "_ML" + str(idx+1) + ".pdb")
            out_name = str(os.path.join(self.cwd, self.new_strc_name) + "_ML" + str(idx+1))
            name = str(self.new_strc_name) + "_ML" + str(idx+1)
            # self.receptor_file_path = os.path.join(self.cwd, str(self.strc + ".pdb"))

            # Save the main file from PyMOL, its is done to be sure to work with the last modified object in PyMOL
            cmd.save(save_to_pdb, self.strc, format = 'pdb', state = idx+1)

            if self.is_receptor:

                self.generate_receptor_pdbqt(self.strc)

                #self.generate_receptor_pdbqt(strc)

            else: # is ligand

                self.generate_ligand_pdbqt(name, out_name)

                #self.generate_ligand_pdbqt(strc)

                old_name = str(name + ".pdbqt")
                new_name = str(self.new_strc_name + ".pdbqt")
                os.rename(old_name, new_name)


    def generate_checked_object(self):

        new_name = self.create_new_name(self.strc)

        index = str((len(self.prepared_objects_list)+1))
        self.new_strc_name = str("0" + index + "_" + new_name + "_" + self.docking_program_name) # Name of the prepared receptor
        self.new_receptor_file_path = str(os.path.join(self.cwd, self.new_strc_name) + str("." + self.format))

        # Save the main file from PyMOL, its is done to be sure to work with the last modified object in PyMOL
        #cmd.h_add(strc)
        cmd.save(self.new_receptor_file_path, self.strc, format = str(self.format), state = 0)

        ##
        ## For those programs that need the generation of the pdbqt file (Vina, Autodock etc ...)
        ##
        if self.generate_pdbqt:

            self.prepare_receptor_path = os.path.join(self.tab.docking_programs_child_tabs.docking_programs.config_path, "prepare_receptor4.py")
            self.prepare_ligand_path = os.path.join(self.tab.docking_programs_child_tabs.docking_programs.config_path, "prepare_ligand4.py")
            self.receptor_file_path = os.path.join(self.cwd, str(self.strc + ".pdb"))

            # Save the main file from PyMOL, its is done to be sure to work with the last modified object in PyMOL
            cmd.save(self.receptor_file_path, self.strc, format = 'pdb', state = 0)

            if self.is_receptor:

                self.generate_receptor_pdbqt(self.strc)

            else: # is ligand

                self.generate_ligand_pdbqt(self.strc, self.strc)

                old_name = str(self.strc + ".pdbqt")
                new_name = str(self.new_strc_name + ".pdbqt")
                os.rename(old_name, new_name)



    def generate_receptor_pdbqt(self, strc):

        self.receptors_settings = ["python", self.prepare_receptor_path, "-r", str(strc + ".pdb"), "-o", str(self.new_strc_name + ".pdbqt"), "-v"]

        if self.tab.structure_frame.add_h.isChecked():
            self.receptors_settings.extend(["-A", "hydrogens"])
            cmd.h_add(strc)

        if self.tab.structure_frame.remove_nonstd.isChecked():
            self.receptors_settings.extend(["-e" "'True'"])

        elif self.tab.structure_frame.remove_water.isChecked():
            self.receptors_settings.extend(["-U", "'waters'"])

        if sys.platform == "win32":

            try:
                subprocess.run(self.receptors_settings,
                shell = True)
                self.generated_object = True

            except subprocess.CalledProcessError as error:
                print(error)
                self.generated_object = False

        else:

            try:
                subprocess.run(self.receptors_settings)
                self.generated_object = True

            except subprocess.CalledProcessError as error:
                print(error)
                self.generated_object = False


    def generate_ligand_pdbqt(self, strc, output_name):


        self.ligands_settings = ["python", self.prepare_ligand_path, "-l", str(strc + ".pdb"), "-v"]

        if self.tab.structure_frame.add_h.isChecked():
            self.ligands_settings.extend(["-A", "hydrogens"])
            cmd.h_add(self.strc)

        if self.tab.structure_frame.none_torsions.isChecked():
            self.ligands_settings.extend(["-Z"])

        elif self.tab.structure_frame.all_torsions.isChecked():
            self.ligands_settings.extend(["-B", "'amide'", "-B" "'guanidinium'"])

        if sys.platform == "win32":

            try:

                process = subprocess.run(self.ligands_settings,
                shell = True)

                # os.rename(str(strc + ".pdbqt"), output_name)

                self.generated_object = True

            except subprocess.CalledProcessError as error:
                print(error)
                self.generated_object = False

        else:
            try:

                process = subprocess.run(self.ligands_settings)

                # os.rename(str(strc + ".pdbqt"), output_name)

                self.generated_object = True

            except subprocess.CalledProcessError as error:
                print(error)
                self.generated_object = False



    def create_new_name(self, strc):

        self.new_name = strc.replace("_", "-")

        return self.new_name


class Split_pdbqt_ligand_files():

    def __init__(self, tab,
    file_name,
    path_to_file,
    ):

        self.file_name = file_name
        self.path_to_file = path_to_file
