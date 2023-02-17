# Copyright 2022 by Serena Rosignoli. All rights reserved.
# This code is part of the DockingPie package and governed by its license. Please
# see the LICENSE file.


import os
import subprocess
import re

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
from lib.docking_program_main.docking_program_gui.new_windows import NewWindow
from lib.docking_program_main.Functions.handle_widgets import HandleWidgets

def check_pymol_object_chemical_type(obj):
    pass


def get_current_sele(tab, main, widget,
                     x_scroll,
                     y_scroll,
                     z_scroll,
                     x_scroll_vis,
                     y_scroll_vis,
                     z_scroll_vis,
                     spacing_scroll_vis):

    # Get the current object
    sel = widget.currentText()

    empty = HandleWidgets.combobox_check_if_empty(self = main,
    widgets_list = [widget])

    if empty:
        pass
    else:
        try:
            tab.show_box_func(main, widget,
                             x_scroll,
                             y_scroll,
                             z_scroll,
                             x_scroll_vis,
                             y_scroll_vis,
                             z_scroll_vis,
                             spacing_scroll_vis)
        except:
            pass


def update_widget_with_pymol_object(main,
                                    widget,
                                    clear_widget = True,
                                    polymer = True,
                                    small_molecule = True,
                                    update_name_if_states = True,
                                    selections = False):

    # List all objects and selection in PyMOL
    if selections:
        selections_list = [str(obj) for obj in cmd.get_names("objects") + cmd.get_names("selections")]
    else:
        selections_list = [str(obj) for obj in cmd.get_names("objects")]

    # Check for importable objects
    if selections_list == []:
        QtWidgets.QMessageBox.warning(main, "PyMOL is empty", str("There isn't any object to import"))
        widget.clear()
    # If importable objects are present, update the ComboBox
    else:
        if clear_widget:
            widget.clear()

        if selections:
            if clear_widget:
                widget.clear()
            for i in selections_list:
                type = cmd.get_type(i)
                if type == str("object:molecule") or type == "selection":
                    widget.addItem(i)

        else:
            polymer_list = []
            small_molecule_list = []
            for i in selections_list:
                type = cmd.get_type(i)
                if polymer:
                    if type == str("object:molecule") and not re.search("Run_", i) and cmd.count_atoms(i) > 200:
                        widget.addItem(i)
                        polymer_list.append(i)

                if small_molecule:
                    if type == str("object:molecule") and not re.search("Run_", i) and cmd.count_atoms(i) < 200:
                        if cmd.count_states(i) > 1:
                            if update_name_if_states:
                                new_name = i + "-(" + str(cmd.count_states(i)) + ")"
                                widget.addItem(new_name)
                            else:
                                widget.addItem(i)
                        else:
                            widget.addItem(i)

                        small_molecule_list.append(i)

            if polymer and not polymer_list:
                QtWidgets.QMessageBox.warning(main, "PyMOL is empty", str("There isn't any Receptor to import"))
                widget.clear()

            if small_molecule and not small_molecule_list:
                QtWidgets.QMessageBox.warning(main, "PyMOL is empty", str("There isn't any Ligand to import"))
                widget.clear()


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

        protein = False
        rna = False
        dna = False
        chain_dict = {}

        # Counts the states of the PyMOL Object
        self.num_states = cmd.count_states(self.file_name)

        # Open the Receptor file
        parsed_file_handle = open(self.file_path, "r")
        # Creates a biopython 'Structure' object and starts to take informations from it.
        self.parsed_biopython_structure = PDBParser(PERMISSIVE=1, QUIET=True).get_structure(self.file_name, parsed_file_handle)
        # Close the Receptor File
        parsed_file_handle.close()

        ### Store residues and heteroresidues information ###

        rec_type_list = []

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

                        if resname in self.aa_list:
                            protein = True
                            string = "PROTEIN"
                        else:
                            if len(resname) == 2:
                                dna = True
                                string = "DNA"
                            else:
                                rna = True
                                string = "RNA"

                rec_type_list.append(str(chain.get_id()) + ": " + string)

        self.receptor_type = '\n'.join(rec_type_list)
        ### Get the receptor type ###

        # found = False
        #
        # for aa in self.aa_list:
        #     if str(self.residues_list[0]) == str(aa):
        #         self.receptor_type = "PROTEIN"
        #         found = True
        #         if found:
        #             break
        #
        # if not found:
        #     self.receptor_type = "POLYNUCLEOTIDE"


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
    tmp_path,
    docking_program_name,
    generate_pdbqt = False,
    pdbqt_dict = {},
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

        # Path to the temp directory of the specific docking program
        self.tmp_path = tmp_path

        # Dictionary of the pdbqt options
        self.pdbqt_dict = pdbqt_dict

        self.generate_receptor()


    def generate_receptor(self):

        self.generated_receptor = False

        # For each object that is loaded in DockingPie
        for self.strc in self.dict:

            # If the object is checked
            obj = self.dict[self.strc]["frame"].strc_checkbox
            if obj.isChecked():

                # Count the number of states
                states = cmd.count_states(obj.text().split()[0])

                # TODO DEAL WITH MULTIPLE STATES

                self.generate_checked_object()


    def generate_checked_object(self):

        # Remove underscores in the original name, to avoid name conflicts
        new_name = self.create_new_name(self.strc)

        # Get the new index from the list of prepared objects
        index = str((len(self.prepared_objects_list)+1))
        # Create the new name of type "index-<name>-<docking_program>" (e.g. 01_1ol5_Vina)
        self.new_strc_name = str("0" + index + "_" + new_name + "_" + self.docking_program_name)
        # Create new name with format (e.g. 01_1ol5_Vina.pdb)
        self.new_strc_name_format = self.new_strc_name + str("." + self.format)
        # Create path to receptor in tmp directory (e.g. <tmp_path>/01_1ol5_Vina.pdb)
        self.new_receptor_file_path = os.path.join(self.tmp_path, self.new_strc_name_format)

        ## For those programs that need the generation of the pdbqt file (Vina, Autodock etc ...)
        if self.generate_pdbqt:

            self.prepare_receptor_path = os.path.join(self.tab.docking_programs_child_tabs.docking_programs.config_path, "prepare_receptor4.py")
            self.prepare_ligand_path = os.path.join(self.tab.docking_programs_child_tabs.docking_programs.config_path, "prepare_ligand4.py")

            if self.is_receptor:

                self.generate_receptor_pdbqt(self.strc)

            else: # is ligand

                self.generate_ligand_pdbqt(self.strc, self.strc)

        else:

            # TODO: adjsust parameters also for programs out of Vina, Smina and ADFR

            # Save the main file from PyMOL, its is done to be sure to work with the last modified object in PyMOL \
            #  PyMOL saves "self.strc" obj (e.g. 1ol5) as "<tmp_path>/01_1ol5_Vina.pdb". Only the current state is saved
            cmd.save(self.new_receptor_file_path, self.strc, format = str(self.format), state = -1)



    def generate_receptor_pdbqt(self, strc):

        self.new_strc_name_format = str(self.new_strc_name + ".pdbqt")

        self.new_receptor_file_path = os.path.join(self.tmp_path, self.new_strc_name_format)

        # Log PDBQT
        self.pdbqt_log_file_name = str(self.new_strc_name + "_PDBQT_LOG.txt")

        self.receptors_settings = ["python", self.prepare_receptor_path, "-r", str(strc + ".pdb"), "-o", self.new_strc_name_format, "-v"]

        # If the User chose to add Hydrogens; default = None
        """
            [-A]  type(s) of repairs to make:
             'bonds_hydrogens': build bonds and add hydrogens
             'bonds': build a single bond from each atom with no bonds to its closest neighbor
             'hydrogens': add hydrogens
             'checkhydrogens': add hydrogens only if there are none already
             'None': do not make any repairs
             (default is 'None')
        """

        if self.pdbqt_dict["add_h"] and self.pdbqt_dict["bonds"]:
            self.receptors_settings.extend(["-A", "bonds_hydrogens"])

        elif self.pdbqt_dict["add_h"]:
            self.receptors_settings.extend(["-A", "hydrogens"])

        elif self.pdbqt_dict["bonds"]:
            self.receptors_settings.extend(["-A", "bonds"])

        # If the User chose to delete any non-standard residue from any chain
        """
            [-e]  delete every nonstd residue from any chain
              'True': any residue whose name is not in this list:
                      ['CYS','ILE','SER','VAL','GLN','LYS','ASN',
                      'PRO','THR','PHE','ALA','HIS','GLY','ASP',
                      'LEU', 'ARG', 'TRP', 'GLU', 'TYR','MET',
                      'HID', 'HSP', 'HIE', 'HIP', 'CYX', 'CSS']
              will be deleted from any chain.
              NB: there are no  nucleic acid residue names at all
              in the list and no metals.
             (default is False which means not to do this)
        """
        if self.pdbqt_dict["remove_nonstd"]:
            self.receptors_settings.extend(["-e", "True"])

        ### -C option: if the user do not chose to add Gasteiger charges
        """
        [-C]  preserve all input charges ie do not add new charges "
        (default is addition of gasteiger charges)"
        """

        if not self.pdbqt_dict["add_gast"]:
            self.receptors_settings.extend(["-C"])

        ### -U options
        """
            [-U]  cleanup type:
             'nphs': merge charges and remove non-polar hydrogens
             'lps': merge charges and remove lone pairs
             'waters': remove water residues
             'nonstdres': remove chains composed entirely of residues of
                      types other than the standard 20 amino acids
             'deleteAltB': remove XX@B atoms and rename XX@A atoms->XX
             (default is 'nphs_lps_waters_nonstdres')
        """
        U_options = []

        # If the User chose to remove Water molecules
        if self.pdbqt_dict["remove_water"]:
            U_options.append('waters')

        if self.pdbqt_dict["remove_lone_pairs"]:
            U_options.append('lps')

        if self.pdbqt_dict["remove_non_polar_H"]:
            U_options.append('nphs')

        if self.pdbqt_dict["remove_non_protein"]:
            U_options.append('nonstdres')

        U_options_string = '_'.join(U_options)

        if U_options:
            self.receptors_settings.extend(["-U", U_options_string])
        else:
            U_options_string = ' '
            self.receptors_settings.extend(["-U", U_options_string])

        if sys.platform == "win32":

            with open(self.pdbqt_log_file_name, "w") as pdbqt_file:
                try:
                    subprocess.run(self.receptors_settings, stdout=pdbqt_file)
                except subprocess.CalledProcessError as error:
                    print(error)

        else:

            with open(self.pdbqt_log_file_name, "w") as pdbqt_file:
                try:
                    subprocess.run(self.receptors_settings, stdout=pdbqt_file)
                except subprocess.CalledProcessError as error:
                    print(error)


    def generate_ligand_pdbqt(self, strc, output_name):

        self.new_strc_name_format = str(self.new_strc_name + ".pdbqt")

        self.new_receptor_file_path = os.path.join(self.tmp_path, self.new_strc_name_format)

        # Log PDBQT
        self.pdbqt_log_file_name = str(self.new_strc_name + "_PDBQT_LOG.txt")

        self.ligands_settings = ["python", self.prepare_ligand_path, "-l", str(strc + ".pdb"), "-o", self.new_strc_name_format, "-v"]

        if self.pdbqt_dict["add_h"]:
            self.ligands_settings.extend(["-A", "hydrogens"])
            cmd.h_add(self.strc)

        if self.pdbqt_dict["none_torsions"]:
            self.ligands_settings.extend(["-Z"])

        elif self.pdbqt_dict["all_torsions"]:
            self.ligands_settings.extend(["-B", "'amide'", "-B" "'guanidinium'"])

        if sys.platform == "win32":

            with open(self.pdbqt_log_file_name, "w") as pdbqt_file:
                try:
                    subprocess.run(self.ligands_settings, stdout=pdbqt_file)
                except subprocess.CalledProcessError as error:
                    print(error)

        else:

            with open(self.pdbqt_log_file_name, "w") as pdbqt_file:
                try:
                    subprocess.run(self.ligands_settings, stdout=pdbqt_file)
                except subprocess.CalledProcessError as error:
                    print(error)



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
