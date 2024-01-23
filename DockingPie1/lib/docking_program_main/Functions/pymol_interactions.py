# Copyright 2022 by Serena Rosignoli. All rights reserved.
# This code is part of the DockingPie package and governed by its license. Please
# see the LICENSE file.


import os
import sys
import shutil
import re
import json
import datetime
import warnings
import math
import subprocess
import statistics
import time

# PyMOL.
from pymol import cmd
from pymol import cmd
from pymol.cgo import *
from pymol.vfont import plain
from pymol import stored
from pymol.Qt import QtWidgets, QtCore, QtGui

# Functionality modules
from lib.docking_program_main.docking_program_gui.new_windows import NewWindow, Import_from_pymol_window_qt
from lib.docking_program_main.Functions.general_docking_func import ObjectParser

# Biopython modules
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.PDBIO import Select


class HideEverythingPyMOL():

    def __init__(self,
    tab,
    to_show = []):

        self.list_of_pymol_representations = ["lines", "spheres", "mesh", "ribbon", "cartoon", "sticks", "dots", "surface", "labels", "extent", "nonbonded", "nb_spheres", "slice", "extent", "slice", "dashes", "angles", "dihedrals", "cgo", "cell", "callback"]
        #self.list_of_pymol_representations = ["lines", "spheres", "ribbon", "cartoon", "sticks", "dots", "surface", "labels", "extent", "nonbonded", "nb_spheres", "slice", "extent", "slice", "dashes", "angles", "dihedrals", "callback"]
        self.list_to_show = to_show

        self.hide_everything_pymol()

    def hide_everything_pymol(self):

        for representation in self.list_of_pymol_representations:
            cmd.hide(representation, "all")

        for objs in self.list_to_show:
            cmd.show("cartoon", objs)
            cmd.show("sticks", objs)


class PyMOL_v25_bug():

    def __init__(self, tab,
    object_not_found):

        self.check_object_not_found_to_load(object_not_found)


    def check_object_not_found_to_load(self, object_not_found):

        type = cmd.get_type(object_not_found)
        if type == str("object:molecule"):
            cmd.delete(object_not_found)
            path_to_obj = os.path.join(self.tab.tab.parent().docking_programs.tmp_dir_path, str(object_not_found + ".pdb"))
            cmd.load(path_to_obj, object_not_found)


class PyMOL_Zoom_Orient_Show_Hide(PyMOL_v25_bug):

    def __init__(self, tab,
    obj,
    show_only = True):

        self.tab = tab
        self.show_only = show_only

        try:

            self.pymol_zoom_orient_show_hide(obj)

        except Exception as e:
                print("PyMOL v 2.5 bug")

                self.check_object_not_found_to_load(obj)

                self.pymol_zoom_orient_show_hide(obj)


    def pymol_zoom_orient_show_hide(self, obj):

        type = PyMOL_get_type(self, obj)

        if self.show_only:
            HideEverythingPyMOL(self)
            cmd.show("sticks", obj)
            cmd.orient(obj)

        else:
            if type.type == "object:map": #map object does not support orient function in PyMOL
                cmd.zoom(obj)
            else:
                cmd.orient(obj)


class PyMOL_Save_Action(PyMOL_v25_bug):

    def __init__(self, tab,
    obj,
    format = None,
    path = None,
    state = 0):

        self.tab = tab
        self.format = format
        self.obj = obj
        self.state = state
        self.path = path

        try:
            self.pymol_save_action()

        except Exception as e:
                print("PyMOL v 2.5 bug")
                self.check_object_not_found_to_load(obj)
                self.pymol_save_action()


    def pymol_save_action(self):

        cmd.save(self.path, self.obj, format = self.format, state = self.state)



class PyMOL_count_atoms(PyMOL_v25_bug):

    def __init__(self, tab,
    obj,
    format = None,
    path = None,
    state = 1):

        self.tab = tab

        try:
            self.atom_number = cmd.count_atoms(obj, state = state)

        except Exception as e:
                print("PyMOL v 2.5 bug")
                self.check_object_not_found_to_load(obj)
                self.atom_number = cmd.count_atoms(obj, state = state)



class PyMOL_get_type(PyMOL_v25_bug):

    def __init__(self, tab,
    obj,
    format = None,
    path = None,
    state = 1):

        self.tab = tab

        try:
            self.type = cmd.get_type(obj)

        except Exception as e:
                print("PyMOL v 2.5 bug")

                self.check_object_not_found_to_load(obj)
                self.type = cmd.get_type(obj)


class Import_from_Pymol():

    def __init__(self, tab,
    current_tab = None,
    is_receptor = False):

        self.tab = tab
        self.current_tab = current_tab
        self.is_receptor = is_receptor

        self.list_of_pymol_objects = []
        self.scrolledlist_items = []

        # Create a first list of all the objects and selections loaded in PyMOL
        self.list_pymol_objects(self.list_of_pymol_objects)

        if self.list_of_pymol_objects == []:
            QtWidgets.QMessageBox.warning(self.tab, "PyMOL Warning", str("There isn't any object to import"))
        else:
            # If the user is working from the Receptor tab
            if self.is_receptor:

                # # Check for multiple states
                # for i in self.list_of_pymol_objects:
                #     # check for multiple states of the object
                #     self.check_for_multiple_states(i)
                #     if self.has_multiple_state:
                #         cmd.split_states(i)
                #         cmd.delete(i)

                self.list_of_pymol_objects = []

                # If multiple states were present and splitted, the names of the objects have changed, hence a new list of PyMOL objects is needed
                self.list_pymol_objects(self.list_of_pymol_objects)

                for i in self.list_of_pymol_objects:
                    # Get PyMOL object type
                    type = cmd.get_type(i)
                    # If the PyMOL object is a molecule
                    if type == str("object:molecule"):
                        # count atoms
                        atom_number = cmd.count_atoms(i, state = 1)
                        if atom_number >= 150:
                            # check for multiple staes of the object
                            self.scrolledlist_items.append(i)


            else:
                # # Check for multiple states
                # for i in self.list_of_pymol_objects:
                #     # check for multiple states of the object
                #     self.check_for_multiple_states(i)
                #     if self.has_multiple_state:
                #         cmd.split_states(i)
                #         cmd.delete(i)

                self.list_of_pymol_objects = []

                # If multiple states were present and splitted, the names of the objects have changed, hence a new list of PyMOL objects is needed
                self.list_pymol_objects(self.list_of_pymol_objects)

                # If multiple states were present and splitted, the names of the objects have changed, hence a new list of PyMOL objects is needed
                for i in self.list_of_pymol_objects:
                    # Get PyMOL object type
                    type = cmd.get_type(i)
                    # If the PyMOL object is a molecule
                    if type == str("object:molecule"):
                        # count atoms
                        atom_number = cmd.count_atoms(i, state = 1)
                        if atom_number < 150:
                            # check for multiple staes of the object
                            self.scrolledlist_items.append(i)

            # New list of PyMOl objects, eventually with the splitted chains
            #self.tmp_list = [str(obj) for obj in cmd.get_names("objects") + cmd.get_names("selections")]

            if self.scrolledlist_items == []:
                QtWidgets.QMessageBox.warning(self.tab, "PyMOL Warning", str("There isn't any object to import"))
            else:
                # Builds a new window with the importable objects, if present
                self.import_from_pymol_window = Import_from_pymol_window_qt(self.tab,
                    title="Import from PyMOL",
                    upper_frame_title="List of Importable Objects",
                    submit_command=self.import_selected_pymol_object,
                    selections_list=self.scrolledlist_items)
                self.import_from_pymol_window.show()


    def list_pymol_objects(self, list):

        for obj in cmd.get_names("objects") + cmd.get_names("selections"):
            if re.search("Run_", str(obj)):
                pass
            else:
                list.append(str(obj))

    def import_selected_pymol_object(self): # When Submit button of the importable objects dialog is pressed

        self.is_vina_tab = False
        self.is_rxdock_tab = False
        self.is_smina_tab = False
        self.is_adfr_tab = False

        if self.current_tab == "ADFR":
            self.is_adfr_tab = True

        if self.current_tab == "Smina":
            self.is_smina_tab = True

        if self.current_tab == "Vina":
            self.is_vina_tab = True

        if self.current_tab == "RxDock":
            self.is_rxdock_tab = True


        if self.is_receptor: # if imported from Receptor tab
            self.selections_to_import = self.import_from_pymol_window.get_objects_to_import() # names list of checked options
            self.import_selected_pymol_receptors()

        else: # if imported from Ligand tab
            self.selections_to_import = self.import_from_pymol_window.get_objects_to_import() # names list of checked options
            self.import_selected_pymol_ligands()

        self.import_from_pymol_window.destroy()


    def import_selected_pymol_receptors(self):

        for pymol_obj in self.selections_to_import:

            tmp_path_name = os.path.join(self.tab.parent().docking_programs.tmp_dir_path, str(pymol_obj + ".pdb"))

            # This save-delete-load-copy-sleep-setname cycle fa veramente cagare i cani  is done to avoid errors while parsing of the files. In this way, despite the type of file that was previously loaded in PyMOL, now it is loaded as a PDB file.
            cmd.save(tmp_path_name, pymol_obj, format = 'pdb')
            cmd.delete(pymol_obj)
            cmd.load(tmp_path_name, pymol_obj)
            # cmd.copy("temp_obj", pymol_obj)
            # cmd.delete(pymol_obj)
            # cmd.set_name("temp_obj", pymol_obj)

            if self.is_vina_tab:
                already_loaded = Load_Object.check_already_loaded_structure(self = self, pymol_obj = pymol_obj, objects_dict = self.tab.parent().docking_programs.vina_receptors_dict)

                if already_loaded:
                    QtWidgets.QMessageBox.warning(self.tab, "Already Loaded File", str("The file '" + pymol_obj + "' is already loaded"))

                else:
                    os.chdir(self.tab.parent().docking_programs.vina_tmp_dir)
                    Load_Object.load_checked_structures(self = self,
                    pymol_obj = pymol_obj,
                    format = "pdb",
                    objects_dict = self.tab.parent().docking_programs.vina_receptors_dict,
                    is_receptor = True)

            elif self.is_rxdock_tab:
                already_loaded = Load_Object.check_already_loaded_structure(self = self, pymol_obj = pymol_obj, objects_dict = self.tab.parent().docking_programs.rxdock_receptors_dict)

                if already_loaded:
                    QtWidgets.QMessageBox.warning(self.tab, "Already Loaded File", str("The file '" + pymol_obj + "' is already loaded"))
                else:
                    os.chdir(self.tab.parent().docking_programs.rxdock_tmp_dir)
                    Load_Object.load_checked_structures(self = self,
                    pymol_obj = pymol_obj,
                    format = "pdb",
                    objects_dict = self.tab.parent().docking_programs.rxdock_receptors_dict,
                    is_receptor = True)

            elif self.is_smina_tab:

                already_loaded = Load_Object.check_already_loaded_structure(self = self, pymol_obj = pymol_obj, objects_dict = self.tab.parent().docking_programs.smina_receptors_dict)

                if already_loaded:
                    QtWidgets.QMessageBox.warning(self.tab, "Already Loaded File", str("The file '" + pymol_obj + "' is already loaded"))
                else:
                    os.chdir(self.tab.parent().docking_programs.smina_tmp_dir)
                    Load_Object.load_checked_structures(self = self,
                    pymol_obj = pymol_obj,
                    format = "pdb",
                    objects_dict = self.tab.parent().docking_programs.smina_receptors_dict,
                    is_receptor = True)

            elif self.is_adfr_tab:

                already_loaded = Load_Object.check_already_loaded_structure(self = self, pymol_obj = pymol_obj, objects_dict = self.tab.parent().docking_programs.adfr_receptors_dict)

                if already_loaded:
                    QtWidgets.QMessageBox.warning(self.tab, "Already Loaded File", str("The file '" + pymol_obj + "' is already loaded"))
                else:
                    os.chdir(self.tab.parent().docking_programs.adfr_tmp_dir)
                    Load_Object.load_checked_structures(self = self,
                    pymol_obj = pymol_obj,
                    format = "pdb",
                    objects_dict = self.tab.parent().docking_programs.adfr_receptors_dict,
                    is_receptor = True)



    def import_selected_pymol_ligands(self):

        for pymol_lig in self.selections_to_import:

            tmp_path_name = os.path.join(self.tab.parent().docking_programs.tmp_dir_path, str(pymol_lig + ".pdb"))

            cmd.save(tmp_path_name, pymol_lig, state = 0, format = 'pdb')
            cmd.delete(pymol_lig)
            cmd.load(tmp_path_name, pymol_lig)
            cmd.copy("temp_obj", pymol_lig)
            cmd.delete(pymol_lig)
            cmd.set_name("temp_obj", pymol_lig)


            if self.is_vina_tab:
                already_loaded = Load_Object.check_already_loaded_structure(self = self, pymol_obj = pymol_lig, objects_dict = self.tab.parent().docking_programs.vina_ligands_dict)

                if already_loaded:
                    QtWidgets.QMessageBox.warning(self.tab, "Already Loaded File", str("The file '" + pymol_lig + "' is already loaded"))
                else:
                    Load_Object.load_checked_structures(self = self,
                    pymol_obj = pymol_lig,
                    format = "pdb",
                    objects_dict = self.tab.parent().docking_programs.vina_ligands_dict)


            if self.is_rxdock_tab:
                already_loaded = Load_Object.check_already_loaded_structure(self = self, pymol_obj = pymol_lig, objects_dict = self.tab.parent().docking_programs.rxdock_ligands_dict)

                if already_loaded:
                    QtWidgets.QMessageBox.warning(self.tab, "Already Loaded File", str("The file '" + pymol_lig + "' is already loaded"))
                else:
                    Load_Object.load_checked_structures(self = self,
                    pymol_obj = pymol_lig,
                    format = "pdb",
                    objects_dict = self.tab.parent().docking_programs.rxdock_ligands_dict)


            if self.is_smina_tab:
                already_loaded = Load_Object.check_already_loaded_structure(self = self, pymol_obj = pymol_lig, objects_dict = self.tab.parent().docking_programs.smina_ligands_dict)

                if already_loaded:
                    QtWidgets.QMessageBox.warning(self.tab, "Already Loaded File", str("The file '" + pymol_lig + "' is already loaded"))
                else:
                    Load_Object.load_checked_structures(self = self,
                    pymol_obj = pymol_lig,
                    format = "pdb",
                    objects_dict = self.tab.parent().docking_programs.smina_ligands_dict)


            if self.is_adfr_tab:

                already_loaded = Load_Object.check_already_loaded_structure(self = self, pymol_obj = pymol_lig, objects_dict = self.tab.parent().docking_programs.adfr_ligands_dict)

                if already_loaded:
                    QtWidgets.QMessageBox.warning(self.tab, "Already Loaded File", str("The file '" + pymol_lig + "' is already loaded"))
                else:
                    os.chdir(self.tab.parent().docking_programs.adfr_tmp_dir)
                    Load_Object.load_checked_structures(self = self,
                    pymol_obj = pymol_lig,
                    format = "pdb",
                    objects_dict = self.tab.parent().docking_programs.adfr_ligands_dict)


    # def split_states_dialog(self, pymol_obj):
    #
    #     qm = QtWidgets.QMessageBox.question(self.tab,'', str("'" + pymol_obj + "' has multiple states. \nDo you want to split them into single objects?"), QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
    #
    #     if qm == QtWidgets.QMessageBox.Yes:
    #         cmd.split_states(pymol_obj)
    #         cmd.delete(pymol_obj)
    #
    #     else:
    #         pass
    #

    def check_for_multiple_states(self, pymol_obj):

        self.has_multiple_state = False

        states = cmd.count_states(pymol_obj)

        if states > 1:
            self.has_multiple_state = True



class PyMOLInteractions:


    """
    A class just to store some functions to interact with PyMOL
    """


    # A part of the following code has been adapted from the "grid settings functions" of Autodock/Vina plugin

    '''
    # Autodock/Vina plugin  Copyright Notice
    # ============================
    #
    # The Autodock/Vina plugin source code is copyrighted, but you can freely use and
    # copy it as long as you don't change or remove any of the copyright
    # notices.
    #
    # ----------------------------------------------------------------------
    # Autodock/Vina plugin is Copyright (C) 2009 by Daniel Seeliger
    #
    #                        All Rights Reserved
    #
    # Permission to use, copy, modify, distribute, and distribute modified
    # versions of this software and its documentation for any purpose and
    # without fee is hereby granted, provided that the above copyright
    # notice appear in all copies and that both the copyright notice and
    # this permission notice appear in supporting documentation, and that
    # the name of Daniel Seeliger not be used in advertising or publicity
    # pertaining to distribution of the software without specific, written
    # prior permission.
    #
    # DANIEL SEELIGER DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
    # SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
    # FITNESS.  IN NO EVENT SHALL DANIEL SEELIGER BE LIABLE FOR ANY
    # SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
    # RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
    # CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
    # CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
    # ----------------------------------------------------------------------
    '''

    def calculate_center(object):

        stored.xyz = []
        cmd.iterate_state(1, self.sel,"stored.xyz.append([x,y,z])")
        xx = statistics.mean(map(lambda a: a[0], stored.xyz))
        yy = statistics.mean(map(lambda a: a[1], stored.xyz))
        zz = statistics.mean(map(lambda a: a[2], stored.xyz))


    def calculate_box_changed(self, x, y, z, spacing, x_vis, y_vis, z_vis):

        x = float(x)
        y = float(y)
        z = float(z)

        xpts = int(x_vis)
        ypts = int(y_vis)
        zpts = int(z_vis)
        spacing = float(spacing)

        cylinder_size = float(0.2)

        size = [xpts*spacing, ypts*spacing, zpts*spacing]
        xmax = x + size[0]/2.
        xmin = x - size[0]/2.
        ymax = y + size[1]/2.
        ymin = y - size[1]/2.
        zmax = z + size[2]/2.
        zmin = z - size[2]/2.
        box_edge_x = [xmin,xmax]
        box_edge_y = [ymin,ymax]
        box_edge_z = [zmin,zmax]
        self.box_coords  = [box_edge_x,box_edge_y,box_edge_z]
        cmd.delete('box')
        self.display_box_changed(self.box_coords, cylinder_size)

        try:
            self.display_wire_box_changed(self.box_coords, spacing)
        except:
            pass


    def display_box_changed(self, box, cylinder_size):
        view = cmd.get_view()
        name = "box"
        obj = []
        # build cgo object
        color = [1.,1.,1.]
        for i in range(2):
            for k in range (2):
                for j in range(2):
                    if i != 1:
                        obj.append(CYLINDER)
                        obj.extend([box[0][i],box[1][j],box[2][k]])
                        obj.extend([box[0][i+1],box[1][j],box[2][k]])
                        obj.append(cylinder_size)
                        obj.extend(color)
                        obj.extend(color)
                        obj.append(COLOR)
                        obj.extend(color)
                        obj.append(SPHERE)
                        obj.extend([box[0][i],box[1][j],box[2][k],cylinder_size])

                    if j != 1:
                        obj.append(CYLINDER)
                        obj.extend([box[0][i],box[1][j],box[2][k]])
                        obj.extend([box[0][i],box[1][j+1],box[2][k]])
                        obj.append(cylinder_size)
                        obj.extend(color)
                        obj.extend(color)
                        obj.append(COLOR)
                        obj.extend(color)
                        obj.append(SPHERE)
                        obj.extend([box[0][i],box[1][j+1],box[2][k],cylinder_size])
                    if k != 1:
                        obj.append(CYLINDER)
                        obj.extend([box[0][i],box[1][j],box[2][k]])
                        obj.extend([box[0][i],box[1][j],box[2][k+1]])
                        obj.append(cylinder_size)
                        obj.extend(color)
                        obj.extend(color)
                        obj.append(COLOR)
                        obj.extend(color)
                        obj.append(SPHERE)
                        obj.extend([box[0][i],box[1][j],box[2][k+1],cylinder_size])
        axes = [[2.0,0.0,0.0],[0.0,2.0,0.0],[0.0,0.0,2.0]]
        xpos = [box[0][1]+(box[0][1]-box[0][0])/5.,box[1][0],box[2][0]]
        cyl_text(obj,plain,xpos,'X',0.10,axes=axes)
        ypos = [box[0][0],box[1][1]+(box[1][1]-box[1][0])/5,box[2][0]]
        cyl_text(obj,plain,ypos,'Y',0.10,axes=axes)
        zpos = [box[0][0],box[1][0],box[2][1]+(box[2][1]-box[2][0])/5]
        cyl_text(obj,plain,zpos,'Z',0.10,axes=axes)
        cmd.load_cgo(obj,name)
        cmd.set_view(view)
        cmd.zoom("box")


    def display_wire_box_changed(self, box, spacing):

        cmd.delete("wirebox")
        color = [1.0,1.0,1.0]
        view = cmd.get_view()
        spacing = float(spacing)
        lwidth = float(1.0)
        xpts = int(round((box[0][1]-box[0][0])/spacing))+1
        ypts = int(round((box[1][1]-box[1][0])/spacing))+1
        zpts = int(round((box[2][1]-box[2][0])/spacing))+1
        obj = []

        for i in range(xpts):
            for k in range(ypts):
                obj.append(BEGIN)
                obj.append(LINE_STRIP)
                obj.append(COLOR)
                obj.extend(color)

                for j in range(zpts):

                    obj.append(VERTEX)
                    obj.extend([box[0][0]+spacing*i,box[1][0]+spacing*k,\
                                box[2][0]+spacing*j])

                obj.append(END)
        for i in range(xpts):
            for j in range (zpts):
                obj.append(BEGIN)
                obj.append(LINE_STRIP)
                obj.append(COLOR)
                obj.extend(color)
                for k in range(ypts):
                    obj.append(VERTEX)
                    obj.extend([box[0][0]+spacing*i,box[1][0]+spacing*k,\
                                box[2][0]+spacing*j])
                obj.append(END)
        for j in range(zpts):
            for i in range (xpts):
                obj.append(BEGIN)
                obj.append(LINE_STRIP)
                obj.append(COLOR)
                obj.extend(color)
                for k in range(ypts):
                    obj.append(VERTEX)
                    obj.extend([box[0][0]+spacing*i,box[1][0]+spacing*k,\
                                box[2][0]+spacing*j])
                obj.append(END)
        for j in range(zpts):
            for k in range (ypts):
                obj.append(BEGIN)
                obj.append(LINE_STRIP)
                obj.append(COLOR)
                obj.extend(color)
                for i in range(xpts):
                    obj.append(VERTEX)
                    obj.extend([box[0][0]+spacing*i,box[1][0]+spacing*k,\
                                box[2][0]+spacing*j])
                obj.append(END)

        cmd.load_cgo(obj,"wirebox")
        cmd.set("cgo_line_width",lwidth)
        cmd.set_view(view)
        cmd.zoom("wirebox")


    def show_crisscross_changed(self, x, y, z):

        center = [float(x), float(y), float(z)]
        cmd.delete("grid_center")
        self.crisscross(center[0], center[1], center[2], 0.5, "grid_center")
        cmd.zoom("grid_center")


    def show_crisscross(self, key):

        center = [float(self.main_window.grid_center[key][0]),
                  float(self.main_window.grid_center[key][1]),
                  float(self.main_window.grid_center[key][2])
                  ]

        cmd.delete("grid_center")
        self.crisscross(center[0], center[1], center[2], 0.5, "grid_center")
        cmd.zoom("grid_center")


    def crisscross(self,x,y,z,d,name="crisscross"):

        obj = [
            LINEWIDTH, 3,

            BEGIN, LINE_STRIP,
            VERTEX, float(x-d), float(y), float(z),
            VERTEX, float(x+d), float(y), float(z),
            END,

            BEGIN, LINE_STRIP,
            VERTEX, float(x), float(y-d), float(z),
            VERTEX, float(x), float(y+d), float(z),
            END,

            BEGIN, LINE_STRIP,
            VERTEX, float(x), float(y), float(z-d),
            VERTEX, float(x), float(y), float(z+d),
            END

            ]
        view = cmd.get_view()
        cmd.load_cgo(obj,name)
        cmd.set_view(view)


    def show_box_func(self, object):

        stored.xyz = []
        cmd.iterate_state(1, self.sel,"stored.xyz.append([x,y,z])")
        self.xx = statistics.mean(map(lambda a: a[0], stored.xyz))
        self.yy = statistics.mean(map(lambda a: a[1], stored.xyz))
        self.zz = statistics.mean(map(lambda a: a[2], stored.xyz))

        return(self.xx, self.yy, self.zz)


    def load_element_in_pymol(self, file_to_load, pymol_object_name):

        cmd.load(file_to_load, pymol_object_name)
        cmd.util.cbac(pymol_object_name)
        cmd.zoom(pymol_object_name)
        cmd.center(pymol_object_name)



#Creating a new class: AtomSelect which calls the class Select
class AtomSelect(Select):
    #Using the method accept_atom available for the Select class
    def accept_atom(self, atom):
        #Selecting all not-disordered atoms and those disordered with the alternative positions labelled 'A'
        if (not atom.is_disordered()) or atom.get_altloc() == "A":
            # Eliminating alt location ID
            atom.set_altloc(" ")
            return True
        else:
            return False


class Load_Object:


    def load_checked_structures(self, pymol_obj, format, objects_dict, is_receptor = False):

            self.is_receptor = is_receptor

            cwd = os.getcwd()

            # create a file of the structure in tmp dir
            self.file_path = os.path.join(cwd, str(pymol_obj  + "." + format))
            self.file_name = pymol_obj
            cmd.save(self.file_path, pymol_obj, format = format)

            self.object = ObjectParser(file_name = self.file_name,
            file_path = self.file_path,
            is_receptor = self.is_receptor
            )

            parsed_file_handle = open(self.file_path, "r")
            self.parsed_biopython_structure = PDBParser(PERMISSIVE=1, QUIET=True).get_structure(self.file_name, parsed_file_handle)
            parsed_file_handle.close()
            io = PDBIO()
            io.set_structure(self.parsed_biopython_structure)
            io.save(self.file_path, select=AtomSelect())

            cmd.delete(pymol_obj)
            cmd.load(self.file_path, pymol_obj)

            # Create a nested dictionary for each loaded object with its information
            objects_dict[self.file_name] = {}
            objects_dict[self.file_name]["parsed_object"] = self.object

            # Build a frame for every single Object
            self.tab.build_structures_frame(self.file_name, objects_dict)

            self.tab.parent().docking_programs.main_window.statusBar().showMessage(str("Loaded file '" + pymol_obj + "'"), 3000)


    def check_already_loaded_structure(self, pymol_obj, objects_dict):

        self.is_already_loaded = False

        for i in objects_dict:
            if str(i) == str(pymol_obj):
                self.is_already_loaded = True

        return self.is_already_loaded
