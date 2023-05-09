# Copyright 2022 by Serena Rosignoli. All rights reserved.
# This code is part of the GenyDock package and governed by its license. Please
# see the LICENSE file.


import os
import shutil

from lib.docking_program_main.Functions.vina_functions import *
from lib.docking_program_main.Functions.rxdock_functions import RxDock_Functions
from lib.docking_program_main.docking_program_gui.dialogs import *
from pymol.Qt import QtWidgets, QtCore, QtGui
import csv


def write_log_titles(file, title, big = True):

    '''
    A function that writes titles inside ###
    n.b. the input must be an opened file, with the method file = 'open("", "w+")'
    '''

    file.write('################################################################')
    file.write('##########  {}  ######################'.format(title))

    if big:
        file.write('################################################################')

    file.write('                                                                    \n')

def check_configuration(self, docking_programs):

    if sys.platform == "linux":
        dir_name = "external_tools_linux"
    if sys.platform == "darwin":
        dir_name = "external_tools_macOS"
    if sys.platform == "win32":
        dir_name = "external_tools_windows"

    ext_tools_path = os.path.join(docking_programs.config_path, dir_name)

    if os.path.isdir(ext_tools_path):
        is_configured = True
    else:
        is_configured = False

    return is_configured


class SelectAll():


    def __init__(self, tab, all_cb, list_of_cb = []):

        self.tab = tab
        self.list_of_cb = list_of_cb

        self.get_state(all_cb)


    def get_state(self, all_cb):

        if all_cb.isChecked():
            self.all_check = True
        else:
            self.all_check = False

        self.all_func()


    def all_func(self):

        for checkbox in self.list_of_cb:
            if self.all_check:
                if checkbox.isChecked():
                    pass
                else:
                    checkbox.setChecked(True)
            else:
                if checkbox.isChecked():
                    checkbox.setChecked(False)
                else:
                    pass



class Check_current_tab():


    def check_docking_program_current_tab(self):

        a = self.parent().docking_programs.docking_programs_tabs.currentIndex()
        current_tab = self.parent().docking_programs.docking_programs_tabs.tabText(a)

        self.is_vina_tab = False
        self.is_rxdock_tab = False
        self.is_adfr_tab = False
        self.is_smina_tab = False

        if current_tab == "Vina":
            self.is_vina_tab = True

        if current_tab == "RxDock":
            self.is_rxdock_tab = True

        if current_tab == "Smina":
            self.is_smina_tab = True

        if current_tab == "ADFR":
            self.is_adfr_tab = True



class Save_to_Csv():

# to transform in general Save to File

    def __init__(self, tab, table):

        """
        Saves the table data to a .csv file.
        """

        self.tab =  tab
        self.table = table

        self.save_to_csv_event(table = self.table)


    def save_to_csv_event(self, table):

        # Let the user select the filepath.
        filepath = asksaveasfile_qt("Save CSV file", name_filter="*.csv")

        if not filepath:
            return None

        try:
            # Writes a .csv file on that path.
            with open(filepath, 'w') as csv_fh:

                writer = csv.writer(csv_fh, delimiter=',', quoting=csv.QUOTE_MINIMAL)

                if table.row_labels is not None:
                    writer.writerow([" "] + table.column_labels)
                else:
                    writer.writerow(table.column_labels)
                for row_idx, row in enumerate(table.data):
                    if table.row_labels is not None:
                        writer.writerow([table.row_labels[row_idx]] + [str(v) for v in row])
                    else:
                        writer.writerow([str(v) for v in row])

        except Exception as e:
            print("- WARNING: could not write a csv file: %s" % str(e))




class OpenFromFile():


    """
    Open from file class.
    This class represent the Opened Files
    """


    def __init__(self, tab, file_type_string):

        self.tab = tab
        self.valid_config = False

        # Get the selected file path
        self.file_path = QtWidgets.QFileDialog.getOpenFileName(self.tab, 'Open file',
            '\\', str(file_type_string))

        if self.file_path == "  ":
            QtWidgets.QFileDialog.close()

        else:
            extension = os.path.splitext(self.file_path[0])[1].replace(".","").lower()

            # Get the selected file name
            self.file_name = os.path.splitext(os.path.basename(self.file_path[0]))[0]
            self.file_name_ext = str(self.file_name + "." + extension)

            self.open_file(extension)


    def open_file(self, extension):

        self.is_valid = False
        self.is_already_loaded = False

        # Check the validity of the file and whether the file is already present in tmp_dirs
        if extension == "grd":
            RxDock_Functions.is_valid_grd_file(self, self.file_path[0])
            self.check_if_already_exist(tmp_dir = self.tab.docking_programs_child_tabs.docking_programs.rxdock_tmp_dir,
            extension = extension)

        elif extension == "as":
            RxDock_Functions.is_valid_as_file(self, self.file_path[0])
            self.check_if_already_exist(tmp_dir = self.tab.docking_programs_child_tabs.docking_programs.rxdock_tmp_dir,
            extension = extension)

        elif extension == "txt":
            self.is_valid_config_file(self.file_path[0])
            self.check_if_already_exist(tmp_dir = self.tab.docking_programs_child_tabs.docking_programs.vina_tmp_dir,
            extension = extension)

        elif extension == "const":
            RxDock_Functions.is_valid_constrains_file(self, self.file_path[0])
            self.check_if_already_exist(tmp_dir = self.tab.docking_programs_child_tabs.docking_programs.rxdock_tmp_dir,
            extension = extension)


    def check_if_already_exist(self, tmp_dir, extension):

        # Estabilish which would be the new path for the config File
        self.new_file_path = os.path.join(tmp_dir, str(self.file_name + "." + extension))

        if os.path.exists(self.new_file_path):
            qm = QtWidgets.QMessageBox.question(self.tab,'Warning', str("The file '" + self.file_name + "' already exist \nDo you want to replace it?"), QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

            if qm == QtWidgets.QMessageBox.Yes:
                try:
                    shutil.copy(self.file_path[0], tmp_dir)
                except:
                    pass

            elif qm == QtWidgets.QMessageBox.No:
                self.is_already_loaded = True

        else:
            shutil.copy(self.file_path[0], tmp_dir)



    def is_valid_config_file(self, file_name):

        self.valid_config = False

        # file_handler = open(file_name, "r")
        #
        # for line in file_handler.readlines():
        #     if line.startswith("center_x") or line.startswith("center_y") or line.startswith("center_z"):
        #             x,y,z = float(line[11:15]), float(line[11:15]), float(line[11:15])
        #             self.valid_config = True
        #
        # file_handler.close()

        self.config_parameters_list = ["center_x", "center_y", "center_z", "size_x", "size_z", "size_y", "exhaustiveness"]

        file_handler = open(file_name, "r")

        for line in file_handler:
            for parameter in self.config_parameters_list:
                if line.startswith(parameter):
                    value = re.findall('[-+]?([0-9]*\.[0-9]+|[0-9]+)', line)
                    if value:
                        pass
                    else:
                        self.config_parameters_list.remove(parameter)

        if self.config_parameters_list:

            self.valid_config = True

        file_handler.close()

        if not self.valid_config:
            QtWidgets.QMessageBox.warning(self.tab, "File Validity Error", str("The file '" + file_name + "' is not valid"))
            self.valid_config = False

        return self.valid_config
