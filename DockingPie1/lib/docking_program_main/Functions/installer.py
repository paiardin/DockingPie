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
import time
import zipfile
from pathlib import Path
import io as io

# PyMOL.
import pymol
from pymol import cmd
from pymol.cgo import *
from pymol.vfont import plain
from pymol.Qt import QtWidgets, QtCore, QtGui
from pymol import stored
import math
import statistics

import urllib.request
#
from lib.docking_program_main.Functions.threads import Protocol_exec_dialog
from lib.docking_program_main.docking_program_gui.new_windows import NewWindow


def catch_errors_installer_threads(function):
    """
    Catches errors in the installer threads and closes the dialogs.
    """

    def wrapper(self, *args, **kwargs):
        try:
            return function(self, *args, **kwargs)
        except Exception as e:
            self.emit_exception.emit(e)

    return wrapper


class Installer_dialog_mixin:
    """
    Mixin class to be incorporated in all the installation process dialogs.
    """

    def keyPressEvent(self, event):
        """
        By overriding this method, the dialog will not close when pressing the "esc" key.
        """
        if event.key() == QtCore.Qt.Key_Escape:
            pass
        else:
            QtWidgets.QDialog.keyPressEvent(self, event)


class External_components_dialog(Installer_dialog_mixin, QtWidgets.QDialog):
    """
    A dialog to download from GitHub an archive with binary files for the external tools of PyMod.
    """

    def __init__(self, tab, url, os_arch, dir_name, local_mode=False):

        self.tab = tab


        QtWidgets.QDialog.__init__(self, parent=self.tab)

        download_path = os.path.join(self.tab.main_window.docking_programs.config_path, "external_tools_linux.zip")
        download_dir = self.tab.main_window.docking_programs.config_path
        self.local_mode = local_mode

        # Performs "local" installation.
        self.backup_download_path = None
        if self.local_mode:
            backup_download_path = askopenfile_qt("Select archive file", name_filter="*.zip", parent=self.get_qt_parent())
            if backup_download_path == "" or not os.path.isfile(backup_download_path):
                self.pymod.main_window.show_warning_message("Input Warning", "Invalid archive file.")
            else:
                self.backup_download_path = backup_download_path

        self.mb_label = "0"
        self.tot_label = "0"
        self.total_size = None
        self.complete_status = False
        self.initUI()


        # Initializes the download thread.
        self.dl_thread = External_tools_download_thread()
        self.dl_thread.update_progressbar.connect(self.on_update_progressbar)
        self.dl_thread.set_params(url, download_path, dir_name = dir_name)
        self.dl_thread.get_current_size.connect(self.on_get_current_size)
        self.dl_thread.get_total_size.connect(self.on_get_total_size)
        self.dl_thread.emit_exception.connect(self.on_emit_exception)
        self.dl_thread.signal_start_install.connect(self.start_install_thread)

        # Initializes the installation thread.
        self.inst_thread = External_tools_installation_thread()
        # self.inst_thread.update_progressbar.connect(self.on_update_progressbar)
        self.inst_thread.set_params(archive_filename=download_path,
                                    archive_dst_directory = download_dir,
                                    dir_name = dir_name,
                                    from_main_window = True)
        self.inst_thread.complete_installation.connect(self.on_complete_installation)
        self.inst_thread.emit_exception.connect(self.on_emit_exception)


    def initUI(self):

        self.setWindowTitle('Download External Tools')

        vertical_layout = QtWidgets.QVBoxLayout()
        self.download_progressbar = QtWidgets.QProgressBar(self)
        # self.progress.setGeometry(0, 0, 340, 25)
        self.download_progressbar.setMaximum(100)
        if not self.local_mode:
            progressbar_text = "Click the button to start components download"
        else:
            if self.backup_download_path is not None:
                progressbar_text = "Click the button to start components installation"
            else:
                progressbar_text = "Invalid archive file"

        self.download_progressbar.setFormat("")
        self.download_progressbar.setValue(0)
        vertical_layout.addWidget(self.download_progressbar)
        self.progressbar_label = QtWidgets.QLabel(progressbar_text)
        vertical_layout.addWidget(self.progressbar_label)

        # Button for starting the installation.
        horizontal_layout = QtWidgets.QHBoxLayout()
        if not self.local_mode:
            start_button_text = 'Start Download'
        else:
            if self.backup_download_path is not None:
                start_button_text = 'Start Installation'
            else:
                start_button_text = 'Exit Installation'
        self.start_button = QtWidgets.QPushButton(start_button_text, self)
        # self.start_button.setStyleSheet(label_style_2)
        self.start_button.clicked.connect(self.on_button_click)
        horizontal_layout.addWidget(self.start_button)

        horizontal_layout.addStretch(1)

        # Button for canceling the installation.
        self.cancel_button = QtWidgets.QPushButton('Cancel', self)
        # self.cancel_button.setStyleSheet(label_style_2)
        self.cancel_button.clicked.connect(self.on_cancel_button_click)
        horizontal_layout.addWidget(self.cancel_button)

        vertical_layout.addLayout(horizontal_layout)

        self.setLayout(vertical_layout)


    # Interactions with the buttons.
    def on_button_click(self):
        # The installation has not been launched yet.
        if not self.complete_status:
            # Remote installation.
            if not self.local_mode:
                self.dl_thread.start()
                self.download_progressbar.setFormat("Connecting...")
                self.progressbar_label.setText("")
                self.start_button.setEnabled(False)
            # Local installation.
            else:
                if self.backup_download_path is not None:
                    self.start_install_thread()
                    self.start_button.setEnabled(False)
                else:
                    self.close()
        # The installation has been already launched.
        else:
            self.close()

    def on_cancel_button_click(self):
        self._terminate_threads()
        self.close()

    def closeEvent(self, evnt):
        self._terminate_threads()

    def _terminate_threads(self):
        if self.dl_thread.isRunning():
            self.dl_thread.terminate()
        if self.inst_thread.isRunning():
            self.inst_thread.terminate()


    # Interactions with the download thread.
    def on_get_total_size(self, value):
        """
        Gets the total size of the archive to download.
        """
        self.total_size = value
        self.tot_label = self._convert_to_mb_string(value)
        if self.total_size == 0:
            pass

    def on_update_progressbar(self, value):
        """
        Updates the percentage in the progressbar.
        """
        self.download_progressbar.setValue(value)

    def on_get_current_size(self, value, state):
        """
        Updates the MB values in the progressbar.
        """
        if not state in ("update", "pass", "completed"):
            raise KeyError(state)

        if state in ("update", "completed"):
            self.mb_label = self._convert_to_mb_string(value)

        if self.total_size != 0:
            if state == "completed":
                self.download_progressbar.setFormat(
                    "Download Completed: %p% (" + self.mb_label + " of " + self.tot_label + " MB)")
            else:
                self.download_progressbar.setFormat(
                    "Download Progress: %p% (" + self.mb_label + " of " + self.tot_label + " MB)")
        else:
            if state == "completed":
                self.download_progressbar.setFormat("Download Completed: " + self.mb_label + " MB")
            else:
                self.download_progressbar.setFormat("Download Progress: " + self.mb_label + " MB of ???")

    def _convert_to_mb_string(self, value):
        if 0 <= value < 100000:
            round_val = 3
        elif 100000 <= value:
            round_val = 1
        else:
            ValueError(value)
        return str(round(value/1000000.0, round_val))


    # Interactions with the installation thread.
    def start_install_thread(self, signal=0):
        self.inst_thread.start()
        # self.download_progressbar.setRange(0, 0)
        # self.download_progressbar.reset()
        self.download_progressbar.setValue(0)
        self.download_progressbar.setFormat("Unzipping (please wait...)")
        self.start_button.setText("Finish Download")

    def on_complete_installation(self):
        self.download_progressbar.setValue(100)
        self.download_progressbar.setFormat("Download Completed Successfully")
        self.start_button.setEnabled(True)
        self.cancel_button.setEnabled(False)
        self.complete_status = True

    def on_emit_exception(self, e):
        """
        Quit the threads and close the installation dialog.
        """
        self._terminate_threads()

        message = "There was an error: %s" % str(e)
        if hasattr(e, "url"):
            message += " (%s)." % e.url
        else:
            message += "."
        message += " Quitting the installation process."
        print(message)
        ##self.tab.show_error_message("Installation Error", message)

        self.close()



class External_tools_download_thread(QtCore.QThread):
    """
    Runs a download thread to download PyMod components.
    """

    # Signals.
    update_progressbar = QtCore.pyqtSignal(int)
    get_current_size = QtCore.pyqtSignal(int, str)
    get_total_size = QtCore.pyqtSignal(int)
    emit_exception = QtCore.pyqtSignal(Exception)
    signal_start_install = QtCore.pyqtSignal(int)


    def set_params(self, url, download_path, dir_name):
        self.url = url
        self.download_path = download_path
        self.block_size = 8192
        self.write_mode = "wb"

    @catch_errors_installer_threads
    def run(self):

        # Connects with the web server.
        with urllib.request.urlopen(self.url) as http_resp:

            file_size = http_resp.headers["Content-Length"]

            if file_size is None:
                file_size = 0
            else:
                file_size = int(file_size)
            self.get_total_size.emit(file_size)

            with open(self.download_path, self.write_mode) as o_fh:

                file_size_dl = 0
                chunks_count = 0

                while True:

                    buffer_data = http_resp.read(self.block_size)
                    file_size_dl += len(buffer_data)

                    # Updates the fraction in the progressbar.
                    if file_size != 0:
                        # Can compute a fraction of the total size of the file to download.
                        frac = file_size_dl/file_size*100
                        self.update_progressbar.emit(int(frac))

                    # Updates the MB values in the progressbar.
                    if chunks_count % 100 == 0:
                        self.get_current_size.emit(file_size_dl, "update")
                    else:
                        self.get_current_size.emit(0, "pass")

                    if not buffer_data:
                        break

                    o_fh.write(buffer_data)
                    chunks_count += 1

                # Completes and updates the GUI.
                self.get_current_size.emit(file_size_dl, "completed")
                self.update_progressbar.emit(100)
                # Wait a little bit of time.
                time.sleep(0.5)
                self.signal_start_install.emit(0)



class External_tools_installation_thread(QtCore.QThread):
    """
    Runs a installation thread to install PyMod components.
    """

    # Signals.
    complete_installation = QtCore.pyqtSignal(int)
    emit_exception = QtCore.pyqtSignal(Exception)

    def set_params(self, archive_filename,
                   archive_dst_directory,
                   dir_name,
                   from_main_window=True):

        self.archive_filename = archive_filename
        self.archive_dst_directory = archive_dst_directory

        self.dir_name = dir_name

        self.from_main_window = from_main_window

    @catch_errors_installer_threads
    def run(self):
        """
        Unzips the archive containing the files of PyMod external tools.
        """

        #-----------------------------------------
        # Check if the file is a valid zip file. -
        #-----------------------------------------

        if not zipfile.is_zipfile(self.archive_filename):
            raise TypeError("The '%s' file is not a zip file." % self.archive_filename)

        #---------------------------------------------
        # Extract the file to a temporary directory. -
        #---------------------------------------------

        shutil.unpack_archive(self.archive_filename, self.archive_dst_directory, format="zip")

        if os.path.isdir(self.archive_dst_directory):
            os.remove(self.archive_filename)

        for file in ["prepare_flexreceptor4.py", "prepare_ligand4.py", "prepare_receptor4.py", "AutoDockTools"]:
            shutil.move(os.path.join(self.archive_dst_directory, self.dir_name, file), self.archive_dst_directory)

        #------------------------------------------------------------------
        # Updates the GUI so that the user can complete the installation. -
        #------------------------------------------------------------------

        if self.from_main_window:
            self.complete_installation.emit(0)




class Installation():

    def __init__(self, tab, program_to_install):

        self.tab = tab
        self.program_to_install = program_to_install

        # Check if Conda is integrated within PyMOL (Open-source PyMOL does not support Conda)

        try:
            import conda
            self.has_conda = True
        except ImportError:
            self.has_conda = False


        self.installation_window()


    def installation_window(self):

        qm = QtWidgets.QMessageBox.question(self.tab,'Installation Dialog',
        str("Do you want to proceed with " + self.program_to_install + " installation?"), QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)

        if qm == QtWidgets.QMessageBox.Yes:

            self.start_installation()

        else:
            pass


    def start_installation(self):

        self.installation_completed = True

        if self.program_to_install == "Openbabel":

            p_dialog = Protocol_exec_dialog(app=self.tab, docking_pie=self.tab,
                                            function=self.start_openbabel_installation,
                                            args=(),
                                            wait_start=0.4, wait_end=0.4,
                                            lock=True,
                                            stdout_silence=False,
                                            title="Installation",
                                            label_text="Installing " + self.program_to_install + ". Please wait.")
            p_dialog.exec_()


        elif self.program_to_install == "RxDock":

            p_dialog = Protocol_exec_dialog(app=self.tab, docking_pie=self.tab,
                                            function=self.start_rxdock_installation,
                                            args=(),
                                            wait_start=0.4, wait_end=0.4,
                                            lock=True,
                                            stdout_silence=False,
                                            title="Installation",
                                            label_text="Installing " + self.program_to_install + ". Please wait.")
            p_dialog.exec_()


        elif self.program_to_install == "Vina":

            p_dialog = Protocol_exec_dialog(app=self.tab, docking_pie=self.tab,
                                            function=self.start_vina_installation,
                                            args=(),
                                            wait_start=0.4, wait_end=0.4,
                                            lock=True,
                                            stdout_silence=False,
                                            title="Installation",
                                            label_text="Installing " + self.program_to_install + ". Please wait.")
            p_dialog.exec_()


        elif self.program_to_install == "Smina":

            p_dialog = Protocol_exec_dialog(app=self.tab, docking_pie=self.tab,
                                            function=self.start_smina_installation,
                                            args=(),
                                            wait_start=0.4, wait_end=0.4,
                                            lock=True,
                                            stdout_silence=False,
                                            title="Installation",
                                            label_text="Installing " + self.program_to_install + ". Please wait.")
            p_dialog.exec_()


        elif self.program_to_install == "sPyRMSD":

            p_dialog = Protocol_exec_dialog(app=self.tab, docking_pie=self.tab,
                                            function=self.start_spyrmsd_installation,
                                            args=(),
                                            wait_start=0.4, wait_end=0.4,
                                            lock=True,
                                            stdout_silence=False,
                                            title="Installation",
                                            label_text="Installing " + self.program_to_install + ". Please wait.")
            p_dialog.exec_()

        elif self.program_to_install == "ADFR":

            p_dialog = Protocol_exec_dialog(app=self.tab, docking_pie=self.tab,
                                            function=self.start_adfr_installation,
                                            args=(),
                                            wait_start=0.4, wait_end=0.4,
                                            lock=True,
                                            stdout_silence=False,
                                            title="Installation",
                                            label_text="Installing " + self.program_to_install + ". Please wait.")
            p_dialog.exec_()

        elif self.program_to_install == "sdsorter":

            p_dialog = Protocol_exec_dialog(app=self.tab, docking_pie=self.tab,
                                            function=self.start_sdsorter_installation,
                                            args=(),
                                            wait_start=0.4, wait_end=0.4,
                                            lock=True,
                                            stdout_silence=False,
                                            title="Installation",
                                            label_text="Installing " + self.program_to_install + ". Please wait.")
            p_dialog.exec_()


        if self.installation_completed == False and self.conda_installation == True and self.has_conda:

            self.about_text_area = QtWidgets.QPlainTextEdit()
            self.about_text_area.setReadOnly(True)

            additional_text = " "

            if re.search('EnvironmentNotWritableError', str(self.s)):
                additional_text = "\n!! Please Read Below !!\n \nPossible solutions:\n- Run PyMOL with administrator privileges\n- Install PyMOL without administrator priviliges"

            self.about_text_area.setPlainText(str(self.s) + additional_text)

            self.detailed_info_window = NewWindow(parent = self.tab,
            title = "Warning", upper_frame_title = "Conda Install Error",
            submit_command = None, with_scroll = True)

            self.detailed_info_window.middle_layout_type.addWidget(self.about_text_area, 0, 0)

            self.detailed_info_window.show()

        elif self.installation_completed == False and self.has_conda:

            self.about_text_area = QtWidgets.QPlainTextEdit()
            self.about_text_area.setReadOnly(True)

            self.about_text_area.setPlainText(str(self.s))

            self.detailed_info_window = NewWindow(parent = self.tab,
            title = "Warning", upper_frame_title = "Install Error",
            submit_command = None, with_scroll = True)

            self.detailed_info_window.middle_layout_type.addWidget(self.about_text_area, 0, 0)

            self.detailed_info_window.show()

        elif self.installation_completed == False and self.has_conda == False:
            QtWidgets.QMessageBox.warning(self.tab, "", str("Conda is not integrated in this PyMOL release. \n" + self.program_to_install + " can not be automatically installed"))

        else:
            if self.program_to_install == "ADFR":
                QtWidgets.QMessageBox.about(self.tab, "DockingPie", str(self.program_to_install) + " installation completed\n Please Restart PyMOL")
            else:
                QtWidgets.QMessageBox.about(self.tab, "DockingPie", str(self.program_to_install) + " installation completed")

        self.tab.check_installation()


    def start_sdsorter_installation(self):

        self.conda_installation = False

        if sys.platform == "linux":
            path_to_sdsorter = (os.path.join(self.tab.docking_programs.config_path, "external_tools_linux", "sdsorter_linux", "bin"))
            name = "sdsorter.static"

        if sys.platform == "darwin":
            path_to_sdsorter = (os.path.join(self.tab.docking_programs.config_path, "external_tools_macOS", "sdsorter_darwin", "bin"))
            name = "sdsorter.osx"

        os.chdir(path_to_sdsorter)
        subprocess.run(["chmod", "755", name])

        path_to_sdsorter = os.path.join(self.tab.docking_programs.path_to_sdsorter)
        try:
            result = subprocess.run([path_to_sdsorter], stdout=subprocess.PIPE, stderr = subprocess.PIPE)
        except PermissionError as e:
            self.s = str(e)
            self.installation_completed = False


    def start_spyrmsd_installation(self):

        self.rosconda("install -c conda-forge spyrmsd")

        # try:
        #
        #     pymol.externing.conda("install -c conda-forge spyrmsd")
        #
        #     self.installation_completed = True
        #
        # except Exception as e:
        #     print(e)
        #     self.installation_completed = False


    def start_openbabel_installation(self):

        self.rosconda("install -c conda-forge openbabel")

        # try:
        #
        #     pymol.externing.conda("install -c conda-forge openbabel")
        #     self.installation_completed = True
        #
        # except Exception as e:
        #     print(e)
        #     self.installation_completed = False


    def start_rxdock_installation(self):

        self.rosconda("install -c bioconda rxdock")

        # try:
        #     pymol.externing.conda("install -c bioconda rxdock")
        #     self.installation_completed = True
        #
        # except Exception as e:
        #     print(e)
        #     self.installation_completed = False
        #

    def start_adfr_installation(self):

        self.conda_installation = False

        if sys.platform == "win32":
            os.chdir(self.tab.docking_programs.path_to_ADFR)
            returned_value = subprocess.call("ADFRsuite_win32_1.0_Setup.exe", shell=True)

        else:
            if sys.platform == "linux":
                tmp_path_bin = (os.path.join(self.tab.docking_programs.config_path, "external_tools_linux", "ADFRsuite_x86_64Linux_1.0", "bin"))
                tmp_path_adfr = (os.path.join(self.tab.docking_programs.config_path, "external_tools_linux", "ADFRsuite_x86_64Linux_1.0"))
                archos = (os.path.join(self.tab.docking_programs.config_path, "external_tools_linux", "ADFRsuite_x86_64Linux_1.0", "bin", "archosv"))

            if sys.platform == "darwin":
                tmp_path_bin = (os.path.join(self.tab.docking_programs.config_path, "external_tools_macOS", "ADFRsuite_x86_64Darwin_1.0", "bin"))
                tmp_path_adfr = (os.path.join(self.tab.docking_programs.config_path, "external_tools_macOS", "ADFRsuite_x86_64Darwin_1.0"))
                archos = (os.path.join(self.tab.docking_programs.config_path, "external_tools_macOS", "ADFRsuite_x86_64Darwin_1.0", "bin", "archosv"))

            os.chdir(tmp_path_adfr)

            subprocess.run(["chmod", "755", "install.sh"], check = True)
            subprocess.run(["./install.sh", "-d", tmp_path_adfr], check =  True)

            subprocess.run(["chmod", "755", archos], check = True)

            path_to_adfr = os.path.join(self.tab.docking_programs.path_to_ADFR)

            try:
                result = subprocess.run([path_to_adfr], stdout=subprocess.PIPE, stderr = subprocess.PIPE)
            except PermissionError as e:
                self.s = str(e)
                self.installation_completed = False


    def start_vina_installation(self):
        self.conda_installation = False

        if sys.platform == "linux":
            path_to_vina = (os.path.join(self.tab.docking_programs.config_path, "external_tools_linux", "vina_linux", "bin"))
            name = "vina"

        if sys.platform == "darwin":
            path_to_vina = (os.path.join(self.tab.docking_programs.config_path, "external_tools_macOS", "vina_darwin", "bin"))
            name = "vina"

        if sys.platform == "win32":
            path_to_vina = (os.path.join(self.tab.docking_programs.config_path, "external_tools_windows", "vina_win32", "bin"))
            name = "vina.exe"

        os.chdir(path_to_vina)
        subprocess.run(["chmod", "755", name])

        path_to_vina = os.path.join(self.tab.docking_programs.path_to_vina)

        try:
            result = subprocess.run([path_to_vina], stdout=subprocess.PIPE, stderr = subprocess.PIPE)
        except PermissionError as e:
            self.s = str(e)
            self.installation_completed = False

        # try:
        #     result = subprocess.run([path_to_vina], stdout=subprocess.PIPE, stderr = subprocess.PIPE)
        #     self.installation_completed = True
        #
        # except Exception as e:
        #     self.installation_completed = False
        #

    def start_smina_installation(self):
        self.conda_installation = False

        if sys.platform == "linux":
            path_to_smina = (os.path.join(self.tab.docking_programs.config_path, "external_tools_linux", "smina_linux", "bin"))
            name = "smina.static"

        if sys.platform == "darwin":
            path_to_smina = (os.path.join(self.tab.docking_programs.config_path, "external_tools_macOS", "smina_darwin", "bin"))
            name = "smina.osx"

        os.chdir(path_to_smina)
        subprocess.run(["chmod", "755", name])

        path_to_smina = os.path.join(self.tab.docking_programs.path_to_smina)
        try:
            result = subprocess.run([path_to_smina], stdout=subprocess.PIPE, stderr = subprocess.PIPE)
        except PermissionError as e:
            self.s = str(e)
            self.installation_completed = False

        #subprocess.run([path_to_smina], stdout=subprocess.PIPE, stderr = subprocess.PIPE)

        # try:
        #     subprocess.run([path_to_smina], stdout=subprocess.PIPE, stderr = subprocess.PIPE)
        #     self.installation_completed = True
        #
        # except Exception as e:
        #     self.installation_completed = False
        #     print(self.installation_completed)


    def rosconda(self, command):


        self.conda_installation = True

        if self.has_conda:

            try:
                import conda.cli
            except ImportError:
                raise pymol.CmdException('conda wrapper only functional in PyMOL bundles')

            import sys, shlex

            args = shlex.split(command)

            if args:
                if args[0] in ('install', 'update'):
                    args[1:1] = ['--yes', '--prefix', sys.prefix]

            from contextlib import redirect_stderr
            #with open('conda_filename.log', 'w') as stderr, redirect_stderr(stderr):
            versione = cmd.get_version()

            if str(versione[0]) == "2.5.3":
                with io.StringIO() as stderr, redirect_stderr(stderr):
                    r = conda.cli.main(*args)
                    self.s = stderr.getvalue()
                if not r:
                    print(' conda finished with success')
                else:
                    self.installation_completed = False

            else:
                with io.StringIO() as stderr, redirect_stderr(stderr):
                    r = conda.cli.main('conda', *args)
                    self.s = stderr.getvalue()
                if not r:
                    print(' conda finished with success')
                else:
                    self.installation_completed = False

        else:
            self.installation_completed = False
