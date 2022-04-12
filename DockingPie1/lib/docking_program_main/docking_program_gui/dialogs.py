# Copyright 2022 by Serena Rosignoli. All rights reserved.
# This code is part of the DockingPie package and governed by its license. Please
# see the LICENSE file.


import os

from pymol.Qt import QtWidgets, QtCore, QtGui


###############################################################################
# Dialogs.                                                                    #
###############################################################################

def askyesno_qt(title, message, parent=None, buttons_text=None):
    """
    Wrapper to a Yes/no dialog in PyQt. If 'buttons_text' is 'None', the default
    "Yes" and "No" buttons will be used. If 'buttons_text' is a list with two
    strings, the first string will be the text of the "Yes" button and the second
    one will be the text of the "No" button.
    """

    # Use yes and no buttons.
    if buttons_text is None:
        answer = QtWidgets.QMessageBox.question(parent, title, message,
                                                QtWidgets.QMessageBox.Yes,
                                                QtWidgets.QMessageBox.No)
        return answer == QtWidgets.QMessageBox.Yes

    # Set custom text on the buttons.
    else:
        dialog = QtWidgets.QMessageBox(parent)
        dialog.setWindowTitle(title)
        dialog.setText(message)
        yesbutton = dialog.addButton(buttons_text[0], QtWidgets.QMessageBox.YesRole)
        nobutton = dialog.addButton(buttons_text[1], QtWidgets.QMessageBox.NoRole)
        answer = dialog.exec_()
        return dialog.clickedButton() is yesbutton


def askopenfile_qt(title, parent=None, initialdir="", initialfile=None, name_filter=""):
    """
    Wrapper to a show a "pick a file to open" dialog in PyQt.
    """
    askfile_dialog = QtWidgets.QFileDialog()
    if initialdir and os.path.isdir(initialdir):
        _initialdir = initialdir
    else:
        _initialdir = ""

    askfile_dialog.setFileMode(QtWidgets.QFileDialog.ExistingFile)

    filepath = askfile_dialog.getOpenFileName(parent, title, _initialdir, name_filter)
    if isinstance(filepath, (tuple, list)):
        filepath = filepath[0]
    else:
        filepath = str(filepath)

    return filepath


def askopenfiles_qt(title, parent=None, initialdir="", name_filter=""):
    """
    Wrapper to a show a "pick multiple files to open" dialog in PyQt.
    """
    askfile_dialog = QtWidgets.QFileDialog()
    if initialdir and os.path.isdir(initialdir):
        _initialdir = initialdir
    else:
        _initialdir = ""

    askfile_dialog.setFileMode(QtWidgets.QFileDialog.ExistingFile)

    _filepaths = askfile_dialog.getOpenFileNames(parent, title, _initialdir, name_filter)
    if isinstance(_filepaths, (tuple, list)):
        filepaths = _filepaths[0]
    else:
        filepaths = _filepaths

    return filepaths


def askdirectory_qt(title, parent=None, initialdir=""):
    """
    Wrapper to a show a "pick a directory" dialog in PyQt.
    """
    askdirectory_dialog = QtWidgets.QFileDialog()
    if initialdir and os.path.isdir(initialdir):
        _initialdir = initialdir
    else:
        _initialdir = ""
    flags = QtWidgets.QFileDialog.ShowDirsOnly
    dirpath = str(askdirectory_dialog.getExistingDirectory(parent, title, _initialdir, flags))
    return dirpath


def asksaveasfile_qt(title, parent=None, initialdir="", name_filter="", check_existent=True):
    """
    Wrapper to a show a "pick a file to save" dialog in PyQt.
    """
    askfile_dialog = QtWidgets.QFileDialog()
    if initialdir and os.path.isdir(initialdir):
        _initialdir = initialdir
    else:
        _initialdir = ""

    _filepath = askfile_dialog.getSaveFileName(parent, title, _initialdir, name_filter)
    if isinstance(_filepath, (tuple, list)):
        filepath = _filepath[0]
        sel_filter = _filepath[1]
    else:
        filepath = str(_filepath)
        sel_filter = ""

    if not filepath:
        return filepath

    if name_filter and check_existent:
        if sel_filter:
            extension = sel_filter[1:]
        else:
            return filepath

        # PyQt has already checked for the existance of the file.
        if filepath.endswith(extension):
            return filepath
        # PyQt may not have seen the file with the extension.
        else:
            if os.path.isfile(filepath + extension):
                title = "Save As"
                message = "File '%s' already exists. Do you want to replace it?" % (filepath + extension)
                choice = askyesno_qt(title, message, parent=parent)
                if choice:
                    return filepath + extension
                else:
                    return ""
            else:
                return filepath + extension
    else:
        return filepath
