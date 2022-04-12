# Copyright 2022 by Serena Rosignoli. All rights reserved.
# This code is part of the DockingPie package and governed by its license. Please
# see the LICENSE file.

"""
Exceptions to be used within DockingPie.
"""

class DockingPieInvalidFile(Exception):
    """
    Used when a sequence or structure file containing some error is opened.
    """
    pass


class DockingPieUnknownFile(Exception):
    """
    Used when a sequence or structure file containing some error is opened.
    """
    pass


class DockingPieMissingStructure(Exception):
    """
    Used when trying to access the 3D structure data of an element which lacks a
    structure.
    """
    pass


class DockingPieSequenceConflict(Exception):
    """
    Used when updating the amino acid sequence of an element with a sequence with
    different amino acids.
    """
    pass


class DockingPieInterruptedProtocol(Exception):
    """
    Used when interrupting a running protocol.
    """
    pass


def catch_protocol_exception(function):
    """
    Function used as a decorator to catch the exceptions raised when running methods
    of protocols.
    """

    def wrapper(self, *args, **kwargs):

        # Launches the method of the protocol.
        try:
            return function(self, *args, **kwargs)

        # The protocol is quitted by the user while the method is running, do not show any meesage.
        except DockingPieInterruptedProtocol as e:
            self.quit_protocol()
            return None

        # The protocol stops because an exception was raised in the method, show an error message.
        except Exception as e:
            # Shows an error message.
            title = "%s Error" % self.protocol_name
            message = "%s stopped because of the following error: %s" % (self.protocol_name, str(e))
            self.DockingPie.main_window.show_error_message(title, message)
            self.quit_protocol()
            return None

    return wrapper
