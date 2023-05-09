
__author__    = "Andre Merzky, Ole Weidner"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


''' Provides a parser class for the file transfer specification as it is
    defined in GFD.90, sction 4.1.3.
'''

import saga.exceptions as se

# 4.1.3 File Transfer Specifications (GFD90 p 176-177)
#
# The syntax of a file transfer directive for the job description is modeled on 
# the LSF syntax (LSF stands for Load Sharing Facility, a commercial job 
# scheduler by Platform Computing), and has the general syntax:
# local_file operator remote_file
# Both the local_file and the remote_file can be URLs. If they are not URLs,
# but full or relative pathnames, then the local_file is relative to the host 
# where the submission is executed, and the remote_file is evaluated on the 
# execution host of the job. The operator is one of the following four:
#
# '>'  copies the local file to the remote file before the job starts. 
#      Overwrites the remote file if it exists.
# '>>' copies the local file to the remote file before the job starts. 
#      Appends to the remote file if it exists.
# '<'  copies the remote file to the local file after the job finishes. 
#      Overwrites the local file if it exists.
# '<<' copies the remote file to the local file after the job finishes. 
#      Appends to the local file if it exists.

class TransferDirectives(object):

    def __init__(self, directives_list):

        self._in_overwrite = dict()
        self._in_append = dict()
        self._out_overwrite = dict()
        self._out_append = dict()

        # each line in directives_list should contain one directive
        for directive in directives_list:
            if (directive.count('>') > 2) or (directive.count('<') > 2):
                msg = "'%s' is not a valid transfer directive string."
                raise se.BadParameter(msg)
            elif '<<' in directive:
                (local, remote) = directive.split('<<')
                self._out_append[local.strip()] = remote.strip()
            elif '>>' in directive:
                (local, remote) = directive.split('>>')
                self._in_append[local.strip()] = remote.strip()
            elif '<' in directive:
                (local, remote) = directive.split('<')
                self._out_overwrite[local.strip()] = remote.strip()
            elif '>' in directive:
                (local, remote) = directive.split('>')
                self._in_overwrite[local.strip()] = remote.strip()
            else:
                msg = "'%s' is not a valid transfer directive string."
                raise se.BadParameter(msg)

    def _dicts_to_string_list(self):
        slist = list()
        for (local, remote) in self._in_overwrite.iteritems():
            slist.append('%s > %s' % (local, remote))
        for (local, remote) in self._in_append.iteritems():
            slist.append('%s >> %s' % (local, remote))
        for (remote, local) in self._out_overwrite.iteritems():
            slist.append('%s < %s' % (local, remote))
        for (remote, local) in self._out_append.iteritems():
            slist.append('%s << %s' % (local, remote))
        return slist

    def __str__(self):
        """ String representation.
        """
        return str(self._dicts_to_string_list())

    @property
    def in_overwrite_dict(self):
        return self._in_overwrite

    @property
    def in_append_dict(self):
        return self._in_append

    @property
    def out_overwrite_dict(self):
        return self._out_overwrite

    @property
    def out_append_dict(self):
        return self._out_append

    @property
    def string_list(self):
        return self._dicts_to_string_list()


def _test_():
    tdp = TransferDirectives(["ab","a>c", "c>>d","f<a","g<<h"])
    print tdp.in_append_dict
    print tdp.in_overwrite_dict
    print tdp.out_append_dict
    print tdp.out_overwrite_dict





