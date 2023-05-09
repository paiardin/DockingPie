
__author__    = "Radical.Utils Development Team (Andre Merzky, Ole Weidner)"
__copyright__ = "Copyright 2013, RADICAL@Rutgers"
__license__   = "MIT"


import os

def which(program):
    '''taken from:
        http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    '''
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)

    if  fpath:
        if  is_exe (program) :
            return program
    else:
        for path in os.environ["PATH"].split (os.pathsep) :
            exe_file = os.path.join (path, program)
            if  is_exe (exe_file) :
                return exe_file

    return None



