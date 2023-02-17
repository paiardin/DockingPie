
__author__    = "Radical.Utils Development Team (Andre Merzky, Ole Weidner)"
__copyright__ = "Copyright 2013, RADICAL@Rutgers"
__license__   = "MIT"


import radical.utils.testing as rut


# ------------------------------------------------------------------------------
#
class TestConfig (rut.TestConfig): 

    #---------------------------------------------------------------------------
    # 
    def __init__ (self, cfg_file):

        # initialize configuration.  We only use the 'radical.utils.tests' 
        # category from the config file.
        rut.TestConfig.__init__ (self, cfg_file, 'radical.utils.tests')


# ------------------------------------------------------------------------------


