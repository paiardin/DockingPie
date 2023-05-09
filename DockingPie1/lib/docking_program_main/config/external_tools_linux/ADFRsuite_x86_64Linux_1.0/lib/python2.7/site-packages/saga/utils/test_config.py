
__author__    = "Andre Merzky, Ole Weidner"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import copy

import radical.utils.testing  as rut
import radical.utils.logger   as rul

import saga.exceptions        as se

import saga

# ------------------------------------------------------------------------------
#
def add_tc_params_to_jd (tc, jd):

    if 'job_walltime_limit' in tc and tc['job_walltime_limit'] != "": 
        jd.wall_time_limit = tc['job_walltime_limit'] 
    if 'job_project' in tc and tc['job_project'] != "": 
        jd.project = tc['job_project']
    if 'job_queue' in tc and tc['job_queue'] != "":
        jd.queue = tc['job_queue']
    if 'job_total_cpu_count' in tc and tc['job_total_cpu_count'] != "":
        jd.total_cpu_count = tc['job_total_cpu_count']
    if 'job_spmd_variation' in tc and tc['job_spmd_variation'] != "": 
        jd.spmd_variation = tc['job_spmd_variation'] 

    return jd

# ------------------------------------------------------------------------------
#
class TestConfig (rut.TestConfig): 

    #-----------------------------------------------------------------
    # 
    def __init__ (self, cfg_file):

        # initialize configuration.  We only use the 'saga.tests' category from
        # the config file.
        rut.TestConfig.__init__ (self, cfg_file, 'saga.tests')

        # setup a saga session for the tests
        # don't populate session with default contexts...
        self.session = saga.Session (default=False) 

        # attempt to create a context from the test config
        if  self.context_type :

            c = saga.Context (self.context_type)

            c.user_id    = self.context_user_id
            c.user_pass  = self.context_user_pass
            c.user_cert  = self.context_user_cert
            c.user_proxy = self.context_user_proxy

            # add it to the session
            self.session.add_context (c)


# ------------------------------------------------------------------------------


