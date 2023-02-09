
__author__    = "Andre Merzky, Ole Weidner"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


"""
List of all registered SAGA adaptor modules.

This registry is used to locate and load adaptors.  The entries must be
formatted in dotted python module notation (e.g. as
'saga.adaptors.context.x509'), and the Python module search path has to be
configured so that the listed modules can be loaded.

Note that a module listed in the registry is not guaranteed to be available:
while the engine will attempt to load the module, it may be disabled by
configure options, or due to missing pre-requisites.
"""

adaptor_registry = [
                    "saga.adaptors.context.myproxy",
                    "saga.adaptors.context.x509",
                    "saga.adaptors.context.ssh",
                    "saga.adaptors.context.userpass",
                    "saga.adaptors.shell.shell_job",
                    "saga.adaptors.shell.shell_file",
                    "saga.adaptors.shell.shell_resource",
                    "saga.adaptors.redis.redis_advert",
                    "saga.adaptors.sge.sgejob",
                    "saga.adaptors.pbs.pbsjob",
                    "saga.adaptors.lsf.lsfjob",
                    # "saga.adaptors.irods.irods_replica",
                    "saga.adaptors.condor.condorjob",
                    "saga.adaptors.slurm.slurm_job",
                    "saga.adaptors.http.http_file",
                    "saga.adaptors.aws.ec2_resource",
                    "saga.adaptors.loadl.loadljob",
                    "saga.adaptors.globus_online.go_file",
                    "saga.adaptors.torque.torquejob",
                    "saga.adaptors.pbspro.pbsprojob"
                   ]
