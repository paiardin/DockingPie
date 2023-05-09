
__author__    = "Andre Merzky, Ashley Z, Ole Weidner"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


""" shell based resource adaptor implementation """

import saga.utils.pty_shell

import saga.adaptors.cpi.base
import saga.adaptors.cpi.resource

from   saga.resource.constants import *
ANY = COMPUTE | STORAGE

import re
import os
import time

SYNC_CALL  = saga.adaptors.cpi.decorators.SYNC_CALL
ASYNC_CALL = saga.adaptors.cpi.decorators.ASYNC_CALL


# --------------------------------------------------------------------
# the adaptor info
#
_ADAPTOR_NAME          = "saga.adaptor.shell_resource"
_ADAPTOR_SCHEMAS       = ["local", "shell"]
_ADAPTOR_OPTIONS       = []

# --------------------------------------------------------------------
# the adaptor capabilities & supported attributes
#
_ADAPTOR_CAPABILITIES  = {
    "rdes_attributes"  : [saga.resource.RTYPE         ,
                          saga.resource.MACHINE_OS    ,
                          saga.resource.MACHINE_ARCH  ,
                          saga.resource.SIZE          ,
                          saga.resource.MEMORY        ,
                          saga.resource.ACCESS       ],
    "res_attributes"   : [saga.resource.RTYPE         ,
                          saga.resource.MACHINE_OS    ,
                          saga.resource.MACHINE_ARCH  ,
                          saga.resource.SIZE          ,
                          saga.resource.MEMORY        ,
                          saga.resource.ACCESS       ],    
    "metrics"          : [saga.resource.STATE, 
                          saga.resource.STATE_DETAIL],
    "contexts"         : {"ssh"      : "public/private keypair",
                          "x509"     : "X509 proxy for gsissh",
                          "userpass" : "username/password pair for ssh"}
}

# --------------------------------------------------------------------
# the adaptor documentation
#
_ADAPTOR_DOC           = {
    "name"             : _ADAPTOR_NAME,
    "cfg_options"      : _ADAPTOR_OPTIONS, 
    "capabilities"     : _ADAPTOR_CAPABILITIES,
    "description"      : """ 
        The Shell resource adaptor. This adaptor attempts to determine what job
        submission endpoint and file system resources are available for a given
        host, and provides the respective access URLs.
        """,
    "example": "examples/jobs/localresource.py",
    "schemas"          : {"shell"  : "find access URLs for shell-job/file adaptors"}
}

# --------------------------------------------------------------------
# the adaptor info is used to register the adaptor with SAGA

_ADAPTOR_INFO          = {
    "name"             : _ADAPTOR_NAME,
    "version"          : "v0.1.beta",
    "schemas"          : _ADAPTOR_SCHEMAS,
    "cpis"             : [
        { 
        "type"         : "saga.resource.Manager",
        "class"        : "ShellResourceManager"
        }, 
        { 
        "type"         : "saga.resource.Compute",
        "class"        : "ShellResourceCompute"
        },
    ]
}

###############################################################################
# The adaptor class

class Adaptor (saga.adaptors.base.Base):
    """ 
    This is the actual adaptor class, which gets loaded by SAGA (i.e. by the
    SAGA engine), and which registers the CPI implementation classes which
    provide the adaptor's functionality.
    """


    # ----------------------------------------------------------------
    #
    def __init__ (self) :

        self.base = super  (Adaptor, self)
        self.base.__init__ (_ADAPTOR_INFO, _ADAPTOR_OPTIONS)

        # for id parsing
        self.id_re = re.compile ('^\[(.*)\]-\[(.*?)\]$')


    # ----------------------------------------------------------------
    #
    def sanity_check (self) :

        # FIXME: also check for gsissh

        pass


    # ----------------------------------------------------------------
    #
    def parse_id (self, id) :
        # split the id '[manager-url]-[resource-url]' in its parts, and return them.

        print id
        match = self.id_re.match (id)

        if  not match or len (match.groups()) != 2 :
            raise saga.BadParameter ("Cannot parse resource id '%s'" % id)

        return (saga.Url (match.group(1)), saga.Url (match.group (2)))




###############################################################################
#
class ShellResourceManager (saga.adaptors.cpi.resource.Manager) :

    # ----------------------------------------------------------------
    #
    def __init__ (self, api, adaptor) :

        self._cpi_base = super  (ShellResourceManager, self)
        self._cpi_base.__init__ (api, adaptor)



    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def init_instance (self, adaptor_state, url, session) :

        self.url     = saga.Url (url)  # deep copy
        self.session = session
        self.access  = {}
        self.access[COMPUTE] = []
        self.access[STORAGE] = []
        self.access[ANY]     = []

        # check for compute entry points
        for schema in ['fork', 'ssh', 'gsissh'] :
            tmp_url = saga.Url (self.url)  # deep copy
            tmp_url.schema = schema

            shell = saga.utils.pty_shell.PTYShell (tmp_url, self.session, self._logger)

            if  shell.alive () :
                self.access[COMPUTE].append (tmp_url)
                self.access[ANY]    .append (tmp_url)
                shell.finalize (True)


        # check for storage entry points
        for schema in ['file', 'sftp', 'gsisftp'] :
            tmp_url = saga.Url (self.url)  # deep copy
            tmp_url.schema = schema

            shell = saga.utils.pty_shell.PTYShell (tmp_url, self.session, self._logger)

            if  shell.alive () :
                self.access[STORAGE].append (tmp_url)
                self.access[ANY]    .append (tmp_url)
                shell.finalize (True)


    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def acquire (self, rd) :

        if  not rd :
            raise saga.BadParameter._log (self._logger, "acquire needs a resource description")

        if  rd.rtype != COMPUTE and \
            rd.rtype != STORAGE     :
            raise saga.BadParameter._log (self._logger, "can only acquire compute and storage resources.")


        # check that only supported attributes are provided
        for attribute in rd.list_attributes():
            if attribute not in _ADAPTOR_CAPABILITIES["rdes_attributes"]:
                msg = "'resource.Description.%s' is not supported by this adaptor" % attribute
                raise saga.BadParameter._log (self._logger, msg)

        if  rd.access :
            access_url = saga.Url (rd.access) 
            if  not access_url in self.access[rd.rtype] :
                msg = "access '%s' is not supported by this backend" % rd.access
                raise saga.BadParameter._log (self._logger, msg)

        if  not len (self.access[rd.rtype]) :
            raise saga.BadParameter._log (self._logger, "resource type is not supported by this backend")

        
        # this dict is passed on to the job adaptor class -- use it to pass any
        # state information you need there.
        adaptor_state = { "resource_access_url"     : self.access[rd.rtype][0], 
                          "resource_type"           : rd.rtype, 
                          "resource_manager"        : self.get_api (), 
                          "resource_manager_url"    : self.url, 
                          "resource_schema"         : self.url.schema }

        if rd.rtype == COMPUTE :
            return saga.resource.Compute (_adaptor=self._adaptor, _adaptor_state=adaptor_state)

        if rd.rtype == STORAGE :
            return saga.resource.Storage (_adaptor=self._adaptor, _adaptor_state=adaptor_state)


    # ----------------------------------------------------------------
    @SYNC_CALL
    def get_url (self) :

        return self.url


    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def list (self, rtype):

        return self.access[rtype]
   
   
    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def release (self, id):

        return # hahaha

   
    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def list_templates (self, rtype) :

        return [] # no templates

   
    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_template (self, name) :

        raise saga.BadParameter ("unknown template %s" % name)


###############################################################################
#
class ShellResourceCompute (saga.adaptors.cpi.resource.Compute) :

    # ----------------------------------------------------------------
    #
    def __init__ (self, api, adaptor) :

        self._cpi_base = super  (ShellResourceCompute, self)
        self._cpi_base.__init__ (api, adaptor)

        self.state       = ACTIVE
        self.rtype       = None
        self.manager     = None
        self.manager_url = None
        self.access_url  = None


    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def init_instance (self, adaptor_info, id, session):

        # eval id if given
        if  id :
            self.manager_url, self.access_url = self._adaptor.parse_id (id)
            self.manager = saga.resource.Manager (self.manager_url)

            if  self.access_url.scheme in self.manager.list (COMPUTE) :
                self.rtype = COMPUTE

            elif self.access_url.scheme in self.manager.list (STORAGE) :
                self.rtype = STORAGE

            else :
                raise saga.BadParameter ("Cannot handle resource type for %s", id)

        # no id -- grab info from adaptor_info
        elif adaptor_info :

            if  not 'resource_access_url'  in adaptor_info or \
                not 'resource_type'        in adaptor_info or \
                not 'resource_manager'     in adaptor_info or \
                not 'resource_manager_url' in adaptor_info    :
                raise saga.BadParameter ("Cannot acquiure resource, insufficient information")

            self.access_url  = adaptor_info['resource_access_url']
            self.rtype       = adaptor_info['resource_type']
            self.manager     = adaptor_info['resource_manager']
            self.manager_url = adaptor_info['resource_manager_url']

            self.id = "[%s]-[%s]" % (self.manager_url, self.access_url)

        else :
            raise saga.BadParameter ("Cannot acquire resource, no contact information")


        return self.get_api ()


    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def get_id           (self) : return self.id

    @SYNC_CALL
    def get_rtype        (self) : return self.rtype

    @SYNC_CALL
    def get_state        (self) : return self.state

    @SYNC_CALL
    def get_state_detail (self) : return None

    @SYNC_CALL
    def get_access       (self) : return self.access_url

    @SYNC_CALL
    def get_manager      (self) : return self.manager

    @SYNC_CALL
    def get_description  (self) : return { ACCESS : self.access_url }

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def reconfig (self):
        raise saga.NotImplemented ("This backend cannot reconfigre resources")


    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def release (self):

        return # hahahah

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def wait (self, state, timeout) : 
        # trick is, we *never* change state...

        if  state == self.state :
            return

        if  timeout >= 0 :
            time.sleep (timeout)
            return

        if  timeout < 0 :
            while True :
                time.sleep (10)
        
        # we never get here...




