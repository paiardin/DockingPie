
__author__    = "Andre Merzky, Ole Weidner"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import os
import subprocess

import saga.context
import saga.adaptors.base
import saga.adaptors.cpi.decorators
import saga.adaptors.cpi.context

SYNC_CALL  = saga.adaptors.cpi.decorators.SYNC_CALL
ASYNC_CALL = saga.adaptors.cpi.decorators.ASYNC_CALL

######################################################################
#
# adaptor meta data
#
_ADAPTOR_NAME          = 'saga.adaptor.myproxy'
_ADAPTOR_SCHEMAS       = ['MyProxy']
_ADAPTOR_OPTIONS       = [
    {
    'category'         : 'saga.adaptor.myproxy',
    'name'             : 'base_workdir',
    'type'             : str,
    'default'          : "$HOME/.saga/adaptors/proxies/",
    'documentation'    : '''The adaptor stores proxies information on the
                          filesystem on the target resource.  This parameter
                          specified what location should be used.''',
    'env_variable'     : None
    }
]

_ADAPTOR_CAPABILITIES  = {
    'attributes'       : [saga.context.TYPE,
                          saga.context.SERVER,
                          saga.context.USER_ID,
                          saga.context.USER_PASS,
                          saga.context.LIFE_TIME]
}

_ADAPTOR_DOC           = {
    'name'             : _ADAPTOR_NAME,
    'cfg_options'      : _ADAPTOR_OPTIONS, 
    'capabilities'     : _ADAPTOR_CAPABILITIES,
    'description'      : """This adaptor fetches an X509 proxy from
                            MyProxy when it is added to a saga.Session.""",
    'schemas'          : {'myproxy' : 'this adaptor can only interact with myproxy backends'},
}

_ADAPTOR_INFO          = {
    'name'             : _ADAPTOR_NAME,
    'version'          : 'v0.1',
    'schemas'          : _ADAPTOR_SCHEMAS,
    'cpis'             : [{ 
        'type'         : 'saga.Context',
        'class'        : 'ContextMyProxy'
        }
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

    def __init__ (self) :

        saga.adaptors.base.Base.__init__ (self, _ADAPTOR_INFO, _ADAPTOR_OPTIONS)

        # there are no default myproxy contexts
        self._default_contexts = []
        self.base_workdir = self.opts['base_workdir'].get_value()


    def sanity_check (self) :
        pass


    def _get_default_contexts (self) :

        return self._default_contexts



######################################################################
#
# job adaptor class
#
class ContextMyProxy (saga.adaptors.cpi.context.Context) :

    def __init__ (self, api, adaptor) :

        _cpi_base = super  (ContextMyProxy, self)
        _cpi_base.__init__ (api, adaptor)
        self.base_workdir = adaptor.base_workdir


    @SYNC_CALL
    def init_instance (self, adaptor_state, type) :

        if not type.lower () in (schema.lower() for schema in _ADAPTOR_SCHEMAS) :
            raise saga.exceptions.BadParameter \
                    ("the MyProxy context adaptor only handles MyProxy contexts - duh!")

        self.get_api ().type = type

        return self.get_api ()


    @SYNC_CALL
    def _initialize (self, session) :

        # make sure we have server, username, password
        api = self.get_api ()


        # set up the myproxy command
        cmd = ""
        
        if api.user_pass :
            cmd = "echo %s | myproxy-logon --stdin_pass"  %  api.user_pass
        else :
            cmd = "myproxy-logon --stdin_pass"

        if api.server :
            if ':' in api.server:
                (server, port) = api.server.split (':', 2)
            else:
                server = api.server
                port = "7512"
            if server    : cmd += " --pshost %s"          %  server
            if port      : cmd += " --psport %s"          %  port

        if  api.user_id : 
            cmd += " --username %s"        %  api.user_id
        
        if  api.life_time and api.life_time > 0 : 
            cmd += " --proxy_lifetime %s"  %  api.life_time

        # store the proxy in a private location
        proxy_store    = "%s" % self.base_workdir
        proxy_location = "%s/myproxy_%d.x509"  %  (proxy_store, id(self))

        if not os.path.exists (proxy_store) :
            try :
                os.makedirs (proxy_store)
            except OSError as e :
                raise saga.exceptions.NoSuccess ("could not create myproxy store: %s"  %  str(e))

        cmd += " --out %s"  %  proxy_location

        expected_result  = "A credential has been received for user %s in %s.\n" \
                         %  (api.user_id, proxy_location)

        process          = subprocess.Popen (cmd, shell=True,
                                             stdout=subprocess.PIPE,
                                             stderr=subprocess.PIPE)
        (stdout, stderr) = process.communicate ()

        if expected_result == stdout :
            self._logger.info (stdout)
        else :
            self._logger.info (stderr)
            raise saga.exceptions.BadParameter ("could not evaluate myproxy context: %s"  %  stderr)

        new_ctx = saga.Context ('X509')

        new_ctx.user_proxy = proxy_location
        new_ctx.life_time  = api.life_time

        session.add_context (new_ctx)




