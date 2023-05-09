
__author__    = "Andre Merzky, Ole Weidner"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import saga.context
import saga.adaptors.base
import saga.adaptors.cpi.context

SYNC_CALL  = saga.adaptors.cpi.decorators.SYNC_CALL
ASYNC_CALL = saga.adaptors.cpi.decorators.ASYNC_CALL

######################################################################
#
# adaptor meta data
#
_ADAPTOR_NAME          = 'saga.adaptor.userpass'
_ADAPTOR_SCHEMAS       = ['UserPass']
_ADAPTOR_OPTIONS       = []

_ADAPTOR_CAPABILITIES  = {
    'attributes'       : [saga.context.TYPE,
                          saga.context.USER_ID,
                          saga.context.USER_PASS]
}

_ADAPTOR_DOC           = {
    'name'             : _ADAPTOR_NAME,
    'cfg_options'      : _ADAPTOR_OPTIONS, 
    'capabilities'     : _ADAPTOR_CAPABILITIES,
    'description'      : """This adaptor stores user_id and user_pass tokens, to
                            be used for backend connections.""",
    'schemas'          : {'userpass' : 'this adaptor can only store username/password pairs.'},
}

_ADAPTOR_INFO          = {
    'name'             : _ADAPTOR_NAME,
    'version'          : 'v0.1',
    'schemas'          : _ADAPTOR_SCHEMAS,
    'cpis'             : [{ 
        'type'         : 'saga.Context',
        'class'        : 'ContextUserPass'
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


    def sanity_check (self) :
        pass


    def _get_default_contexts (self) :

        return self._default_contexts



######################################################################
#
# job adaptor class
#
class ContextUserPass (saga.adaptors.cpi.context.Context) :

    def __init__ (self, api, adaptor) :

        _cpi_base = super  (ContextUserPass, self)
        _cpi_base.__init__ (api, adaptor)


    @SYNC_CALL
    def init_instance (self, adaptor_state, type) :
        
        if not type.lower () in (schema.lower() for schema in _ADAPTOR_SCHEMAS) :
            raise saga.exceptions.BadParameter \
                    ("the UserPass context adaptor only handles UserPass contexts - duh!")

        self._type = type

        return self


    @SYNC_CALL
    def _initialize (self, session) :
        pass





