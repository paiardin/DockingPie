
__author__    = "Andre Merzky, Ole Weidner"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


""" the adaptor base class. """

import radical.utils         as ru
import radical.utils.config  as ruc
import radical.utils.logger  as rul


from   saga.exceptions import *


# ------------------------------------------------------------------------------
# adaptor base class
#
class Base (ruc.Configurable) :

    # We only need one instance of this adaptor per process (actually per
    # engine, but engine is a singleton, too...) -- the engine will though
    # create new CPI implementation instances as needed (one per SAGA API
    # object).
    __metaclass__ = ru.Singleton

    
    # --------------------------------------------------------------------------
    #
    # FIXME: adaptor_options type...
    #
    def __init__ (self, adaptor_info, adaptor_options=[]) :

        self._info    = adaptor_info
        self._opts    = adaptor_options
        self._name    = adaptor_info['name']
        self._schemas = adaptor_info['schemas']

        self._lock    = ru.RLock      (self._name)
        self._logger  = ru.get_logger('radical.saga.api')

        has_enabled = False
        for option in self._opts :
            if option['name'] == 'enabled' :
                has_enabled = True

        if not has_enabled :
            # *every* adaptor needs an 'enabled' option!
            self._opts.append ({ 
                'category'         : self._name,
                'name'             : 'enabled', 
                'type'             : bool, 
                'default'          : True, 
                'valid_options'    : [True, False],
                'documentation'    : "Enable / disable loading of the adaptor",
                'env_variable'     : None
                }
            )


        ruc.Configurable.__init__       (self, 'saga')
        ruc.Configurable.config_options (self, self._name, self._opts)


    # --------------------------------------------------------------------------
    #
    # if sanity_check() is commented out here, then we will only load adaptors
    # which implement the method themselves.
    #
    def sanity_check (self) :
        """ This method can be overloaded by adaptors to check runtime
            conditions on adaptor load time.  The adaptor should raise an
            exception if it will not be able to function properly in the given
            environment, e.g. due to missing dependencies etc.
        """
        raise BadParameter ("Adaptor %s does not implement sanity_check()"  \
                         % self._name)


    # --------------------------------------------------------------------------
    #
    def register (self) :
        """ Adaptor registration function. The engine calls this during startup
            to retrieve the adaptor information.
        """

        return self._info


    # --------------------------------------------------------------------------
    #
    def get_name (self) :
        return self._name


    # --------------------------------------------------------------------------
    #
    def get_schemas (self) :
        return self._schemas


    # --------------------------------------------------------------------------
    #
    def get_info (self) :
        return self._info





