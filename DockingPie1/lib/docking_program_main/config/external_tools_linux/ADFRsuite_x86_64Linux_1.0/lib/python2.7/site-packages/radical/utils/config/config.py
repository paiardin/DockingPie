
__author__    = "Radical.Utils Development Team (Andre Merzky, Ole Weidner)"
__copyright__ = "Copyright 2013, RADICAL@Rutgers"
__license__   = "MIT"


''' Provides API handles for configuration systems. '''

import os
import copy
import threading

from ..singleton import Singleton

from   configfile import ConfigFileReader


# ------------------------------------------------------------------------------
#
def getConfig(name):
    """ 
    Returns a handle to the global configuration object.
    """
    return Configuration(name) 


# ------------------------------------------------------------------------------
#
class ConfigOption(object):
    """ 
    Represent a (mutable) configuration option.
    """

    # --------------------------------------------------------------------------
    #
    def __init__(self, category, name, val_type, default_val, valid_options,
                 documentation, env_var):

        self._category      = category
        self._name          = name
        self._val_type      = val_type
        self._default_val   = default_val
        self._valid_options = valid_options
        self._env_var       = env_var
        self._documentation = documentation
        self._value         = None


    # --------------------------------------------------------------------------
    #
    def __str__(self):
        return str({'name':self._name, 'value':self._value})


    # --------------------------------------------------------------------------
    #
    def as_dict(self):
        return {self._name:self._value}


    # --------------------------------------------------------------------------
    #
    @property
    def category(self):
        return self._category


    # --------------------------------------------------------------------------
    #
    def set_value(self, value):
        # make sure we got the right value type
        if type(value) != self._val_type:
            raise ValueTypeError(self._category, self._name, 
              type(value), self._val_type)

        self._value = value


    # --------------------------------------------------------------------------
    #
    def get_value(self):
        return self._value


# ------------------------------------------------------------------------------
#
# a singleton class which keeps all configurations
#
class _Configurations (object) :
    """
    This singleton class maintains references to all global configurations.
    """
    __metaclass__ = Singleton
    _rlock = threading.RLock ()

    # --------------------------------------------------------------------------
    #
    def __init__ (self) :


        with self._rlock :
            self._configs = {}


    # --------------------------------------------------------------------------
    #
    def put (self, name, config) :

        with self._rlock :
            self._configs[name] = dict (config)
          # self._configs[name] = copy.deepcopy (config)


    # --------------------------------------------------------------------------
    #
    def get (self, name) :

        with self._rlock :
            if  not name in self._configs :
                self._configs[name] = dict()

            return dict (self._configs[name])
          # return copy.deepcopy (self._configs[name])


# ------------------------------------------------------------------------------
#
_configurations = _Configurations ()


# ------------------------------------------------------------------------------
#
class Configuration (object): 
    """ 
    Represents a global configuration for a given namespace.

    The Configuration class can be used to introspect and modify configuration
    options.  It uses a 'Singleton' data store, which means that  multiple
    instances all point to the same global configuration, for a specific given
    namespace.
    """    

    _rlock = threading.RLock ()


    # --------------------------------------------------------------------------
    #
    def __init__ (self, name):

        with self._rlock :

            self._name = name

            # initialize and sync options
            self._initialize()


    # --------------------------------------------------------------------------
    #
    def _put_state (self):

        with self._rlock :

            # push state back into the store
            state = {}
            state['all_valid_options'] = self._all_valid_options
            state['master_config']     = self._master_config

            _configurations.put (self._name, state)


    # --------------------------------------------------------------------------
    #
    def _get_state (self):

        with self._rlock :

            # populate this instance with previously derived config options
            state = _configurations.get (self._name)

            # _all_valid_options contains the raw option sets
            # added to configuration by 'Configurable' classes
            if 'all_valid_options' in state :
                self._all_valid_options = state['all_valid_options']
            else :
                self._all_valid_options = list()

            # _master config holds the parsed options, including 
            # values fetched from configuartion files and env variables
            if 'master_config' in state :
                self._master_config = state['master_config']
            else :
                self._master_config = dict()



    # --------------------------------------------------------------------------
    #
    def _update (self, namespace, valid_options):

        with self._rlock :

            # sync for singletonism
            self._get_state ()

            # add the new options to the global dictionary
            self._all_valid_options += valid_options

            # sync for singletonism
            self._put_state ()

            # re-_initialize the configuration object.
            self._initialize()

    # --------------------------------------------------------------------------
    #
    def _initialize(self, inject_cfg_file=None, add_cfg_file=None):
        """ Initialize the global configuration.

            :param inject_cfg_file: is used *only* for testing purposes 
             and overwrites / ignores the regular config file locations 
             /etc/<name>.cfg & $HOME/.<name>.cfg

            :param add_cfg_file: is used *only* for testing purposes 
             and adds a specific config file to the list of evaluated files.  
        """
        with self._rlock :

            # get state from cache
            self._get_state ()

            cfg_files = list()
            if inject_cfg_file is not None:
                cfg_files.append(inject_cfg_file)
            else:
                # check for the existence of regular configuration files
                sys_cfg = '/etc/%s.cfg' % self._name
                if os.path.exists(sys_cfg):
                    cfg_files.append(sys_cfg)

                uconf = "%s_CONFIG" % self._name.upper ()
                uconf = uconf.replace ('.', '_')
                uconf = uconf.replace ('-', '_')
                if uconf in os.environ :
                    usr_cfg = os.environ[uconf]
                    if  not os.path.exists (usr_cfg) :
                        print "WARNING: %s set to %s, but file does not exist" % (uconf, usr_cfg)
                else:
                    usr_cfg = '%s/.%s.cfg' % (os.path.expanduser("~"), self._name)

                if os.path.exists(usr_cfg):
                    cfg_files.append(usr_cfg)

            if add_cfg_file is not None:
                cfg_files.append(add_cfg_file)

            cfr = ConfigFileReader(cfg_files)
    
            # cfg_file_dict holds all configuration variables that 
            # were read from either a system-wide or user configuration file.
            cfg_file_dict = cfr.get_config_dict()


            # load valid options and add them to the configuration
            for option in self._all_valid_options:
                cat = option['category']
                if cat not in self._master_config:
                    # first occurrence - add new category key
                    self._master_config[cat] = dict()

                ev = None
                if option['env_variable'] :
                    ev = os.environ.get(option['env_variable'])

                if (option['category'] in cfg_file_dict) and \
                    option['name']     in cfg_file_dict[option['category']]:
                    # found entry in configuration file -- use it
                    tmp_value = cfg_file_dict[option['category']][option['name']]
                    # some list types need type conversion
                    if option['type'] == list:
                        value = tmp_value.split(",")
                    elif option['type'] == bool:
                        if tmp_value.lower() == 'true':
                          value = True
                        elif tmp_value.lower() == 'false':
                          value = False 
                        else:
                          raise ValueTypeError(option['category'], option['name'],
                              tmp_value, option['type'])
                    elif option['type'] == int:
                        value = int(tmp_value)
                    elif option['type'] == float:
                        value = float(tmp_value)
                    else:
                        value = str(tmp_value)

                elif ev is not None:
                    tmp_value = ev
                    if option['type'] == float:
                        value = float(tmp_value)
                    elif option['type'] == int:
                        value = int(tmp_value)
                    elif option['type'] == list:
                        value = tmp_value.split(",")
                    elif option['type'] == bool:
                        if tmp_value.lower() == 'true':
                          value = True
                        elif tmp_value.lower() == 'false':
                          value = False 
                        elif tmp_value == '1':
                          value = True
                        elif tmp_value == '0':
                          value = False
                        else:
                          raise ValueTypeError(option['category'], option['name'],
                              tmp_value, option['type'])
                    elif option['type'] == int:
                        value = int(tmp_value)
                    elif option['type'] == float:
                        value = float(tmp_value)
                    else:
                        value = tmp_value
                else:
                    value = option['default']

                if not 'valid_options' in option :
                    option['valid_options'] = None

                self._master_config[cat][option['name']] = ConfigOption(
                    option['category'],
                    option['name'],
                    option['type'],
                    option['default'],
                    option['valid_options'],
                    option['documentation'],
                    option['env_variable'])

                self._master_config[cat][option['name']].set_value(value) 


            # now walk through the cfg_file_dict -- for all entries not yet handled
            # (i.e. which have not been registered before), use default
            # ConfigOptions:
            #       category       : from cfg_file_dict
            #       name           : from cfg_file_dict
            #       value          : from cfg_file_dict
            #       type           : string
            #       default        : ""
            #       valid_options  : None
            #       documentation  : ""
            #       env_variable   : None
            #
            # If later initialize_ is called again, and that option has been
            # registered meanwhile, this entry will be overwritten.

            for section in cfg_file_dict.keys () :
                for name in cfg_file_dict[section].keys () :
                    value = cfg_file_dict[section][name]

                    # check if this is a registered entry
                    exists = False
                    if section in self._master_config :
                        if name in self._master_config[section] :
                            exists = True
                    else :
                        # section does not exist - create it, and add option
                        self._master_config[section] = {}

                    # if not, add it with default ConfigOptions
                    if not exists :
                        self._master_config[section][name] = ConfigOption(
                            section   , # category
                            name      , # name
                            str       , # type
                            ""        , # default
                            None      , # valid_options
                            ""        , # documentation
                            None      ) # env_variable

                        self._master_config[section][name].set_value(value)


            # sync cache with newly configured state
            self._put_state ()

    # --------------------------------------------------------------------------
    #
    def has_category(self, category_name):
        """ Check for a specific configuration category.
        """
        with self._rlock :

            # sync for singletonism
            self._get_state ()

            if category_name not in self._master_config:
                return False
            else:
                return True

    # --------------------------------------------------------------------------
    #
    def get_category(self, category_name):
        """ Return a specific configuration category.
        """
        with self._rlock :

            # sync for singletonism
            self._get_state ()

            if category_name not in self._master_config:
                raise CategoryNotFound(category_name)
            else:
                return self._master_config[category_name]

    # --------------------------------------------------------------------------
    #
    def has_option(self, category_name, option_name):

        with self._rlock :

            # sync for singletonism
            self._get_state ()

            if category_name not in self._master_config:
                return False
            else:
                if option_name not in self._master_config[category_name]:
                    return False
                else:
                    return True


    # --------------------------------------------------------------------------
    #
    def get_option(self, category_name, option_name):

        with self._rlock :

            # sync for singletonism
            self._get_state ()

            if category_name not in self._master_config:
                raise CategoryNotFound(category_name)
            else:
                if option_name not in self._master_config[category_name]:
                    raise OptionNotFound(category_name, option_name)
                else:
                    return self._master_config[category_name][option_name]


    # --------------------------------------------------------------------------
    #
    def as_dict (self, cn = None) :

        with self._rlock :

            # sync for singletonism
            self._get_state ()

            ret = {}

            if  cn : 
                for on in self._master_config[cn] :
                    ret[on] = self.get_option (cn, on).get_value ()

            else :
                for cn in self._master_config :
                    ret[cn] = {}
  
                    for on in self._master_config[cn] :
                        ret[cn][on] = self.get_option (cn, on).get_value ()

            return ret



# ------------------------------------------------------------------------------
#
class Configurable (object) :
    """ 
    This class provides an interface for all configurable objects.  
    """

    # --------------------------------------------------------------------------
    #
    def __init__(self, name):

        self._name = name

    # --------------------------------------------------------------------------
    #
    def config_options (self, namespace, valid_options):

        ## register a new 'configurable' object
        getConfig(self._name)._update(namespace, valid_options)


    # --------------------------------------------------------------------------
    #
    def get_config(self, namespace):
        return getConfig(self._name).get_category(namespace)



    # --------------------------------------------------------------------------
    #
    def get_config_as_dict (self): return getConfig(self._name).as_dict ()
    def            as_dict (self): return getConfig(self._name).as_dict ()


# ------------------------------------------------------------------------------
#
class CategoryNotFound(LookupError):
    # --------------------------------------------------------------------------
    #
    def __init__(self, name):
        msg = "A category with name '%s' could not be found." % name
        LookupError.__init__(self, msg)


# ------------------------------------------------------------------------------
#
class OptionNotFound(LookupError):
    # --------------------------------------------------------------------------
    #
    def __init__(self, category_name, option_name):
        name = "%s.%s" % (category_name, option_name)
        msg  = "An option with name '%s' could not be found." % (name)
        LookupError.__init__(self, msg)


# ------------------------------------------------------------------------------
#
class ValueTypeError(TypeError):
    # --------------------------------------------------------------------------
    #
    def __init__(self, category_name, option_name, value_type, required_type):
        name = "%s.%s" % (category_name, option_name)
        msg  = "Option %s requires value of type '%s' but got '%s'." \
             % (name, required_type, value_type)
        TypeError.__init__ (self, msg)


# ------------------------------------------------------------------------------
#


