
__author__    = "Radical.Utils Development Team (Andre Merzky, Ole Weidner)"
__copyright__ = "Copyright 2013, RADICAL@Rutgers"
__license__   = "MIT"


''' Provides utilities to work with configuration files and query strings.  '''

from copy import deepcopy
import ConfigParser, os

class ConfigFileReader(object):

    def __init__(self, filenames):
        config = ConfigParser.RawConfigParser()
        config.read(filenames)

        self.config_dict = dict()
        
        for section in config.sections():
            self.config_dict[section] = dict()
            for option in config.options(section):
                self.config_dict[section][option] = config.get(section, option)

    def get_config_dict(self):
        return self.config_dict

class BaseConfig(object):
    
    def __init__(self):
        self._query_string = ""
        self._dictionary = {}

    def _re_eval(self, query_string=None):
        if query_string is not None:
            # initial eval from query string
            self._query_string = query_string
            for param in self._query_string.split('&'):
                try:
                    (key, value) = param.split("=")
                    self._dictionary[key] = value
                except ValueError:
                    if param != '':
                        self._dictionary[param] = ""
        else:
            # re-genereate query string
            new_query_string = ''
            for (key, value) in self._dictionary.iteritems():
                if value is None:
                    new_query_string += "&%s" % key
                else:
                    new_query_string += "&%s=%s" % (key, value)
            if new_query_string.startswith("&"):
                new_query_string = new_query_string[1:]
            self._query_string = new_query_string
            

    def as_query_string(self, subsection_name=None):
        ''' Return the configuration (or a subsection of it) as a URL 
            query string.
        '''
        self._re_eval()
        return deepcopy(self._query_string)

    def as_config_file(self, subsection_name=None):
        ''' Return the configuration (or a subsection of it) as a  
            configuration file string.
        '''
        self._re_eval()

    def as_dict(self, subsection_name=None):
        ''' Return the configuration (or a subsection of it) as a
            Python dictionary.
        '''
        self._re_eval()
        return deepcopy(self._dictionary)

    def set_option(self, key, value=None):
        self._dictionary[key] = value
        self._re_eval

    def remove_option(self, key):
        del self._dictionary[key]
        self._re_eval

class ConfigFile(BaseConfig):
    ''' Read and parse a configuration from a local file. 
    '''
    def __init__(self, file_path):
        Configuration.__init__(self)

class ConfigQuery(BaseConfig):
    ''' Parse a configuration from a URL query string.
    '''
    def __init__(self, query_string):
        Configuration.__init__(self)
        self._re_eval(query_string=query_string)


def _test_():
    cq = ConfigQuery("WhenToTransferOutput=ON_EXIT&should_transfer_files=YES&notifications=Always")
    print cq.as_dict()
    cq.set_option('VERBOSE')
    print cq.as_query_string()
    cq.set_option('VERBOSE',False)
    print cq.as_query_string()
    cq.remove_option('VERBOSE')
    print cq.as_dict()




