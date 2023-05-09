
__author__    = "Radical.Utils Development Team (Andre Merzky)"
__copyright__ = "Copyright 2013, RADICAL@Rutgers"
__license__   = "MIT"


import os
import imp
import sys
import glob

import singleton

from .logger import get_logger


# ------------------------------------------------------------------------------
#
class _PluginRegistry (dict) :
    """
    The plugin registry helper class avoids that plugins are loaded twice.
    """

    __metaclass__ = singleton.Singleton


    # --------------------------------------------------------------------------
    #
    def __init__ (self) :

        self._registry = dict ()


    # --------------------------------------------------------------------------
    #
    def register (self, namespace, plugins) :

        if  not namespace in self._registry :
            self._registry[namespace] = plugins


    # --------------------------------------------------------------------------
    #
    def retrieve (self, namespace) :

        if  namespace in self._registry :
            return self._registry[namespace]

        return None


# ------------------------------------------------------------------------------
#
class PluginManager (object) :
    """ 
    The RADICAL plugin management and loading utility.

    The plugin manager allows to manage plugins of a specific types.  For those
    types, the manager can search for installed plugins, list and describe
    plugins found, load plugins, and instantiate the plugin for further use.

    Example::

        # try to load the 'echo' plugin from the 'radical' namespace
        plugin_type = 'echo'

        pm = radical.utils.PluginManager ('radical')

        for plugin_name in pm.list (plugin_type) :
            print plugin_name
            print pm.describe (plugin_type, plugin_name)

        default_plugin = pm.load ('echo', 'default')

        default_plugin.init_plugin ("world")
        default_plugin.run ()  # prints "hello default world"


    The plugins are expected to follow a specific naming and coding schema to be
    recognized by the plugin manager.  The naming schema is:

        [namespace].plugins.[ptype].plugin_[ptype]_[pname].py

    i.e. for the example above: `radical.plugins.echo.plugin_echo_default.py`

    The plugin code consists of two parts:  a plugin description, and a plugin
    class.  The description is a module level dictionary named
    `PLUGIN_DESCRIPTION`, the plugin class must be named `PLUGIN_CLASS`, and
    must have a class constructor `__init__(*args, **kwargs)` to create plugin
    instances for further use.

    At this point, we leave the definition of the exact plugin signatures open,
    but expect that to be more strictly defined per plugin type in the future.

    Note that the PluginManager construction is, at this point, not considered
    thread safe.
    """


    #---------------------------------------------------------------------------
    # 
    def __init__ (self, namespace) :
        """
        namespace: name of module (plugins are expected in namespace/plugins/)
        """

        # import here to avoid circular imports
        import radical.utils.logger as logger

        self._namespace = namespace
        self._logger    = get_logger('radical.utils')
        self._registry  = _PluginRegistry ()  # singleton
        self._plugins   = self._registry.retrieve (self._namespace)

        # load adaptors if registry didn't have any registered, yet
        if  not self._plugins :
            self._load_plugins ()
            self._registry.register (self._namespace, self._plugins)


    #---------------------------------------------------------------------------
    # 
    def _load_plugins (self) :
        """ 
        Load all plugins for the given namespace.  Previously loaded plugins
        are overloaded.
        """

        # start wish a fresh plugin registry
        self._plugins = dict () 

        self._logger.info ('loading plugins for namespace %s' % self._namespace)

        # avoid to load plugins twice in case of redundant sys paths
        seen = list() 

        # search for plugins in all system module paths
        for spath in sys.path :

            # we only load plugins installed under the namespace hierarchy
            npath = self._namespace.replace ('.', '/')
            ppath = "%s/%s/plugins/"  %  (spath, npath)
            pglob = "*/plugin_*.py"  

            # make sure the 'plugins' dir exists
            if  not os.path.isdir (ppath) :
                continue

            # we assume that all python sources in that location are
            # suitable plugins
            pfiles = glob.glob (ppath + pglob)

            if  not pfiles :
                continue

            for pfile in pfiles :

                # from the full plugin file name, derive a short name for more
                # useful logging, duplication checks etc. -- simply remove the
                # namespace path portion...
                if  pfile.startswith (spath) :
                    pshort = pfile[len(spath):]
                else :
                    pshort = pfile

              # print "pfile : %s" % pfile
              # print "npath : %s" % npath
              # print "spath : %s" % spath
              # print "pshort: %s" % pshort

                # check for duplication
                if  pshort in seen :
                    continue
                else :
                    seen.append (pshort)

                # modname for 'load_source' needs to be unique, otherwise global
                # vars in the plugin file (such as, aehm, PLUGIN_DESCRIPTION)
                # will be overwritten by the next plugin load.
                pmodname = "%s.plugins.%s" % (self._namespace, 
                                              os.path.basename(os.path.dirname(pfile)))
                modname  = "%s.%s"         % (pmodname,
                                              os.path.splitext(os.path.basename(pfile))[0])

                try :
                    # load and register the plugin

                    # modname is unique and correct -- but load_source raises
                    # a RuntimeWarning, because, apparently, the plugin's parent
                    # module cannot be found.  In fact, the parent is a proper
                    # python module, and it loads fine via 'import' -- but it is
                    # is not imported before, load_source doesn't like that.
                    # Well, the parent module actually SHOULD NOT be imported --
                    # but we do it here anyways to silence the warning.  Thanks
                    # python...
                    __import__ (pmodname)

                    # now load the plugin proper
                    plugin = imp.load_source (modname, pfile)

                    # get plugin details from description
                    ptype  = plugin.PLUGIN_DESCRIPTION.get ('type',        None)
                    pname  = plugin.PLUGIN_DESCRIPTION.get ('name',        None)
                    pvers  = plugin.PLUGIN_DESCRIPTION.get ('version',     None)
                    pdescr = plugin.PLUGIN_DESCRIPTION.get ('description', None)

                    # make sure details are complete
                    if  not ptype  : 
                        self._logger.error ('no plugin type in %s' % pshort)
                        continue

                    if  not pname  : 
                        self._logger.error ('no plugin name in %s' % pshort)
                        continue

                    if  not pvers  : 
                        self._logger.error ('no plugin version in %s' % pshort)
                        continue

                    if  not pdescr : 
                        self._logger.error ('no plugin description in %s' % pshort)
                        continue

                    # now put the plugin and plugin info into the plugin
                    # registry.  Duh!
                    if  not ptype in self._plugins :
                        self._plugins[ptype] = {}

                    if  pname in self._plugins[ptype] :
                        self._logger.warn ('overloading plugin %s' % pshort)

                    self._plugins[ptype][pname] = {
                        'class'       : plugin.PLUGIN_CLASS,
                        'type'        : ptype, 
                        'name'        : pname, 
                        'version'     : pvers, 
                        'description' : pdescr,
                        'instance'    : None
                    }

                    self._logger.debug ('loading plugin %s' % pfile)
                    self._logger.info  ('loading plugin %s' % pshort)

                except Exception as e :
                    self._logger.error ('loading plugin %s failed: %s' % (pshort, e))


    #---------------------------------------------------------------------------
    # 
    def list_types (self) :
        """
        return a list of loaded plugin types
        """
        return self._plugins.keys ()


    #---------------------------------------------------------------------------
    # 
    def list (self, ptype) :
        """
        return a list of loaded plugins for a given plugin type
        """
        if  not ptype in self._plugins :
            self._logger.debug (self.dump_str())
            raise LookupError ("No such plugin type %s in %s" \
                    % (ptype, self._plugins.keys()))

        return self._plugins[ptype].keys ()


    #---------------------------------------------------------------------------
    # 
    def describe (self, ptype, pname) :
        """
        return a plugin details for a given plugin type / name pair
        """
        if  not ptype in self._plugins :
            self._logger.debug (self.dump_str())
            raise LookupError ("No such plugin type %s in %s" \
                    % (ptype, self._plugins.keys()))

        if  not pname in self._plugins[ptype] :
            self._logger.debug (self.dump_str())
            raise LookupError ("No such plugin name %s (type: %s) in %s" \
                    % (pname, ptype, self._plugins[ptype].keys()))

        return self._plugins[ptype][pname]


    #---------------------------------------------------------------------------
    # 
    def load (self, ptype, pname) :
        """
        check if a plugin with given type and name was loaded, if so, instantiate its
        plugin class and return it.
        """

        if  not ptype in self._plugins :
            self._logger.debug (self.dump_str())
            raise LookupError ("No such plugin type %s in %s" \
                    % (ptype, self._plugins.keys()))

        if  not pname in self._plugins[ptype] :
            self._logger.debug (self.dump_str())
            raise LookupError ("No such plugin name %s (type: %s) in %s" \
                    % (pname, ptype, self._plugins[ptype].keys()))

        # create new plugin instance
        return self._plugins[ptype][pname]['class']()


    #---------------------------------------------------------------------------
    # 
    def dump (self) :

        import pprint
        pprint.pprint (self._plugins)


    #---------------------------------------------------------------------------
    # 
    def dump_str (self) :

        import pprint
        return "\n%s" % pprint.pformat (self._plugins)


# ------------------------------------------------------------------------------

