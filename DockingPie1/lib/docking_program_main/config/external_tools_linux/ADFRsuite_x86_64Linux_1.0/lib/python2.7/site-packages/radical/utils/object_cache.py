
__author__    = "Radical.Utils Development Team (Andre Merzky)"
__copyright__ = "Copyright 2013, RADICAL@Rutgers"
__license__   = "MIT"


import threading
import lockable
import singleton


# ------------------------------------------------------------------------------
#
# default timeout for delayed object removal.
#
_TIMEOUT    = 10
_DEFAULT_NS = 'global'


# ------------------------------------------------------------------------------
#
@lockable.Lockable
class ObjectCache (object) :

    """ 
    This is a singleton object caching class -- it maintains a reference
    counted registry of existing objects.
    """
    
    # TODO: we should introduce namespaces -- this is a singleton, but we may want
    # to use it in several places, thus need to make sure to not use colliding
    # names...
    #
    # FIXME: this class is not thread-safe!

    __metaclass__ = singleton.Singleton

    # --------------------------------------------------------------------------
    #
    def __init__ (self, timeout=_TIMEOUT) :
        """
        Make sure the object cache dict is initialized, exactly once.

        If timeout is 0 or smaller, the objects are removed immediately --
        otherwise removal is delayed by the specified timeout in seconds, to
        avoid thrashing on frequent removal/creation.
        """

        self._timeout = timeout
        self._cache   = dict()


    # --------------------------------------------------------------------------
    #
    def get_obj (self, oid, creator, ns=_DEFAULT_NS) :
        """
        For a given object id, attempt to retrieve an existing object.  If that
        object exists, increase the reference counter, as there is now one more
        user for that object.  
        
        If that object does not exist, call the given creator, then register and
        return the object thusly created.

        oid     : id of the object to get from the cache.  
        creator : method to use to create a new object instance

                  Example:
                      def creator () :
                          return get_logger (name)

                      ret = object_cache.get_object (name, creator)
        """

        with self :

            if  not ns in self._cache :
                self._cache[ns] = dict()
            ns_cache = self._cache[ns]


            oid = str(oid)

            if  not oid in ns_cache :

                obj = creator ()

                ns_cache [oid]        = {}
                ns_cache [oid]['cnt'] = 0
                ns_cache [oid]['obj'] = obj

            ns_cache [oid]['cnt'] += 1

            return ns_cache [oid]['obj']


    # --------------------------------------------------------------------------
    #
    def rem_obj (self, obj, ns=_DEFAULT_NS) :
        """
        For a given objects instance, decrease the refcounter as the caller
        stops using that object.  Once the ref counter is '0', remove all traces
        of the object -- this should make that object eligable for Python's
        garbage collection.  Returns 'True' if the given object was indeed
        registered, 'False' otherwise.

        The removal of the object is actually time-delayed.  That way, we will
        keep the object around *just* a little longer, which provides caching
        semantics in the case of frequent creation/dstruction cycles.
        """

        with self :

            if  not ns in self._cache :
                self._cache[ns] = dict()
            ns_cache = self._cache[ns]


            for oid in ns_cache.keys () :

                if  obj == ns_cache [oid]['obj'] :

                    if  self._timeout > 0 :
                        # delay actual removeal by _timeout seconds
                        threading.Timer (self._timeout, self._rem_obj, [oid]).start ()

                    else :
                        # immediate removeal 
                        self._rem_obj (oid)
                        
                    return True

            return False  # obj not found


    # --------------------------------------------------------------------------
    #
    def _rem_obj (self, oid, ns=_DEFAULT_NS) :
        """
        actual removal of an object (identified by oid) from the cache -- see
        :func:`rem_obj()` for details.
        """

        with self :

            if  not ns in self._cache :
                self._cache[ns] = dict()
            ns_cache = self._cache[ns]


            ns_cache [oid]['cnt'] -= 1

            if  ns_cache [oid]['cnt'] == 0 :
                ns_cache [oid]['obj'] = None  # free the obj reference
                ns_cache.pop (oid, None)      # remove the cache entry


# ------------------------------------------------------------------------------

