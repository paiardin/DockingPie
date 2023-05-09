
__author__    = "Radical.Utils Development Team (Andre Merzky)"
__copyright__ = "Copyright 2013, RADICAL@Rutgers"
__license__   = "MIT"


import time
import threading
import contextlib
import weakref
import singleton as rus


# ------------------------------------------------------------------------------
#
TIMEOUT   = 0.1 # sleep between lock checks
READONLY  = 'ReadOnly'
READWRITE = 'ReadWrite'


# ------------------------------------------------------------------------------
#
class _Registry (object) :
    """
    This singleton registry class maintains a set of object instances, by some
    entity ID.  Consumers can acquire and release those instances, for
    `READONLY` or `READWRITE` activities.  An unlimited number of READONLY
    activities can be active at any time, on any given instance -- but if any
    `READWRITE` acquisition will block until no other instances are leased
    anymore, and any `READONLY` acquisition will block until an eventual release
    of a `READWRITE` lease::

        my_obj = MyClass (eid='my_id')
        print my_obj.id   # prints my_id
        ...

        import raducal.utils as  ru
        ro_instance_1 = ru.Registry.acquire ('my_id', ru.READONLY)
        ro_instance_2 = ru.Registry.acquire ('my_id', ru.READONLY)
        ...

        rw_instance_1 = ru.Registry.acquire ('my_id', ru.READWRITE)

        # the last call would block until the following two calls have been
        # issued:
        ru.Registry.release ('my_id')
        ru.Registry.release ('my_id')

    The correct semantic functionality of this class depends on the careful
    release of unused instances.  Using a `with` statement is encouraged, and
    a `lease` context manager is provided for that purpose::

        with ru.Registry.lease (eid, ru.READONLY) as ro_instance :
            ro_instance.do_something ()
            ro_instance.do_nothing   ()
            ro_instance.do_something ()

        with ru.Registry.lease (eid, ru.READWRITE) as rw_instance :
            rw_instance.change_something

    The registry will thus ensure that a consumer will always see instances
    which are not changed by a third party over the scope of a lease.

    *Requirements*

    Registered object instances must fulfill two requirements:
    
    * they must be `lockable`, i.e. during `aquire` we can call `entity.lock()`
      on them to activate a `threading.RLock`, and during `release()` we can
      call `entity.unlock()` to deactivate it.  Instances will thus guaranteed
      to be locked during lease.

    * they must be `identifiable`, i.e. they must have an `id` property.
      Alternatively, an entity ID can be passed as optional parameter during
      registration -- all followup methods will require that eid.  


    *Note* (for Troy as consumer of this registry, but applicable in general):

    It is in general not a favourable design to have large, all-visible
    registries in a multi-component software stack, as ownership of state
    transitions can easily become blurry, and as a registry can also become
    a performance bottleneck on frequent queries -- so why are we doing this?

    First of all, state transitions in Troy *are* blurry to some extent, as Troy
    aims to support a variety of state diagram transitions, and thus the order
    of transitions is not pre-defined (first derive overlay then schedule CUs,
    or vice versa?).  Also, the ownership of state transitions for a workload is
    not easily scoped (the :class:`Planner` will move a workload from `DESCRIBED` 
    to `PLANNED`, the :class:`WorkloadManager` will move it from `PLANNED` to
    `TRANSLATED`, etc.  And, finally, we want to allow for re-scheduling,
    re-translation, re-planning etc, which would require us to pass control of
    a workload back and forth between different modules.  Finally, this seems to
    be useful for inspection and introspection of Troy activities related to
    specific workload instances.

    In that context, a registry seems the (much) lesser of two devils: The
    registry class will allow to register workload and overlay instances, and to
    acquire/release control over them.  The module which acquires control needs
    to ascertain that the respective workload and overlay instances are in
    a usable state for that module -- the registry is neither interpreting nor
    enforcing any state model on the managed instances -- that is up to the
    leasing module.  Neither will the registry perform garbage collection --
    active unregistration will remove instances from the registry, but not
    enforce a deletion.

    This is a singleton class.  We assume that Workload *and* Overlay IDs are
    unique.
    """
    __metaclass__ = rus.Singleton


    # --------------------------------------------------------------------------
    #
    def __init__ (self) :
        """
        Create a new Registry instance.  
        """

        # make this instance lockable
        self.lock = threading.RLock ()

        self._registry = dict()
        self._session  = None  # this will be set by Troy.submit_workload()


    # ------------------------------------------------------------------------------
    #
    @contextlib.contextmanager
    def lease (self, oid, mode=READWRITE) :
        try :
            yield self.acquire (oid, mode)
        finally :
            self.release (oid)

    # --------------------------------------------------------------------------
    #
    def register (self, entity, eid=None) :
        """
        register a new object instance.
        """

        # lock manager before checking/manipulating the registry
        with self.lock :

            # check if instance is lockable
            if  not hasattr (entity, '__enter__') or \
                not hasattr (entity, '__exit__' ) or \
                not hasattr (entity, 'lock'     ) or \
                not hasattr (entity, 'unlock'   ) or \
                not hasattr (entity, 'locked'   ) or \
                not hasattr (entity, '_locked'  ) or \
                not hasattr (entity, '_rlock'   ) :
                raise TypeError ("Registry only manages lockables")

            # check if instance is identifiable
            if  not eid :
                if  not hasattr (entity, 'id') :
                    raise TypeError ("Registry only manages identifiables")
                eid = entity.id

            if  eid in self._registry :
                raise ValueError ("'%s' is already registered" % eid)

            self._registry[eid] = {}
            self._registry[eid]['ro_leases'] = 0  # not leased
            self._registry[eid]['rw_leases'] = 0  # not leased
            self._registry[eid]['entity']    = weakref.ref (entity)


    # --------------------------------------------------------------------------
    #
    def acquire (self, eid, mode) :
        """
        temporarily relinquish control over the referenced identity to the
        caller.
        """

        # sanity check
        if  not eid in self._registry :
            raise KeyError ("'%s' is not registered" % eid)

        # wait for the entity to be fee for the expected usage
        while True :

            # lock manager before checking/manipulating the registry
            with  self.lock :

                if  mode == READONLY :
                    # make sure we have no active write lease

                    if  self._registry[eid]['rw_leases'] :
                        # need to wait
                        time.sleep (TIMEOUT)
                        continue

                    else :
                        # free for READONLY use
                        self._registry[eid]['ro_leases'] += 1
                        break

                if  mode == READWRITE :
                    # make sure we have no active read or write lease

                    if  self._registry[eid]['rw_leases'] or \
                        self._registry[eid]['ro_leases'] :
                        # need to wait
                        time.sleep (TIMEOUT)
                        continue

                    else :
                        # free for READWRITE use
                        self._registry[eid]['rw_leases'] += 1
                        break

        # acquire entity lock
        entity = self._registry[eid]['entity']()
        
        if  None == entity :
            raise KeyError ("'%s' was deallocated" % eid)
            
        else :
            # all is well...
            entity.lock ()
            return entity


    # --------------------------------------------------------------------------
    #
    def release (self, eid) :
        """
        relinquish the control over the referenced entity
        """

        # sanity check
        if  not eid in self._registry :
            raise KeyError ("'%s' is not registered" % eid)

        # lock manager before checking/manipulating the registry
        with  self.lock :

            if  not self._registry[eid]['ro_leases'] and \
                not self._registry[eid]['rw_leases'] :
                raise ValueError ("'%s' was not acquired" % eid)

            # release entity lease
            if  self._registry[eid]['ro_leases'] :
                self._registry[eid]['ro_leases'] -= 1
            elif self._registry[eid]['rw_leases'] :
                self._registry[eid]['rw_leases'] -= 1

            # release entity lock
            entity = self._registry[eid]['entity']()

            if  entity == None :
                raise KeyError ("'%s' was deallocated" % eid)

            else :
                # all is well...
                entity.unlock ()


    # --------------------------------------------------------------------------
    #
    def unregister (self, eid) :
        """
        remove the reference entity from the registry, but do not explicitly
        call the entity's destructor.  This will unlock the entity.
        """

        # sanity check
        if  not eid in self._registry :
            raise KeyError ("'%s' is not registered" % eid)

        # lock manager before checking/manipulating the registry
        with  self.lock :

            entity = self._registry[eid]['entity']()

            if  entity == None :
                raise KeyError ("'%s' was deallocated" % eid)


            # unlock entity
            while self._registry[eid]['ro_leases'] :
                entity.unlock ()
                self._registry[eid]['ro_leases'] -= 1

            while self._registry[eid]['rw_leases'] :
                entity.unlock ()
                self._registry[eid]['rw_leases'] -= 1


            # remove entity from registry, w/o a trace...
            del self._registry[eid]


# ------------------------------------------------------------------------------
#
# create a global registry instance
#
Registry = _Registry ()


# ------------------------------------------------------------------------------

