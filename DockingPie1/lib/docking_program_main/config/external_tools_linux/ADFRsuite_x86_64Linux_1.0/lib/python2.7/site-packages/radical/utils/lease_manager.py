
__author__    = "Radical.Utils Development Team (Andre Merzky)"
__copyright__ = "Copyright 2013, RADICAL@Rutgers"
__license__   = "MIT"


import os
import time
import lockable
import singleton
import threading

from .logger import get_logger

# default settings for lease manager
MAX_POOL_SIZE = 15     # unlimited
MAX_POOL_WAIT = 60     # seconds
MAX_OBJ_AGE   = 60*10  # 10 minutes


# ------------------------------------------------------------------------------
#
class _LeaseObject (object) :

    _uid = 0

    # --------------------------------------------------------------------------
    def __init__ (self, lm, logger, creator, args) :

        self.lm         = lm
        self.used       = False
        self.log        = logger
        self.obj        = creator (*args)
        self.uid        = 'lo.%04d' % _LeaseObject._uid
        self.t_created  = time.time()
        self.t_leased   = None
        self.t_released = time.time() # we take control *now*

        _LeaseObject._uid += 1

    # --------------------------------------------------------------------------
    def __cmp__ (self, other) :

        return (other == self.obj)

    # --------------------------------------------------------------------------
    def __enter__ (self) :

        return self.obj

    # --------------------------------------------------------------------------
    def __exit__ (self, *args) :

        self.lm.release (self)

    # --------------------------------------------------------------------------
    def lease (self, *args) :

        if  self.used :
            raise RuntimeError ("LeaseObject is already leased: %s" % self)

        self.used     = True
        self.t_leased = time.time()

      # self.log.debug ("%s was unused for %6.1fs" % (self.uid, self.t_leased-self.t_released))

    # --------------------------------------------------------------------------
    def release (self, *args) :

        if  not self.used :
            raise RuntimeError ("LeaseObject is not leased: %s" % self)

        self.used       = False
        self.t_released = time.time ()

      # self.log.debug ("%s was leased for %6.1fs" % (self.uid, self.t_released-self.t_leased))

    # --------------------------------------------------------------------------
    def is_leased (self, *args) :

        return self.used
        
# ------------------------------------------------------------------------------
#
@lockable.Lockable
class LeaseManager (object) :

    """ 
    This is a lease manager -- it creates resource instances on demand and hands
    out leases to them.  If for a given ID no object instance exists, one is
    created, locked and a lease is returned.  If for an ID an object instance
    exists and is not leased, it is locked and returned.  If one or more
    instances exist but all are in use already (leased), a new instance is
    created (up to MAX_POOL_SIZE -- can be overwritten in the lease call).
    If that limit is reached, no objects are returned, and instead the lease
    call blocks until one of the existing objects gets released.
    """
    

    # --------------------------------------------------------------------------
    #
    def __init__ (self, max_pool_size=MAX_POOL_SIZE, max_pool_wait=MAX_POOL_WAIT, max_obj_age=MAX_OBJ_AGE) :
        """
        Make sure the object dict is initialized, exactly once.
        """

        import radical.utils.logger as rul

        self._log = get_logger('radical.utils')
      # self._log.setLevel ('DEBUG')

        self._log.debug ('lm new manager')
        self._pools = dict()
        self._max_pool_size = max_pool_size
        self._max_pool_wait = max_pool_wait 
        self._max_obj_age   = max_obj_age


    # --------------------------------------------------------------------------
    #
    def _initialize_pool (self, pool_id) :
        """
        set up a new pool, but do not create any instances, yet.
        """

        with self :

          # self._log.debug ('lm check   pool (%s)' % self._pools.keys())

            if  not pool_id in self._pools :

                self._log.debug ('lm create  pool   for %s (%s) (%s)' \
                        % (pool_id, type(pool_id), self))

                self._pools[pool_id] = dict()
                self._pools[pool_id]['objects'] = list()
                self._pools[pool_id]['freed']   = None
                self._pools[pool_id]['event']   = threading.Event ()

              # import pprint
              # pprint.pprint (self._pools)

            return self._pools[pool_id]


    # --------------------------------------------------------------------------
    #
    def _create_object (self, pool_id, creator, args) :
        """
        a new instance is needed -- create one, unless max_pool_size is reached.
        If that is the case, return `None`, otherwise return the created object
        (which is locked before return).
        """

        with self :

            self._log.debug ('lm create  object for %s' % pool_id)

            if  pool_id not in self._pools :
                raise RuntimeError ("internal error: no pool for '%s'!" % pool_id)

            pool = self._pools[pool_id]

            # check if a poolsize cap is set and reched
            if  self._max_pool_size >  0 and \
                self._max_pool_size <= len(pool['objects']) :

                # no more space...
                return None

            # poolsize cap not reached -- increase pool.  If creating a new
            # object does not work for any reason, return None.
            obj = None
            try :
                obj = _LeaseObject (self, self._log, creator, args)
                obj.lease ()
                pool['objects'].append (obj)

            except Exception as e :
                # this exception needs to fall through -- we can't wait
                # for object creation problems to fix themself over time...
                self._log.exception ("Could not create lease object")
                raise

            return obj



    # --------------------------------------------------------------------------
    #
    def _remove_object (self, pool_id, obj) :
        """
        a new instance is needed -- create one, unless max_pool_size is reached.
        If that is the case, return `None`, otherwise return the created object
        (which is locked before return).
        """

        with self :

          # self._log.debug ('lm remove  object for %s' % pool_id)

            if  pool_id not in self._pools :
                raise RuntimeError ("internal error: no pool for '%s'!" % pool_id)

            pool = self._pools[pool_id]

            # poolsize cap not reached -- increase pool.  If creating a new
            # object does not work for any reason, return None.
            obj = None
            try :
                pool['objects'].remove (obj)

            except Exception as e :
                self._log.exception ("Could not remove lease object")

            return obj



    # --------------------------------------------------------------------------
    #
    def lease (self, pool_id, creator, args=[]) :
        """

        For a given object identified, attempt to retrieve an existing object
        from the pool.  If such a free (released) object is found, lock and
        refturn it.  If such objects exist but all are in use, create a new one
        up to max_pool_size (default: 10).  
        used, block untill it is freed.  If that object does not exist, create
        it and proceed per above.
        return the object thusly created.

        pool_id       : id of the pool to lease from.  The pool ID essentially
                        determines the scope of validity for the managed objects
                        -- they form a namespace, where objects are shared when
                        living under the same name space entry (the pool_id).

        creator       : method to use to create a new object instance

                        Example:
                            def creator () :
                                return get_logger (name)

                            ret = lease_manager.lease (name, creator)
        """

        pool_id = str(pool_id)

        if  not isinstance (args, list) :
            args = [args]

        with self :

            # make sure the pool exists
            pool = self._initialize_pool (pool_id)

          # self._log.debug ('lm lease   object for %s (%s)' \
          #               % (pool_id, len(pool['objects'])))
          # self._log.debug (pool['objects'])

            # find an unlocked object instance in the pool
            # NOTE: we iterate over a copy of the list, as an eventual object 
            # removeal would screw up the index...
            for obj in pool['objects'][:] :

              # self._log.debug ('lm lease   object %s use: %s' % (obj, obj.is_leased()))

                if  not obj.is_leased () :

                    # check age
                    age = time.time() - obj.t_created
                    if  age > self._max_obj_age :
                        # too old -- remove and continue to search for a younger unleased object
                      # self._log.debug ('lm retire  object %s (%6.2fs)' % (obj, age))
                        self._remove_object (pool_id, obj)
                        continue

                    # found one -- lease/lock and return it
                  # self._log.debug ('lm lease   object %s use: ok!' % obj)
                    obj.lease ()
                    return obj

            # no unlocked object found -- create a new one 
            obj = self._create_object (pool_id, creator, args)

            # FIXME: we could try_catch the above error, and then check if the
            #        pool has anything worth to wait on.  That might be useful
            #        for creation errors which are transient.  Alas, we don't
            #        have any means to distinguish them from non-transient
            #        errors, so we don't do that at this point...

            # check if we got an object
            if  obj is not None :

                # we got a locked object -- return
              # self._log.debug ('lm lease   object %s use: %s -- new!' \
              #               % (obj, obj.is_leased()))
                return obj


        # pool is full, nothing is free -- we need to wait for an event on
        # the pool to free an instance.
        # Now, this is where deadlocks will happen: any application leasing
        # too many instances in the same thread (with 'too many' meaning
        # more than max_pool_size) will lock that thread here, thus having no
        # chance to release other instances.  We thus will print a log error
        # here, and will raise a timeout exception after MAX_POOL_WAIT
        # seconds.
        # Not that we release our lock here, to give other threads the
        # chance to release objects.
        self._log.warn ('lm lease   object: pool is full')
        timer_start = time.time()
        timer_now   = time.time()

        while (timer_now - timer_start) < self._max_pool_wait :

            # wait for any release activity on the pool
            timer_left = self._max_pool_wait - (timer_now - timer_start)
            pool['event'].wait (timer_left)

            with self :

                if  pool['event'].is_set () :

                    # make sure we don't lock up
                    pool['event'].clear ()

                    # we won the race!  now get the freed object
                    obj = pool['freed']
                    pool['freed'] = None

                    if  obj is None :

                        # object was deleted -- we should have space to create
                        # a new one!
                        obj = self._create_object (pool_id, creator, args)

                        # check if we got an object
                        if  obj is not None :

                            # we got a locked object -- return
                            return obj

                        else :

                            # find did not find a freed object, did not get to
                            # create a new one -- is there a free one in the
                            # pool?
                            for obj in pool['objects'] :

                                if  not obj.is_leased () :

                                    # found one -- lease/lock and return it
                                    obj.lease ()
                                    return obj

                            # apparently not - so we give up for now...


                    else :
                        # we got a freed object -- lock and return
                        obj.lease ()
                        return obj


                # none free, none created - or we lost the rase for handling the 
                # event.  Wait again.
                timer_now = time.time ()


        # at this point we give up: we can't create a new object, can't find
        # a free one, and we are running out of wait time...
        raise LookupError ("stop waiting on object lease")


    # --------------------------------------------------------------------------
    #
    def release (self, instance, delete=False) :
        """
        the given object is not needed right now -- unlock it so that somebody
        else can lease it.  This will not delete the object,
        """

      # self._log.debug ('lm release object')

        with self :

          # import pprint
          # pprint.pprint (self._pools)

            for pool_id in self._pools :

                for obj in self._pools [pool_id]['objects'] :
                    if  instance is not obj :
                        # this is not the object you are looking for.
                        continue

                  # self._log.debug ('lm release object %s for %s' % (obj, pool_id))

                    # remove the lease lock on the object
                    obj.release ()

                    if  delete :
                        # remove the object from the pool (decreasing its 
                        # ref counter and thus making it eligible for 
                        # garbage collection).  
                        self._pools [pool_id]['objects'].remove (obj)
                        self._pools [pool_id]['freed'] = None

                    else :
                        # mark the object as freed for lease.  
                        self._pools [pool_id]['freed'] = obj

                    # notify waiting threads about the lease or creation
                    # opportunity.
                    self._pools [pool_id]['event'].set ()

                    # object has been released
                    return 

            # FIXME: log warning
            # for now we ignore double-frees
            # raise RuntimeError ("cannot release object -- not managed")


# ------------------------------------------------------------------------------

