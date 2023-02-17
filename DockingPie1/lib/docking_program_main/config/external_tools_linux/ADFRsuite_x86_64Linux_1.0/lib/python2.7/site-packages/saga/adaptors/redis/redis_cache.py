
__author__    = "Andre Merzky"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import time

import threading
import redis_ordered_dict    as rod

CACHE_DEFAULT_SIZE = 10000
CACHE_DEFAULT_TTL  = 1.0    # 1 second

VAL = 'val'
TTL = 'ttl'

######################################################################
#
class Cache :

    # ----------------------------------------------------------------
    #
    def __init__ (self, logger, size=CACHE_DEFAULT_SIZE, ttl=CACHE_DEFAULT_TTL) :

        if int (size) < 1 :
            raise AttributeError ('size < 1 or not a number')

        if int (ttl) < 0 :
            raise AttributeError ('ttl < 0 or not a number')

        self.size   = size
        self.ttl    = ttl
        self.dict   = rod.OrderedDict ()
        self.lock   = threading.RLock ()
        self.logger = logger
        self.hit    = 0
        self.miss   = 0

        # start a thread which, with low priority, cleans out the dict now and
        # then (pops items until a live one is found


    # ----------------------------------------------------------------
    #
    def _dump (self) :
        print " ---------------------------------------------- "
        print " CACHE STATISTICS : "
        print " size: %5d" % len(self.dict)
        print " hit : %5d" % self.hit 
        print " miss: %5d" % self.miss
        print self.dict.keys()
        print " ---------------------------------------------- "


    # ----------------------------------------------------------------
    #
    def get (self, key) :

        self.logger.debug ("redis_cache_get %s", key)

        with self.lock:

            # check if we have a live entry
            if key in self.dict :

                now = time.time ()

                if self.ttl and self.dict[key][TTL] > now :
                    # if yes, cache hit!
                    # return data -- doh!
                    self.hit += 1
                    return self.dict[key][VAL]

                else :
                    # entry timed out
                    del self.dict[key]

            # cache entry not found, or timed out
            self.miss += 1
            raise AttributeError ("cache miss for '%s' " % key)


    # ----------------------------------------------------------------
    #
    def set (self, key, value) :

        with self.lock :

            # remove superfluous(?) entries
            while len (self.dict) >= self.size :
                self.dict.popitem (last=False)

            self.dict[key]      = {}
            self.dict[key][VAL] = value
            self.dict[key][TTL] = time.time () + self.ttl


    # ----------------------------------------------------------------
    #
    def delete (self, key) :

        with self.lock :
            del self.dict[key]


    # ----------------------------------------------------------------




