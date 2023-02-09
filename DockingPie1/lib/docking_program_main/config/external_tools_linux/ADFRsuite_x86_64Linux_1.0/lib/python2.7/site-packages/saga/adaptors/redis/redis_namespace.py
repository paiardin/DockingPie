
__author__    = "Andre Merzky"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


""" A collection of utilities which maps a namespace structure to redis nosql
keys.  It mirrors the Python ``os`` API, to some extent.

The redis key layout is like the following::

    # namespace tree
    /node:/                  : { type : DIR,   ... }
    /node:/etc/              : { type : DIR,   ... }
    /node:/etc/passwd        : { type : ENTRY, ... }
    /node:/etc/pass          : { type : LINK,  tgt : /etc/passwd, ... }
    ...

    # attribute storage for ns entries, aka adverts (incl. object storage)
    /data:/etc/passwd        : { _obj : <serial>, key_1 : val_1, key_2 : val_2, ... }
    ...

    # index for ls type ops
    /kids:/                  : [etc, usr, lib, ...]
    /kids:/etc/              : [init.d, rc.0, rc.1, ...
    ...

    # the two below are only enabled on 'attribute_indexing = True
    # index for attribute key lookups
    /keys:/etc/passwd:key_1  : [/vals/etc/passwd, ...]
    /keys:/etc/passwd:key_2  : [/vals/etc/passwd, ...]

    # index for attribute value lookups
    /vals:/etc/passwd:val_1  : [/keys/etc/passwd, ...]
    /vals:/etc/passwd:val_2  : [/keys/etc/passwd, ...]

all wildcard lookup versions will be slow -- only solution (iiuc) would be to
blow the indexes to cover for wildcard expansion - which is always incomplete
anyway...

All parts of the above structure are set/get in a single pipelined transactional
'MULTI' block.  One could consider to move the ops into a lua script, if
performance is insufficient.  

TODO:
    - use locks to make thread safe(r)

"""

import re
import os
import time
import string
import redis

import radical.utils         as ru
import radical.utils.logger  as rul
import radical.utils.threads as rut

import redis_cache

from   saga.exceptions       import *
from   saga.advert.constants import *

import saga.utils.misc       as sumisc


TYPE   = 'type'
OBJECT = 'object'

DIR    = 'dir'
ENTRY  = 'entry'
LINK   = 'link'

NODE   = 'node'
DATA   = 'data'
KIDS   = 'kids'
KEYS   = 'keys'
VALS   = 'vals'

MON    = 'saga-advert-events'

# --------------------------------------------------------------------
#
# for some reason, POSIX allows two leading slashes, but we need path names to
# be unique (they are used as keys...)
def redis_ns_parent (path) :
    if path == '/' or path == '//' : return '/'
    return os.path.split (path)[0]

def redis_ns_name (path) :
    if path == '/' or path == '//' : return '/'
    return os.path.split (path)[1]


# --------------------------------------------------------------------
#
class redis_ns_monitor (ru.Thread) :

    # ----------------------------------------------------------------
    #
    def __init__ (self, r, pub) :

        self.r      = r
        self.pub    = pub
        self.logger = r.logger

        self.pat = {}
        self.pat['ATTRIBUTE'] = re.compile ('\s*\[(?P<key>[^=]+)=(?P<val>.+)]\s*')

        rut.Thread.__init__ (self, self.work)
        self.setDaemon (True)


    # ----------------------------------------------------------------
    #
    def work (self) :

        try :
        
            callbacks = self.r.callbacks
            sub       = self.pub.listen ()

            while sub :

                info = sub.next ()
                data = info['data']

                if not type (data) == type ("") :
                    self.logger.warn ("ignoring event : %s"  %  data)
                    continue

                # FIXME: need proper regex parsing
                elems = data.split ()
                if not len (elems) == 3 :
                    self.logger.warn ("ignoring event args : %s"  %  data)
                    continue


                self.logger.debug ("sub %s"  %  data)
                event = elems[0]
                path  = elems[1]
                args  = elems[2:]

                if path in callbacks :
                    
                    if event == 'ATTRIBUTE' :

                        for arg in args :
                            
                            # args are formatted like '[key=val]'
                            match = self.pat[event].match (string.join (args, ' '))

                            if  not match :
                                self.logger.warn ("event parse error for %s" % args)
                                continue

                            # FIXME: error check
                            key     = match.group ('key')
                            val     = match.group ('val')

                            if key in callbacks[path] :
                                for idx in callbacks[path][key] :

                                    cb  = callbacks[path][key][idx][0]
                                    obj = callbacks[path][key][idx][1]

                                    obj.set_attribute (key, val, obj._UP)

                    if event == 'ATTRIBUTES' :
                        self.logger.warn ("unknown event type %s" % event)
                        pass

        except Exception as e :
            self.logger.critical ("redis monitoring thread crashed - disable callback handling (%s)") % str(e)
            return


# --------------------------------------------------------------------
#
class redis_ns_server (redis.Redis) :

    def __init__ (self, url) :

        if url.scheme != 'redis' :
            raise BadParameter ("scheme in url is not supported (%s != redis://...)" %  url)

        self.url        = url
        self.host       = 'localhost'
        self.port       = 6379
        self.db         = 0
        self.password   = None
        self.errors     = 'strict'

        if url.host     : self.host     = url.host
        if url.port     : self.port     = url.port
        if url.username : self.username = url.username
        if url.password : self.password = url.password

        # create redis client 
        t1 = time.time ()
        redis.Redis.__init__   (self, 
                                host      = self.host,
                                port      = self.port,
                                db        = self.db,
                                password  = self.password,
                                 errors    = self.errors)
        t2 = time.time ()

        # add a logger 
        self.logger = ru.get_logger ('radical.saga')

        # create a cache dict and attach to redis client instance.  Cache
        # lifetime is set to 10 times the redis-connect latency.
        self.cache = redis_cache.Cache (logger=self.logger, ttl=((t2-t1)*10))

        # create a second client to manage the (blocking) 
        # pubsub communication for event notifications
        self.r2 = redis.Redis  (host      = self.host,
                                port      = self.port,
                                db        = self.db,
                                password  = self.password,
                                errors    = self.errors)

        # set up pubsub endpoint, and start a thread to monitor channels
        self.callbacks = {}
        self.pub = self.r2.pubsub ()
        self.pub.subscribe (MON)
        # FIXME: create one pubsub channel per path (for paths which have
        #        callbacks registered)

        self.monitor = redis_ns_monitor (self, self.pub)
        self.monitor.start ()



    def __del__ (self) :

        if self.pub :
            self.pub.unsubscribe (MON)


# --------------------------------------------------------------------
#
class redis_ns_entry :

    # ----------------------------------------------------------------
    #
    def __init__ (self, r, path) :

        self.r         = r
        self.path      = path
        self.node      = {TYPE : None}
        self.data      = {}
        self.kids      = []
        self.valid     = False # not initialized
        self.logger    = r.logger
        self.cache     = r.cache
        self.callbacks = r.callbacks


    # ----------------------------------------------------------------
    #
    def _dump (self) :

        print "self.r         : %s" % str(self.r        )
        print "self.path      : %s" % str(self.path     )
        print "self.node      : %s" % str(self.node     )
        print "self.data      : %s" % str(self.data     )
        print "self.kids      : %s" % str(self.kids     )
        print "self.valid     : %s" % str(self.valid    )
        print "self.logger    : %s" % str(self.logger   )
        print "self.cache     : %s" % str(self.cache    )
        print "self.callbacks : %s" % str(self.callbacks)


    # ----------------------------------------------------------------
    #
    @classmethod
    def opendir (self, r, path, flags) :
    
        r.logger.debug ("redis_ns_entry.opendir %s" % path)
    
        e = redis_ns_entry (r, path)

        try :
            e.fetch () 

        except Exception as e :
    
            if  CREATE         & flags or \
                CREATE_PARENTS & flags    :
                e = redis_ns_entry (r, path)
                e.mkdir (flags)
            
            else :
                raise BadParameter ("Cannot open %s (no such directory)" % path)

    
        if not e.is_dir () :
            raise BadParameter ("Cannot open %s (not a directory)" % path)

        return e
            
    
    # ----------------------------------------------------------------
    #
    @classmethod
    def open (self, r, path, flags) :
        # FIXME: the checks below make open quite slow, as we travel down the
        #        path components.  This should be done in a single pipeline.

    
        r.logger.debug ("redis_ns_entry.open %s" % path)
    
        # make sure parent dir exists
        try :
            parent_path = redis_ns_parent (path)

            # go down the rabbit hole
            parent = self.opendir (r, parent_path, flags)


        except Exception as e :

            raise BadParameter ("Cannot open parent %s (%s)" % (parent, e))


        # try to open entry itself
        e = redis_ns_entry (r, path)
        try :
            e.fetch () 

        except Exception :
    
            if  CREATE & flags :
                e.node[TYPE] = ENTRY
                e.create (flags)
            
            else :
                raise BadParameter ("Cannot open %s (no such entry)" % path)

    
        return e
            
    
    # ----------------------------------------------------------------
    #
    def mkdir (self, flags) :
        """
        Don't call this  directly -- to create a dir, call opendir with
        'create'/'create_parents' and 'exclusive'.  If called, assumes that
        entry is invalid (and thus does not yet exist).
        """

        path = self.path

        if  self.valid :
            raise IncorrectState ("mkdir on %s fails, entry already exists" %  path)

        self.logger.debug ("redis_ns_entry.mkdir %s" % path)
    
        # if / does not exist, we always create it - no need for checks
        if path != '/' :
    
          # if CREATE_PARENTS is set, we need to check all parents, and need to create
          # them as needed.  We go top down, and terminate once a parent is found
          if CREATE_PARENTS & flags :
    
              # this will recursively travel down the chimney hole, and stop
              # whenever it finds an existing  directory
              parent = redis_ns_parent (path)
              pe     = redis_ns_entry.opendir (self.r, parent, CREATE_PARENTS)
    
          else :
          # if 'CREATE_PARENTS is not set, parent must exist.
              parent = redis_ns_parent (path)
              pe     = None
                
              try :
                  pe = redis_ns_entry (self.r, parent)
                  pe.fetch ()

              except Exception as e :
                  raise BadParameter ("mkdir %s fails, parent does not exist: %s" %  (path, e))

              if not pe.is_dir () :
                  raise BadParameter ("mkdir %s fails, parent is no directory: %s" %  (path, parent))

        self.node[TYPE] = DIR
        self.create (flags)
    

    # ----------------------------------------------------------------
    #
    def create (self, flags=0) :
        """
        This assumes that the target entry does not exist.  If flags contains
        CREATE though, we should not raise an exception if it in fact does not.
        """

        path = self.path

        if  self.valid :
            raise IncorrectState ("mkdir on %s fails, entry already exists" % path)

        # FIXME: need to ensure this via a WATCH call.
    
        self.logger.debug ("redis_ns_entry.create %s" % path)
    
        name   = redis_ns_name   (path)
        parent = redis_ns_parent (path)
    
        now    = time.time()
    
        self.node['mtime'] = now
        self.node['ctime'] = now
    
        p = self.r.pipeline ()
        # FIXME: add guard
        if len(self.node) : p.hmset  (NODE+':'+path, self.node)
        if len(self.data) : p.hmset  (DATA+':'+path, self.data)
        if len(self.kids) : p.hmset  (KIDS+':'+path, self.kids)
    
        # add entry as  kid to parent
        # FIXME: avoid duplicated entries!
        if path != '/' :
            p.sadd (KIDS+':'+parent, path)
    
        # add new index entries
        for key in self.data :
            val = self.data[key]
            p.sadd (KEYS+':'+str(key), path)
            p.sadd (VALS+':'+str(val), path)
    
    
        # FIXME: eval vals
        p.execute ()
    
        # issue notification about entry creation to parent dir
        self.logger.debug ("pub CREATE %s [%s]"  %  (parent, name))
        self.r.publish   (MON, "CREATE %s [%s]"  %  (parent, name))
    
        # refresh cache state
        self.cache.set (NODE+':'+path, self.node)
        self.cache.set (DATA+':'+path, self.data)
        self.cache.set (KIDS+':'+path, self.kids)

        self.valid = True
    

    # ----------------------------------------------------------------
    #
    def __str__ (self) :

        return "[%-4s] %-25s %s %s : %s" \
             % (self.node[TYPE], self.path, self.node, self.kids, self.data) 


    # ----------------------------------------------------------------
    #
    def is_dir (self) :
        
        if self.node[TYPE] == DIR :
            return True
        else :
            return False


    # ----------------------------------------------------------------
    #
    def list (self) :
        
        if  not self.node[TYPE] == DIR :
            raise IncorrectState ("'list()' is only supported on directories")
        
        self.fetch ()

        return self.kids


    # ----------------------------------------------------------------
    #
    def fetch (self) :

        self.logger.debug ("redis_ns_entry.fetch %s" % self.path)

        path = self.path

        try :
            self.node = self.cache.get (NODE+':'+path)
            self.data = self.cache.get (DATA+':'+path)

            if self.node[TYPE] == DIR :
                self.kids = self.cache.get (KIDS+':'+path)

            self.valid = True
            return

        except Exception :
            # some cache ops failed, so we need to properly fetch data.  We simply
            # fetch all of it
            pass


        try :
            p = self.r.pipeline ()
            p.hgetall  (NODE+':'+path)
            p.hgetall  (DATA+':'+path)
            p.smembers (KIDS+':'+path)
            values = p.execute ()
        
            if len (values) != 3 :
                self.valid = False
                return

            self.valid = True

            # FIXME: check val types
            self.node = values[0]
            self.data = values[1]
            self.kids = values[2] # will be 'None' for non-DIR entries

            if len (self.node) == 0 :
                self.valid = False
                raise IncorrectState ("backend entry seems to be gone or corrupted")

            # cache our newly found entries
            self.cache.set (NODE+':'+path, self.node)
            self.cache.set (DATA+':'+path, self.data)
            self.cache.set (KIDS+':'+path, self.kids)

            # fetched from redis ok
            self.valid = True

        except Exception as e :
            self.valid = False
            raise IncorrectState ("backend entry is gone or corrupted: %s" % str(e))



    # ----------------------------------------------------------------
    #
    def get_data (self) :

        self.logger.debug ("redis_ns_entry.get_data %s" % self.path)
    
        self.fetch () # refresh cache/state as needed

        return self.data


    # ----------------------------------------------------------------
    #
    def get_object (self) :
    
        self.logger.debug ("redis_ns_entry.get_object %s" % self.path)
    
        self.fetch () # refresh cache/state as needed

        # FIXME: dig '_object' out of the data hash, and de-serialize it to the
        # respective SAGA API object (if possible)
    
        return None


    # ----------------------------------------------------------------
    #
    def get_key (self, key) :
    
        self.logger.debug ("redis_ns_entry.get_key %s" % (key))
    
        self.fetch () # refresh cache/state as needed

        if not key in self.data :
            raise BadParameter ("no such attribute (%s)" %  key)

        val = self.data[key]

        if key == 'foo' and val == 'stop' :
            self.cache._dump ()

        return val


    # ----------------------------------------------------------------
    #
    def set_key (self, key, val) :
    
        path = self.path
        self.logger.debug ("set_key %s: %s" % (path, key))
    
        self.fetch () # refresh cache/state as needed

        # FIXME: we only fetch() for the indexes - we should optimize that again
        # by moving index consolidation into a separate thread (p.srem below)

        if key in self.data and self.data[key] == val :

            # nothing changed - so just trigger the set event
            self.logger.debug ("Pub ATTRIBUTE %s [%s=%s]"  %  (path, key, val))
            self.r.publish   (MON, "ATTRIBUTE %s [%s=%s]"  %  (path, key, val))

            # nothing else to do
            return


        # need to set the key, and update the key/val indexes
        p = self.r.pipeline ()
        # FIXME: add guard
        p.hmset  (NODE+':'+path, {'mtime': time.time()})
        p.hmset  (DATA+':'+path, {key    : val})
    
        # delete old invalid index entry
        if key in self.data :
            # we keep the key index entry around
            p.srem (VALS+':'+str(self.data[key]), path)
        else :
            # new key: add new key index entry
            p.sadd (KEYS+':'+str(key), path)
    
        # always add new value index entry
        p.sadd (VALS+':'+str(val), path)
    
        # FIXME: eval return types / values
        vals = p.execute ()
    
        # issue notification about key creation/update
        self.logger.debug ("PUB ATTRIBUTE %s [%s=%s]"  %  (path, key, val))
        self.r.publish   (MON, "ATTRIBUTE %s [%s=%s]"  %  (path, key, val))
    
        # update cache
        self.data[key] = val
        self.cache.set (DATA+':'+path, self.data)


    # ----------------------------------------------------------------
    #
    def manage_callback (self, key, id, cb, obj) :
        # FIXME: this needs a shared lock with the monitor thread!
    
        self.logger.debug ("redis_ns_entry.manage__callback %s : %s" % (self.path, key))
    
        if not self.path in self.callbacks :
            self.callbacks[self.path] = {}
            if not key in self.callbacks[self.path] :
                self.callbacks[self.path][key] = {}
    
        if id == None :
            del self.callbacks[self.path][key]
    
        if cb :
            self.callbacks[self.path][key][id] = [cb, obj]
        else :
            # cb == None: remove that callback
            # FIXME: if we have no callbacks for key, remove [key]
            # FIXME: if we have no keys for path, remove [path]
            # FIXME: if we remove a [path], unsubscribe for its notifications
            del self.callbacks[self.path][key][id]

  



