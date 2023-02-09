
__author__    = "Andre Merzky"
__copyright__ = "Copyright 2013, The SAGA Project"
__license__   = "MIT"


import redis
from   pprint import pprint

r = redis.Redis (host='localhost', password='securedis')

keep = True
keep = False

print "-----------------------------------------------"

keys = r.keys ("*")
keys.sort()

for k in keys :
    t = r.type (k)

    if  t == 'hash' :
        d = r.hgetall (k)
        if not 'url' in d :
            if not keep: r.delete (k)
        print "%-35s [%-6s] %s " % (k, t, d)
        if not keep : r.delete (k)

    elif t == 'list' :
        print "%-35s [%-6s] %s " % (k, t, r.lrange (k, 0, -1))
        if not keep: r.delete (k)

    elif t == 'set' :
        print "%-35s [%-6s] %s " % (k, t, r.smembers (k))
        if not keep: r.delete (k)

    else :
        print "%-35s [%-6s] %s " % (k, t, r.get (k))
        if not keep: r.delete (k)


print "-----------------------------------------------"

