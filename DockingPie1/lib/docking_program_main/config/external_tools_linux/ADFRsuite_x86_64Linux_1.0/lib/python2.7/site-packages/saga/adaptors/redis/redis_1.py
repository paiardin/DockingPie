
__author__    = "Andre Merzky"
__copyright__ = "Copyright 2013, The SAGA Project"
__license__   = "MIT"


import redis

r = redis.Redis (host='localhost', password='securedis')

print "------------------ set / get ------------------"
print r.set ("name", "DeGizmo")
print r.get ("name")

print "------------------ set / inc / decr / get -----"
print r.set  ("hit_counter", 1)
print r.incr ("hit_counter")
print r.get  ("hit_counter")
print r.decr ("hit_counter")
print r.get  ("hit_counter")

print "------------------ rpush / lrange / ... -------"
print r.rpush  ("members", "Adam")
print r.rpush  ("members", "Bob")
print r.rpush  ("members", "Carol")
print r.lrange ("members", 0, -1)
print r.llen   ("members")
print r.lindex ("members", 1)

print "------------------ dict set -------------------"

print r.hmset ("key1", {'11' : 'ONE', '1' : 'one'})
print r.hmset ("key2", {'22' : 'TWO', '2' : 'two'})
print r.hmset ("key3", {'33' : 'TRE', '3' : 'tre'})

print "------------------ pipeline dict get ----------"
pipe = r.pipeline()

for key in ['key1', 'key2', 'key3'] :
    pipe.hgetall (key)

for val in pipe.execute ():
    print val

print "------------------ list keys ------------------"
print r.keys ("*")


print "------------------ lua script -----------------"
lua = """
local value  = redis.call ('GET', KEYS[1])
value        = tonumber   (value)
return value * ARGV[1]"""

multiply = r.register_script (lua)
print r.set ('foo', 2)
print multiply (keys=['foo'], args=[5])

print "-----------------------------------------------"

