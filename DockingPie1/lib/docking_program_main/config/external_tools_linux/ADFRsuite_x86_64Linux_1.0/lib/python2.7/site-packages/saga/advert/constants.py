
__author__    = "Andre Merzky"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import saga.namespace.constants as ns

# filesystem flags enum:
OVERWRITE      = ns.OVERWRITE       #      1
RECURSIVE      = ns.RECURSIVE       #      2
DEREFERENCE    = ns.DEREFERENCE     #      4
CREATE         = ns.CREATE          #      8
EXCLUSIVE      = ns.EXCLUSIVE       #     16
LOCK           = ns.LOCK            #     32
CREATE_PARENTS = ns.CREATE_PARENTS  #     64
TRUNCATE       =                         128
# APPEND       = reserved           #    256
READ           =                         512
WRITE          =                        1024
READ_WRITE     =                        1536
# BINARY       = reserved           #   2048


# attributes
ATTRIBUTE      = 'Attribute'
OBJECT         = 'Object'
TTL            = 'TTL'
CHANGE         = 'Change'
NEW            = 'New'
DELETE         = 'Delete'




