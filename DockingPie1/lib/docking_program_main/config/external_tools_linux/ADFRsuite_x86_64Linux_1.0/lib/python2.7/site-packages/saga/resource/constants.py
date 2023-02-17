
__author__    = "Andre Merzky"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


# FIXME: OS enums, ARCH enums

# resource type enum
COMPUTE      =  1;            """ resource accepting jobs """
STORAGE      =  2;            """ storage resource (duh) """
NETWORK      =  4;            """ connects compute and storage resources """

# resource state enum """
UNKNOWN      =  None;         """ wut? """
NEW          =  1;            """ requsted, not accepting jobs, yet;
                                  initial state """
PENDING      =  2;            """ accepting jobs, will become active eventually """
ACTIVE       =  4;            """ accepting jobs, jobs can run """
CANCELED     =  8;            """ released by user; final state """
EXPIRED      = 16;            """ released by system; final state """
DONE         = EXPIRED;       """ alias """
FAILED       = 32;            """ released unexpectedly by system or internally;
                                  final state """
FINAL        = CANCELED | DONE | FAILED

# resource attributes """
ID           = 'Id';          """ url identifying a resource instance """
RTYPE        = 'Rtype';       """ type enum, identifying the resource type """
STATE        = 'State';       """ state enum, identifying the rsource state """
STATE_DETAIL = 'StateDetail'; """ string, representing the native backend state """
MANAGER      = 'Manager';     """ url, pointing to the resource's manager """
DESCRIPTION  = 'Description'; """ dict, containing resource attributes  """

# generic resource description attributes """
TEMPLATE     = 'Template';    """ string, template to which the resource 
                                  was created"""
IMAGE        = 'Image';       """ FIXME: """

# resource lifetime attributes """
DYNAMIC      = 'Dynamic';     """ bool, enables/disables on-demand 
                                  resource 
                                  resizing """
START        = 'Start';       """ time, expected time at which resource 
                                  becomes ACTIVE """
END          = 'End';         """ time, expected time at which resource 
                                  will EXPIRE """
DURATION     = 'Duration';    """ time, expected time span between ACTIVE 
                                  and EXPIRED """

# resource type specific (non-generic) attributes """
MACHINE_OS   = 'MachineOS';   """ enum, identifying the resource's operating 
                                  system """
MACHINE_ARCH = 'MachineArch'; """ enum, identifying the machine architecture """
SIZE         = 'Size';        """ int, identifying No. of slots / size in 
                                  MB / No. of IPs """
MEMORY       = 'Memory';      """ int, identifying memory size in MegaByte """
ACCESS       = 'Access';      """ string, identifying the hostname/ip, mount
                                  point or provisioning URL to access the
                                  resource """



