
__author__    = "Andre Merzky, Ole Weidner"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import radical.utils.signatures as rus

import saga.attributes       as sa
import saga.exceptions       as se
import constants             as const

#-------------------------------------------------------------------------------
#
class Description (sa.Attributes) :
    """ 
    The resource description class. 

    Resource descriptions are used for two purposes:

      * an application can pass a description instances to a
        :class:`saga.resource.Manager` instance, to request control 
        over the resource slice described in the description; 

      * an application can request a resource's description for 
        inspection of resource properties.

    There are three specific types of descriptions: 
      
      * :class:`saga.resource.ComputeDescription` for the description of
        resources with compute capabilities;

      * :class:`saga.resource.StorageDescription` for the description of
        resources with data storage capabilities;

      * :class:`saga.resource.NetworkDescription` for the description of
        resources with communication capabilities.

    There is at this point no notion of resources which combine different
    capabilities.

    For all these capabilities, the following attributes are supported:

      * `RType`      : `Enum`, describing the capabilities of the resource
                       (`COMPUTE`, `STORAGE` or `NETWORK`)

      * `Template`   : `String`, a backend specific resource class with some
                       pre-defined hardware properties to apply to the resource.

      * `Image`      : `String`, a backend specific resource class with some 
                       pre-defined software properties to apply to the resource.

      * `Dynamic`    : `Boolean, if `True` signifies that the resource may
                       dynamically change its properties at runtime

      * `Start`      : `Integer (seconds) since epoch when the resource is 
                       expected to enter / when the resource entered `ACTIVE` 
                       state.

      * `End`        : `Integer (seconds) since epoch when the resource is 
                       expected to enter / when the resource entered a `FINAL` 
                       state.

      * `Duration`   : `Integer`, seconds for which the resource is expected to
                       remain / the resource remained in `ACTIVE` state.

      * `MachineOS`  : `String`, for `COMPUTE` resources, specifies the
                       operating system type running on that resource.

      * `MachineArch : `String`, for `COMPUTE` resources, specifies the
                       machine architecture of that resource.

      * `Size`       : `Integer`, for `COMPUTE` resources, specifies the
                       number of process slots provided, for `STORAGE` resource
                       specifies the number of bytes, of the resource.

      * `Memory`     : `Integer`, for `COMPUTE` resources, specifies the
                       number of bytes provided as memory.

      * `Access`     : `String`, usually an URL, which specifies the contact
                       point for the resource capability interface / service
                       interface.
    """

    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Description', 
                  rus.optional (dict))
    @rus.returns (rus.nothing)
    def __init__ (self, d=None):
        """
        __init__()

        Create a new Description instance.
        """

        # set attribute interface properties

        self._attributes_extensible  (False)
        self._attributes_camelcasing (True)

        # register properties with the attribute interface

        self._attributes_register  (const.RTYPE       , None , sa.ENUM  , sa.SCALAR, sa.WRITEABLE) 
        self._attributes_register  (const.TEMPLATE    , None , sa.STRING, sa.SCALAR, sa.WRITEABLE) 
        self._attributes_register  (const.IMAGE       , None , sa.STRING, sa.SCALAR, sa.WRITEABLE) 
        self._attributes_register  (const.DYNAMIC     , False, sa.BOOL  , sa.SCALAR, sa.WRITEABLE) 
        self._attributes_register  (const.START       , None , sa.TIME  , sa.SCALAR, sa.WRITEABLE) 
        self._attributes_register  (const.END         , None , sa.TIME  , sa.SCALAR, sa.WRITEABLE) 
        self._attributes_register  (const.DURATION    , None , sa.TIME  , sa.SCALAR, sa.WRITEABLE) 
        self._attributes_register  (const.MACHINE_OS  , None , sa.STRING, sa.SCALAR, sa.WRITEABLE) 
        self._attributes_register  (const.MACHINE_ARCH, None , sa.STRING, sa.SCALAR, sa.WRITEABLE) 
        self._attributes_register  (const.SIZE        , 1    , sa.INT   , sa.SCALAR, sa.WRITEABLE) 
        self._attributes_register  (const.MEMORY      , None , sa.STRING, sa.SCALAR, sa.WRITEABLE) 
        self._attributes_register  (const.ACCESS      , None , sa.STRING, sa.SCALAR, sa.WRITEABLE) 

        self._attributes_set_enums (const.RTYPE, [const.COMPUTE ,
                                                  const.STORAGE ,
                                                  const.NETWORK ])

        # FIXME: initialization should be done in Attributes: initialization
        # from dict or from other attributable
        #
        if  d :
            for key in d.list_attributes () :
                self.set_attribute (key, d.get_attribute (key)) 


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Description', 
                  ('Description', dict))
    @rus.returns ('Description')
    def __deepcopy__ (self, memo) :
        """
        An alias for `clone()`.
        """
        return self.clone()


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Description', 
                  rus.optional ('Description'))
    @rus.returns ('Description')
    def clone (self, other=None) :
        """ 
        clone()

        Implements deep copy. 

        Unlike the default python assignment (copy object reference),
        a deep copy will create a new object instance with the same state --
        after a deep copy, a change on one instance will not affect the other.
        """

        # a job description only has attributes - so create a new instance,
        # clone the attribs, and done.
        if not other :
            other = Description ()

        return self._attributes_deep_copy (other)


# ------------------------------------------------------------------------------
#
class ComputeDescription (Description) : 
    """
    A `ComputeDescription` is a specific description for a resource which
    provides compute capabilities, i.e. the ability to run application
    executables / code.
    """

    # --------------------------------------------------------------------------
    #
    @rus.takes   ('ComputeDescription', 
                  rus.optional (dict))
    @rus.returns (rus.nothing)
    def __init__ (self, d=None) :
        
        if  d :
            if  const.RTYPE in d and d[const.RTYPE] != const.COMPUTE :
                raise se.BadParameter ("Cannot create ComputeResource with type '%s'" \
                                    % d[const.RTYPE])

        self._descr = super  (ComputeDescription, self)
        self._descr.__init__ (d)

        self.rtype = const.COMPUTE


# ------------------------------------------------------------------------------
#
class StorageDescription (Description) :
    """
    A `StorageDescription` is a specific description for a resource which
    provides storage capabilities, i.e. the ability to persistently store,
    organize and retrieve data files.
    """

    # --------------------------------------------------------------------------
    #
    @rus.takes   ('StorageDescription', 
                  rus.optional (dict))
    @rus.returns (rus.nothing)
    def __init__ (self, d=None) :
        
        if  d :
            if  const.RTYPE in d and d[const.RTYPE] != const.STORAGE :
                raise se.BadParameter ("Cannot create StorageResource with type '%s'" \
                                    % d[const.RTYPE])

        self._descr = super  (StorageDescription, self)
        self._descr.__init__ (d)

        self.rtype = const.STORAGE


# ------------------------------------------------------------------------------
#
class NetworkDescription (Description) :
    """
    A `NetworkDescription` is a specific description for a resource which
    provides network capabilities.
    """

    # --------------------------------------------------------------------------
    #
    @rus.takes   ('NetworkDescription', 
                  rus.optional (dict))
    @rus.returns (rus.nothing)
    def __init__ (self, d=None) :
        
        if  d :
            if  const.RTYPE in d and d[const.RTYPE] != const.NETWORK :
                raise se.BadParameter ("Cannot create NetworkResource with type '%s'" \
                                    % d[const.RTYPE])

        self._descr = super  (NetworkDescription, self)
        self._descr.__init__ ()

        self.rtype = const.NETWORK





