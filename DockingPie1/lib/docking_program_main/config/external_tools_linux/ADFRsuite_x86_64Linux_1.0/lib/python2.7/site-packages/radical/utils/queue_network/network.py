
""" 
Provide an abstraction for distributed queue network.  

The network nodes are instances of queue_network::Node.  A node can exist either
as thread, as a local process, or as a remote process, the latter spawned via
SAGA.

The communication between nodes is always via queues, were
groups of nodes feed into the same queues, and/or are fed from the same queues.
A node can be fed from multiple queues, where queues are prioritized according
to some policy.  Depending if two communicating nodes are co-threads,
co-processes, or remote processes, different queue implementations (with the
same semantic) are used.

A network is created according to a specific state model.  That model describes
the state of entities passing through the queue network.  Each active state is
enacted by a network node, and each pending state is represented by the entity
being passed through a queue.  A queue connection thus exists for every valid
state transition.

Nodes also support the observer patter, ie. they can notify observers about
enacted state transitions. 
"""

# ------------------------------------------------------------------------------
#
# http://stackoverflow.com/questions/9539052/python-dynamically-changing-base-classes-at-runtime-how-to
#
# Depending on agent architecture (which is specific to the resource type it
# runs on) can switch between different component types: using threaded (when
# running on the same node), multiprocessing (also for running on the same node,
# but avoiding python's threading problems, for the prices of slower queues),
# and remote processes (for running components on different nodes, using zeromq
# queues for communication).
#
# We do some trickery to keep the actual components independent from the actual
# schema:
#
#   - we wrap the different queue types into a rpu.Queue object
#   - we change the base class of the component dynamically to the respective type
#
# This requires components to adhere to the following restrictions:
#
#   - *only* communicate over queues -- no shared data with other components or
#     component instances.  Note that this also holds for example for the
#     scheduler!
#   - no shared data between the component class and it's run() method.  That
#     includes no sharing of queues.
#   - components inherit from base_component, and the constructor needs to
#     register all required component-internal and -external queues with that
#     base class -- the run() method can then transparently retrieve them from
#     there.
#


# ------------------------------------------------------------------------------
#
class network (object):

    # --------------------------------------------------------------------------
    #
    def __init__ (self, state_model={}, network_description={}):
        """
        span the network according to state model and network description
        """

        self.state_model = state_model


    # --------------------------------------------------------------------------
    #
    def feed (self, entities):
        """
        feed entities into the network
        """

        if not isinstance(entities, list):
            entities = [entities]

        for entity in entities:
            assert (entity.state == self.state_model['initial'])

        self.feeder_queue.push (entities)


    # --------------------------------------------------------------------------
    #
    def subscribe (self, states, callback):
        """
        subscribe for callback notifications for state transitions *into* the
        specified states.  If callbacks return 'False' or 'None', they are
        unregistered.
        """

        if not isinstance(states, list):
            states = [states]




