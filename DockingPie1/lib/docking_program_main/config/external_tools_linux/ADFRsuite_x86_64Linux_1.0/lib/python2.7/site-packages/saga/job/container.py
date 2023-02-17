
__author__    = "Andre Merzky, Ole Weidner"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import saga.task as st


# ------------------------------------------------------------------------------
#
# 'forward' declaration of job.Container
#
class Container (st.Container):
    """ :todo: document me

        .. py:attribute:: jobs

           The (read-only) `jobs` property returns a list of all job objects in the container.

           :rtype: :class:`saga.job.Job` list

        .. py:attribute:: size

           The (read-only) `size` property returns the number of job objectis in the container.

           :rtype: int

        .. py:attribute:: states

           The (read-only) `states` property returns a list of states that represent the 
           states of the individual jobs in the container. 

           :rtype: list
    """

    # --------------------------------------------------------------------------
    #
    def __init__ (self) :

        self._task_container = super  (Container, self)
        self._task_container.__init__ ()

        import saga.attributes as sa

        self._attributes_register   ("Jobs",   [], sa.ANY, sa.VECTOR, sa.READONLY)
        self._attributes_set_getter ("Jobs",   self.get_jobs)


    # --------------------------------------------------------------------------
    #
    def get_jobs (self) :
        """ This is similar to get_tasks(), but returns only Job typed entries
        from the container.
        """

        tasks = self.get_tasks ()
        jobs  = []

        import saga.job  as sjob

        for task in tasks :
            if isinstance (task, sjob.job.Job) :
                jobs.append (task)

        return jobs




