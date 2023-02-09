
__author__    = "SAGA Development Team"
__copyright__ = "Copyright 2013, RADICAL"
__license__   = "MIT"


import os
import radical.utils        as ru
import radical.utils.logger as rul


# ------------------------------------------------------------------------------
#
from   saga.constants      import *

from   saga.task           import Task, Container
from   saga.attributes     import Attributes, Callback
from   saga.session        import Session, DefaultSession
from   saga.context        import Context
from   saga.url            import Url

from   saga.exceptions     import SagaException
from   saga.exceptions     import NotImplemented
from   saga.exceptions     import IncorrectURL
from   saga.exceptions     import BadParameter
from   saga.exceptions     import AlreadyExists
from   saga.exceptions     import DoesNotExist
from   saga.exceptions     import IncorrectState
from   saga.exceptions     import PermissionDenied
from   saga.exceptions     import AuthorizationFailed
from   saga.exceptions     import AuthenticationFailed
from   saga.exceptions     import Timeout
from   saga.exceptions     import NoSuccess

import saga.job
import saga.filesystem
import saga.replica
import saga.advert
import saga.resource


# ------------------------------------------------------------------------------

pwd     = os.path.dirname (__file__)
root    = "%s/.." % pwd
version, version_detail, version_branch, sdist_name, sdist_path = ru.get_version ([root, pwd])

# FIXME: the logger init will require a 'classical' ini based config, which is
# different from the json based config we use now.   May need updating once the
# radical configuration system has changed to json
_logger = ru.get_logger('radical.saga')
_logger.info ('radical.saga         version: %s' % version_detail)


# ------------------------------------------------------------------------------

