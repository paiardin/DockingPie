

import re
import os
import sys
import subprocess    as sp


# ------------------------------------------------------------------------------
#
# versioning mechanism:
#
#   - short_version:  1.2.3                   - is used for installation
#   - long_version:   1.2.3-9-g0684b06-devel  - is used as runtime (ru.version)
#   - both are derived from the last git tag and branch information
#   - VERSION files are created on install, with the long_version
#
def get_version (paths=None):
    """
    paths:
        a VERSION file containing the long version is checked for in every
        directory listed in paths.   When we find a VERSION file, we also look
        for an SDIST file, and return the found name and location as absolute
        path to the sdist.
    """

    if not paths:
        # by default, get version for myself
        pwd   = os.path.dirname (__file__)
        root  = "%s/.." % pwd
        paths = [root, pwd]
    
    if not isinstance (paths, list):
        paths = [paths]
    
    long_version  = None
    short_version = None
    branch_name   = None
    sdist_name    = None
    sdist_path    = None
    version_path  = None
    
    
    # if in any of the paths a VERSION file exists, we use the long version
    # in there.
    for path in paths:
    
        try:
    
            version_path = "%s/VERSION" % path
    
            with open (version_path) as f:
                line = f.readline()
                line.strip()
                pattern = re.compile ('^\s*(?P<long>(?P<short>[^-@]+?)(-[^@]+?)?(?P<branch>@.+?)?)\s*$')
                match   = pattern.search (line)
    
                if match:
                    long_version  = match.group ('long')
                    short_version = match.group ('short')
                    branch_name   = match.group ('branch')
                  # print 'reading  %s' % version_path
                    break
    
        except Exception as e:
            # ignore missing VERSION file -- this is caught below
            pass


    if long_version:
        # check if there is also an SDIST near the version_path
        sdist_path = version_path.replace ('/VERSION', '/SDIST')
        try:
            sdist_name = open(sdist_path).read().strip()
        except Exception as e:
            # ignore missing SDIST file
            pass
        
        sdist_path = version_path.replace ('/VERSION', '/%s' % sdist_name)

    # check if any one worked ok
    if long_version:
        return short_version, long_version, branch_name, sdist_name, sdist_path
    else:
        raise RuntimeError ("Cannot determine version from %s" % paths)
    

# ------------------------------------------------------------------------------

