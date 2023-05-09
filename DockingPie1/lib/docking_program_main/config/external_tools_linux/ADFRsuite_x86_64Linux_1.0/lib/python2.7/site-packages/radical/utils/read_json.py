
__author__    = "Radical.Utils Development Team (Andre Merzky)"
__copyright__ = "Copyright 2013, RADICAL@Rutgers"
__license__   = "MIT"


import re
import json


# ------------------------------------------------------------------------------
#
def read_json (filename) :
    """
    Comments in the form of
        # rest of line
    are stripped from json before parsing

    use like this::

        import pprint
        pprint.pprint (read_json (sys.argv[1]))

    """

    with open (filename) as f:

        content = ''

        # weed out comments
        for line in f.readlines () :
            content += re.sub (r'^\s*#.*', '', line)

        return json.loads (content)


# ------------------------------------------------------------------------------
#
def read_json_str (filename) :
    """
    same as read_json, but converts unicode strings to simple strings.
    This could trivially done with YAML, but that introduces a compile
    dependency -- so we convert manually...
    """

    return _unicode_to_strings (read_json (filename))


# ------------------------------------------------------------------------------
#
def write_json (data,  filename) :
    """
    thin wrapper around python's json write, for consistency of interface

    """

    # thanks to
    # http://stackoverflow.com/questions/956867/#13105359
    def _byteify(input):
        if isinstance(input, dict):
            return {_byteify(key):_byteify(value) for key,value in input.iteritems()}
        elif isinstance(input, list):
            return [_byteify(element) for element in input]
        elif isinstance(input, unicode):
            return input.encode('utf-8')
        else:
            return input


    with open (filename, 'w') as f :
        json.dump (_byteify(data), f, sort_keys=True, indent=4, ensure_ascii=False)
        f.write ("\n")


# ------------------------------------------------------------------------------
#
def parse_json (json_str, filter_comments=True) :
    """
    Comments in the form of
        # rest of line
    are stripped from json before parsing
    """

    if not filter_comments :
        return json.loads (json_str)

    else :
        content = ''
        for line in json_str.split ('\n') :
            content += re.sub (r'^\s*#', '', line)
            content += '\n'

        return json.loads (content)


# ------------------------------------------------------------------------------
#
def parse_json_str (json_str) :
    """
    same as parse_json, but converts unicode strings to simple strings
    """

    return _unicode_to_strings (parse_json (json_str))


# ------------------------------------------------------------------------------
#
def _unicode_to_strings (unicode_data) :

    if isinstance (unicode_data, dict) :
        ret = dict ()
        for key, value in unicode_data.iteritems() :
            ret[_unicode_to_strings(key)] = _unicode_to_strings(value) 
        return ret

    elif isinstance (unicode_data, list) :
        ret = list ()
        for element in unicode_data :
            ret.append (_unicode_to_strings (element))
        return ret

    elif isinstance (unicode_data, unicode) :
        return unicode_data.encode ('utf-8')

    else:
        return unicode_data


# ------------------------------------------------------------------------------
