
__author__    = "Radical.Utils Development Team (Andre Merzky, Matteo Turilli)"
__copyright__ = "Copyright 2013, RADICAL@Rutgers"
__license__   = "MIT"


import sys
import string
import singleton
# import colorama as c

# ------------------------------------------------------------------------------
#
class Reporter(object):

    # we want reporter style to be consistent in the scope of an application
    __metaclass__ = singleton.Singleton

    # COLORS = {'white'       : c.Style.BRIGHT    + c.Fore.WHITE   ,
    #           'yellow'      : c.Style.BRIGHT    + c.Fore.YELLOW  ,
    #           'green'       : c.Style.BRIGHT    + c.Fore.GREEN   ,
    #           'blue'        : c.Style.BRIGHT    + c.Fore.BLUE    ,
    #           'cyan'        : c.Style.BRIGHT    + c.Fore.CYAN    ,
    #           'red'         : c.Style.BRIGHT    + c.Fore.RED     ,
    #           'magenta'     : c.Style.BRIGHT    + c.Fore.MAGENTA ,
    #           'black'       : c.Style.BRIGHT    + c.Fore.BLACK   ,
    #           'darkwhite'   : c.Style.DIM       + c.Fore.WHITE   ,
    #           'darkyellow'  : c.Style.DIM       + c.Fore.YELLOW  ,
    #           'darkgreen'   : c.Style.DIM       + c.Fore.GREEN   ,
    #           'darkblue'    : c.Style.DIM       + c.Fore.BLUE    ,
    #           'darkcyan'    : c.Style.DIM       + c.Fore.CYAN    ,
    #           'darkred'     : c.Style.DIM       + c.Fore.RED     ,
    #           'darkmagenta' : c.Style.DIM       + c.Fore.MAGENTA ,
    #           'darkblack'   : c.Style.DIM       + c.Fore.BLACK   ,
    #           'reset'       : c.Style.RESET_ALL + c.Fore.RESET
    #       }

    COLORS     = {'reset'        : '\033[39m',
                  'black'        : '\033[30m',
                  'red'          : '\033[31m',
                  'green'        : '\033[32m',
                  'yellow'       : '\033[33m',
                  'blue'         : '\033[34m',
                  'magenta'      : '\033[35m',
                  'cyan'         : '\033[36m',
                  'lightgray'    : '\033[37m',
                  'darkgray'     : '\033[90m',
                  'lightred'     : '\033[91m',
                  'lightgreen'   : '\033[92m',
                  'lightyellow'  : '\033[93m',
                  'lightblue'    : '\033[94m',
                  'lightmagenta' : '\033[95m',
                  'lightcyan'    : '\033[96m',
                  'white'        : '\033[97m'}

    COLOR_MODS = {'reset'        : '\033[0m',
                  'bold'         : '\033[1m',
                  'underline'    : '\033[4m',
                  'blink'        : '\033[5m',
                  'inverse'      : '\033[7m',
                  ''             : ''}


    # Define terminal colors for the reporter
    TITLE    = 'bold lightblue'
    HEADER   = 'bold lightyellow'
    INFO     = 'lightblue'
    IDLE     = 'lightwhite'
    PROGRESS = 'lightwhite'
    OK       = 'lightgreen'
    WARN     = 'lightyellow'
    ERROR    = 'lightred'

    EMPTY   = ''
    DOTTED  = '.'
    SINGLE  = '-'
    DOUBLE  = '='
    HASHED  = '#'

    LINE_LENGTH = 80

    # --------------------------------------------------------------------------
    #
    def __init__(self, title=None, targets=['stdout']):

        '''
        settings.style:
          E : empty line
          T : tabulator
          L : line of line segments
          M : text to report
        '''

        self._title    = title
        self._pos      = 0
        self._targets  = targets
        self._settings = {
                'title' : {
                    'color'   : self.TITLE,
                    'style'   : 'ELMLE',
                    'segment' : self.DOUBLE
                    },
                'header' : {
                    'color'   : self.HEADER,
                    'style'   : 'ELME',
                    'segment' : self.SINGLE
                    },
                'info' : {
                    'color'   : self.INFO,
                    'style'   : 'M',
                    'segment' : self.EMPTY
                    },
                'idle' : {
                    'color'   : self.IDLE,
                    'style'   : 'M',
                    'segment' : self.EMPTY
                    },
                'progress' : {
                    'color'   : self.PROGRESS,
                    'style'   : 'M',
                    'segment' : self.EMPTY
                    },
                'ok' : {
                    'color'   : self.OK,
                    'style'   : 'M',
                    'segment' : self.EMPTY
                    },
                'warn' : {
                    'color'   : self.WARN,
                    'style'   : 'M',
                    'segment' : self.EMPTY
                    },
                'error' : {
                    'color'   : self.ERROR,
                    'style'   : 'M',
                    'segment' : self.EMPTY
                    },
                'plain' : {
                    'color'   : '',
                    'style'   : 'M',
                    'segment' : self.EMPTY,
                    }
                }
        self._idle_sequence = '/-\\|'
        self._idle_pos      = dict()

        # set up the output target streams
        self._color_streams = list()
        self._streams       = list()

        for tgt in self._targets:

            if  tgt.lower() == 'stdout':
                self._color_streams.append(sys.stdout)

            elif tgt.lower() == 'stderr':
                self._color_streams.append(sys.stderr)

            else:
                # >>&     color stream in append    mode
                # >&      color stream in overwrite mode (default)
                # >>  non-color stream in append    mode
                # >   non-color stream in overwrite mode

                if   tgt.startswith('>>&'): self._color_streams.append(open(tgt[3:], 'a'))
                elif tgt.startswith('>&') : self._color_streams.append(open(tgt[2:], 'w'))
                elif tgt.startswith('>>') : self._streams.append      (open(tgt[2:], 'a'))
                elif tgt.startswith('>')  : self._streams.append      (open(tgt[1:], 'w'))
                else                      : self._color_streams.append(open(tgt,     'w'))

        # and send the title to all streams
        if self._title:
            self.title(self._title)


    # --------------------------------------------------------------------------
    #
    def _out(self, color, msg):

        color_mod = ''
        if  ' ' in color:
            color_mod, color = color.split(' ', 2)

        color     = self.COLORS.get(color.lower(), '')
        color_mod = self.COLOR_MODS.get(color_mod.lower(), '')

        color += color_mod


        # make sure we count tab length on line start correctly
        msg = msg.replace('\n\t', '\n        ')

        # make sure we don't extent a long line further
        if self._pos >= (self.LINE_LENGTH) and msg and msg[0] != '\n':
            while msg[0] == '\b':
                msg = msg[1:]
            msg = '\n        %s' % msg

        # special control characters:
        #
        #   * '>>' will, at it's place, insert sufficient spaces to make the
        #     remainder of the string right-aligned.  Only one >> is
        #     interpreted, linebreaks before it are ignored
        #
        #   * '<<' will insert a line break if the position is not already on
        #     the beginning of a line
        #
      # print "[%s:%s]" % (self.__hash__(), self._pos),
        slash_f  = msg.find('>>')
        if slash_f >= 0:
            copy   = msg[slash_f+1:].strip()
            spaces = self.LINE_LENGTH - self._pos - len(copy) + 1 # '>>'
            if spaces < 0:
                spaces = 0
            msg = msg.replace('>>', spaces * ' ')

        slash_cr = msg.find('<<')
        if slash_cr >= 0:
          # print "{%s:%s}" % (self.__hash__(), self._pos),
            if self._pos + slash_cr > 0:
                spaces = self.LINE_LENGTH - self._pos - 1
                msg = msg.replace('<<', '%s\\\n' % (spaces * ' '))
            else:
                msg = msg.replace('<<', '')

        mlen  = len(filter(lambda x: x in string.printable, msg))
        mlen -= msg.count('\b')

      # print "<%s>" % (self._pos),
        # find the last \n and then count how many chars we are writing after it
        slash_n = msg.rfind('\n')
        if slash_n >= 0:
          # print "(%s" % (self._pos),
            self._pos = mlen - slash_n - 1
          # print ": %s)" % (self._pos),
        else:
          # print "'%s'[%s" % (msg, self._pos),
            self._pos += mlen
          # print ": %s]" % (self._pos),


        for stream in self._color_streams:
            stream.write(color)
            stream.write(msg)
            stream.write(self.COLORS['reset'])
            stream.write(self.COLOR_MODS['reset'])
            try:
                stream.flush()
            except Exception as e:
                pass

        for stream in self._streams:
            stream.write(msg)
            try:
                stream.flush()
            except Exception as e:
                pass

    # --------------------------------------------------------------------------
    #
    def _format(self, msg, settings=None):

        if not msg:
            msg = ''

        if not settings:
            settings = {}

        color   = settings.get('color',   '')
        style   = settings.get('style',   'M')
        segment = settings.get('segment', '')

        for c in style:

            if  c == 'M':
                self._out(color, "%s" % msg)

            if  c == 'T':
                self._out(color, "\t")

            elif c == 'E':
                self._out(color, "\n")

            elif c == 'L':
                if segment:
                    self._out(color, "%s\n" % (self.LINE_LENGTH * segment))


    # --------------------------------------------------------------------------
    #
    def set_style(self, which, color=None, style=None, segment=None):

        if which not in self._settings:
            raise LookupError('reporter does not support style "%s"' % which)

        settings = self._settings[which]

        if color  : settings['color']   = color
        if style  : settings['style']   = style
        if segment: settings['segment'] = segment


    # --------------------------------------------------------------------------
    #
    def title(self, title=''):

        if not title:
            title = self._title

        if title:
            fmt   = " %%-%ds\n" % (self.LINE_LENGTH-1)
            title = fmt % title

        self._format(title, self._settings['title'])


    # --------------------------------------------------------------------------
    #
    def header(self, msg=''):

        if msg:
            fmt = "%%-%ds\n" % self.LINE_LENGTH
            msg = fmt % msg

        self._format(msg, self._settings['header'])


    # --------------------------------------------------------------------------
    #
    def info(self, msg=''):

        self._format(msg, self._settings['info'])


    # --------------------------------------------------------------------------
    #
    def idle(self, c=None, mode=None, color=None, idle_id=None):

        if not idle_id:
            idle_id = 'default'

        if color: col = self._settings[color]['color']
        else    : col = self._settings['idle']['color']

        idx = 0
        if   mode == 'start': self._out(col, 'O')
        elif mode == 'stop' : self._out(col, '\b ')
        else:
            if not c:
                idx  = self._idle_pos.get(idle_id, 0)
                c    = self._idle_sequence[idx % len(self._idle_sequence)]
                idx += 1
                self._out(col, '\b%s' % c)
            else:
                idx += 1
                self._out(col, '\b%s|' % c)

        self._idle_pos[idle_id] = idx


    # --------------------------------------------------------------------------
    #
    def progress(self, msg=''):

        if not msg:
            msg = '.'
        self._format(msg, self._settings['progress'])


    # --------------------------------------------------------------------------
    #
    def ok(self, msg=''):

        self._format(msg, self._settings['ok'])


    # --------------------------------------------------------------------------
    #
    def warn(self, msg=''):

        self._format(msg, self._settings['warn'])


    # --------------------------------------------------------------------------
    #
    def error(self, msg=''):

        self._format(msg, self._settings['error'])


    # --------------------------------------------------------------------------
    #
    def exit(self, msg='', exit_code=0):

        self.error(msg)
        sys.exit(exit_code)


    # --------------------------------------------------------------------------
    #
    def plain(self, msg=''):

        self._format(msg)


# ------------------------------------------------------------------------------

if __name__ == "__main__":

    import radical.utils as ru

    r = ru.Reporter(title='test')

    r.header  ('header  \n')
    r.info    ('info    \n')
    r.progress('progress\n')
    r.ok      ('ok      \n')
    r.warn    ('warn    \n')
    r.error   ('error   \n')
    r.plain   ('plain   \n')

    i = 0
    j = 0
    for cname,col in r.COLORS.items():
        if cname == 'reset':
            continue
        i+=1
        for mname,mod in r.COLOR_MODS.items():
            if mname == 'reset':
                continue
            j+=1
            sys.stdout.write("%s%s[%12s-%12s] " % (col, mod, cname, mname))
            sys.stdout.write("%s%s" % (r.COLORS['reset'], r.COLOR_MODS['reset']))
        sys.stdout.write("\n")
        j = 0

    import time
  # r.info('test idler:')
  # r.idle(mode='start')
  # for i in range(3):
  #     r.idle()
  #     time.sleep(0.3)
  # r.idle(color='ok', c='.')
  # r.idle(color='error', c='.')
  # for i in range(3):
  #     r.idle()
  #     time.sleep(0.1)
  #
  # r.idle(mode='stop')
  # r.ok('>>done\n')

    r.info('idle test\n')
    r.info('1234567891         2         3         4         5         6         7         8\n\t')
    r.info('.0.........0.........0.........0.........0.........0.........0.........0')
    r.idle(mode='start')
    for i in range(200):
        r.idle(); time.sleep(0.01)
        r.idle(); time.sleep(0.01)
        r.idle(); time.sleep(0.01)
        r.idle(); time.sleep(0.01)
        r.idle(color='ok', c="+")
    r.idle(mode='stop')

    r.set_style('error', color='yellow', style='ELTTMLE', segment='X')
    r.error('error ')

    r.exit    ('exit    \n', 1)

