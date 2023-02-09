################################################################################
##
## This library is free software; you can redistribute it and/or
## modify it under the terms of the GNU Lesser General Public
## License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## 
## This library is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## Lesser General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public
## License along with this library; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
##
## (C) Copyrights Dr. Michel F. Sanner and TSRI 2016
##
################################################################################

import os

class RepeatPrinter:
    """
    The class allows printing repeating messages on the last line of a terminal (Unix only). This can be used to display a counter that does not scroll the terminal

    example:
        printer = RepeatPrinter()
        printer.prompt('Hello World')
        for i in range(20):
            printer.update(str(i))

    thsi code will display 'Hello World' followed by a counter
    """

    def __init__(self):
        # use tcup lines to find out how many lines in the terminal
        self.nbLines = int(os.popen("tput lines").readlines()[0])
        self.repeatingPosition = None


    def prompt(self, prompt):
        # print a message leading the repeating messages

        #os.system("tput sc")
        os.system("tput cup %d 0"%(self.nbLines-1))
        self.repeatingPosition = len(prompt) + 1  
        print prompt
        
    def update(self, msg):
        os.system("tput cup %d %d"%(self.nbLines-2, self.repeatingPosition))
        print msg



if __name__ == '__main__':
    from time import sleep

    printer = RepeatPrinter()
    printer.prompt('Hello World')

    sleep(0.5)
    for i in range(20):
        printer.update(str(i))
        sleep(0.2)

