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

import sys
from mglutil.regression import testplus
ed = None
withThreads = 1 # default is: multi-threading on

# allow for additional user input
if len(sys.argv):
    for myArg in sys.argv[1:]:
        if myArg[:11] == 'withThreads':
            withThreads = int(string.strip(myArg)[-1])


def startEditor():
    global ed
    from Vision.VPE import VisualProgramingEnvironment
    ed = VisualProgramingEnvironment(name='Vision')
    ed.root.update_idletasks()
    ed.configure(withThreads=withThreads)


def quitEditor():
    ed.master.after(1000, ed.exit )
    ed.mainloop()


def pause(sleepTime=0.4):
    from time import sleep
    ed.master.update()
    sleep(sleepTime)


def test_loadDejaVuLib():
    from DejaVu.VisionInterface.DejaVuNodes import vizlib
    ed.addLibrary(vizlib)


def test_allDejaVuNodes():
    # test the symserv nodes
    libs = ed.libraries
    posx = 150
    posy = 150
    for lib in libs.keys():
        ed.ModulePages.selectpage(lib)
        ed.root.update_idletasks()
        for cat in libs[lib].libraryDescr.keys():
            for node in libs[lib].libraryDescr[cat]['nodes']:
                klass = node.nodeClass
                kw = node.kw
                args = node.args
                netNode = apply( klass, args, kw )
                print 'testing: '+node.name # begin node test
                #add node to canvas
                ed.currentNetwork.addNode(netNode,posx,posy)
                # show widget in node if available:
                widgetsInNode = netNode.getWidgetsForMaster('Node')
                if len(widgetsInNode.items()) is not None:
                    for port,widget in widgetsInNode.items():
                        netNode.createOneWidget(port)
                        ed.root.update_idletasks()
                    # and then hide it
                    for port,widget in widgetsInNode.items():
                        netNode.hideInNodeWidget(port.widget)
                        ed.root.update_idletasks()

                # show widgets in param panel if available:
                widgetsInPanel = netNode.getWidgetsForMaster('ParamPanel')
                if len(widgetsInPanel.items()):
                    netNode.paramPanel.show()
                    ed.root.update_idletasks()

                    #and then hide it
                    netNode.paramPanel.hide()
                    ed.root.update_idletasks()
                        
                # and now delete the node
                ed.currentNetwork.deleteNodes([netNode])
                ed.root.update_idletasks()
                
                print 'passed: '+node.name # end node test



harness = testplus.TestHarness( "dejavu",
                                connect = (startEditor, (), {}),
                                funs = testplus.testcollect( globals()),
                                disconnect = quitEditor
                                )

if __name__ == '__main__':
    testplus.chdir()
    print harness
    sys.exit( len( harness))
