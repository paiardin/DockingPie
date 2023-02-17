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

########################################################################
#
#    Vision Macro - Python source code - file generated by vision
#    Thursday 22 December 2005 12:10:49 
#    
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Daniel Stoffler, Michel Sanner and TSRI
#   
# revision: Guillaume Vareille
#  
#########################################################################
#
# $Header: /mnt/raid/services/cvs/DejaVu2/VisionInterface/StereoSep.py,v 1.1.1.1.4.1 2017/07/13 22:20:08 annao Exp $
#
# $Id: StereoSep.py,v 1.1.1.1.4.1 2017/07/13 22:20:08 annao Exp $
#

from NetworkEditor.macros import MacroNode
class StereoSep(MacroNode):

    def __init__(self, constrkw={}, name='StereoSep', **kw):
        kw['name'] = name
        apply( MacroNode.__init__, (self,), kw)

    def beforeAddingToNetwork(self, net):
        MacroNode.beforeAddingToNetwork(self, net)
        ## loading libraries ##
        from DejaVu2.VisionInterface.DejaVu2Nodes import vizlib
        net.editor.addLibraryInstance(vizlib,"DejaVu2.VisionInterface.DejaVu2Nodes", "vizlib")

        from Vision.StandardNodes import stdlib
        net.editor.addLibraryInstance(stdlib,"Vision.StandardNodes", "stdlib")


    def afterAddingToNetwork(self):
        from NetworkEditor.macros import MacroNode
        MacroNode.afterAddingToNetwork(self)
        ## loading libraries ##
        from DejaVu2.VisionInterface.DejaVu2Nodes import vizlib
        from Vision.StandardNodes import stdlib
        ## building macro network ##
        StereoSep_0 = self
        from traceback import print_exc

        ## loading libraries ##
        from DejaVu2.VisionInterface.DejaVu2Nodes import vizlib
        self.macroNetwork.getEditor().addLibraryInstance(vizlib,"DejaVu2.VisionInterface.DejaVu2Nodes", "vizlib")

        from Vision.StandardNodes import stdlib
        self.macroNetwork.getEditor().addLibraryInstance(stdlib,"Vision.StandardNodes", "stdlib")

        try:

            ## saving node input Ports ##
            input_Ports_1 = self.macroNetwork.ipNode
            input_Ports_1.move(176, 6)
        except:
            print "WARNING: failed to restore MacroInputNode named input Ports in network self.macroNetwork"
            print_exc()
            input_Ports_1=None
        try:

            ## saving node output Ports ##
            output_Ports_2 = self.macroNetwork.opNode
            output_Ports_2.move(194, 314)
        except:
            print "WARNING: failed to restore MacroOutputNode named output Ports in network self.macroNetwork"
            print_exc()
            output_Ports_2=None
        try:

            ## saving node Get cameras ##
            from Vision.StandardNodes import GetAttr
            Get_cameras_3 = GetAttr(constrkw = {}, name='Get cameras', library=stdlib)
            self.macroNetwork.addNode(Get_cameras_3,275,67)
            apply(Get_cameras_3.getInputPortByName('objects').configure, (), {'datatype': 'viewer'})
            apply(Get_cameras_3.getInputPortByName('attr').widget.configure, (), {'choices': ('cameras',)})
            Get_cameras_3.getInputPortByName("attr").widget.set("cameras")
        except:
            print "WARNING: failed to restore GetAttr named Get cameras in network self.macroNetwork"
            print_exc()
            Get_cameras_3=None
        try:

            ## saving node Slice Data ##
            from Vision.StandardNodes import SliceData
            Slice_Data_4 = SliceData(constrkw = {}, name='Slice Data', library=stdlib)
            self.macroNetwork.addNode(Slice_Data_4,301,167)
            apply(Slice_Data_4.getInputPortByName('data').configure, (), {'datatype': 'list'})
            Slice_Data_4.getInputPortByName("_slice").widget.set("[0]")
        except:
            print "WARNING: failed to restore SliceData named Slice Data in network self.macroNetwork"
            print_exc()
            Slice_Data_4=None
        try:

            ## saving node Dial ##
            from Vision.StandardNodes import DialNE
            Dial_5 = DialNE(constrkw = {}, name='Dial', library=stdlib)
            self.macroNetwork.addNode(Dial_5,497,16)
            Dial_5.getInputPortByName("dial").widget.set(0.159390991469)
        except:
            print "WARNING: failed to restore DialNE named Dial in network self.macroNetwork"
            print_exc()
            Dial_5=None
        try:

            ## saving node Redraw ##
            from DejaVu2.VisionInterface.DejaVu2Nodes import Redraw
            Redraw_6 = Redraw(constrkw = {}, name='Redraw', library=vizlib)
            self.macroNetwork.addNode(Redraw_6,199,232)
        except:
            print "WARNING: failed to restore Redraw named Redraw in network self.macroNetwork"
            print_exc()
            Redraw_6=None
        try:

            ## saving node call method ##
            from Vision.StandardNodes import CallMethod
            call_method_7 = CallMethod(constrkw = {}, name='call method', library=stdlib)
            self.macroNetwork.addNode(call_method_7,353,228)
            apply(call_method_7.addInputPort, (), {'datatype': 'float', 'width': 12, 'required': False, 'name': 'sideBySideTranslation', 'height': 12})
            call_method_7.getInputPortByName("signature").widget.set("Set sideBySideTranslation")
        except:
            print "WARNING: failed to restore CallMethod named call method in network self.macroNetwork"
            print_exc()
            call_method_7=None
        self.macroNetwork.freeze()

        ## saving connections for network StereoSep ##
        if Get_cameras_3 is not None and Slice_Data_4 is not None:
            self.macroNetwork.connectNodes(
                Get_cameras_3, Slice_Data_4, "attrs", "data", blocking=True)
        if Dial_5 is not None and call_method_7 is not None:
            self.macroNetwork.connectNodes(
                Dial_5, call_method_7, "value", "sideBySideTranslation", blocking=True)
        if Slice_Data_4 is not None and call_method_7 is not None:
            self.macroNetwork.connectNodes(
                Slice_Data_4, call_method_7, "data", "objects", blocking=True)
        if call_method_7 is not None and Redraw_6 is not None:
            self.macroNetwork.connectNodes(
                call_method_7, Redraw_6, "objects", "trigger", blocking=True)
        input_Ports_1 = self.macroNetwork.ipNode
        if input_Ports_1 is not None and Redraw_6 is not None:
            self.macroNetwork.connectNodes(
                input_Ports_1, Redraw_6, "new", "viewer", blocking=True)
        if input_Ports_1 is not None and Get_cameras_3 is not None:
            self.macroNetwork.connectNodes(
                input_Ports_1, Get_cameras_3, "Redraw_viewer", "objects", blocking=True)
        self.macroNetwork.unfreeze()

        StereoSep_0.shrink()

        ## reset modifications ##
        StereoSep_0.resetTags()
        StereoSep_0.buildOriginalList()
