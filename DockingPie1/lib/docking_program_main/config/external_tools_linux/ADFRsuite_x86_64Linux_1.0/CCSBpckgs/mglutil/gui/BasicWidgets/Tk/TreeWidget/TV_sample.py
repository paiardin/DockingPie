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

from TreeWidget.tree import TreeView

if __name__ == '__main__':
    tv = TreeView()
    # addNode(nodename, parentname = None)
    tv.addNode('protein_1')
    tv.addNode('residue_11',    parent='protein_1')
    tv.addNode('AminoAcid',     parent='protein_1|residue_11')
    tv.addNode('A',             parent='protein_1|residue_11|AminoAcid')
    tv.addNode('H',             parent='protein_1|residue_11|AminoAcid')
   

    tv.addNode('protein_2')
    tv.addNode('protein_3')


    tv.addNode('residue_21',    parent='protein_2')
    tv.addNode('residue_25',    parent='protein_2')
    tv.addNode('basdfe',        parent='protein_2|residue_21')
    tv.addNode('AminoAcid',     parent='protein_2|residue_21')
   
 
    tv.addNode('etc',       parent='protein_1|residue_11')
    tv.addNode('color',     parent='protein_1|residue_11|etc')
    tv.addNode('density',   parent='protein_1|residue_11|etc')
    tv.addNode('residue_12',parent='protein_1')
    
  
    for a in range(1):
        name = 'AA' + str(a)
        tv.addNode(name, parent='protein_1|residue_11|AminoAcid')


    tv.addNode('2', parent='protein_2|residue_21')
    tv.addNode('3', parent='protein_2|residue_21')
    tv.addNode('4', parent='protein_2|residue_21')

    tv.addNode('L', parent='protein_2|residue_21|AminoAcid')
    tv.addNode('S', parent='protein_2|residue_21|AminoAcid')
 
    
    for a in range(10):
        name = 'A' + str(a)
        tv.addNode(name,  parent='protein_2|residue_21|AminoAcid')


    tv.addNode('protein_4')
    tv.addNode('residue_22', parent='protein_2')

##  to delete a node:
#   tv.deleteNode(nodename, parentname)
# e.g. >>> tv.deleteNode('residue_21', 'protein_2')
# e.g. >>> tv.deleteNode('protein_2')
