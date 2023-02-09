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

#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner and TSRI 2015
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/ADFR/torTreeFromPDBQT.py,v 1.3 2017/07/08 02:10:18 sanner Exp $
#
# $Id: torTreeFromPDBQT.py,v 1.3 2017/07/08 02:10:18 sanner Exp $
#
import weakref

class TreeNode:
    def __init__(self, nodeNum, parent=None):
        self.num = nodeNum
        if parent is not None:
            parent = weakref.ref(parent)
        self.parent = parent # weakref to parent or None
        self.atoms = [] # atom indices
        self.bond = [] # atom indices of rotatable bond leading to this node
        self.distanceFromLeaf = 0
        self.children = [] # list of TreeNodes

class TorTreeFromPDBQT:
    """
    build torsion tree from the pdbqt source code
    """

    def _addTorsionInfo(self, root, child):
        ## recursive function to create node.torsionIndices and
        ## node.torsionAtomNames. The indices are corresponding to
        ## the molecule's atom indices
        i1, i2 = [root.serialToIndex[x] for x in child.bond]
        a1, a2 = root.mol._ag[[i1,i2]]
        for b1 in a1.iterBonded():
            if b1!=a2:
                break
        for b2 in a2.iterBonded():
            if b2!=a1:
                break
        #print b1.getName(),a1.getName(),a2.getName(),b2.getName()
        child.torsionIndices = (b1.getIndex(),a1.getIndex(),a2.getIndex(),b2.getIndex())
        child.torsionAtomNames = (b1.getName(),a1.getName(),a2.getName(),b2.getName())
        if len(child.children): # not a leaf
            for c in child.children:
                self._addTorsionInfo(root, c)
                if child.distanceFromLeaf < c.distanceFromLeaf + 1:
                    child.distanceFromLeaf = c.distanceFromLeaf + 1
        else:
            child.distanceFromLeaf=0

    def addTorsionInfo(self, root):
        for child in root.children:
            self._addTorsionInfo(root, child)

    def __call__(self, filename, mol=None):

        if mol is None:
            from MolKit2 import Read
            mol = Read(filename)

        f = open(filename)
        lines = f.readlines()
        f.close()
        serialToIndex = {}
        atCount = 0
        hasModels = False
        for line in lines:
            if line.startswith('MODEL'):
                hasModels = True
                continue
            if hasModels and line.startswith('ENDMDL'):
                break
            if line.startswith('ROOT'):
                root = TreeNode(0)
                root.mol = mol
                nodeNum = 0
                currentNode = root

            elif line.startswith('BRANCH'):
                nodeNum += 1
                newNode = TreeNode(nodeNum, parent=currentNode)
                w = line.split()
                bat1 = int(w[1])
                bat2 = int(w[2])
                newNode.bond = (bat1, bat2)
                currentNode.children.append(newNode)
                currentNode = newNode

            elif line.startswith('ENDBRANCH'):
                currentNode = currentNode.parent()

            elif line.startswith('ATOM') or line.startswith('HETATM'):
                serial = int(line.split()[1])
                serialToIndex[serial] = atCount
                currentNode.atoms.append(serial)
                atCount += 1

            elif line.startswith('TORSDOF'):
                root.torsdof = int(line.split()[1])
        root.serialToIndex = serialToIndex
        self.addTorsionInfo(root)
        
        return root

def weedBonds(mol, torTree):
    """
    removed <- weedBonds(mol, torTree)

    tree - ligandTorTree

    return a list (i,j) atom indices for bonds between atoms in the same
    rigid body
    """
    # define a recursive function that will walk the tree
    # and remove add to removedPairs all pairs of atoms
    # located in the same rigid body (i.e. node in the troTree)
    # as well as any inteaction between atoms in the node and
    # the 2 atoms of the bond traversed to reach this node

    def _handleRigidBody(root, node):#, removedPairs):
        # first extend this node with the second atom
        # from each bond connecting to a child node
        atoms = [root.serialToIndex[x] for x in node.atoms]
        for c in node.children:
            atoms.append(root.serialToIndex[c.bond[1]])
        if node.bond: # handle bond through which we came here from parent
            atoms.append(root.serialToIndex[node.bond[0]])
        for i, ai in enumerate(atoms):
            for bi in atoms[i+1:]:
                removedPairs.append((ai,bi))

        for c in node.children:
            _handleRigidBody(root, c)

    # call recursive function to populate self.removedPairs
    removedPairs = []
    _handleRigidBody(torTree, torTree)

    return removedPairs
