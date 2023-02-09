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
## Copyright (c) MGL TSRI 2016
##
################################################################################

# adfr.py
# Our build/install process builds our C++ code into a Python
# extension module, available via ADFRcc/adfrcc.py.
# However, you should now no longer directly "import ADFRcc.adfrcc".
# Instead, import via this file, e.g. as follows:
#     from ADFRcc.adfr import *
# You'll observe that adfr has classes with the same names as the
# Python-visible C++ classes from adfrcc:
# --In some cases, all that's here is a line of form "MyClass = CPP.MyClass".
#   This is just to make it possible to have the safe rule that all
#   C++ classes must be imported from here, not directly from adfrcc.
# --In other cases, we define MyClass here as a subclass of CPP.MyClass.
#   This is needed because Python does not know about C++ pointers,
#   so would free objects it had created but no longer had a
#   direct reference to, oblivious to the fact that other C++ objects
#   were still using these objects. Here we fix this problem by
#   maintaining Python references that parallel the underlying
#   cross-references between the C++ objects.
# Whenever you create a new Python-visible C++ class, remember to add it
# to this file -- either as a simple "MyClass = CPP.MyClass" line or,
# if needed, as a same-named subclass.

import ADFRcc.adfrcc as CPP

class AtomSet(CPP.AtomSet):
    def __init__(self, aAtomSetStatic, aName=""):
        self._atomSetStatic = aAtomSetStatic
        super(AtomSet,self).__init__(aAtomSetStatic, aName)

AtomSetStatic = CPP.AtomSetStatic

AtomType = CPP.AtomType

Element = CPP.Element

class FlexibleReceptorScorer(CPP.FlexibleReceptorScorer):
    def addGridMap(self, aGridMap):
        if not hasattr(self, "_gridMaps"): self._gridMaps = []
        self._gridMaps.append(aGridMap)
        super(FlexibleReceptorScorer,self).addGridMap(aGridMap)
    def setLigandAtomSet(self, aLigandAtomSet):
        self._ligandAtomSet = aLigandAtomSet
        super(FlexibleReceptorScorer,self).setLigandAtomSet(aLigandAtomSet)
    def setRrAtomSet(self, aRrAtomSet):
        self._rrAtomSet = aRrAtomSet
        super(FlexibleReceptorScorer,self).setRrAtomSet(aRrAtomSet)
    def setFrAtomSet(self, aFrAtomSet):
        self._frAtomSet = aFrAtomSet
        super(FlexibleReceptorScorer,self).setFrAtomSet(aFrAtomSet)

# For the FT* classes (Flexibility Tree nodes), getting a little fancier
# to avoid needless code replication. Each CPP.FT* class inherits (unchanged)
# the CPP.FTBase methods addChild() and removeChild(). We need to add
# some safe retaining to these methods. Multiple inheritance of each of
# our new subclasses from _FTAddChild should accomplish that:
class _FTAddChild:
    def addChild(self, child):
        if not hasattr(self, "_children"): self._children = set()
        self._children.add(child)
        child._parent = self
        CPP.FTBase.addChild(self, child)
    def removeChild(self, child):
        self._children.remove(child)
        child._parent = None
        CPP.FTBase.removeChild(self, child)

class FTAtomsNode(_FTAddChild, CPP.FTAtomsNode): pass

class FTBase(_FTAddChild, CPP.FTBase): pass

class FTDeclareAtomSet(_FTAddChild, CPP.FTDeclareAtomSet):
    def setAtomSet(self, aAtomSet):
        self._declaredAtomSet = aAtomSet
        super(FTDeclareAtomSet,self).setAtomSet(aAtomSet)

class FTDiscreteConformation(_FTAddChild, CPP.FTDiscreteConformation): pass

class FTDiscreteTransformation(_FTAddChild, CPP.FTDiscreteTransformation): pass

class FTDiscreteTranslation(_FTAddChild, CPP.FTDiscreteTranslation): pass

class FTRotationAboutBond(_FTAddChild, CPP.FTRotationAboutBond): pass

class FTRotationAboutPointQuat(_FTAddChild, CPP.FTRotationAboutPointQuat): pass

class FTTorsion(_FTAddChild, CPP.FTTorsion):
    def setDihedralAtomIndexes(self, aAtomIndex1, aAtomIndex2, aAtomIndex3, aAtomIndex4, aAtomSet=None):
        self._declaredAtomSet = aAtomSet
        super(FTTorsion,self).setDihedralAtomIndexes(aAtomIndex1, aAtomIndex2, aAtomIndex3, aAtomIndex4, aAtomSet)

class Genome(CPP.Genome):
    def __init__(self, aMotionNodes, aFtRoot, aScorer):
        self._motionNodes = aMotionNodes
        self._ftRoot = aFtRoot
        self._scorer = aScorer
        super(Genome,self).__init__(aMotionNodes, aFtRoot, aScorer)

GridMap = CPP.GridMap

class GridScorer(CPP.GridScorer):
    def initialize(self, aAtomSet, aName=""):
        self._atomSet = aAtomSet;
        super(GridScorer,self).initialize(aAtomSet, aName)
    def addGridMap(self, aGridMap):
        if not hasattr(self, "_gridMaps"): self._gridMaps = []
        self._gridMaps.append(aGridMap)
        super(GridScorer,self).addGridMap(aGridMap)

class Individual(CPP.Individual):
    # C++ Individual has two constructors: standard constructor (aArg of class Genome,
    # aInitialGenes possibly non-null) and copy constructor (aArg of class Individual,
    # aInitialGenes null):
    def __init__(self, aArg, aInitialGenes=None):
        self._genome = aArg.getGenome() if isinstance(aArg, CPP.Individual) else aArg
        super(Individual,self).__init__(aArg, aInitialGenes)

class PairwiseScorer(CPP.PairwiseScorer):
    def initialize(self, aAtomSet1, aAtomSet2, aName=""):
        self._atomSet1 = aAtomSet1;
        self._atomSet2 = aAtomSet2;
        super(PairwiseScorer,self).initialize(aAtomSet1, aAtomSet2, aName)

Logger = CPP.Logger

MathTablesManager = CPP.MathTablesManager

Parameters = CPP.Parameters

Randomize = CPP.Randomize

class RigidReceptorScorer(CPP.RigidReceptorScorer):
    def addGridMap(self, aGridMap):
        if not hasattr(self, "_gridMaps"): self._gridMaps = []
        self._gridMaps.append(aGridMap)
        super(RigidReceptorScorer,self).addGridMap(aGridMap)
    def setLigandAtomSet(self, aLigandAtomSet):
        self._ligandAtomSet = aLigandAtomSet
        super(RigidReceptorScorer,self).setLigandAtomSet(aLigandAtomSet)
    def setRrAtomSet(self, aRrAtomSet):
        self._rrAtomSet = aRrAtomSet
        super(RigidReceptorScorer,self).setRrAtomSet(aRrAtomSet)

Scorer = CPP.Scorer

SolisWets = CPP.SolisWets
