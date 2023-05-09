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

# test.py
# This file is a Python equivalent of ADFRcc/src/ADFRTest.h,cpp.
# It can be used to help test that the ADFR C++ code, and the
# swig interface to that code, have been built into a correctly
# functioning Python extension module, ADFRcc.adfrcc.

# NEW: Rather than importing directly from ADFRcc.adfrcc, we now import
# via adfr. Also note that, since test.py is in the same directory as adfr.py,
# I use "from adfr import *" to import from this local/original version of adfr.py;
# this is convenient when adfr.py is being frequently revised. More general code
# would use "from ADFRcc.adfr import *", which references the copy of
# adfr.py that is in the ADFRcc package's installation directory:
from adfr import *
import numpy
import os

# Set FTYPE to 'd' if using double-precision numpy floats here,
# or 'f' for single-precision floats.
# (It should preferably be possible to use either 'd' or 'f' here,
# independent of whether USE_64_BIT_FLOATS is true or false in the
# C++ code.  But in fact, the combination of FTYPE='d' with
# USE_64_BIT_FLOATS=false currently does not work.)
FTYPE = 'd'

# Loads the default atomic parameters, then logs a (lengthy)
# description of them.
def testParameters():
    print "\nBeginning test of the atomic parameters:"
    parameters = Parameters.getParameters()
    parameters.printDebugDescription()

# Allocates, initializes, and returns a very simple AtomSetStatic
# suitable for testing purposes.  To keep things really simple,
# there are no hbonding atoms in this AtomSetStatic.
def createAtomSetStaticTypeA():
    numAtoms = 4
    atomTypeNames = ["H", "C", "N", "Zn"]
    charges = numpy.array([0.5, 0.3, -1., 0.], FTYPE)
    origCoords =  numpy.array([
        [1., 1., 1.],
        [3., 1., 1.],
        [10., 10., 10.],
        [11., 11., 11.]], FTYPE)
    # The simplest case: none of the atoms are covalently bonded to each other:
    covalentBonds = numpy.array([-1], 'i')
    
    atomSetStatic = AtomSetStatic(numAtoms,  "simple_a")
    atomSetStatic.setAtomTypes(atomTypeNames)
    atomSetStatic.setCharges(charges)
    atomSetStatic.setOrigCoords(origCoords)
    atomSetStatic.setCovalentBonds(covalentBonds)
    return atomSetStatic

# Similar to createAtomSetStaticTypeA(), except the new AtomSetStatic
# is a little less trivial.  In particular, it includes several hbond
# donors and acceptors.
def createAtomSetStaticTypeB():
    numAtoms = 14
    # * Atom 0 (NX) is a directional donor; its hydrogen neighbor Atom 1 (HD)
    #   determines its directionality.  And similarly for Atom 4 and Atom 5,
    #   and for Atom 9 and Atom 10.
    # * Atom 2 (NA) is a unidirectional acceptor; its covalent neighbor Atom 3 (C)
    #   determines its directionality.
    # * Atom 6 (OA) is a bidirectional acceptor.  Its directionality is
    #   determined by its (one) covalent neighbor Atom 7 (C),
    #   and by Atom 7's neighbor Atom 8 (C).
    # * Atom 11 (OA) is a bidirectional acceptor.  Its directionality is
    #   determined by its (two) neighbors, Atom 12 (C) and Atom 13 (C):
    atomTypeNames = ["NX", "HD", "NA", "C", "NX", "HD", "OA", "C", "C", "NX", "HD", "OA", "C", "C"]
    charges = numpy.array([-0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.4, 0.4, 0, -0.5, 0.5, -0.6, 0.3, 0.3], 'd')
    origCoords = numpy.array([
        [0., 0., 0.],
        [1., 0., 0.],
        [3.5, 0.5, 0.],
        [5., 0.5, 0.],
        [0., 5., 0.],
        [1., 5., 0.],
        [3.5, 5.5, 0.],
        [5, 5.5, 0.],
        [6.5, 6., 0.],
        [0., 10., 0.],
        [1., 10., 0.],
        [3.5, 10.5, 0.],
        [5., 9.5, 0.],
        [5., 12., 0.]], FTYPE)
    # Declares covalent bonds between Atom 0 and Atom 1, between Atom 2 and Atom 3, etc.:
    covalentBonds = numpy.array([0, 1, 2, 3, 4, 5, 6, 7, 7, 8, 9, 10, 11, 12, 11, 13, -1], 'i')
    
    atomSetStatic = AtomSetStatic(numAtoms, "simple_b")
    atomSetStatic.setAtomTypes(atomTypeNames)
    atomSetStatic.setCharges(charges)
    atomSetStatic.setOrigCoords(origCoords)
    atomSetStatic.setCovalentBonds(covalentBonds)
    # Not really needed, as 3 is the default value; but for reference:
    atomSetStatic.setNoScoreHops(3)
    return atomSetStatic

# Similar to createAtomSetStaticTypeA().  Creates and returns a very simple
# AtomSetStatic that I use to test the GridScorer.  Coordinates and atom-types
# are designed to match the grid from the 1ppc test-data.
def createAtomSetStaticTypeC():
    numAtoms = 4
    atomTypeNames = ["N", "C", "N", "S"]
    charges = numpy.array([-0.5, 0.3, -1., 0.2], FTYPE)
    origCoords = numpy.array([
        [-13.6, -1.4, 15.9],
        [-5.1, -0.5, 24.9],
        [-3.8, -1.6, 34.9],
        [4.5, -3.8, 25.8]], FTYPE)
    covalentBonds = numpy.array([-1], 'i')
    
    atomSetStatic = AtomSetStatic(numAtoms, "simple_c")
    atomSetStatic.setAtomTypes(atomTypeNames)
    atomSetStatic.setCharges(charges)
    atomSetStatic.setOrigCoords(origCoords)
    atomSetStatic.setCovalentBonds(covalentBonds)
    return atomSetStatic

# Similar to createAtomSetStaticTypeC().  Creates and returns a very simple
# AtomSetStatic suitable for GridScoring against the grid from the
# 1ppc test-data.
def createAtomSetStaticTypeD():
    numAtoms = 4
    atomTypeNames = ["NA", "C", "NA", "S"]
    charges = numpy.array([-0.5, 0.3, -1., 0.2], FTYPE)
    origCoords = numpy.array([
        [-8.6, 3.6, 20.9],
        [-3.1, 1.5, 22.9],
        [-1.8, 0.4, 32.9],
        [2.5, -1.8, 23.8]], FTYPE)
    covalentBonds = numpy.array([0, 1, 2, 3, -1], 'i')
    
    atomSetStatic = AtomSetStatic(numAtoms, "simple_d")
    atomSetStatic.setAtomTypes(atomTypeNames)
    atomSetStatic.setCharges(charges)
    atomSetStatic.setOrigCoords(origCoords)
    atomSetStatic.setCovalentBonds(covalentBonds)
    return atomSetStatic

# Another simple test AtomSetStatic.
def createAtomSetStaticTypeE():
    numAtoms = 2
    atomTypeNames = ["NX", "HD"]
    charges = numpy.array([-0.5, 0.3], FTYPE)
    origCoords = numpy.array([
        [-6.6, 1.6, 18.9],
        [-7.6, 2.6, 19.7]], FTYPE)
    covalentBonds = numpy.array([0, 1, -1], 'i')
    
    atomSetStatic = AtomSetStatic(numAtoms, "simple_e")
    atomSetStatic.setAtomTypes(atomTypeNames)
    atomSetStatic.setCharges(charges)
    atomSetStatic.setOrigCoords(origCoords)
    atomSetStatic.setCovalentBonds(covalentBonds)
    return atomSetStatic

# Another simple test AtomSetStatic, this one used in particular
# to test Flexibility Tree manipulations. The data here is inspired
# by an ASP residue, but I have made no effort to set realistic
# atom-types or charges, and no hydrogens have been added.
def createAtomSetStaticTypeF():
    numAtoms = 8
    # Atom indexes:          0    1    2    3    4    5    6     7
    # PDB-style atom-names:  N    CA   C    O    CB   CG   OD1   OD2
    atomTypeNames =        ["N", "C", "C", "O", "C", "C", "OA", "OA"]
    charges = numpy.array([-0.2, 0.1, 0.1, -0.2, 0.1, 0.1, -0.3, -0.3], FTYPE)
    origCoords = numpy.array([
        [0.903, 3.406, 5.025],
        [1.455, 2.603, 3.941],
        [2.946, 2.929, 3.788],
        [3.365, 4.076, 3.950],
        [0.713, 2.918, 2.640],
        [1.111, 2.001, 1.503],
        [0.774, 0.800, 1.544],
        [1.766, 2.478, 0.563]], FTYPE)
    covalentBonds = numpy.array([0,1, 1,2, 1,4, 2,3, 4,5, 5,6, 5,7, -1], 'i')
    
    atomSetStatic = AtomSetStatic(numAtoms, "simple_f")
    atomSetStatic.setAtomTypes(atomTypeNames)
    atomSetStatic.setCharges(charges)
    atomSetStatic.setOrigCoords(origCoords)
    atomSetStatic.setCovalentBonds(covalentBonds)
    return atomSetStatic

# Creates an AtomSet and logs a description of it.
# You may instead want to use testPairwiseScorer(), which does all
# that this function does, and more.
def testAtomSet():
    print "\nBeginning test of creating an AtomSet:"
    atomSetStatic = createAtomSetStaticTypeB()
    atomSet = AtomSet(atomSetStatic)
    atomSet.printDebugDescription()

# Creates two AtomSets, and logs a description of each one.
# Then scores the first against itself, the first against the second,
# and the second against itself.  Logs a description of each of
# the three scoring results, as well as several additional tests
# of the PairwiseScorer get*() score-info accessors.
def testPairwiseScorer():
    print "\nBeginning test of pairwise scoring:"
    atomSetA = AtomSet(createAtomSetStaticTypeA())
    atomSetA.printDebugDescription()
    
    atomSetB = AtomSet(createAtomSetStaticTypeB())
    atomSetB.printDebugDescription()
    
    # Score A against itself:
    scorerAA = PairwiseScorer()
    scorerAA.initialize(atomSetA, atomSetA)
    scorerAA.calculateScores() 
    scorerAA.printDebugDescription()
    
    # Score A against B:
    scorerAB = PairwiseScorer()
    scorerAB.initialize(atomSetA, atomSetB)
    scorerAB.calculateScores()
    scorerAB.printDebugDescription()
    
    # Score B against itself:
    scorerBB = PairwiseScorer()
    scorerBB.initialize(atomSetB, atomSetB)
    scorerBB.calculateScores()
    scorerBB.printDebugDescription()
    
    # Also a few tests of the PairwiseScorer get*() methods:
    scorer = scorerBB
    scorerName = scorer.getName()
    print "Testing some of the PairwiseScorer get*() methods (on the scoring-results from scorer '%s'):" % (scorerName,)
    print "isScoringDoneYet:",  scorer.getIsScoringDoneYet()
    print "isAbortedDueToCollisions:", scorer.getIsAbortedDueToCollisions()
    print "numCollisions:", scorer.getNumCollisions()
    print "collisionAtomIndex1:", scorer.getCollisionAtomIndex1()
    print "collisionAtomIndex2:", scorer.getCollisionAtomIndex2()
    print "totalScore:", scorer.getTotalScore()
    print "totalEnergy:", scorer.getTotalEnergy()
    print "totalVdwEnergy:", scorer.getTotalVdwEnergy()
    print "totalHbEnergy:", scorer.getTotalHbEnergy()
    print "totalEstatEnergy:", scorer.getTotalEstatEnergy()
    print "totalSolvEnergy:", scorer.getTotalSolvEnergy()
    print "atomEnergy(0):", scorer.getAtomEnergy(0)
    print "atomVdwEnergy(0):", scorer.getAtomVdwEnergy(0)
    print "atomHbEnergy(0):", scorer.getAtomHbEnergy(0)
    print "atomEstatEnergy(0):", scorer.getAtomEstatEnergy(0)
    print "atomSolvEnergy(0):", scorer.getAtomSolvEnergy(0)
    print "isHbondingPair(0,2):", scorer.getIsHbondingPair(0, 2)
    print "distance(0,2):", scorer.getDistance(0, 2)
    print "energy(0,2):", scorer.getEnergy(0, 2)
    print "vdwEnergy(0,2):", scorer.getVdwEnergy(0, 2)
    print "hbEnergy(0,2):", scorer.getHbEnergy(0, 2)
    print "estatEnergy(0,2):", scorer.getEstatEnergy(0, 2)
    print "solvEnergy(0,2):", scorer.getSolvEnergy(0, 2)    
    vdwEnergyMatrix = scorer.getVdwEnergyMatrix()
    print "For the scoring-results from scorer '%s', the complete vdw energy-matrix is:" % (scorerName,)
    print vdwEnergyMatrix

# Data paths/filenames (needed by testGridMaps() and testGridScorer()):
# Absolute path of the (installed) ADFRcc/ directory:
import ADFRcc
ADFRcc_DIR = os.path.abspath(ADFRcc.__path__[0])
# Absolute path of the dir holding our test-data (e.g., .map files):
TEST_DATA_PATH = os.path.join(ADFRcc_DIR, "Tests/Data/1ppc")
# Each member of MAP_FILES is a tuple of two strings:
# * First, a map-file type, which must be "ELECTROSTATIC", "DESOLVATION",
#   or an AD atom-type name.
# * Second, a corresponding map-file name.  This must be the filename
#   of an AutoGrid-style .map file (relative to TEST_DATA_PATH):
MAP_FILES = [
    ("ELECTROSTATIC", "1ppc_rec.e.map"),
    ("DESOLVATION", "1ppc_rec.d.map"),
    ("A", "1ppc_rec.A.map"),
    ("C", "1ppc_rec.C.map"),
    ("HD", "1ppc_rec.HD.map"),
    ("N", "1ppc_rec.N.map"),
    ("NA", "1ppc_rec.NA.map"),
    ("OA", "1ppc_rec.OA.map"),
    ("S", "1ppc_rec.S.map"),
]

# Loads a specified list of AutoGrid-style .map files, and logs
# a description of one of them.  Also, the .map files are supposed
# to be for the same rigid receptor, so should all have the same
# bounding-box etc.; the function tests this, and asserts if it is not true.
# WARNING: The logged description of the .map file includes the map's
# full data, so is very long.  You may instead want to use testGridScorer().
def testGridMaps():
    print "\nBeginning test of loading GridMaps:"
    maps = []
    for (mapType,mapFileName) in MAP_FILES:
        map = GridMap()
        map.loadFromMapFile(mapType, TEST_DATA_PATH, mapFileName)
        if len(maps) == 0:
            map.printDebugDescription()
        else:
            map.checkConsistency(maps[0])
        maps.append(map)

# Creates an AtomSet, and logs a description of it.
# Also loads a specified list of AutoGrid-style .map files.
# Then scores the AtomSet against the grid, and logs the scoring results.
# The current test has been designed to work correctly against the .map files
# in the 1ppc test-data (which is now included in the distribution).
def testGridScorer():
    print "\nBeginning test of grid-scoring:"
    atomSet = AtomSet(createAtomSetStaticTypeC())
    atomSet.printDebugDescription()
    
    gridScorer = GridScorer()
    gridScorer.initialize(atomSet)
    
    for (mapType,mapFileName) in MAP_FILES:
        map = GridMap()
        map.loadFromMapFile(mapType, TEST_DATA_PATH, mapFileName)
        gridScorer.addGridMap(map)
    
    gridScorer.calculateScores()
    gridScorer.printDebugDescription()

# Uses a FlexibleReceptorScorer to test-score a system consisting
# of a ligand AtomSet, a flexible receptor AtomSet, and a
# rigid receptor represented by grid-maps.
# The current test has been designed to work correctly against the .map files
# in the 1ppc test-data (which is now included in the distribution).
def testFlexibleReceptorScorer():
    print "\nBeginning test of a FlexibleReceptorScorer:"
    
    ligandAtomSet = AtomSet(createAtomSetStaticTypeC())
    print "Using this as the ligand AtomSet:"
    ligandAtomSet.printDebugDescription()
    
    frAtomSet = AtomSet(createAtomSetStaticTypeD())
    print "Using this as the flexible receptor AtomSet:"
    frAtomSet.printDebugDescription()
    
    scorer = FlexibleReceptorScorer()
    scorer.setLigandAtomSet(ligandAtomSet)
    scorer.setFrAtomSet(frAtomSet)
    for (mapType,mapFileName) in MAP_FILES:
        map = GridMap()
        map.loadFromMapFile(mapType, TEST_DATA_PATH, mapFileName)
        scorer.addGridMap(map)
    
    scorer.calculateScores()
    scorer.printDebugDescription()

# Uses a FlexibleReceptorScorer to test-score a system consisting
# of a ligand AtomSet, a flexible receptor AtomSet, and a
# rigid receptor represented by both grid-maps and an AtomSet.
# It is assumed that the grid-maps are "vdw only"; hence each
# scoring against the rigid receptor is done using both a GridScorer
# (for vdw/electrostatic/desolvation) and a PairwiseScorer (for hb).
# WARNING: Results will be unrealistic, because (a) our grid-maps
# used here for testing are not actually "vdw only", and
# (b) the rigid receptor AtomSet used here is an arbitrary one
# that bears no relation to the atom-set that generated the grid-maps.
def testMixedScoring():
    print "\nBeginning test of \"mixed\" scoring:"
    
    ligandAtomSet = AtomSet(createAtomSetStaticTypeC())
    print "Using this as the ligand AtomSet:"
    ligandAtomSet.printDebugDescription()
    
    frAtomSet = AtomSet(createAtomSetStaticTypeD())
    print "Using this as the flexible receptor AtomSet:"
    frAtomSet.printDebugDescription()
    
    rrAtomSet = AtomSet(createAtomSetStaticTypeE())
    print "Using this (arbitrarily) as the rigid receptor AtomSet:"
    rrAtomSet.printDebugDescription()
    
    scorer = FlexibleReceptorScorer()
    scorer.setLigandAtomSet(ligandAtomSet)
    scorer.setFrAtomSet(frAtomSet)
    scorer.setRrAtomSet(rrAtomSet)
    for (mapType,mapFileName) in MAP_FILES:
        map = GridMap()
        map.loadFromMapFile(mapType, TEST_DATA_PATH, mapFileName)
        scorer.addGridMap(map)
    
    # If the grid-maps really are "vdw only" (as indicated by the .map files
    # having a header line reading "VDW_ONLY 1"), the scorer will detect this
    # and automatically use mixed scoring.  Since we don't really have any
    # "vdw only" .map files at the moment, we need the following line to
    # force testing of mixed scoring:
    scorer.forceScoreMixed()
    
    scorer.calculateScores()
    scorer.printDebugDescription()

# Creates a Flexibility Tree, and tests that it can transform
# an example AtomSet.  As our example AtomSet, consider
# an AtomSet created per createAtomSetStaticTypeF(), which is
# similar to an ASP residue. Let's also allow for a torsion of
# the dihedral angle defined by atoms 0,1,4,5, and another torsion
# of the dihedral angle defined by atoms 1,4,5,6. Then (per my setup;
# it could be done in several ways) atoms 0,1,2,3,4
# are the "body" that is subject only to the overall translation and rotation;
# atom 5 is a "branch" that is subject to the first torsion as well;
# and atoms 6,7 are a "sub-branch" that is subject to both torsions.
# Here we create a Flexibility Tree capable of representing the above
# rules, then verify that it can carry out example transformations
# of the AtomSet.
# NEW: Test extended to make use of Genome/Individual/SolisWets.
# NEW: Test extended to test FTDiscreteConformation/FTDiscreteTransformation.
def testFlexibilityTree():
    print "\nBeginning test of the Flexibility Tree code:"
    
    # Create an AtomSet to be manipulated:
    atomSet = AtomSet(createAtomSetStaticTypeF())
    print "We will use the following AtomSet to test the Flexibility Tree code:"
    atomSet.printDebugDescription()
    
    # Create the Flexibility Tree:
    ftRoot = FTDeclareAtomSet("Root")
    ftRoot.setAtomSet(atomSet)
    
    ftTrans = FTDiscreteTranslation("Translation")
    ftRoot.addChild(ftTrans)
    ftTrans.setBounds(numpy.array([-5.,-5.,-5.], FTYPE), numpy.array([10.,10.,10.], FTYPE))
    ftTrans.setPreferredPoints(numpy.array([[0.,0.,0.], [1.,1.,1.]], FTYPE))
    
    ftRotQuat = FTRotationAboutPointQuat("PointQuat")
    ftTrans.addChild(ftRotQuat)
    ftRotQuat.setRotPoint(1., 2., 1.)
    
    ftAtoms_Body = FTAtomsNode("Body")
    ftRotQuat.addChild(ftAtoms_Body)
    ftAtoms_Body.setAtomIndexes(numpy.array([0,1,2,3,4], 'i'))
    
    ftTorsion1 = FTTorsion("Torsion1")
    ftAtoms_Body.addChild(ftTorsion1)
    ftTorsion1.setDihedralAtomIndexes(0, 1, 4, 5, atomSet)
    
    ftAtoms_Branch = FTAtomsNode("Branch")
    ftTorsion1.addChild(ftAtoms_Branch)
    ftAtoms_Branch.setAtomIndexes(numpy.array([5], 'i'))
    
    ftTorsion2 = FTTorsion("Torsion2")
    ftAtoms_Branch.addChild(ftTorsion2)
    ftTorsion2.setDihedralAtomIndexes(1, 4, 5, 6, atomSet)
    
    ftAtoms_SubBranch = FTAtomsNode("SubBranch")
    ftTorsion2.addChild(ftAtoms_SubBranch)
    ftAtoms_SubBranch.setAtomIndexes(numpy.array([6,7], 'i'))
    
    # Next create a Scorer. For my rather trivial testing here, I just want
    # to score atomSet against itself. However, the Genome code is generally
    # used with a RigidReceptorScorer or FlexibleReceptorScorer; so as a
    # compromise, I am using a RigidReceptorScorer with only ligand-ligand
    # scoring enabled:
    scorer = RigidReceptorScorer()
    scorer.setLigandAtomSet(atomSet)
    scorer.setLrrScoringEnabled(False)
    
    # Next, create a Genome that wraps the Flexibility Tree's motionNodes:
    motionNodes = [ftTrans, ftRotQuat, ftTorsion1, ftTorsion2]
    genome = Genome(motionNodes, ftRoot, scorer)
    
    # Next, create an Individual based on this Genome:
    individual = Individual(genome)
    
    # First test: we constructed the Individual without specifying initial
    # gene values, which is equivalent to calling initialize() on it.
    # Hence the FTDiscreteTranslation will be set to the identity translation;
    # the FTRotationAboutPointQuat will do no rotation (though it will still
    # translate the rotPoint to the origin); the FTTorsion nodes will leave
    # the origAngles unchanged. Let's score to verify this is true:
    score = individual.calculateScores()
    genes = individual.getGenesPy()
    print  "After initialize(), the genes are %s, while the AtomSet scores %f against itself and is now as follows (coords should be like origCoords, except for translation of rotPoint to the origin):" % (genes,score)
    atomSet.printDebugDescription()
    
    # Second test: let's set the genes to (arbitrary) new values, and score again:
    arbGenes = numpy.array([0.4666666667, 0.6, 0.6, 0.575, 0.625, 0.7031009601, 0.9330127019, 0.6031591627, 0.0816952176], FTYPE)
    assert(len(arbGenes) == individual.getNumVariables(), "test.py: testFlexibilityTree(): Adjust code if needed; length of arbGenes no longer matches individual's number of motion-control variables (= %d)." % (individual.getNumVariables(),))
    individual.setGenes(arbGenes)
    score = individual.calculateScores()
    genes = individual.getGenesPy()  # should equal arbGenes
    print "After an arbitrary transformation, the genes are %s, while the AtomSet scores %f against itself and is now:" % (genes,score)
    atomSet.printDebugDescription()
    
    # Third test: let's randomize the genes, and score again:
    Randomize.setRandSeed(25)
    individual.randomize()
    score = individual.calculateScores()
    genes = individual.getGenesPy()
    print "After randomize(), the genes are %s, while the AtomSet scores %f against itself and is now:" % (genes,score)
    atomSet.printDebugDescription()
    
    # Fourth test: Starting from the above randomized state, let's do
    # a SolisWets minimization:
    sw = SolisWets.getSingleton()
    sw.setEnableDebugOutput(True)
    sw.minimize(individual, 2, 2, 100, 10, 0.3)
    score = individual.getTotalScore()  # calculateScores() isn't needed after sw.minimize()
    genes = individual.getGenesPy()
    print "After SolisWets minimization, the genes are %s, while the AtomSet scores %f against itself and is now:" % (genes,score)
    atomSet.printDebugDescription()
    
    # NEW: As an (artificial) test of FTDiscreteConformation and
    # FTDiscreteTransformation, we're going to do surgery on the Flexibility Tree,
    # removing ftTorsion1 and ftAtoms_Branch, and replacing them with
    # an FTDiscreteConformation and FTDiscreteTransformation,
    # which we will set up so as to replicate two of the above test conformations.
    
    # Step 1: Start surgery on the Flexibility Tree. Recall that it is currently
    # ftRoot >> ftTrans >> ftRotQuat >> ftAtoms_Body >> ftTorsion1 >> ftAtoms_Branch >>
    # ftTorsion2 >> ftAtoms_SubBranch. After we finish our surgery (Step 4 below), it will be:
    # ftRoot >> ftTrans >> ftRotQuat >> ftAtoms_Body >> ftConformation >> ftTransformation >>
    # ftTorsion2 >> ftAtoms_SubBranch:
    ftAtoms_Branch.removeChild(ftTorsion2)
    ftAtoms_Body.removeChild(ftTorsion1)
    ftTorsion1 = None
    ftAtoms_Branch = None
    
    # Step 2: Create ftConformation, an FTDiscreteConformation node holding two possible
    # conformations for the "Branch" atoms (which are in fact a single atom, atom 5).
    # The first-conformation coords for atom 5 have been set to match
    # the result of applying ftTorsion1's localMatrix from the 'arbitrary transformation'
    # test above to atom 5's orig coords; the second-conformation coords have similarly
    # been set to match the 'after initialize()' test above:
    ftConformation = FTDiscreteConformation("DiscreteConformation")
    ftAtoms_Body.addChild(ftConformation)
    ftConformation.setAtomIndexes(numpy.array([5], 'i'))
    ftConformation.setCoords(numpy.array([[[0.512936,1.695214,1.770053]], [[1.111,2.001,1.503]]], FTYPE))
    # Uncomment to test the feature that an FTDiscreteConformation may specify
    # per-conformation offsets to be added to the computed energy of the node's atomSet:
    #ftConformation.setEnergyOffsets(numpy.array([-100., 200.], FTYPE))
    
    # Step 3: Create ftTransformation, an FTDiscreteTransformation node holding
    # two matrices, coordinating with ftConformation's choice of two conformations.
    # The first matrix has been set to match ftTorsion1's localMatrix from the
    # 'arbitrary transformation' test above; the second matrix has similarly been set
    # to match ftTorsion1's localMatrix from the 'after initialize()' test above:
    ftTransformation = FTDiscreteTransformation("DiscreteTransformation")
    ftConformation.addChild(ftTransformation)
    ftTransformation.setMatricesData(numpy.array(
        [[[0.897516,0.411659,0.158121,-1.54559],
          [-0.438396,0.871701,0.218966,0.10888],
          [-0.0476951,-0.265846,0.962835,0.907859],
          [0,0,0,1]],
         [[1,0,0,0],
          [0,1,0,0],
          [0,0,1,0],
          [0,0,0,1]]], FTYPE))
    
    # Step 4: Finishing our surgery on the Flexibility Tree:
    ftTransformation.addChild(ftTorsion2)
    
    # Step 5: Re-create genome and individual to match the modified Flexibility Tree.
    # Note that numMotionNodes and numVariables haven't changed:
    motionNodes[2] = ftConformation  # has replaced ftTorsion1
    genome = Genome(motionNodes, ftRoot, scorer)
    individual = Individual(genome)
    assert(len(arbGenes) == individual.getNumVariables(), "test.py: testFlexibilityTree(): Code assumes that numVariables hasn't changed; but in fact, numVariables appears to have changed from %d to %d." % (len(arbGenes), individual.getNumVariables()));
    
    # Step 6: First FTDiscreteConformation test, whose results should match
    # (up to round-off error) the 'arbitrary transformation' test above.
    # We do this by setting the same arbGenes on the individual, except that
    # arbGenes[7] (previously the ftTorsion1 angle) is now the choice of conformation
    # (gene-value 0. corresponds to conformationChoice = 1):
    arbGenes[7] = 0.
    individual.setGenes(arbGenes)
    score = individual.calculateScores()
    genes = individual.getGenesPy()  # should equal arbGenes
    print "Replicating the 'arbitrary transformation' test above using FTDiscreteConformation. The genes are %s, while the AtomSet scores %f against itself and is now:" % (genes,score)
    atomSet.printDebugDescription()
    
    # Step 7: Second FTDiscreteConformation test, whose results should match
    # (up to round-off error) the 'after initialize()' test above.
    # We do this by calling initialize() on the individual, then changing
    # gene [7] to 1. (which translates to conformationChoice = 2):
    individual.initialize()
    tempGenes = individual.getGenesPy()
    tempGenes[7] = 1.
    individual.setGenes(tempGenes)
    score = individual.calculateScores()
    genes = individual.getGenesPy()  # should equal tempGenes
    print "Replicating the 'after initialize()' test above using FTDiscreteConformation. The genes are %s, while the AtomSet scores %f against itself and is now:" % (genes,score)
    atomSet.printDebugDescription()

def masterTest():
    print "Starting ADFRcc test from Python."
    #testParameters()
    #testAtomSet()
    #testPairwiseScorer()
    #testGridMaps()
    #testGridScorer()
    #testFlexibleReceptorScorer()
    #testMixedScoring()
    testFlexibilityTree()
    print "Completed ADFRcc test from Python."

if __name__ == "__main__":
    masterTest()
