from prody.tests import TestCase

from numpy import array, log, zeros, char
from numpy.testing import assert_array_equal, assert_array_almost_equal

from prody.tests.datafiles import *

from prody import LOGGER, refineMSA, parseMSA, calcMSAOccupancy, mergeMSA
from prody import uniqueSequences

LOGGER.verbosity = None

FASTA = parseMSA(pathDatafile('msa_Cys_knot.fasta'))
FASTA_ALPHA = char.isalpha(FASTA._msa)
NUMSEQ = FASTA.numSequences() * 1.

class TestRefinement(TestCase):


    def testLabel(self):
        label = 'FSHB_BOVIN'
        index = FASTA.getIndex(label)
        refined = refineMSA(FASTA, label=label)._getArray()

        expected = FASTA._getArray().take(FASTA_ALPHA[index].nonzero()[0], 1)

        assert_array_equal(refined, expected)


    def testRowocc(self):

        refined = refineMSA(FASTA, rowocc=0.9)._getArray()
        expected = FASTA._getArray()[FASTA_ALPHA.sum(1) / 112. >= 0.9,:]

        assert_array_equal(refined, expected)

    def testColocc(self):

        refined = refineMSA(FASTA, colocc=0.9)._getArray()
        expected = FASTA._getArray()[:, FASTA_ALPHA.sum(0) / NUMSEQ >= 0.9]

        assert_array_equal(refined, expected)

    def testRowCol(self):

        rowocc = 0.9
        colocc = 1.0
        refined = refineMSA(FASTA, rowocc=rowocc,
                            colocc=colocc)._getArray()
        rows = FASTA_ALPHA.sum(1) / 112. >= rowocc
        expected = FASTA._getArray()[rows]
        cols = char.isalpha(expected).sum(0,
                    dtype=float) / expected.shape[0] >= colocc

        expected = expected.take(cols.nonzero()[0], 1)
        assert_array_equal(refined, expected)

    def testSeqid(self):

        seqid = 0.50
        label = 'FSHB_BOVIN'
        unique = uniqueSequences(FASTA, seqid)
        refined = refineMSA(FASTA, seqid=seqid)
        assert_array_equal(refined._getArray(), FASTA._getArray()[unique])

    def testSeqidLabel(self):

        seqid = 0.50
        label = 'FSHB_BOVIN'
        labeled = refineMSA(FASTA, label=label)
        unique = uniqueSequences(labeled, seqid)
        unique[FASTA.getIndex(label)] = True
        refined = refineMSA(FASTA, label=label, seqid=seqid)
        assert_array_equal(refined._getArray(), labeled._getArray()[unique])


    def testAll(self):

        rowocc = 0.9
        colocc = 0.9
        seqid = 0.98
        label = 'FSHB_BOVIN'
        refined = refineMSA(FASTA, label=label, seqid=seqid,
                            rowocc=rowocc, colocc=colocc)

        index = FASTA.getIndex(label)
        which = FASTA_ALPHA[index].nonzero()[0]
        expected = FASTA._getArray().take(which, 1)

        expected = expected[uniqueSequences(expected, seqid)]

        expected = expected[calcMSAOccupancy(expected, 'row') >= rowocc]

        which = (calcMSAOccupancy(expected) >= colocc).nonzero()[0]
        expected = expected.take(which, 1)

        assert_array_equal(refined._getArray(), expected)


class TestMerging(TestCase):


    def testMerge(self):

        merged = mergeMSA(FASTA, FASTA)
        length = FASTA.numResidues()
        self.assertEqual(merged[:, :length], merged[:, length:])
