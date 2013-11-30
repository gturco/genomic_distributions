import unittest
import sys
from flatfeature import Bed
from pyfasta import Fasta
#sys.path.append("../scripts")
from perfect_target_region import rel_pos

class TestPerfectTargetRegion(unittest.TestCase):
    def setUp(self):
        self.gene_name = "Os01g02110"
        self.bed = Bed("ricetest.bed")
        self.fasta = Fasta("ricetest.fasta")
        self.gene = self.bed.accn(self.gene_name)
        self.exons = self.gene['locs']


    def test_rel_pos(self):

        self.assertEqual((376,486),rel_pos(self.gene,self.exons[0]))
        self.assertEqual((1289,1789),rel_pos(self.gene,self.exons[-1]))

    def test_fasta(self):
        exon = self.exons[-1]
        seq = self.fasta[self.gene_name][:]
        self.assertTrue(1789 <= len(seq))

    def test_pattern(self):
        e = exons[-1]
        start, stop = rel_pos(self.gene,e)
        for exon in self.exons:

if __name__ == '__main__':
    unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSequenceFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)
