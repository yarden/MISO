
import unittest
import pysplicing

class TestMISO(unittest.TestCase):
    def test_miso1(self):
        gene=pysplicing.createGene( ((1,100), (201,300), (401,500)),
                                    ((0,1), (0,2), (0,1,2)) )
        pysplicing.noIso(gene)
        pysplicing.isoLength(gene)
        reads=pysplicing.simulateReads(gene, 0L, (0.2,0.3,0.5), 2000L, 33L)
        est=pysplicing.MISO(gene, 0L, reads[1], reads[2], 33L, 5000L, 500L, 
                            10L, (1.0,1.0,1.0))
        [ sum(e)/len(e) for e in est[0] ]

if __name__ == '__main__':
    unittest.main()
