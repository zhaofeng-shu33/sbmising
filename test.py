import unittest
from sbmising import sbm_graph, SIBM
from utility import get_ground_truth, compare

class TestGenerator(unittest.TestCase):
    def test_sbm(self):
        with self.assertRaises(ValueError):
            sbm_graph(50, 2, 16, 4)
        with self.assertRaises(ValueError):
            sbm_graph(100, 2, 3, 4)
        with self.assertRaises(ValueError):
            sbm_graph(101, 2, 16, 4)

class TestIsing(unittest.TestCase):
    def test_ising_2(self):
        G = sbm_graph(100, 2, 16, 4)
        results = SIBM(G, k=2)
        print(results)
        labels_true = get_ground_truth(G)
        self.assertAlmostEqual(compare(results, labels_true), 1.0)        


if __name__ == '__main__':
    unittest.main()