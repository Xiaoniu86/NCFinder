import unittest
import pandas as pd
from src.bidirectional_detect import bidirect_detect_chromosome

class TestBidirectionalDetect(unittest.TestCase):
    def test_bidirectional_detection(self):
        cl_df = pd.DataFrame({'chr': ['chrI'], 'strand': ['+'], 'start': [100], 'end': [200]})
        gene_df = pd.DataFrame({'chromosome': ['chrI'], 'strand': ['-'], 'gene_start': [150], 'gene_end': [250]})
        result = bidirect_detect_chromosome(cl_df, gene_df, "chrI", 50, 50)
        self.assertEqual(len(result), 1)

if __name__ == '__main__':
    unittest.main()
