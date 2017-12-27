from unittest import TestCase
from Bio.Seq import Seq
import pandas as pd

class Test_orf_overlap(TestCase):
    def test_create_orf_df(self):
        from orfcoverage.orf_overlap import create_orf_df

        seq = Seq("NNNTAANNNTAGNNN")
        orf_df = create_orf_df(seq)
        print(orf_df)

        seq2 = seq.reverse_complement()
        orf_df = create_orf_df(seq2)
        print(orf_df)
        self.fail() #TODO use from pandas.util.testing import assert_frame_equal for testing

    def test_count_overlap(self):
        from orfcoverage.orf_overlap import count_overlap

        orf_df = pd.DataFrame({"lane": [0, 3], "start": [0, 3], "end": [6, 9]})
        seqLength = 12

        overlapCount_arr = count_overlap(orf_df, seqLength)

        self.assertEqual(overlapCount_arr[0], 3)
        self.assertEqual(overlapCount_arr[1], 6)
        self.assertEqual(overlapCount_arr[2], 3)
        for i in range(3, 7):
            self.assertEqual(overlapCount_arr[i], 0)
