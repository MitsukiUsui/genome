from unittest import TestCase

class Test_relative_overlap(TestCase):
    def test_count_relative_overlap(self):
        from orfcoverage.relative_overlap import count_relative_overlap
        from orfcoverage.relative_overlap import Interval

        cdsInterval_lst = [Interval("cds1", 3, 9, "-0")]
        orfInterval_lst = [Interval("orf1", 0, 6, "+0"),
                           Interval("orf2", 3, 9, "-0"),
                           Interval("orf3", 7, 10, "-1")]
        relativeOverlap_dct = count_relative_overlap(cdsInterval_lst, orfInterval_lst)

        self.assertEqual(relativeOverlap_dct["+0"], 6)
        self.assertEqual(relativeOverlap_dct["+1"], 0)
        self.assertEqual(relativeOverlap_dct["+2"], 2)
        self.assertEqual(relativeOverlap_dct["-0"], 3)
        self.assertEqual(relativeOverlap_dct["-1"], 0)
        self.assertEqual(relativeOverlap_dct["-2"], 0)
        self.assertEqual(relativeOverlap_dct["intergenic"], 4)


