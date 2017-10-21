from unittest import TestCase
from orfcoverage.phase import Phase

class Test_phase(TestCase):

    def test_operate(self):
        phase = Phase()
        self.assertEqual(phase.operate("+1", "-1"), "-2")

    def test_revopes(self):
        """
        this is the definition of revops
        """

        phase = Phase()
        for p in phase.phase_lst:
            for ops in phase.phase_lst:
                revops = phase.revops(ops)
                self.assertEqual(phase.operate(phase.operate(p, ops), revops), p)
