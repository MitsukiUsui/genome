from unittest import TestCase
from orfcoverage.phase import PhaseController

class Test_phase(TestCase):

    def test_operate(self):
        pc = PhaseController()
        self.assertEqual(pc.operate("+1", "-1"), "-2")

    def test_revopes(self):
        """
        this is the definition of revops
        """

        pc = PhaseController()
        for p in PhaseController.phase_lst:
            for ops in PhaseController.phase_lst:
                revops = pc.revops(ops)
                self.assertEqual(pc.operate(pc.operate(p, ops), revops), p)
