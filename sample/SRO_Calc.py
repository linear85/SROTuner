from SROTuner import SRO
from RHEA import RHEA
import os
import numpy as np

def testSRO_Check(testFile: str, WCP: np, element: int, structure: str,  messges: str, tolerance = 0.05) -> bool:
    alloy = SRO(testFile, structure, element)
    calc_WCP = alloy.WCP()
    diff = np.abs(calc_WCP - WCP)
    diff[diff <= tolerance] = 0
    assert np.sum(diff) <= 0, f"Do not pass test: the calculated WCPs is different with the given WCPs: {messges}"
    print(f"Pass the test of WCPs (SRO) calculation for {messges}")


def testSRO_BCC_Ternary_Equiatomic_Random():
    baseDir = os.path.dirname(os.path.realpath(__file__))
    testSRO_Check_File = os.path.join(baseDir, "SRO_Calc", "BCC", "ternary", "MoNbTa_Random.lmp")
    Actual_WCP  = np.zeros((3, 3))
    testSRO_Check(testSRO_Check_File, Actual_WCP, 3, "BCC", "BCC Ternary Equiatomic Random")


def testSRO_BCC_Ternary_Equiatomic_Order():
    baseDir = os.path.dirname(os.path.realpath(__file__))
    testSRO_Check_File = os.path.join(baseDir, "SRO_Calc", "BCC", "ternary", "MoNbTa_Order.lmp")
    Actual_WCP  = np.array([[ 0.64082326,  0.06680182, -0.70809454],
                            [ 0.06680182, -0.40116626,  0.3343155 ],
                            [-0.70809454,  0.3343155,   0.37429779]])
    testSRO_Check(testSRO_Check_File, Actual_WCP, 3, "BCC", "BCC Ternary Equiatomic Order")


def testSRO_BCC_Quaternary_Equiatomic_Random():
    baseDir = os.path.dirname(os.path.realpath(__file__))
    testSRO_Check_File = os.path.join(baseDir, "SRO_Calc", "BCC", "quaternary", "Mo25Nb25Ta25W25_Random.dump")
    Actual_WCP  = np.zeros((4, 4))
    testSRO_Check(testSRO_Check_File, Actual_WCP, 4, "BCC", "BCC Quaternary Equiatomic Random")


def testSRO_BCC_Quaternary_Equiatomic_Order():
    baseDir = os.path.dirname(os.path.realpath(__file__))
    testSRO_Check_File = os.path.join(baseDir, "SRO_Calc", "BCC", "quaternary", "Mo25Nb25Ta25W25_Order.dump")
    Actual_WCP = np.array(RHEA().getWCP_300K()["25_25_25_25"])
    testSRO_Check(testSRO_Check_File, Actual_WCP, 4, "BCC", "BCC Quaternary Equiatomic Order")


def testSRO_FCC_Ternary_Equiatomic_Random():
    baseDir = os.path.dirname(os.path.realpath(__file__))
    testSRO_Check_File = os.path.join(baseDir, "SRO_Calc", "FCC", "ternary", "FeCoNi_random.lmp")
    Actual_WCP  = np.zeros((3, 3))
    testSRO_Check(testSRO_Check_File, Actual_WCP, 3, "FCC", "FCC Ternary Equiatomic Random")


def testSRO_FCC_Ternary_Equiatomic_Order():
    baseDir = os.path.dirname(os.path.realpath(__file__))
    testSRO_Check_File = os.path.join(baseDir, "SRO_Calc", "FCC", "ternary", "FeCoNi_order.lmp")
    Actual_WCP = np.array([[ 0.3,  0.1, -0.4],
                           [ 0.1, -0.2,  0.1],
                           [-0.4,  0.1,  0.3]])
    testSRO_Check(testSRO_Check_File, Actual_WCP, 3, "FCC", "FCC Ternary Equiatomic Order")


def testSRO_FCC_Quaternary_Equiatomic_Random():
    baseDir = os.path.dirname(os.path.realpath(__file__))
    testSRO_Check_File = os.path.join(baseDir, "SRO_Calc", "FCC", "quaternary", "FeCoCrNi_Random.lmp")
    Actual_WCP  = np.zeros((4, 4))
    testSRO_Check(testSRO_Check_File, Actual_WCP, 4, "FCC", "FCC Quaternary Equiatomic Random")


def testSRO_FCC_Quaternary_Equiatomic_Order():
    baseDir = os.path.dirname(os.path.realpath(__file__))
    testSRO_Check_File = os.path.join(baseDir, "SRO_Calc", "FCC", "quaternary", "FeCoCrNi_order.lmp")
    Actual_WCP  = np.array([[0.1, 	-0.1,  -0.1,	0.1],
                            [-0.1,	-0.2, 	0.2,    0.1],
                            [-0.1,	 0.2,	0.1,   -0.2],
                            [0.1,    0.1,  -0.2,	0]])
    testSRO_Check(testSRO_Check_File, Actual_WCP, 4, "FCC", "FCC Quaternary Equiatomic Order")


testSRO_BCC_Ternary_Equiatomic_Random()
testSRO_BCC_Ternary_Equiatomic_Order()
testSRO_BCC_Quaternary_Equiatomic_Random()
testSRO_BCC_Quaternary_Equiatomic_Order()
testSRO_FCC_Ternary_Equiatomic_Random()
testSRO_FCC_Ternary_Equiatomic_Order()
testSRO_FCC_Quaternary_Equiatomic_Random()
testSRO_FCC_Quaternary_Equiatomic_Order()