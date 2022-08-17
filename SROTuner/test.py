import numpy as np
import os

from src.checker import SRO
from src.tuner import Quaternary
from src.tuner import Ternary
np.set_printoptions(suppress=True)

def testSRO_Check(testFile: str, WCP: np, element: int, structure: str,  messges: str, tolerance = 0.05) -> bool:
    alloy = SRO(testFile, structure, element)
    calc_WCP = alloy.WCP()
    # print(calc_WCP)
    diff = np.abs(calc_WCP - WCP)
    diff[diff <= tolerance] = 0
    assert np.sum(diff) <= 0, f"Do not pass test: the calculated WCPs is different with the given WCPs: {messges}"
    print(f"Pass the test of WCPs (SRO) calculation for {messges}")


def testSRO_Tuner(testFile: str, WCP: np, element: int, structure: str, tolerance=20):
    baseDir = os.getcwd()
    outputStr = os.path.join(baseDir, "tmp.lmp")
    match element:
        case 3: tuner = Ternary(testFile, WCP, outputStr, structure, tolerance) 
        case 4: tuner = Quaternary(testFile, WCP, outputStr, structure, tolerance) 
    tuner.adjustSRO()
    WCP2 = SRO(outputStr, structure, element).WCP()
    # print(WCP2)
    # print(WCP)
    os.remove(outputStr)
    return


def testSRO_BCC_Ternary_Equiatomic_Random():
    baseDir = os.getcwd()
    testSRO_Check_File = os.path.join(baseDir, "SROTuner", "testData", "SRO_Check", "BCC", "ternary", "MoNbTa_Random.lmp")
    Actual_WCP  = np.zeros((3, 3))
    testSRO_Check(testSRO_Check_File, Actual_WCP, 3, "BCC", "BCC Ternary Equiatomic Random")


def testSRO_BCC_Ternary_Equiatomic_Order():
    baseDir = os.getcwd()
    testSRO_Check_File = os.path.join(baseDir, "SROTuner", "testData", "SRO_Check", "BCC", "ternary", "MoNbTa_Order.lmp")
    Actual_WCP  = np.array([[ 0.64082326,  0.06680182, -0.70809454],
                            [ 0.06680182, -0.40116626,  0.3343155 ],
                            [-0.70809454,  0.3343155,   0.37429779]])
    testSRO_Check(testSRO_Check_File, Actual_WCP, 3, "BCC", "BCC Ternary Equiatomic Order")


def testSRO_BCC_Quaternary_Equiatomic_Random():
    baseDir = os.getcwd()
    testSRO_Check_File = os.path.join(baseDir, "SROTuner", "testData", "SRO_Check", "BCC", "quaternary", "Mo25Nb25Ta25W25_Random.dump")
    Actual_WCP  = np.zeros((4, 4))
    testSRO_Check(testSRO_Check_File, Actual_WCP, 4, "BCC", "BCC Quaternary Equiatomic Random")


def testSRO_BCC_Quaternary_Equiatomic_Order():
    baseDir = os.getcwd()
    testSRO_Check_File = os.path.join(baseDir, "SROTuner", "testData", "SRO_Check", "BCC", "quaternary", "Mo25Nb25Ta25W25_Order.dump")
    Actual_WCP  = np.array([[0.845214844, 	-0.266357422,	-0.943603516,	0.364746094],
                            [-0.266357422,	-0.043212891, 	0.470703125,	-0.161132813],
                            [-0.943359375,	0.471435547,	0.708984375,   	-0.237060547],
                            [0.364746094,  	-0.161132813,	-0.237060547,	0.033447266]])
    testSRO_Check(testSRO_Check_File, Actual_WCP, 4, "BCC", "BCC Quaternary Equiatomic Order")


def testSRO_FCC_Ternary_Equiatomic_Random():
    baseDir = os.getcwd()
    testSRO_Check_File = os.path.join(baseDir, "SROTuner", "testData", "SRO_Check", "FCC", "ternary", "FeCoNi_random.lmp")
    Actual_WCP  = np.zeros((3, 3))
    testSRO_Check(testSRO_Check_File, Actual_WCP, 3, "FCC", "FCC Ternary Equiatomic Random")

def testSRO_FCC_Quaternary_Equiatomic_Random():
    baseDir = os.getcwd()
    testSRO_Check_File = os.path.join(baseDir, "SROTuner", "testData", "SRO_Check", "FCC", "quaternary", "FeCoCrNi_Random.lmp")
    Actual_WCP  = np.zeros((4, 4))
    testSRO_Check(testSRO_Check_File, Actual_WCP, 4, "FCC", "FCC Quaternary Equiatomic Random")


def testTuner_BCC_Ternary_Equiatomic():
    baseDir = os.getcwd()
    inputFile = os.path.join(baseDir, "SROTuner", "testData", "SRO_Tune", "BCC", "ternary", "MoNbTa_Random.lmp")
    WCP  = np.array([[ 0.64082326,  0.06680182, -0.70809454],
                     [ 0.06680182, -0.40116626,  0.3343155 ],
                     [-0.70809454,  0.3343155,   0.37429779]])
    testSRO_Tuner(inputFile, WCP, 3, "BCC")


def testTuner_BCC_Quaternary_Equiatomic():
    baseDir = os.getcwd()
    inputFile = os.path.join(baseDir, "SROTuner", "testData", "SRO_Tune", "BCC", "quaternary", "Mo25Nb25Ta25W25_Random.lmp")
    WCP  = np.array([[0.845214844, 	-0.266357422,	-0.943603516,	0.364746094],
                     [-0.266357422,	-0.043212891, 	0.470703125,	-0.161132813],
                     [-0.943359375,	0.471435547,	0.708984375,   	-0.237060547],
                     [0.364746094,  -0.161132813,	-0.237060547,	0.033447266]])
    testSRO_Tuner(inputFile, WCP, 4, "BCC")

testSRO_BCC_Ternary_Equiatomic_Random()
testSRO_BCC_Ternary_Equiatomic_Order()
testSRO_BCC_Quaternary_Equiatomic_Random()
testSRO_BCC_Quaternary_Equiatomic_Order()
testTuner_BCC_Ternary_Equiatomic()
testTuner_BCC_Quaternary_Equiatomic()

testSRO_FCC_Ternary_Equiatomic_Random()
testSRO_FCC_Quaternary_Equiatomic_Random()

