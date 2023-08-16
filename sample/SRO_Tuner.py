import numpy as np
import os

from RHEA import RHEA
from SROTuner import SRO
from SROTuner import Quaternary
from SROTuner import Ternary
np.set_printoptions(suppress=True)


def testSRO_Tuner(testFile: str, WCP: np, element: int, structure: str, tolerance=20):
    baseDir = os.getcwd()
    outputStr = os.path.join(baseDir, "tmp.lmp")
    match element:
        case 3: tuner = Ternary(testFile, WCP, outputStr, structure, tolerance) 
        case 4: tuner = Quaternary(testFile, WCP, outputStr, structure, tolerance)
    tuner.adjustSRO()
    WCP2 = SRO(outputStr, structure, element).WCP()
    print(WCP2)
    # print(WCP)
    os.remove(outputStr)
    return



def testTuner_BCC_Ternary_Equiatomic():
    baseDir = os.getcwd()
    inputFile = os.path.join(baseDir, "SRO_Tuner", "BCC", "ternary", "MoNbTa_Random.lmp")
    WCP  = np.array([[ 0.64082326,  0.06680182, -0.70809454],
                     [ 0.06680182, -0.40116626,  0.3343155 ],
                     [-0.70809454,  0.3343155,   0.37429779]])
    testSRO_Tuner(inputFile, WCP, 3, "BCC")


def testTuner_BCC_Quaternary_Equiatomic():
    baseDir = os.getcwd()
    inputFile = os.path.join(baseDir, "SRO_Tuner", "BCC", "quaternary", "Mo25Nb25Ta25W25_Random.lmp")
    WCP = RHEA().getWCP_300K()["25_25_25_25"]
    testSRO_Tuner(inputFile, WCP, 4, "BCC")


def testTuner_FCC_Ternary_Equiatomic():
    baseDir = os.getcwd()
    inputFile = os.path.join(baseDir, "SRO_Tuner", "FCC", "ternary", "FeCoNi_random.lmp")
    WCP  = np.array([[ 0.3,  0.1, -0.4],
                     [ 0.1, -0.2,  0.1],
                     [-0.4,  0.1,  0.3]])
    testSRO_Tuner(inputFile, WCP, 3, "FCC")


def testTuner_FCC_Quaternary_Equiatomic():
    baseDir = os.getcwd()
    inputFile = os.path.join(baseDir, "SRO_Tuner", "FCC", "quaternary", "FeCoCrNi_Random.lmp")
    WCP  = np.array([[0.1, 	-0.1,  -0.1,	0.1],
                     [-0.1,	-0.2, 	0.2,    0.1],
                     [-0.1,	 0.2,	0.1,   -0.2],
                     [0.1,   0.1,  -0.2,	0  ]])
    testSRO_Tuner(inputFile, WCP, 4, "FCC")




testTuner_BCC_Ternary_Equiatomic()
testTuner_BCC_Quaternary_Equiatomic()
testTuner_FCC_Ternary_Equiatomic()
testTuner_FCC_Quaternary_Equiatomic()


