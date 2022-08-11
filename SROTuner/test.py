import numpy as np

from checker import Alloy
np.set_printoptions(suppress=True)

testFile = "testData/Mo25Nb25Ta25W25_Random.dump"
# testFile = "testData/Mo25Nb25Ta25W25_Order.dump"
RHEA25 = Alloy(testFile, "BCC")
print(RHEA25.WCP())


