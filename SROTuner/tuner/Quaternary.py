from tuner.Alloy import Alloy
import random
import numpy as np
from tqdm import tqdm

class Quaternary(Alloy):
    def __init__(self, strcuture_file: str, WCPs: np, saved_path: str or bool, structure: str, tolerance=30):
        super().__init__(strcuture_file, WCPs, saved_path, structure, tolerance)
    
    def getWantWCPs(self, WCPs:np) -> np:
        compArray = np.array([[self.comp[0], self.comp[1], self.comp[2], self.comp[3]],
                              [self.comp[0], self.comp[1], self.comp[2], self.comp[3]],
                              [self.comp[0], self.comp[1], self.comp[2], self.comp[3]],
                              [self.comp[0], self.comp[1], self.comp[2], self.comp[3]]]) * self.N
        WCPs = (1 - WCPs) * compArray
        WCPs[0, :] *= self.comp[0]*self.totalNumber
        WCPs[1, :] *= self.comp[1]*self.totalNumber 
        WCPs[2, :] *= self.comp[2]*self.totalNumber 
        WCPs[3, :] *= self.comp[3]*self.totalNumber
        return WCPs
    
    def getCurWCPs(self) -> np:
        WCPs    = np.zeros((4,4))
        df1     = self.data.df_All.loc[self.data.df_All["type"].isin([1])]
        df2     = self.data.df_All.loc[self.data.df_All["type"].isin([2])]
        df3     = self.data.df_All.loc[self.data.df_All["type"].isin([3])]
        df4     = self.data.df_All.loc[self.data.df_All["type"].isin([4])]
        WCPs[0] = df1.sum().tolist()[1:]
        WCPs[1] = df2.sum().tolist()[1:]
        WCPs[2] = df3.sum().tolist()[1:]
        WCPs[3] = df4.sum().tolist()[1:]
        return WCPs
    
    def statistics(self):
        self.randomAssign()
        for index in tqdm(range(self.totalNumber)):
            if (index in self.E1):
                self.count(index, 1)
            elif (index in self.E2):
                self.count(index, 2)
            elif (index in self.E3):
                self.count(index, 3)
            else:
                self.count(index, 4)
        print("Done with creating dataframe")
        return
    
    def randomAssign(self):
        self.E_all  = set([i for i in range(self.totalNumber)])
        self.E1     = set(random.sample([i for i in range(self.totalNumber)], int(self.totalNumber*self.comp[0])))
        self.E_all  = self.E_all - self.E1
        self.E2     = set(random.sample(list(self.E_all), int(self.totalNumber*self.comp[1])))
        self.E_all  = self.E_all - self.E2
        self.E3     = set(random.sample(list(self.E_all), int(self.totalNumber*self.comp[2])))
        self.E4     = self.E_all - self.E3
        print("Done with random selection")
    
    def count(self, index:int, atom_type:int):
        count_1, count_2, count_3, count_4 = 0, 0, 0, 0
        for nei in self.data.graph[index]:
            if (nei in self.E1):
                count_1 += 1
            if (nei in self.E2):
                count_2 += 1
            if (nei in self.E3):
                count_3 += 1
            if (nei in self.E4):
                count_4 += 1
        self.data.df_All["n_1"].iloc[index]  = count_1
        self.data.df_All["n_2"].iloc[index]  = count_2
        self.data.df_All["n_3"].iloc[index]  = count_3
        self.data.df_All["n_4"].iloc[index]  = count_4
        self.data.df_All["type"].iloc[index] = atom_type
        return


    def removeAtom(self, index:int, prev_type:int):
        match prev_type:
            case 1: self.data.df_All["n_1"].iloc[index] -= 1
            case 2: self.data.df_All["n_2"].iloc[index] -= 1
            case 3: self.data.df_All["n_3"].iloc[index] -= 1
            case 4: self.data.df_All["n_4"].iloc[index] -= 1
            case _: raise "This is a Quaternary alloy, only have four elements"
        return
    
    # add the type of atom from neighbor atom
    def addAtom(self, index:int, next_type:int):
        match next_type:
            case 1: self.data.df_All["n_1"].iloc[index] += 1
            case 2: self.data.df_All["n_2"].iloc[index] += 1
            case 3: self.data.df_All["n_3"].iloc[index] += 1
            case 4: self.data.df_All["n_4"].iloc[index] += 1
            case _: raise "This is a Quaternary alloy, only have four elements"
        return
   
   # calculate the change of WCPs due to the swap of the two atoms
    def getChange(self, atom_1:int, atom_2:int) -> np:
        diff   = np.zeros((4,4))
        L1     = self.data.df_All.iloc[atom_1].tolist()
        L2     = self.data.df_All.iloc[atom_2].tolist()
        t1, t2 = L1.pop(0)-1, L2.pop(0)-1
        d1, d2 = np.array(L2) - np.array(L1), np.array(L1) - np.array(L2)
        diff[t1, :] += d1
        diff[:, t1] += d1
        diff[t2, :] += d2
        diff[:, t2] += d2
        return diff