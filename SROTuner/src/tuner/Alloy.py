from abc import abstractclassmethod
from tuner.SRO_Data import SRO_Data
import random
import pandas as pd
import numpy as np
import os
import re
import time
np.set_printoptions(suppress=True)

# This class is to construct atomic configure from given SRO level (WCPs) for quaternary alloys
class Alloy():
    def __init__(self, strcuture_file:str, WCPs:np, saved_path:str|bool, structure: str, tolerance = 30):
        '''
        Args:
            strcuture_file (string):                the lammps data about the alloy
            WCPs (np.array):                        the Warren-Cowley parameters
            saved_path (string or False):           the path and filename to save the created atomic configuration / Do not save
            structure (str):                        the atomic number in the first shell (indicate the crystal structure, for example BCC: 8)
            tolerance (int):                        the tolerance (difference between real and want WCP) to finish the tuning
        '''
        assert (structure == "BCC" or structure == 'FCC'), "only support BCC and FCC structure"
        self.N           = 8 if structure == "BCC" else 12
        self.readFile    = strcuture_file
        # self.data        = SRO_Data(strcuture_file, self.N)
        self.data        = self.create_SRO_Data()
        print("Finish reading")
        self.totalNumber = self.data.totalNumber
        self.range       = [i for i in range(self.totalNumber)]
        self.comp        = self.data.composition
        self.wantWCPs    = self.getWantWCPs(WCPs)
        self.inputPath   = strcuture_file
        self.statistics()
        print("Finish random assign")
        self.curWCPs     = self.getCurWCPs()
        # print(self.curWCPs)
        # print(self.wantWCPs)
        self.name = saved_path
        self.tolerance = tolerance
    
    @abstractclassmethod
    def create_SRO_Data(self):
        pass
    
    @abstractclassmethod
    def getWantWCPs(self):
        pass

    @abstractclassmethod
    def getCurWCPs(self):
        pass

    # main function to randomly assign the atoms and count the distribution
    @abstractclassmethod
    def statistics(self):
        pass

    # Helper function to randomly assign the atoms (randomly select the sites for different elements)
    @abstractclassmethod
    def randomAssign(self):
        pass

    # Helper function to count the distribution of the given atom and type
    @abstractclassmethod
    def count(self, index:int, atom_type:int):
        pass
    
    # remove the type of atom from neighbor atom
    @abstractclassmethod
    def removeAtom(self, index:int, prev_type:int):
        pass

    # add the type of atom from neighbor atom
    @abstractclassmethod
    def addAtom(self, index:int, next_type:int):
        pass

    # calculate the change of WCPs due to the swap of the two atoms
    @abstractclassmethod
    def getChange(self, atom_1:int, atom_2:int) -> np:
        pass

    # main function to adjust the SRO level and save it
    def adjustSRO(self) -> None:
        print("Start tuning")
        t_D = self.wantWCPs - self.curWCPs
        while True:
            if (self.checkDiff(np.copy(t_D))):
                print("Done with tuning")
                if (not self.name):
                    return
                self.__toTxt()
                self.__addPreifx()
                return
            atom_1, atom_2 = self.randomTwoAtoms()
            l_D = self.getChange(atom_1, atom_2)
            if (self.checkAccept(np.copy(t_D), np.copy(l_D))):
                t_D = t_D - l_D
                # print(l_D)
                self.swapAtoms(atom_1, atom_2)
                tmp_sum = np.sum(np.absolute(t_D))
                print(int(tmp_sum), end='\r')
                
                # print(t_D)
    
    # Helper function to check the different between actual SRO and idea SRO is smaller than the tolerance or not
    def checkDiff(self, diff:np):
        diff = np.abs(diff)
        diff[diff <= self.tolerance] = 0
        return np.sum(diff) <= 0
    
    # Helper function to swap two atoms and update the distributions
    def swapAtoms(self, atom_1:int, atom_2:int) -> None:
        t1, t2 = self.data.df_All["type"].iloc[atom_1], self.data.df_All["type"].iloc[atom_2]
        for nei in self.data.graph[atom_1]:
            self.update_Neighbor(nei, t1, t2)
        for nei in self.data.graph[atom_2]:
            self.update_Neighbor(nei, t2, t1)
        self.data.df_All["type"].iloc[atom_1] = t2
        self.data.df_All["type"].iloc[atom_2] = t1
        return

    # Helper function to update the distribution of neighbor atom
    def update_Neighbor(self, index:int, prev_type:int, next_type:int) -> None:
        self.removeAtom(index, prev_type)
        self.addAtom(index, next_type)

    # random select two atoms with different types
    def randomTwoAtoms(self):
        while True:
            atom_1 = random.choice(self.range)
            atom_2 = random.choice(self.range)
            if (atom_2 != atom_1 and self.data.df_All['type'].iloc[atom_1] != self.data.df_All['type'].iloc[atom_2]):
                return atom_1, atom_2


    
    # check the diff is on the right direction or not
    @staticmethod
    def checkAccept(t_D:np, l_D:np) -> bool:
        prev = np.sum(np.absolute(t_D))
        late = np.sum(np.absolute(t_D - l_D))
        if (random.uniform(0, 1) < 0.0001):
             return late <= prev+6
        return late <= prev
    
    # write the atomic configuration into a temp file
    def __toTxt(self):
        saved_df = self.data.df_ind
        saved_df[["type"]] = self.data.df_All[["type"]]
        saved_df[['x', 'y', 'z']] = self.data.df_pos
        self.tmp_path ='modified.dump'
        saved_df.to_csv(self.tmp_path, header=None, index=None, sep=' ')
    
    # transform the temp file into a lammps readable file (add the prefix into it) and save it, and delete the temp file
    def __addPreifx(self):
        data_1 = self.__getPrefix()
        with open(self.tmp_path) as fp:
            data2 = fp.read()
        data_1 += "\n"
        data_1 += "\n"
        data_1 += data2
        with open (self.name, 'w') as fp:
            fp.write(data_1)
        os.remove(self.tmp_path)
    
    # get the header (first part) of the dump file
    def __getPrefix(self) -> str:
        f_in = open(self.inputPath, 'r')
        script = f_in.read()
        f_in.close()
        index = "Atoms # atomic"
        idx = script.index(index)
        return script[:idx+len(index)]




