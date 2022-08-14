from ovito.io import import_file
from ovito.data import NearestNeighborFinder
import pandas as pd
from collections import defaultdict
import numpy as np

# This class to create some DataFrame and graph to save the information of the crystal structure
class SRO_Data():
    def __init__(self, file_path:str, N:int):
        pipeline            = import_file(file_path)
        data                = pipeline.compute()
        self.totalNumber    = data.particles.count
        self.composition    = self.getComposition(data)
        finder              = NearestNeighborFinder(N, data)
        self.graph          = self.__toGraph(finder)
        self.df_All         = self.__toDF_All()
        self.df_pos         = self.__position_DF(data)
        self.df_ind         = self.__index_DF(data)
                

    def __toGraph(self, finder):
        graph = defaultdict(set)
        for index in range(self.totalNumber):
            for neigh in finder.find(index):
                graph[index].add(neigh.index)
        return graph
    
    @staticmethod
    def __position_DF(data):
        df_pos = pd.DataFrame(data.particles.position, columns=['x', 'y', 'z'])
        return df_pos
    
    @staticmethod
    def __index_DF(data):
        df_index = pd.DataFrame(data.particles.identifiers, columns=['id'])
        return df_index

    def __toDF_All(self):
        all_data = [[0 for _ in range(self.elements+1)] for _ in range(self.totalNumber)]
        columns=["type", "n_1", "n_2", "n_3", "n_4", "n_5"]
        df_All   = pd.DataFrame(all_data, columns=columns[:self.elements+1])
        return df_All
    
    # test the atomic configuration is perfect or not
    def isPerfectCrystal(self):
        for i in range(self.totalNumber):
            for nei in self.graph[i]:
                assert(i in self.graph[nei])
        print("Pass the test!!!")
    
    def getComposition(self, data:object) -> np:
        ptypes = data.particles.particle_types
        self.elements = ptypes.type + 1
        tmp = np.array([0 for _ in range(self.elements)])
        for i in ptypes:
            tmp[i-1] += 1
        return tmp/self.totalNumber




