from typing import Tuple
from ovito.io import import_file
from ovito.data import NearestNeighborFinder
import numpy as np
from tqdm import tqdm

class SRO:
    def __init__(self, dataFile: str, structure: str, elements: int, composition = None) -> None:
        pipeline = import_file(dataFile)
        self.data = pipeline.compute()
        self.eleTypes = self.data.particles.particle_types
        self.elements = elements
        self.composition = self.getComposition() if (not composition) else composition
        assert self.elements == len(self.composition), f"the given composition shows {len(self.composition)} elements\n, but ovito detects {self.elements} elements"
        self.SRO_WCP = [[] for _ in range(self.elements**2)]
        assert (structure == "BCC" or structure == 'FCC'), "only support BCC and FCC structure"
        self.neighCount = 8 if structure == 'BCC' else 12
        
    # get the composition (np)
    def getComposition(self) -> np.ndarray:
        tmp = np.array([0 for _ in range(self.elements)])
        for i in self.eleTypes:
            tmp[i-1] += 1
        return tmp/self.data.particles.count
    
    # calculate the WCP
    def WCP(self) -> np.ndarray:
        finder = NearestNeighborFinder(self.neighCount, self.data)
        for idx in tqdm(range(self.data.particles.count)):
            eleType, counts = self.localAtomicEnv(idx, finder)
            self.single_atom_WCP(eleType, self.neighCount, counts)
        SRO_para = self.reshapeWCP()
        return SRO_para
        
    def reshapeWCP(self) -> np.ndarray:
        SRO_para = []
        for sro in self.SRO_WCP:
            tmp = sum(sro)/len(sro)
            SRO_para.append(tmp)
        SRO_para = np.array(SRO_para)
        SRO_para = SRO_para.reshape(self.elements, self.elements)
        return SRO_para

    # get the type of central atom and count the number of each element in first shell
    def localAtomicEnv(self, idx: int, finder: NearestNeighborFinder) -> Tuple[int, list[int]]:
        eleType = self.eleTypes[idx]
        counts  = [0 for _ in range(self.elements)]
        for neigh in finder.find(idx):
            neighType = self.eleTypes[neigh.index]
            counts[neighType-1] += 1
        return eleType, counts
    
    # calculate the local WCP for a specific atom
    def single_atom_WCP(self, eleType:int, N:int, counts:list[int]) -> None:
        base = (eleType-1)*self.elements
        i = 0
        for num, comp in zip(counts, self.composition):             
            tmp_SRO  = 1 - num/(comp*N)
            self.SRO_WCP[base+i].append(tmp_SRO)
            i += 1
        return

