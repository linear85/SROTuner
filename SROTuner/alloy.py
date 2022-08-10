from ovito.io import import_file
import numpy as np

class alloy:
    def __init__(self, dataFile: str,  composition = None) -> None:
        pipeline = import_file(dataFile)
        self.data = pipeline.compute()
        self.eleTypes = self.data.particles.particle_types
        self.elements = self.eleTypes.type
        self.composition = self.getComposition() if (not composition) else composition
        assert self.elements == len(self.composition), f"the given composition shows {len(self.composition)} elements\n, but ovito detects {self.elements} elements"
        

    # get the composition (np)
    def getComposition(self) -> np:
        tmp = np.array([0 for _ in range(self.elements)])
        for i in self.eleTypes:
            tmp[i-1] += 1
        return tmp/self.data.particles.count
    
    # calculate the WCP
    def WCP(self):
        return
    
    def localAtomicEnv(self, idx: int, finder: object):
        eleType = self.eleTypes[idx]
        counts  = [0 for _ in range(self.elements)]
        for neigh in finder.find(idx):
            neighType = self.eleTypes[neigh.index]
            counts[neighType-1] += 1
        return eleType, counts
    
    # calculate the local WCP for a specific atom
    def single_atom_WCP(self, eleType:int, N:int, counts:list[int], SRO: np) -> None:
        base = (eleType-1)*eleType
        i = 0
        for num, comp in zip(counts, self.composition):             
            tmp_SRO  = 1 - num/(comp*N)
            SRO[base+i].append(tmp_SRO)
            i += 1
        return SRO

    

