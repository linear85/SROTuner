from laegr import GNNData, TrainEdge
from typing import List, Set, Dict
from ovito.io import import_file
import torch
from torch import Tensor
from tqdm import tqdm
import random
import pandas as pd
import os


def create_output_model(input_szie: int) -> torch.nn:
    model = torch.nn.Sequential(
        torch.nn.Linear(input_szie, 512),
        torch.nn.ReLU(),
        torch.nn.Linear(512, 64),
        torch.nn.ReLU(),
        torch.nn.Linear(64, 1)
    )
    return model



class EdgeModel(torch.nn.Module):
    def __init__(self, output_size: int) -> None:
        super().__init__()
        self.layer = torch.nn.Sequential(
            torch.nn.Linear(3, 512),
            torch.nn.ReLU(),
            torch.nn.Linear(512, 64),
            torch.nn.ReLU(),
            torch.nn.Linear(64, output_size)
        )

    def forward(self, x: Tensor) -> Tensor:
        return self.layer(x)



class GNN_MC:
    def __init__(self, dump_path: str, model_path: str, output_path: str) -> None:
        self.graph_data = GNNData(dump_path, criteria=26, feature="Onehot3D").toGraph()
        self.connectivity = self.getConnectivity(self.graph_data.edge_index)
        self.model = self.readModel(model_path)
        self.cur_PE = self.model_prediction()

    def run(self, steps: int) -> None:
        for _ in tqdm(range(steps)):
            atom_1, atom_2 = self.randomTwoAtoms()
            self.swap(atom_1, atom_2)
            self.accept(atom_1, atom_2)
    
    def saveStructure(self, output_template: str, output_path: str) -> None:
        types = self.__toAtomicType()
        self.changeComp(output_template, output_path, types)

    @staticmethod
    def readModel(path: str) -> torch.nn:
        agent = TrainEdge("NNConv")
        edge_model_1 = EdgeModel(100*4)
        # edge_model_2 = EdgeModel(100*100)
        gnn_model = agent.build_gnn_model(4, 100, 1, [edge_model_1])
        output_model = create_output_model(100)
        model = agent.build_model(gnn_model, output_model)
        model.load_state_dict(torch.load(path))
        return model

    def model_prediction(self) -> torch.tensor:
        nodes = torch.tensor(self.graph_data.nodes, dtype=torch.float)
        edge_index = torch.tensor(self.graph_data.edge_index, dtype=torch.long)
        edge_index=edge_index.t().contiguous()
        edge_attr = torch.tensor(self.graph_data.edge_attr, dtype=torch.float)
        pred = self.model(nodes, edge_index, edge_attr)
        return torch.sum(pred)
        
    @staticmethod
    def getConnectivity(edge_index: List[List[int]]) -> Dict[int, Set[int]]:
        connectivity = {}
        for [x, y] in edge_index:
            if x in connectivity:
                connectivity[x].add(y)
            else:
                connectivity[x] = {y}
            if y in connectivity:
                connectivity[y].add(x)
            else:
                connectivity[y] = {x}
        return connectivity
    
    def accept(self, atom_1: int, atom_2: int) -> None:
        next_PE = self.model_prediction()
        if (next_PE < self.cur_PE):
            self.cur_PE = next_PE
        else:
            self.swap(atom_1, atom_2)
    
    def swap(self, atom_1: int, atom_2: int) -> None:
        self.graph_data.nodes[atom_1], self.graph_data.nodes[atom_2] = self.graph_data.nodes[atom_2], self.graph_data.nodes[atom_1]
    
    def randomTwoAtoms(self) -> None:
        while True:
            atom_1 = random.randint(0, len(self.connectivity))
            atom_2 = random.randint(0, len(self.connectivity))
            if (atom_2 != atom_1 and self.graph_data.nodes[atom_1] != self.graph_data.nodes[atom_2]):
                return atom_1, atom_2

    def __toAtomicType(self) -> List[int]:
        types = []
        for i in self.graph_data.nodes:
            cur_type = i.index(1)
            types.append(cur_type+1)
        return types

    def changeComp(self, template: str, output_path: str, types: List[int]) -> None:
        s = open(template, 'r')
        First_Part = ""
        while True:
            line = s.readline()
            if ("Masses" in line):
                for _ in range(5):
                    s.readline()
                continue
            else:
                First_Part = First_Part + line
            if ("Atoms # atomic" in line):
                First_Part = First_Part + '\n'
                break
        df = pd.read_csv(s, header=None, delimiter=' +', engine='python')
        s.close()
        df.iloc[:, 1] = types
        df.to_csv("tmp.lmp", header=None, index=None, sep=' ', mode='a')
        f_tmp = open("tmp.lmp", 'r')
        text_Second_part = f_tmp.read()
        f_tmp.close()
        f_write = open(output_path, 'w')
        text = First_Part + text_Second_part
        f_write.write(text)
        f_write.close()
        os.remove("tmp.lmp")


# dump_path = r"E:\GNN_data\data\HEA_PE_PerfectStructure\MD_Simulation\BCC_MoNbTaW\25_25_25_25_0K\mc.0.dump"
# output_templet = "random.data"
# model_path = r"C:\Users\89721\Desktop\L3_GNN_Model_AllData"
# output_path = "test.dump"
# steps = 10
# model = GNN_MC(dump_path, model_path, output_path)
# print(model.cur_PE)
# model.run(steps)
# print(model.cur_PE)
# model.saveStructure(output_templet, output_path)