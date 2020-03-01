from rdkit import Chem
import sys
from rdkit import Chem
from rdkit.Chem import AllChem

from concurrent import futures

from qmkit.base import BaseInterface
from qmkit.api import GeometoryOptimizer



if __name__ == "__main__":
    mol = Chem.MolFromSmiles("CC(=O)OCCC(/C)=C\C[C@H](C(C)=C)CCC=C")
    optimizer = GeometoryOptimizer(mol, n_jobs=2)
    optimizer.n_gen = 100
    optimizer.run()
