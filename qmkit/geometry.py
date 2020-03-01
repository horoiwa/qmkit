from tqdm import tqdm
import psi4
from rdkit.Chem import AllChem


def generate_conformers(mol, confgen=1000, rms=2.0, n_jobs=1):
    """ コンフォーマ生成
        confgen個の３D conformerを生成後、MMFFでエネルギー計算
    """
    conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=confgen,
                                          pruneRmsThresh=rms,
                                          numThreads=n_jobs)

    energy = []
    prop = AllChem.MMFFGetMoleculeProperties(mol)
    for cid in tqdm(conf_ids):
        mmff = AllChem.MMFFGetMoleculeForceField(mol, prop, confId=cid)
        mmff.Minimize()
        energy.append((mmff.CalcEnergy(), cid))

    # 4. エネルギーをソートし，相対エネルギーとIDをリストに格納
    energy.sort()
    mmenergy_conf_ids = [(i-energy[0][0], j) for i, j in energy]

    return mol, mmenergy_conf_ids


def mm_optimize(mol):
    return mol


def qm_optimize(mol):
    return mol


def td_dft(mol, n):
    log = ""
    return log
