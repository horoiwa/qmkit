from rdkit import Chem
from rdkit.Chem import AllChem


def generate_conformers(mol, n_conf=1000, rms=1.0, n_jobs=1):
    conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=n_conf,
                                          pruneRmsThresh=rms,
                                          numThreads=n_jobs)
    return mol, list(conf_ids)


def mm_optimize(mol):
    return mol


def qm_optimize(mol):
    return mol


def td_dft(mol, n):
    log = ""
    return log
