import uuid
from pathlib import Path
from logging import DEBUG, Formatter, StreamHandler, getLogger

from rdkit import Chem


def get_logger():
    logger = getLogger(str(uuid.uuid4()))
    logger.setLevel(DEBUG)
    stream_handler = StreamHandler()
    handler_format = Formatter(
        '(%(levelname)s)[%(asctime)s]\n%(message)s')
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)
    return logger


def to_sdf_by_confIds(mol, confIds, filename):
    writer = Chem.SDWriter(filename)
    for confID in confIds:
        writer.write(mol, confId=confID)
    writer.close()


def from_multisdf(filepath):
    if not Path(filepath).exists():
        raise FileNotFoundError(f"{str(filepath)}")

    suppl = Chem.SDMolSupplier(filepath)
    mols = [mol for mol in suppl if mol is not None]
    mols = list(map(Chem.AddHs, mols))
    return mols
