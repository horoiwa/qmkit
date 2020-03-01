import uuid
from pathlib import Path
from dataclasses import asdict
from pprint import pprint

import rdkit
from rdkit import Chem

from qmkit.util import get_logger
from qmkit.config import CONFIG_LIST


class BaseInterface:

    def __init__(self, mol, save_dir=None, config=None):

        if not isinstance(mol, rdkit.Chem.rdchem.Mol):
            raise TypeError(f"mol must be rdchem.Mol object, given {type(mol)}")

        self.mol = mol

        if save_dir is None:
            self.savedir = (Path.home() / ".qmkit").resolve()
        elif Path(save_dir).exists():
            self.savedir = self.rename_savedir(save_dir).resolve()
        else:
            self.savedir = Path(save_dir).resolve()

        self.tmpfile = str(self.savedir / f"{uuid.uuid4()}")

        if config:
            self.config = asdict(CONFIG_LIST[config]())

        self.log = {"in": self.mol, "out": None}

        self.logger = get_logger()

        self._post_init()

    def _post_init(self):
        self.savedir.mkdir(parents=True, exist_ok=True)

        Chem.AddHs(self.mol)

        Chem.AllChem.Compute2DCoords(self.mol)

    def __del__(self):
        """一時ファイル置き場を使用しているなら削除
        """
        if self.savedir == (Path.home() / ".qmkit").resolve():
            name = Path(self.tmpfile).stem
            tmpfiles = list(self.savedir.glob(f"{name}*"))
            for file in tmpfiles:
                file.unlink()

    def show_config(self):
        pprint(self.config)

    def get_result(self):
        raise NotImplementedError()

    @staticmethod
    def rename_savedir(save_dir: Path):
        n = 1
        while True:
            tmp = Path(str(save_dir) + f"_{n}")
            if not tmp.exists():
                save_dir = tmp
                break
            else:
                n += 1

        return save_dir

    @classmethod
    def from_file(cls, file_path, project_dir):

        file_path = Path(file_path)
        if not file_path.exists():
            raise FileNotFoundError()

        mol = cls.mol_from_file(file_path)
        return cls(mol, project_dir)

    @staticmethod
    def mol_from_file(file_path):
        suffix = file_path.suffix
        if suffix == ".mol":
            raise NotImplementedError()
        elif suffix == ".mol2":
            raise NotImplementedError()
        elif suffix == ".sdf":
            raise NotImplementedError()
        else:
            raise NotImplementedError(f'Read {suffix} is NotImplemented')
