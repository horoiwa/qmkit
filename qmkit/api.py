from qmkit.base import BaseInterface
from qmkit.geometry import generate_conformers
from qmkit import util


class GeometoryOptimizer(BaseInterface):

    def __init__(self, mol, save_dir=None, n_jobs=1):
        super().__init__(mol, save_dir, n_jobs)

        self.n_gen = 1000

        self.rms = 1.0

        self.forcefield = ""

        self.n_qm = 5

    def run(self):
        mol = self.mol
        tmpfile = self.tmpfile + ".sdf"

        mol, confIds = generate_conformers(mol, self.n_gen,
                                           self.rms, self.n_jobs)

        self.logger.info(f"{len(confIds)} conformers generated")
        util.to_sdf_by_confIds(mol, confIds, tmpfile)
        mols = util.from_multisdf(tmpfile)
        print(len(mols))

