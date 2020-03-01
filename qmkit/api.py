from qmkit.base import BaseInterface
from qmkit.geometry import generate_conformers
from qmkit import util


class GeometoryOptimizer(BaseInterface):

    def __init__(self, mol, save_dir=None, n_jobs=1):

        super().__init__(mol, save_dir=save_dir, config="OptConfig")

        self.show_config()

    def run(self, n_jobs=1):
        mol = self.mol
        tmpfile = self.tmpfile + "_confs.sdf"

        self.logger.info("Start conformer generation")
        mol, mmenergy_confIds = generate_conformers(mol,
                                                    self.config["confgen"],
                                                    self.config["rms"],
                                                    n_jobs)

        self.logger.info(f"{len(mmenergy_confIds)} conformers generated")

        confIds = [mmenergy_confIds[i][1] for i in range(self.config["n_qm"])]
        util.to_sdf_by_confIds(mol, confIds, tmpfile)
        mols = util.from_multisdf(tmpfile)

        self.logger.info(f"QM-optimization: {len(mols)} conformers ")
