from gpaw import GPAW, PW, FermiDirac
from ase.io import read
from ase.parallel import parprint
from ase.optimize.lbfgs import LBFGS
from ase.constraints import StrainFilter
from gpaw.xc.libvdwxc import vdw_df_cx
from ase.constraints import FixAtoms
import numpy as np

atoms = read("Li3N_opt_110_opt_Li.cif")

constraints = FixAtoms(indices=np.ndarray.tolist(np.arange(0,124,1)))
atoms.set_constraint(constraints)

# atoms = read("GDY_Li_top_opt.traj")

calc = GPAW(mode=PW(450), kpts=(3,2,1), xc='PBE',
            convergence={'eigenstates': 1.e-8},
            random=True, occupations=FermiDirac(0.01), txt='Li3N_opt.txt')

atoms.calc = calc

#sf = StrainFilter(atoms)
relax = LBFGS(atoms=atoms, logfile='Li3N_cell_opt.log', trajectory='Li3N_cell_opt.traj')
relax.run(fmax=0.05)
