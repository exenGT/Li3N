from catlearn.optimize.mlneb import MLNEB

from gpaw import GPAW, PW, FermiDirac
from ase.io import read
from ase.parallel import parprint
from ase.optimize.lbfgs import LBFGS
from ase.constraints import StrainFilter
from gpaw.xc.libvdwxc import vdw_df_cx
from ase.constraints import FixAtoms
import copy


calc = GPAW(mode=PW(450), kpts=(5,5,1), xc='PBE',
            convergence={'eigenstates': 1.e-7, 'density': 2.0e-4},
            random=True, occupations=FermiDirac(0.05), txt='neb_opt.txt',
            parallel={'augment_grids': True})


neb_catlearn = MLNEB(start='initial_opt.traj', # Initial end-point.
                     end='final_opt.traj', # Final end-point.
                     ase_calc=calc, # Calculator, it must be the same as the one used for the optimizations.
                     n_images=7, # Number of images (interger or float, see above).
                     interpolation='idpp', # Choose between linear or idpp interpolation (as implemented in ASE).
                     restart=False)

neb_catlearn.run(fmax=0.05, trajectory='ML-NEB.traj', full_output=True)
