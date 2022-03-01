from gpaw import GPAW, PW, FermiDirac
from ase.io import read
from ase.parallel import parprint
from ase.optimize.bfgs import BFGS
from ase.constraints import StrainFilter
from gpaw.xc.libvdwxc import vdw_df_cx


#atoms = read("Li3N.cif")

atoms = read("Li3N_cell_opt.traj")

calc = GPAW(mode=PW(500), kpts=(11,11,11), xc='PBE',
            convergence={'eigenstates': 1.e-8},
            random=True, occupations=FermiDirac(0.01), txt='Li3N_opt_2.txt')

atoms.calc = calc

sf = StrainFilter(atoms)
relax = BFGS(sf, logfile='Li3N_cell_opt_2.log', trajectory='Li3N_cell_opt_2.traj')
relax.run(fmax=0.02)
