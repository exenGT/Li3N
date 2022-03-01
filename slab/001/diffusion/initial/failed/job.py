from gpaw import GPAW, PW, FermiDirac
from ase.io import read
from ase.parallel import parprint
from ase.optimize.lbfgs import LBFGS
from ase.constraints import StrainFilter
from gpaw.xc.libvdwxc import vdw_df_cx
from ase.constraints import FixAtoms


atoms = read("initial_opt_2.traj")

# atoms = read("GDY_Li_top_opt.traj")

calc = GPAW(mode=PW(500), kpts=(5,6,1), xc=vdw_df_cx(),
            convergence={'eigenstates': 1.e-7},
            random=True, occupations=FermiDirac(0.03), txt='initial_opt_3.txt',
            parallel={'augment_grids': True})

atoms.calc = calc

Z = atoms.get_positions()[:, 2]
indices = [i for i, z in enumerate(Z) if z < 3.0]
constraint = FixAtoms(indices=indices)
atoms.set_constraint(constraint)

#sf = StrainFilter(atoms)
relax = LBFGS(atoms=atoms, logfile='initial_opt_3.log', trajectory='initial_opt_3.traj')
relax.run(fmax=0.05)

