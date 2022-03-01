from gpaw import GPAW, PW, FermiDirac
from ase.io import read
from ase.parallel import parprint
from ase.optimize.lbfgs import LBFGS
from ase.constraints import StrainFilter
from gpaw.xc.libvdwxc import vdw_df_cx
from ase.constraints import FixAtoms


atoms = read("initial.cif")

# atoms = read("GDY_Li_top_opt.traj")

calc = GPAW(mode=PW(500), kpts=(5,5,1), xc='PBE',
            convergence={'eigenstates': 1.e-7, 'density': 1.0e-4},
            random=True, occupations=FermiDirac(0.05), txt='initial_opt.txt',
            parallel={'augment_grids': True})

atoms.calc = calc

Z = atoms.get_positions()[:, 2]
indices = [i for i, z in enumerate(Z) if z < 8.5]
constraint = FixAtoms(indices=indices)
atoms.set_constraint(constraint)

#sf = StrainFilter(atoms)
relax = LBFGS(atoms=atoms, logfile='initial_opt.log', trajectory='initial_opt.traj')
relax.run(fmax=0.05)

