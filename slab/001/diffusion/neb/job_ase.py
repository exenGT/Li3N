from ase.io import read
from ase.constraints import FixAtoms
from ase.neb import NEB
from ase.optimize.lbfgs import LBFGS
from gpaw.mpi import rank, size

from gpaw import GPAW, PW, FermiDirac

initial = read('initial_opt.traj')
final = read('final_opt.traj')

Z = initial.get_positions()[:, 2]
indices = [i for i, z in enumerate(Z) if z < 8.5]
constraint = FixAtoms(indices=indices)

n = size // 6      # number of cpu's per image
j = 1 + rank // n  # my image number
assert 6 * n == size

images = [initial]

for i in range(6):
    ranks = range(i * n, (i + 1) * n)
    image = initial.copy()

    if rank in ranks:

        calc = GPAW(mode=PW(500), kpts=(5, 5, 1), xc='PBE',
                    convergence={'eigenstates': 1.e-7, 'density': 1.0e-4},
                    occupations=FermiDirac(0.05),
                    parallel={'augment_grids': True},
                    txt=f'neb_{j}.txt',
                    communicator=ranks)

        image.calc = calc

    image.set_constraint(constraint)
    images.append(image)

images.append(final)

neb = NEB(images, parallel=True, climb=True)
neb.interpolate(method='idpp')

qn = LBFGS(neb, logfile='neb.log', trajectory='neb.traj')
qn.run(fmax=0.05)

