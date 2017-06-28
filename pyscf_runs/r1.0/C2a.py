#!/usr/bin/python
import numpy
import math
from pyscf import gto, scf, ao2mo, mcscf, tools, fci
from pyscf.future.shciscf import shci, settings

mol = gto.M(
    atom = 'C 0 0 0; C 0 0 1.0',
    basis = 'cc-pvdz',
    verbose=4,
    symmetry=1,
    spin = 0)
myhf = scf.RHF(mol)
myhf.kernel()
myhf.analyze()
