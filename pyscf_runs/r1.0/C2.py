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

mycas = mcscf.CASSCF(myhf, 26, 8)
mycas.fcisolver = shci.SHCI(mol)
mycas.fcisolver.useExtraSymm = False
mycas.fcisolver.sweep_iter = [ 0  ]
mycas.fcisolver.sweep_epsilon = [ 1e9 ]
mycas.max_cycle_macro = 0
e_noPT = mycas.mc2step()[0]
exit(0)
