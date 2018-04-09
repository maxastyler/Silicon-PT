import matplotlib.pyplot as plt
import numpy as np
from glob import glob

def load_energy_ecut(ecut):
    with open("si.ecut.{}.out".format(ecut)) as f:
        lines = f.readlines()
        return float(next(filter(lambda x: "!" in x, lines)).split()[-2])

def load_pressure_ecut(ecut):
    with open("si.ecut.{}.out".format(ecut)) as f:
        lines = f.readlines()
        return float(next(filter(lambda x: "kbar" in x, lines)).split()[-1].strip('P='))

def load_energy_k(k):
    with open("si.k.{}.out".format(k)) as f:
        lines = f.readlines()
        return float(next(filter(lambda x: "!" in x, lines)).split()[-2])

def load_pressure_k(k):
    with open("si.k.{}.out".format(k)) as f:
        lines = f.readlines()
        return float(next(filter(lambda x: "kbar" in x, lines)).split()[-1].strip('P='))

ecuts = list(filter(lambda x: x<200, sorted(map(lambda x: int(x.split('.')[-2]), glob("*.ecut.*.out")))))
energies = [load_energy_ecut(ecut)*13605.693 for ecut in ecuts]
pressures = [load_pressure_ecut(ecut) for ecut in ecuts]

plt.xlabel(r'ecutwfc $(Ry)$')
plt.ylabel(r'$|\Delta E|$ $(meV)/atom$')
plt.semilogy(ecuts, [abs(e-energies[-1]) for e in energies], color='black')
plt.gcf().set_size_inches(4, 3.5)
plt.tight_layout()
plt.savefig('../../Sodium-DFT-Project/project_report/figures/silicon_convergence_e_ecut.svg')
#plt.show()
plt.clf()
plt.xlabel(r'ecutwfc $(Ry)$')
plt.ylabel(r'$|\Delta P|$ $(kbar)$')
plt.semilogy(ecuts, [abs(p-pressures[-1]) for p in pressures], color='black')
plt.gcf().set_size_inches(4, 3.5)
plt.tight_layout()
plt.savefig('../../Sodium-DFT-Project/project_report/figures/silicon_convergence_p_ecut.svg')
#plt.show()
plt.clf()

ks = list(filter(lambda x: x<200, sorted(map(lambda x: int(x.split('.')[-2]), glob("*.k.*.out")))))
energies = [load_energy_k(k)*13605.693 for k in ks]
pressures = [load_pressure_k(k) for k in ks]

plt.xlabel(r'k grid')
plt.ylabel(r'$|\Delta E|$ $(meV)/atom$')
plt.gcf().set_size_inches(4, 3.5)
plt.tight_layout()
plt.semilogy(ks, [abs(e-energies[-1]) for e in energies], color='black')
plt.savefig('../../Sodium-DFT-Project/project_report/figures/silicon_convergence_e_k.svg')
#plt.show()
plt.clf()
plt.xlabel(r'k grid')
plt.ylabel(r'$|\Delta P|$ $(kbar)$')
plt.semilogy(ks, [abs(p-pressures[-1]) for p in pressures], color='black')
plt.gcf().set_size_inches(4, 3.5)
plt.tight_layout()
plt.savefig('../../Sodium-DFT-Project/project_report/figures/silicon_convergence_p_k.svg')
#plt.show()
