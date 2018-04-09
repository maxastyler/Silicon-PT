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
energies = [load_energy_ecut(ecut) for ecut in ecuts]
pressures = [load_pressure_ecut(ecut) for ecut in ecuts]

plt.plot(ecuts, [abs(e-energies[-1]) for e in energies])
#plt.plot(ecuts, [abs(p-pressures[-1]) for p in pressures])
plt.show()

ks = list(filter(lambda x: x<200, sorted(map(lambda x: int(x.split('.')[-2]), glob("*.k.*.out")))))
energies = [load_energy_k(k) for k in ks]
pressures = [load_pressure_k(k) for k in ks]

plt.plot(ks, [abs(e-energies[-1]) for e in energies])
#plt.plot(ks, [abs(p-pressures[-1]) for p in pressures])
plt.show()
