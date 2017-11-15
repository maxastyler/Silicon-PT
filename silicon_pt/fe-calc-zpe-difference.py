import numpy as np
import matplotlib.pyplot as plt

def load_zpes(lats):
    return [get_zpe(i) for i in lats]

def get_zpe(lat_string):
    with open("si.{lat}/si.{lat}.fe".format(lat=lat_string)) as f:
        return(float(f.readlines()[0].split()[1]))

def get_electrons_energy(lat_string):
    with open("si.{lat}/si.scf.{lat}.out".format(lat=lat_string)) as f:
        for line in f.readlines(): 
            if "!" in line: 
                return float(line.split()[-2])

def load_electron_energies(lats):
    return [get_electrons_energy(i) for i in lats]

lats = np.arange(10.28, 10.38, 0.002)
lat_strings = []
for i in lats: 
    lat_strings.append("{lat:.3f}".format(lat=i))
#plot_electron_energy(lats, lat_strings)
zpes=load_zpes(lat_strings)
es=load_electron_energies(lat_strings)
e_sums=[es[i]+zpes[i] for i in range(len(es))]
fig, ax1 = plt.subplots()
x1 = plt.subplots()
ax1.plot(lats, es, 'b-')
ax1.set_ylabel('Free energy with no ZPE', color='b')
ax2 = ax1.twinx()
ax2.plot(lats, e_sums, 'r-')
ax2.set_ylabel('Free energy with ZPE', color='r')
plt.show()
