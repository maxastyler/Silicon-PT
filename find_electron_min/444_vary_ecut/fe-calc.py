#TODO: need to include electron energy
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit

ADD_ELECTRON_ENERGY=False

#This way of calculating gives the free energies as FREE_ENERGIES[x][y]
#The first index, x specifies the temperature and the second specifies the lattice parameter

#Load lattice parameters by a list of strings which are the lattice parameter eg. 10.210
def load_lats(lats, lat_strings):
    fe_lat_temp = []
    electron_energy=[]
    for i in lat_strings:
        fe_str = "si.{lat}/si.{lat}.fe".format(lat=i)
        f_es = np.loadtxt(fe_str).transpose().tolist()
        if ADD_ELECTRON_ENERGY: 
            lat_energy = get_electrons_energy(i)
            electron_energy.append(lat_energy)
            for j in range(len(f_es[1])):
                f_es[1][j] += lat_energy
        fe_lat_temp.append(f_es)
    
    LAT, TEMP = np.meshgrid(lats, fe_lat_temp[0][0])
    FREE_ENERGIES = np.array([fe_lat_temp[i][1] for i in range(len(lats))]).transpose()
    #fig = plt.figure()
    #ax = fig.gca(projection = '3d')
    #surf = ax.plot_surface(LAT, TEMP, FREE_ENERGIES)
     
    plt.plot(LAT[0], FREE_ENERGIES[0])
    #plt.plot(LAT[0], electron_energy)
#    plt.contourf(LAT, TEMP, FREE_ENERGIES, 100)
    plt.show()

def plot_electron_energy(lats, lat_strings):
    es = [get_electrons_energy(i) for i in lat_strings]
    plt.plot(lats, es) 
    plt.show()

def plot_electron_energy_difference(lats, lat_strings):
    es = [get_electrons_energy(i) for i in lat_strings]
    e_diffs = [abs(i - es[-1]) for i in es]
    plt.semilogy(lats, e_diffs) 
    plt.show()

def plot_electron_pressure(lats, lat_strings):
    ps = [get_electrons_pressure(i) for i in lat_strings]
    plt.plot(lats, ps) 
    plt.show()

def plot_electron_pressure_difference(lats, lat_strings):
    ps = [get_electrons_pressure(i) for i in lat_strings]
    p_diffs = [abs(i - ps[-1]) for i in ps]
    plt.semilogy(lats, p_diffs) 
    plt.show()

def plot_electron_energy_pressure_difference(lats, lat_strings):
    es = [get_electrons_energy(i) for i in lat_strings]
    e_diffs = [abs(i - es[-1]) for i in es]
    plt.subplot(121)
    plt.title("Energy convergence vs ecutwfc")
    plt.semilogy(lats, e_diffs) 
    ps = [get_electrons_pressure(i) for i in lat_strings]
    p_diffs = [abs(i - ps[-1]) for i in ps]
    plt.subplot(122)
    plt.title("Pressure convergence vs ecutwfc")
    plt.semilogy(lats, p_diffs) 
    plt.show()
def get_electrons_energy(lat_string):
    with open("si.{lat}/si.scf.{lat}.out".format(lat=lat_string)) as f:
        for line in f.readlines(): 
            if "!" in line: 
                return float(line.split()[-2])

def get_electrons_pressure(lat_string):
    with open("si.{lat}/si.scf.{lat}.out".format(lat=lat_string)) as f:
        for line in f.readlines(): 
            if "total   stress  (Ry/bohr**3)                   (kbar)" in line:
                return float(line.split()[-1])

def load_closer_to_lat_param():
    lats = np.arange(10.32, 10.36, 0.001)
    lat_strings = []
    for i in lats:
        lat_strings.append("{lat:.3f}".format(lat=i))
    load_lats(lats, lat_strings)

def load_more_spaced_out():
    lats = np.arange(10.2, 10.37, 0.01)
    lat_strings = []
    for i in lats: 
        lat_strings.append("{lat:.2f}".format(lat=i))
    load_lats(lats, lat_strings)

#load_closer_to_lat_param()
#load_more_spaced_out()

ecuts = np.arange(10, 60, 2)
ecut_strings = []
for i in ecuts: 
    ecut_strings.append("{ecut}".format(ecut=i))
plot_electron_energy_pressure_difference(ecuts, ecut_strings)

