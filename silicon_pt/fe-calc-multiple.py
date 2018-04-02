#TODO: need to include electron energy
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D

ADD_ELECTRON_ENERGY=True

#This way of calculating gives the free energies as FREE_ENERGIES[x][y]
#The first index, x specifies the temperature and the second specifies the lattice parameter

#Load lattice parameters by a list of strings which are the lattice parameter eg. 10.210
def load_lats(lats, lat_strings):
    fe_lat_temp = []
    electron_energy=[]
    for i in lat_strings:
        fe_str = "si.{lat}/si.{lat}.fe".format(lat=i)
        f_es_np = np.loadtxt(fe_str).transpose()
        if f_es_np.ndim==1: print(i)
        f_es = f_es_np.tolist()
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

    #calculate the fitted free energies to give: fitted_energies[temp][lat]
    fitted_energies=[]
    for temp_index in range(len(LAT)):
        p_fit = curve_fit(lambda x, a, b, c: a*(x-b)**2+c, LAT[temp_index], FREE_ENERGIES[temp_index])
        a=p_fit[0]
        if temp_index==0:
            print(a)
        fitted_energies.append([a[0]*(x-a[1])**2+a[2] for x in LAT[temp_index]])

    fig, ax1 = plt.subplots()
    #x1 = plt.subplots()
    #plt.title(r'Free Energy vs Lattice Parameter in Silicon at T=0K, T=500K')
    ax1.plot(LAT[0], FREE_ENERGIES[0], 'b-')
    ax1.plot(LAT[0], fitted_energies[0], 'b.')
    ax1.set_xlabel('Lattice Paramter (Bohr)')
    # Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel(r'Free energy at $T=0K$', color='b')
    ax1.tick_params('y', colors='b')
    
    ax2 = ax1.twinx()
    ax2.plot(LAT[25], FREE_ENERGIES[25], 'r-')
    ax2.plot(LAT[25], fitted_energies[25], 'r.')
    ax2.set_ylabel(r'Free energy at $T=500K$', color='r')
    ax2.tick_params('y', colors='r')
    
    fig.tight_layout()
    fig.set_size_inches(6, 5)
    fig.savefig('../../Sodium-DFT-Project/project_presentation/silicon_0_500_comparison.png')
    #plt.show()


def plot_electron_energy(lats, lat_strings):
    es = [get_electrons_energy(i) for i in lat_strings]
    plt.plot(lats, es) 
    plt.show()

def get_electrons_energy(lat_string):
    with open("si.{lat}/si.scf.{lat}.out".format(lat=lat_string)) as f:
        for line in f.readlines(): 
            if "!" in line: 
                return float(line.split()[-2])

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

lats = np.arange(10.28, 10.38, 0.002)
lat_strings = []
for i in lats: 
    lat_strings.append("{lat:.3f}".format(lat=i))
#plot_electron_energy(lats, lat_strings)
load_lats(lats, lat_strings)
