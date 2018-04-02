#TODO: need to include electron energy
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D

ADD_ELECTRON_ENERGY=True

det_cell = np.linalg.det(np.array([[-0.5, 0, 0.5], [0, 0.5, 0.5], [-0.5, 0.5, 0]]))

def murnaghan(v, v0, e0, b0, db0):
    eta = ((v/v0)**(1/3))
    bovo=b0*v0
    bracket = eta**(-3*db0)/(db0-1)+1
    last = bovo/(db0-1)
    middle = (bovo*eta**3/db0)*bracket
    return e0+middle-last

def birch_murnaghan(v, v0, e0, b0, db0):
    eta = ((v/v0)**(1/3))
    b_first = ((eta**(-2)-1)**3)*db0
    b_last = ((eta**(-2)-1)**2)*(6-4*(eta**(-2)))
    return e0 + (9/16)*b0*v0*(b_first+b_last)

def vinet(v, v0, e0, b0, db0):
    eta = ((v/v0)**(1/3))
    bovo_eta = b0*v0/((db0-1)**2)
    exp_br = (3/2)*(db0-1)*(1-eta)
    other_br = 3*(db0-1)*(1-eta)-2
    return e0 + 4*bovo_eta + 2*bovo_eta*np.exp(exp_br)*other_br

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
    v_0 = [[], []]
    for (temp_index, _) in enumerate(LAT):
        print(temp_index)
        p_fit = curve_fit(vinet, [i**3*det_cell for i in LAT[temp_index]], FREE_ENERGIES[temp_index], p0=[300, -21, 0.005, 4])
        v_0[0].append(p_fit[0][0])
        v_0[1].append(p_fit[0][1])
        a=p_fit[0]
        #if temp_index==0:
        print(a)
        fitted_energies.append([vinet(x, *a) for x in LAT[temp_index]])

    #fig1 = plt.subplot(211)
    #plt.figure().set_size_inches(5, 5)
    plt.title(r"Free Energy vs Lattice Parameter in Silicon with Isotherms")
    for i in range(0, len(LAT), 2):
        plt.plot(LAT[i], FREE_ENERGIES[i], 'b-')
        #fig1.plot(LAT[0], fitted_energies[0], 'b.')
    #print(v_0[0])
    l_vals = [(v/det_cell)**(1/3) for v in v_0[0]]
    #print(l_vals)
    plt.plot([(v/det_cell)**(1/3) for v in v_0[0]], v_0[1], c = 'r')
    plt.xlabel(r'Lattice Paramter (Bohr)')
        # Make the y-axis label, ticks and tick labels match the line color.
    plt.ylabel(r'Free energy $(Ry)$')
    plt.tick_params('y')
    plt.tight_layout()
    plt.savefig('../../Sodium-DFT-Project/project_presentation/silicon_isotherms.png')
    #plt.subplot(212)
    plt.clf()
    #plt.figure().set_size_inches(4, 6)
    plt.title(r"Silicon Volume vs Temperature")
    ts = [i for i in range(0, 701, 20)]
    plt.plot(ts, v_0[0], c='r')
    plt.xlabel(r"Temp $(K)$")
    plt.ylabel(r"Volume (Bohr$^3$)")
    plt.tight_layout()
    plt.savefig('../../Sodium-DFT-Project/project_presentation/silicon_min_volume.png')
    #plt.subplot(313)
    #plt.plot(*mid_gradient(ts, v_0[0]))
    
    #ax2 = ax1.twinx()
    #ax2.plot(LAT[25], FREE_ENERGIES[25], 'r-')
    #ax2.plot(LAT[25], fitted_energies[25], 'r.')
    #ax2.set_ylabel('Free energy at T=500K', color='r')
    #ax2.tick_params('y', colors='r')
    #
    #fig.tight_layout()
    #plt.show()

def mid_gradient(x, y):
    mids = list(map(lambda a: (a[1]-a[0])/2+a[0], zip(x, x[1:])))
    dxs = list(map(lambda a: (a[1]-a[0])/2, zip(x, x[1:])))
    dys = list(map(lambda a: (a[1]-a[0])/2, zip(y, y[1:])))
    grads = list(map(lambda a: a[1]/a[0], zip(dxs, dys)))
    
    return (mids, grads)

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
