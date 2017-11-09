import numpy as np
import matplotlib.pyplot as plt

lats = np.arange(10.2, 10.37, 0.01)
fe_lat_temp = []
for i in lats: 
    fe_str = "si.{lat:.2f}/si.{lat:.2f}.fe".format(lat=i)
    fe_lat_temp.append(np.loadtxt(fe_str).transpose().tolist())

LAT, TEMP = np.meshgrid(lats, fe_lat_temp[0][0])
FREE_ENERGIES = np.array([fe_lat_temp[i][1] for i in range(len(lats))]).transpose()
print(np.argmin(FREE_ENERGIES, 1))

plt.contourf(LAT, TEMP, FREE_ENERGIES, 100)
print(lats[14])
plt.show()
