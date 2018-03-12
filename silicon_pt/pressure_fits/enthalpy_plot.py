import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def load_pressure(T):
    with open("si.{}.pressure".format(T)) as f:
        lines=f.readlines()[9:]
        pressures = [float(i.split()[4]) for i in lines]
        enthalpies = [float(i.split()[2]) for i in lines]
    return list(reversed(pressures)), list(reversed(enthalpies))
p1, e1 = load_pressure(0)
p2, e2 = load_pressure(600)
print(p1, e1)

plt.plot(p1, e1)
plt.show()
