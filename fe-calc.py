import numpy as np
import matplotlib.pyplot as plt

lats = np.arange(10.2, 10.5, 0.01)
l_str = lambda x: "si.{lat:.2f}".format(lat=x)
for i in lats: 
    print(l_str(i))
