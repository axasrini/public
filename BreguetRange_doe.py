"""
Battery kWh vs range for 21700 NCM cell 
Varying Cl/Cd for different ground effects
WTO based on half a ton aircraft with varying Wpay
Average efficiencies for propulsion
"""

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

CL_CD_ratio = np.array(
    [15, 20, 25, 30, 35, 40, 45, 50], dtype=float
)

MTOW = 500;  #max taekoff weight init  
Wbatt_ratio = np.array(        #nasa electric x57 is 0.28 https://www.nasa.gov/centers/armstrong/news/FactSheets/FS-109.html
	[0.05, 0.1, 0.2], dtype=float
)

few = 0.4    # baseline from https://link.springer.com/article/10.1007/s13272-021-00530-w
""""""
np.array(
	[0.35,0.4,0.5], dtype=float
)
empty_weight=few*MTOW
cb=540000  #J/kg for good density LFP battery pack 150Wh/kg
range_electric = np.empty(shape=(8,3))
""""""
# Other constants
AIR_DENSITY = 1.1671
g= 9.81
n_i = 0.96
n_m = 0.94
n_p = 0.9

class nf(float):
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()

fmt = '%r'

def calc_range_breguet():
    """
    Function calculates Range from breguet equation
    :param Wbatt_ratio: ratio of battery pack weight to take off weight of plane
    :return:numpy array with electric range
    """
    for i in range(len(CL_CD_ratio)):
        CL_CD = CL_CD_ratio[i]
        for j in range(len(Wbatt_ratio)):
            Wbatt_Wto = Wbatt_ratio[j]
            range_electric[i, j] = cb / g * CL_CD * Wbatt_Wto * n_i * n_m * n_p / 1000

    return range_electric


if __name__ == '__main__':
    calc_range_breguet()
    print(range_electric)

    #assuming a half ton aircraft calculate the payload and battery estimates based on Wbatt_Wto ratio
    Wbatt = Wbatt_ratio * MTOW
    Wpay = MTOW - few*MTOW - Wbatt*1.1  #assuming it takes 10% more for take off and climb
    batt_cap = cb * Wbatt / 3600 /1000   #kWh
    print("Payloads are ", Wpay)
    print("Battery weights are ", Wbatt)
    print("Empty weight is ", few*MTOW)

#plot some stuff

#z3 = griddata( CL_CD_ratio, Wbatt, (xi, yi), method='linear')
plt.figure()
levels_range = np.linspace(0,500,10)
c3 = plt.contour(Wpay, CL_CD_ratio, range_electric, levels=levels_range, linewidths=1.2,colors='k',label='Cl/Cd')
c3.levels = [nf(val) for val in c3.levels]
plt.clabel(c3,c3.levels,inline=True,fmt=fmt,fontsize=10)
plt.xlabel(r'Payload (kg)')
plt.ylabel(r'Cl_Cd ratio')
plt.annotate('Range in km ',(1000,100),color='k')
plt.title('Battery Range curve for seaflight')

plt.figure()
levels_range = np.linspace(0,500,10)
c4 = plt.contour(batt_cap*1.1, CL_CD_ratio, range_electric, levels=levels_range, linewidths=1.2,colors='k',label='Cl/Cd')
c4.levels = [nf(val) for val in c4.levels]
plt.clabel(c4,c4.levels,inline=True,fmt=fmt,fontsize=10)
plt.xlabel(r'Battery capacity (kWh)')
plt.ylabel(r'Cl_Cd ratio')
plt.annotate('Range in km ',(1000,100),color='k')
plt.title('Battery Range curve for seaflight')
plt.show()


