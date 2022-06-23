"""
Battery kWh vs range for 21700 NCM cell 
Varying Cl/Cd for different ground effects
WTO based on half a ton aircraft with varying Wpay
Average efficiencies for propulsion
"""

from cmath import sqrt
from email.headerregistry import ContentDispositionHeader
from unittest.util import _MIN_COMMON_LEN
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

CL_CD_ratio = np.array(
    [15, 20, 25, 30, 35, 40, 45, 50], dtype=float
)

MTOW = 500;  #kg max taekoff weight init  
Wbatt_ratio = np.array(        #nasa electric x57 is 0.28 https://www.nasa.gov/centers/armstrong/news/FactSheets/FS-109.html
	[0.05, 0.1, 0.2], dtype=float
)

few = 0.5    # baseline from https://link.springer.com/article/10.1007/s13272-021-00530-w
""""""
np.array(
	[0.35,0.4,0.5], dtype=float
)
""""""
empty_weight=few*MTOW
cb=540000  #J/kg for good density LFP battery pack 150Wh/kg
range_electric = np.empty(shape=(8,3))

# Other constants
rho_sl = 1.1671
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
    Wpay = MTOW - few*MTOW - Wbatt
    batt_cap = cb * Wbatt / 3600 /1000   #kWh
    print("Payloads are ", Wpay)
    print("Battery weights are ", Wbatt)
    print("Empty weight is ", few*MTOW)

#plot some stuff

#z3 = griddata( CL_CD_ratio, Wbatt, (xi, yi), method='linear')
plt.figure()
levels_range = np.linspace(100,1000,10)
c3 = plt.contour(Wbatt, CL_CD_ratio, range_electric, levels=levels_range, linewidths=1.2,colors='k',label='Cl/Cd')
c3.levels = [nf(val) for val in c3.levels]
plt.clabel(c3,c3.levels,inline=True,fmt=fmt,fontsize=10)
plt.xlabel(r'Battery weight (kg)')
plt.ylabel(r'Cl_Cd ratio')
plt.annotate('Range in km ',(1000,100),color='k')
plt.title('Battery Range curve')

plt.figure()
levels_range = np.linspace(100,1000,10)
c4 = plt.contour(batt_cap, CL_CD_ratio, range_electric, levels=levels_range, linewidths=1.2,colors='k',label='Cl/Cd')
c4.levels = [nf(val) for val in c4.levels]
plt.clabel(c4,c4.levels,inline=True,fmt=fmt,fontsize=10)
plt.xlabel(r'Battery capacity (kWh)')
plt.ylabel(r'Cl_Cd ratio')
plt.annotate('Range in km ',(1000,100),color='k')
plt.title('Battery Range curve')
#plt.show()


## Calculate range based on performance 

#inputs for sea level drag 
# rho_sl = 0.002377  #sl/ft3

rho_sl = 1.225 #kg/m3
#MTOW
Cdo = 0.028
AR = 4
e=0.9
pi=3.14
K = 1 / (pi * AR * e)
S = 4*1.875   #m2
#alt = 400
#rho_alt = 
V = np.arange(5,150,10)  #mps at sea level
print("Velocity is ",V)

drag = np.empty(shape=(len(V)))
P = np.empty(shape=(len(V)))
energy_flight = np.empty(shape=(len(V)))
flight_time = 160934 / V  #hawaiian trip in meters/mps for flight time in seconds

def calc_drag_force():
    for i in range(len(V)):
        drag[i] = (0.5 * rho_sl * np.square(V[i]) * S)* (Cdo + K * np.square((MTOW*9.81)/ (0.5*rho_sl* np.square(V[i]) * S))) #N
        P[i] = drag[i] * V[i] /1000 #kW at sea level
        energy_flight[i] = P[i] * 1000 * flight_time[i] / 3600 #Wh

    return drag,P,energy_flight

if __name__ == '__main__':
    calc_drag_force()
    print("Power kW is ",P)
    print("Drag force is ", drag)
    print("Energy consumed for flight is", energy_flight)

# V_ALT = V * sqrt(rho_alt/ rho_sl)
# drag_ALT = Cdo * 0.5 * rho_sl * V_ALT^2 * S + (2* K * (MTOW*2.204)^2)/ rho_sl* V_ALT^2 * S
# P_ALT = drag_ALT * V_ALT

# calculate design velocity
design_vel = V[np.argmin(P)]*3.6 #kph
print("Velocity for minimum power consumption is ", design_vel)  
design_energy = energy_flight[np.argmin(P)]/1000  #kWh
print("Energy at design velocity is ", design_energy)
design_flight_time = flight_time[np.argmin(P)]
min_P = np.min(P)
print("Minimum power drawn in kW is: ", min_P)
#calculate V to achieve L/D max another way to compare
Cle = np.sqrt(Cdo/K)
V_E = np.sqrt((2*MTOW*9.81)/(rho_sl*S*Cle))*3.6
print("Speed for min range and L/Dmax in kph is: ", V_E)

plt.figure(3)
plt.plot(V* 3.6, drag, label="Drag newtons")
plt.legend(loc="upper right")
plt.xlabel(r'Velocity kph')
plt.ylabel(r'Drag/Energy consumed')
plt.title('drag for straight and level flight - 100 mi')

plt.figure(4)
plt.plot(V* 3.6, P, label="Power kW")
plt.legend(loc="upper right")
plt.xlabel(r'Velocity kph')
plt.ylabel(r'Power kW')
plt.title('Power for straight and level flight - 100 mi')

plt.figure(5)
plt.plot(V* 3.6, energy_flight, label="Energy Wh")
plt.legend(loc="upper right")
plt.xlabel(r'Velocity kph')
plt.ylabel(r'energy Wh')
plt.title('energy consumed for straight and level flight in Hawaii - 100 mi')
plt.plot(design_vel, design_energy*1000, marker="o", markersize=10) 
plt.plot(V_E, design_energy*1000, marker="x", markersize=10) 
plt.text(0.65, 500000, "Vel in kph using min power method: %s" % design_vel)
plt.text(0.65, 450000, "Energy in kWh: %s" % design_energy)
plt.text(0.65, 400000, "Flight time in s: %s" % design_flight_time)
plt.text(0.65, 350000, "Vel in kph using L/Dmax method: %s" % V_E)
plt.show()


# check against eviation alice plane

#9 minutes
#150 ktas 

""""
V = np.arange(20,150,20)
Cdo = 0.029
e = 0.83
S = 28.9
MTOW = 7491
W = MTOW*9.81 #newtons
Rho_10000ft = 0.904 
rho_sl = Rho_10000ft
AR = 12.7
calc_drag_force()

design_vel = V[np.argmin(P)]*3.6 #kph
print("Velocity for minimum power consumption in kph is ", design_vel)  
design_energy = energy_flight[np.argmin(P)]/1000  #kWh
print("Energy at design velocity is ", design_energy)
design_flight_time = flight_time[np.argmin(P)]

plt.figure()
plt.plot(V* 3.6, drag, label="Drag newtons")
plt.legend(loc="upper right")
plt.xlabel(r'Velocity kph')
plt.ylabel(r'Drag/Energy consumed')
plt.title('drag for straight and level flight in Hawaii - 100 mi')

plt.figure()
plt.plot(V* 3.6, P, label="Power kW")
plt.legend(loc="upper right")
plt.xlabel(r'Velocity kph')
plt.ylabel(r'Power kW')
plt.title('Power for straight and level flight in Hawaii - 100 mi')

plt.figure()
plt.plot(V* 3.6, energy_flight, label="Energy Wh")
plt.legend(loc="upper right")
plt.xlabel(r'Velocity kph')
plt.ylabel(r'energy Wh')
plt.title('energy consumed for straight and level flight in Hawaii - 100 mi')
plt.plot(design_vel, design_energy*1000, marker="o", markersize=10) 
plt.text(0.65, 50000, "Vel in kph: %s" % design_vel)
plt.text(0.65, 47500, "Energy in kWh: %s" % design_energy)
plt.text(0.65, 45000, "Flight time in s: %s" % design_flight_time)

plt.show()

K = 1 / (pi * AR * e)

Cle = np.sqrt(Cdo/ K)
Ve = np.sqrt((2*MTOW*9.81)/(rho_sl*S*Cle))*3.6
print("Best range speed in kph is: ",Ve)

"""
