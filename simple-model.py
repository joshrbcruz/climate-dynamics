'''
Dynamical system tracking average global temperature based off of http://archive.dimacs.rutgers.edu/MPE/Energy/DIMACS-EBM.pdf
Describes earth's climate in an equilibrium state
Using albedo varying with time, there are 3 equilibria
Two temperatures are stable equilibria: around 233K and 288K
There is an unstable equilibrium at 265K
Without additional human input, this model assumes that the Earth's average global temperature remains stable at either 233K (snowball earth) or 288K

Human input can be expressed by varying the epsilon parameter to lower/raise greenhouse gas effect
'''

import math
import numpy as np
from scipy import integrate 
import matplotlib.pyplot as plt

#Input constants
r = 2.912 #Average heat capacity of the Earth/atmosphere system (W-yr/m^2K)
q = 342 #Annual global mean incoming solar radiation (W/m^2)
alpha = 0.3 #Albedo (Dimensionless)
sigma = 5.67 * math.pow(10, -8) #Stefan-Boltzmann Constant (W/m^2K^4)
epsilon = 0.6137212149814136 #Greenhouse effect (0<epsilon<1). Models effect of greenhouse gases
#This value of epsilon results in equilibrium, where average temp tends to 288 K (15 deg. celsius)
#Calculated by calculating the zero of (1-alpha)Q - epsilon*sigma*288^4

tempInit = 288 

#Range of years
t = np.linspace(0, 10)

def albedo(temp):
	return 0.7-0.4*(math.e**((temp-265)/5)/(1+math.e**((temp-265)/5)))

#Based off section 4 in http://archive.dimacs.rutgers.edu/MPE/Energy/DIMACS-EBM.pdf
def temperatureDynamics(temperature, t):
	dTemperature = (q*(1-albedo(temperature)) - sigma*epsilon*math.pow(temperature, 4))/r
	return dTemperature

temperatures = integrate.odeint(temperatureDynamics, tempInit, t)	

plt.plot(t, temperatures)
plt.xlabel('time (years)')
plt.ylabel('temperature')
plt.show()