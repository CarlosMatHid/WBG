"""
COMPARING WBG with CovGalileons and LCDM

"""

# import libraries
import matplotlib.pyplot as plt
import numpy as np


# load data

data = np.loadtxt('QuarticGalileon_right_values_cache_FRW.dat')
dataO = np.loadtxt('QuarticGalileon_right_values_cache_BackgroundEFT.dat')
# read 
a, Hconf = data[:,0], data[:,3]*2.99792458*10**5/70.9
aO, Omega, Lambda, c = dataO[:,0], dataO[:,3]-1, dataO[:,9]/((70.9/(2.99792458*10**5))**2), dataO[:,7]

'Hconf(a)/H0'
# OmR = 0
#plt.plot(a, 0.00223607*( 1/a**2 + 29300/a + (  1 + 58600* a + 8.5849*10**8*a**2 + 2.82796*10**10*a**8)**(1/2)/a**2)**(1/2), label='G4 CovGal - WBG(background), OmR = 0', color="purple",marker='o', linestyle='None',markersize=2)
# OmR = e-5
#plt.plot(a, 0.00223606797749979*(1/a**2 + 29300/a + (  1 + 58600* a + 8.5849*10**8* a**2 + 2.82796*10**10 *a**8)**(1/2)/a**2)**(1/2), label='G4 CovGal - WBG(background), OmR = e-5', color="blue",marker='o', linestyle='None',markersize=2)
#OmR = 8.32e-5
plt.plot(a, 7.83741*10**(-9)*( 6.77249*10**11/a**2 + 2.38502*10**15/a + (  114123*(3.52169*10**13 + 2.48042*10**17*a + 4.36757*10**20*a**2 + 1.43858*10**22*a**8)**(1/2))/a**2)**(1/2), label='G4 CovGal - WBG(background), OmR = 8.32e-5', color="blue",marker='o', linestyle='None',markersize=2)
# EFTCAMB:
plt.plot(a,Hconf, label='G4 CovGal EFTCAMB',marker='o', linestyle='None', color="orange",markersize=2)

#add details

plt.xscale('log')
plt.yscale('log')
plt.xlabel('a')
plt.ylabel('H(a)')
plt.title('H(a)')
plt.legend()
#plt.figure(figsize=(25,10))
plt.figure() #figsize=(25,10)

'Omega'
#EFTCAMB:
plt.plot(aO ,Omega, label='G4 CovGal EFTCAMB',marker='o', linestyle='None', color="orange",markersize=2)
#CovGalileon:
plt.plot(aO,-1 - (3.17958*10**31*a**8)/(6.77249*10**11 + 2.38502*10**15*a +114123*(3.52169*10**13 + 2.48042*10**17* a + 4.36757*10**20 *a**2 + 1.43858*10**22 *a**8)**(1/2))**2, label='G4 CovGal - WBG(background)', color="blue",marker='o', linestyle='None',markersize=2)
#add details
plt.xlabel('a')
plt.ylabel('Omega(a)')
plt.title('Omega(a)')
plt.legend()
plt.figure() #figsize=(25,10)

'plot Lambda'
#EFTCAMB
plt.plot(aO,Lambda, label='G4 CovGal EFTCAMB',marker='o', linestyle='None', color="orange",markersize=2)
#CovGalileon:

#add details
plt.xlabel('a')
plt.ylabel('Lambda(a)')
plt.title('Lambda(a)')
plt.legend()

'plot Lambda'
plt.figure()
#EFTCAMB
plt.plot(aO, c, label='G4 CovGal EFTCAMB',marker='o', linestyle='None', color="orange",markersize=2)
#CovGalileon:

#add details
plt.xlabel('a')
plt.ylabel('c(a)')
plt.title('c(a)')
plt.legend()