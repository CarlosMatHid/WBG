"""
COMPARING WBG with CovGalileons and LCDM

"""

# import libraries
import matplotlib.pyplot as plt
import numpy as np


# load data
#data=np.loadtxt('H_QuarticCG_cache_FRW.dat')
data=np.loadtxt('QuarticCG_cache_FRW.dat')
data3=np.loadtxt('QuarticCGnew_cache_FRW.dat')
data2 = np.loadtxt('QuarticCG_cache_BackgroundEFT.dat')
data4 = np.loadtxt('QuarticCGtestOm_m0p5_cache_FRW.dat')
data5 = np.loadtxt('QuarticCGtestCRAZY_cache_FRW.dat')
# read 
a, Hconf, Hdot= data[:,0], data[:,3]*3*10**5/100, data[:,4]
a_2, Omega, Lambda = data2[:,0], data2[:,3]-1, data2[:,9]
a3, Hconf3, Hdot3= data3[:,0], data3[:,3]*3*10**5/100, data3[:,4]
a4, Hconf4, Hdot4= data4[:,0], data4[:,3]*3*10**5/100, data4[:,4]
a5, Hconf5 = data5[:,0], data5[:,3]*3*10**5/100
#xi = 2.5
#c4 = -1/9*xi**(-2)+2/3*0.7*xi**(-4)
#c5 = 0
#Omega = a**4/2/H^4 *xi**4 *(c4) #-6 *c5 *xi*(1-Hdot/(H**2))

#plt.plot(a,Hconf, color="black")
#0.00223607(1/(a**2) + 31110/a + (  1 + 62220*a + 9.67832*10**8 *a**2 + 2.75556*10**10 *a**8)**(1/2)/(a**2))**(1/2)
#plt.plot(a,Hconf, label='G4 CovGal EFTCAMB',marker='o', linestyle='None', color="green",markersize=2)
plt.plot(a3,Hconf3, label='G4 CovGal new EFTCAMB',marker='o', linestyle='None', color="red",markersize=2)
plt.plot(a5,Hconf5, label='G4 CovGalcrazy new EFTCAMB',marker='o', linestyle='None', color="orange",markersize=2)
##plt.plot(a, 0.00223606797749979*(a**(-2) + 30000/a + (1 + 60000*a + 9*10**8*a**2 + 2.79996*10**10*a**8)**(1/2)*a**(-2))**(1/2), label='G4 CovGal - WBG(background)', color="blue",marker='o', linestyle='None')
#plt.plot(a, 0.00223607*(1/(a**2) + 31110/a + (  1 + 62220*a + 9.67832*10**8*a**2 + 2.75556*10**10 *a**8)**(1/2)/(a**2))**(1/2), label='G4 CovGal - WBG(background)', color="blue",marker='o', linestyle='None',markersize=2)
#plt.plot(a4, 0.00223606797749979*( 1/a4**2 + 26400/a4 + (1+ 52800*a4 + 6.9696*10**8*a4**2 + 2.94396*10**10*a4**8)**(1/2)/a4**2)**(1/2), label='G4 CovGal new - WBG(background)', color="purple",marker='o', linestyle='None',markersize=2)

#######MALplt.plot(a5, 0.00223607*(1/a5**2 + 29300/a5 + (1+ 58600*a5 + 8.5849*10**8*a5**2 + 2.82796*10**10*a5**8)**(1/2)/a5**2)**(1/2), label='G4 CovGal new - WBG(background)', color="purple",marker='o', linestyle='None',markersize=2)
#plt.scatter(a,a**4/2/Hconf**4 *xi**4 *(c4), label='QuarticGalileon', color="red")
#add details

plt.xscale('log')
plt.yscale('log')
plt.xlabel('a')
plt.ylabel('H(a)')
plt.title('H(a)')
plt.legend()
#plt.figure(figsize=(25,10))
plt.figure() #figsize=(25,10)

'plot Omega'
plt.plot(a_2,Omega, label='G4 CovGal EFTCAMB',marker='o', linestyle='None', color="green",markersize=2)
plt.plot(a, -1 - (5.2065* a**8)/(0.000032144 + 1 *a +  0.000032144 *(1 + 62220 *a + 9.67832*10**8* a**2 + 2.75556*10**10 *a**8)**(1/2))**2 , label='G4 CovGal - WBG(background)', color="blue",marker='o', linestyle='None',markersize=2)
plt.xlabel('a')
plt.ylabel('Omega(a)')
plt.title('Omega(a)')
plt.legend()
plt.figure() #figsize=(25,10)

'plot Lambda'
plt.plot(a_2,Lambda/(a_2**2), label='G4 CovGal EFTCAMB',marker='o', linestyle='None', color="green",markersize=2)
#plt.plot(a, -1 - (5.2065* a**8)/(0.000032144 + 1 *a +  0.000032144 *(1 + 62220 *a + 9.67832*10**8* a**2 + 2.75556*10**10 *a**8)**(1/2))**2 , label='G4 CovGal - WBG(background)', color="blue",marker='o', linestyle='None',markersize=2)
plt.xlabel('a')
plt.ylabel('Lambda(a)')
plt.title('Lambda(a)')
plt.legend()