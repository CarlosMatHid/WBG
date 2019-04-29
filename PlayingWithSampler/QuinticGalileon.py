
"""
    QUINTIC GALILEON
    
Created on Mon Feb 25 11:51:10 2019
@author: mateos
"""

# import libraries
import matplotlib.pyplot as plt
import numpy as np

# load data
data=np.loadtxt('G5_Stability_Space.dat')

# read 
x, y = data[:,1], data[:,2]
s = data[:,3]
x_stable = []
y_stable = []
x_unstable = []
y_unstable = []

for i in range(len(s)):
    if s[i]==0:
        x_unstable.append(x[i])
        y_unstable.append(y[i]) 
    else:
        x_stable.append(x[i])
        y_stable.append(y[i]) 

#plot
plt.figure(figsize=(10,10))

plt.scatter(x_stable, y_stable, color='blue', label='stable')
plt.scatter(x_unstable, y_unstable, color='white')
#add details
plt.legend()

plt.xlabel('xi')
plt.ylabel('c3')
plt.title('G5 Stability')
plt.legend()
plt.savefig('G5_Stability')


dataTT = np.loadtxt('G5_scalCls.dat') 
dataMPS = np.loadtxt('G5_matterpower.dat') 
# read data
xTT, yTT = dataTT[:,0], dataTT[:,1]
xMPS, yMPS = dataMPS[:,0], dataMPS[:,1]

# Plot data
plt.figure(figsize=(15,8))
plt.plot(xTT,yTT, label='TT', color="red")
plt.xscale('log')
plt.xlabel('l')
plt.ylabel('C_l')
plt.title('CMB')
plt.legend()
plt.savefig('G5_CMB')


plt.figure(figsize=(15,8))
plt.xscale('log')
plt.xlabel('k')
plt.ylabel('P')
plt.title('Matter PS')
plt.legend()
plt.savefig('G5_MatterPS')
plt.plot(xMPS,yMPS, label='MPS', color="blue")