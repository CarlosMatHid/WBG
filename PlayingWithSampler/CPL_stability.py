#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 17:21:20 2019

@author: mateos
"""

# import libraries
import matplotlib.pyplot as plt
import numpy as np

# load data
data=np.loadtxt('SampleCPL_Stability_Space.dat')

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

plt.scatter(x_stable, y_stable, color='blue')
plt.scatter(x_unstable, y_unstable, color='red')
#add details
plt.legend()

plt.xlabel('w0')
plt.ylabel('wa')
plt.title('CPL Stability')

plt.savefig('CPL_Stability')