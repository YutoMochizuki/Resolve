#!/usr/local/bin/python
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'serif'
import numpy as np
import scipy as sp
import scipy.optimize as so
import scipy.special
import sys



file="linefitting_result_all.txt"
data=np.loadtxt(file)
species=data[:,0]
fwhm=data[:,1]
fwhme=data[:,2]
gain=data[:,3]
gaine=data[:,4]

#plt.ylim(0,15)
plt.xlabel('EPI (eV)')
plt.ylabel('EPI shift (eV)')
plt.errorbar(species,gain*species,gaine*species,capsize=3, fmt='o', markersize=3, label='all')
plt.legend()
plt.savefig('../figures/gaintrend_all.png')
plt.clf()
plt.close()




