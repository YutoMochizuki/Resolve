#!/usr/local/bin/python
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'serif'
import numpy as np
import scipy as sp
import scipy.optimize as so
import scipy.special
import sys


data_recom=np.loadtxt('./CalLines/recommended_TC4use.txt', dtype='unicode')

def lineDB(line_num):
  if line_num < len(data_recom):
    lineDBname=data_recom[line_num].split('_')[0] + '_' + data_recom[line_num].split('_')[1][:2]
    filename="./CalLines/" + data_recom[line_num] + ".dat"
    data=np.loadtxt(filename, comments='#')
    energyDB=data[:,0]
    lgammaDB=data[:,1]
    ampDB=data[:,2]
    return lineDBname, energyDB, lgammaDB, ampDB


file="linefitting_result.txt"
data=np.loadtxt(file)
pixel=data[:,0]
species=data[:,1]
fwhm=data[:,2]
fwhme=data[:,3]
gain=data[:,4]
gaine=data[:,5]


def data_field(pixelnum):
  species_field=species[pixel==pixelnum]
  fwhm_field=fwhm[pixel==pixelnum]
  fwhme_field=fwhme[pixel==pixelnum]
  centerE=[]
  for species_use in species_field:
    linename, energy, lgamma, amp = lineDB(int(species_use))
    centerE.append(energy[0])
  return centerE, fwhm_field, fwhme_field


plt.ylim(0,15)
plt.xlabel('EPI (eV)')
plt.ylabel('FWHM (eV)')
for pixelnum in range(9):
  E,fw,fwe=data_field(pixelnum)
  plt.errorbar(E,fw,fwe,capsize=3, fmt='o', markersize=3, label='pix%d' % pixelnum)
plt.legend()
plt.savefig('../figures/FWHMtrend_A0.png')
plt.clf()
plt.close()


plt.ylim(0,15)
plt.xlabel('EPI (eV)')
plt.ylabel('FWHM (eV)')
for pixelnum in range(9,18):
  if pixelnum != 12 :
    E,fw,fwe=data_field(pixelnum)
    plt.errorbar(E,fw,fwe,capsize=3, fmt='o', markersize=3, label='pix%d' % pixelnum)
plt.legend()
plt.savefig('../figures/FWHMtrend_A1.png')
plt.clf()
plt.close()


plt.ylim(0,15)
plt.xlabel('EPI (eV)')
plt.ylabel('FWHM (eV)')
for pixelnum in range(18,27):
    E,fw,fwe=data_field(pixelnum)
    plt.errorbar(E,fw,fwe,capsize=3, fmt='o', markersize=3, label='pix%d' % pixelnum)
plt.legend()
plt.savefig('../figures/FWHMtrend_B0.png')
plt.clf()
plt.close()


plt.ylim(0,15)
plt.xlabel('EPI (eV)')
plt.ylabel('FWHM (eV)')
for pixelnum in range(27,36):
    E,fw,fwe=data_field(pixelnum)
    plt.errorbar(E,fw,fwe,capsize=3, fmt='o', markersize=3, label='pix%d' % pixelnum)
plt.legend()
plt.savefig('../figures/FWHMtrend_B1.png')
plt.clf()
plt.close()



file="linefitting_result_all.txt"
data=np.loadtxt(file)
species=data[:,0]
fwhm=data[:,1]
fwhme=data[:,2]
gain=data[:,3]
gaine=data[:,4]

def data_field_all():
  species_field=species
  fwhm_field=fwhm
  fwhme_field=fwhme
  centerE=[]
  for species_use in species_field:
    linename, energy, lgamma, amp = lineDB(int(species_use))
    centerE.append(energy[0])
  return centerE, fwhm_field, fwhme_field


plt.ylim(0,15)
plt.xlabel('EPI (eV)')
plt.ylabel('FWHM (eV)')
E,fw,fwe=data_field_all()
plt.errorbar(E,fw,fwe,capsize=3, fmt='o', markersize=3, label='all')
plt.legend()
plt.savefig('../figures/FWHMtrend_all.png')
plt.clf()
plt.close()




