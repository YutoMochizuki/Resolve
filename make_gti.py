#!/usr/local/anaconda2/bin/python2.7

#gtiを作るプログラム


############################################################################
# [Function]
# Generate GTI FITS file from a text file.
#
# [Usage]
# make_sxs_gti.py (infile) (outfile)
#
############################################################################

########################################################
# Imports
########################################################
from astropy.io import fits
import pandas as pd
import numpy as np
import sys

########################################################
# User-defined parameters
########################################################

########################################################
# Functions
########################################################

########################################################
# Main routine
########################################################
if __name__ == '__main__':

    if len(sys.argv)!=3:
        sys.exit(1)

    # Read data.
    infile = sys.argv[1]
    outfile = sys.argv[2]
    data = pd.read_csv(infile, delim_whitespace=True, header=None, dtype=object)
    data.columns=['START_TIME','END_TIME']

    start=np.array([])
    stop =np.array([])
    
    for index, row in data.iterrows():
        start = np.append(start, row.START_TIME)
        stop  = np.append(stop,  row.END_TIME)
    
    col1 = fits.Column(name='START', format='E', array=start)
    col2 = fits.Column(name='STOP', format='E', array=stop)
    cols = fits.ColDefs([col1, col2])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    prihdr = tbhdu.header
    prihdr.set('EXTNAME', 'STDGTI')
    prihdr.set('HDUCLASS', 'OGIP')
    prihdr.set('HDUCLAS1', 'GTI')
    prihdr.set('HDUCLAS2', 'STANDARD')    
    prihdr.set('MJDREF', '59507.16667')    
    prihdr.set('TIMEZERO', '0')    
    
    #os.remove('./stt.gti')
    tbhdu.writeto(outfile)
