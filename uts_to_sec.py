#!/bin/sh

#utcから秒に直すプログラム

import os
import datetime

from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objs as go

#secに直す関数
def utc_to_time(y,m,d,h):
    'Get UTC from TIME corrected for the leap seconds.'

    utc0 = datetime.datetime(2014,1,1)
    utc1 = datetime.datetime(y,m,d,h)
    
    t0 = Time(utc0, scale='utc')
    t1 = Time(utc1, scale='utc')

    leap_sec0=(t0.tai.value-t0.utc.value).total_seconds()
    leap_sec1=(t1.tai.value-t1.utc.value).total_seconds()
    
    leap_sec=leap_sec1-leap_sec0
    
    td = utc1 - utc0 - datetime.timedelta(seconds=-round(leap_sec))
    
    sec = td.total_seconds()

    return sec

a = utc_to_time(2022,6,1,0)
print("sec =",a)
