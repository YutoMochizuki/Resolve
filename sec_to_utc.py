#!/bin/sh

#秒からUTCに直すプログラム

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

#UTCに直す関数
def time_to_utc(time):
    'Get UTC from TIME corrected for the leap seconds.'

    utc0 = datetime.datetime(1980,1,1)
    utc1 = utc0 + datetime.timedelta(seconds=time)

    t0 = Time(utc0, scale='utc')
    t1 = Time(utc1, scale='utc')

    leap_sec0=(t0.tai.value-t0.utc.value).total_seconds()
    leap_sec1=(t1.tai.value-t1.utc.value).total_seconds()

    leap_sec=leap_sec1-leap_sec0

    utc1 = utc1 + datetime.timedelta(seconds=-round(leap_sec))

    return utc1

b = time_to_utc(1318754066.875)
print(b)
