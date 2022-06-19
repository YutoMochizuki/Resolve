#!/bin/sh

#rslmxsgtiを手動で動かすために、mxsのLEDのgtiを作るプログラム


import os
import datetime
import math

from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objs as go

#%matplotlib inline


#hk01を入力
r_dir = "./"
r_file = "xa097081220rsl_a0.hk1.gz"

#image_fileの中にeventfiileを代入
image_file = os.path.join(r_dir, r_file)

#fitsのヘッダーを代入
hdu_list = fits.open(image_file)

#header=3列目を代入
data = Table.read(os.path.join(r_dir, r_file), hdu=3)

#tableをasciiに変換してtxtをpandasに読み込む
#LED1
ascii.write([data["TIME"], data["FWE_TI_LED1_ON"],data["FWE_TI_LED1_OFF"],data["FWE_LED1_PLS_LEN_CAL"], data["FWE_LED1_PLS_SPC_CAL"]], "a.txt", delimiter='|', overwrite=True)
df = pd.read_table('a.txt', header=0,delimiter="|")

#一時ファイルのdelete
os.remove('./a.txt')

#sに直す
df["FWE_TI_LED1_ON"] = df["FWE_TI_LED1_ON"]/64.0
df["FWE_TI_LED1_OFF"] = df["FWE_TI_LED1_OFF"]/64.0

#UTCに直す関数
def time_to_utc(time):
    'Get UTC from TIME corrected for the leap seconds.'

    utc0 = datetime.datetime(1980,1,6)
    utc1 = utc0 + datetime.timedelta(seconds=time)

    t0 = Time(utc0, scale='utc')
    t1 = Time(utc1, scale='utc')

    leap_sec0=(t0.tai.value-t0.utc.value).total_seconds()
    leap_sec1=(t1.tai.value-t1.utc.value).total_seconds()

    leap_sec=leap_sec1-leap_sec0

    utc1 = utc1 + datetime.timedelta(seconds=-round(leap_sec))

    return utc1

#timeに直す関数
def utc_to_time(y,m,d,h,mi,s,ms):
    'Get UTC from TIME corrected for the leap seconds.'

    utc0 = datetime.datetime(2014,1,1)
    utc1 = datetime.datetime(y,m,d,h,mi,s)

    t0 = Time(utc0, scale='utc')
    t1 = Time(utc1, scale='utc')

    leap_sec0=(t0.tai.value-t0.utc.value).total_seconds()
    leap_sec1=(t1.tai.value-t1.utc.value).total_seconds()

    leap_sec=leap_sec1-leap_sec0

    td = utc1 - utc0 - datetime.timedelta(seconds=-round(leap_sec))

    sec = td.total_seconds()

    return sec

b = time_to_utc(df.iloc[len(df) - 1][1])

#bのLEDが始まった時間
a = utc_to_time(b.year,b.month,b.day,b.hour,b.minute,b.second,b.microsecond)
#!!!動作確認必要!!!
#prin(a)


#何回LEDをON/OFFしたのかどうか
times = (df.iloc[len(df) - 1][0] - a)/df.iloc[len(df) - 1][4]/0.001

start = []
stop =[]
start.append(a)
for i in range(1, times):
    stop.append(start[i-1]+df.iloc[len(df) - 1][3]*0.001)
    if i == times:
        break
    else:
        start.append(stop[i-1]+(df.iloc[len(df) - 1][4]-df.iloc[len(df) - 1][3])*0.001)


#startとstopのgti
dim2list = [[start[i], stop[i]] for i in range(len(stop))]

dd = pd.DataFrame(dim2list, columns=['start', 'stop'])

dd.to_csv('mxs_LED1_gti.txt', sep='\t', index=False, header=None)
