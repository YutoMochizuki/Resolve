#!/bin/sh

#clean eventからlightcurveを作るプログラム

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

#%matplotlib inline

#UTCに直す関数
def time_to_utc(time):
    'Get UTC from TIME corrected for the leap seconds.'

    utc0 = datetime.datetime(2014,1,1)
    utc1 = utc0 + datetime.timedelta(seconds=time)

    t0 = Time(utc0, scale='utc')
    t1 = Time(utc1, scale='utc')

    leap_sec0=(t0.tai.value-t0.utc.value).total_seconds()
    leap_sec1=(t1.tai.value-t1.utc.value).total_seconds()

    leap_sec=leap_sec1-leap_sec0

    utc1 = utc1 + datetime.timedelta(seconds=-round(leap_sec))

    return utc1

#timeに直す関数
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

r_dir = "./../evt/"
r_file = "xa097081200rsl_p0px1010_cl.evt"

#image_fileの中にeventfiileを代入
image_file = os.path.join(r_dir, r_file)

#fitsのヘッダーを代入
hdu_list = fits.open(image_file)
header = hdu_list['EVENTS'].header

#tableに代入
data = Table.read(os.path.join(r_dir, r_file), hdu=1)

#tableをasciiに変換してtxtをpandasに読み込む
ascii.write([data["PIXEL"], data["ITYPE"], data["EPI2"], data["TIME"]], "a.txt", delimiter='|', overwrite=True)
df = pd.read_table('a.txt', header=0,delimiter="|")

os.remove('./a.txt')

#data整理
#Hp

#light curve
tmax_count = max(data["TIME"])
tmin_count = min(data["TIME"])

#tmax_obs = header["TSTOP"]
#tmin_obs = header["TSTART"]
#一時的にこうしてるだけで本当はTSTARTを使う
tmax_obs = utc_to_time(2021,10,20,4)
tmin_obs = utc_to_time(2021,10,20,0)


print(time_to_utc(tmax_obs))
print(time_to_utc(tmin_obs))

print(time_to_utc(tmax_count))
print(time_to_utc(tmin_count))

tbin = 100
nbins = int((tmax_obs - tmin_obs) / tbin)

y_all, edge = np.histogram(df["TIME"], bins=nbins, range=(tmin_count, tmax_count))
x_all = (edge[:-1] + edge[1:]) / 2.

x_len = len(x_all)

y_Hp = np.zeros((36, x_len))
x_Hp = np.zeros((36, x_len))
for i in range(0, 36):
        pi = df[(df["PIXEL"]==i) & (df["ITYPE"]==0)]
        y_Hp[i], edge = np.histogram(pi["TIME"], bins=nbins, range=(tmin_count,tmax_count))
        x_Hp[i] = (edge[:-1] + edge[1:]) / 2.
 
time = [time_to_utc(t) for t in x_Hp[0]]

scatters = []
scatters.append(
    go.Scatter(
        x=time, y=y_all/float(tbin), name="All"
    )
)
for i in range(36):
    pixel_str = str(i).rjust(2, "0")
    scatters.append(
        go.Scatter(
            x=time, y=y_Hp[i]/float(tbin), name=f"Hp({pixel_str})"
        )
    )

fig = go.Figure(data=scatters)
#fig.update_layout(xaxis_type="linear", yaxis_type="linear")
fig.update_xaxes(title="Time (s)") # X軸タイトルを指定
fig.update_yaxes(title="Counts /sec") # Y軸タイトルを指定
#fig.update_yaxes(scaleanchor="x", scaleratio=1) # Y軸のスケールをX軸と同じに（plt.axis("equal")）
fig.update_layout(title="lightcurve (2021/10/22 00:00:00 - 2021/10/22 04:00:00)") # グラフタイトルを設定
fig.update_layout(font={"family":"Meiryo", "size":20}) # フォントファミリとフォントサイズを指定
fig.update_layout(showlegend=True) # 凡例を強制的に表示（デフォルトでは複数系列あると表示）
#fig.update_layout(width=800, height=600) # 図の高さを幅を指定
fig.update_layout(template="plotly_white") # 白背景のテーマに変更

fig.write_image("./lc.png")

fig.show()
