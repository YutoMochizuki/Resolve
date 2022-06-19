#!/bin/sh

#ピクセルごとのイメージを作るプログラム

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
ascii.write([data["PIXEL"], data["ITYPE"], data["EPI2"], data["TIME"]], "a.txt" ,delimiter='|', overwrite=True)
df = pd.read_table('a.txt', header=0, delimiter="|")

os.remove('./a.txt')

#countの整理
#count開始時間の最大最小
tmax_count = max(data["TIME"])
tmin_count = min(data["TIME"])

#tmax_obs = header["TSTOP"]
#tmin_obs = header["TSTART"]
#一時的にこうしてるだけで本当はTSTARTを使う
tmax_obs = utc_to_time(2021,10,20,4)
tmin_obs = utc_to_time(2021,10,20,0)

#ビン
tbin = 100
nbins = int( (tmax_obs - tmin_obs)/tbin )

#全ピクセル全グレードのヒストグラムの作成
y_all, edge = np.histogram(df["TIME"], bins=nbins, range=(tmin_count, tmax_count))
x_all = (edge[:-1] + edge[1:]) / 2.

x_len = len(x_all)

#ピクセルごとのcountの作成
y_Hp = np.zeros((36, x_len))
x_Hp = np.zeros((36, x_len))
counts = [0]*36
for i in range(0, 36):
        #a = df[(df["PIXEL"]==i) & (df["ITYPE"]==0)]
        pi = df[(df["PIXEL"]==i)]
        y_Hp[i], edge = np.histogram(pi["TIME"], bins=nbins, range=(tmin_count, tmax_count))
        x_Hp[i] = (edge[:-1] + edge[1:]) / 2.
        for j in range(len(x_Hp[0])):
            counts[i]  += y_Hp[i][j]/(tmax_obs - tmin_obs)

z=np.zeros((7,7))
z[:,:] = np.nan

z[0,0] = counts[12]
z[1,2] = counts[14]
z[1,3] = counts[16]
z[1,4] = counts[8]
z[1,5] = counts[6]
z[1,6] = counts[5]

z[2,1] = counts[11]
z[2,2] = counts[13]
z[2,3] = counts[15]
z[2,4] = counts[7]
z[2,5] = counts[4]
z[2,6] = counts[3]

z[3,1] = counts[9]
z[3,2] = counts[10]
z[3,3] = counts[17]
z[3,4] = counts[0]
z[3,5] = counts[2]
z[3,6] = counts[1]

z[4,1] = counts[19]
z[4,2] = counts[20]
z[4,3] = counts[18]
z[4,4] = counts[35]
z[4,5] = counts[28]
z[4,6] = counts[27]

z[5,1] = counts[21]
z[5,2] = counts[22]
z[5,3] = counts[25]
z[5,4] = counts[33]
z[5,5] = counts[31]
z[5,6] = counts[29]

z[6,1] = counts[23]
z[6,2] = counts[24]
z[6,3] = counts[26]
z[6,4] = counts[34]
z[6,5] = counts[32]
z[6,6] = counts[30]

#heatmapのxとy
x = []
y = []

l=0
a=0

#data = np.zeros((10000,3))

#for i in range(0, 7):
#    for j in range(0, 7):
#        data.append = (i,j,z[i,j])



#for i in range(0, 7):
#    for j in range(0, 7):
#        if i==1 and j==1:
#            l+=1
#            #print("l=",l)
#            continue
#        if math.isnan(z[i,j]): #nanを0に変換
#            z[i,j] = 0
#            #print("i=",i)
#            #print("j=",j)
#        p = int(z[i,j])
#        if p == 0:
#            a+=1
#            #print("a=",a)
#            continue
#        for k in range(p):
#            x.append(i)
#            y.append(j)

#print("z[0,0]=",z[0,0])
#print("z[1,1]=",z[5,5])
#print("x=",x)
#x.append(0)
#y.append(0)

#plotly
plot = []

#plot.append(
#go.Histogram2d(
#x = x,
#y = y,
#colorbar = dict(len=0.8,  # カラーバーの長さを0.8に（デフォルトは1）
#outlinecolor='black',  # カラーバーの枠線の色
#outlinewidth=2,  # カラーバーの枠線の太さ
#bordercolor='gray',  # カラーバーとラベルを含むカラーバー自体の枠線の色
#borderwidth=1,  # カラーバーとラベルを含むカラーバー自体の枠線の太さ
#title=dict(
#text='bar',
#side='right',  # カラーバーのタイトルをつける位置（デフォルトはtop）
#),
#),
#colorscale='Spectral',  # カラーバーの色合いをjetに変更
#zmin=0,  # カラーバーの最小値
#zmax=5,  # カラーバーの最大値
#))



for i in range(0, 7):
    for j in range(0, 7):
        plot.append(
	go.Heatmap(
	x = [i],
	y = [j],
	z = [z[i,j]],
	colorbar = dict(len=0.8,  # カラーバーの長さを0.8に（デフォルトは1）
	outlinecolor='black',  # カラーバーの枠線の色
	#outlinewidth=2,  # カラーバーの枠線の太さ
	#bordercolor='gray',  # カラーバーとラベルを含むカラーバー自体の枠線の色
	borderwidth=1,  # カラーバーとラベルを含むカラーバー自体の枠線の太さ
	title=dict(
	#text='bar',
	#side='right',  # カラーバーのタイトルをつける位置（デフォルトはtop）
	),
	),
	colorscale='Spectral',  # カラーバーの色合いをjetに変更
	zmin=0,  # カラーバーの最小値
	zmax=1,  # カラーバーの最大値
	)
        )


#plot.append(
#go.Histogram2d(
#z = [z[0,0]],
#x = [0],
#y = [0],
#colorbar = dict(len=0.8,  # カラーバーの長さを0.8に（デフォルトは1）
#outlinecolor='black',  # カラーバーの枠線の色
#outlinewidth=2,  # カラーバーの枠線の太さ
#bordercolor='gray',  # カラーバーとラベルを含むカラーバー自体の枠線の色
#borderwidth=1,  # カラーバーとラベルを含むカラーバー自体の枠線の太さ
#title=dict(
#text='bar',
#side='right',  # カラーバーのタイトルをつける位置（デフォルトはtop）
#),
#),
#colorscale='Spectral',  # カラーバーの色合いをjetに変更
#min=0,  # カラーバーの最小値
#zmax=5,  # カラーバーの最大値
#))

#plot.append(
#go.Histogram2d(
#z = [z[1,1]],
#x = [1],
#y = [1],
#colorbar = dict(len=0.8,  # カラーバーの長さを0.8に（デフォルトは1）
#outlinecolor='black',  # カラーバーの枠線の色
#outlinewidth=2,  # カラーバーの枠線の太さ
#bordercolor='gray',  # カラーバーとラベルを含むカラーバー自体の枠線の色
#borderwidth=1,  # カラーバーとラベルを含むカラーバー自体の枠線の太さ
#title=dict(
#text='bar',
#side='right',  # カラーバーのタイトルをつける位置（デフォルトはtop）
#),
#),
#colorscale='Spectral',  # カラーバーの色合いをjetに変更
#zmin=0,  # カラーバーの最小値
#zmax=5,  # カラーバーの最大値
#))


# Define layout
layout = go.Layout(
        title='Rate(1/s) of HP+MP+MS+LP+LS(2021/10/22 00:00:00 - 2021/10/22 04:00:00)',
        xaxis=dict(
            title = 'X',
            type = 'linear',
            autorange = False,
            range = [-0.5, 6.5],
            showgrid=False,
            ),
        yaxis=dict(
            title = 'Y',
            type = 'linear',
            autorange = False,
            range = [-0.5, 6.5],
            showgrid=False,
            ),
        width=800,
        height=800,
        showlegend=False,
        #auto_open=False
)

fig = go.Figure(
    plot, layout)

fig.write_image("./image.png")

fig.show()

