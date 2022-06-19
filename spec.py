#!/bin/sh

#eventからスペクトルを作成するプログラム

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
#%matplotlib inline
from astropy.io import fits
from astropy.table import Table
from astropy.io import ascii
import plotly.graph_objs as go
import plotly.io as pio


#cleaneventを入力
r_dir = "./../evt/"
r_file = "xa097081200rsl_p0px1010_cl.evt"

#image_fileの中にeventfiileを代入
image_file =os.path.join(r_dir, r_file)

#fitsのヘッダーを代入
hdu_list = fits.open(image_file)
header = hdu_list['EVENTS'].header

#tableに代入
data = Table.read(os.path.join(r_dir, r_file),hdu=1)

#tableをasciiに変換してtxtをpandasに読み込む
ascii.write([data["PIXEL"],data["ITYPE"],data["EPI2"],data["TIME"]],"a.txt",delimiter='|',overwrite=True)
df_a = pd.read_table('a.txt',header=0,delimiter="|")

os.remove('./a.txt')


#spectra
smax = max(data["EPI2"])
smin = min(data["EPI2"])
sbin = 100
nbins = int( (smax - smin)/sbin )

y_Hp = [0]*37
x_Hp = [0]*37
for i in range(0,36):
    a = df_a[(df_a["PIXEL"] == i) & (df_a["ITYPE"] == 0)]
    y_Hp[i], edge = np.histogram(a["EPI2"], bins=nbins, range=(smin,smax))
    x_Hp[i] = (edge[:-1] + edge[1:]) / 2.

y_all, edge = np.histogram(df_a["EPI2"], bins=nbins, range=(smin,smax))
x_all = (edge[:-1] + edge[1:]) / 2.

fig = go.Figure(data=[
        go.Scatter(
        x=x_all, y=y_all/float(sbin),name="All",
        ),
    
    go.Scatter(
        x=x_Hp[0], y=y_Hp[0]/float(sbin),name="H(px00)",
        ),

    go.Scatter(
        x=x_Hp[1], y=y_Hp[1]/float(sbin),name="H(px01)",
        ),
    
    go.Scatter(
        x=x_Hp[2], y=y_Hp[2]/float(sbin),name="H(px02)",
        ),
    
    go.Scatter(
        x=x_Hp[3], y=y_Hp[3]/float(sbin),name="H(px03)",
        ),
    
    go.Scatter(
        x=x_Hp[4], y=y_Hp[4]/float(sbin),name="H(px04)",
        ),
    
    go.Scatter(
        x=x_Hp[5], y=y_Hp[5]/float(sbin),name="H(px05)",
        ),
    
    go.Scatter(
        x=x_Hp[6], y=y_Hp[6]/float(sbin),name="H(px06)",
        ),
    
    go.Scatter(
        x=x_Hp[7], y=y_Hp[7]/float(sbin),name="H(px07)",
        ),
    
    go.Scatter(
        x=x_Hp[8], y=y_Hp[8]/float(sbin),name="H(px08)",
        ),
    
    go.Scatter(
        x=x_Hp[9], y=y_Hp[9]/float(sbin),name="H(px09)",
        ),
    
    go.Scatter(
        x=x_Hp[10], y=y_Hp[10]/float(sbin),name="H(px10)",
        ),
    
    go.Scatter(
        x=x_Hp[11], y=y_Hp[11]/float(sbin),name="H(px11)",
        ),
    
    go.Scatter(
        x=x_Hp[12], y=y_Hp[12]/float(sbin),name="H(px12)",
        ),
    
    go.Scatter(
        x=x_Hp[13], y=y_Hp[13]/float(sbin),name="H(px13)",
        ),
    
    go.Scatter(
        x=x_Hp[14], y=y_Hp[14]/float(sbin),name="H(px14)",
        ),
    
    go.Scatter(
        x=x_Hp[15], y=y_Hp[15]/float(sbin),name="H(px15)",
        ),
    
    go.Scatter(
        x=x_Hp[16], y=y_Hp[16]/float(sbin),name="H(px16)",
        ),
    
    go.Scatter(
        x=x_Hp[17], y=y_Hp[17]/float(sbin),name="H(px17)",
        ),

    go.Scatter(
        x=x_Hp[18], y=y_Hp[18]/float(sbin),name="H(px18)",
        ),
    
    go.Scatter(
        x=x_Hp[19], y=y_Hp[19]/float(sbin),name="H(px19)",
        ),
    
    go.Scatter(
        x=x_Hp[20], y=y_Hp[20]/float(sbin),name="H(px20)",
        ),
    
    go.Scatter(
        x=x_Hp[21], y=y_Hp[21]/float(sbin),name="H(px21)",
        ),
    
    go.Scatter(
        x=x_Hp[22], y=y_Hp[22]/float(sbin),name="H(px22)",
        ),
    
    go.Scatter(
        x=x_Hp[23], y=y_Hp[23]/float(sbin),name="H(px23)",
        ),
    
    go.Scatter(
        x=x_Hp[24], y=y_Hp[24]/float(sbin),name="H(px24)",
        ),
    
    go.Scatter(
        x=x_Hp[25], y=y_Hp[25]/float(sbin),name="H(px25)",
        ),
    
    go.Scatter(
        x=x_Hp[25], y=y_Hp[25]/float(sbin),name="H(px25)",
        ),
    
    go.Scatter(
        x=x_Hp[26], y=y_Hp[26]/float(sbin),name="H(px26)",
        ),
    
    go.Scatter(
        x=x_Hp[27], y=y_Hp[27]/float(sbin),name="H(px27)",
        ),
    
    go.Scatter(
        x=x_Hp[28], y=y_Hp[28]/float(sbin),name="H(px28)",
        ),
    
    go.Scatter(
        x=x_Hp[29], y=y_Hp[29]/float(sbin),name="H(px29)",
        ),
    
    go.Scatter(
        x=x_Hp[30], y=y_Hp[30]/float(sbin),name="H(px30)",
        ),
    
    go.Scatter(
        x=x_Hp[31], y=y_Hp[31]/float(sbin),name="H(px31)",
        ),
    
    go.Scatter(
        x=x_Hp[32], y=y_Hp[32]/float(sbin),name="H(px32)",
        ),
    
    go.Scatter(
        x=x_Hp[33], y=y_Hp[33]/float(sbin),name="H(px33)",
        ),
    
    go.Scatter(
        x=x_Hp[34], y=y_Hp[34]/float(sbin),name="H(px34)",
        ),
    
    go.Scatter(
        x=x_Hp[35], y=y_Hp[35]/float(sbin),name="H(px35)",
        )

])
fig.update_layout(xaxis_type="log", yaxis_type="log")
fig.update_xaxes(title="Energy(eV)") # X軸タイトルを指定
fig.update_yaxes(title="Counts / eV") # Y軸タイトルを指定
#fig.update_yaxes(scaleanchor="x", scaleratio=1) # Y軸のスケールをX軸と同じに（plt.axis("equal")）
fig.update_layout(title="spec (2021/10/22 00:00:00 - 2021/10/22 04:00:00)") # グラフタイトルを設定
fig.update_layout(font={"family":"Meiryo", "size":20}) # フォントファミリとフォントサイズを指定
fig.update_layout(showlegend=True) # 凡例を強制的に表示（デフォルトでは複数系列あると表示）
fig.update_layout(width=800, height=600) # 図の高さを幅を指定
fig.update_layout(template="plotly_white") # 白背景のテーマに変更

fig.write_image('./spec.png')

fig.show()

