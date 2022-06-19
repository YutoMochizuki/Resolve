#!/bin/sh

#gain historyからその妥当性を確かめるプログラム


#  Col  Name             Format[Units](Range)      Comment
#    1 TIME               1D [sec]             Seconds from 01 Jan 2019 00:00:00
#    2 PIXEL              1I                   Pixel Number range 0-35
#    3 COR_FIT            1D                   Ratio of fitted energy to profile
#    4 COR_AVE            1D                   Ratio of bin average energy to profile
#    5 CHISQ              1D                   Reduced Chi^2 of fit
#    6 AVGUNBIN           1D                   Unbinned average energy
#    7 AVGBIN             1D                   Binned average energy
#    8 AVGFIT             1D                   Energy from fit
#    9 SHIFT              1D                   Fitted energy shift
#   10 SCALE              1D                   Fitted count scale factor
#   11 BGRND              1D                   Fitted background counts
#   12 SLOPE              1D                   Fitted background slope
#   13 WIDTH              1D                   Fitted Gaussian convolution width
#   14 TELAPSE            1D                   Seconds between first and last events
#   15 EXPOSURE           1D                   Exposure seconds between first and last events
#   16 NEVENT             1I                   Number of events in fitting group
#   17 BINMESH            309D                 Energy bins
#   18 SPECTRUM           309D                 Counts spectrum
#   19 FITPROF            309D                 Fitted profile
#   20 TEMP_FIT           1D                   Temperature from fit
#   21 TEMP_AVE           1D                   Temperature from bin average

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

r_dir = "/Users/mochizuki/Desktop/"
r_file = "gain_out_run.fits"

image_file =os.path.join(r_dir, r_file)

hdu_list = fits.open(image_file)
#hdu_list.info()

data = Table.read(os.path.join(r_dir, r_file), hdu=1)

fig = go.Figure(data=[
        go.Scatter(
        x=data['TIME'], y=data["COR_FIT"],name="COR_FIT"),
            go.Scatter(
        x=data['TIME'], y=data["BGRND"],name="BGRND"),
    go.Scatter(
        x=data['TIME'], y=data["CHISQ"],name="CHISQ"),
    go.Scatter(
        x=data['TIME'], y=data["WIDTH"],name="WIDTH"),
    go.Scatter(
        x=data['TIME'], y=data["SHIFT"],name="SHIFT"),
    go.Scatter(
        x=data['TIME'], y=data["SCALE"],name="SCALE"),
    go.Scatter(
        x=data['TIME'], y=data["SLOPE"],name="SLOPE"
            )
]
               )

#fig.update_layout(xaxis_type="linear", yaxis_type="linear")
fig.update_xaxes(title="Time (s)") # X軸タイトルを指定
fig.update_yaxes(title="Counts /sec") # Y軸タイトルを指定
#fig.update_yaxes(scaleanchor="x", scaleratio=1) # Y軸のスケールをX軸と同じに（plt.axis("equal")）
#fig.update_layout(title="lightcurve (2021/10/22 00:00:00 - 2021/10/22 04:00:00)") # グラフタイトルを設定
fig.update_layout(font={"family":"Meiryo", "size":20}) # フォントファミリとフォントサイズを指定
fig.update_layout(showlegend=True) # 凡例を強制的に表示（デフォルトでは複数系列あると表示）
#fig.update_layout(width=800, height=600) # 図の高さを幅を指定
fig.update_layout(template="plotly_white") # 白背景のテーマに変更

fig.show()
#fig.write_image("./lc.png")
