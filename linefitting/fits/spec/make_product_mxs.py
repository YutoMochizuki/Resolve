#make products from event fits file
#for linefitting

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

r_dir = "/Users/mochizuki/Desktop/Resolve_linefitting/fits/spec"
r_file = "event_out_mxs.fits"

image_file =os.path.join(r_dir, r_file)

hdu_list = fits.open(image_file)
#hdu_list.info()

data = Table.read(os.path.join(r_dir, r_file), hdu=1)

ascii.write([data["EPI2"],data["TIME"]],"c.txt",delimiter='|',overwrite=True)

df_a = pd.read_table('c.txt', header=0,delimiter="|")
#print(df_a)

os.remove('./c.txt')

df_a = df_a.dropna(how='any')
#df_s = df_a.sort_values('EPI2')
#print(df_s)

#a_df = df_s.values
#print(a_df)


###########################
#make lightcurve
###########################


#a = int(len(a_df)) - 1
#tmax = int(a_df[a][1])
#tmin = int(a_df[0][1])

tmax = max(data["TIME"])
tmin = min(data["TIME"])

tbin = 10
nbins = int( (tmax - tmin)/tbin )

lc_src_y, edge = np.histogram(df_a['TIME'],range=(tmin,tmax),bins=nbins)
lc_src_x = (edge[:-1] + edge[1:]) / 2.

fig = go.Figure(data=[
        go.Scatter(
        x=lc_src_x, y=lc_src_y/float(tbin),name="All",
#    error_y = dict(
 #           type = 'data', # value of error bar given in data coordinates
  #          array=np.sqrt(lc_src_y)/float(tbin),
   #         visible=True)
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


###########################
#make spec
###########################



#a = int(len(a_df)) - 1
#smax = int(a_df[a][0])
#smin = int(a_df[0][0])

smax = max(df_a["EPI2"])
smin = min(df_a["EPI2"])

sbin = 0.5
nbins = int( (smax - smin)/sbin )

spec_y, edge = np.histogram(df_a["EPI2"], 
	bins=nbins, range=(smin,smax))
spec_x = (edge[:-1] + edge[1:]) / 2.
spec_x_err = (edge[:-1] - edge[1:]) / 2.
spec_y_err = np.sqrt(spec_y/float(sbin))

fig = go.Figure(data=[
        go.Scatter(
        x=spec_x, y=spec_y/float(sbin),name="All",
#    error_y = dict(
 #           type = 'data', # value of error bar given in data coordinates
  #          array=spec_y_err,
   #         visible=True)
            )
]
               )

fig.update_layout(xaxis_type="log", yaxis_type="log")
fig.update_xaxes(title="Energy(eV)") # X軸タイトルを指定
fig.update_yaxes(title="Counts / eV") # Y軸タイトルを指定
#fig.update_yaxes(scaleanchor="x", scaleratio=1) # Y軸のスケールをX軸と同じに（plt.axis("equal")）
#fig.update_layout(title="spec (2021/10/22 00:00:00 - 2021/10/22 04:00:00)") # グラフタイトルを設定
#fig.update_layout(font={"family":"Meiryo", "size":20}) # フォントファミリとフォントサイズを指定
#fig.update_layout(showlegend=True) # 凡例を強制的に表示（デフォルトでは複数系列あると表示）
fig.update_layout(width=800, height=600) # 図の高さを幅を指定
fig.update_layout(template="plotly_white") # 白背景のテーマに変更

fig.show()

#fig.write_image("./spec.png")


#############################################
#make spec_all.txt
#############################################


df1 = pd.DataFrame(
    data={'x': spec_x, 
          'y': spec_y,
          }
)
#print(df1)
#20keV以上は切り捨て
df1 = df1[df1["x"] < 20000]

#df1 = df1[df1["x"] > 5000]

df1.to_csv('/Users/mochizuki/Desktop/Resolve_linefitting/fits/spec/spec_all_mxs.txt', sep=' ', index=None, header=None)


















