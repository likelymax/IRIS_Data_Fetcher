#%%
import os
import numpy as np
from obspy import read
import sys
import argparse
import pandas as pd
import datetime as dt

#%%
# Create an ArgumentParser object
parser = argparse.ArgumentParser(description='Example command line arguments')

# %%
# Add arguments to the parser
parser.add_argument('-R', type=str, help='RF file', required=True)
#%%
# Parse the arguments from the command line
args = parser.parse_args()

# %% read the file
filename = args.R
f = read(filename)
net = f[0].stats.network
stn = f[0].stats.station
stlo = f[0].stats.sac.stlo
stla = f[0].stats.sac.stla
year = f[0].stats.sac.nzyear
day = np.int64(f[0].stats.sac.nzjday)
day = int(day)
hour = f[0].stats.sac.nzhour
minut = f[0].stats.sac.nzmin
sec = int(f[0].stats.sac.nzmsec/1000)
date = dt.datetime(year, 1, 1) + dt.timedelta(day - 1)
year = date.year
mon = date.month
day = date.day
tmptime = dt.datetime(year, mon, day, hour, minut, int(sec/1000))
endtime = tmptime + dt.timedelta(hours = 0.5)
startime = endtime - dt.timedelta(hours = 1)
format_startime = startime.strftime("%Y-%m-%dT%H:%M:%S")
format_endtime = endtime.strftime("%Y-%m-%dT%H:%M:%S")
with open("out.txt", "w") as f:
    f.write("%s %s %f %f %s %s\n" % (net, stn, stlo, stla, format_startime, format_endtime))

