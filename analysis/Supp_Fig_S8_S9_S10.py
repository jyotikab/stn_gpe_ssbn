

#----------------------------------------------------------------------------------------------------------------------
# Reproduces Supp Fig S8, S9, S10
# call the script with: 
# "non-osc" 0.4 0.3
# "osc" 0.4 0.3
# "border" 0.4 0.3
#--------------------------------------------------------------------------------------------------------------------


import numpy as np
import itertools
import sys
home_directory = ""
sys.path.append(home_directory+'/common/')
import analyze_data_ssbn as adata
import analysis_funcs_wo_framework as anal_wo
reload(adata)
import sys
from numpy import *
reload(adata)
import pdb
import pylab as pl
import matplotlib.cm as cm
import plotly.plotly as py
import os
import cPickle as pickle
import matplotlib.pyplot as plt
import pdb
import seaborn as sns

import matplotlib.gridspec as gridspec	
num_threads = 8

start_stim = 500
stop_stim = 5500.

def get_sim_3d(match_pars):
	print path1
	res = adata.analyze_data(prefix1,data_path=path1)
	return res


def concatFiles(path,prefix):
	fhand = file(path+prefix+'.gdf', 'a+')
	for x in xrange(num_threads):
	   filename = path+"/"+prefix+"-{}.gdf".format(x)
	   if os.stat(filename).st_size == 0:
	       continue
	   temp = np.loadtxt(filename)
	   np.savetxt(fhand, temp)
	fhand.close()


#"Enter arguments: 1) regime - "osc", "non-osc" or "border", 2) Gpe bursting ratio (0-1.0), STN bursting ratio (0-1.0)

#print "Enter arguments: 1) regime - osc, non-osc or border, 2) Gpe bursting ratio (0-1.0), STN bursting ratio (0-1.0)"

regime = sys.argv[1]
gpe_bs = sys.argv[2]
stn_bs = sys.argv[3]

if regime == "border": 
	gpe_ip = np.array([900.])[0]
	stn_ip = np.array([1400.])[0]
elif regime == "osc":
	gpe_ip = np.array([500.])[0]
	stn_ip = np.array([1000.])[0]
elif regime == "non-osc":
	gpe_ip = np.array([1300.])[0]
	stn_ip = np.array([1400.])[0]


gpe_rat = gpe_bs 
stn_rat = stn_bs
pre = 200
post = 400
gpe_chg_time = 1500
stn_chg_time = 3500
simtime = 7500

ticks_step = 200
binsize=10.
path1 = '../ssbn_effect/'


fig = pl.figure(figsize=(12,12))
t1 = fig.add_subplot(121)
t2 = fig.add_subplot(122)

prefix = "psdbEffect_midway_change_inh_ip_"+str(gpe_ip)+"_exc_ip_"+str(stn_ip)+"_inh_ratio_"+str(gpe_rat)+"_exc_ratio_"+str(stn_rat)+"_spikes"
data = np.load("output/"+prefix+".npy")
fname = "output/inh_ip_"+str(gpe_ip)+"_exc_ip_"+str(stn_ip)+"_inh_ratio_"+str(gpe_rat)+"_exc_ratio_"+str(stn_rat)+"_spikes.gdf"
np.savetxt(fname,data,fmt = '%.6f')

stn_ids = np.arange(3,1000+3)
gpe_ids = np.arange(1000+3,3000+3)

stn_ids_b = np.arange(3,int(1000*float(stn_bs)+3))
stn_ids_nb = np.arange(int(1000*float(stn_bs))+3,1000+3)
 
gpe_ids_b = np.arange(1000+3,1000+int(2000*float(gpe_bs))+3)
gpe_ids_nb = np.arange(1003+int(2000*float(gpe_bs)),3003)

ids_gpe_b = [ i  for i,x in enumerate(data[:,0]) if x in gpe_ids_b]
ids_gpe_nb = [ i  for i,x in enumerate(data[:,0]) if x in gpe_ids_nb]

ids_stn_b = [ i  for i,x in enumerate(data[:,0]) if x in stn_ids_b]
ids_stn_nb = [ i  for i,x in enumerate(data[:,0]) if x in stn_ids_nb]
ids_stn_b1_phase = np.random.random_integers(np.min(stn_ids_b), np.max(stn_ids_b),10)
 


def psth_rescale(data_stn_b,data_stn_nb,data_gpe_b,data_gpe_nb,bins,y_min,y_max):
	a_gpe_b,b_gpe_b = np.histogram(data_gpe_b,bins=bins)
	a_gpe_nb,b_gpe_nb = np.histogram(data_gpe_nb,bins=bins)	
        a_stn_b,b_stn_b = np.histogram(data_stn_b,bins=bins)
        a_stn_nb,b_stn_nb = np.histogram(data_stn_nb,bins=bins)

	max_a = np.max(np.hstack((a_gpe_b, a_gpe_nb, a_stn_b, a_stn_nb)))
	print max_a
	#y_min - stn_b, stn_nb, gpe_b, gpe_nb
	a_stn_b_scaled = a_stn_b +y_min[0]
	a_stn_nb_scaled = a_stn_nb +y_min[1]
	a_gpe_b_scaled = a_gpe_b + y_min[2]
	a_gpe_nb_scaled = a_gpe_nb + y_min[3]

	return a_stn_b_scaled, a_stn_nb_scaled, a_gpe_b_scaled, a_gpe_nb_scaled, b_gpe_b, b_stn_b
	

# Gpe change
step = 1
ms = 2.5 
#ms = 8.5 
t1.plot(data[ids_stn_b[::step],1],data[ids_stn_b[::step],0],'.',color='steelblue',alpha=0.8,markersize=ms)
t1.plot(data[ids_stn_nb[::step],1],data[ids_stn_nb[::step],0],'.',color='darkslategrey',alpha=0.8,markersize=ms)
t1.plot(data[ids_gpe_b[::step],1],data[ids_gpe_b[::step],0],'.',color='darkorange',alpha=0.8,markersize=ms)
t1.plot(data[ids_gpe_nb[::step],1],data[ids_gpe_nb[::step],0],'.',color='sienna',alpha=0.8,markersize=ms)

bins = np.arange(gpe_chg_time-pre,gpe_chg_time+post+180,binsize)

a_stn_b_scaled, a_stn_nb_scaled, a_gpe_b_scaled, a_gpe_nb_scaled, b_gpe_b_s, b_stn_b_s = psth_rescale(data[ids_stn_b,1], data[ids_stn_nb,1], data[ids_gpe_b,1], data[ids_gpe_nb,1], bins, [np.min(data[ids_stn_b,0]), np.min(data[ids_stn_nb,0]), np.min(data[ids_gpe_b,0]), np.min(data[ids_gpe_nb,0])], [np.max(data[ids_stn_b,0]), np.max(data[ids_stn_nb,0]), np.max(data[ids_gpe_b,0]), np.max(data[ids_gpe_nb,0])])




t1.plot(b_stn_b_s[:-1],a_stn_b_scaled,'k-',linewidth=2.0)
t1.plot(b_stn_b_s[:-1],a_stn_nb_scaled,'k-',linewidth=2.0)

t1.plot(b_gpe_b_s[:-1],a_gpe_b_scaled,'k-',linewidth=2.0)
t1.plot(b_gpe_b_s[:-1],a_gpe_nb_scaled,'k-',linewidth=2.0)
t1.set_xticks(np.arange(gpe_chg_time-pre,gpe_chg_time+post+ticks_step,ticks_step))
t1.set_xticklabels([ str(x) for x in np.arange(gpe_chg_time-pre,gpe_chg_time+post+ticks_step,ticks_step)])
for i,x in enumerate(t1.get_xticklabels()):
		x.set_fontsize(14)
		x.set_fontweight('bold')
t1.set_title(str(int(float(gpe_rat)*100))+"%"+" GPe bursty at 1500ms",fontsize=14,fontweight='bold')
t1.vlines(1500.0,3,3003,'k',linewidth=3.0,linestyles='dashed')
t1.set_xlim(gpe_chg_time-pre,gpe_chg_time+post)
t1.set_ylim(0,3003)

# STN change
 
t2.plot(data[ids_stn_b[::step],1],data[ids_stn_b[::step],0],'.',color='steelblue',alpha=0.8,markersize=ms)
t2.plot(data[ids_stn_nb[::step],1],data[ids_stn_nb[::step],0],'.',color='darkslategrey',alpha=0.8,markersize=ms)
t2.plot(data[ids_gpe_b[::step],1],data[ids_gpe_b[::step],0],'.',color='darkorange',alpha=0.8,markersize=ms)
t2.plot(data[ids_gpe_nb[::step],1],data[ids_gpe_nb[::step],0],'.',color='sienna',alpha=0.8,markersize=ms)

bins_isi = np.linspace(1,50,30)
isi_gpe_b=[]
isi_gpe_nb=[]
isi_stn_b=[]
isi_stn_nb=[]
for x in ids_gpe_b[::5]:
	t_b = data[np.where(data[:,0]==x),1]
	t_b = t_b[t_b>=stn_chg_time]
	if len(t_b) > 0:
		isi_gpe_b.append(np.diff(t_b))


for x in ids_gpe_nb[::5]:
	t_b = data[np.where(data[:,0]==x),1]
	t_b = t_b[t_b>=stn_chg_time]
	if len(t_b) > 0:
		isi_gpe_nb.append(np.diff(t_b))

for x in ids_stn_nb[::5]:
	t_b = data[np.where(data[:,0]==x),1]
	t_b = t_b[t_b>=stn_chg_time]
	if len(t_b) > 0:
		isi_stn_nb.append(np.diff(t_b))

for x in ids_stn_b[::5]:
	t_b = data[np.where(data[:,0]==x),1]
	t_b = t_b[t_b>=stn_chg_time]
	if len(t_b)> 0:
		isi_stn_b.append(np.diff(t_b))

bins = np.arange(stn_chg_time-pre,simtime,binsize)
a_stn_b_scaled, a_stn_nb_scaled, a_gpe_b_scaled, a_gpe_nb_scaled, b_gpe_b_s, b_stn_b_s = psth_rescale(data[ids_stn_b,1], data[ids_stn_nb,1], data[ids_gpe_b,1], data[ids_gpe_nb,1], bins, [np.min(data[ids_stn_b,0]), np.min(data[ids_stn_nb,0]), np.min(data[ids_gpe_b,0]), np.min(data[ids_gpe_nb,0])], [np.max(data[ids_stn_b,0]), np.max(data[ids_stn_nb,0]), np.max(data[ids_gpe_b,0]), np.max(data[ids_gpe_nb,0])])


t2.plot(b_stn_b_s[:-1],a_stn_b_scaled,'k-',linewidth=2.0)
t2.plot(b_stn_b_s[:-1],a_stn_nb_scaled,'k-',linewidth=2.0)

t2.plot(b_gpe_b_s[:-1],a_gpe_b_scaled,'k-',linewidth=2.0)
t2.plot(b_gpe_b_s[:-1],a_gpe_nb_scaled,'k-',linewidth=2.0)
t2.set_xticks(np.arange(stn_chg_time-pre,stn_chg_time+post+ticks_step,ticks_step))
t2.set_xticklabels([ str(x) for x in np.arange(stn_chg_time-pre,stn_chg_time+post+ticks_step,ticks_step)])

for i,x in enumerate(t2.get_xticklabels()):
	x.set_visible(True)
	x.set_fontsize(14)
	x.set_fontweight('bold')


t2.set_title(str(int(float(stn_rat)*100))+"%"+" STN bursty at 3500ms",fontsize=14,fontweight='bold')
t2.vlines(3500.0,3,3003,'k',linewidth=3.0,linestyles='dashed')
t2.set_xlim(stn_chg_time-pre,stn_chg_time+post)
t2.set_ylim(0,3003)

fig.savefig("figs/"+prefix+"_"+regime+"_raster"+".png")


pl.show()
