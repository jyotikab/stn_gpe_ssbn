

#----------------------------------------------------------------------------------------------
# Reproduces Figure 6
# Requires the file: "burst_data_all_border.pickle"
# In order to generate burst_data_all_border.pickle, run calc_burst_length_all.py
#
#-----------------------------------------------------------------------------------------------





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
from plotly.graph_objs import *
import os
import cPickle as pickle
import matplotlib.pyplot as plt
import pdb
import seaborn as sns

import matplotlib.gridspec as gridspec	
import scipy.signal as sp_sig
from scipy.interpolate import interp2d, interp1d
from itertools import groupby
import scipy.signal as sp_sig
import scipy as sp



pl.ion()
regime = sys.argv[1]
gpe_bs = sys.argv[2]
stn_bs = sys.argv[3]

if regime == "border":

	gpe_ip = np.array([900.])[0]
	stn_ip = np.array([1400.])[0]
	max_pow = 1 # also serves as threshold 
	min_pow = 0


bin_width=5 #ms
num_trials = 5

fs = 200 # 1000/5ms
nperseg=40 # segment length = 40*5 = 200ms
noverlap=10 # overlap = 10*5 = 50ms
start_time=1.5/0.02 # t returned by spectrogram
start_time_psth = 2000/5 
time_scale = 1000
burst_bins = np.arange(0,10.,0.1)

step = 1



def smooth_spectrogram(t,f,Sxx):
	func = interp2d(t, f, Sxx, kind='cubic')	
	t_new = np.arange(np.min(t),np.max(t),0.02)
	f_new = np.arange(np.min(f),np.max(f),1)
	Sxx_new = func(t_new,f_new)
	t_N, f_N = np.meshgrid(t_new,f_new)
	
	return t_N, f_N, Sxx_new


subfig_hands = []
smoothen_time_unit = 20. # ms, had to compare the figure time axis and had to figure this out 


burst_data = pickle.load(open("../ssbn_effect/burst_data_all_border.pickle","r"))
max_pow = burst_data["threshold"]
#max_pow = 0.7


fig = plt.figure(figsize=(16.5,12))
ax0 = plt.subplot2grid((3,3), (0,0),colspan=2) # PSTH
ax0_g = plt.subplot2grid((3,3), (1,0),colspan=2) # PSTH
ax1 = plt.subplot2grid((3,3), (0,2)) # colormap for bl
ax1_g = plt.subplot2grid((3,3), (1,2)) # colormap for ba
ax2_1 = plt.subplot2grid((3,3), (2,0)) # Burst width vs burst amplitude for a point in colormap
ax2_2 = plt.subplot2grid((3,3), (2,1)) # Burst width vs burst amplitude for a point in colormap
ax2_3 = plt.subplot2grid((3,3), (2,2)) # Burst width vs burst amplitude for a point in colormap

gpe_eg_bs = [gpe_bs, 0.4, 0.1]
stn_eg_bs = [stn_bs, 0.4, 0.8]

hands = [ax2_1, ax2_2, ax2_3]
labels = ["E", "F", "G"]
colors = ["cyan","coral", "teal"]

prefix = "psdbEffect_inh_ip_"+str(gpe_ip)+"_exc_ip_"+str(stn_ip)+"_inh_ratio_"+str(gpe_bs)+"_exc_ratio_"+str(stn_bs)
data = pickle.load(open("output/"+prefix+".pickle","r"))

j = 1

print "trial no:",j

# Read and plot PSTH
psth_stn = data["psth_stn"][0][j]
psth_gpe = data["psth_gpe"][0][j]
bins = 	data["psth_stn"][1]

# Generate psth that folows poisson distribution of the same rate as this combination
lam = np.mean(psth_stn)
size_poss = len(psth_stn)
psth_stn_poss = np.random.poisson(lam=lam,size=size_poss)


N  = 2    # Filter order
B, A = sp_sig.butter(N, [0.15,0.2], btype='band')
psth_stn_new = psth_stn	
psth_gpe_new = psth_gpe




# Calculate spectrogram from psth
f_stn,t_stn,Sxx_stn = sp_sig.spectrogram(psth_stn - np.mean(psth_stn),fs,nperseg=nperseg,noverlap=noverlap)
# Smoothing the spectrogram
t_stn_new,f_stn_new,Sxx_stn_new = smooth_spectrogram(t_stn,f_stn,Sxx_stn)	

# Calculate spectrogram from poisson psth
f_poss,t_poss,Sxx_poss = sp_sig.spectrogram(psth_stn_poss - np.mean(psth_stn_poss),fs,nperseg=nperseg,noverlap=noverlap)
# Smoothing the spectrogram
t_poss_new,f_poss_new,Sxx_poss_new = smooth_spectrogram(t_poss,f_poss,Sxx_poss)


psth_stn_filt = np.abs(sp_sig.filtfilt(B,A,psth_stn))

psth_stn_poss_filt = np.abs(sp_sig.filtfilt(B,A,psth_stn_poss))


ax0_g.plot(bins[:-1][start_time_psth::step]/1000.,psth_stn_filt[start_time_psth::step],'-',color='blue',linewidth=1.0,label='STN')
ax0_g.plot(bins[:-1][start_time_psth::step]/1000.,psth_stn_poss_filt[start_time_psth::step],'-',color='orange',linewidth=1.0,label='Poisson',alpha=0.5)
ax0_g.set_title("Filtered PSTH",fontsize=15,fontweight='bold')
for x in ax0_g.get_xticklabels():
	x.set_fontsize(10)
	x.set_fontweight('bold')
for x in ax0_g.get_yticklabels():
	x.set_fontsize(10)
	x.set_fontweight('bold')
ax0_g.legend(prop={'size':12,'weight':'bold'})
ax0_g.set_xlim(2.0,7.0)	

y_lab = 4.1
ax0_g.text(1.6,y_lab,"B",fontweight='bold',fontsize=15)


avg_STN_env = np.mean([ Sxx_stn_new[f_stn_new==x] for x in np.arange(15,20,1)],0)
avg_poss_env = np.mean([ Sxx_poss_new[f_poss_new==x] for x in np.arange(15,20,1)],0)

 
ax0.plot(t_stn_new[f_stn_new==20.][start_time::step],avg_STN_env[start_time::step],'-',color='steelblue',linewidth=2.0,label='STN')

ax0.plot(t_poss_new[f_poss_new==20.][start_time::step],avg_poss_env[start_time::step],'-',color='sienna',linewidth=2.0,label='Poisson',alpha=0.5)

print "max_pow",max_pow
ax0.hlines(xmin=np.min(t_poss_new[f_poss_new==20.][start_time::step]) ,xmax=np.max(t_poss_new[f_poss_new==20.][start_time::step]),y=max_pow,color='r',linestyles='dashed')	
ax0.legend(prop={'size':12,'weight':'bold'})
ax0.set_xlim(2.0,7.0)
ax0.set_title("Spectogram envelope",fontsize=15,fontweight='bold')
for x in ax0.get_xticklabels():
	x.set_fontsize(10)
	x.set_fontweight('bold')
for x in ax0.get_yticklabels():
	x.set_fontsize(10)
	x.set_fontweight('bold')

y_lab = 1.21
ax0.text(1.6,y_lab ,"A",fontweight='bold',fontsize=15)


gpe_rat = np.arange(0.,1.1,0.1)
stn_rat = np.arange(0.,1.1,0.1)

ax1.set_title("Burst length (s)",fontsize=15,fontweight='bold')
im1 = ax1.pcolor(gpe_rat,stn_rat,burst_data["bl_mean"],cmap=cm.inferno,vmin=0.0)
ax1.set_ylabel("GPe bursting ratio",fontsize=15,fontweight='bold')
for x in ax1.get_xticklabels():
	x.set_fontsize(10)
	x.set_fontweight('bold')
for x in ax1.get_yticklabels():
	x.set_fontsize(10)
	x.set_fontweight('bold')
ax1.plot(float(stn_bs)+0.05,float(gpe_bs)+0.05,'*',color=colors[0],markersize=15)
ax1.plot(float(stn_eg_bs[1])+0.05,float(gpe_eg_bs[1])+0.05,'*',color=colors[1],markersize=15)
ax1.plot(float(stn_eg_bs[2])+0.05,float(gpe_eg_bs[2])+0.05,'*',color=colors[2],markersize=15)
ax1.set_xticks(stn_rat)
ax1.set_yticks(stn_rat)

ax1.text(-0.15,1.05,"C",fontweight='bold',fontsize=15)
ax1.grid(color='k',which='both',linestyle='solid',linewidth=1.0)
fig.colorbar(im1,ax=ax1)

ax1_g.set_title("Burst amplitude",fontsize=15,fontweight='bold')
im1_g = ax1_g.pcolor(gpe_rat,stn_rat,burst_data["ba_mean"],cmap=cm.inferno,vmin=0.0)
ax1_g.set_xlabel("STN bursting ratio",fontsize=15,fontweight='bold')
ax1_g.set_ylabel("GPe bursting ratio",fontsize=15,fontweight='bold')
for x in ax1_g.get_xticklabels():
	x.set_fontsize(10)
	x.set_fontweight('bold')
for x in ax1_g.get_yticklabels():
	x.set_fontsize(10)
	x.set_fontweight('bold')

ax1_g.plot(float(stn_bs)+0.05,float(gpe_bs)+0.05,'*',color='cyan',markersize=15)
ax1_g.plot(float(stn_eg_bs[1])+0.05,float(gpe_eg_bs[1])+0.05,'*',color=colors[1],markersize=15)
ax1_g.plot(float(stn_eg_bs[2])+0.05,float(gpe_eg_bs[2])+0.05,'*',color=colors[2],markersize=15)

ax1_g.set_xticks(stn_rat)
ax1_g.set_yticks(stn_rat)

ax1_g.text(-0.15,1.05,"D",fontweight='bold',fontsize=15)
ax1_g.grid(color='k',which='both',linestyle='solid',linewidth=1.0)
fig.colorbar(im1_g,ax=ax1_g)



text_loc = [(0.9,0.585),(3.0,0.76), (1.4,0.7)]
text_loc1 = [(0.9,0.565),(3.0,0.725), (1.4,0.67)]
for i,(g1,s1) in enumerate(zip(gpe_eg_bs,stn_eg_bs)):

	gpe_ind = np.where(np.round(gpe_rat,1)==float(g1))[0]
	stn_ind = np.where(np.round(stn_rat,1)==float(s1))[0]

	bl = burst_data["bl_all"][gpe_ind][stn_ind]
	ba = burst_data["ba_all"][gpe_ind][stn_ind]
	avg_freq = burst_data["burst_freq_median"][gpe_ind,stn_ind][0]
	avg_bl = np.round(np.median(np.hstack(bl)),2)
	
	for x,y in zip(bl,ba):
		hands[i].plot(x,y,'.',markersize=10,color=colors[i])
	hands[i].set_xlabel("Burst length (s)",fontsize=15,fontweight='bold')

	if i == 0:
		hands[i].set_ylabel("Burst amplitude",fontsize=15,fontweight='bold')
	for x in hands[i].get_xticklabels():
		x.set_fontsize(10)
		x.set_fontweight('bold')
	for x in hands[i].get_yticklabels():
		x.set_fontsize(10)
		x.set_fontweight('bold')


	slope, intercept, r_value, p_value, std_err = sp.stats.linregress(np.hstack(bl), np.hstack(ba))
	print p_value
	z_check = slope*np.arange(0,np.max(np.hstack(bl)),0.1)+intercept


	hands[i].text(text_loc[i][0],text_loc[i][1],"Avg. burst length: "+str(avg_bl)+"s",fontsize=10,fontweight='bold',color='k')
	hands[i].text(text_loc[i][0],text_loc1[i][1],"Avg. intraburst freq: "+str(avg_freq)+"Hz",fontsize=10,fontweight='bold',color='k')

	hands[i].plot(np.arange(0,np.max(np.hstack(bl)),0.1),z_check,'-',color=colors[i],linewidth=2.0)
	if i == 0:
		hands[i].text(np.median(np.hstack(bl))+0.25,np.median(np.hstack(ba))-0.05,'r='+str(np.round(r_value,2))+', p='+str(np.round(p_value,4)),color='k',fontweight='bold')
	else:
		hands[i].text(np.mean(np.hstack(bl))+0.5,np.median(np.hstack(ba))-0.01,'r='+str(np.round(r_value,2))+', p='+str(np.round(p_value,4)),color='k',fontweight='bold')


	y_lab = hands[i].get_ylim()[1]+0.02*hands[i].get_ylim()[1]
	x_lab = hands[i].get_xlim()[1]*0.1

	if i == 0:
		hands[i].text(-x_lab*2.5,y_lab,labels[i],fontweight='bold',fontsize=15)
	else:
		hands[i].text(-x_lab*1.2,y_lab,labels[i],fontweight='bold',fontsize=15)
	if i == 0:
		hands[i].set_xlim(0,2.0)
	elif i == 2:
		hands[i].set_xlim(0,3.0)
	else:
		hands[i].set_xlim(0,7.0)


fig.subplots_adjust(hspace=0.2,left=0.06,right=0.96,bottom=0.05,top=0.92,wspace=0.2)	
fig.savefig("figs/burst_length_summary.pdf")
