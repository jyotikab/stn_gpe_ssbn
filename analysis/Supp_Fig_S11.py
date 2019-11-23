
#----------------------------------------------------------------------------------------------------------------------------------------------
# Reproduces the Supp Fig S11 in the manuscript, and additional single trial spectograms
#--------------------------------------------------------------------------------------------------------------------------------------------
#
#


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
from scipy.interpolate import interp2d

pl.ion()
which = sys.argv[1]  #"midway_change_sinstn"
regime = sys.argv[2] #"border"
gpe_bs = sys.argv[3] # "0.4"

if regime == "border": 
	gpe_ip = np.array([900.])[0]
	stn_ip = np.array([1400.])[0]
	max_pow = 5.0
	min_pow = 0

elif regime == "osc":
	gpe_ip = np.array([500.])[0]
	stn_ip = np.array([1000.])[0]
	max_pow = 40
	min_pow = 0
elif regime == "non-osc":
	gpe_ip = np.array([1300.])[0]
	stn_ip = np.array([1400.])[0]
	max_pow = 2.0
	min_pow = 0


gpe_rat = gpe_bs 
bin_width=5
stn_bs = np.arange(0,1.1,0.1)
num_trials = 5
fs = 200 # 1000/5ms
nperseg=40 # segment length = 40*5 = 200ms
noverlap=10 # overlap = 10*5 = 50ms
start_time=1000/bin_width
time_scale = 1000


sxx_nb_all=[]
sxx_b_all=[]


def smooth_spectrogram(t,f,Sxx):
	func = interp2d(t, f, Sxx, kind='cubic')	
	t_new = np.arange(np.min(t),np.max(t),0.02)
	f_new = np.arange(np.min(f),np.max(f),1)
	Sxx_new = func(t_new,f_new)
	t_N, f_N = np.meshgrid(t_new,f_new)
	
	return t_N, f_N, Sxx_new




for j in xrange(0,num_trials):
	fig = pl.figure(figsize=(12,12))
	subfig_hands_b = []

	sxx_b_trial = []
	for i,stn_rat in enumerate(stn_bs):
		 
		prefix = "psdbEffect_"+which+"_inh_ip_"+str(gpe_ip)+"_exc_ip_"+str(stn_ip)+"_inh_ratio_"+str(gpe_rat)+"_exc_ratio_"+str(stn_rat)

		data = pickle.load(open("output/"+prefix+".pickle","r"))
		subfig_hands_b.append(fig.add_subplot(6,2,i+1))
			
		psth_stn_b = data["psth_stn"][0][j]
		f_b,t_b,Sxx_b = sp_sig.spectrogram(psth_stn_b[start_time:] - np.mean(psth_stn_b[start_time:]),fs,nperseg=nperseg,noverlap=noverlap)
		# Smoothing the spectrogram
		t_b_new,f_b_new,Sxx_b_new = smooth_spectrogram(t_b,f_b,Sxx_b)	

		im = subfig_hands_b[-1].pcolormesh(t_b_new*time_scale+(start_time*bin_width),f_b_new,Sxx_b_new,shading='gouraud',vmin=min_pow,vmax=max_pow)
		if i ==5:
			subfig_hands_b[-1].set_ylabel("Frequency (Hz)",fontsize=15,fontweight='bold')
		subfig_hands_b[-1].set_title("STN osc:"+str(stn_rat*100)+"%",fontsize=10,fontweight='bold')
		
		sxx_b_trial.append(Sxx_b_new)
		subfig_hands_b[-1].set_yticks([5,20,35])

		for x in subfig_hands_b[-1].get_yticklabels():
			x.set_fontsize(10)
			x.set_fontweight("bold")


		if i <10:
			subfig_hands_b[-1].set_xticklabels([],visible=False)
		else:
			for x in subfig_hands_b[-1].get_xticklabels():
				x.set_fontsize(10)
				x.set_fontweight('bold')	
			subfig_hands_b[-1].set_xlabel("Time",fontsize=15,fontweight='bold')
	
		subfig_hands_b[-1].set_ylim(5,40)
		subfig_hands_b[-1].set_xlim((start_time*bin_width+t_b[0]*time_scale),np.max(t_b*time_scale))
		subfig_hands_b[-1].set_aspect(40)
	
	sxx_b_all.append(sxx_b_trial)
	cbar_ax = fig.add_axes([0.6, 0.1, 0.3, 0.03])
	fig.colorbar(im, cax=cbar_ax,orientation='horizontal')	
	for x in cbar_ax.get_xticklabels():
		x.set_fontsize(10)
		x.set_fontweight('bold')		
	plt.tight_layout()
	fig.savefig("figs/"+which+"_"+regime+"spectrogram_gpe_rat_"+str(gpe_rat)+"_trial_"+str(j)+".png")




fig = pl.figure(figsize=(12,12)) 
for i,stn_rat in enumerate(stn_bs):
	ft_b = fig.add_subplot(6,2,i+1)

	mean_sxx_b = np.mean(np.array(sxx_b_all)[:,i,:,:],axis=0) # Because for 0th element, there was no element in bursting
	
	im = ft_b.pcolormesh(t_b_new*time_scale+(start_time*bin_width),f_b_new,mean_sxx_b,shading='gouraud',vmin=min_pow,vmax=max_pow)
	ft_b.set_title("STN osc:"+str(stn_rat*100)+"%",fontsize=12,fontweight='bold')

	ft_b.set_yticks([5,20,35])

	for x in ft_b.get_yticklabels():
		x.set_fontsize(12)
		x.set_fontweight("bold")

			
	if i <10:
		ft_b.set_xticklabels([],visible=False)
	else:
		for x in ft_b.get_xticklabels():
			x.set_fontsize(12)
			x.set_fontweight('bold')	
		ft_b.set_xlabel("Time",fontsize=15,fontweight='bold')
	if i ==5:
		ft_b.set_ylabel("Frequency (Hz)",fontsize=15,fontweight='bold')

	ft_b.set_ylim(5,40)
	ft_b.set_xlim((start_time*bin_width+t_b[0]*time_scale),np.max(t_b*time_scale))

	ft_b.set_aspect(40)

cbar_ax = fig.add_axes([0.6, 0.1, 0.3, 0.03])
fig.colorbar(im, cax=cbar_ax,orientation='horizontal')	

for x in cbar_ax.get_xticklabels()[::2]:
	x.set_fontsize(12)
	x.set_fontweight('bold')

for x in cbar_ax.get_xticklabels()[1::2]:
	x.set_visible(False)

fig.suptitle("Mean spectrogram over 5 trials",fontsize=15,fontweight='bold')
fig.subplots_adjust(hspace=0.2,left=0.05,right=0.96,bottom=0.05,top=0.94,wspace=0.2)	
fig.savefig("figs/"+which+"_"+regime+"spectrogram_gpe_rat_"+str(gpe_rat)+"_mean.png")

pl.show()


