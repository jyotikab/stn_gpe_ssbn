

#-----------------------------------------------------------------------------------------
# Pre requisite for Figure 6 and Fig6.py
# Generates the file burst_data_all_border.pickle and additional figures
#
#-------------------------------------------------------------------------------------------




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

# "border", "osc" or "non-osc"
regime = sys.argv[1]


if regime == "border":
	gpe_ip = np.array([900.])[0]
	stn_ip = np.array([1400.])[0]
	max_pow = 1 # also serves as threshold 
	min_pow = 0


gpe_bs = np.arange(0.0,1.1,0.1)
stn_bs = np.arange(0.0,1.1,0.1)

bin_width=5 #ms
num_trials = 5

fs = 200 # 1000/5ms
nperseg=40 # segment length = 40*5 = 200ms
noverlap=10 # overlap = 10*5 = 50ms
start_time=1000/bin_width
time_scale = 1000
burst_bins = np.arange(0,10.,0.1) # in seconds
freq_bins = np.arange(5.,40.,1.) # in Hz



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


burst_data_all = dict()

burst_length_mean = np.zeros((len(gpe_bs),len(stn_bs)))
burst_freq_mean = np.zeros((len(gpe_bs),len(stn_bs)))
burst_length_median = np.zeros((len(gpe_bs),len(stn_bs)))
burst_amplitude_mean = np.zeros((len(gpe_bs),len(stn_bs)))
burst_amplitude_median = np.zeros((len(gpe_bs),len(stn_bs)))



per_poss_trials=[]
for i,gpe_rat in enumerate(gpe_bs):
	for j,stn_rat in enumerate(stn_bs):
		prefix = "psdbEffect_inh_ip_"+str(gpe_ip)+"_exc_ip_"+str(stn_ip)+"_inh_ratio_"+str(gpe_rat)+"_exc_ratio_"+str(stn_rat)
		data = pickle.load(open("output/"+prefix+".pickle","r"))
		for k in xrange(0,num_trials):
			psth_stn = data["psth_stn"][0][k]

			# Generate psth that folows poisson distribution of the same rate as this combination
			lam = np.mean(psth_stn)
			size_poss = len(psth_stn)
			psth_stn_poss = np.random.poisson(lam=lam,size=size_poss)

			# Calculate spectrogram from poisson psth
			f_poss,t_poss,Sxx_poss = sp_sig.spectrogram(psth_stn_poss - np.mean(psth_stn_poss),fs,nperseg=nperseg,noverlap=noverlap)
			# Smoothing the spectrogram
			t_poss_new,f_poss_new,Sxx_poss_new = smooth_spectrogram(t_poss,f_poss,Sxx_poss)	

			Sxx_poss_new_beta_band = np.mean([ Sxx_poss_new[f_poss_new == x] for x in np.arange(15.,20.,1.)],0)
			per_max = np.max(Sxx_poss_new_beta_band[50:]) # Because the dominant frequency changes
			per_poss_trials.append(per_max)


max_pow = np.mean(per_poss_trials)
burst_data_all["threshold"] = max_pow
burst_data_all["per_poss_trials"] = per_poss_trials

burst_lengths_all=[]
burst_amplitudes_all =[]
for i,gpe_rat in enumerate(gpe_bs):

	fig = plt.figure(figsize=(16,12))
	dominant_freq = []
	burst_amplitudes = []	
	burst_lengths = []
	fig.suptitle("Regime:"+regime+" GPe bs:"+str(gpe_rat*100)+"%",fontsize=15,fontweight='bold')
	ax1_hands=[]
	ax2_hands=[]
	ax3_hands=[]
	for j,stn_rat in enumerate(stn_bs):
		prefix = "psdbEffect_inh_ip_"+str(gpe_ip)+"_exc_ip_"+str(stn_ip)+"_inh_ratio_"+str(gpe_rat)+"_exc_ratio_"+str(stn_rat)
		data = pickle.load(open("../ssbn_effect/"+prefix+".pickle","r"))
		
		dominant_freq_stn_bs = []
		burst_amplitudes_stn_bs = []
		burst_lengths_stn_bs = []
		
		ax1 = fig.add_subplot(11,3,j*3+1) # Burst amplitude vs burst length
		ax2 = fig.add_subplot(11,3,j*3+2) # Burst length hist
		ax3 = fig.add_subplot(11,3,j*3+3) # Dominant frequency hist
		ax1_hands.append(ax1)
		ax2_hands.append(ax2)
		ax3_hands.append(ax3)
		
		for k in xrange(0,num_trials):
			psth_stn = data["psth_stn"][0][k]
			psth_gpe = data["psth_gpe"][0][k]
			bins = 	data["psth_stn"][1]
			
 			# Generate psth that folows poisson distribution of the same rate as this combination
			lam = np.mean(psth_stn)
			size_poss = len(psth_stn)
			psth_stn_poss = np.random.poisson(lam=lam,size=size_poss)


			# Calculate spectrogram from psth
			f_stn,t_stn,Sxx_stn = sp_sig.spectrogram(psth_stn - np.mean(psth_stn),fs,nperseg=nperseg,noverlap=noverlap)
			# Smoothing the spectrogram
			t_stn_new,f_stn_new,Sxx_stn_new = smooth_spectrogram(t_stn,f_stn,Sxx_stn)	
			Sxx_stn_new_beta_band = np.mean([ Sxx_stn_new[f_stn_new == x] for x in np.arange(15.,20.,1.)],0)
			f_stn_new_beta_band = np.mean([ f_stn_new[f_stn_new == x] for x in np.arange(15.,20.,1.)],0)
			ind = np.where(Sxx_stn_new_beta_band > max_pow)	
			
			ind1 = np.where(Sxx_stn_new > max_pow) 
			
			# ignoring the initial transient bit
			ind_mid = ind[0][ind[0]>50]
			ind_mid1 = ind1[1][ind1[1]>50]
			ind_freq = ind1[0][ind1[1]>50]
			
			# similar time slots - at different frequencies, 19,20,21, only frequency is enough
			
			grp_freq = [(k, sum(1 for i in g)) for k,g in groupby(ind_freq)]
			if len(grp_freq) > 0:
				dom_freq = np.array(grp_freq)[np.argmax(np.array(grp_freq)[:,1]),0]
				ind_mid_dom = ind_mid1[ind_freq==dom_freq]
				
				dominant_freq_stn_bs.append(dom_freq)

				diff_ind = np.diff(ind_mid_dom)
				
				bursts = [(k, sum(1 for i in g)) for k,g in groupby(diff_ind) if k == 1]
				grp_ind = [(k, sum(1 for i in g)) for k,g in groupby(diff_ind)]
				
				burst_amp =[ np.mean(Sxx_stn_new[ind_freq[m:m+x[1]],ind_mid1[m:m+x[1]]])	for m,x in enumerate(grp_ind) if x[0] ==1 and x[1] > 1]	

				burst_length=[]
				for x in bursts:
					if x[1] > 1:
						burst_length.append(x[1]*smoothen_time_unit)			
				if len(burst_amp) > 0:
					burst_amplitudes_stn_bs.append(burst_amp)
					burst_lengths_stn_bs.append(np.array(burst_length)/1000.)
				else:
					burst_amplitudes_stn_bs.append([0.0])
					burst_lengths_stn_bs.append([0.0])
				
			else:
				burst_amplitudes_stn_bs.append([0.0])
				burst_lengths_stn_bs.append([0.0])
			

		burst_amplitude_mean[i][j] = np.mean(np.hstack(burst_amplitudes_stn_bs))
		burst_amplitude_median[i][j] = np.median(np.hstack(burst_amplitudes_stn_bs))
		burst_length_mean[i][j] = np.mean(np.hstack(burst_lengths_stn_bs))
		burst_length_median[i][j] = np.median(np.hstack(burst_lengths_stn_bs))
		burst_freq_mean[i][j] = np.median(dominant_freq_stn_bs)	

		for x,y in zip(burst_lengths_stn_bs,burst_amplitudes_stn_bs):
			ax1.plot(x,y,'.',markersize=10,color='darkslategray')
		if j == 10:
			ax1.set_xlabel("Burst length (s)",fontsize=15,fontweight='bold')
		else:
			for x in ax1.get_xticklabels():
				x.set_visible(False)				
		if j == 5:
			ax1.set_ylabel("Burst amplitude",fontsize=15,fontweight='bold')	

		if j == 0 or j == 10:
			ax1.set_ylabel("STN bs:"+str(stn_rat*100)+"%",fontsize=10,fontweight='bold')


		bl,_,_ = ax2.hist(burst_lengths_stn_bs,bins=burst_bins,fc='darkolivegreen',alpha=0.5)
		min_bl = np.min(bl)
		max_bl = np.max(bl)
		ax2.set_yticks(np.round(np.linspace(min_bl,max_bl,3),1))
		ax2.set_yticklabels([ str(np.round(x,1)) for x in np.linspace(min_bl,max_bl,3)])
		mean = np.round(np.mean(np.hstack(burst_lengths_stn_bs)),1)
		ax2.text(5,(min_bl+max_bl)/2.,"Mean:"+str(mean),fontsize=10,fontweight='bold')

		if j == 10:
			ax2.set_xlabel("Burst length (s)",fontsize=15,fontweight='bold')
		else:
			for x in ax2.get_xticklabels():
				x.set_visible(False)				
			

		fr,_,_ = ax3.hist(dominant_freq_stn_bs,bins=freq_bins,fc='dodgerblue',alpha=0.5)
		min_fr = np.min(fr)
		max_fr = np.max(fr)
		ax3.set_yticks(np.round(np.linspace(min_fr,max_fr,3),1))
		ax3.set_yticklabels([ str(np.round(x,1)) for x in np.linspace(min_fr,max_fr,3)])
		mean = np.round(np.mean(np.hstack(dominant_freq_stn_bs)),1)
		ax3.text(30,(min_fr+max_fr)/2.,"Mean:"+str(mean),fontsize=10,fontweight='bold')


		if j == 10:
			ax3.set_xlabel("Dominant frequency (Hz)",fontsize=15,fontweight='bold')
		else:
			for x in ax3.get_xticklabels():
				x.set_visible(False)				
		
		
		burst_amplitudes.append(burst_amplitudes_stn_bs)
		burst_lengths.append(burst_lengths_stn_bs)
		dominant_freq.append(dominant_freq_stn_bs)
	if type(np.min(burst_amplitudes)) == list:
		min_amp = np.min(burst_amplitudes)[0]
	else:
		min_amp = np.min(burst_amplitudes)
	if type(np.max(burst_amplitudes)) == list:
		max_amp = np.max(burst_amplitudes)[0]
	else:
		max_amp = np.max(burst_amplitudes)


	for x in ax1_hands:
		x.set_ylim(min_amp,max_amp)
		x.set_yticks(np.round(np.linspace(min_amp,max_amp,3),1))
		x.set_yticklabels([ str(np.round(x1,1)) for x1 in np.linspace(min_amp,max_amp,3)])
		x.set_xlim(0.0,6.5)


	
	burst_lengths_all.append(burst_lengths)
	burst_amplitudes_all.append(burst_amplitudes)

	fig.subplots_adjust(hspace=0.2,left=0.08,right=0.96,bottom=0.05,top=0.92)	
	fig.savefig("figs/"+regime+"gpe_bs_"+str(gpe_rat)+"_burst_length_hist.png")



burst_data_all["bl_mean"] = burst_length_mean
burst_data_all["bl_median"] = burst_length_median
burst_data_all["ba_mean"] = burst_amplitude_mean
burst_data_all["ba_median"] = burst_amplitude_median
burst_data_all["burst_freq_median"] = burst_freq_mean
burst_data_all["bl_all"] = burst_lengths_all
burst_data_all["ba_all"] = burst_amplitudes_all
burst_data_all["corr_r_val"] = np.zeros((10,10))
burst_data_all["corr_p_val"] = np.zeros((10,10))



for i in np.arange(0,10): # GPe_bs
	for j in np.arange(0,10):
	        bl = burst_data_all["bl_all"][i][j]
	        ba = burst_data_all["ba_all"][i][j]

	        slope, intercept, r_value, p_value, std_err = sp.stats.linregress(np.hstack(bl), np.hstack(ba))	
		burst_data_all["corr_r_val"][i][j] = r_value
		burst_data_all["corr_p_val"][i][j] = p_value




for i in np.arange(0,10): # GPe_bs
	for j in np.arange(0,10):
		test=0
		if burst_data_all["bl_median"][i][j]>=0.1 and burst_data_all["bl_median"][i][j]<=0.3:
			test = test+1
		if burst_data_all["burst_freq_median"][i][j] >=18 and burst_data_all["burst_freq_median"][i][j]<=20:
			test = test+1
		if burst_data_all["corr_p_val"][i][j] < 0.05:
			test = test+1
		if test > 2:
			print "All conditions met"
			print burst_data_all["bl_median"][i][j]
			print burst_data_all["burst_freq_median"][i][j]
			print burst_data_all["corr_p_val"][i][j]
			print "GPe bs",i*10
			print "STN bs",j*10
			

pickle.dump(burst_data_all,open("../ssbn_effect/burst_data_all_"+regime+".pickle","w"))

