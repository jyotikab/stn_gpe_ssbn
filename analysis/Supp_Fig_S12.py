
#--------------------------------------------------------------------------------------------------------------------
# Reproduces Supp Figure S12
#  Extra figure, average power spectrum for two phases: GPe turns bursty, STN turns bursty
# Call script with "border" 0.4 0.9
#
#-----------------------------------------------------------------------------------------------------------------------


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
import scipy.signal as sp_sig
from scipy.interpolate import interp2d

pl.ion()


def smooth_spectrogram(t,f,Sxx):
	func = interp2d(t, f, Sxx, kind='cubic')	
	t_new = np.arange(np.min(t),np.max(t),0.02)
	f_new = np.arange(np.min(f),np.max(f),1)
	Sxx_new = func(t_new,f_new)
	t_N, f_N = np.meshgrid(t_new,f_new)
	
	return t_N, f_N, Sxx_new

def smooth_psth(a,box):
	box_kern = np.ones(box)/box
	a_new = np.convolve(a,box_kern,mode='same')
	return a_new
	


regime = sys.argv[1]  # "border"
gpe_bs = sys.argv[2]  # 0.4
stn_bs = sys.argv[3]  # 0.9

if regime == "border": 
	gpe_ip = np.array([900.])[0]
	stn_ip = np.array([1400.])[0]
	max_pow = 0.09
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


binwidth = 10.
gpe_rat = gpe_bs 
stn_rat = stn_bs
pre = 200
post = 400
simtime = 7500
fs = 1000/binwidth # 1000/5ms
nperseg=50 # 
noverlap=10 # 
start_time=1000/binwidth
time_scale = 1000
ticks_step = 200
binsize=10.
path1 = '../ssbn_effect/'


prefix = "psdbEffect_midway_change_inh_ip_"+str(gpe_ip)+"_exc_ip_"+str(stn_ip)+"_inh_ratio_"+str(gpe_rat)+"_exc_ratio_"+str(stn_rat)+"_spikes"
data = np.load("/"+prefix+".npy")

#seed = np.random.randint(0,9999999,1)[0]
seed = 7527479
print "seed",seed
np.random.seed(seed)

def plot_rate_dist(ax,data,bins,txt):
	a,b,c = ax.hist(data,bins=bins,facecolor='darkolivegreen',normed=True)
	mean = np.round(np.mean(data),2)
	ax.text(60,np.max(a),"mean:"+str(mean)+" spks/s",color='k',fontsize=12,fontweight='bold')
	ax.set_title(txt,fontsize=15,fontweight='bold')
	return a

def pick_plot_neuron_ids(ids,data1,which,ax,y_ax):
	mean = np.mean(data1)
	print "mean",mean
	std = np.std(data1)
	print "std",std
	if which == "less":
		ind = np.where(data1<=np.abs(mean-std ))[0]
		np.random.shuffle(ind)
		
		ids_sub = ids[ind[0]]
		print data1[ind[0]]
		ax.text(data1[ind[0]],y_ax*0.75,"#1",color='k',fontsize=12,fontweight='bold')
		ax.vlines(x=data1[ind[0]],ymin=0,ymax=y_ax*0.7,color='k',linestyles='dashed',linewidth=2.0)
		return ids_sub
	elif which == "more":
		ind = np.where(data1>(mean+std))[0]
		np.random.shuffle(ind)
		
		ids_sub = ids[ind[0]]
		print data1[ind[0]]
		ax.text(data1[ind[0]],y_ax*0.2,"#3",color='k',fontsize=12,fontweight='bold')
		ax.vlines(x=data1[ind[0]],ymin=0,ymax=y_ax*0.15,color='k',linestyles='dashed',linewidth=2.0)
		return ids_sub
	else:
		ind = np.where(np.abs(data1-mean)<=1.0)[0]
		np.random.shuffle(ind)
		
		ids_sub = ids[ind[0]]
		print data1[ind[0]]
		ax.text(data1[ind[0]],y_ax*0.4,"#2",color='k',fontsize=12,fontweight='bold')
		ax.vlines(x=data1[ind[0]],ymin=0,ymax=y_ax*0.35,color='k',linestyles='dashed',linewidth=2.0)
		return ids_sub


stn_ids_b = np.arange(3,int(1000*float(stn_bs)+3))
stn_ids_nb = np.arange(int(1000*float(stn_bs))+3,1000+3)

max_rate = 400

stn_ids_b_rates = np.array([   len(np.where(data[:,0]==x)[0])/7.5  for x in stn_ids_b   ])
stn_ids_nb_rates = np.array([ len(np.where(data[:,0]==x)[0])/7.5  for x in stn_ids_nb ])

 
gpe_ids_b = np.arange(1000+3,1000+int(2000*float(gpe_bs))+3)
gpe_ids_nb = np.arange(1003+int(2000*float(gpe_bs)),3003)

gpe_ids_b_rates = np.array([ len(np.where(data[:,0]==x)[0])/7.5  for x in gpe_ids_b ])
gpe_ids_nb_rates = np.array([ len(np.where(data[:,0]==x)[0])/7.5  for x in gpe_ids_nb ])



bins_rates = np.linspace(0,120,50)

fig3 = pl.figure(3,figsize=(12,12))
t1 = fig3.add_subplot(221)
t2 = fig3.add_subplot(222)
t3 = fig3.add_subplot(223)
t4 = fig3.add_subplot(224)

a3 = plot_rate_dist(t1,gpe_ids_b_rates,bins_rates,"GPe bursty")
a4 = plot_rate_dist(t2,gpe_ids_nb_rates,bins_rates,"GPe non-bursty")
a1 = plot_rate_dist(t3,stn_ids_b_rates,bins_rates,"STN bursty")
a2 = plot_rate_dist(t4,stn_ids_nb_rates,bins_rates,"STN non-bursty")

gpe_b_3 = [pick_plot_neuron_ids(gpe_ids_b,gpe_ids_b_rates,"less",t1,np.max(a3)),pick_plot_neuron_ids(gpe_ids_b,gpe_ids_b_rates,"equal",t1,np.max(a3)),pick_plot_neuron_ids(gpe_ids_b,gpe_ids_b_rates,"more",t1,np.max(a3))]
gpe_nb_3 = [pick_plot_neuron_ids(gpe_ids_nb,gpe_ids_nb_rates,"less",t2,np.max(a4)),pick_plot_neuron_ids(gpe_ids_nb,gpe_ids_nb_rates,"equal",t2,np.max(a4)),pick_plot_neuron_ids(gpe_ids_nb,gpe_ids_nb_rates,"more",t2,np.max(a4))]
stn_b_3 = [pick_plot_neuron_ids(stn_ids_b,stn_ids_b_rates,"less",t3,np.max(a1)),pick_plot_neuron_ids(stn_ids_b,stn_ids_b_rates,"equal",t3,np.max(a1)),pick_plot_neuron_ids(stn_ids_b,stn_ids_b_rates,"more",t3,np.max(a1))]
stn_nb_3 = [pick_plot_neuron_ids(stn_ids_nb,stn_ids_nb_rates,"less",t4,np.max(a2)),pick_plot_neuron_ids(stn_ids_nb,stn_ids_nb_rates,"equal",t4,np.max(a2)),pick_plot_neuron_ids(stn_ids_nb,stn_ids_nb_rates,"more",t4,np.max(a2))]

fig3.savefig("figs/single_neuron_firing_rates_histogram.png")

# Pick 3 GPe bursy, 3 Gpe non bursty, 3 stn bursty, 3 STn non bursty


print "gpe_b_3",gpe_b_3
print "gpe_nb_3",gpe_nb_3


print "stn_b_3",stn_b_3
print "stn_nb_3",stn_nb_3



bins = np.arange(0,simtime,binwidth)

def plot_spectrogram(ax,ts,tit):
	a,b = np.histogram(ts,bins=bins)
	f,t,Sxx = sp_sig.spectrogram(a - np.mean(a),fs,nperseg=nperseg,noverlap=noverlap)
	t_new,f_new,Sxx_new = smooth_spectrogram(t,f,Sxx)
	im = ax.pcolormesh(t_new*time_scale+(start_time*binwidth),f_new,Sxx_new,shading='gouraud',vmin=min_pow,vmax=max_pow)	
	ax.set_title(tit,fontsize=15,fontweight='bold')
	ax.set_xlim((start_time*binwidth+t[0]*time_scale),np.max(t*time_scale))

	return f,t,Sxx


def calc_fft(psth,binw):
	fftfreq = np.fft.fftfreq(len(psth),d= binw)
	fft = np.fft.fft(psth-np.mean(psth))
	fft_y = np.abs(fft)[:len(psth)/2]
	# Smooth fft
	box_kern = np.ones(5)/5.
	fft_y_new = np.convolve(fft_y,box_kern,mode='same')

	return fftfreq,fft_y_new


def plot_fft(ax,fftfreq,fft_y_new,std_fft,color,label_flag,label,binw,tit):

	if label_flag == True:
		ax.plot(fftfreq[:len(fft_y_new)]*1000,fft_y_new,'-',color=color,alpha=0.5,linewidth=2.0,label=label)
		ax.fill_between(fftfreq[:len(fft_y_new)]*1000,fft_y_new-std_fft,fft_y_new+std_fft,color=color,alpha=0.1)
	else:
		ax.plot(fftfreq[:len(fft_y_new)]*1000,fft_y_new,'-',color=color,alpha=0.5,linewidth=2.0)
		ax.fill_between(fftfreq[:len(fft_y_new)]*1000,fft_y_new-std_fft,fft_y_new+std_fft,color=color,alpha=0.1)

	ax.set_title(tit,fontsize=15,fontweight='bold')
	if "GPe" in tit:
		ax.set_xticklabels([])
	else:
		for x in ax.get_xticklabels():
			x.set_fontweight('bold')
		ax.set_xlabel("Frequency (Hz)",fontweight='bold',fontsize=15)
	for x in ax.get_yticklabels():
		x.set_fontweight('bold')

	ax.legend()

fig = pl.figure(1,figsize=(12,20))
fig1 = pl.figure(2,figsize=(12,12))
subhands = []
sub_gpe_b = fig1.add_subplot(2,2,1)
sub_gpe_nb = fig1.add_subplot(2,2,2)
sub_stn_b = fig1.add_subplot(2,2,3)
sub_stn_nb = fig1.add_subplot(2,2,4)
binw = 10.
before_1500_colors = cm.get_cmap('Blues',5)
before_3500_colors = cm.get_cmap('Reds',5)
before_7500_colors = cm.get_cmap('Greens',5)

for i in np.arange(0,12):
	subhands.append(fig.add_subplot(4,3,i+1))
	if i < 3:
		if i == 0:
			b4_3500_fft = []
			b4_3500_freq = []
			b4_7500_fft = []
			b4_7500_freq = []

		ts = data[np.where(data[:,0]==gpe_b_3[i]),1]
		f,t,Sxx = plot_spectrogram(subhands[-1],ts,"GPe bursty#"+str(i+1))
		a1,b1 = np.histogram(ts[np.where(ts<=3500)],bins=np.arange(1000,3500,binw))
		t1,t2 = calc_fft(a1,binw)
		b4_3500_freq.append(t1)
		b4_3500_fft.append(t2)
	
		a2,b2 = np.histogram(ts[np.where(ts>3500)],bins=np.arange(3500,7500,binw))

		t1,t2 = calc_fft(a2,binw)
		b4_7500_freq.append(t1)
		b4_7500_fft.append(t2)
		if i == 2:
			plot_fft(sub_gpe_b,np.mean(b4_3500_freq,axis=0),np.mean(b4_3500_fft,axis=0),np.std(b4_3500_fft,axis=0),before_3500_colors(i+2),True,"1500-3500ms",binw,"GPe bursty")	
			plot_fft(sub_gpe_b,np.mean(b4_7500_freq,axis=0),np.mean(b4_7500_fft,axis=0),np.std(b4_7500_fft,axis=0),before_7500_colors(i+2),True,"3500-7500ms",binw,"GPe bursty")	

	if i > 2 and i < 6:	
		if i == 3:
			b4_3500_fft = []
			b4_3500_freq = []
			b4_7500_fft = []
			b4_7500_freq = []

		ts = data[np.where(data[:,0]==gpe_nb_3[i-3]),1]
		f,t,Sxx = plot_spectrogram(subhands[-1],ts,"GPe non-bursty#"+str(i-3+1))
		a1,b1 = np.histogram(ts[np.where(ts<=3500)],bins=np.arange(1000,3500,binw))
		t1,t2 = calc_fft(a1,binw)
		b4_3500_freq.append(t1)
		b4_3500_fft.append(t2)

		a2,b2 = np.histogram(ts[np.where(ts>3500)],bins=np.arange(3500,7500,binw))
		t1,t2 = calc_fft(a2,binw)
		b4_7500_freq.append(t1)
		b4_7500_fft.append(t2)
		if i == 5:
			plot_fft(sub_gpe_nb,np.mean(b4_3500_freq,axis=0),np.mean(b4_3500_fft,axis=0),np.std(b4_3500_fft,axis=0),before_3500_colors(i+2),True,"1500-3500ms",binw,"GPe non-bursty")	
			plot_fft(sub_gpe_nb,np.mean(b4_7500_freq,axis=0),np.mean(b4_7500_fft,axis=0),np.std(b4_7500_fft,axis=0),before_7500_colors(i+2),True,"3500-7500ms",binw,"GPe non-bursty")	


	if i > 5 and i < 9:	
		ts = data[np.where(data[:,0]==stn_b_3[i-6]),1]
		f,t,Sxx = plot_spectrogram(subhands[-1],ts,"STN bursty#"+str(i-6+1))
		if i == 6:
			b4_3500_fft = []
			b4_3500_freq = []
			b4_7500_fft = []
			b4_7500_freq = []
		a1,b1 = np.histogram(ts[np.where(ts<=3500)],bins=np.arange(1000,3500,binw))
		t1,t2 = calc_fft(a1,binw)
		b4_3500_freq.append(t1)
		b4_3500_fft.append(t2)

		a2,b2 = np.histogram(ts[np.where(ts>3500)],bins=np.arange(3500,7500,binw))
		t1,t2 = calc_fft(a2,binw)
		b4_7500_freq.append(t1)
		b4_7500_fft.append(t2)
		if i == 8:
			plot_fft(sub_stn_b,np.mean(b4_3500_freq,axis=0),np.mean(b4_3500_fft,axis=0),np.std(b4_3500_fft,axis=0),before_3500_colors(i+2),True,"1500-3500ms",binw,"STN bursty")	
			plot_fft(sub_stn_b,np.mean(b4_7500_freq,axis=0),np.mean(b4_7500_fft,axis=0),np.std(b4_7500_fft,axis=0),before_7500_colors(i+2),True,"3500-7500ms",binw,"STN bursty")	

	if i > 8 and i < 12:	
		ts = data[np.where(data[:,0]==stn_nb_3[i-9]),1]
		f,t,Sxx = plot_spectrogram(subhands[-1],ts,"STN non-bursty#"+str(i-9+1))
		if i == 9:
			b4_3500_fft = []
			b4_3500_freq = []
			b4_7500_fft = []
			b4_7500_freq = []
		a1,b1 = np.histogram(ts[np.where(ts<=3500)],bins=np.arange(1000,3500,binw))
		t1,t2 = calc_fft(a1,binw)
		b4_3500_freq.append(t1)
		b4_3500_fft.append(t2)

		a2,b2 = np.histogram(ts[np.where(ts>3500)],bins=np.arange(3500,7500,binw))
		t1,t2 = calc_fft(a2,binw)
		b4_7500_freq.append(t1)
		b4_7500_fft.append(t2)
		if i == 11:
			plot_fft(sub_stn_nb,np.mean(b4_3500_freq,axis=0),np.mean(b4_3500_fft,axis=0),np.std(b4_3500_fft,axis=0),before_3500_colors(i+2),True,"1500-3500ms",binw,"STN non-bursty")	
			plot_fft(sub_stn_nb,np.mean(b4_7500_freq,axis=0),np.mean(b4_7500_fft,axis=0),np.std(b4_7500_fft,axis=0),before_7500_colors(i+2),True,"3500-7500ms",binw,"STN non-bursty")	


	for x in subhands[-1].get_yticklabels():
		x.set_fontsize(10)
		x.set_fontweight("bold")


	if i <9:
		subhands[-1].set_xticklabels([],visible=False)
	else:
		for x in subhands[-1].get_xticklabels():
			x.set_fontsize(10)
			x.set_fontweight('bold')	
		subhands[-1].set_xlabel("Time",fontsize=15,fontweight='bold')
plt.figure(1)
plt.tight_layout()

plt.figure(2)
plt.tight_layout()

fig.savefig("figs/single_neuron_spectograms.png")
fig1.savefig("figs/single_neuron_ffts.png")


