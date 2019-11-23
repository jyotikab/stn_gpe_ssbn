
import numpy as np
import itertools
import sys
home_directory = ""
sys.path.append(home_directory+'/common/')
import analyze_data_ssbn as adata
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
import pickle




def comp_mean_rate(act,total_time,numNeurons,time_range=[]):
	evs = act[:,0]
	ts = act[:,1]
	binsize = 5.
	if time_range!=[]:
	    idx = (ts>time_range[0]) & (ts<=time_range[1])
	    spikes = ts[idx]
	if time_range==[]:
	  total_time =total_time 
	else:
	  total_time = time_range[1] - time_range[0]
	mean_rate = 1.*len(spikes)/(numNeurons*total_time*1e-3)
	return mean_rate
	#psth,xx = np.histogram(ts,bins = np.arange(0,total_time,binsize))
	#return np.mean(psth/((binsize/1000.)*numNeurons))

def psd(act,total_time, bin_w = 5.,time_range = []):
	evs = act[:,0]
	ts = act[:,1]

	if time_range!=[]:
	    idx = (ts>time_range[0]) & (ts<=time_range[1])
	    spikes = ts[idx]
	if time_range==[]:
	  total_time =total_time 
	else:
	  total_time = time_range[1] - time_range[0]

	if len(spikes) == 0:
	  print 'psd: spike array is empty'
	  return np.nan, np.nan, np.nan, np.nan
	  
	ids = np.unique(evs)
	nr_neurons = len(ids)
	#psd, max_value, freq,h = misc2.psd_sp(spikes[:,1],nr_bins,nr_neurons)
	bins = np.arange(time_range[0],time_range[1],bin_w)
	a,b = np.histogram(spikes, bins)
	ff = abs(np.fft.fft(a- np.mean(a)))**2
	Fs = 1./(bin_w*0.001)
	freq2 = np.fft.fftfreq(len(bins))[0:len(bins/2)+1]
	freq = np.linspace(0,Fs/2,len(ff)/2+1)
	px = ff[0:len(ff)/2+1]
	max_px = np.max(px[1:])
	idx = px == max_px
	corr_freq = freq[pl.find(idx)]
	new_px = px
	max_pow = new_px[pl.find(idx)]
	return new_px,freq, freq2, corr_freq[0], max_pow


def spec_entropy(act,total_time,bin_w = 5.,time_range=[],freq_range = []):
	'''Function to calculate the spectral entropy'''
	power,freq,dummy,dummy,dummy = psd(bin_w = bin_w,time_range = time_range,act=act,total_time=total_time)
	if freq_range != []:
	  power = power[(freq>freq_range[0]) & (freq < freq_range[1])]
	  freq = freq[(freq>freq_range[0]) & (freq < freq_range[1])]	
	k = len(freq)
	power = power/sum(power)
	sum_power = 0
	for ii in range(k):
	  sum_power += (power[ii]*np.log(power[ii]))
	spec_ent = -(sum_power/np.log(k))  
	return spec_ent


