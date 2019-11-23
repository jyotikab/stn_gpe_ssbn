import sys
import numpy as np
import pylab as pl
import sys
home_directory = ""
sys.path.append(home_directory+'/common/')
import params_d_ssbn as params_d
import cPickle as cp
import misc2
import scipy.stats as st
from scipy.stats import norm
import math
from collections import Counter
import random
import pdb



# Class with functions to analyze simulated data

class analyze_data(object):
	    
    pars = params_d.Parameters([])

    def __init__(self,prefix,netPars,seed,data_path=''):	
    	if not data_path=='':
	    #self.data_path = pars.data_path
	    self.data_path = data_path
	self.pars.prefix = prefix
	# Also copy the GIDs of excitatory and inhibitory pops, to plot them in separate colors
	if "pops_exc" in dir(netPars):
		self.pops_exc = netPars.pops_exc
	if "pops_inh" in dir(netPars):
		self.pops_inh = netPars.pops_inh
	self.pars.T_total = netPars.T_total
	self.seed = seed
	self.events = {}

    def get_pop_spikes(self,spikes,nmax,pop_id='all'):
	
	if pop_id == 'all':
	    #neuron_nr = self.pars['neurons_tot']
	    neuron_nr = len(self.pars.pops_exc)+len(self.pars.pops_inh)
	else:
	    if pop_id==[]:
		print 'population not found'
		return [],[],[]
	    else:
		#pop_id = pop_id[0]
		if nmax==[]:
		    idx = (spikes[:,0]>=self.pars[pop_id][0]) & (spikes[:,0]<=self.pars[pop_id][-1])		    
		else: 
		    idx = (spikes[:,0]>=pop_id[0]) & (spikes[:,0]<pop_id[0]+nmax)
		    neuron_nr = nmax
		neuron_nr = min(pop_id[-1] - pop_id[0]+1,nmax)
		spikes = spikes[idx]	
		
	return pop_id,spikes, neuron_nr
	    
    def load_info(self,prefix,data_path=[]):
	fname = self.data_path + prefix + '.info'
	fh = open(fname,'r')
	info1 = cp.load(fh)
	return(info1)		    
	
	
    def load_spikes(self):
	if not 'spikes' in self.events.keys():
	    fname = self.pars.data_path +self.pars.prefix+'_spikes.npy'
	    print "fname",fname
	    spikes = np.load(fname)	
	    self.events['spikes'] = spikes
	return 
		

    def load_vm(self):
	if not 'vm' in self.events.keys():
	    fname = self.data_path + self.pars['prefix'] + '_vm.npy'
	    vm = np.load(fname)	
	    self.events['vm'] = vm
	    if len(vm)>0:
		print 'Imported Vm data'
	    else:
		print 'Vm array is empty !'
	return 
	
    def plot_raster(self,nmax=100,kernel_w=1):
	''' plot raster of spikes for all neurons from 1-nr'''
	self.load_spikes()
	spikes = self.events['spikes']
	if len(spikes)==0:
	    print 'Plot raster: spike array is empty !'
	    return	
	# JB: Separate the spikes of the two populations 
	# Assuming the GIDs assigned are contigous, if thats not the case, this will fail

	idxpopE =np.logical_and(spikes[:,0]>=self.pops_exc[0],spikes[:,0]<=self.pops_exc[-1])
	idxpopI =np.logical_and(spikes[:,0]>=self.pops_inh[0],spikes[:,0]<=self.pops_inh[-1])
	
	
	pl.plot(spikes[idxpopE,1],spikes[idxpopE,0],'r.',label='STN')
	pl.plot(spikes[idxpopI,1],spikes[idxpopI,0],'b.',label='GPe')
	pl.legend()
	pl.show()    
    
    def plot_vm(self,nid=[]):
	self.load_vm()
	self.load_spikes()	
	vm = self.events['vm']
	spikes = self.events['spikes']
	if nid!=[]:
	    idx1 = vm[:,0] ==nid
	    idx2 = spikes[:,0] ==nid
	else:
	    idx1 = vm[:,0]>0
	    idx2 = idx1
	offset = np.max(vm[idx1,2])*0.8
	pl.plot(vm[idx1,1],vm[idx1,2])
	pl.plot(spikes[idx2,1],spikes[idx2,0]*0+offset,'ro')		
	pl.show()
    
    #def comp_psth1(self,fig,pop_id='all',nmax=100,form = '-',kernel_w=1,res=0.1,plot_fl=1, label = 'one',kernel = 'normal', time_range = [], color=np.array([1.,1.,1.])/255.,binsize=100.):	
    def comp_psth1(self,pop_id='all',nmax=100,form = '-',kernel_w=1,res=0.1,plot_fl=1, label = 'one',kernel = 'normal', time_range = [], color=np.array([1.,1.,1.])/255.,binsize=50.):	
	''' compute psth'''
	self.load_spikes()	
	spikes = self.events['spikes']
	if len(spikes)==0:
	    print 'Comp psth: spike array is empty!'
	    return [],[]		
	
	pop_id,spikes,neuron_nr = self.get_pop_spikes(spikes,nmax,pop_id)
	print neuron_nr

	if spikes == []:
	    return [],[]	    		
	if time_range!=[]:
	    idx = (spikes[:,1]>time_range[0]) & (spikes[:,1]<=time_range[1])
	    spikes = spikes[idx]
	    sim_time = time_range[1] - time_range[0]
	
	else:
	    sim_time = self.pars.T_total	

	psth,xx = np.histogram(spikes[:,1],bins = np.arange(0,sim_time,binsize))
	
	#fig.plot(xx[:-1],psth/((binsize/1000.)*neuron_nr),'-',linewidth=2.0,color=color)
	return psth/((binsize/1000.)*neuron_nr),xx
	

    def comp_psth(self,pop_id='all',nmax=100,form = '-',kernel_w=1,res=0.1,plot_fl=1, label = 'one',kernel = 'normal', time_range = [], color=np.array([1.,1.,1.])/255.):
	''' compute psth'''
	self.load_spikes()	
	spikes = self.events['spikes']
	if len(spikes)==0:
	    print 'Comp psth: spike array is empty!'
	    return [],[]		
	
	pop_id,spikes,neuron_nr = self.get_pop_spikes(spikes,nmax,pop_id)
	print neuron_nr

	if spikes == []:
	    return [],[]	    		
	if time_range!=[]:
	    idx = (spikes[:,1]>time_range[0]) & (spikes[:,1]<=time_range[1])
	    spikes = spikes[idx]
	    sim_time = time_range[1] - time_range[0]
	
	else:
	    sim_time = self.pars.T_total	

	
	
	if plot_fl:
	    pl.plot(xx,psth,form,lw = 3.,label = label, color = color)
	    pl.xlim([0,self.pars.T_sim + self.pars.T_wup])
	return(psth,xx)
	
            	
    def comp_mean_rate(self,time_range=[],pop_id='all',nmax=100):
	self.load_spikes()	
	spikes = self.events['spikes']
	if len(spikes)==0:
	    print 'Comp mean rate: spike array is empty !'
	    return np.nan
	pop_id,spikes,neuron_nr = self.get_pop_spikes(spikes,nmax,pop_id)
	if len(spikes) == 0:
	    return np.nan
	    
	if time_range!=[]:
	    idx = (spikes[:,1]>time_range[0]) & (spikes[:,1]<=time_range[1])
	    spikes = spikes[idx]
	if time_range==[]:
	  total_time = self.pars['T_total']
	else:
	  total_time = time_range[1] - time_range[0]
	mean_rate = 1.*len(spikes)/neuron_nr/total_time*1e3
	return mean_rate
	    
      
    # To compute the power spectral density and return the two biggest values 
    def psd(self, bin_w = 100., nmax = 4000, time_range = [], pop_id = 'all'):
      self.load_spikes()
      spikes = self.events['spikes']
      pop_id,spikes,neuron_nr = self.get_pop_spikes(spikes,nmax,pop_id)
      if len(spikes) == 0:
	print 'psd: spike array is empty'
	return np.nan, np.nan, np.nan, np.nan
      if time_range!=[]:
	    idx = (spikes[:,1]>time_range[0]) & (spikes[:,1]<=time_range[1])
	    spikes = spikes[idx]
      if time_range==[]:
	total_time = self.pars.T_total
      else:
	total_time = time_range[1] - time_range[0]
      
      if len(spikes) == 0:
	print 'psd: spike array is empty'
	return np.nan, np.nan, np.nan, np.nan
	
      ids = np.unique(spikes[:,0])[:nmax]
      nr_neurons = len(ids)
      #psd, max_value, freq,h = misc2.psd_sp(spikes[:,1],nr_bins,nr_neurons)
      if time_range==[]:
		bins = np.arange(self.pars.T_wup,self.pars.T_total,bin_w)	
      else:
		bins = np.arange(time_range[0],time_range[1],bin_w)
      a,b = np.histogram(spikes[:,1], bins)
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
      
	
    def spec_entropy(self,nmax,bin_w = 5.,pop_id = 'pops_exc',time_range=[],freq_range = []):
      '''Function to calculate the spectral entropy'''
      power,freq,freq2,dummy,dummy = self.psd(pop_id = pop_id,bin_w = bin_w,time_range = time_range,nmax=nmax)
      #print freq	
      if freq_range != []:
	power = power[(freq>freq_range[0]) & (freq < freq_range[1])]
	freq = freq[(freq>freq_range[0]) & (freq < freq_range[1])]	
      k = len(freq)
      power = power/sum(power)
      sum_power = 0
      for ii in range(k):
	sum_power += (power[ii]*np.log(power[ii]))
      spec_ent = -(sum_power/np.log(k))  
      return spec_ent,power,freq,freq2
      
      
      
	  
	



	
            

	
	    
