import numpy as np
import itertools
import sys
import sys
import pdb
import pylab as pl
import matplotlib.cm as cm
import os
import pickle


# Had to change function spectral entropy to accept FFT spectrum instead of spike trains
def spec_entropy(power,freq,time_range=[],freq_range = []):
	'''Function to calculate the spectral entropy'''
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

# Spectral entropy vs noise
dt = 0.1
time = np.arange(0,100,dt)
sig = np.sin(2*np.pi*0.2*time)
norm_noise_sig = np.arange(0,5.0,1)

fig = pl.figure(figsize=(16,16))
t11 = fig.add_subplot(221)
t21 = fig.add_subplot(222)

cmap_choice = 'Dark2'

cmap1 = cm.get_cmap(cmap_choice,6)
norm_noise_sig[0] = 0.0001
for i,x in enumerate(norm_noise_sig):
    sig_noisy = sig + np.random.normal(0,x,len(sig))

    fft_freq = np.fft.fftfreq(len(sig_noisy),d=dt)
    fft = np.fft.fft(sig_noisy-np.mean(sig_noisy))

    se = spec_entropy(power=np.abs(fft)**2,freq=fft_freq,freq_range=[0.10,0.35])


    t11.plot(time/dt,sig_noisy,'.-',color=cmap1(i),label="se ="+str(np.round(se,2)),alpha=0.35,linewidth=2.0)
    fft_norm = np.abs(fft[:int(len(fft_freq)/2)])/np.sum(np.abs(fft[:int(len(fft_freq)/2)]))
    t21.plot(fft_freq[:int(len(fft_freq)/2)]*100,fft_norm,'.-',color=cmap1(i),label="se ="+str(np.round(se,2)),alpha=0.5)

t11.set_xlim(0,300)
t11.set_ylabel("Amplitude (au)",fontsize=15,fontweight='bold')
t21.set_xlim(0,40)
t21.set_xlabel("Frequency (Hz)",fontsize=15,fontweight='bold')
t21.set_ylabel("Normalized power",fontsize=15,fontweight='bold')


t11.legend(prop={'size':10,'weight':'bold'})
t21.legend(prop={'size':10,'weight':'bold'})
t11.set_title("Spectral entropy vs amount of noise",fontsize=15,fontweight='bold')
#fig.savefig("HS_noise.png")

#cmap1 = cm.get_cmap(cmap_choice,6)
# Spectral entropy vs second peak
sig2 = np.sin(2*np.pi*0.18*time) 
sig2_amp_ratio = np.arange(0,1.1,0.2)

#fig1 = pl.figure(figsize=(12,12))
t12 = fig.add_subplot(223)
t22 = fig.add_subplot(224)


for i,x in enumerate(sig2_amp_ratio):
    sig_noisy = sig2*x +(1-x)*sig 

    fft_freq = np.fft.fftfreq(len(sig_noisy),d=dt)
    fft = np.fft.fft(sig_noisy-np.mean(sig_noisy))

    se = spec_entropy(power=np.abs(fft)**2,freq=fft_freq,freq_range=[0.10,0.35])


    t12.plot(time/dt,sig_noisy,'.-',color=cmap1(i),label="se ="+str(np.round(se,2)),alpha=0.35,linewidth=2.0)
    fft_norm = np.abs(fft[:int(len(fft_freq)/2)])/np.sum(np.abs(fft[:int(len(fft_freq)/2)]))
    t22.plot(fft_freq[:int(len(fft_freq)/2)]*100,fft_norm,'.-',color=cmap1(i),label="se ="+str(np.round(se,2)),alpha=0.5)

t12.set_xlim(0,300)
t12.set_xlabel("Time (ms)",fontsize=15,fontweight='bold')
t12.set_ylabel("Amplitude (au)",fontsize=15,fontweight='bold')
t22.set_xlim(0,40)
t22.set_xlabel("Frequency (Hz)",fontsize=15,fontweight='bold')
t22.set_ylabel("Normalized power",fontsize=15,fontweight='bold')


t12.legend(prop={'size':10,'weight':'bold'})
t22.legend(prop={'size':10,'weight':'bold'})
t12.set_title("Spectral entropy vs multiple peaks",fontsize=15,fontweight='bold')
fig.savefig("HS_example.png")


pl.show()



