#To test the SSBN neuron model


import numpy as np
import nest
import pylab as pl

labelsize = 14.
ticksize = 14.
tsim = 500000.
# Custom chosen colors
col_lis = [np.array([63.,25.,255.])/255., np.array([204.,99.,20.])/255.,np.array([178.,120.,76.])/255., np.array([200.,0.,0.])/255.,np.array([153.,88.,61.])/255., np.array([204.,187.,20.])/255.]

# Range of input for the FI curve
amp_arr = np.arange(5000.,10000.,1000.)
sup_fr_lis =[]
fig = pl.figure(1,(8,5))
ax1 = fig.add_subplot(111)
# Range of burst lengths (1 to 5)
for kk,jj in enumerate(np.arange(5.)):
  fr_lis = []
  seed = np.random.randint(0,9999999,1)
  print seed
  for ii in amp_arr:
    nest.ResetKernel()

    nest.SetKernelStatus({'resolution':0.1,'grng_seed':seed[0]}) # Setting a new random seed each time
    # Poisson generator with rate as amp_arr
    dc_gen = nest.Create("poisson_generator",params = {'rate':ii})
    # Create ssbn neuron
    aa = nest.Create("ssbn",params = {'spb':jj+1})
    
    sd = nest.Create('spike_detector')
    # Connect ssbn to spike detector
    nest.Connect(aa,sd)
    # Poisson generator to ssbn 
    nest.Connect(dc_gen,aa)
    # Simulate
    nest.Simulate(tsim)
    # Read spikes
    spikes = nest.GetStatus(sd,'events')[0]
    num_spikes = len(spikes['senders'])
    print num_spikes
    # Calculate rate
    f_rate = (num_spikes*1000)/tsim
    fr_lis.append(f_rate)
  # Plot FI
  ax1.plot(amp_arr, fr_lis,lw = 5., alpha = 0.7, color = col_lis[kk], label = str(jj+1)) 
    
pl.legend(loc = 'best', prop = {'size':10.})
#pl.xlim(300,500)
ax1.set_ylim(0,40)
ax1.set_ylabel("Firing rate (Hz)",size = labelsize)
ax1.set_xticks(amp_arr[::2])
ax1.set_yticks(np.arange(0,40,10))
for tl in ax1.get_xticklabels():
  tl.set_fontsize(ticksize)
  ax1.xaxis.get_major_formatter().set_powerlimits((0,1))
ax1.set_xlabel("Poisson rate (Hz)",fontsize = labelsize)
ax1.text(6.2e4,40,'B',fontsize = labelsize + 3, style = 'normal')
for tl in ax1.get_yticklabels():
  tl.set_fontsize(ticksize)
  
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.tick_params(axis = 'both', which = 'both',top = 'off', right = 'off')
pl.savefig("FIs_psdb.png")
pl.show()  
    
