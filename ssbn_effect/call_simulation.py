import sys

# Set the home directory here
home_directory = ""

sys.path.append(home_directory+'common/')
import misc2
import analyze_data_ssbn as adata
reload(adata)
import main_simulation as simd
reload(simd)
import pylab as pl
import numpy as np
import cPickle as pickle
import pdb
numTrials = 5

# Set the path where the seeds are stored
seed_path = ""
seeds = pickle.load(open(seed_path+"seeds.pickle","r"))
def runSim(pars,prefix): 
	psth_exc_all=[]
	psth_inh_all=[]
	se_exc_all =[]
	freq_exc_all =[]
	se_inh_all= []
	freq_inh_all= []
	# For 5 trials
	for num in np.arange(0,numTrials):
		pars["seed"] = [seeds[num]]
		net1 = simd.Simulation(prefix,new_pars=pars)
		net1.sim()

		
		ad2 = adata.analyze_data(prefix,net1.pars,net1.seed)
		ad2.pars.pops_exc = net1.pars.pops_exc
		ad2.pars.pops_inh = net1.pars.pops_inh
		ad2.pars.T_sim = net1.pars.T_sim
		ad2.pars.t_ref = net1.pars.t_ref
		ad2.pars.perts = net1.pars.perts
		print "len(ad2.pars.pops_exc)",len(ad2.pars.pops_exc)
		print "len(ad2.pars.pops_inh)",len(ad2.pars.pops_inh)
		print net1.events['spikes']

		# Calculate PSTH of STN and GPe neurons
		psth_exc,xx_exc = ad2.comp_psth1(nmax=len(ad2.pars.pops_exc),pop_id=ad2.pars.pops_exc,color='r',binsize=5.)
		psth_inh,xx_inh = ad2.comp_psth1(nmax=len(ad2.pars.pops_inh),pop_id=ad2.pars.pops_inh,color='r',binsize=5.)

		# Calculate spectral entropy
		se_exc,power_exc,freq_exc,freq2_exc = ad2.spec_entropy(pop_id = ad2.pars.pops_exc,nmax=len(ad2.pars.pops_exc),freq_range=[10.,35.])
		se_inh,power_inh,freq_inh,freq2_inh = ad2.spec_entropy(pop_id = ad2.pars.pops_inh,nmax=len(ad2.pars.pops_inh),freq_range=[10.,35.])

		psth_exc_all.append(psth_exc)
		psth_inh_all.append(psth_inh)
		se_exc_all.append(se_exc)
		se_inh_all.append(se_inh)
		freq_exc_all.append(freq2_exc*100)
		freq_inh_all.append(freq2_inh*100)
		print "se_exc_all",se_exc_all
		print "se_inh_all",se_inh_all



	data=dict()
	data["psth_stn"] = [psth_exc_all,xx_exc]
	data["psth_gpe"] = [psth_inh_all,xx_inh]
	data["gids_stn"] = ad2.pars.pops_exc
	data["gids_gpe"] = ad2.pars.pops_inh
	data["freqSpec_stn"] = [se_exc_all,freq_exc_all]
	data["freqSpec_gpe"] = [se_inh_all,freq_inh_all]
	print "data",data
	pickle.dump(data,open("output/"+prefix+".pickle","w"),protocol=2)





