import numpy as np
import itertools
import shutil
import os
import pickle
import sys
import time
# Set the home directory
home_directory = ""

sys.path.append(home_directory+'/common/')
import misc2
import analyze_data_ssbn as adata
reload(adata)
import main_simulation as simd
reload(simd)

# If midway give sinoidal stimulation to stn neurons call import call_simulation_midway_sinstn 
import call_simulation as call_sim
#import call_simulation_midway_sinstn as call_sim

#storage_home = os.environ.get('STORAGE_HOME')
storage_home = home_directory+'/ssbn_effect/'

path1 = storage_home+'/scripts/' # Path to store the results
path2 = storage_home+'/Dist/' # Path to jdf files
path3 = storage_home+'/jdf/'
path5 = storage_home+'/output/' # Path for output 


poi_rate_bkg_stn = np.arange(1000.,2000.,200.)
poi_rate_bkg_gpe = np.arange(300.,1500.,200.)


# Refractory time can be given as input, default is 5.0ms
ref = sys.argv[1]
# Is this is the robustness analysis ? If "y", the connection probability, synaptic strength etc is not read from default params file but a separate file
perts = sys.argv[2]

# % of neurons that are made bursty
gpe_ratio_psdb = np.arange(0.,1.1,0.1)
stn_ratio_psdb = np.arange(0.,1.1,0.1)

for gpe_ip in poi_rate_bkg_gpe:
	for stn_ip in poi_rate_bkg_stn:
		for gpe_ratio in gpe_ratio_psdb:
			for stn_ratio in stn_ratio_psdb:
				sim_name = "psdbEffect_inh_ip_"+str(gpe_ip)+"_exc_ip_"+str(stn_ip)+"_inh_ratio_"+str(gpe_ratio)+"_exc_ratio_"+str(stn_ratio)+"ref_"+ref+"perts_"+perts
				#sim_name = "psdbEffect_inh_ip_no_beta_band_"+str(gpe_ip)+"_exc_ip_"+str(stn_ip)+"_inh_ratio_"+str(gpe_ratio)+"_exc_ratio_"+str(stn_ratio)
				params = dict()
				params['model_type'] = 'ssbn'
				params['add_spikes'] = 0
				params['T_sim'] = 9000.
				params["t_ref"] = float(ref)
				params["perts"] = perts
                                params["gpe_rat"] = gpe_ratio# % of bursting neurons in inh pop
                                params["stn_rat"] = stn_ratio# % of bursting neurons in exc pop

				params['pg_rate_exc'] = stn_ip
				params['pg_rate_inh'] = gpe_ip
				params['num_spikes_burst_stn'] = 4.
				params['num_spikes_burst_gpe'] = 4.
				params['add_spikes']= 0
				params['rnd_seeds'] = 1
				# The whole simulation with same % of neurons as bursty. If 1, neurons are changed after some time in simulation from non-bursty to bursty 

                                params['change_type'] = 1       # Changes the neuron midway
                                params['num_chg_stn'] = int(pars.Parameters.N[0]*pars.Parameters.order*stn_ratio)
                                params['chg_time_stn'] = 3500. # Change STN into bursty after 3500ms
                                params['num_chg_gpe'] = int(pars.Parameters.N[1]*pars.Parameters.order*gpe_ratio)
                                params['chg_time_gpe'] = 1500. # Change GPe into bursty after 1500ms

				
				params['record_vm'] =[]	# For now switch off the multimeter
		

				call_sim(params,sim_name)	

				
