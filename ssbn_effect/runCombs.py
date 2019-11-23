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
import call_simulation as call_sim
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
				params['change_type'] = 0	# Changes the neuron midway
				
				params['record_vm'] =[]	# For now switch off the multimeter
		

				call_sim(params,sim_name)	

				'''
				# To run on SLURM cluster
				fh = open(path3 + '%s.py'%(sim_name),'w')
				fh.write('import sys\n')
				fh.write("sys.path.insert(0,'/home/j.bahuguna/PSDB/psdb_effect/')\n")
				fh.write('import call_simulation as single_sim\n')
				fh.write('single_sim.runSim(%s,%s)\n'% (params,repr(sim_name)))
				fh.close()
				
				print sim_name

				fh = open(path3 + '%s.jdf'%(sim_name),'w')
				content = [
				'#!/bin/bash\n',
				#   '#PBS -t 0-%d\n'%(len(comb_pars)-1),
				#'#SBATCH -o /home/j.bahuguna/homology/vAModel/output/test_job_.out',
				'#SBATCH --output=/home/j.bahuguna/PSDB/psdb_effect/output/PDF_%s.out'%(sim_name),
				'#SBATCH --error=/home/j.bahuguna/PSDB/psdb_effect/output/PDF_%s.err'%(sim_name),
				#'#PBS -j oe\n',
				'#SBATCH --job-name=%s\n'%(sim_name),
				'#SBATCH --mem-per-cpu=3500mb\n',
				'#SBATCH --time 3:00:00\n',
				#'#SBATCH -p long \n',
				#'#SBATCH --output=%s%s_%s.txt\n'%(path3,sim_name,postfix),
				#'export PYTHONPATH=/clustersw/cns/nest-2.2.2/lib/python2.7/dist-packages/:$PYTHONPATH\n',
				'python /home/j.bahuguna/PSDB/psdb_effect/jdf/%s.py\n'%(sim_name),
				#'python /bcfgrid/data/padmanabhan/scripts/levcomp/batch/%s_%s.py'%(sim_name,str(nr))
				]
				fh.writelines(content)
				fh.close()
				filename = path3 + '%s.jdf'%(sim_name)	
				os.system('sbatch  %s'%filename )
				'''		
				
