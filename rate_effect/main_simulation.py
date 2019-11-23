import numpy as np
import nest
import scipy
import scipy.io as sio
import os.path
import sys
import os

# Set the home directory here
home_directory = ""

# For the common param files
sys.path.append(home_directory+"../common/")

import params_d_ssbn as params_d

neuron_param = {'V_th':-54.0, 'V_reset': -70.0, 't_ref': 5.0, 'g_L':10.0,'C_m':200.0, 'E_ex': 0.0, 'E_in': -80.0, 'tau_syn_ex':5.0,'tau_syn_in': 10.0,'tau_minus':20.0}


connect_stn_gpe = True
connect_poisson_bkg = True

poi_rate_bkg_gpe = np.arange(300.,1500.,200.)

pars = params_d.Parameters() # The param files object with all common parameters, like simtimes, number of neurons in GPe/STN, connectivity, delays etc

simtime = pars.T_total # 
Ngpe = pars.order*pars.N[1]
Nstn = pars.order*pars.N[0]
print "Ngpe",Ngpe
print "Nstn",Nstn

def runSim(params):
	stn_inp_rate = params["stn_inp"]
	seed = params["seed"]
	for ii in range(len(poi_rate_bkg_gpe)):

		nest.ResetKernel()
		path = os.getcwd()
		path = path+"/data/"
		if os.path.isdir(path+str(seed)+"/") == False:
			os.mkdir(path+str(seed)+"/")
		nest.SetStatus([0],{'data_path':'./data/'+str(seed),'overwrite_files': True,'grng_seed':seed})

		gp_neurons = nest.Create('ssbn',Ngpe,neuron_param)
		st_neurons = nest.Create('ssbn',Nstn,neuron_param)

		# Spike detectors
		gp_sd = nest.Create("spike_detector", 1)
		st_sd = nest.Create("spike_detector", 1)
		
		f_name_gp = 'GP_gp_' + str(poi_rate_bkg_gpe[ii]) + '_stn_' + str(stn_inp_rate)
		f_name_st = 'ST_gp_' + str(poi_rate_bkg_gpe[ii]) + '_stn_' + str(stn_inp_rate)
		
		
		nest.SetStatus(gp_sd,[{"label":f_name_gp,"withtime": True,"withgid": True,'to_memory':False,'to_file':True}])
		nest.SetStatus(st_sd,[{"label": f_name_st,"withtime": True,"withgid": True,'to_memory':False,'to_file':True}])

		if connect_poisson_bkg:
			print 'Connecting Poisson'
			#PG to GPE
			pg_gen_gpe = nest.Create('poisson_generator',1,{'rate':poi_rate_bkg_gpe[ii]})
			weights = np.random.uniform(low=0.5,high=1.5,size=Ngpe)
			delays = np.ones(Ngpe)
			nest.Connect(pg_gen_gpe,gp_neurons,syn_spec={"weight":[ [x] for x in weights],"delay":[[x] for x in delays]})
			# PG TO STN
			pg_gen_stn = nest.Create('poisson_generator',1,{'rate':stn_inp_rate})
			weights = np.random.uniform(low=0.5,high=1.5,size=Nstn)
			delays = np.ones(Nstn)
			nest.Connect(pg_gen_stn,st_neurons,syn_spec={'weight':[ [x] for x in weights],'delay':[[x] for x in delays]})


		if connect_stn_gpe:
			print 'STN GPE Connect'
			# random connectivity, synapse numbers
			syn_stn_gpe = {'rule':'fixed_outdegree','outdegree':int(Ngpe*pars.epsilon_stn_gpe)}
			syn_gpe_gpe = {'rule':'fixed_outdegree','outdegree':int(Ngpe*pars.epsilon_gpe_gpe)}
			syn_gpe_stn = {'rule':'fixed_outdegree','outdegree':int(Nstn*pars.epsilon_gpe_stn)}
			
			print syn_stn_gpe
			print syn_gpe_gpe
			print syn_gpe_stn
			# STN-STN == No connections
			print 'Connect STN-GPE'
			nest.Connect(st_neurons,gp_neurons,conn_spec=syn_stn_gpe,syn_spec={'weight':pars.J_stn_gpe,'delay':pars.del_stn_gpe})
			print 'Connect GPE-STN'                                                                                 
			nest.Connect(gp_neurons,st_neurons,conn_spec=syn_gpe_stn,syn_spec={'weight':pars.J_gpe_stn,'delay':pars.del_gpe_stn})
			print 'Connect GPE-GPE'                                                                                 
			nest.Connect(gp_neurons,gp_neurons,conn_spec=syn_gpe_gpe,syn_spec={'weight':pars.J_gpe_gpe,'delay':pars.del_gpe_gpe})

		# record spikes
		nest.Connect(gp_neurons,gp_sd)
		nest.Connect(st_neurons,st_sd)

		# simulate
		print 'Simulating '
		nest.Simulate(simtime)
		print 'Simulation finished'


