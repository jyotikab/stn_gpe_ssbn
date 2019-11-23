import numpy as np
import nest
import scipy
import scipy.io as sio
import os.path
import sys
import pickle
import os
home_directory = ""
sys.path.append(home_directory+'/common/')
import params_d_ssbn as params_d
import analysis_funcs_wo_framework as anal_wo
neuron_param = {'V_th':-54.0, 'V_reset': -70.0, 't_ref': 5.0, 'g_L':10.0,'C_m':200.0, 'E_ex': 0.0, 'E_in': -80.0, 'tau_syn_ex':5.0,'tau_syn_in': 10.0,'tau_minus':20.0,'spb':1.0}

connect_stn_gpe = True
connect_poisson_bkg = True


# lesser resolution for firing rates to STN and GPe
poi_rate_bkg_gpe = np.arange(300.,1500.,200.)
poi_rate_bkg_stn = np.arange(1000.,2000.,200.)

pars = params_d.Parameters() # Also get the connectivity, delays and synaptic strengths from here
pars.T_sim = 1000
simtime = pars.T_sim + pars.T_wup
Ngpe = pars.order*pars.N[1]
Nstn = pars.order*pars.N[0]
print "Ngpe",Ngpe
print "Nstn",Nstn

seed_path = ""
seeds = pickle.load(open(seed_path+"/seeds.pickle","r"))[:3]


def runSim(samp=False,new_params=[]):


	stn_fr_inps_all = np.zeros((len(seeds),len(poi_rate_bkg_gpe),len(poi_rate_bkg_stn)))
	gpe_fr_inps_all = np.zeros((len(seeds),len(poi_rate_bkg_gpe),len(poi_rate_bkg_stn)))
	stn_se_inps_all = np.zeros((len(seeds),len(poi_rate_bkg_gpe),len(poi_rate_bkg_stn)))

	for kk in range(len(seeds)):
		for jj in range(len(poi_rate_bkg_stn)):
			stn_inp_rate = poi_rate_bkg_stn[jj]
			for ii in range(len(poi_rate_bkg_gpe)):
			    #for jj in range(len(poi_rate_bkg_stn)):
				print jj
				print ii
				nest.ResetKernel()

				path = os.getcwd()
				path = path+"/../rate_effect/"
				nest.SetStatus([0],{'grng_seed':seeds[kk],'local_num_threads':4})
				gp_neurons = nest.Create('psdb',Ngpe,neuron_param)
				st_neurons = nest.Create('psdb',Nstn,neuron_param)

				# Spike detectors
				gp_sd = nest.Create("spike_detector", 1)
				st_sd = nest.Create("spike_detector", 1)
				
				f_name_gp = 'GP_gp_' + str(poi_rate_bkg_gpe[ii]) + '_stn_' + str(stn_inp_rate) + "_samp_"+str(samp)
				f_name_st = 'ST_gp_' + str(poi_rate_bkg_gpe[ii]) + '_stn_' + str(stn_inp_rate) + "_samp_"+str(samp)
				
				
				nest.SetStatus(gp_sd,[{"label":f_name_gp,"withtime": True,"withgid": True,'to_memory':True,'to_file':False}])
				nest.SetStatus(st_sd,[{"label": f_name_st,"withtime": True,"withgid": True,'to_memory':True,'to_file':False}])

				if connect_poisson_bkg:
					print 'Connecting Poisson'
					#PG to GPE
					pg_gen_gpe = nest.Create('poisson_generator',1,{'rate':poi_rate_bkg_gpe[ii]})
					weights = np.random.uniform(low=0.5,high=1.5,size=Ngpe)
					delays = np.ones(Ngpe)
					nest.Connect(pg_gen_gpe,gp_neurons,syn_spec={'weight':[ [x] for x in weights],'delay':[[x] for x in delays]})
					# PG TO STN
					pg_gen_stn = nest.Create('poisson_generator',1,{'rate':stn_inp_rate})
					weights = np.random.uniform(low=0.5,high=1.5,size=Nstn)
					delays = np.ones(Nstn)
					nest.Connect(pg_gen_stn,st_neurons,syn_spec={'weight':[ [x] for x in weights],'delay':[[x] for x in delays]})


				if connect_stn_gpe:
					print 'STN GPE Connect'
					if samp == True:
						pars.epsilon_stn_gpe = new_params["epsilon_stn_gpe"]
						pars.epsilon_gpe_gpe = new_params["epsilon_gpe_gpe"]
						pars.epsilon_gpe_stn = new_params["epsilon_gpe_stn"]
						pars.J_stn_gpe = new_params["J_stn_gpe"]
						pars.J_gpe_stn = new_params["J_gpe_stn"]
						pars.J_gpe_gpe = new_params["J_gpe_gpe"]
						pars.del_stn_gpe = new_params["del_stn_gpe"]
						pars.del_gpe_gpe = new_params["del_gpe_gpe"]
						pars.del_gpe_stn = new_params["del_gpe_stn"]


					print "pars.J_stn_gpe",pars.J_stn_gpe	
					print "pars.J_gpe_stn",pars.J_gpe_stn
					print "pars.epsilon_gpe_stn",pars.epsilon_gpe_stn
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

				gpeAct = np.vstack((nest.GetStatus(gp_sd)[0]['events']['senders'], nest.GetStatus(gp_sd)[0]['events']['times'])).T 
				stnAct = np.vstack((nest.GetStatus(st_sd)[0]['events']['senders'],nest.GetStatus(st_sd)[0]['events']['times'])).T


				
				# Calculate avg firing rates and spectral entropies
				stn_fr_inps_all[kk][ii][jj] = anal_wo.comp_mean_rate(time_range = [0.0,pars.T_wup+ pars.T_sim],act=stnAct,total_time=simtime,numNeurons=Nstn) 
				gpe_fr_inps_all[kk][ii][jj] = anal_wo.comp_mean_rate(time_range = [0.0,pars.T_wup+ pars.T_sim],act=gpeAct,total_time=simtime,numNeurons=Ngpe) 

				stn_se_inps_all[kk][ii][jj] = anal_wo.spec_entropy(time_range = [pars.T_wup,pars.T_wup+ pars.T_sim], freq_range = [10.,35.],act=gpeAct,total_time=pars.T_total)

	return np.nanmean(stn_fr_inps_all,axis=0), np.nanmean(gpe_fr_inps_all,axis=0),np.nanmean(stn_se_inps_all,axis=0)

			

