
import numpy as np
import sys
import pdb
import socket
import os


home_directory = ""
class Parameters(object):
	    	data_path = home_directory+'/ssbn_effect/output/' # Store results
    		nest_path_tmp = home_directory+'/ssbn_effect/tmp/' # nest throws up files
	

		print_time = False


		# GENERAL PARAMETERS    
		T_sim = 3000.		# (ms) main simulation time

		T_wup = 500.		# (ms) warmup time
		T_cdown = 0.		# (ms) cool-down time	# inputs are off in this last period
		T_total = T_sim + T_wup + T_cdown


		#NETWORK PARAMETERS   
		order = 1000    # Number os neurons
		N = [1,2]	# list with number of neurons per sub-population; negative indicates inhibitory population
		inh2_ratio = 0.	# ratio of 2nd inhibitory population, 0: no 2nd inh population # ratio of busrty neurons
		exc2_ratio = 0.      # ratio of 2nd excitatory population, 0: no 2nd population

		model_type = 'ssbn'
		min_del = 0.1
		max_del = 2.



		NSTN = N[0]*order
		NGPe = N[1]*order

		#Parameters for changing neuron types midway
		change_type = 0 # 1 changes type and 0 does not
		chg_time = 472. #Time at which the simulations need to be stopped and neuron types changed.s
		num_chg = 500   #Number of neurons for which parameters need to be changed
		chg_spk_bs = 4. #Number of spikes per burst for the changing neurons


		  

		# NEST SIMULATION PARAMETERS    
		owr_files = False
		dt = 0.1    	       	# ms simulation resolution (precise times are used, thus sort of irrelevant)

		#np.random.seed()		# different seed each time parameter object is created	

		rnd_seeds = 1
		#print rnd_seeds

		#record_spikes = 'all' 
		#record_spikes = list(np.arange(self.pops,5004))
		record_vm = [100]

		add_spikes = 0 # 0 if no additional inh. spikes to be added and one if additional spikes to be added.
		extra_spk = 'inh'
		C_m_exc = 200.
		C_m = 250.
		tau_syn_in = 10.0     # (ms) synaptic time constant
		tau_syn_ex = 5.0
		tau_m  = 20.0		# (ms) membrane time constant

		delay_inter = 5.                 # (ms) synaptic transmission delay

		delay_intra = 2.
				

		rnd_dist = False		# True if distributions are used for parameters


		dc1_pars = {			# external DC generator
		'start':400.,
		'stop':0.,
		'amplitude':0.,
		}

		dc2_pars = {			# external DC generator
		'start':600.,
		'stop':0.,
		'amplitude':0.,
		}


		pg_rate_exc = 11500.
		pg_rate_inh = 11500.
		extra_inh = 0.

		num_spikes_burst_gpe = 4.0
		num_spikes_burst_stn = 4.0
		J_ext = 1.0

		J_gpe_gpe =-0.725
		J_gpe_stn = -0.8	# Amin 
		J_stn_gpe = 1.2

		epsilon_gpe_gpe = 0.02
		epsilon_stn_gpe = 0.023
		epsilon_gpe_stn = 0.035

		del_stn_gpe = 6.
		del_gpe_gpe = 3.
		del_gpe_stn = 6.

		# Example parameter combination from robustness analysis
		'''
                J_gpe_gpe =-0.67
                J_gpe_stn = -1.0       
                J_stn_gpe = 1.04

                epsilon_gpe_gpe = 0.02
                epsilon_stn_gpe = 0.02
                epsilon_gpe_stn = 0.03

                del_stn_gpe = 5.96
                del_gpe_gpe = 3.14
                del_gpe_stn = 5.34
		'''



		def __init__(self,new_vals=[]):
			if new_vals!=[]: 	# overwrite pars
					for ii in new_vals.keys():	 
							if ii=='dc1_pars_amp':
								self.dc1_pars['amplitude'] = new_vals[ii]		
							elif isinstance(new_vals[ii],(float,int)):
								ss = 'self.%s = %f'%(ii,new_vals[ii])
							elif isinstance(new_vals[ii],str):
								ss = 'self.%s = "%s"'%(ii,new_vals[ii])
							exec(ss)	    

			self.N_exc = [(1-self.exc2_ratio)*self.N[0],self.exc2_ratio*self.N[0]]	
			self.K_exc = self.N_exc                             # for use in simulation_d to avoid non existent population
			# remove one population if it is zero size
			self.spiking_mode_exc = ['no_bursting','bursting']
			  

			if self.N_exc[0]==0:
				try:
					self.N_exc.pop(0)
					self.spiking_mode_exc.pop(0)
				except:
					print 'some error with removing zero population'		    
			elif len(self.N_exc)==2:
				if self.N_exc[1]==0:
						try:
							self.N_exc.pop(1)
							self.spiking_mode_exc.pop(1)
						except:
							print 'some error with removing zero population'
				
			self.N_inh = [(1-self.inh2_ratio)*self.N[1],self.inh2_ratio*self.N[1]]	
			self.K_inh = self.N_inh                             # for use in simulation_d to avoid non existent population
			self.spiking_mode_inh = ['no_bursting','bursting']
			  
				
			if self.N_inh[0]==0:
				try:
					self.N_inh.pop(0)
					self.spiking_mode_inh.pop(0)
				except:
					print 'some error with removing zero population'		    
			elif len(self.N_inh)==2:
				if self.N_inh[1]==0:
						try:
							self.N_inh.pop(1)
							self.spiking_mode_inh.pop(1)
						except:
							print 'some error with removing zero population'

			self.N_exc = np.array(self.N_exc) * int(self.order)
			self.N_exc = np.array(self.N_exc.round(),int)

			self.N_inh = np.array(self.N_inh) * int(self.order)
			self.N_inh = np.array(self.N_inh.round(),int)

			self.C_ext = self.epsilon*self.N[0]	# number of external connections per neuron


			self._set_neuron_pars()	#set single neuron parameters

			
			
		def _set_neuron_pars(self):
				self.neuron_params_exc = range(len(self.N_exc))

					
				if self.model_type == 'ssbn':
				  for ntypes in np.arange(0,len(self.N_exc)):
					if self.spiking_mode_exc[ntypes] == 'bursting':
					  self.neuron_params_exc[ntypes]={
							  'spb':self.num_spikes_burst_stn,
							  #'tau_m':self.tau_m,
							  'C_m':self.C_m_exc,
							  'tau_syn_in':self.tau_syn_in,
							  'tau_syn_ex':self.tau_syn_ex,
							'V_th':-54.0,
							'V_reset':-70.0,
							't_ref':5.0,
							'g_L':10.0,
							'E_ex':0.0,
							'E_in':-80.0
 
							}
							
					elif self.spiking_mode_exc[ntypes] == 'excitatory':
					  self.neuron_params_exc[ntypes]={
											 'spb':1.0,
							  'C_m':self.C_m_exc,
							  'tau_syn_in':self.tau_syn_in,
							  'tau_syn_ex':self.tau_syn_ex,
							'V_th':-54.0,
							'V_reset':-70.0,
							't_ref':5.0,
							'g_L':10.0,
							'E_ex':0.0,
							'E_in':-80.0
							}
							  


					elif self.spiking_mode_exc[ntypes] == 'no_bursting':
					  self.neuron_params_exc[ntypes]={
							'spb': 1.0,
							  'C_m':self.C_m_exc,
							  'tau_syn_in':self.tau_syn_in,
							  'tau_syn_ex':self.tau_syn_ex,
							'V_th':-54.0,
							'V_reset':-70.0,
							't_ref':5.0,
							'g_L':10.0,
							'E_ex':0.0,
							'E_in':-80.0
							}
							
				self.neuron_params_inh = range(len(self.N_inh))


					
				if self.model_type == 'ssbn':
				  for ntypes in np.arange(0,len(self.N_inh)):
					if self.spiking_mode_inh[ntypes] == 'bursting':
					  self.neuron_params_inh[ntypes]={
							  'spb':self.num_spikes_burst_gpe,
							  'C_m':self.C_m_exc,
							  'tau_syn_in':self.tau_syn_in,
							  'tau_syn_ex':self.tau_syn_ex,
							'V_th':-54.0,
							'V_reset':-70.0,
							't_ref':5.0,
							'g_L':10.0,
							'E_ex':0.0,
							'E_in':-80.0
							}
							
					elif self.spiking_mode_inh[ntypes] == 'excitatory':
					  self.neuron_params_inh[ntypes]={
											 'spb':1.0,
							  'C_m':self.C_m_exc,
							  'tau_syn_in':self.tau_syn_in,
							  'tau_syn_ex':self.tau_syn_ex,
							'V_th':-54.0,
							'V_reset':-70.0,
							't_ref':5.0,
							'g_L':10.0,
							'E_ex':0.0,
							'E_in':-80.0
							}
					elif self.spiking_mode_inh[ntypes] == 'no_bursting':
					  self.neuron_params_inh[ntypes]={
							'spb': 1.0,
							  'C_m':self.C_m_exc,
							  'tau_syn_in':self.tau_syn_in,
							  'tau_syn_ex':self.tau_syn_ex,
							'V_th':-54.0,
							'V_reset':-70.0,
							't_ref':5.0,
							'g_L':10.0,
							'E_ex':0.0,
							'E_in':-80.0
							}


