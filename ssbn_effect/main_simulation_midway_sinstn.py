import sys
#sys.path.insert(0,'/home/bahuguna/Utils/nest/lib/python2.7/site-packages/')
import nest
import numpy as np
import pylab as pl
import misc2
import cPickle as cp
import time
import os
import shutil
from datetime import datetime
import imp
import tempfile
import glob
import random
import pdb

home_directory = ""

class Simulation(object):
	    
    def __init__(self,prefix,new_pars=[],pars_file=[]):
        '''Creates and simulates a network in NEST'''
	#Temporary way to run without the selected pars file
        pars_file = []
	print new_pars
        start_build_net = time.time()

        if pars_file==[]:	# import generic params_d file
            import sys
	    sys.path.append(home_directory+"/common/") 		
	    import params_d_ssbn
	    reload(params_d_ssbn)
	    pars = params_d_ssbn.Parameters(new_pars)		
	else: 			# import specific params_d file
	    fobj,pathname,description = imp.find_module(pars_file)
	    params_d_sp = imp.load_module(pars_file,fobj,pathname,description)
	    pars = params_d_sp.Parameters(new_pars)
	    print "else"

	print pars.__dict__
	print "pars.num_chg_gpe",pars.num_chg_gpe
	print "pars.num_chg_stn",pars.num_chg_stn
	self.T_sim = pars.T_sim + pars.T_wup + pars.T_cdown
	self.record_vm = pars.record_vm
	self.recorders = {}
	self.events = {'spikes':[],'vm':[]}	
	self.pars = pars
	self.pars.prefix = prefix	
	self.pars.seed = new_pars["seed"]
	seed = self.pars.seed

	self.seed = seed
	print "seed",self.seed
	# INITIALIZE NETWORK -----------------------------------------------------------------------	        
        nest_path_tmp = tempfile.mktemp(prefix=pars.nest_path_tmp)
        os.mkdir(nest_path_tmp)
        nest.ResetKernel()
	print seed[0]
	#nest.SetKernelStatus({"resolution":0.1})
        shutil.rmtree(nest.GetStatus([0],'data_path')[0],ignore_errors=True)
        nest.SetKernelStatus({'resolution': pars.dt, 'print_time': pars.print_time,
        'overwrite_files':pars.owr_files, 'grng_seed':seed[0],
        'data_path':nest_path_tmp})
        
        #print '\nBuilding network...'
        
        # CREATE SOURCES ----------------------------------------------------------------------------
	# For the sine wave
        self.pg_exc_part_sin = nest.Create('sinusoidal_poisson_generator', 1)
        self.pg_exc_part_norm = nest.Create('poisson_generator', 1)
        self.pg_exc = nest.Create('poisson_generator', 1)
        self.pg_inh = nest.Create('poisson_generator', 1)

        nest.SetStatus(self.pg_exc, {'rate': pars.pg_rate_exc,'stop':self.pars.chg_time_stn})
	# Sine wave of frequency 15 Hz	
        nest.SetStatus(self.pg_exc_part_sin, {'rate': pars.pg_rate_exc*0.5,'start':self.pars.chg_time_stn,'frequency':15.,'amplitude':pars.pg_rate_exc*0.5})
        nest.SetStatus(self.pg_exc_part_norm, {'rate': pars.pg_rate_exc,'start':self.pars.chg_time_stn})
        nest.SetStatus(self.pg_inh, {'rate': pars.pg_rate_inh})
                
        
         # CREATE POPULATIONS -----------------------------------------------------------------------
	#print 'Creating populations...\n'
	
	# STN	
	neurons_exc = []
	self.pops_exc = range(len(pars.N_exc))
	for ii,nr in enumerate(pars.N_exc):
	  self.pops_exc[ii] = nest.Create(pars.model_type, abs(nr))
	  neurons_exc.extend(self.pops_exc[ii])
	  
	#  set neuron parameters	for every population independently              
	for ntypes in range(len(pars.N_exc)):
	  nest.SetStatus(self.pops_exc[ntypes], pars.neuron_params_exc[ntypes])	
	  
	
	if pars.rnd_dist:
	  nest.SetStatus(neurons_inh,'tau_m',pars.tau_m_rnd)
	  
	neurons_inh = []
	self.pops_inh = range(len(pars.N_inh))
	for ii,nr in enumerate(pars.N_inh):
	  self.pops_inh[ii] = nest.Create(pars.model_type, abs(nr))
	  neurons_inh.extend(self.pops_inh[ii])
	  
	#  set neuron parameters	for every population independently              
	for ntypes in range(len(pars.N_inh)):
	  nest.SetStatus(self.pops_inh[ntypes], pars.neuron_params_inh[ntypes])
	  
	
	if pars.rnd_dist:	# False by default
	  nest.SetStatus(neurons_inh,'tau_m',pars.tau_m_rnd)
	  

	if pars.change_type:	# Changes neuron type midway, 0 by default
	  # Changing order, first Gpe, then STN 
	  self.time_lis = [0,pars.chg_time_gpe,pars.chg_time_stn]
	  
	  
	self.pops = self.pops_exc + self.pops_inh
	
	self.pops_exc = [item for sublist in self.pops_exc for item in sublist]
	self.pops_inh = [item for sublist in self.pops_inh for item in sublist]
	print "len(self.pops_exc)",len(self.pops_exc)
	print "len(self.pops_inh)",len(self.pops_inh)
 
	 # Make connections -------------------------------------------------------------------------
	self.pars.neurons_tot = len(self.pops_exc) + len(self.pops_inh) 

	self.pars.pops_exc = self.pops_exc
	self.pars.pops_inh = self.pops_inh
	

	weights_exc = np.random.uniform(low=0.5,high=1.5,size=len(self.pops_exc))
	delays_exc = np.ones(len(self.pops_exc))

	nest.Connect(self.pg_exc, self.pops_exc,syn_spec={'weight':[ [x] for x in weights_exc],'delay':[[x] for x in delays_exc]})
 
	weights_inh = np.random.uniform(low=0.5,high=1.5,size=len(self.pops_inh))
	delays_inh = np.ones(len(self.pops_inh))

	self.weights_exc = weights_exc
	self.weights_inh = weights_inh
	self.delays_inh = delays_inh
	self.delays_exc = delays_exc

	nest.Connect(self.pg_inh, self.pops_inh,syn_spec={'weight':[ [x] for x in weights_inh],'delay':[[x] for x in delays_inh]})
	print "STN params"
	print nest.GetStatus([self.pops_exc[-1]])

	print "GPe params"
	print nest.GetStatus([self.pops_inh[-1]])

         #STN connections
	num_stn_gpe = int(pars.epsilon_stn_gpe * len(self.pops_inh))
	nest.Connect(self.pops_exc,self.pops_inh,conn_spec={'rule':'fixed_outdegree','outdegree':num_stn_gpe,'autapses':False,'multapses':False}, syn_spec={'weight':pars.J_stn_gpe, 'delay': pars.del_stn_gpe}) 


	#GPE connections
	num_gpe_gpe = int(pars.epsilon_gpe_gpe * len(self.pops_inh))
	nest.Connect(self.pops_inh,self.pops_inh, conn_spec={'rule':'fixed_outdegree','outdegree':num_gpe_gpe,'autapses':False,'multapses':False}, syn_spec={'weight':pars.J_gpe_gpe, 'delay': pars.del_gpe_gpe})

	num_gpe_stn = int(pars.epsilon_gpe_stn* len(self.pops_exc))
	nest.Connect(self.pops_inh,self.pops_exc, conn_spec={'rule':'fixed_outdegree','outdegree':num_gpe_stn,'autapses':False,'multapses':False}, syn_spec={'weight':pars.J_gpe_stn, 'delay': pars.del_gpe_stn}) 

	self.record_spikes = self.pops_exc+self.pops_inh
		
	sd = nest.Create('spike_detector',1)
	nest.SetStatus(sd,{'to_file':True,'to_memory':True,'withgid':True})
	nest.Connect(self.record_spikes,sd)
	self.recorders['sd'] = sd
	
	
	if self.pars.record_vm != []:
	    vm = nest.Create('voltmeter',1)
	    #print 'Id of vm recorder: ',vm
	    nest.SetStatus(vm,{'withtime':True,'withgid':True,'to_file':True,'to_memory':False})
	    nest.Connect(vm,self.pars.record_vm)
	    nest.SetStatus(self.pars.record_vm,{'V_th':1000.}) # record free Vm
	    self.recorders['vm'] = vm
	    
	self.build_net_time = time.time()-start_build_net
	
    def sim(self):
	start_sim_time = time.time()
	if self.pars.change_type == 1:
	  for ii in self.time_lis:
	    if ii == self.pars.chg_time_gpe and ii < self.pars.chg_time_stn: # 
	      print 'changing type - gpe'
	      nest.SetStatus(self.pops[1][:int(self.pars.num_chg_gpe)], params = {'spb':self.pars.chg_spk_bs})
	      nest.Simulate(self.pars.chg_time_stn - self.pars.chg_time_gpe) 
	    elif ii == self.pars.chg_time_stn and ii < self.T_sim:	
	      print 'changing type - stn'
	      if int(self.pars.num_chg_stn) == 0:
              	nest.Connect(self.pg_exc_part_norm, self.pops[0][int(self.pars.num_chg_stn):],syn_spec={'weight':[ [x] for x in self.weights_exc[int(self.pars.num_chg_stn):]],'delay':[[x] for x in self.delays_exc[int(self.pars.num_chg_stn):]]})
	      elif int(self.pars.num_chg_stn) == 1000:
	      	nest.Connect(self.pg_exc_part_sin, self.pops[0][:int(self.pars.num_chg_stn)],syn_spec={'weight':[ [x] for x in self.weights_exc[:int(self.pars.num_chg_stn)]],'delay':[[x] for x in self.delays_exc[:int(self.pars.num_chg_stn)]]})
	      else:
	      	nest.Connect(self.pg_exc_part_sin, self.pops[0][:int(self.pars.num_chg_stn)],syn_spec={'weight':[ [x] for x in self.weights_exc[:int(self.pars.num_chg_stn)]],'delay':[[x] for x in self.delays_exc[:int(self.pars.num_chg_stn)]]})
		nest.Connect(self.pg_exc_part_norm, self.pops[0][int(self.pars.num_chg_stn):],syn_spec={'weight':[ [x] for x in self.weights_exc[int(self.pars.num_chg_stn):]],'delay':[[x] for x in self.delays_exc[int(self.pars.num_chg_stn):]]})
	      nest.Simulate(self.T_sim - self.pars.chg_time_stn) 
	    else:
	      nest.Simulate(self.pars.chg_time_gpe)
	else:
	  nest.Simulate(self.T_sim)	
	  print self.T_sim
	self.sim_real_time = time.time() - start_sim_time
	if not os.path.exists(self.pars.data_path):
	    os.makedirs(self.pars.data_path)
	    
	self.save_spikes()
	self.save_info()
	# important to remove any tmp folders after simulation
	tmp_dir = nest.GetStatus([0],'data_path')[0]	
	tmp_files = glob.glob(tmp_dir+'/*')
	for ff in tmp_files:	   
	    print not 'nfs' in ff
	    if not 'nfs' in ff:
		os.remove(ff)
		#print 'Removed file %s'%ff
	shutil.rmtree(tmp_dir,ignore_errors=True)
	#print 'Removed dir: ',tmp_dir
	    
    def get_spikes(self):
	  events = nest.GetStatus(self.recorders['sd'],'events')[0]
	  spikes = np.zeros((len(events['senders']),2))
	  spikes[:,0] = events['senders']
	  spikes[:,1] = events['times']
	  self.events['spikes'].append(spikes)
	  #self.events['spikes'] =spikes
	    
    def get_vm(self):
	  events = np.loadtxt(nest.GetStatus(self.recorders['vm'])[0]['filenames'][0])
	  potentials = np.zeros((len(events),3))
	  potentials[:,0] = events[:,0]
	  potentials[:,1] = events[:,1]
	  potentials[:,2] = events[:,2]
	  
	  self.events['vm'].append(potentials)
	    
    def save_spikes(self,data_path=[]):
      '''merge spikes from all populations into one vector
      and save it as a numpy array'''
      
      if data_path == []:
	  data_path = self.pars.data_path
      #if self.events['spikes']==[]:
      self.get_spikes()
      spikes = self.events['spikes']
      print "spikes",np.shape(spikes)	
      spikes_tot = np.array([0,0],ndmin=2)
      fname = data_path + self.pars.prefix +'_spikes.npy'
      for ii,sp_pop in enumerate(spikes[0:]):
	  if len(sp_pop)>0:
	      spikes_tot = np.concatenate((spikes_tot,sp_pop))
      spikes_tot = spikes_tot[1:]	# ignore first row used for concatenation
      np.save(fname,spikes_tot)
	
    def save_vm(self,data_path=[]):
	if data_path == []:
	    data_path = self.pars.data_path
	if self.events['vm']==[]:
	    vm = self.get_vm()
	vm = self.events['vm']
	vm_tot = vm[0]
	fname = data_path + self.pars.prefix + '_vm.npy'
	for ii,vm_pop in enumerate(vm[1:]):	    
	    vm_tot = np.concatenate((vm_tot,vm_pop))
	np.save(fname,vm_tot)
	
    def save_info(self,data_path=[]):
      if data_path == []:
	  data_path = self.pars.data_path	
      info1 = dict()
      pars0 = self.pars	
      #pars = dict()
      # convert to dict, because cPickle has problems with some objects
      pars = dict([(name,getattr(pars0,name))  for name in dir(pars0) if not name.startswith('_')])

      info1['pars'] = pars	
      tstamp = datetime.now()
      info1['tstamp'] = tstamp.strftime("%d-%m-%Y %H:%M:%S")
      info1['time_build'] = self.build_net_time
      info1['time_sim_real'] = self.sim_real_time

      pops = []
      for pop in self.pops:
	  cells = np.array(pop,dtype=int)
	  pops.append((cells[0],cells[-1]))
      info1['pops'] = pops
      #print '\n printing info:\n',info1,'\n'
      
      fname = data_path + self.pars.prefix + '.info'
      fh = open(fname,'w')
      cp.dump(info1,fh)	
      fh.close()
      #print 'Saved info-file in %s'%fname
      

	
