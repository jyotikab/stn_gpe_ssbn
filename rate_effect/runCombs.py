import numpy as np
import itertools
import shutil
import os
import pickle

import main_simulation as main_sim


poi_rate_bkg_stn = np.arange(1000.,2000.,200.)

#seeds = np.random.randint(0,9999999,5) # Either generate the seeds everytime, or generate them once, save them and resuse for them for reproducibility

#pickle.dump(seeds,open(path5+"seeds.pickle","w"))
# Set the seed path
seed_path = ""
seeds = pickle.load(open(seed_path+"seeds.pickle","r"))

for seed in seeds:
	for ip in poi_rate_bkg_stn:


		sim_name = "rateEffect"
		params = dict()
		params["stn_inp"] = ip
		params["seed"] = seed
		sim_name = sim_name+"_"+str(ip)+"_"+str(seed)
		main_sim.runSim(params)



