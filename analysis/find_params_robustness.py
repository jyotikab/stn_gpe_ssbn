
import matplotlib as mpl
mpl.use('Agg')

import pylab as pl
import numpy as np
import matplotlib.cm as cm
from pylab import *
import scipy.sparse.linalg as sp
import itertools
import marshal
import pickle
import logging
import gc
from scipy.integrate import odeint
import sys
import glob
from operator import attrgetter

import copy
from deap import base, creator
from deap import tools
import deap
import random
import numpy as np
import nest
import scipy
import scipy.io as sio
import os.path
import sys
import main_simulation_robustness_sampling as stn_gpe_samp
import os
from sklearn import svm
import pdb
import itertools as it

home_directory = ""
sys.path.append(home_directory+'/common/')
import params_d_ssbn as params_d

path = "../rate_effect/"
pars = params_d.Parameters()




def fitness_function(sols):

	print "sols",sols
		
	params = dict()
	params["J_gpe_gpe"] = sols[0] # J_gpe_gpe, mean =-0.725
	params["J_gpe_stn"] = sols[1] # J_gpe_stn , mean = -0.8
	params["J_stn_gpe"] = sols[2] # J_stn_gpe, mean = 1.2
	params["epsilon_gpe_gpe"] = sols[3] # epsilon_gpe_gpe , mean = 0.02
	params["epsilon_stn_gpe"] = sols[4] # epsilon_stn_gpe, mean = 0.023
	params["epsilon_gpe_stn"] = sols[5] # epsilon_gpe_stn, mean = 0.035
	params["del_gpe_gpe"] = sols[6]
	params["del_gpe_stn"] = sols[7]
	params["del_stn_gpe"] = sols[8]

	
	err = 0

	stn_fr_mat, gpe_fr_mat, stn_se_mat =  stn_gpe_samp.runSim(samp=True,new_params=params)

	print stn_se_mat
	ind_osc = np.where(stn_se_mat <= 0.4)
	ind_non_osc = np.where(stn_se_mat >= 0.55)
	XX,YY = np.meshgrid(poi_rate_bkg_gpe,poi_rate_bkg_stn)

	xy =  np.vstack([XX.ravel(), YY.ravel()]).T
	
	if len(ind_osc[0]) == 0 or len(ind_non_osc[0]) == 0: #All oscillatory or all non-oscillatory 
		return 0
	labels = np.zeros(np.shape(stn_se_mat))
	labels[ind_osc] = 1
	labels = (labels.T).flatten()

	ind_osc = np.where(labels ==1)
	ind_non_osc = np.where(labels==0)


	clf = svm.SVC(kernel='linear',C=1.0)
	clf.fit(xy,labels)
	score = clf.score(xy,labels)

	w = clf.coef_[0]
	m = -w[0]/w[1]  # Shoudl be positive
	c = clf.intercept_[0]/w[1]
	y_dec = m*np.unique(xy[:,0])-c


	
	print "score ",score

	if score > 0.99 and m > 0: # linearly separable and slope > 0
		print "found"
		ans["sols"].append(sols)
		ans["score"].append(score)	
		ans["stn_fr_mat"].append(stn_fr_mat)
		ans["gpe_fr_mat"].append(gpe_fr_mat)
		ans["stn_se_mat"].append(stn_se_mat)
		print ans
	
		pickle.dump(ans,open(path+"Combinations_normal_dist_new_"+str(seed1)+".pickle","w"))
	return score,


def sample_normal_dist(mu=0,var=1,samp_num=10):
	val = np.random.normal(mu,var,samp_num)
	if mu < 0:
		val = val[np.where(val<0)]
	else:
		val = val[np.where(val>0)]
	return val	


def paramSearch(params,anti=0):

	seed = params["seed"]
	print "in Paramsearch"

	logging.basicConfig(level=logging.DEBUG)
	leak = -0.1

	global seed1
	seed1 = seed
	np.random.seed(seed1)

	
	global ans
	ans = dict()
	ans["sols"] = []
	ans["score"] = []
	ans["stn_fr_mat"] = []
	ans["gpe_fr_mat"] = []
	ans["stn_se_mat"] = []


	num_samp = 5
	j_gpe_gpe_dist = sample_normal_dist(mu = pars.J_gpe_gpe,var = np.abs(pars.J_gpe_gpe*0.1),samp_num=num_samp)
	j_gpe_stn_dist = sample_normal_dist(mu = pars.J_gpe_stn,var = np.abs(pars.J_gpe_stn*0.1),samp_num=num_samp)
	j_stn_gpe_dist = sample_normal_dist(mu = pars.J_stn_gpe,var = np.abs(pars.J_stn_gpe*0.1),samp_num=num_samp)
	
	eps_gpe_gpe = sample_normal_dist(mu = pars.epsilon_gpe_gpe,var = np.abs(pars.epsilon_gpe_gpe*0.1),samp_num=num_samp)
	eps_stn_gpe = sample_normal_dist(mu = pars.epsilon_stn_gpe,var = np.abs(pars.epsilon_stn_gpe*0.1),samp_num=num_samp)
	eps_gpe_stn = sample_normal_dist(mu = pars.epsilon_gpe_stn,var = np.abs(pars.epsilon_gpe_stn*0.1),samp_num=num_samp)

	delay_gpe_gpe = sample_normal_dist(mu = pars.del_gpe_gpe,var = np.abs(pars.del_gpe_gpe*0.1),samp_num=num_samp)
	delay_gpe_stn = sample_normal_dist(mu = pars.del_gpe_stn,var = np.abs(pars.del_gpe_stn*0.1),samp_num=num_samp)
	delay_stn_gpe = sample_normal_dist(mu = pars.del_stn_gpe,var = np.abs(pars.del_stn_gpe*0.1),samp_num=num_samp)

		

	combs = list(it.product(j_gpe_gpe_dist,j_gpe_stn_dist,j_stn_gpe_dist,eps_gpe_gpe,eps_stn_gpe,eps_gpe_stn,delay_gpe_gpe,delay_gpe_stn,delay_stn_gpe))
	
	np.random.shuffle(combs)
			
	for comb in combs:
		fitness_function(comb)










