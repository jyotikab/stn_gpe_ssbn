
#-----------------------------------------------------------------------
# Reproduces Supp Fig S2
# Should be run after "find_params_robustness.py"
# If run_everything = "y", then also plots Fig 3A for all alternate solutions 
#----------------------------------------------------------------------


import numpy as np
import itertools
import sys
home_directory = ""
sys.path.append(home_directory+'/common/')
import analyze_data_ssbn as adata
import params_d_ssbn as params_d
import analysis_funcs_wo_framework as anal_wo
reload(adata)
import sys
from numpy import *
reload(adata)
import pdb
import pylab as pl
import matplotlib.cm as cm
import plotly.plotly as py
from plotly.graph_objs import *
import os
import pickle
import matplotlib.pyplot as plt
import matplotlib
from sklearn import svm
import os
import glob 
import seaborn as sns
import pandas as pd

poi_rate_bkg_stn = np.arange(1000.,2000.,200.)
poi_rate_bkg_gpe = np.arange(300.,1500.,200.)


#
run_everything = sys.argv[1]
params = ["j_gpe_gpe","j_gpe_stn","j_stn_gpe","eps_gpe_gpe","eps_stn_gpe","eps_gpe_stn","del_gpe_gpe","del_gpe_stn","del_stn_gpe"]

sols =[]
sols_final =[]
stn_frs =[]
stn_frs_final =[]
gpe_frs = []
gpe_frs_final = []
stn_se = []
stn_se_final = []

#all_files = glob.glob("output/old/C*[0-9]*.pickle")
all_files = glob.glob("../rate_effect/C*[0-9]*.pickle")
for fi in all_files:
	comb = pickle.load(open(fi,"r"))
	for x in comb["sols"]:
		sols.append(x)
	for x in comb["stn_fr_mat"]:
		stn_frs.append(x)
	for x in comb["gpe_fr_mat"]:
		gpe_frs.append(x)
	for x in comb["stn_se_mat"]:
		stn_se.append(x)


XX,YY = np.meshgrid(poi_rate_bkg_gpe,poi_rate_bkg_stn)
xy =  np.vstack([XX.ravel(), YY.ravel()]).T

if run_everything == "y":

	for i,x in enumerate(stn_se):
		labels = np.zeros(np.shape(x))
		ind_osc = np.where(x <= 0.4)
		ind_non_osc = np.where(x >= 0.55)
		# Check that there is at least one oscillatory and non-oscillatory combination
		if len(ind_osc[0]) ==0 or len(ind_non_osc[0])==0:
			continue

		labels[ind_osc] = 1
		labels = (labels.T).flatten()
		print "labels",labels	
		clf = svm.SVC(kernel='linear',C=1.0)
		clf.fit(xy,labels)
		score = clf.score(xy,labels)

		w = clf.coef_[0]
		m = -w[0]/w[1]

		c = clf.intercept_[0]/w[1]
		'''	
		print "score",score
		print "m",m
		print "w",w
		print "c",c
		'''
		sol = sols[i]
		tit = ''
		for j,s in enumerate(sol):
			tit = tit+params[j]+":"+str(np.round(s,2))+","

				
		if score > 0.95 and m > 0:
			print "score",score
			print "m",m
			print "w",w
			print "c",c
			
			sols_final.append(sols[i])
			stn_frs_final.append(stn_frs[i])
			gpe_frs_final.append(gpe_frs[i])
			stn_se_final.append(x)

			pl.figure(figsize=(16,16))
			pl.pcolor(np.append(poi_rate_bkg_gpe,1500),np.append(poi_rate_bkg_stn,2000.),x.T,cmap=cm.Greens,vmin=0.4,vmax=0.6)	
			y_dec = m*XX[0]-c

			pl.plot(XX[0][:3],y_dec[:3],'k--')
			pl.xlim(np.min(poi_rate_bkg_gpe),np.max(poi_rate_bkg_gpe))
			pl.ylim(np.min(poi_rate_bkg_stn),np.max(poi_rate_bkg_stn))

			pl.suptitle(tit,fontsize=10,fontweight='bold')
			pl.colorbar()
			pl.savefig("alternate_sols_"+str(i)+".png")
			#pl.show() 

	all_final = dict()
	all_final["sols"] = sols_final
	all_final["stn_fr_mat"] = stn_frs_final
	all_final["gpe_fr_mat"] = gpe_frs_final
	all_final["stn_se_mat"] = stn_se_final

	pickle.dump(all_final,open("../rate_effect/final_sols.pickle","w"))
	
else:
	sols_final = pickle.load(open("../rate_effect/final_sols.pickle","r"))["sols"]


def plot_violinplot(ax,data,color,ticklabels,title,orig):

	sur_data = [ np.random.normal(x,np.abs(x)*0.1,2000) for x in orig]

	min_data = [np.min(x) for x in sur_data]
	max_data = [np.max(x) for x in sur_data]

	print"min_data",min_data
	print"max_data",max_data


	df1 = pd.DataFrame()
	df1["Values"] = np.vstack((data[:,0],data[:,1],data[:,2])).flatten()
	df1["Variable name"] = np.hstack([[ticklabels[0]]*len(data[:,0]),[ticklabels[1]]*len(data[:,1]),[ticklabels[2]]*len(data[:,2])])
	df1["Distributions"] = "Result"


	df2 = pd.DataFrame()
	df2["Values"] = np.vstack((sur_data[0],sur_data[1],sur_data[2])).flatten()
	df2["Variable name"] =  np.hstack([[ticklabels[0]]*len(sur_data[0]),[ticklabels[1]]*len(sur_data[1]),[ticklabels[2]]*len(sur_data[2])])
	df2["Distributions"] = "Sampled"

	frames= [df1,df2]
	result = pd.concat(frames)		
	sns.set(font_scale=2)
	vio_plot = sns.violinplot(x="Variable name",y="Values",hue="Distributions",data=result,split=True,ax=ax,palette={"Sampled":"slategray","Result":"sienna"},legend=False)

	if "strengths" in title:
		ylabel = "nS"
	elif "delays" in title:
		ylabel = "ms"
	else:
		ylabel =""

	vio_plot.set_ylabel(ylabel)
	vio_plot.set_xlabel("")
	vio_plot.axes.legend(loc=2)	
	ax.set_title(title,fontsize=25,fontweight='bold')	

plt.rcParams["text.usetex"] =True
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

fig = pl.figure(figsize=(20,11))
t1 = fig.add_subplot(131)
t2 = fig.add_subplot(132)
t3 = fig.add_subplot(133)
plot_violinplot(t1,np.array(sols_final)[:,:3],"darkolivegreen",[r'$J_{gpe \rightarrow gpe}$', r'$J_{gpe \rightarrow stn}$',r'$J_{stn \rightarrow gpe}$'],"Synaptic strengths",[-0.725,-0.8,1.2])
plot_violinplot(t2,np.array(sols_final)[:,3:6],"darkolivegreen",[r'$\epsilon_{gpe \rightarrow gpe}$', r'$\epsilon_{stn \rightarrow gpe}$', r'$\epsilon_{gpe \rightarrow stn}$' ],"Connection probability",[0.02,0.023,0.035] )
plot_violinplot(t3,np.array(sols_final)[:,6:],"darkolivegreen",[r'$\tau_{gpe \rightarrow gpe}$', r'$\tau_{gpe \rightarrow stn}$', r'$\tau_{stn \rightarrow gpe}$' ], "Synaptic delays",[3.0,6.0,6.0])


fig.subplots_adjust(left=0.06,wspace=0.25,right=0.96)

fig.savefig("Param_dists_violin.png")
fig.savefig("Param_dists_violin.pdf")






         	
	                
	        	
                	
                	
                	
                	
                	
                        
                	
                	
                        
                	
                	
                
