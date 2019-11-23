
#------------------------------------------------------------------------------------------
# This code reproduces Fig 2 and Supp Fig S4 in the manuscript
#------------------------------------------------------------------------------------------
#
#

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
import os
import pickle
import matplotlib.pyplot as plt
num_threads = 8


run_everything = sys.argv[1] #"y" - It will calculate the avg firing rate and spectral entropies from the raw data files and saves it in a dictionary called Data_<ref>.pickle. Should be ran at least once with "y". After that can be run with "n" where it reads the Data_<ref>.pickle and plots the figure

ref = sys.argv[2] # refractory period. This means there should be raw data files from simulations for the particular refractory period.

pars = params_d.Parameters()
def get_sim_3d(match_pars,seed):
	res = adata.analyze_data(prefix1,data_path=path1,netPars=pars,seed=seed)
	return res


def concatFiles(path,prefix):
	fhand = file(path+prefix+'.gdf', 'a+')
	for x in xrange(num_threads):
	   filename = path+"/"+prefix+"-{}.gdf".format(x)
	   if os.stat(filename).st_size == 0:
	       continue
	   temp = np.loadtxt(filename)
	   np.savetxt(fhand, temp)
	fhand.close()


poi_rate_bkg_stn = np.arange(1000.,1800.,100.)
poi_rate_bkg_gpe = np.arange(300.,1500.,100.)
exc_range = [0] 
inh_range = [0]
ext_rate = poi_rate_bkg_stn
add_inh_rate = poi_rate_bkg_gpe
# Set the seed path
seed_path = ""
seeds = pickle.load(open(seed_path+"/seeds.pickle","r"))
path1 = os.getcwd()+'/data/'
if run_everything == 'y':
	p1 = []
	lab_names = []
	legend_dict = dict()

	sim_id=1	
	big_gpe_ff_lis_mult = []
	big_gpe_fr_lis_mult = []
	big_stn_ff_lis_mult = []
	big_stn_fr_lis_mult = []
	big_stn_freq_lis_mult = []
	big_stn_spec_lis_mult = []
	big_gpe_freq_lis_mult = []
	big_gpe_spec_lis_mult = []

	for seed in seeds:
		if ref == "":
			path2 = path1+str(seed)+"/"	
		else:
			path2 = path1+str(seed)+"_ref_"+str(ref)+"/"	
		big_gpe_ff_lis = []
		big_gpe_fr_lis = []
		big_stn_ff_lis = []
		big_stn_fr_lis = []
		big_stn_freq_lis = []
		big_stn_spec_lis = []
		big_gpe_freq_lis = []
		big_gpe_spec_lis = []

		for exc_count, exc_val in enumerate(ext_rate):
			stn_ff_lis = []
			stn_fr_lis = []
			gpe_ff_lis = []
			gpe_fr_lis = []
			gpe_spec_lis = []
			gpe_freq_lis = []
			stn_spec_lis = []
			stn_freq_lis = []
			for inh_count, inh_val in enumerate(add_inh_rate):
				prefix1 = "GP_gp_"+str(inh_val)+"_stn_"+str(exc_val)
				prefix2 = "ST_gp_"+str(inh_val)+"_stn_"+str(exc_val)	
				GPact = np.loadtxt(path2+"GP_gp_"+str(inh_val)+"_stn_"+str(exc_val)+"-3001-0.gdf")
				STNact = np.loadtxt(path2+"ST_gp_"+str(inh_val)+"_stn_"+str(exc_val)+"-3002-0.gdf")
 
				numGP = pars.NGPe
				numSTN = pars.NSTN
				ad1 = get_sim_3d({'pg_rate_exc':exc_val, 'pg_rate_inh': inh_val, 'exc2_ratio':0., 'inh2_ratio':0.},seed)
				stn_fr_val = anal_wo.comp_mean_rate(time_range = [0.0,ad1.pars.T_wup+ ad1.pars.T_sim],act=STNact,total_time=ad1.pars.T_total,numNeurons=numSTN)			
				gpe_fr_val = anal_wo.comp_mean_rate(time_range = [0.0,ad1.pars.T_wup+ ad1.pars.T_sim],act=GPact,total_time=ad1.pars.T_total,numNeurons=numGP)

				xx,xx,xx,peak_freq_gpe, peak_pow_gpe = anal_wo.psd(time_range = [ad1.pars.T_wup,ad1.pars.T_wup+ ad1.pars.T_sim],act=GPact,total_time=ad1.pars.T_total)
				gpe_freq_val = peak_freq_gpe
				gpe_spec_val = anal_wo.spec_entropy(time_range = [ad1.pars.T_wup,ad1.pars.T_wup+ ad1.pars.T_sim], freq_range = [10.,35.],act=GPact,total_time=ad1.pars.T_total)

				xx,xx,xx,peak_freq_stn, peak_pow_stn = anal_wo.psd(time_range = [ad1.pars.T_wup,ad1.pars.T_wup+ ad1.pars.T_sim],act=STNact,total_time=ad1.pars.T_total)
				stn_freq_val = peak_freq_stn
				stn_spec_val = anal_wo.spec_entropy(time_range = [ad1.pars.T_wup,ad1.pars.T_wup+ ad1.pars.T_sim], freq_range = [10.,35.],act=STNact,total_time=ad1.pars.T_total)
				
				stn_fr_lis.append(stn_fr_val)
				gpe_fr_lis.append(gpe_fr_val)
				gpe_spec_lis.append(gpe_spec_val)
				gpe_freq_lis.append(gpe_freq_val)
				stn_spec_lis.append(stn_spec_val)
				stn_freq_lis.append(stn_freq_val)
				print exc_count,inh_count
			big_gpe_fr_lis.append(gpe_fr_lis)
			big_stn_fr_lis.append(stn_fr_lis)
			big_stn_spec_lis.append(stn_spec_lis)
			big_gpe_spec_lis.append(gpe_spec_lis)

		big_gpe_fr_lis_mult.append(big_gpe_fr_lis)
		big_stn_fr_lis_mult.append(big_stn_fr_lis)
		big_stn_spec_lis_mult.append(big_stn_spec_lis)
		big_gpe_spec_lis_mult.append(big_gpe_spec_lis)

		
	data = dict()
	data["bck_rate_stn"] = poi_rate_bkg_stn
	data["bck_rate_gpe"] = poi_rate_bkg_gpe
	data["gpe_rate"] = big_gpe_fr_lis_mult
	data["stn_rate"] = big_stn_fr_lis_mult
	data["gpe_spec"] = big_gpe_spec_lis_mult
	data["stn_spec"] = big_stn_spec_lis_mult
	if ref == "":
		pickle.dump(data,open(path1+"Data.pickle","w"))
	else:
		pickle.dump(data,open(path1+"Data_"+ref+".pickle","w"))

else:
	if ref == "":
		data = pickle.load(open(path1+"Data.pickle","r"))
	else:
		data = pickle.load(open(path1+"Data_"+ref+".pickle","r"))
	poi_rate_bkg_stn=data["bck_rate_stn"]  
	poi_rate_bkg_gpe=data["bck_rate_gpe"]  
	big_gpe_fr_lis_mult  =data["gpe_rate"]  
	big_stn_fr_lis_mult  =data["stn_rate"]  
	big_gpe_spec_lis_mult=data["gpe_spec"]  
	big_stn_spec_lis_mult=data["stn_spec"]  


add_inh_rate = np.append(add_inh_rate,1500.)
ext_rate = np.append(ext_rate,1800.)


fontDict = {'size':14, 'weight':'bold','color':'red'}

fig = pl.figure(1,figsize=(9.5,13.5))
ax1 = fig.add_subplot(421)
im1 = pl.pcolor(add_inh_rate, ext_rate, np.mean(array(big_gpe_fr_lis_mult),0), cmap = cm.BuGn)#,vmin=30.,vmax=45.)
ax1.set_xlim(add_inh_rate[0], add_inh_rate[-1])
ax1.set_ylim(ext_rate[0], ext_rate[-1])
ax1.set_xticks([])
ax1.set_yticks(ext_rate[::3])
ax1.set_ylabel('STN input',fontsize=15,fontweight='bold')
for x in ax1.get_yticklabels():
	x.set_fontsize(12)
	x.set_fontweight('bold')

ax1.set_title('GPe firing rate (spks/s)',fontsize=12,fontweight='bold')
ax1.text(150,1808,"A",fontsize=25,fontweight='bold')
cb1 = fig.colorbar(im1,ax=ax1)
for x in cb1.ax.yaxis.get_ticklabels():
	x.set_fontweight('bold')


ax1 = fig.add_subplot(422)
im2 = pl.pcolor(add_inh_rate, ext_rate, np.mean(array(big_stn_fr_lis_mult),0), cmap = cm.BuGn)#,vmin=30.,vmax=45.)
ax1.set_xlim(add_inh_rate[0], add_inh_rate[-1])
ax1.set_ylim(ext_rate[0], ext_rate[-1])
ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_title('STN firing rate (spks/s)',fontsize=12,fontweight='bold')
ax1.text(150,1808,"B",fontsize=25,fontweight='bold')
cb2 = fig.colorbar(im2,ax=ax1)
for x in cb2.ax.yaxis.get_ticklabels():
	x.set_fontweight('bold')


ax1 = fig.add_subplot(423)
im3 = pl.pcolor(add_inh_rate, ext_rate, np.nanmean(array(big_gpe_spec_lis_mult),0), cmap = cm.BuGn,vmin=0.3,vmax=0.7)
ax1.set_xlim(add_inh_rate[0], add_inh_rate[-1])
ax1.set_ylim(ext_rate[0], ext_rate[-1])
ax1.set_xticks(add_inh_rate[::3])
ax1.set_yticks(ext_rate[::3])
ax1.set_ylabel('STN input',fontsize=15,fontweight='bold')
ax1.set_xlabel('GPe input',fontsize=15,fontweight='bold')
ax1.text(150,1808,"C",fontsize=25,fontweight='bold')
for x in ax1.get_yticklabels():
	x.set_fontsize(12)
	x.set_fontweight('bold')
for x in ax1.get_xticklabels():
	x.set_fontsize(12)
	x.set_fontweight('bold')



ax1.set_title('GPe Spectral Entropy',fontsize=12,fontweight='bold')
x1 = 650.;x2=1150.;y1=1000.;y2=1800.
m = (y2-y1)/(x2-x1)
c = y1-m*x1
x_gen = np.arange(x1,x2+100,100)
y_gen = m*x_gen+c


ax1.plot(x_gen,y_gen,'r--',linewidth=2.5)
ax1.text(320,1450, "Oscillatory",fontdict=fontDict)
ax1.text(900,1150, "Non-",fontdict=fontDict)
ax1.text(850,1050, "Oscillatory",fontdict=fontDict)

cb3 = fig.colorbar(im3,ax=ax1)
for x in cb3.ax.yaxis.get_ticklabels():
	x.set_fontweight('bold')



ax1 = fig.add_subplot(424)
ind = np.where(np.isnan(np.array(big_stn_spec_lis_mult))==True)
np.array(big_stn_spec_lis_mult)[ind] = 1.0
im4 = pl.pcolor(add_inh_rate, ext_rate, np.nanmean(array(big_stn_spec_lis_mult),0), cmap = cm.BuGn,vmin=0.3,vmax=0.7)
ax1.set_xlim(add_inh_rate[0], add_inh_rate[-1])
ax1.set_ylim(ext_rate[0], ext_rate[-1])
ax1.set_xticks(add_inh_rate[::3])
ax1.set_yticks([])
ax1.set_xlabel('GPe input',fontsize=15,fontweight='bold')
ax1.set_title('STN Spectral Entropy',fontsize=12,fontweight='bold')
for x in ax1.get_xticklabels():
	x.set_fontsize(12)
	x.set_fontweight('bold')

ax1.text(320,1450, "Oscillatory",fontdict=fontDict)
ax1.text(900,1150, "Non-",fontdict=fontDict)
ax1.text(850,1050, "Oscillatory",fontdict=fontDict)

ax1.text(150,1808,"D",fontsize=25,fontweight='bold')

ax1.plot(x_gen,y_gen,'r--',linewidth=2.5)

cb4 = fig.colorbar(im4,ax=ax1)
for x in cb4.ax.yaxis.get_ticklabels():
	x.set_fontweight('bold')


ax2 = plt.subplot2grid((4,2),(2,0),colspan=2,rowspan=2)
im5 = ax2.pcolor(np.mean(array(big_gpe_fr_lis_mult),0),np.mean(array(big_stn_fr_lis_mult),0),np.nanmean(array(big_stn_spec_lis_mult),0),cmap=cm.BuGn,vmin=0.3,vmax=0.7,edgecolor='k')
ax2.set_xlabel("GPe firing rate (spks/s)",fontsize=15,fontweight='bold')
ax2.set_ylabel("STN firing rate (spks/s)",fontsize=15,fontweight='bold')
for x in ax2.get_xticklabels():
	x.set_fontsize(15)
	x.set_fontweight('bold')
for x in ax2.get_yticklabels():
	x.set_fontsize(15)
	x.set_fontweight('bold')

cb5 = fig.colorbar(im5,ax=ax2)
for x in cb5.ax.yaxis.get_ticklabels():
	x.set_fontweight('bold')

cb5.ax.yaxis.set_label_text("Spectral Entropy (SE)",fontsize=15,fontweight='bold')
ax2.text(23,42,"E",fontsize=25,fontweight='bold')
#rateVsSE = np.zeros(len(

plt.subplots_adjust(hspace=0.35,wspace=0.3,top=0.95,bottom=0.05)


pl.savefig("rate_effect_5trials_"+ref+".pdf")
pl.savefig("rate_effect_5trials_"+ref+".jpg")



#--------------------------------------------------- This plots the Supp Fig 4 -----------------------------------------------------------------



fig1 = pl.figure(2,figsize=(9,16))
t1 = fig1.add_subplot(211)
for i,x in enumerate(big_gpe_fr_lis_mult):
	t1.plot(x[0],big_gpe_spec_lis_mult[i][0],'.',color=cm.BuGn(float(i+1)/len(big_gpe_fr_lis_mult)),markersize=13,label="Trial "+str(i+1))
	t1.plot(x,big_gpe_spec_lis_mult[i],'.',color=cm.BuGn(float(i+1)/len(big_gpe_fr_lis_mult)),markersize=13)

t1.set_title("GPe",fontsize=25,fontweight='bold')
t1.set_ylabel("SE",fontsize=25,fontweight='bold')
t1.legend(prop={"size":10,"weight":"bold"})
for x in t1.get_yticklabels():
	x.set_fontsize(15)
	x.set_fontweight('bold')
for x in t1.get_xticklabels():
	x.set_fontsize(15)
	x.set_fontweight('bold')


t2 = fig1.add_subplot(212)
for i,x in enumerate(big_stn_fr_lis_mult):
	t2.plot(x[0],big_stn_spec_lis_mult[i][0],'.',color=cm.BuGn(float(i+1)/len(big_stn_fr_lis_mult)),markersize=13,label="Trial "+str(i+1))
	t2.plot(x,big_stn_spec_lis_mult[i],'.',color=cm.BuGn(float(i+1)/len(big_stn_fr_lis_mult)),markersize=13)

t2.set_title("STN",fontsize=25,fontweight='bold')
t2.set_ylabel("SE",fontsize=25,fontweight='bold')
t2.legend(prop={"size":10,"weight":"bold"})
t2.set_xlabel("Firing rates",fontsize=25,fontweight='bold')
for x in t2.get_yticklabels():
	x.set_fontsize(15)
	x.set_fontweight('bold')
for x in t2.get_xticklabels():
	x.set_fontsize(15)
	x.set_fontweight('bold')

plt.subplots_adjust(hspace=0.15,top=0.95,bottom=0.05)
fig1.savefig("RateVsSE_new_"+ref+".jpg")
fig1.savefig("RateVsSE_new_"+ref+".pdf")





