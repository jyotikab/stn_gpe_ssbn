
#--------------------------------------------------------------------------------------------
# Reproduces Figure 5  
#--------------------------------------------------------------------------------------------



import numpy as np
import itertools
home_directory = ""
sys.path.append(home_directory+'/common/')
import analyze_data_psd as adata
import params_d_psdb as params_d
import analysis_funcs_wo_framework as anal_wo
import simulation_nest_psdb as simd
#import common.ff_analyze_data as adata
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
num_threads = 8
import matplotlib.patches as patches
from mpl_toolkits.axes_grid.inset_locator import inset_axes


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

gpe_ratio_psdb = np.arange(0.,1.1,0.1)
stn_ratio_psdb = np.arange(0.,1.1,0.1)

exc_range = [0] 
inh_range = [0]
ext_rate = poi_rate_bkg_stn
add_inh_rate = poi_rate_bkg_gpe
# Set the seed path
seed_path = ""
seeds = pickle.load(open(seed_path+"/seeds.pickle","r"))
# Firing rate effect, for the light background
path1 = os.getcwd()+'/rate_effect/'
data = pickle.load(open(path1+"Data.pickle","r"))
poi_rate_bkg_stn=data["bck_rate_stn"]  
poi_rate_bkg_gpe=data["bck_rate_gpe"]  
big_gpe_fr_lis_mult  =data["gpe_rate"]  
big_stn_fr_lis_mult  =data["stn_rate"]  
big_gpe_spec_lis_mult=data["gpe_spec"]  
big_stn_spec_lis_mult=data["stn_spec"]  

fontDict = {'size':20, 'weight':'bold','color':'red'}

fig = pl.figure(1,figsize=(16,16))

ax2 = fig.add_subplot(111)
im5 = ax2.pcolor(np.mean(array(big_gpe_fr_lis_mult),0),np.mean(array(big_stn_fr_lis_mult),0),np.nanmean(array(big_stn_spec_lis_mult),0),cmap=cm.BuGn,vmin=0.3,vmax=0.7,edgecolor='k',alpha=0.4,linewidth=0.2)
ax2.set_xlabel("GPe firing rate (spks/s)",fontsize=20,fontweight='bold')
ax2.set_ylabel("STN firing rate (spks/s)",fontsize=20,fontweight='bold')
for x in ax2.get_xticklabels():
	x.set_fontsize(15)
	x.set_fontweight('bold')
for x in ax2.get_yticklabels():
	x.set_fontsize(15)
	x.set_fontweight('bold')


def plot_bursting_regime_wise(gpe_rate,stn_rate,subplot_ax,gpe_rate_list,stn_rate_list,se_mat,regime_txt,txt_offset_x,txt_offset_y):
	ind = np.where(np.logical_and(poi_rate_bkg_gpe == gpe_rate,poi_rate_bkg_stn == stn_rate)==True)[0]

	im3 = subplot_ax.pcolor(np.array(gpe_rate_list)[ind,ind][0],np.array(stn_rate_list)[ind,ind][0],np.array(se_mat)[ind,ind][0],cmap=cm.BuGn,vmin=0.3,vmax=0.7,edgecolor='k',alpha=0.7)

	subplot_ax.text(np.mean(np.array(gpe_rate_list)[ind,ind][0])+txt_offset_x,np.mean(np.array(stn_rate_list)[ind,ind][0])+txt_offset_y,regime_txt,color='r',fontweight='bold',fontsize=15)

	return im3



GPeSTNpairs = np.array([[500.0,1000.0 ],[900.0,1400.0],[1300.0,1400.0],[700.0,1000.0],[500.0,1400.0],[1100.0,1200.0]])	# 0th is non-monotonic, 1st is all oscillatory, 2nd is non-oscillatory regimes
poi_rate_bkg_gpe = GPeSTNpairs[:,0]
poi_rate_bkg_stn = GPeSTNpairs[:,1]
big_gpe_fr_lis = np.zeros((len(poi_rate_bkg_gpe),len(poi_rate_bkg_stn),len(gpe_ratio_psdb),len(stn_ratio_psdb)))

big_stn_fr_lis = np.zeros((len(poi_rate_bkg_gpe),len(poi_rate_bkg_stn),len(gpe_ratio_psdb),len(stn_ratio_psdb)))
big_stn_spec_lis = np.zeros((len(poi_rate_bkg_gpe),len(poi_rate_bkg_stn),len(gpe_ratio_psdb),len(stn_ratio_psdb)))
big_gpe_spec_lis = np.zeros((len(poi_rate_bkg_gpe),len(poi_rate_bkg_stn),len(gpe_ratio_psdb),len(stn_ratio_psdb)))

# ssbn_effect
path2 = '../ssbn_effect/output/'
for i,gpe_ip in enumerate(poi_rate_bkg_gpe):
	for j,stn_ip in enumerate(poi_rate_bkg_stn):
		for k,gpe_rat in enumerate(gpe_ratio_psdb):
			for l,stn_rat in enumerate(stn_ratio_psdb):
				prefix = "psdbEffect_inh_ip_"+str(gpe_ip)+"_exc_ip_"+str(stn_ip)+"_inh_ratio_"+str(gpe_rat)+"_exc_ratio_"+str(stn_rat)
				data1 = pickle.load(open(path2+prefix+".pickle","r"))
				big_gpe_spec_lis[i][j][k][l] = np.nanmean(data1["freqSpec_gpe"][0])
				big_stn_spec_lis[i][j][k][l] = np.nanmean(data1["freqSpec_stn"][0])
				mean_gpe = np.array(np.mean(data1["psth_gpe"][0],0))
				mean_stn = np.mean(data1["psth_stn"][0],0)
				big_gpe_fr_lis[i][j][k][l] = np.mean(mean_gpe[1:])
				big_stn_fr_lis[i][j][k][l] = np.mean(mean_stn[1:])

im3 = plot_bursting_regime_wise(900.0,1400.0,ax2,big_gpe_fr_lis,big_stn_fr_lis,big_stn_spec_lis,'2',2,3) # border
im3 = plot_bursting_regime_wise(700.0,1000.0,ax2,big_gpe_fr_lis,big_stn_fr_lis,big_stn_spec_lis,"2'",2,2.5) # border1
im3 = plot_bursting_regime_wise(500.0,1000.0,ax2,big_gpe_fr_lis,big_stn_fr_lis,big_stn_spec_lis,"1",2,3) # osc
im3 = plot_bursting_regime_wise(500.0,1400.0,ax2,big_gpe_fr_lis,big_stn_fr_lis,big_stn_spec_lis,"1'",2,3) # osc1
im3 = plot_bursting_regime_wise(1300.0,1400.0,ax2,big_gpe_fr_lis,big_stn_fr_lis,big_stn_spec_lis,"3",2,3) # non-osc
im3 = plot_bursting_regime_wise(1100.0,1200.0,ax2,big_gpe_fr_lis,big_stn_fr_lis,big_stn_spec_lis,"3'",1,1) # non-osc1

inset_axes = inset_axes(ax2, 
                    width="40%", # width = 30% of parent_bbox
                    height=2.5, # height : 1 inch
                    loc=2)

fontDict = {'size':12, 'weight':'bold','color':'red'}
inset_axes.pcolor(np.nanmean(array(big_gpe_spec_lis_mult),0), cmap = cm.BuGn,vmin=0.3,vmax=0.7,alpha=0.8)
x1 = 3.5;x2=8.5;y1=0.;y2=8.
m = (y2-y1)/(x2-x1)
c = y1-m*x1
x_gen = np.arange(x1,x2+0.1,0.1)
y_gen = m*x_gen+c
inset_axes.plot(x_gen,y_gen,'r--',linewidth=2.5)
inset_axes.text(0.5,3, "Oscillatory",fontdict=fontDict)
inset_axes.text(8,3.5, "Non-",fontdict=fontDict)
inset_axes.text(7.4,0.5, "Oscillatory",fontdict=fontDict)
inset_axes.text(2.3,0.1, "1",fontdict=fontDict)
inset_axes.text(2.3,5.1, "1'",fontdict=fontDict)
inset_axes.text(6.5,4.1, "2",fontdict=fontDict)
inset_axes.text(4.3,0.2, "2'",fontdict=fontDict)
inset_axes.text(10.3,4.1, "3",fontdict=fontDict)
inset_axes.text(8.3,2.1, "3'",fontdict=fontDict)
inset_axes.set_ylim(0,8)
inset_axes.set_xticks([])
inset_axes.set_yticks([])
boxcolor='black'
lw = 2.0
inset_axes.add_patch(patches.Rectangle((2, 0),1,1, fill=False, edgecolor=boxcolor, linewidth=lw, alpha=1.0))
inset_axes.add_patch(patches.Rectangle((2, 5),1,1, fill=False, edgecolor=boxcolor, linewidth=lw, alpha=1.0))
inset_axes.add_patch(patches.Rectangle((6, 4),1,1, fill=False, edgecolor=boxcolor, linewidth=lw, alpha=1.0))
inset_axes.add_patch(patches.Rectangle((4, 0),1,1, fill=False, edgecolor=boxcolor, linewidth=lw, alpha=1.0))
inset_axes.add_patch(patches.Rectangle((10, 4),1,1, fill=False, edgecolor=boxcolor, linewidth=lw, alpha=1.0))
inset_axes.add_patch(patches.Rectangle((8, 2),1,1, fill=False, edgecolor=boxcolor, linewidth=lw, alpha=1.0))






cb5 = fig.colorbar(im3,ax=ax2)
for x in cb5.ax.yaxis.get_ticklabels():
	x.set_fontweight('bold')

cb5.ax.yaxis.set_label_text("Spectral Entropy (SE)",fontsize=20,fontweight='bold')


ax2.set_xlim(23,52)
fig.subplots_adjust(right = 0.97,left=0.08,top=0.97,bottom=0.08)

fig.savefig("figs/SE_graph_all.pdf",dpi=600)
fig.savefig("figs/SE_graph_all.png",dpi=600)


pl.show()


