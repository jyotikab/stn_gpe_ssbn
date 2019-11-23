
#--------------------------------------------------------------------------------------------
# This script reproduces Figure 3 in the manuscript, expects "psdbEffect_inh_ip_"+prefix1+str(gpe_ip)+"_exc_ip_"+str(stn_ip)+"_inh_ratio_"+str(gpe_rat)+"_exc_ratio_"+str(stn_rat), files from simulations different proprtion of STn and GPe bursty neurons (gpe_rat:0-1.0, stn_rat:0-1.0) and at least 3 firing rate combinations:np.array([[500.0,1000.0 ],[900.0,1400.0],[1300.0,1400.0]])
#--------------------------------------------------------------------------------------------

import numpy as np
import itertools
import sys
home_directory = ""
sys.path.append(home_directory+'/common/')
import analyze_data_ssbn as adata
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
import cPickle as pickle
import matplotlib.pyplot as plt
import pdb
import seaborn as sns
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FormatStrFormatter
import matplotlib.patches as patches

num_threads = 8

# "" for Fig 3
prefix1 = sys.argv[1]


def get_sim_3d(match_pars):
	print path1
	res = adata.analyze_data(prefix1,data_path=path1)
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

gpe_ratio_psdb = np.arange(0.,1.1,0.1)
stn_ratio_psdb = np.arange(0.,1.1,0.1)


path1 = os.getcwd()+'/output/'


GPeSTNpairs = np.array([[500.0,1000.0 ],[900.0,1400.0],[1300.0,1400.0]])	# 0th is non-monotonic, 1st is all oscillatory, 2nd is non-oscillatory regimes
gpe_ratio_psdb = np.arange(0.,1.1,0.1)
stn_ratio_psdb = np.arange(0.,1.1,0.1)
poi_rate_bkg_gpe = GPeSTNpairs[:,0]
poi_rate_bkg_stn = GPeSTNpairs[:,1]
big_gpe_fr_lis = np.zeros((len(poi_rate_bkg_gpe),len(gpe_ratio_psdb),len(stn_ratio_psdb)))

big_stn_fr_lis = np.zeros((len(poi_rate_bkg_gpe),len(gpe_ratio_psdb),len(stn_ratio_psdb)))
big_stn_spec_lis = np.zeros((len(poi_rate_bkg_gpe),len(gpe_ratio_psdb),len(stn_ratio_psdb)))
big_gpe_spec_lis = np.zeros((len(poi_rate_bkg_gpe),len(gpe_ratio_psdb),len(stn_ratio_psdb)))
big_stn_freq_lis = np.zeros((len(poi_rate_bkg_gpe),len(gpe_ratio_psdb),len(stn_ratio_psdb)))
big_gpe_freq_lis = np.zeros((len(poi_rate_bkg_gpe),len(gpe_ratio_psdb),len(stn_ratio_psdb)))




data = pickle.load(open("../rate_effect/Data.pickle","r"))
poi_rate_bkg_stn_all=data["bck_rate_stn"]  
poi_rate_bkg_gpe_all=data["bck_rate_gpe"]  
big_gpe_fr_lis_mult  =data["gpe_rate"]  
big_stn_fr_lis_mult  =data["stn_rate"]  
big_gpe_spec_lis_mult=data["gpe_spec"]  
big_stn_spec_lis_mult=data["stn_spec"]  

fontDict = {'size':15, 'weight':'bold','color':'red'}

fig = plt.figure(figsize=(18,17))


top = 0.7
bottom = 0.05
gs0 = GridSpec(1,3)
gridWidth = 0.65
gridSpace = 0.05
left0 = 0.18
hspace=0.05
gs0.update(left=left0,right=left0+gridWidth,top=top+0.27,bottom=top+2.5*hspace,wspace=0.05,hspace=hspace)
gs0Hands =[]
gs0Hands.append(plt.subplot(gs0[0,1]))
poi_rate_bkg_gpe_all = np.append(poi_rate_bkg_gpe_all,1500.)
poi_rate_bkg_stn_all = np.append(poi_rate_bkg_stn_all,1800.)

im3 = gs0Hands[0].pcolor(poi_rate_bkg_gpe_all, poi_rate_bkg_stn_all, np.nanmean(array(big_gpe_spec_lis_mult),0), cmap = cm.BuGn,vmin=0.3,vmax=0.7,alpha=0.8)
gs0Hands[0].set_xlim(poi_rate_bkg_gpe_all[0], poi_rate_bkg_gpe_all[-1])
gs0Hands[0].set_ylim(poi_rate_bkg_stn_all[0], poi_rate_bkg_stn_all[-1])
gs0Hands[0].set_xticks(poi_rate_bkg_gpe_all[::3])
gs0Hands[0].set_yticks(poi_rate_bkg_stn_all[::3])
gs0Hands[0].set_ylabel('STN input',fontsize=15,fontweight='bold')
gs0Hands[0].set_xlabel('GPe input',fontsize=15,fontweight='bold')
for x in gs0Hands[0].get_yticklabels():
	x.set_fontsize(12)
	x.set_fontweight('bold')
for x in gs0Hands[0].get_xticklabels():
	x.set_fontsize(12)
	x.set_fontweight('bold')

x1 = 650.;x2=1150.;y1=1000.;y2=1800.
m = (y2-y1)/(x2-x1)
c = y1-m*x1
x_gen = np.arange(x1,x2+100,100)
y_gen = m*x_gen+c
gs0Hands[0].plot(x_gen,y_gen,'r--',linewidth=2.5)
gs0Hands[0].text(320,1450, "Oscillatory",fontdict=fontDict)
gs0Hands[0].text(950,1250, "Non-",fontdict=fontDict)
gs0Hands[0].text(850,1150, "Oscillatory",fontdict=fontDict)
step = 20
  
gs0Hands[0].text(500+step,1000+step, "1",fontdict=fontDict)
gs0Hands[0].text(900+step*2,1400+step*0.5, "2",fontdict=fontDict)
gs0Hands[0].text(1300+step,1400+step, "3",fontdict=fontDict)
gs0Hands[0].text(180,1850,"A",fontweight='bold',fontsize=18)
boxcolor='black'
lw = 2.0
gs0Hands[0].add_patch(patches.Rectangle((500, 1000),100,100, fill=False, edgecolor=boxcolor, linewidth=lw, alpha=1.0))
gs0Hands[0].add_patch(patches.Rectangle((900, 1400),100,100, fill=False, edgecolor=boxcolor, linewidth=lw, alpha=1.0))
gs0Hands[0].add_patch(patches.Rectangle((1300, 1400),100,100, fill=False, edgecolor=boxcolor, linewidth=lw, alpha=1.0))



gs1 = GridSpec(2,2)
left = 0.08
gridWidth = 0.25
gridSpace = 0.05
gs1.update(left=left, right=left+gridWidth,top=top,bottom=bottom, wspace=0.05)
gs1Hands =[]				
gs1Hands.append(plt.subplot(gs1[0,0]))	# All oscillatory
gs1Hands.append(plt.subplot(gs1[0,1]))
gs1Hands.append(plt.subplot(gs1[1,0]))
gs1Hands.append(plt.subplot(gs1[1,1]))

left2 = left+gridWidth+gridSpace
gs2 = GridSpec(2,2)
gs2.update(left=left2, right=left2+gridWidth,top=top,bottom=bottom, wspace=0.05)
gs2Hands =[]
gs2Hands.append(plt.subplot(gs2[0,0])) # Non-monotonic
gs2Hands.append(plt.subplot(gs2[0,1]))
gs2Hands.append(plt.subplot(gs2[1,0]))
gs2Hands.append(plt.subplot(gs2[1,1]))

left3=left2+gridWidth+gridSpace
gs3 = GridSpec(2,2)
gs3.update(left=left3, right=left3+gridWidth,top=top,bottom=bottom, wspace=0.05)
gs3Hands =[]
gs3Hands.append(plt.subplot(gs3[0,0])) # All non-oscillatory
gs3Hands.append(plt.subplot(gs3[0,1]))
gs3Hands.append(plt.subplot(gs3[1,0]))
gs3Hands.append(plt.subplot(gs3[1,1]))


cbar1s = [ fig.add_axes([left+0.015,top+0.02,0.09,0.02]),  fig.add_axes([left2+0.015,top+0.02,0.09,0.02 ]), fig.add_axes([left3+0.015,top+0.02,0.09,0.02 ])]
cbar12s = [ fig.add_axes([(left+gridWidth/2.)+0.015,top+0.02,0.09,0.02]), fig.add_axes([(left2+gridWidth/2.)+0.015,top+0.02,0.09,0.02]), fig.add_axes([(left3+gridWidth/2.)+0.015,top+0.02,0.09,0.02])]
i=0
for gpe_ip,stn_ip in zip(poi_rate_bkg_gpe,poi_rate_bkg_stn) :
		for k,gpe_rat in enumerate(gpe_ratio_psdb):
			for l,stn_rat in enumerate(stn_ratio_psdb):
				prefix = "psdbEffect_inh_ip_"+prefix1+str(gpe_ip)+"_exc_ip_"+str(stn_ip)+"_inh_ratio_"+str(gpe_rat)+"_exc_ratio_"+str(stn_rat)
				data = pickle.load(open("output/"+prefix+".pickle","r"))
				big_gpe_spec_lis[i][k][l] = np.nanmean(data["freqSpec_gpe"][0])
				big_stn_spec_lis[i][k][l] = np.nanmean(data["freqSpec_stn"][0])
					
				mean_gpe = np.array(np.mean(data["psth_gpe"][0],0))
				mean_stn = np.mean(data["psth_stn"][0],0)
				big_gpe_fr_lis[i][k][l] = np.mean(mean_gpe[1:])
				big_stn_fr_lis[i][k][l] = np.mean(mean_stn[1:])

				big_gpe_freq_lis[i][k][l] = np.mean(data["freqSpec_gpe"][1:])
				big_stn_freq_lis[i][k][l] = np.mean(data["freqSpec_gpe"][1:])


		print i
		if i == 0:
			gsHands = np.copy(gs1Hands)
			figH = gs1
		elif i == 1:
			gsHands = np.copy(gs2Hands)
			figH = gs2
		else:
			gsHands = np.copy(gs3Hands)
			figH = gs3

		gsHands[2].pcolor(stn_ratio_psdb,gpe_ratio_psdb,big_gpe_spec_lis[i],cmap=cm.BuGn,vmin=0.3,vmax=0.7)
		gsHands[2].set_xlabel("STN bursting ratio",fontsize=15,fontweight='bold')
		gsHands[2].set_xlim(0,1.0)
		gsHands[2].set_ylim(0,1.0)
		if i == 0:
			for y in gsHands[2].get_yticklabels():
				y.set_visible(True)
				y.set_fontweight('bold')
				y.set_fontsize(12)

			gsHands[2].set_ylabel("GPe bursting ratio",fontsize=15,fontweight='bold')
		else:
			for y in gsHands[2].get_yticklabels():
				y.set_visible(False)
	
		for y in gsHands[2].get_xticklabels():
			y.set_visible(False)
		for y in gsHands[2].get_xticklabels()[1::2]:
			y.set_visible(True)
			y.set_fontweight('bold')
			y.set_fontsize(12)
			
		gsHands[2].set_title("GPe SE",fontsize=15,fontweight='bold')


		img4=gsHands[3].pcolor(stn_ratio_psdb,gpe_ratio_psdb,big_stn_spec_lis[i],cmap=cm.BuGn,vmin=0.3,vmax=0.7)
		gsHands[3].set_xlim(0,1.0)
		gsHands[3].set_ylim(0,1.0)
		for y in gsHands[3].get_yticklabels():
			y.set_visible(False)
		for y in gsHands[3].get_xticklabels():
			y.set_visible(False)
		gsHands[3].set_title("STN SE",fontsize=15,fontweight='bold')
		
		if i == 2:
			cbar_ax3 = fig.add_axes([left3+gridWidth+0.015, 0.05, 0.015, top*0.4])
			cb1 = fig.colorbar(img4,cax=cbar_ax3)
			cb1.ax.set_yticklabels(cb1.ax.get_yticklabels(),fontsize=14,fontweight='bold')
			

		minGPe = np.min(big_gpe_fr_lis[i])
		maxGPe = np.max(big_gpe_fr_lis[i])
		print "gpe_ip",gpe_ip
		print "stn_ip",stn_ip
		print "minGPe",minGPe
		print "maxGPe",maxGPe

		img3 = gsHands[0].pcolor(stn_ratio_psdb,gpe_ratio_psdb,big_gpe_fr_lis[i],cmap=cm.BuGn)#,vmin=20,vmax=60)

		img3.set_clim(np.round([minGPe,maxGPe],2))
		for y in gsHands[0].get_xticklabels():
			y.set_visible(False)
		cbar1 =  cbar1s[i]

		midGPe = ((maxGPe-minGPe)/2.)+minGPe
		cbarlabels = [round(x,2) for x in list([minGPe,midGPe,maxGPe])]
		print cbarlabels
		fig.colorbar(img3,cax=cbar1,orientation='horizontal',ticks=cbarlabels,format="%.1f")

		cbar1.set_xticklabels([str(round(x,1)) for x in cbarlabels] , fontweight='bold',fontsize=12)
		print [x.get_text() for x in cbar1.get_xticklabels()]
		if i == 0:
			for y in gsHands[0].get_yticklabels():
				y.set_visible(True)
				y.set_fontweight('bold')
				y.set_fontsize(12)

			gsHands[0].set_ylabel("GPe bursting ratio",fontweight='bold',fontsize=15)
		else:
			for y in gsHands[0].get_yticklabels():
				y.set_visible(False)


		gsHands[0].set_title("GPe rates",fontsize=15,fontweight='bold',y=1.18)	
		if i == 0:
			gsHands[0].text(-0.3,1.2,"B",fontweight='bold',fontsize=18)
		if i == 1:
			gsHands[0].text(-0.3,1.2,"C",fontweight='bold',fontsize=18)
		if i == 2:
			gsHands[0].text(-0.3,1.2,"D",fontweight='bold',fontsize=18)
		minSTN = np.min(big_stn_fr_lis[i])
		maxSTN = np.max(big_stn_fr_lis[i])

		print "minSTN",minSTN
		print "maxSTN",maxSTN

		img32 = gsHands[1].pcolor(stn_ratio_psdb,gpe_ratio_psdb,big_stn_fr_lis[i],cmap=cm.BuGn)#,vmin=20,vmax=60)

		img32.set_clim(np.round([minSTN,maxSTN],2))
		for y in gsHands[1].get_yticklabels():
			y.set_visible(False)
		for y in gsHands[1].get_xticklabels():
			y.set_visible(False)
		cbar12 = cbar12s[i]

		midSTN = ((maxSTN-minSTN)/2.)+minSTN
		cbarlabels12 =[ round(x,2)  for x in   list([minSTN,midSTN,maxSTN])]
		print cbarlabels12

		fig.colorbar(img32,cax=cbar12,orientation='horizontal',ticks=cbarlabels12,format="%.1f")
		cbar12.xaxis.set_ticklabels( [ str(round(x,1)) for x in  cbarlabels12], fontweight='bold')
		for y in cbar12.get_xticklabels():
			y.set_fontsize(12)
			#cbar12.set_xticklabels(str(np.round(float(y.get_item()),1)))


		gsHands[1].set_title("STN rates",fontsize=15,fontweight='bold',y=1.18)	


		i+=1

fig.savefig("figs/Fig5"+prefix1+".pdf")
fig.savefig("figs/Fig5"+prefix1+".png")
