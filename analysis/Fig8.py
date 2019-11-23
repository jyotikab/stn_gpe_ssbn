
#--------------------------------------------------------------------------------------
# Reproduces Figure 8
# Only plots the background, with legend etc
# The arrows were added by inkscape
#--------------------------------------------------------------------------------------

import numpy as np
import itertools
import sys
home_directory = ""
sys.path.append(home_directory+'/common/')
import analyze_data_ssbn as adata
import analysis_funcs_wo_framework as anal_wo
reload(adata)

import pylab as pl
import matplotlib.cm as cm
import os
import cPickle as pickle
import matplotlib as mpl
import pdb
import seaborn as sns
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FormatStrFormatter
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerPatch

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['mathtext.fontset'] = 'stix'
gpe_ratio_psdb = np.arange(0.,1.1,0.1)
stn_ratio_psdb = np.arange(0.,1.1,0.1)
pl.rc('axes', linewidth=4.5)

path1 = os.getcwd()+'/output/'
osc_color="olive"
nonosc_color="lightskyblue"

GPeSTNpairs = np.array([[500.0,1000.0 ],[900.0,1400.0],[1300.0,1400.0]])	# 0th is non-monotonic, 1st is all oscillatory, 2nd is non-oscillatory regimes
gpe_ratio_psdb = np.arange(0.,1.1,0.1)
stn_ratio_psdb = np.arange(0.,1.1,0.1)
poi_rate_bkg_gpe = GPeSTNpairs[:,0]
poi_rate_bkg_stn = GPeSTNpairs[:,1]
big_gpe_fr_lis = np.zeros((len(poi_rate_bkg_gpe),len(gpe_ratio_psdb),len(stn_ratio_psdb)))

big_stn_fr_lis = np.zeros((len(poi_rate_bkg_gpe),len(gpe_ratio_psdb),len(stn_ratio_psdb)))
big_stn_fr_lis = np.zeros((len(poi_rate_bkg_gpe),len(gpe_ratio_psdb),len(stn_ratio_psdb)))
big_gpe_fr_lis = np.zeros((len(poi_rate_bkg_gpe),len(gpe_ratio_psdb),len(stn_ratio_psdb)))
big_stn_freq_lis = np.zeros((len(poi_rate_bkg_gpe),len(gpe_ratio_psdb),len(stn_ratio_psdb)))
big_gpe_freq_lis = np.zeros((len(poi_rate_bkg_gpe),len(gpe_ratio_psdb),len(stn_ratio_psdb)))


def draw_arrow(arrow_color,x1,y1,x2,y2,radius,nuc="GPe",scale=2.5):
	style="Simple,tail_width=5,head_width=18,head_length=25"
	kw = dict(arrowstyle=style, color=arrow_color)
	if nuc == "STN":
		a3 = patches.FancyArrowPatch((x1,y1), (x2,y2),connectionstyle="arc3,rad="+str(radius*scale), **kw)
	else:
		a3 = patches.FancyArrowPatch((x1,y1), (x2,y2),connectionstyle="arc3,rad=0", **kw)
	return a3

def plot_ei_balance(data_regime,filename,zoom_in=False,zoom_lims=[],t3=[],israte=0,xticks_vis=True,yticks_vis=True,subplot_label='a'):

	class HandlerEllipse(HandlerPatch):
	    def create_artists(self, legend, orig_handle,
			       xdescent, ydescent, width, height, fontsize, trans):
		center = 0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent
		p = mpatches.Ellipse(xy=center, width=width + xdescent,
				     height=height + ydescent)
		self.update_prop(p, orig_handle, legend)
		p.set_transform(trans)
		return [p]



	def plot_ei_regime(ax,data,JEI,JII,filename,fig,zoom_in,zoom_lims,scales):
		JII_b = np.array(data["points_fr_gpe"])*0.725*40*10
		JEI_b = np.array(data["points_fr_stn"])*1.2*23*5
		

		ax1 = ax.twinx().twiny()

		#if zoom_in == True:
		#	ax1.margins(x=zoom_lims[0],y=zoom_lims[1])
		for k,(x_row,y_row) in enumerate(zip(JEI_b,JII_b)):
			arr_col = cm.inferno(float(k)/len(JEI_b))
			for l,(x_col,y_col) in enumerate(zip(x_row,y_row)):
				if l < len(x_row)-1: 
					ax1.arrow(x_col,y_col,x_row[l+1]-x_col,y_row[l+1]-y_col,head_width=50,length_includes_head=True,fc = arr_col, ec = arr_col,linestyle='-',linewidth=1.0)
			if k%2 == 0:
				ax1.plot(x_col,y_col,marker='o',markersize=12,color=arr_col,label='gpe ratio:'+str(gpe_rat[k]),alpha=0.8)
			else:
				ax1.plot(x_col,y_col,marker='o',markersize=12,color=arr_col,alpha=0.8)
			ax1.plot(x_row[0],y_row[0],marker='*',markersize=12,color=arr_col,alpha = 0.8)

		handles, labels = ax1.get_legend_handles_labels()	
		ax1.set_xlim(ax.get_xlim())
		ax1.set_ylim(ax.get_ylim())
		ax1.set_xticklabels([])
		ax1.set_yticklabels([])
				
		
		

	fontDict = {'size':18, 'weight':'bold','color':'red'}	 
	start_point = data_A["start_point"]
	points_fr_gpe = data_A["points_fr_gpe"]
	points_fr_stn = data_A["points_fr_stn"]
	gpe_rat = data_A["gpe_rat"]
	stn_rat = data_A["stn_rat"]

	data = pickle.load(open("../Arvind/data/Data.pickle","r"))
	poi_rate_bkg_stn_all=data["bck_rate_stn"]  
	poi_rate_bkg_gpe_all=data["bck_rate_gpe"]  
	big_gpe_fr_lis_mult  =data["gpe_rate"]  
	big_stn_fr_lis_mult  =data["stn_rate"]  
	big_gpe_spec_lis_mult  =data["gpe_spec"]  
	big_stn_spec_lis_mult  =data["stn_spec"]  
	

	
	JII = np.mean(np.array(big_gpe_fr_lis_mult)*0.725*40*10,0)
	JEI = np.mean(np.array(big_stn_fr_lis_mult)*1.2*23*5,0)

	osc_ind = np.where(np.mean(np.array(big_gpe_spec_lis_mult),0) <=0.4)
	non_osc_ind = np.where(np.mean(np.array(big_gpe_spec_lis_mult),0) > 0.4)

	a_ind_g,a_ind_s = np.where(poi_rate_bkg_gpe_all == 500.)[0],np.where(poi_rate_bkg_stn_all == 1000.)[0]
	b_ind_g,b_ind_s = np.where(poi_rate_bkg_gpe_all == 900.)[0],np.where(poi_rate_bkg_stn_all == 1400.)[0]
	c_ind_g,c_ind_s = np.where(poi_rate_bkg_gpe_all == 1300.)[0],np.where(poi_rate_bkg_stn_all == 1400.)[0]
	
	print np.mean(np.array(big_gpe_fr_lis_mult),0)[a_ind_s,a_ind_g]
	t3.set_ylim(8000,15000)
	
	t3.plot(JEI[non_osc_ind],JII[non_osc_ind],'o',color=nonosc_color,markersize=10,alpha=0.5,markeredgecolor=nonosc_color,markeredgewidth=1.5)
	t3.plot(JEI[osc_ind],JII[osc_ind],'o',color=osc_color,markersize=10,alpha=0.5,markeredgecolor=osc_color,markeredgewidth=1.5)

	t3.plot(JEI[a_ind_s,a_ind_g],JII[a_ind_s,a_ind_g],'*',color='r',markersize=10)
	t3.text(JEI[a_ind_s,a_ind_g]+100,JII[a_ind_s,a_ind_g],"1",fontdict = fontDict)

	t3.plot(JEI[b_ind_s,b_ind_g],JII[b_ind_s,b_ind_g],'*',color='r',markersize=10)
	t3.text(JEI[b_ind_s,b_ind_g]+100,JII[b_ind_s,b_ind_g],"2",fontdict = fontDict)

	t3.plot(JEI[c_ind_s,c_ind_g],JII[c_ind_s,c_ind_g],'*',color='r',markersize=10)
	t3.text(JEI[c_ind_s,c_ind_g]+100,JII[c_ind_s,c_ind_g]-50,"3",fontdict = fontDict)
	t3.set_xlabel('Effective excitation to a GPe neuron',fontsize=25,fontweight='bold',labelpad=20)
	t3.set_ylabel('Effective inibition to a GPe neuron',fontsize=25,fontweight='bold',labelpad=20)

	for x in t3.get_xticklabels():
		x.set_visible(False)
	for x in t3.get_yticklabels():
		x.set_visible(False)

	t3.arrow(1300,13000,800,0,head_width=150, fc='k', ec='k',linewidth=4.5)	
	t3.text(1100,13200,"Increase in STN firing rate",color='k',fontsize=20,fontweight='bold')	

	t3.text(3500,14000,"Oscillatory",color='k',fontsize=20,fontweight='bold')
	t3.text(500,14000,"Non-oscillatory",color='k',fontsize=20,fontweight='bold')


	arrow_colors = [ cm.brg(x)   for x in np.linspace(0,1,10) ]
	print arrow_colors
	for x in arrow_colors:
		print mpl.colors.rgb2hex(x)	


	gpe_b_ind = 2
	stn_b_ind = 6
	rect = patches.Rectangle((310,12250),100,120,linewidth=4,edgecolor=arrow_colors[gpe_b_ind],facecolor='none',joinstyle='round')
	t3.add_patch(rect)
	rect = patches.Rectangle((650,13100),100,120,linewidth=4,edgecolor=arrow_colors[stn_b_ind],facecolor='none',joinstyle='round')
	t3.add_patch(rect)


	rect = patches.Rectangle((1350,11000),100,120,linewidth=4,edgecolor=arrow_colors[gpe_b_ind],facecolor='none',joinstyle='round')
	t3.add_patch(rect)
	elli = patches.Ellipse((1740,12000),100,120,linewidth=4,edgecolor=arrow_colors[stn_b_ind],facecolor='none',joinstyle='round')
	t3.add_patch(elli)

	elli = patches.Ellipse((1850,8800),100,120,linewidth=4,edgecolor=arrow_colors[gpe_b_ind],facecolor='none',joinstyle='round')
	t3.add_patch(elli)
	elli = patches.Ellipse((2300,9800),100,120,linewidth=4,edgecolor=arrow_colors[stn_b_ind],facecolor='none',joinstyle='round')
	t3.add_patch(elli)


	circ = mpatches.Circle((0, 0), 1, facecolor="none", edgecolor='k', linewidth=4)

	legend_elements = [Line2D([0], [0], color=arrow_colors[gpe_b_ind], lw=4, label='GPe bursting'), Line2D([0], [0], color=arrow_colors[stn_b_ind], lw=4, label='STN bursting'),  patches.Patch(facecolor='none', edgecolor='k',label='Non-oscillatory state',joinstyle='round',linewidth=4),circ ]


	t3.fill([0,1100,1500,2000,0,0],[8000,8000,11000,15000,15000,8000],color=nonosc_color,alpha=0.15)
	t3.fill([1100,2000,5000,5000],[8000,15000,15000,8000],color=osc_color,alpha=0.15)


	t3.set_xlim(np.min(JEI),np.max(JEI))
	t3.set_ylim(np.min(JII),np.max(JII))

	for x in t3.get_xticklabels():
		x.set_fontsize(12)
		x.set_fontweight('bold')

	for x in t3.get_yticklabels():
		x.set_fontsize(12)
		x.set_fontweight('bold')


	if xticks_vis == False:
		for x in t3.get_xticklabels():
			x.set_visible(False)
		t3.set_xlabel("")
	if yticks_vis == False:
		for x in t3.get_yticklabels():
			x.set_visible(False)
		t3.set_ylabel("")
	

	plt.legend(handles=legend_elements,labels=["GPe bursting", "STN bursting", "Non-oscillatory state", "Oscillatory state"], loc='lower right', prop={'size':20,"weight":"bold"},handler_map={mpatches.Circle: HandlerEllipse()} )
	if israte == 0:
		JII_b = np.array(data_regime["points_fr_gpe"])*0.725*40*10
		JEI_b = np.array(data_regime["points_fr_stn"])*1.2*23*5


		scale_osc = JEI[a_ind_s,a_ind_g]/JEI_b[a_ind_s,a_ind_g], JII[a_ind_s,a_ind_g]/JII_b[a_ind_s,a_ind_g]
		scale_border = JEI_b[b_ind_s,b_ind_g]/JEI[b_ind_s,b_ind_g], JII_b[b_ind_s,b_ind_g]/JII[b_ind_s,b_ind_g]
		scale_nonosc = JEI_b[c_ind_s,c_ind_g]/JEI[c_ind_s,c_ind_g], JII_b[c_ind_s,c_ind_g]/JII[c_ind_s,c_ind_g]

		plot_ei_regime(t3,data_regime,JEI,JII,filename,fig2,zoom_in,zoom_lims,[scale_osc,scale_nonosc,scale_border])


def read_data():
	i=0

	for gpe_ip,stn_ip in zip(poi_rate_bkg_gpe,poi_rate_bkg_stn) :
			for k,gpe_rat in enumerate(gpe_ratio_psdb):
				for l,stn_rat in enumerate(stn_ratio_psdb):
					prefix = "psdbEffect_inh_ip_"+str(gpe_ip)+"_exc_ip_"+str(stn_ip)+"_inh_ratio_"+str(gpe_rat)+"_exc_ratio_"+str(stn_rat)
					data1 = pickle.load(open("output/"+prefix+".pickle","r"))
					mean_gpe = np.array(np.mean(data1["psth_gpe"][0],0))
					mean_stn = np.mean(data1["psth_stn"][0],0)

					big_gpe_fr_lis[i][k][l] = np.mean(mean_gpe[0:])
					big_stn_fr_lis[i][k][l] = np.mean(mean_stn[0:])
					
					if gpe_rat == 0.0 and stn_rat == 0.0:
						print big_gpe_fr_lis[i][k][l]
						print big_stn_fr_lis[i][k][l]
			i+=1

	return big_gpe_fr_lis, big_stn_fr_lis

fig2 = pl.figure(figsize=(16,12))
t1 = fig2.add_subplot(111)


filename="EI_balance_summary"

big_gpe_fr_lis, big_stn_fr_lis = read_data()

a_ind = 0
b_ind = 1
c_ind = 2


data_A = dict()
data_A["start_point"] = (500.,1000.)
data_A["points_fr_gpe"] = big_gpe_fr_lis[a_ind]
data_A["points_fr_stn"] = big_stn_fr_lis[a_ind]
data_A["gpe_rat"] = gpe_ratio_psdb
data_A["stn_rat"] = stn_ratio_psdb

plot_ei_balance(data_A,"EI_balance_rate",t3=t1,israte=1,subplot_label='A')



fig2.subplots_adjust(left=0.08,right=0.95,bottom=0.05,wspace=0.1,hspace=0.12)
fig2.savefig("figs/"+filename+".png")
fig2.savefig("figs/"+filename+".pdf")
fig2.savefig("figs/"+filename+".svg")


pl.show()
