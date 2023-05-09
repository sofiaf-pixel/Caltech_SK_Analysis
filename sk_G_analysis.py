#!/usr/bin/env python

import matplotlib.pyplot as plt
import scipy.io as sio
import numpy as np
import math
import pickle
from scipy.optimize import curve_fit
import sys, os
import os.path
import shutil
import pylab as pl
import mce_data
#sys.path.insert(0, "/home/cheng/analysis/load_curve/sk_dark/scripts_house")
import calib_SK
import get_LCs as lc
import single_pr_darkrun
import ba40_ModuleMapping
from Gtemplists import gettemplists
templists = gettemplists()

#RUN = 'L1'
#DATE='20230302'

RUN=sys.argv[1]
DATE=sys.argv[2]

fit_beta=False

rnpsat = 0.040
templist = templists[DATE]

def GTcModel(T,G,Tc,beta):
	return G/(beta+1)*(Tc**(beta+1)-T**(beta+1))/Tc**beta

def GTcModel2(T,G,Tc):
	#beta = 2.
	return G/(2+1)*(Tc**(2+1)-T**(2+1))/Tc**2

def main():

	Gc=[]
	Tc=[]
	beta=[]
	RnTi=[]
	Psat308=[]
	G450=[]


	lcdata = {}
	prdata = {}
	colors = pl.cm.jet(np.linspace(0,1,len(templist)))
	rpsat=np.zeros(len(templist))

	for temp in templist:
		#===================================#
		# Loading data
		#===================================#
		in_path  = '%s/'%(DATE)
		filename = 'LC_G_sk_FPU'+str(temp)+'mK_datamode1_source1'
		out_path_main = 'output/%s/'%(DATE)
		if not os.path.isdir(out_path_main):
			os.makedirs(out_path_main)
		#out_path = out_path_main + filename + '/'
		#if not os.path.isdir(out_path):
			#os.makedirs(out_path)
		out_Gplot = out_path_main + 'L%s_Gplot/'%(RUN)
		if not os.path.isdir(out_Gplot):
			os.makedirs(out_Gplot)

		if fit_beta==True:
			if not os.path.isdir(out_Gplot + 'Fit_beta/'):
				os.makedirs(out_Gplot + 'Fit_beta/')

		datafn = in_path + filename + '/' + filename
		#print(datafn)

		biasfn = datafn + '.bias'
		f = mce_data.MCEFile(datafn)
		dname = os.path.split(datafn)[0]

		#===================================#
		# Calib
		#===================================#

		calib = calib_SK.calib_sk()

		# calib.BIAS_CAL = (calib.V_B_MAX/(2^calib.BITS_BIAS))./(calib.R_BIAS + calib.R_WIRE);
		bias = np.loadtxt(biasfn,skiprows=1)
		biasCalib = bias*calib["BIAS_CAL"][1]

		# calib.FB_CAL = (calib.V_FB_MAX/(2^calib.BITS_FB))./(calib.R_FB+calib.R_WIRE) ./ calib.M_FB;
		y = -1.0*f.Read(row_col=True,unfilter='DC').data
		yCalib =  1.0*f.Read(row_col=True,unfilter='DC').data*calib["FB_CAL"][1]

		nr,nc,nt = y.shape
		rows = np.zeros((nc,nr),dtype=np.int64)
		cols = np.zeros((nc,nr),dtype=np.int64)

		row_i=int(sys.argv[3])
		col_i=int(sys.argv[4])

		Gfn = 'sk_G_'+str(DATE)+'_row'+str(row_i)+'_col'+str(col_i)


		#===================================#
		# Cook LC
		#===================================#
		#for col in range(16):
		for col in [col_i]:
			# for row in range(41):
			for row in [row_i]:
				fitrange = {}

				fitrange["rnti_low"] = 6000.00
				fitrange["rnti_hgh"] = 7000.00
				fitrange["sc_low"] = 0.00
				fitrange["sc_hgh"] = 50.00

				ibias, ites, rnti1 = lc.get_LCs(bias, y[row,col], row, col, calib, fitrange)	#if out_path is included in the arguments, the function will produce plot. Otherwise no. out put dic; see definition. (lc_cook func is no longer available. need to be changed to getLCs in script_house)
				# if row==0&col==0:
				# 	print (ibias)
				# 	print (ites)

				temp_i = templist.index(temp)
				rr, pp, rnti2, psat, rpsat_i= lc.get_PRs(ibias, ites, row, col, calib, rnpsat)
				rpsat[temp_i]=rpsat_i

				if rnti2<0:
					rnti=rnti1
				else:
					rnti=rnti2

				lcdata[str(temp)+'_r'+str(row)+'c'+str(col)] = [ibias, ites]
				prdata[str(temp)+'_r'+str(row)+'c'+str(col)] = [rr, pp, rnti, psat]#, 1/(ites*rr)]

	#============================================#
	# Get psat and rn, do Gfit
	#============================================#

	templist.sort()
	iteslist = [0]*len(templist)
	ibiaslist = [0]*len(templist)
	rnpsatlowlist = [0]*len(templist)
	rnpsathighlist = [0]*len(templist)
	psatlist = [0]*len(templist)
	rntilist = [0]*len(templist)

	DetInfos = {}

	for col in [col_i]:
		for row in [row_i]:

			dinfos = {}
			#im, detcol,detrow,detpol = ba40_ModuleMapping.mce2det(col,row)
			im, detcol,detrow,detpol = 0, 0, 0, 0
			for itt in range(len(templist)):

				psatlist[itt] = prdata[str(templist[itt])+'_r'+str(row)+'c'+str(col)][3]*1.00e12
				rntilist[itt] = prdata[str(templist[itt])+'_r'+str(row)+'c'+str(col)][2]*1.00e3

			#	lcIs = lcdata[str(templist[itt])+'_r'+str(row)+'c'+str(col)]
			#	rr, pp, rnti, psat = lc.get_PRs(lcIs[0], lcIs[1], row, col, calib, rnpsat)
			rnti = np.mean(rntilist)
			xdata = np.linspace(templist[0]-1, templist[-1]+1, 100)

			mask=np.asarray(psatlist)>0

			psatlist_good=np.asarray(psatlist)[mask]
			templist_good=np.asarray(templist)[mask]


			fig = pl.figure(figsize=(20,6.5), dpi=80)
			pl.clf()

			ax = pl.subplot(1, 3, 1)
			pl.suptitle('Row %02d'%row + ' Col %02d'%col)
			pl.xlabel('Ibias [uA]', fontsize=15)
			pl.ylabel('Ites [uA]', fontsize=15)
			pl.xlim(0,500)
			#pl.ylim(0,50)
			for itemp,temp in enumerate(templist_good):
			        pl.plot(lcdata[str(temp)+'_r'+str(row)+'c'+str(col)][0]*1.e6,lcdata[str(temp)+'_r'+str(row)+'c'+str(col)][1]*1.e6, color=colors[itemp], label=str(temp)+'mK')
			plt.tick_params(labelsize=14)
			pl.grid()
			pl.legend(loc=1, prop={'size': 14})
			#pl.axvline(x = rnpsat*1.0e3,color='r', linestyle='--')


			ax = pl.subplot(1, 3, 2)
			pl.suptitle('Row %02d'%row + ' Col %02d'%col)
			pl.xlabel('R [mOhms]', fontsize=15)
			pl.ylabel('P [pW]', fontsize=15)
			# if np.isnan(rnti):
			pl.xlim(0,80)
			# else:
			# 	pl.xlim(0,rnti+20)
			pl.ylim(0,30)
			for itemp,temp in enumerate(templist_good):
				pl.plot(prdata[str(temp)+'_r'+str(row)+'c'+str(col)][0]*1.00e3,prdata[str(temp)+'_r'+str(row)+'c'+str(col)][1]*1.00e12,color=colors[itemp],label=str(temp)+'mK')
			plt.tick_params(labelsize=14)
			pl.grid()
			pl.legend(loc=1, prop={'size': 14})
			pl.axvline(x = rpsat[temp_i]*1000,color='r', linestyle='--')

			ax = pl.subplot(1, 3, 3)
			pl.scatter(templist_good,psatlist_good)
			if 'popt' in locals():
				del popt
				del pcov

			if np.isnan(rnti):
				x_cut=0*[]
				for i,temp11 in enumerate(templist_good):

					x_resist=prdata[str(temp11)+'_r'+str(row)+'c'+str(col)][0]
					mask=(x_resist>rnpsat)&(x_resist<rnpsat+0.005)

					y_power=prdata[str(temp11)+'_r'+str(row)+'c'+str(col)][1]
					if np.mean(y_power[mask])>0 and np.mean(y_power[mask])<20:
						x_resist=prdata[str(temp11)+'_r'+str(row)+'c'+str(col)][0]
						x_resist_mO=x_resist*1000.
						y_power=prdata[str(temp11)+'_r'+str(row)+'c'+str(col)][1]
						if len(y_power)<3:
							x_cut.append(-999999999999999999.9)
						else:
							#print(y_power)
							r_index=np.where(y_power>np.mean(y_power[mask])+1e-13)
							x_cut.append(np.min(x_resist_mO[r_index]))
				x_cut=np.asarray(x_cut)
				if len(x_cut[x_cut>0])<1:
					rnti=0.
				else:
					rnti=np.mean(x_cut[x_cut>0.])


			idx=np.argsort(templist_good)
			templist_good_sorted=templist_good[idx]
			psatlist_good_sorted=psatlist_good[idx]
			psatlist_gg=[]
			templist_gg=[]
			for i in range (0,len(psatlist_good_sorted)-2):
				if psatlist_good_sorted[i]>psatlist_good_sorted[i+1]:
					psatlist_gg.append(psatlist_good_sorted[i])
					templist_gg.append(templist_good_sorted[i])
			psatlist_gg.append(psatlist_good_sorted[-2])
			templist_gg.append(templist_good_sorted[-2])
			#pl.scatter(templist_good_sorted, psatlist_good_sorted, c='y', label='before correction')
			pl.scatter(templist_gg, psatlist_gg, c='r')

			try:

				if fit_beta==True:
					popt, pcov = curve_fit(GTcModel, templist_gg, psatlist_gg)
					Gc = popt[0]*1e3
					Tc = popt[1]
					beta = popt[2]
					rnti = rnti
					Psat308 = psatlist_gg[3]
					#pl.plot(xdata, GTcModel(xdata, *popt), 'g--',label='fit: Gc=%5.3f, Tc=%5.3f, beta=%5.3f' %(Gc, Tc, beta))
					pl.plot(xdata, GTcModel(xdata, *popt), 'g--')

				else:
					popt, pcov = curve_fit(GTcModel2, templist_gg, psatlist_gg)
					Gc = popt[0]*1e3
					Tc = popt[1]
					beta = 2
					rnti = rnti
					Psat308 = psatlist_gg[3]
					#pl.plot(xdata, GTcModel(xdata, *popt), 'g--',label='fit: Gc=%5.3f, Tc=%5.3f, beta=%5.3f' %(Gc, Tc, beta))
					pl.plot(xdata, GTcModel2(xdata, *popt), 'g--')

				pl.text(0.61,0.95,'Gc='+str(round(Gc,1))+' pW/K',transform=ax.transAxes, fontsize=14)
				pl.text(0.61,0.9,'Tc='+str(round(Tc))+' mK',transform=ax.transAxes, fontsize=14)
				pl.text(0.61,0.85,'beta='+str(round(beta,2)),transform=ax.transAxes, fontsize=14)
				pl.text(0.61,0.8,'RnTi='+str(round(rnti))+' mOhm',transform=ax.transAxes, fontsize=14)
				pl.text(0.61,0.75,'Psat308='+str(round(Psat308,1))+' pW',transform=ax.transAxes, fontsize=14)

			except:
				Gc = -1
				Tc = -1
				beta = -1
				rnti = -1
				Psat308 = -1


			pl.grid()
			pl.xlim(250, 500)
			pl.ylim(0)
			pl.tick_params(labelsize=14)
			pl.xlabel('T [mK]', fontsize=15)
			pl.ylabel('Psat [pW]', fontsize=15)

			gtdata = {'Gc': Gc, 'Tc':Tc, 'beta':beta, 'RnTi': rnti, 'Psat308':Psat308}

			print('col ', col_i, ' row ', row_i)
			print('gtdata=', gtdata)

			if fit_beta==False:

				fn = os.path.join(out_Gplot,'sk_G_row%d'%(row) + '_col%d'%col)
				pl.suptitle('mcerow %02d mcecol %02d, detRow %02d detCol %02d Pol-%s'%(row,col,detrow,detcol,detpol), fontsize=15)
				#fnd = os.path.join(out_Gplot,'sk_G_detrow%d'%(detrow) + '_detcol%d'%detcol+'_pol%s'%detpol)
				fnpickle = out_Gplot + Gfn +'.pkl'

			else:

				fn = os.path.join(out_Gplot + 'Fit_beta/', 'sk_G_row%d'%(row) + '_col%d'%col)
				pl.suptitle('mcerow %02d mcecol %02d, detRow %02d detCol %02d Pol-%s'%(row,col,detrow,detcol,detpol), fontsize=15)
				#fnd = os.path.join(out_Gplot + 'Fit_beta/', 'sk_G_detrow%d'%(detrow) + '_detcol%d'%detcol+'_pol%s'%detpol)
				fnpickle = out_Gplot + 'Fit_beta/' + Gfn +'.pkl'

			pl.savefig(fn)
			#pl.savefig(fnd)
			pl.close()

			#shutil.copy2(os.path.realpath(__file__), out_Gplot + (os.path.realpath(__file__).split("/")[-1]).replace(".py",".txt"))
			fn=open(fnpickle,'wb')
			pickle.dump(gtdata, fn)
			fn.close()


	'''
	read data: d = pickle.l/ad(/pen(fnpickle,'r'))
		d[0] --> lcdata
		d[1] --> prdata
		d[0]['457_r4c0']['rnti'/Ohms,'calibFB'/A,'calibBIAS'/A,'rntifit','sc_slope']
		d[1]['457_r4c0']['resistance'/Ohms,'power'/W,'psat'/W,'rnti'/Ohms, responsivity/(1/V)]
		d[2] --> GTdata 'r28c1': {'rnti': 85.20458899931745, 'beta': 1.5068964400413496, 'Gc': 15.159194238818376, 'Tc': 482.51052453878344}
		d[3] --> DetInfos 'r20c1': {'pol': 'A', 'mcecol': 1, 'mcerow': 20, 'detcol': 4, 'detrow': 2}

	'''

if __name__=='__main__':
    main()
