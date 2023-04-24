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

RUN = 'L9'
rnpsat = 0.040

#DATE = '20190730'
DATE='20230302'
templist = templists[DATE]

def GTcModel(T,G,Tc,beta):
	return G/(beta+1)*(Tc**(beta+1)-T**(beta+1))/Tc**beta

def GTcModel2(T,G,Tc,beta=2):
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
		#filename = 'LC_G_FPU_'+str(temp)+'mK_datamode1_run1'
		filename = 'LC_G_sk_FPU'+str(temp)+'mK_datamode1_source1'
		out_path_main = 'output/%s/'%(DATE)
		if not os.path.isdir(out_path_main):
                        os.makedirs(out_path_main)
		out_path = out_path_main + filename + '/'
		if not os.path.isdir(out_path):
			os.makedirs(out_path)
		out_Gplot = out_path_main + '%s_Gplot/'%(RUN)
		if not os.path.isdir(out_Gplot):
                        os.makedirs(out_Gplot)

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
		# print(nc,nr)
		# sys.exit()
		rows = np.zeros((nc,nr),dtype=np.int)
		cols = np.zeros((nc,nr),dtype=np.int)

		col_i=int(sys.argv[2])
		row_i=int(sys.argv[1])

		Gfn = 'sk_G_'+str(DATE)+'_row'+str(row_i)+'_col'+str(col_i)


		#===================================#
		# Cook LC
		#===================================#
		#for col in range(16):
		for col in [col_i]:
			# for row in range(41):
			for row in [row_i]:
				fitrange = {}

				fitrange["rnti_low"] = 5000.00
				fitrange["rnti_hgh"] = 5500.00
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
	# GTdata = {}
	DetInfos = {}
	#for col in range (16):
	for col in [col_i]:
		#for row in range(41):
		for row in [row_i]:


			dinfos = {}
			#im, detcol,detrow,detpol = ba40_ModuleMapping.mce2det(col,row)
			im, detcol,detrow,detpol = 0, 0, 0, 0
			for itt in range(len(templist)):
				#rnpsatlowlist[itt] = prdata[str(templist[itt])+'_r'+str(row)+'c'+str(col)][4]
				#rnpsathighlist[itt] = prdata[str(templist[itt])+'_r'+str(row)+'c'+str(col)][5]
				psatlist[itt] = prdata[str(templist[itt])+'_r'+str(row)+'c'+str(col)][3]*1.00e12
				rntilist[itt] = prdata[str(templist[itt])+'_r'+str(row)+'c'+str(col)][2]*1.00e3
			#rnpsat = 0.5*(np.max(rnpsatlowlist) + np.min(rnpsathighlist))
			#rnpsat = 0.075
			#for itt in range(len(templist)):
			#	lcIs = lcdata[str(templist[itt])+'_r'+str(row)+'c'+str(col)]
			#	rr, pp, rnti, psat = lc.get_PRs(lcIs[0], lcIs[1], row, col, calib, rnpsat)
			rnti = np.mean(rntilist)
			xdata = np.linspace(templist[0]-1, templist[-1]+1, 100)

			mask=np.asarray(psatlist)>0

			psatlist_good=np.asarray(psatlist)[mask]
			templist_good=np.asarray(templist)[mask]

			#xdata_good = np.linspace(templist_good[0]-1, templist_good[-1]+1, 100)

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


			if psatlist_good[len(psatlist_good)-1]< psatlist_good[0]:
				psatlist_good=psatlist_good[:-1]
				templist_good=templist_good[:-1]


			try:
				popt, pcov = curve_fit(GTcModel, templist_good, psatlist_good)
				#GoT=GTcModel(xdata, *popt)
				##print('xdata:', xdata)
				#G_450=Got(np.where(xdata>449 and xdata<451))
				pl.plot(xdata, GTcModel(xdata, *popt), 'g--',label='fit: Gc=%5.3f, Tc=%5.3f, beta=%5.3f' % tuple(popt))
				#histo parameters

				Gc = popt[0]*1e3
				Tc = popt[1]
				beta = popt[2]
				rnti = rnti
				Psat308 = psatlist[3]
					#G450.append(G_450*1e3)
			except:
				# popt = np.array([-1,-1,-1])
				# pcov = np.array([-1,-1,-1])
				# pass
				#
				# #pl.xlim(200,510)
				# #pl.ylim(0,6)
				# pl.ylabel('P [pW]', fontsize=15)
				# pl.xlabel('T [mK]', fontsize=15)
				# plt.tick_params(labelsize=14)
				# plt.grid()
				Gc = -1
				Tc = -1
				beta = -1
				rnti = -1
				Psat308 = -1

			##print('G:', Gc)
				#if popt[0]>0 and popt[1]>0 and popt[1]<550:
				try:
					pl.text(0.5,0.9,'Gc='+str(round(popt[0]*1e3,1))+' pW/K',transform=ax.transAxes, fontsize=14)
					#print("Gc Try good")
					pl.text(0.5,0.85,'Tc='+str(round(popt[1]))+' mK',transform=ax.transAxes, fontsize=14)
					#print("Tc Try good")
					pl.text(0.5,0.80,'beta='+str(round(popt[2],2)),transform=ax.transAxes, fontsize=14)
					#print("Beta Try good")
					pl.text(0.5,0.75,'RnTi='+str(round(rnti))+' mOhm',transform=ax.transAxes, fontsize=14)
					#print("RnTi Try good")
					#pl.text(0.5,0.7,'Psat308='+str(Psat308)+' pW',transform=ax.transAxes, fontsize=14)
					##print("Psat Try good")
				# 	gtdata['Gc'] = popt[0]*1e3
				# 	gtdata['Tc'] = popt[1]
				# 	gtdata['beta'] = popt[2]
				# 	gtdata['rnti'] = rnti
				# 	gtdata['Psat308'] = psatlist[3]
				# 	dinfos['mcecol'] = col
				# 	dinfos['mcerow'] = row
				# 	dinfos['detcol'] = detcol
				# 	dinfos['detrow'] = detrow
				# 	dinfos['pol'] = detpol
				except Exception as e:
					print(e)

				#
				# 	#print('1')
				# 	gtdata['Gc'] = float('nan')
				# 	gtdata['Tc'] = float('nan')
				# 	gtdata['beta'] = float('nan')
				# 	gtdata['rnti'] = float('nan')
				# 	gtdata['Psat308'] = float('nan')
				# 	dinfos['mcecol'] = float('nan')
				# 	dinfos['mcerow'] = float('nan')
				# 	dinfos['detcol'] = float('nan')
				# 	dinfos['detrow'] = float('nan')
				# 	dinfos['pol'] = float('nan')
					# # else:
					# 	#print('2')
					# 	gtdata['Gc'] = float('nan')
					# 	gtdata['Tc'] = float('nan')
					# 	gtdata['beta'] = float('nan')
					# 	gtdata['rnti'] = float('nan')
					# 	gtdata['Psat308'] = float('nan')
					# 	dinfos['mcecol'] = float('nan')
					# 	dinfos['mcerow'] = float('nan')
					# 	dinfos['detcol'] = float('nan')
					# 	dinfos['detrow'] = float('nan')
					# 	dinfos['pol'] = float('nan')

			#
			# GTdata['r'+str(row)+'c'+str(col)] = gtdata
			# DetInfos['r'+str(row)+'c'+str(col)] = dinfos
			# del gtdata
			# del dinfos

			gtdata = {'Gc': Gc, 'Tc':Tc, 'beta':beta, 'rnti': rnti, 'Psat308':Psat308}

			print('gtdata=', gtdata)

			fn = os.path.join(out_Gplot,'sk_G_row%d'%(row) + '_col%d'%col)
			pl.suptitle('mcerow %02d mcecol %02d, detRow %02d detCol %02d Pol-%s'%(row,col,detrow,detcol,detpol), fontsize=15)
			fnd = os.path.join(out_Gplot,'sk_G_detrow%d'%(detrow) + '_detcol%d'%detcol+'_pol%s'%detpol)

			pl.savefig(fn)
			pl.savefig(fnd)
			pl.close()



			#shutil.copy2(os.path.realpath(__file__), out_Gplot + (os.path.realpath(__file__).split("/")[-1]).replace(".py",".txt"))
			fnpickle = out_Gplot + Gfn +'.pkl'
			fn=open(fnpickle,'wb')
			pickle.dump(gtdata, fn)
			fn.close()






	#makes par histo

	#print('G:',Gc)

	#
	# fig = plt.figure(figsize=(20,6.5), dpi=80)
	# plt.tick_params(labelsize=18)
	# bb= np.linspace(5, 20, 12)
	# Gc_hist=np.asarray(Gc)[np.asarray(Gc)>0]
	# #print('G:',Gc_hist)
	# plt.hist(np.asarray(Gc_hist), bins=bb)  # arguments are passed to np.histogram
	# #print(np.histogram(Gc_hist, bins=bb, density=True))
	# plt.xlabel('pW/K', fontsize=18)
	# plt.title("Gc_hist -"+"Mean="+str(round(np.mean(Gc_hist),3))+"pW/K", fontsize=18)
	# fn_histo = os.path.join(out_Gplot,'Gc_histo')
	# plt.savefig(fn_histo)
	# #print('Gc_mean=',np.mean(np.asarray(Gc_hist)))
	# #print('Gc_std=',np.std(np.asarray(Gc_hist)))
	#
	#
	# fig = plt.figure(figsize=(20,6.5), dpi=80)
	# plt.tick_params(labelsize=18)
	# bb= np.linspace(480, 510, 13)
	# Tc_hist=np.asarray(Tc)[np.asarray(Tc)>0]
	# #Tc_hist=Tc_hist[Tc_hist<550]
	# plt.hist(np.asarray(Tc_hist), bins=bb)  # arguments are passed to np.histogram
	# #print(np.histogram(Tc_hist, bins=bb, density=True))
	# plt.xlabel('mK',fontsize=18)
	# plt.title("Tc_hist -"+" Mean="+str(round(np.mean(Tc_hist),3))+"mK", fontsize=18)
	# fn_histo = os.path.join(out_Gplot,'Tc_histo')
	# plt.savefig(fn_histo)
	# #print('Tc_mean=',np.mean(np.asarray(Tc_hist)))
	# #print('Tc_std=',np.std(np.asarray(Tc_hist)))
	#
	#
	# fig = plt.figure(figsize=(20,6.5), dpi=80)
	# plt.tick_params(labelsize=18)
	# bb= np.linspace(1, 2.5, 10)
	# beta=np.asarray(beta)[np.asarray(beta)>0]
	# plt.hist(np.asarray(beta), bins=bb)  # arguments are passed to np.histogram
	# #print(np.histogram(beta, bins=bb, density=True))
	# plt.title("Beta_hist -"+" Mean="+str(round(np.mean(beta),3)),fontsize=18)
	# fn_histo = os.path.join(out_Gplot,'beta_histo')
	# plt.savefig(fn_histo)
	#
	# #print('beta_mean=',np.mean(np.asarray(beta)))
	# #print('beta_std=',np.std(np.asarray(beta)))
	#
	#
	# fig = plt.figure(figsize=(20,6.5), dpi=80)
	# plt.tick_params(labelsize=18)
	# bb= np.linspace(0, 150, 10)
	# RnTi=np.asarray(RnTi)[np.asarray(RnTi)>0]
	# RnTi=RnTi[np.asarray(RnTi)<500]
	# plt.hist(np.asarray(RnTi), bins=bb)  # arguments are passed to np.histogram
	# #print(np.histogram(RnTi, bins=bb, density=True))
	# plt.xlabel('mOhm',fontsize=18)
	# plt.title("RnTi -"+" Mean="+str(round(np.mean(RnTi),2))+"mOhm",fontsize=18)
	# fn_histo = os.path.join(out_Gplot,'RnTi_histo')
	# plt.savefig(fn_histo)
	# #print('RnTi_mean=',np.mean(np.asarray(RnTi)))
	# #print('RnTi_std=',np.std(np.asarray(RnTi)))
	#
	#
	# fig = plt.figure(figsize=(20,6.5), dpi=80)
	# plt.tick_params(labelsize=18)
	# bb= np.linspace(0, 5, 15)
	# Psat308=np.asarray(Psat308)[np.asarray(Psat308)>0]
	# plt.hist(np.asarray(Psat308), bins=bb)  # arguments are passed to np.histogram
	# #print(np.histogram(Psat308, bins=bb, density=True))
	# plt.xlabel('pW',fontsize=18)
	# plt.title("Psat308 -"+" Mean="+str(round(np.mean(Psat308),3))+"pW",fontsize=18)
	# fn_histo = os.path.join(out_Gplot,'Psat308_histo')
	# plt.savefig(fn_histo)
	# #print('Psat308_mean=',np.mean(np.asarray(Psat308)))
	# #print('Psat308_std=',np.std(np.asarray(Psat308)))
	# plt.close()
	#



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
