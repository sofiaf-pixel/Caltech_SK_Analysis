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
#rnpsat = 0.040

DATE='20230302'
templist = templists[DATE]



def make_histo(out_Gplot, Gfn, nrows, ncols):

	Gc=[]
	Tc=[]
	beta=[]
	RnTi=[]
	Psat308=[]

	for row in range (nrows):
		for col in range (ncols):
			try:

				fnpickle = out_Gplot + Gfn +'_row'+str(row)+'_col'+str(col)+'.pkl'

				fn=open(fnpickle,'rb')
				gtdata = pickle.load(fn)
				fn.close()

				print(gtdata)

				Gc.append(gtdata['Gc'])
				Tc.append(gtdata['Tc'])
				beta.append(gtdata['beta'])
				RnTi.append(gtdata['rnti'])
				Psat308.append(gtdata['Psat308'])

			except Exception as e:

				print(e)


	print ('Gc = ', Gc)
	print('Tc =', Tc)
	print('beta = ', beta)



	fig = plt.figure(figsize=(20,6.5), dpi=80)
	plt.tick_params(labelsize=18)
	bb= np.linspace(60, 140, 12)

	Gc_array=np.asarray(Gc)
	#print('Gc_array=', Gc_array)
	Gc_mask=np.where(Gc_array>0)[0]
	Gc_hist=Gc_array[Gc_mask]
	#print('G:',Gc_hist)
	plt.hist(np.asarray(Gc_hist), bins=bb, facecolor='r', edgecolor='k')  # arguments are passed to np.histogram
	#print(np.histogram(Gc_hist, bins=bb, density=True))
	plt.xlabel('pW/K', fontsize=18)
	plt.title("Gc_hist -"+"Median="+str(round(np.nanmedian(Gc_hist),3))+"pW/K", fontsize=18)
	fn_histo = os.path.join(out_Gplot,'Gc_histo')
	plt.savefig(fn_histo)
	#print('Gc_Median=',np.mean(np.asarray(Gc_hist)))
	#print('Gc_std=',np.std(np.asarray(Gc_hist)))


	fig = plt.figure(figsize=(20,6.5), dpi=80)
	plt.tick_params(labelsize=18)
	bb= np.linspace(480, 540, 13)
	Tc_hist=np.asarray(Tc)[np.asarray(Tc)>0]
	#Tc_hist=Tc_hist[Tc_hist<550]
	plt.hist(np.asarray(Tc_hist), bins=bb, facecolor='r', edgecolor='k')  # arguments are passed to np.histogram
	#print(np.histogram(Tc_hist, bins=bb, density=True))
	plt.xlabel('mK',fontsize=18)
	plt.title("Tc_hist -"+" Median="+str(round(np.nanmedian(Tc_hist),3))+"mK", fontsize=18)
	fn_histo = os.path.join(out_Gplot,'Tc_histo')
	plt.savefig(fn_histo)
	#print('Tc_mean=',np.mean(np.asarray(Tc_hist)))
	#print('Tc_std=',np.std(np.asarray(Tc_hist)))


	fig = plt.figure(figsize=(20,6.5), dpi=80)
	plt.tick_params(labelsize=18)
	bb= np.linspace(0, 3, 10)
	beta=np.asarray(beta)[np.asarray(beta)>0]
	plt.hist(np.asarray(beta), bins=bb, facecolor='r', edgecolor='k')  # arguments are passed to np.histogram
	#print(np.histogram(beta, bins=bb, density=True))
	plt.title("Beta_hist -"+" Median="+str(round(np.nanmedian(beta),3)),fontsize=18)
	fn_histo = os.path.join(out_Gplot,'beta_histo')
	plt.savefig(fn_histo)

	#print('beta_mean=',np.mean(np.asarray(beta)))
	#print('beta_std=',np.std(np.asarray(beta)))

	fig = plt.figure(figsize=(20,6.5), dpi=80)
	plt.tick_params(labelsize=18)
	bb= np.linspace(40, 120, 10)
	RnTi=np.asarray(RnTi)[np.asarray(RnTi)>0]
	RnTi=RnTi[np.asarray(RnTi)<500]
	plt.hist(np.asarray(RnTi), bins=bb, facecolor='r', edgecolor='k')  # arguments are passed to np.histogram
	#print(np.histogram(RnTi, bins=bb, density=True))
	plt.xlabel('mOhm',fontsize=18)
	plt.title("RnTi -"+" Median="+str(round(np.nanmedian(RnTi),2))+"mOhm",fontsize=18)
	fn_histo = os.path.join(out_Gplot,'RnTi_histo')
	plt.savefig(fn_histo)
	#print('RnTi_mean=',np.mean(np.asarray(RnTi)))
	#print('RnTi_std=',np.std(np.asarray(RnTi)))


	fig = plt.figure(figsize=(20,6.5), dpi=80)
	plt.tick_params(labelsize=18)
	bb= np.linspace(10, 25, 15)
	Psat308=np.asarray(Psat308)[np.asarray(Psat308)>0]
	plt.hist(np.asarray(Psat308), bins=bb, facecolor='r', edgecolor='k')  # arguments are passed to np.histogram
	#print(np.histogram(Psat308, bins=bb, density=True))
	plt.xlabel('pW',fontsize=18)
	plt.title("Psat308 -"+" Median="+str(round(np.nanmedian(Psat308),3))+"pW",fontsize=18)
	fn_histo = os.path.join(out_Gplot,'Psat308_histo')
	plt.savefig(fn_histo)
	#print('Psat308_mean=',np.mean(np.asarray(Psat308)))
	#print('Psat308_std=',np.std(np.asarray(Psat308)))
	plt.show()




def main():

	out_path_main = 'output/%s/'%(DATE)
	out_Gplot= out_path_main + '%s_Gplot/'%(RUN)
	Gfn = 'sk_G_'+str(DATE)

	make_histo(out_Gplot, Gfn, 41, 16)





if __name__=='__main__':
    main()
