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

RUN = 'L1'
DATE='20230302'

fit_beta=False

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

				#print(gtdata)

                Gc.append(gtdata['Gc'])
                Tc.append(gtdata['Tc'])
                beta.append(gtdata['beta'])
                RnTi.append(gtdata['RnTi'])
                Psat308.append(gtdata['Psat308'])

            except Exception as e:
                print(e)

    fig = plt.figure(figsize=(20,6.5), dpi=80)
    plt.tick_params(labelsize=18)
    bb= np.linspace(60, 120, 11)
    Gc_hist=np.asarray(Gc)[np.asarray(Gc)>60]
    Gc_hist=np.asarray(Gc_hist)[np.asarray(Gc_hist)<120]

    #Gc_array=np.asarray(Gc)
	#print('Gc_array=', Gc_array)
    #Gc_mask=np.where(Gc_array>0)[0]
    #Gc_hist=Gc_array[Gc_mask]
	#print('G:',Gc_hist)

    plt.hist(np.asarray(Gc_hist), bins=bb, facecolor='r', edgecolor='k')  # arguments are passed to np.histogram
	#print(np.histogram(Gc_hist, bins=bb, density=True))
    plt.xlabel('pW/K', fontsize=18)
    plt.title("Gc -"+" Median = "+str(round(np.nanmedian(Gc_hist),3))+" pW/K", fontsize=18)
    fn_histo = os.path.join(out_Gplot + 'Histograms/', '%s_Gc_histo'%RUN)
    print(fn_histo)
    plt.savefig(fn_histo)
	#print('Gc_Median=',np.mean(np.asarray(Gc_hist)))
	#print('Gc_std=',np.std(np.asarray(Gc_hist)))


    fig = plt.figure(figsize=(20,6.5), dpi=80)
    plt.tick_params(labelsize=18)
    bb= np.linspace(480, 510, 11)
    #plt.xticks([480,485,490,495,500,505,510])
    Tc_hist=np.asarray(Tc)[np.asarray(Tc)>480]
    Tc_hist=np.asarray(Tc_hist)[np.asarray(Tc_hist)<510]
	#Tc_hist=Tc_hist[Tc_hist<550]
    plt.hist(np.asarray(Tc_hist), bins=bb, facecolor='r', edgecolor='k')  # arguments are passed to np.histogram
	#print(np.histogram(Tc_hist, bins=bb, density=True))
    plt.xlabel('mK',fontsize=18)
    plt.title("Tc -"+" Median = "+str(round(np.nanmedian(Tc_hist),3))+" mK", fontsize=18)
    fn_histo = os.path.join(out_Gplot + 'Histograms/', '%s_Tc_histo'%RUN)
    print(fn_histo)
    plt.savefig(fn_histo)
	#print('Tc_mean=',np.mean(np.asarray(Tc_hist)))
	#print('Tc_std=',np.std(np.asarray(Tc_hist)))


    fig = plt.figure(figsize=(20,6.5), dpi=80)
    plt.tick_params(labelsize=18)
    bb= np.linspace(1, 3, 11)
    plt.xticks([1,1.5,2,2.5,3])
    beta=np.asarray(beta)[np.asarray(beta)>1]
    beta=np.asarray(beta)[np.asarray(beta)<3]
    plt.hist(np.asarray(beta), bins=bb, facecolor='r', edgecolor='k')  # arguments are passed to np.histogram
	#print(np.histogram(beta, bins=bb, density=True))
    plt.title("Beta -"+" Median = "+str(round(np.nanmedian(beta),3)),fontsize=18)
    fn_histo = os.path.join(out_Gplot + 'Histograms/', '%s_beta_histo'%RUN)
    print(fn_histo)
    plt.savefig(fn_histo)
	#print('beta_mean=',np.mean(np.asarray(beta)))
	#print('beta_std=',np.std(np.asarray(beta)))


    fig = plt.figure(figsize=(20,6.5), dpi=80)
    plt.tick_params(labelsize=18)
    bb= np.linspace(45, 95, 11)
    RnTi=np.asarray(RnTi)[np.asarray(RnTi)>45]
    RnTi=RnTi[np.asarray(RnTi)<95]
    plt.hist(np.asarray(RnTi), bins=bb, facecolor='r', edgecolor='k')  # arguments are passed to np.histogram
	#print(np.histogram(RnTi, bins=bb, density=True))
    plt.xlabel('mOhm',fontsize=18)
    plt.title("RnTi -"+" Median = "+str(round(np.nanmedian(RnTi),2))+" mOhm",fontsize=18)
    fn_histo = os.path.join(out_Gplot + 'Histograms/', '%s_RnTi_histo'%RUN)
    print(fn_histo)
    plt.savefig(fn_histo)
	#print('RnTi_mean=',np.mean(np.asarray(RnTi)))
	#print('RnTi_std=',np.std(np.asarray(RnTi)))


    fig = plt.figure(figsize=(20,6.5), dpi=80)
    plt.tick_params(labelsize=18)
    bb= np.linspace(5, 15, 11)
    Psat308=np.asarray(Psat308)[np.asarray(Psat308)>5]
    Psat308=np.asarray(Psat308)[np.asarray(Psat308)<15]
    plt.hist(np.asarray(Psat308), bins=bb, facecolor='r', edgecolor='k')  # arguments are passed to np.histogram
	#print(np.histogram(Psat308, bins=bb, density=True))
    plt.xlabel('pW',fontsize=18)
    plt.title("Psat308 -"+" Median = "+str(round(np.nanmedian(Psat308),3))+" pW",fontsize=18)
    fn_histo = os.path.join(out_Gplot + 'Histograms/', '%s_Psat308_histo'%RUN)
    print(fn_histo)
    plt.savefig(fn_histo)
	#print('Psat308_mean=',np.mean(np.asarray(Psat308)))
	#print('Psat308_std=',np.std(np.asarray(Psat308)))
    #plt.show()




def main():

    out_path_main = 'output/%s/'%(DATE)

    if fit_beta==False:
        out_Gplot= out_path_main + '%s_Gplot/'%(RUN)
        Gfn = 'sk_G_'+str(DATE)

    else:
        out_Gplot= out_path_main + '%s_Gplot/'%(RUN) + 'Fit_beta/'
        Gfn = 'sk_G_'+str(DATE)

    if not os.path.isdir(out_Gplot + 'Histograms/'):
        os.makedirs(out_Gplot + 'Histograms/')

    make_histo(out_Gplot, Gfn, 41, 16)





if __name__=='__main__':
    main()
