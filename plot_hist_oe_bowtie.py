import os,sys
sys.path.insert(0, "../")
import numpy as np
import cPickle as pickle
import histogram as hist
import shutil
import easyplots as esp
import Lx_ModuleMapping as mm
# inputs
#inpath = '/home/czhang/analysispy/tespy/output/20200122/20200122_BA_BA1_plot/'

module = 'L1'
date = '20220328'
run = '%s_SK_%s'%(date,module)
imod = 0 #(0,1) = 0, (8,9) = 1

inpath = '//scr2/sofiafatigoni/plots/%s/%s/'%(date, run)
fn = '%s_dpdt_rnti.pkl'%run

infile = inpath + fn
d = pickle.load(open(infile,'r'))
# output
outpath = inpath+"hists/"
#outpath = inpath.replace('cryo', 'output')
if not os.path.isdir(outpath):
        os.makedirs(outpath)
shutil.copy2(os.path.realpath(__file__), outpath + (os.path.realpath(__file__).split("/")[-1]).replace(".py",".txt"))

rnti = np.full((24, 41), float('nan'))
rnti_L = np.full((24, 41), float('nan'))
rnti_U = np.full((24, 41), float('nan'))
dpdt = np.full((24, 41), float('nan'))
dpdt_L = np.full((24, 41), float('nan'))
dpdt_U = np.full((24, 41), float('nan'))

for col in range(16):
	for row in range(41):
		try:
			#rnti[col, row] = d[1]['r%dc%d'%(row, col)][2]*1e3
                        rnti[col, row] = d[1][col][row]
                        dpdt[col, row] = d[0][col][row]

                        iM,dc,dr,p = mm.mce2det(col,row)
		except:
			continue



#im = 4
for im in range(imod,imod+1):
    rnti_cut = np.array([x for x in rnti.flatten() if (x>0. and x<800.)])
    hist.plot_1Dhist(rnti_cut.flatten(), outpath, 'module%d_rnti_hist'%im,
    		maintitle=r'$R_n(Ti)$, %s'%module,
    	        xlabel=r'$R_n(Ti)$', 
    		xunit=r'm$\Omega$',
    		binrange=[100., 450.])
    
    fig,ax = esp.presetting(7.4,6,lx="channels",ly="Rn, mOhm")
    dpdt_cut = np.array([x for x in dpdt.flatten() if (x>.0001)])
    print(dpdt.flatten())
    hist.plot_1Dhist(dpdt_cut, outpath, 'module%d_dpdt_hist'%im,
                    maintitle=r'$dP/dT$, %s'%module,
                    xlabel=r'$dP/dT$',
                    xunit=r'pW/K',
                    binrange=[.02, .08])

    fig,ax = esp.presetting(7.4,6,lx="channels",ly="dpdt, pW")

