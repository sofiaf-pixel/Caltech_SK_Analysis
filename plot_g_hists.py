import os,sys
sys.path.insert(0, "../")
import shutil
import numpy as np
import cPickle as pickle
import histogram as hist
import Lx_ModuleMapping as minfo
import triangle_mapping_Lx as tri
# inputs
date = '20221017'
module = 'SK_L4_150GHz'
im = 0

inpath  = '//scr2/sofiafatigoni/plots/%s/%s_SK_%s_Gplot/'%(date,date,module)
fn='SK_G_%s_%s.pkl'%(module,date)

infile = inpath + fn
d = pickle.load(open(infile,'r'))

# output
outpath = inpath+"hists/"
if not os.path.isdir(outpath):
        os.makedirs(outpath)
shutil.copy2(os.path.realpath(__file__), outpath + (os.path.realpath(__file__).split("/")[-1]).replace(".py",".txt"))

gc = np.full((32, 41), float('nan'))
g450_low = np.full((32, 41), float('nan'))
g450_high = np.full((32, 41), float('nan'))
g450 = np.full((32, 41), float('nan'))
rnti = np.full((32, 41), float('nan'))
beta = np.full((32, 41), float('nan'))
tc = np.full((32, 41), float('nan'))
psat = np.full((32, 41), float('nan'))
t = 270.
for col in range(32):
	for row in range(41):
                try:
			rnti[col, row] = d[2]['r%dc%d'%(row, col)]['rnti']
			gc[col, row] = d[2]['r%dc%d'%(row, col)]['Gc']
			tc[col, row] = d[2]['r%dc%d'%(row, col)]['Tc']
			beta[col, row] = d[2]['r%dc%d'%(row, col)]['beta']
			g450[col, row] = gc[col, row]*(450./tc[col, row])**beta[col, row]
			psat[col, row] = gc[col, row]/(2.+1.)*(tc[col, row]**(2+1)-t**(2+1))/tc[col, row]**2*1e-3
                        blah, detcol,detrow,detpol= minfo.mce2det(col%16,row)
                        
		except:
			continue

hist.plot_1Dhist((rnti[im*16:im*16+16,:]), outpath, '%s_rnti_hist'%module,
		maintitle=r'$R_n(Ti)$, %s'%module,
		xlabel=r'$R_n(Ti)$', 
		xunit=r'm$\Omega$',
		binrange=[40., 250.])

tri.triange_mapping_single(rnti[im*16:im*16+16,:], itile = 0,
                vmin = 40.,
                vmax = 250.,
                vari = r'$R_n(Ti)$', 
                unit = r'm$\Omega$',
                outpath = outpath + '%s_rnti_tile.png'%module,
                mask_wire = None,
                mapping='SK',
                maintitle= r'$R_n(Ti)$, %s'%module,
                cmap = 'jet')

hist.plot_1Dhist((beta[im*16:im*16+16,:]), outpath, '%s_beta_hist'%module,
		maintitle=r'$\beta$, %s'%module,
		xlabel=r'$\beta$', 
		xunit=r'',
		binrange=[1.2, 2.4])

tri.triange_mapping_single(beta[im*16:im*16+16,:], itile = 0,
                vmin = 1.2,
                vmax = 2.4,
                vari = r'$\beta$', 
                unit = r'',
                outpath = outpath + '%s_beta_tile.png'%module,
                mask_wire = None,
                mapping='SK',
                maintitle= r'$\beta$, %s'%module,
                cmap = 'jet')

hist.plot_1Dhist((gc[im*16:im*16+16, :]), outpath, '%s_gc_hist'%module,
                maintitle=r'$G_c$, %s'%module,
                xlabel=r'$G_c$',
                xunit=r'pW/K',
                binrange=[30., 80.])

tri.triange_mapping_single(gc[im*16:im*16+16,:], itile = 0,
                vmin = 30.,
                vmax = 80.,
                vari = r'$G_c$', 
                unit = r'pW/K',
                outpath = outpath + '%s_gc_tile.png'%module,
                mask_wire = None,
                mapping='SK',
                maintitle= r'$G_c$, %s'%module,
                cmap = 'jet')


hist.plot_1Dhist((g450[im*16:im*16+16, :]), outpath, '%s_g450_hist'%module,
                maintitle=r'$G_{450}$, %s'%module,
                xlabel=r'$G_{450}$',
                xunit=r'pW/K',
                binrange=[30., 55.])

tri.triange_mapping_single(g450[im*16:im*16+16,:], itile = 0,
                vmin = 30.,
                vmax = 55.,
                vari = r'$G_450$', 
                unit = r'pW/K',
                outpath = outpath + '%s_g450_tile.png'%module,
                mask_wire = None,
                mapping='SK',
                maintitle= r'$G_450$, %s'%module,
                cmap = 'jet')


hist.plot_1Dhist((tc[im*16:im*16+16, :]), outpath, '%s_tc_hist'%module,
                maintitle=r'$T_c$, %s'%module,
                xlabel=r'$T_c$',
                xunit=r'mK',
                binrange=[400., 570.])

tri.triange_mapping_single(tc[im*16:im*16+16,:], itile = 0,
                vmin = 400.,
                vmax = 570.,
                vari = r'$T_c$', 
                unit = r'mK',
                outpath = outpath + '%s_tc_tile.png'%module,
                mask_wire = None,
                mapping='SK',
                maintitle= r'$T_c$, %s'%module,
                cmap = 'jet')

hist.plot_1Dhist((psat[im*16:im*16+16, :]), outpath, '%s_psat270_hist'%module,
                maintitle=r'$P_{sat}(270 mK)$, %s'%module,
                xlabel=r'$P_{sat}$',
                xunit=r'pW',
                binrange=[2.0, 10.0])

tri.triange_mapping_single(psat[im*16:im*16+16,:], itile = 0,
                vmin = 2.,
                vmax = 10.,
                vari = r'$P_{sat}(270 mK)$', 
                unit = r'pW',
                outpath = outpath + '%s_psat_tile.png'%module,
                mask_wire = None,
                mapping='SK',
                maintitle= r'$P_{sat}(270 mK)$, %s'%module,
                cmap = 'jet')
