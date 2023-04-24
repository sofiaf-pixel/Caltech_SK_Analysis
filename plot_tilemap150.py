import os,sys
import numpy as np
import pickle
import Lx_ModuleMapping as minfo
import triangle_mapping_Lx as tri
import histogram as hist
import shutil
# inputs
# module = 'L1'
# date = '20220914'
# run = '%s_BA_L1'%date
im = 1
# inpath = '//scr2/sofiafatigoni/plots/%s/%s/'%(date, run)
# fn = run+'_dpdt_rnti.pkl'


DATE='20230302'
RUN = 'L9'

out_path_main = 'output/%s/'%(DATE)
out_Gplot= out_path_main + '%s_Gplot/'%(RUN)

inpath="../Ganalysis/output/"+DATE+"/L9_Gplot/"


cmap = 'jet'


# output
outpath = inpath+"tilemaps/"
print(outpath)
if not os.path.isdir(outpath):
        os.makedirs(outpath)
shutil.copy2(os.path.realpath(__file__), outpath + (os.path.realpath(__file__).split("/")[-1]).replace(".py",".txt"))

rnti = np.array([[np.nan for i2 in range(41)] for i1 in range(24)])
psat = np.array([[np.nan for i2 in range(41)] for i1 in range(24)])
Gc = np.array([[np.nan for i2 in range(41)] for i1 in range(24)])
Tc = np.array([[np.nan for i2 in range(41)] for i1 in range(24)])
beta = np.array([[np.nan for i2 in range(41)] for i1 in range(24)])

for col in range(16):
    for row in range(41):
        try:
            fn = "sk_G_"+DATE+"_row"+str(row)+"_col"+str(col)+".pkl"
            infile = inpath+fn
            f = open(infile,'rb')
            d = pickle.load(f)
            f.close()
            print(d)
            rnti[col, row] = d['rnti']
            psat[col, row] = d['Psat308']
            Gc[col, row] = d['Gc']
            Tc[col, row] = d['Tc']
            beta[col, row] = d['beta']
        except:
            continue

# hist.plot_1Dhist((psat), outpath, '%s_Psat308_hist'%RUN,
# 		maintitle=r'Psat308, %s'%RUN,
# 		xlabel=r'Psat308$',
# 		xunit=r'pW/K$',
# 		binrange=[9, 12])

tri.triange_mapping_single(psat, itile = 0,
                        vmin = 8,
                        vmax = 12,
                        vari = r'Psat308',
                        unit = r'pW',
                        outpath = outpath + '%s_Psat308.png'%RUN,
                        mask_wire = None,
                        mapping='SK',
                        maintitle= r'Psat308, %s'%RUN,
                        cmap = cmap)

tri.triange_mapping_single(Tc, itile = 0,
                        vmin = 480,
                        vmax = 500,
                        vari = r'Tc',
                        unit = r'mK',
                        outpath = outpath + '%s_Tc.png'%RUN,
                        mask_wire = None,
                        mapping='SK',
                        maintitle= r'Tc, %s'%RUN,
                        cmap = cmap)

tri.triange_mapping_single(Gc, itile = 0,
                        vmin =60,
                        vmax = 120,
                        vari = r'Gc',
                        unit = r'pW/K',
                        outpath = outpath + '%s_Gc.png'%RUN,
                        mask_wire = None,
                        mapping='SK',
                        maintitle= r'Gc, %s'%RUN,
                        cmap = cmap)

tri.triange_mapping_single(beta, itile = 0,
                        vmin = 1,
                        vmax = 3,
                        vari = r'beta',
                        unit = r'',
                        outpath = outpath + '%s_beta.png'%RUN,
                        mask_wire = None,
                        mapping='SK',
                        maintitle= r'beta, %s'%RUN,
                        cmap = cmap)
# hist.plot_1Dhist((rnti), outpath, '%s_rnti_hist'%RUN,
# 		maintitle=r'$R_n(Ti)$, %s'%RUN,
# 		xlabel=r'$R_n(Ti)$',
# 		xunit=r'm$\Omega$',
# 		binrange=[45., 85.])

tri.triange_mapping_single(rnti, itile = 0,
			vmin = 85,
			vmax = 45.,
			vari = r'$R_n$',
			unit = r'm$\Omega$',
			outpath = outpath + '%s_rnti.png'%RUN,
			mask_wire = None,
			mapping='SK',
			maintitle= r'$R_n(Ti)$, %s'%RUN,
			cmap = cmap)
