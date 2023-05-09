import os,sys
import numpy as np
import pickle
import Lx_ModuleMapping as minfo
import triangle_mapping_Lx as tri
import histogram as hist
import shutil

RUN = 'L1'
DATE='20230302'

fit_beta=False

main_path = 'output/%s/'%(DATE)

if fit_beta==False:
    inpath = main_path +  '%s_Gplot/'%(RUN)
else:
    inpath = main_path +  '%s_Gplot/'%(RUN) + 'Fit_beta/'

cmap = 'jet'


# output
outpath = inpath + "Tilemaps/"
if not os.path.isdir(outpath):
        os.makedirs(outpath)
print(outpath)
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
            #print(d)
            rnti[col, row] = d['RnTi']
            psat[col, row] = d['Psat308']
            Gc[col, row] = d['Gc']
            Tc[col, row] = d['Tc']
            beta[col, row] = d['beta']
        except Exception as e:
            print(e)
            continue

# hist.plot_1Dhist((psat), outpath, '%s_Psat308_hist'%RUN,
# 		maintitle=r'Psat308, %s'%RUN,
# 		xlabel=r'Psat308$',
# 		xunit=r'pW/K$',
# 		binrange=[9, 12])

#print('beta=',beta)
#print('Gc=', Gc)


tri.triange_mapping_single(psat, itile = 0,
                        vmin = 8,
                        vmax = 13,
                        vari = r'Psat308',
                        unit = r'pW',
                        outpath = outpath + '%s_Psat308_tilemap.png'%RUN,
                        mask_wire = None,
                        mapping='SK',
                        maintitle= r'Psat308, %s'%RUN,
                        cmap = cmap)

tri.triange_mapping_single(Tc, itile = 0,
                        vmin = 480,
                        vmax = 500,
                        vari = r'Tc',
                        unit = r'mK',
                        outpath = outpath + '%s_Tc_tilemap.png'%RUN,
                        mask_wire = None,
                        mapping='SK',
                        maintitle= r'Tc, %s'%RUN,
                        cmap = cmap)

tri.triange_mapping_single(Gc, itile = 0,
                        vmin =70,
                        vmax = 110,
                        vari = r'Gc',
                        unit = r'pW/K',
                        outpath = outpath + '%s_Gc_tilemap.png'%RUN,
                        mask_wire = None,
                        mapping='SK',
                        maintitle= r'Gc, %s'%RUN,
                        cmap = cmap)

tri.triange_mapping_single(beta, itile = 0,
                        vmin = 1.8,
                        vmax = 2.2,
                        vari = r'beta',
                        unit = r'',
                        outpath = outpath + '%s_beta_tilemap.png'%RUN,
                        mask_wire = None,
                        mapping='SK',
                        maintitle= r'beta, %s'%RUN,
                        cmap = cmap)

tri.triange_mapping_single(rnti, itile = 0,
			            vmin = 40,
			            vmax = 85,
			            vari = r'$R_n$',
			            unit = r'm$\Omega$',
			            outpath = outpath + '%s_RnTi_tilemap.png'%RUN,
			            mask_wire = None,
			            mapping='SK',
			            maintitle= r'$R_n(Ti)$, %s'%RUN,
			            cmap = cmap)
