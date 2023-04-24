#Example of Script that makes Tile plots
#(it makes real tile plots for light det and fake tile plots for dark dets)
#This was used for BA1 G Measurements tile plots
#
#SF - 2019/07

import numpy as np
import sys,os
import pickle as pickle
import triangle_mapping_150_SF as tm

# load data
DATE='20230302'

RUN = 'L9'

out_path_main = 'output/%s/'%(DATE)
out_Gplot= out_path_main + '%s_Gplot/'%(RUN)
Gfn = 'sk_G_'+str(DATE)

# fnpickle = "output/"+DATE+"/ba40Gscript_Gplot/sk_G_"+DATE+".pkl"
# d = pickle.load(open(fnpickle,'rb'))
# Gdata = d[-2]
#

Gc = np.full((16,41), float("nan"))
biaslist = np.linspace(0., 2000., 100)

for col in range(16):
	for row in range(41):
		try:
			fnpickle = out_Gplot + Gfn +'_row'+str(row)+'_col'+str(col)+'.pkl'

			fn=open(fnpickle,'rb')
			gtdata = pickle.load(fn)
			fn.close()

			print(gtdata)

			Gc[col,row] =   gtdata['Gc']

		except Exception as e:
			print(e)

# prepare data
Gca = tm.prepare_data(Gc, mapping="BA2")

#For dark detectors
#Gca_dark, Gca_dark_short = tm.prepare_data_dark(Gc, mapping="SK")

#prepare triangle
lld = 5
lld_dark=2
[x,y,triangles] = tm.prepare_triangles(lld)
[x_dark,y_dark,triangles_dark] = tm.prepare_triangles(lld_dark)

#mask
'''mask_wire = [ 0,float('Inf'),float('Inf'),0,0,  #col1A
              0,0,0,float('Inf'),0,     #2A
              0,0,0,0,0,        #3A
              0,0,0,0,float('Inf'),     #4A
              0,0,0,float('Inf'),float('Inf'),  #5A
              0,0,0,float('Inf'),float('Inf'),  #col1B
              0,0,0,0,0,        #2B
              0,float('Inf'),0,0,0,     #3B
              float('Inf'),0,0,0,0,     #4B
              0,0,float('Inf'),0,0]     #5B'''


mask_wire = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,	#1A
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,	#1A
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]



mask_wire_dark = [ 0,0,  #1A
		   0,0,  #2A
		   0,0,  #1B
		   0,0]  #2B

# make plots
outpath = "output/"+DATE+"/TilePlots/"
#tm.plot_a_map(mcecooda,triangles,x,y,mask_wire,0, 0.0, 1, 'mce col',' - ',outpath)


# betamask=beta[np.logical_not(np.isnan(beta))]
Gcamask=Gc[np.logical_not(np.isnan(Gc))]
#
# Psat308amask=Psat308[np.logical_not(np.isnan(Psat308))]
# Psat308amask=Psat308amask[Psat308amask>0]
#
# Tcamask=Tc[np.logical_not(np.isnan(Tc))]
# Tcamask=Tcamask[Tcamask>450]
#
# RnTiamask=RnTi[np.logical_not(np.isnan(RnTi))]
# RnTiamask=RnTiamask[RnTiamask<500]



tm.plot_a_map(Gca,triangles, x, y, mask_wire, 0, 5, 30, 'Gc_L1','pW/K',outpath, 'Gc_L1_Tile', ifshow=True, cmap="jet") #changer 2 and 3 number, lower and upper limit

# tm.plot_a_map(Tca,triangles,x,y,mask_wire,0, 490, 498, 'Tc_T6','mK',outpath, cmap="jet")
# tm.plot_a_map(beta_a,triangles,x,y,mask_wire,0, 1.5, 1.9, 'beta_T6','',outpath, cmap="jet")
# tm.plot_a_map(RnTia,triangles,x,y,mask_wire,0, 50, 150, 'RnTi_T6','mOhm',outpath, cmap="jet")
# tm.plot_a_map(Psat308a,triangles,x,y,mask_wire,0, 1, 4.5, 'Psat308_T6','pW',outpath, cmap="jet")
#
#
# tm.plot_a_map_dark(Gca_dark_short,triangles_dark,x_dark,y_dark,mask_wire_dark,0, np.min(np.asarray(Gca_dark_short)), np.max(np.asarray(Gca_dark_short)), 'Gc_T6_Dark','pW/K',outpath, mapping='SK', cmap="jet")
# tm.plot_a_map_dark(Tca_dark_short,triangles_dark,x_dark,y_dark,mask_wire_dark,0, np.min(np.asarray(Tca_dark_short)), np.max(np.asarray(Tca_dark_short)), 'Tc_T6_Dark','mK',outpath, mapping='SK',cmap="jet")
# tm.plot_a_map_dark(beta_a_dark_short,triangles_dark,x_dark,y_dark,mask_wire_dark,0, np.min(np.asarray(beta_a_dark_short)), np.max(np.asarray(beta_a_dark_short)), 'beta_T6_Dark','',outpath, mapping='SK', cmap="jet")
# tm.plot_a_map_dark(RnTia_dark_short,triangles_dark,x_dark,y_dark,mask_wire_dark,0, np.min(np.asarray(RnTia_dark_short)), np.max(np.asarray(RnTia_dark_short)), 'RnTi_T6_Dark','mOhm',outpath, mapping='SK', cmap="jet")
# tm.plot_a_map_dark(Psat308a_dark_short,triangles_dark,x_dark,y_dark,mask_wire_dark,0, np.min(np.asarray(Psat308a_dark_short)), np.max(np.asarray(Psat308a_dark_short)), 'Psat308_T6_Dark','pW',outpath, mapping='SK', cmap="jet")
