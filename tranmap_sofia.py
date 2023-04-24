import numpy as np
import sys,os
import pickle as pickle
#sys.path.insert(0, "/home/cheng/analysis/load_curve/sk_dark/scripts_house")
import triangle_mapping as tm

# load data
DATE='20230302'
RUN = 'L9'

out_path_main = 'output/%s/'%(DATE)
out_Gplot= out_path_main + '%s_Gplot/'%(RUN)

fnpickle = "output/"+DATE+"/ba40Gscript_Gplot/sk_G_"+DATE+".pkl"

d = pickle.load(open(fnpickle,'rb'))

Gdata = d[-2]

Gc = np.full((2,33), float("nan"))
Tc = np.full((2,33), float("nan"))
beta = np.full((2,33), float("nan"))
RnTi = np.full((2,33), float("nan"))
Psat308 = np.full((2,33), float("nan"))

biaslist = np.linspace(0., 2000., 100)

for col in [0,1]:
	for row in range(33):
		Gc[col,row] =   Gdata["r%dc%d"%(row, col)]["Gc"]
		Tc[col,row] =   Gdata["r%dc%d"%(row, col)]["Tc"]
		beta[col,row] = Gdata["r%dc%d"%(row, col)]["beta"]
		RnTi[col,row] = Gdata["r%dc%d"%(row, col)]["rnti"]
		Psat308[col,row] = Gdata["r%dc%d"%(row, col)]["Psat308"]

print(Psat308)

#col_dark=[0,1]
#row_dark=[8,9,10,11]

#for col in col_dark:
        #for row in row_dark:
                #Gc_dark[col,row] =   Gdata["r%dc%d"%(row, col)]["Gc"]
                #Tc_dark[col,row] =   Gdata["r%dc%d"%(row, col)]["Tc"]
                #beta_dark[col,row] = Gdata["r%dc%d"%(row, col)]["beta"]
                #RnTi_dark[col,row] = Gdata["r%dc%d"%(row, col)]["rnti"]
                #Psat308_dark[col,row] = Gdata["r%dc%d"%(row, col)]["Psat308"]





# prepare data
Gca = tm.prepare_data(Gc, mapping="SK")
Tca = tm.prepare_data(Tc, mapping="SK")
beta_a = tm.prepare_data(beta, mapping="SK")
RnTia = tm.prepare_data(RnTi, mapping="SK")
Psat308a = tm.prepare_data(Psat308, mapping="SK")


#For dark detectors

Gca_dark, Gca_dark_short = tm.prepare_data_dark(Gc, mapping="SK")
Tca_dark, Tca_dark_short = tm.prepare_data_dark(Tc, mapping="SK")
beta_a_dark, beta_a_dark_short = tm.prepare_data_dark(beta, mapping="SK")
RnTia_dark, RnTia_dark_short = tm.prepare_data_dark(RnTi, mapping="SK")
Psat308a_dark, Psat308a_dark_short = tm.prepare_data_dark(Psat308, mapping="SK")



#print(RnTi_a)

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


mask_wire = [ 0,0,0,0,0,  #1A
	      0,0,0,0,0,     #2A
	      0,0,0,0,0,        #3A
	      0,0,0,0,0,     #4A
	      0,0,0,0,0,  #5A
	      0,0,0,0,0,  #1B
	      0,0,0,0,0,        #2B
	      0,0,0,0,0,     #3B
	      0,0,0,0,0,     #4B
	      0,0,0,0,0]     #5B



mask_wire_dark = [ 0,0,  #1A
		   0,0,  #2A
		   0,0,  #1B
		   0,0]  #2B






# make plots
outpath = "output/"+DATE+"/TilePlots/"
#tm.plot_a_map(mcecooda,triangles,x,y,mask_wire,0, 0.0, 1, 'mce col',' - ',outpath)


betamask=beta[np.logical_not(np.isnan(beta))]
Gcamask=Gc[np.logical_not(np.isnan(Gc))]

Psat308amask=Psat308[np.logical_not(np.isnan(Psat308))]
Psat308amask=Psat308amask[Psat308amask>0]

Tcamask=Tc[np.logical_not(np.isnan(Tc))]
Tcamask=Tcamask[Tcamask>450]

RnTiamask=RnTi[np.logical_not(np.isnan(RnTi))]
RnTiamask=RnTiamask[RnTiamask<500]


#tm.plot_a_map(Gca,triangles,x,y,mask_wire,0, np.min(np.asarray(Gcamask)), np.max(np.asarray(Gcamask)), 'Gc T6','pW/K',outpath, cmap="jet") #changer 2 and 3 number, lower and upper limit
#tm.plot_a_map(Tca,triangles,x,y,mask_wire,0, np.min(np.asarray(Tcamask)), np.max(np.asarray(Tcamask)), 'Tc T6','mK',outpath, cmap="jet")
#tm.plot_a_map(beta_a,triangles,x,y,mask_wire,0, np.min(np.asarray(betamask)), np.max(np.asarray(betamask)), 'beta T6','',outpath, cmap="jet")
#tm.plot_a_map(RnTia,triangles,x,y,mask_wire,0, np.min(np.asarray(RnTiamask)), np.max(np.asarray(RnTiamask)), 'RnTi T6','mOhm',outpath, cmap="jet")
#tm.plot_a_map(Psat308a,triangles,x,y,mask_wire,0, np.min(np.asarray(Psat308amask)), np.max(np.asarray(Psat308amask)), 'Psat308 T6','pW',outpath, cmap="jet")


tm.plot_a_map(Gca,triangles,x,y,mask_wire,0, 5, 30, 'Gc_T6','pW/K',outpath, cmap="jet") #changer 2 and 3 number, lower and upper limit
tm.plot_a_map(Tca,triangles,x,y,mask_wire,0, 490, 498, 'Tc_T6','mK',outpath, cmap="jet")
tm.plot_a_map(beta_a,triangles,x,y,mask_wire,0, 1.5, 1.9, 'beta_T6','',outpath, cmap="jet")
tm.plot_a_map(RnTia,triangles,x,y,mask_wire,0, 50, 150, 'RnTi_T6','mOhm',outpath, cmap="jet")
tm.plot_a_map(Psat308a,triangles,x,y,mask_wire,0, 1, 4.5, 'Psat308_T6','pW',outpath, cmap="jet")

print('hello')
print(beta_a_dark)
print(beta_a_dark_short)
#print(np.isfinite(Gca_dark))
#mask=np.isfinite(Gca_dark)
#Gca_dark=np.asarray(Gca_dark)
#print(Gca_dark[mask])
#print(len(Gca_dark))

tm.plot_a_map_dark(Gca_dark_short,triangles_dark,x_dark,y_dark,mask_wire_dark,0, 0, 35, 'Gc_T6_Dark','pW/K',outpath, cmap="jet") #changer 2 and 3 number, lower and upper limit
tm.plot_a_map_dark(Tca_dark_short,triangles_dark,x_dark,y_dark,mask_wire_dark,0, 490, 498, 'Tc_T6_Dark','mK',outpath, cmap="jet")
tm.plot_a_map_dark(beta_a_dark_short,triangles_dark,x_dark,y_dark,mask_wire_dark,0, 1.5, 1.9, 'beta_T6_Dark','',outpath, cmap="jet")
tm.plot_a_map_dark(RnTia_dark_short,triangles_dark,x_dark,y_dark,mask_wire_dark,0, 50, 150, 'RnTi_T6_Dark','mOhm',outpath, cmap="jet")
tm.plot_a_map_dark(Psat308a_dark_short,triangles_dark,x_dark,y_dark,mask_wire_dark,0, 1, 4.5, 'Psat308_T6_Dark','pW',outpath, cmap="jet")
