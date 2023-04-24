#Module containing functions to make tile plots
#modified to include BA2 option
# SF - 2022/04

import matplotlib.pyplot as plt
import matplotlib.tri as tri
import scipy.io as sio
import numpy as np
import math
import sys,os
import pickle as pickle
import ba40_ModuleMapping_BA as bamap
import ba150_ModuleMapping_BA as ba2map
import ba40_ModuleMapping as skmap

Ndet_mce = 16*41
Ndet_mce_dark= 2*4
Ndet_det = 16*41
npixel = 18*18 # pixels in one tile
npixel_dark = 4
ntile = 1
lld = 18 # how many dets in one row pyhsically
lld_dark = 2 #fake map for dark det
w,h = 2*npixel,ntile
w_dark=2*npixel_dark

def prepare_triangles(lld):
	xy = []
	for i in range((lld+1)*(lld+1)):
		xy.append([])
		for j in range(2):
			xy[i].append(0)
	for i in range((lld+1)*(lld+1)):
		xy[i][0] = (i-i%(lld+1))/(lld+1)+1
		xy[i][1] = i%(lld+1)+1
	xy = np.asarray(xy)
	x = xy[:, 0]
	y = xy[:, 1]

	triangles = []
	for i in range((lld)*(2*lld)):
		triangles.append([])
		for j in range(3):
			triangles[i].append(0)
	for i in range(lld):
		for j in range(lld*2):
			if j<lld:
				triangles[(j*lld+i)][0] = j*(lld+1)+i
				triangles[(j*lld+i)][1] = j*(lld+1)+i+1
				triangles[(j*lld+i)][2] = j*(lld+1)+i+lld+1
			else:
				triangles[(j*lld+i)][0] = (j-(lld-1))*(lld+1)+i
				triangles[(j*lld+i)][1] = (j-(lld-1))*(lld+1)+i+1
				triangles[(j*lld+i)][2] = (j-(lld-1))*(lld+1)+i-lld
	triangles = np.asarray(triangles)
	return x,y,triangles



def prepare_data(input_inmce, itile=0, mapping="BA2"):
	data_arr = [float('nan') for x in range(w)]
#da = [0]*h
	DATA = []
	IDet_col = []
	IDet_row = []
	POLAR = []

	if mapping == 'BA2':
		tot_mcerows = 41
		tot_mcecols = 16
	else:
		tot_mcerows = 33
		tot_mcecols = 2

	for icol in range(tot_mcecols):
		for irow in range(tot_mcerows):
			DATA.append(input_inmce[icol][irow])
			if mapping=="BA":
				[im, dc,dr,dp] = bamap.mce2det(2*itile+icol,irow)
			if mapping=="BA2":
				[dc,dr,dp] = ba2map.mce2det(2*itile+icol,irow)
			else:
				[im, dc,dr,dp] = skmap.mce2det(icol,irow)
			IDet_col.append(dc)
			IDet_row.append(dr)
			POLAR.append(dp)
	for idet in range(Ndet_mce):
		if POLAR[idet] == 'A':

			# if IDet_col[idet]<0 or IDet_col[idet]<0:
			# 	continue
			idet_col = int(IDet_col[idet])
			idet_row = int(IDet_row[idet])
			data_arr[(idet_col-1)*lld+(idet_row-1)] = DATA[idet]

		elif POLAR[idet] == 'B':

			# if IDet_col[idet]<0 or IDet_col[idet]<0:
            #                     continue
			idet_col = int(IDet_col[idet])
			idet_row = int(IDet_row[idet])
			data_arr[(idet_col-1)*lld+(idet_row-1)+npixel] = DATA[idet]
		else:
			continue
	#for itile in range(ntile):
	#	da[itile] = np.asarray(data_arr[itile])

	return data_arr



def prepare_data_dark(input_inmce, itile=0, mapping="BA"):

	data_arr = [float('nan') for x in range(w_dark)]
	data_arr_A = []
	data_arr_B = []

	DATA=[]
	IDet_col=[]
	IDet_row=[]
	IDet_col_dark = []
	IDet_row_dark = []
	POLAR = []
	POLAR_dark = []


	for icol in [0,1]:
		for irow in range(33):
			#DATA.append(input_inmce[icol][irow])
			if mapping == "BA":
				[im, dc,dr,dp] = bamap.mce2det(2*itile+icol,irow)
				n=0
			else:
				[im, dc,dr,dp] = skmap.mce2det(icol,irow)
				n=1

			print("n:", n)
			i=0
			if dp=='D':
				DATA.append(input_inmce[icol][irow])
				print(irow)
				print(icol)
				print(input_inmce[icol][irow])
				if icol==0:
					dc=0
					if irow==7+n:
						dr=0
						dp='DB'       #O
						data_arr_A.append(input_inmce[icol][irow])
					if irow==8+n:
						dr=0
						dp='DB' #O
						data_arr_B.append(input_inmce[icol][irow])
					if irow==9+n:
						dr=1
						dp='DA' #I
						data_arr_A.append(input_inmce[icol][irow])
					if irow==10+n:
						dr=1
						dp='DA'  #I
						data_arr_B.append(input_inmce[icol][irow])

				if icol==1:
					dc=1
					if irow==7+n:
						dr=0
						dp='DA' #I
						data_arr_A.append(input_inmce[icol][irow])
					if irow==8+n:
						dr=0
						dp='DA'#I
						data_arr_B.append(input_inmce[icol][irow])
					if irow==9+n:
						dr=1
						dp='DB'#O
						data_arr_A.append(input_inmce[icol][irow])
					if irow==10+n:
						dr=1
						dp='DB' #O
						data_arr_B.append(input_inmce[icol][irow])
				IDet_col_dark.append(dc)
				IDet_row_dark.append(dr)
				POLAR_dark.append(dp)

				i=i+1


			#IDet_col.append(dc)
                        #IDet_row.append(dr)
                        #POLAR.append(dp)

	for idet in range(Ndet_mce_dark):
		if POLAR_dark[idet] == 'DA':
			#if IDet_col_dark[idet]<0 or IDet_col_dark[idet]<0:
				#continue
			idet_col = int(IDet_col_dark[idet])
			idet_row = int(IDet_row_dark[idet])
			data_arr[(idet_col-1)*lld_dark+(idet_row-1)] = DATA[idet]
			#data_arr_A.append(DATA[idet])

		elif POLAR_dark[idet] == 'DB':

			#if IDet_col_dark[idet]<0 or IDet_col_dark[idet]<0:
				#continue
			idet_col = int(IDet_col_dark[idet])
			idet_row = int(IDet_row_dark[idet])
			data_arr[(idet_col-1)*lld_dark+(idet_row-1)+npixel_dark] = DATA[idet]
			#data_arr_B.append(DATA[idet])
		else:
			continue

        #for itile in range(ntile):
        #       da[itile] = np.asarray(data_arr[itile])
	#data_arr = np.isfinite(data_arr)
		data_arr_short=data_arr_A+data_arr_B

	return data_arr, data_arr_short



def prepare_data_alpha(inpath,filename):

	Tc1_arr = [[float('nan') for x in range(w)] for y in range(h)]
	Tc2_arr = [[float('nan') for x in range(w)] for y in range(h)]
	Tc_aver_arr = [[float('nan') for x in range(w)] for y in range(h)]
	r1_arr = [[float('nan') for x in range(w)] for y in range(h)]
	r2_arr = [[float('nan') for x in range(w)] for y in range(h)]
	Tc1a = [0]*h
	Tc2a = [0]*h
	Tcaa = [0]*h
	r1a = [0]*h
	r2a = [0]*h

	infile = inpath+filename
	if os.path.exists(infile):
                d2 = pickle.load(open(infile,'r'))

	ICol = []
	IRow = []
	IDet_col = []
	IDet_row = []
	POLAR = []
	TC1 = []
	TC2 = []
	TCaver = []
	R1 = []
	R2 = []

	for icol in [0,1]:
		for irow in range(33):
			ICol.append(icol)
			IRow.append(irow)
			[im,dc,dr,dp]=bamap.mce2det(icol,irow)
			IDet_col.append(dc)
			IDet_row.append(dr)
			POLAR.append(dp)
			TC1.append(d2['t1'][icol][irow])
			TC2.append(d2['t2'][icol][irow])
			TCaver.append((d2['t1'][icol][irow]+d2['t2'][icol][irow])/2.)
			R1.append(d2['r1'][icol][irow])
			R2.append(d2['r1'][icol][irow]+d2['r2'][icol][irow])

	for idet in range(Ndet_mce):

		if POLAR[idet] == 'A':

			polar = 'A'

			icol = int(ICol[idet])

			irow = int(IRow[idet])


			if IDect_col[idet]<0 or IDet_col[idet]<0:

				continue
			idet_col = int(IDet_col[idet])

			idet_row = int(IDet_row[idet])

			itile = 0

			Tc1_arr[itile][(idet_col-1)*lld+(idet_row-1)] = TC1[idet]

			Tc2_arr[itile][(idet_col-1)*lld+(idet_row-1)] = TC2[idet]

			Tc_aver_arr[itile][(idet_col-1)*lld+(idet_row-1)] = TCaver[idet]

			r1_arr[itile][(idet_col-1)*lld+(idet_row-1)] = R1[idet]

			r2_arr[itile][(idet_col-1)*lld+(idet_row-1)] = R2[idet]

		elif POLAR[idet] == 'B':

			polar = 'B'

			icol = int(ICol[idet])

			irow = int(IRow[idet])

			if IDet_col[idet]<0 or IDet_col[idet]<0:

				continue

			idet_col = int(IDet_col[idet])

			idet_row = int(IDet_row[idet])

			itile = 0

			Tc1_arr[itile][(idet_col-1)*lld+(idet_row-1)+npixel] = TC1[idet]

			Tc2_arr[itile][(idet_col-1)*lld+(idet_row-1)+npixel] = TC2[idet]

			Tc_aver_arr[itile][(idet_col-1)*lld+(idet_row-1)+npixel] = TCaver[idet]

			r1_arr[itile][(idet_col-1)*lld+(idet_row-1)+npixel] = R1[idet]

			r2_arr[itile][(idet_col-1)*lld+(idet_row-1)+npixel] = R2[idet]
		else:
			continue

	for itile in range(ntile):

		Tc1a[itile] = np.asarray(Tc1_arr[itile])

		Tc2a[itile] = np.asarray(Tc2_arr[itile])

		Tcaa[itile] = np.asarray(Tc_aver_arr[itile])

		r1a[itile] = np.asarray(r1_arr[itile])

		r2a[itile] = np.asarray(r2_arr[itile])

	return Tc1a, Tc2a, Tcaa, r1a, r2a

def prepare_data_G(inpath,filename):

	Rn_arr = [[float('nan') for x in range(w)] for y in range(h)]

	Gc_arr = [[float('nan') for x in range(w)] for y in range(h)]

	beta_arr = [[float('nan') for x in range(w)] for y in range(h)]

	Tc_arr = [[float('nan') for x in range(w)] for y in range(h)]

	Ra = [0]*h

	Gca = [0]*h

	Tca = [0]*h

	ba = [0]*h

	infile = inpath+filename

	if os.path.exists(infile):
                d2 = pickle.load(open(infile,'r'))


	ICol = []

	IRow = []

	IDet_col = []

	IDet_row = []

	POLAR = []

	RN = []

	G = []

	BETA = []

	TC = []

	for icol in [0,1]:

		for irow in range(33):

			ICol.append(d2[3]['r'+str(irow)+'c'+str(icol)]['mcecol'])

			IRow.append(d2[3]['r'+str(irow)+'c'+str(icol)]['mcerow'])

			IDet_col.append(d2[3]['r'+str(irow)+'c'+str(icol)]['detcol'])

			IDet_row.append(d2[3]['r'+str(irow)+'c'+str(icol)]['detrow'])

			POLAR.append(d2[3]['r'+str(irow)+'c'+str(icol)]['pol'])

			RN.append(d2[2]['r'+str(irow)+'c'+str(icol)]['rnti'])

			G.append(d2[2]['r'+str(irow)+'c'+str(icol)]['Gc'])

			BETA.append(d2[2]['r'+str(irow)+'c'+str(icol)]['beta'])

			TC.append(d2[2]['r'+str(irow)+'c'+str(icol)]['Tc'])

	for idet in range(Ndet_mce):

		if POLAR[idet] == 'A':

			polar = 'A'

			icol = int(ICol[idet])
                        #if icol == 11:
                        #       continue

			irow = int(IRow[idet])

			if np.isnan(IDet_col[idet]) or np.isnan(IDet_col[idet]):
				continue

			if IDet_col[idet]<0 or IDet_col[idet]<0:
				continue

			idet_col = int(IDet_col[idet])

			idet_row = int(IDet_row[idet])
                        #itile = (icol-icol%4)/4

			itile = 0 # for E14 strange wiring

			Rn_arr[itile][(idet_col-1)*lld+(idet_row-1)] = RN[idet]

			Gc_arr[itile][(idet_col-1)*lld+(idet_row-1)] = G[idet]

			Tc_arr[itile][(idet_col-1)*lld+(idet_row-1)] = TC[idet]

			beta_arr[itile][(idet_col-1)*lld+(idet_row-1)] = BETA[idet]

		elif POLAR[idet] == 'B':

			polar = 'B'

			icol = int(ICol[idet])
                        #if icol == 11:
                        #        continue

			irow = int(IRow[idet])

			if np.isnan(IDet_col[idet])|np.isnan(IDet_col[idet]):

				continue

			if IDet_col[idet]<0 or IDet_col[idet]<0:

				continue

			idet_col = int(IDet_col[idet])

			idet_row = int(IDet_row[idet])
                        #itile = (icol-icol%4)/4

			itile = 0 # for E14 strange wiring

			Rn_arr[itile][(idet_col-1)*lld+(idet_row-1)+npixel] = RN[idet]

			Gc_arr[itile][(idet_col-1)*lld+(idet_row-1)+npixel] = G[idet]

			Tc_arr[itile][(idet_col-1)*lld+(idet_row-1)+npixel] = TC[idet]

			beta_arr[itile][(idet_col-1)*lld+(idet_row-1)+npixel] = BETA[idet]

		else:

			continue


	for itile in range(ntile):

		Ra[itile] = np.asarray(Rn_arr[itile])

		Gca[itile] = np.asarray(Gc_arr[itile])

		ba[itile] = np.asarray(beta_arr[itile])

		Tca[itile] = np.asarray(Tc_arr[itile])

	return Ra,Gca,ba,Tca


def plot_a_map(data_array,triangles,x,y,mask_wire,itile,value_min,value_max,plot_title,plot_unit,inpath, savefn, ifshow=False, cmap="jet"):

	fig = plt.figure()

	plt.gca().set_aspect('equal')



	data_array = np.ma.masked_where(np.isinf(mask_wire),data_array)

	data_array = np.ma.masked_where(np.isnan(data_array),data_array)

	plt.tripcolor(x, y, triangles, facecolors=np.array([float(np.isinf(mm)+1)*0.5 for mm in mask_wire]), cmap=plt.cm.Greys_r, edgecolors='k')

	plt.tripcolor(x, y, triangles, facecolors=data_array, cmap=cmap, vmin=value_min, vmax=value_max , edgecolors='k')

	clb = plt.colorbar()
	clb.ax.set_title(plot_unit)

	if not itile =='None':
		plt.title(plot_title+' Tile{0:d}'.format(itile+1))
	else:
		plt.title(plot_title)

	plt.xlabel('idet_col',fontsize=15)

	plt.ylabel('idet_row',fontsize=15)

	plt.tick_params(axis='both',which='major',labelsize=15)

	plt.text(1.2, 1.2, 'A', fontsize=12)

	plt.text(1.6, 1.6, 'B', fontsize=12)

	fig.savefig(inpath+savefn+'.png')

	if ifshow:
        	plt.show()



def plot_a_map_dark(data_array,triangles,x,y,mask_wire,itile,value_min,value_max,plot_title,plot_unit,inpath, mapping='BA', ifshow=False, cmap="jet"):

	fig = plt.figure()

	plt.gca().set_aspect('equal')

	data_array = np.ma.masked_where(np.isinf(mask_wire),data_array)

	data_array = np.ma.masked_where(np.isnan(data_array),data_array)

	plt.tripcolor(x, y, triangles, facecolors=np.array([float(np.isinf(mm)+1)*0.5 for mm in mask_wire]), cmap=plt.cm.Greys_r, vmin=0.0, vmax=1.0 , edgecolors='k')

	plt.tripcolor(x, y, triangles, facecolors=data_array, cmap=cmap, vmin=value_min, vmax=value_max , edgecolors='k')

	plt.colorbar()

	plt.title(plot_title+' '+plot_unit+' Tile{0:d}'.format(itile+1))

	plt.xlabel('idet_col',fontsize=15)

	plt.ylabel('idet_row',fontsize=15)

	plt.tick_params(axis='both',which='major',labelsize=15)

	plt.text(1.5, 2.9, 'mce col 0', fontsize=10)
	plt.text(2.5, 2.9, 'mce col 1', fontsize=10)

	if mapping=='BA':

		plt.text(1.2, 1.2,'I', fontsize=15)
		plt.text(1.6, 1.6, 'O', fontsize=15)

	if mapping=='SK':

                plt.text(1.2, 1.2,'O', fontsize=15)
                plt.text(1.6, 1.6, 'I', fontsize=15)



	fig.savefig(inpath+plot_title+'_Tile%gid.png'%itile)

	if ifshow:
		plt.show()




def triange_mapping_single(input_inmce,mask_wire,vmin,vmax,vari,unit,outpath):
	[x,y,triangles] = prepare_triangles(lld)
	da = prepare_data(input_inmce)
	print(mask_wire)
	for itile in range(ntile):
                plot_a_map(da,triangles,x,y,mask_wire,itile,vmin,vmax,vari,unit,outpath)


def triangle_mapping():


	RUN = 'ba40M5T6'
	#inpath = '/home/cheng/analysis/DATA/output/20181125/'
	inpath = '/output/20180726/ba40script_G_Gplot/'
	filename = 'sk_G_20190726.pkl'
	[x,y,triangles] = prepare_triangles(lld)
	#[Tc1a, Tc2a, Tcaa,r1a,r2a] = prepare_data_alpha(inpath,filename)
	[Ra,Gca,ba,Tca] = prepare_data_G(inpath,filename)

# Rather than create a Triangulation object, can simply pass x, y and triangles
# arrays to tripcolor directly.  It would be better to use a Triangulation
# object if the same triangulation was to be used more than once to save
# duplicated calculations.
# Can specify one color value per face rather than one per point by using the
# facecolors kwarg.
		#row #1,2,3,4,5		#col

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

	for itile in range(ntile):

		plot_a_map(Gca,triangles,x,y,mask_wire,itile,20.0, 70.0, 'Gc','pW/K',inpath)
		plot_a_map(ba,triangles,x,y,mask_wire,itile,1.8, 2.6, 'beta',' - ',inpath)
		plot_a_map(Tca,triangles,x,y,mask_wire,itile,460,540,'Tc','mK',inpath)
		#plot_a_map(Tc2a,triangles,x,y,mask_wire,itile,0.495,0.510,'Tc2','K',inpath)
		#plot_a_map(Tcaa,triangles,x,y,mask_wire,itile,0.495,0.510,'TTaver','K',inpath)
		plot_a_map(Ra,triangles,x,y,mask_wire,itile,100,250,'RnTi','mOhms',inpath)
		'''
		fig = plt.figure()
		plt.gca().set_aspect('equal')

		Ra[itile] = np.ma.masked_where(np.isinf(mask_wire),Ra[itile])
		Ra[itile] = np.ma.masked_where(np.isnan(Ra[itile]),Ra[itile])
		plt.tripcolor(x, y, triangles, facecolors=np.array([float(np.isinf(mm)+1)*0.5 for mm in mask_wire]), cmap=plt.cm.Greys_r, vmin=0.0, vmax=1.0 , edgecolors='k')
		plt.tripcolor(x, y, triangles, facecolors=Ra[itile], cmap=plt.cm.jet, vmin=100, vmax=250 , edgecolors='k')

		plt.colorbar()
		plt.title('Rn [mOhms] Tile{0:d}'.format(itile+1))
		plt.xlabel('idet_col',fontsize=15)
		plt.ylabel('idet_row',fontsize=15)
		plt.tick_params(axis='both',which='major',labelsize=15)
		plt.text(1.2, 1.2, 'A', fontsize=15)
		plt.text(1.6, 1.6, 'B', fontsize=15)
		fig.savefig(inpath+'Rn_map.png')
		plt.show()

		fig = plt.figure()
                plt.gca().set_aspect('equal')

                Gca[itile] = np.ma.masked_where(np.isinf(mask_wire),Gca[itile])
                Gca[itile] = np.ma.masked_where(np.isnan(Gca[itile]),Gca[itile])
                plt.tripcolor(x, y, triangles, facecolors=np.array([float(np.isinf(mm)+1)*0.5 for mm in mask_wire]), cmap=plt.cm.Greys_r, vmin=0.0, vmax=1.0 , edgecolors='k')
                plt.tripcolor(x, y, triangles, facecolors=Gca[itile], cmap=plt.cm.jet, vmin=30, vmax=60 , edgecolors='k')

                plt.colorbar()
                plt.title('Gc [pW/K] Tile{0:d}'.format(itile+1))
                plt.xlabel('idet_col',fontsize=15)
                plt.ylabel('idet_row',fontsize=15)
                plt.tick_params(axis='both',which='major',labelsize=15)
                plt.text(1.2, 1.2, 'A', fontsize=15)
                plt.text(1.6, 1.6, 'B', fontsize=15)
                fig.savefig(inpath+'G_map.png')
		plt.show()


		fig = plt.figure()
                plt.gca().set_aspect('equal')

                ba[itile] = np.ma.masked_where(np.isinf(mask_wire),ba[itile])
                ba[itile] = np.ma.masked_where(np.isnan(ba[itile]),ba[itile])
                plt.tripcolor(x, y, triangles, facecolors=np.array([float(np.isinf(mm)+1)*0.5 for mm in mask_wire]), cmap=plt.cm.Greys_r, vmin=0.0, vmax=1.0 , edgecolors='k')
                plt.tripcolor(x, y, triangles, facecolors=ba[itile], cmap=plt.cm.jet, vmin=2.0, vmax=2.8 , edgecolors='k')

                plt.colorbar()
                plt.title('beta Tile{0:d}'.format(itile+1))
                plt.xlabel('idet_col',fontsize=15)
                plt.ylabel('idet_row',fontsize=15)
                plt.tick_params(axis='both',which='major',labelsize=15)
                plt.text(1.2, 1.2, 'A', fontsize=15)
                plt.text(1.6, 1.6, 'B', fontsize=15)
                fig.savefig(inpath+'beta_map.png')
		plt.show()

		fig = plt.figure()
                plt.gca().set_aspect('equal')

                Tca[itile] = np.ma.masked_where(np.isinf(mask_wire),Tca[itile])
                Tca[itile] = np.ma.masked_where(np.isnan(Tca[itile]),Tca[itile])
                plt.tripcolor(x, y, triangles, facecolors=np.array([float(np.isinf(mm)+1)*0.5 for mm in mask_wire]), cmap=plt.cm.Greys_r, vmin=0.0, vmax=1.0 , edgecolors='k')
                plt.tripcolor(x, y, triangles, facecolors=Tca[itile], cmap=plt.cm.jet, vmin=495, vmax=510 , edgecolors='k')

                plt.colorbar()
                plt.title('Tc [mK] Tile{0:d}'.format(itile+1))
                plt.xlabel('idet_col',fontsize=15)
                plt.ylabel('idet_row',fontsize=15)
                plt.tick_params(axis='both',which='major',labelsize=15)
                plt.text(1.2, 1.2, 'A', fontsize=15)
                plt.text(1.6, 1.6, 'B', fontsize=15)
                fig.savefig(inpath+'Tc_map.png')
		plt.show()


		fig = plt.figure()
		plt.gca().set_aspect('equal')
		DPDTa[itile] = np.ma.masked_where(np.isnan(DPDTa[itile]),DPDTa[itile])
		plt.tripcolor(x, y, triangles, facecolors=DPDTa[itile], vmin=0.085, vmax=0.1 , edgecolors='k')
		plt.colorbar()
		plt.title('dPdT [pW/K] Tile{0:d}'.format(itile+1))
		plt.xlabel('idet_col', fontsize=15)
		plt.ylabel('idet_row', fontsize=15)
		plt.tick_params(axis='both',which='major',labelsize=15)
		plt.text(1.2, 1.2, 'A', fontsize=15)
		plt.text(1.6, 1.6, 'B', fontsize=15)
		#plt.show()
		fig.savefig('/home/cheng/analysis/DATA/output/20180103/run1/E14_Tile{0:d}_dpdt.png'.format(itile+1))
'''
	#return 0

if __name__ == '__main__':
	#in_path = str(sys.argv[1])
	sys.stdout.write(str(triangle_mapping()))
