import matplotlib.pyplot as plt
import os.path
import shutil
import pylab as pl
import mce_data
#sys.path.insert(0, "/home/cheng/analysis/load_curve/sk_dark/scripts_house")                                                                         
import calib_SK
import get_LCs as lc
import single_pr_darkrun
import ba40_ModuleMapping as mm
import numpy as np

from Gtemplists import gettemplists
templists = gettemplists()

RUN = 'ba40Gscript'
rnpsat = 0.040

DATE = '20190730'
templist = templists[DATE]

def GTcModel(T,G,Tc,beta):
        return G/(beta+1)*(Tc**(beta+1)-T**(beta+1))/Tc**beta

def GTcModel2(T,G,Tc,beta=2):
        #beta = 2.                                                                                                                                   
        return G/(2+1)*(Tc**(2+1)-T**(2+1))/Tc**2

def main():

	Gfn = 'sk_G_%s'%(DATE)
	lcdata = {}
	prdata = {}
	colors = pl.cm.jet(np.linspace(0,1,len(templist)))

	R_n = np.zeros((2,33,len(templist)))
	R_sc=np.zeros((2,33,len(templist)))
	
	i=0
	for temp in templist:
                #===================================#                                                                                                            
                # Loading data                                                                                                                                   
                #===================================#                                                                                                            
		in_path  = './%s/'%(DATE)
		filename = 'LC_G_FPU_'+str(temp)+'mK_datamode1_run1'
		out_path_main = './output/%s/SC_slope'%(DATE)
		if not os.path.isdir(out_path_main):
			os.makedirs(out_path_main)
		out_path = out_path_main + filename + '/'
		if not os.path.isdir(out_path):
			os.makedirs(out_path)
		out_Gplot = out_path_main + '%s_/'%(RUN)
		if not os.path.isdir(out_Gplot):
			os.makedirs(out_Gplot)

		datafn = in_path + filename + '/' + filename

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
		rows = np.zeros((nc,nr),dtype=np.int)
		cols = np.zeros((nc,nr),dtype=np.int)

		#rowmax=33
		rowmax=33

                #===================================#                                                                                                            
                # Cook LC                                                                                                                                        
                #===================================#                                                                                                            
		for col in [0,1]:
			for row in range(rowmax):
				fitrange = {}

				fitrange["rnti_low"] = 5000.00
				fitrange["rnti_hgh"] = 5500.00
				fitrange["sc_low"] = 0.00
				fitrange["sc_hgh"] = 50.00

				ibias, ites, rnti1 = lc.get_LCs(bias, y[row,col], row, col, calib, fitrange)    #if out_path is included in the arguments, the function will produce plot.


				ibias_sc=ibias[len(ibias)-10:len(ibias)]*1e6
				ites_sc=ites[len(ibias)-10:len(ibias)]*1e6
				
				ibias_n=ibias[0:50]*1e6
				ites_n=ites[0:50]*1e6

				m_sc, off_sc=np.polyfit(ibias_sc, ites_sc, 1)
				m_n, off_n=np.polyfit(ibias_n, ites_n, 1)

				R_n[col,row,i]=m_n 
				R_sc[col,row,i]=m_sc

				print(m_n)
				print(R_n[col,row,i])

				fig = plt.figure(figsize=(20,6.5), dpi=80)
					
				plt.plot(ibias*1e6, ites*1e6)
				plt.plot(ibias_sc, ites_sc, c='r', label="Ang_coeff_sc=%.3f"%m_sc)
				plt.plot(ibias_n, ites_n, c='g', label="R_n=%.3f"%m_n)
				plt.suptitle('Row %02d'%row + ' Col %02d'%col)
				plt.title('Temperature %02d'%temp+'mK')
				plt.xlabel('Ibias [uA]', fontsize=15)
				plt.ylabel('Ites [uA]', fontsize=15)
				plt.legend()
				#plt.show()

				rowcol_path='Row'+str(row)+'_Col'+str(col)+'/'
				new_path=out_Gplot

				if not os.path.isdir(new_path):
					os.makedirs(new_path)

				fn = os.path.join(new_path,'Row %02d'%row +'_Col %02d'%col+'_Temp %02d'%temp)
				plt.savefig(fn)

				rr, pp, rnti2, psat = lc.get_PRs(ibias, ites, row, col, calib, rnpsat)

				if rnti2<0:
					rnti=rnti1
				else:
					rnti=rnti2

				lcdata[str(temp)+'_r'+str(row)+'c'+str(col)] = [ibias, ites]
				prdata[str(temp)+'_r'+str(row)+'c'+str(col)] = [rr, pp, rnti, psat, 1/(ites*rr)]

		i=i+1

	for col in [0,1]:
		for row in range(rowmax):
			
			im,detcol,detrow,detpol=mm.mce2det(col,row)
			
			rowcol_path='det_row'+str(detrow)+'_col'+str(detcol)+'/'
			new_path=out_Gplot+rowcol_path

			fig = plt.figure(figsize=(20,6.5), dpi=80)
			plt.tick_params(labelsize=18)
			mask_n=R_n[col,row,:]<1
			bb= np.linspace(R_n[col,row,mask_n].min()-.001, R_n[col,row,mask_n].max()+.001, 100)
			plt.hist(R_n[col,row,:], bins=bb)  # arguments are passed to np.histogram                                                                               
			plt.xlabel('R_n',fontsize=18)
			plt.title("R_n -"+" Mean="+str(round(np.mean(R_n[col,row,:]),3)), fontsize=18)
			fn_histo = os.path.join('slopes/','R_n_'+'Row%02d'%detrow + 'Col%02d'%detcol)
			plt.savefig(fn_histo)


			fig = plt.figure(figsize=(20,6.5), dpi=80)
			plt.tick_params(labelsize=18)
			R_sc_mask=(R_sc[col, row,:]>0.7)&(R_sc[col, row,:]<1.3)
			bb= np.linspace(1.08,1.12 , 8)
			plt.hist(R_sc[col,row,:], bins=bb)  # arguments are passed to np.histogram                                                                                     
			plt.xlabel('Sc_slopes',fontsize=18)
			plt.title("Sc_slopes -"+" Mean="+str(round(np.mean(R_sc[col,row,R_sc_mask]),3)), fontsize=18)
			fn_histo_sc = os.path.join("slopes/"+'sk_slopes_row'+str(row)+'_col'+str(col)+'.png')
			plt.savefig(fn_histo_sc)



if __name__=='__main__':
    main()

