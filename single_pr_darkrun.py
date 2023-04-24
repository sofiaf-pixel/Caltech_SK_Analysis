'''
function to get the P-R plot from dark run LC.
'''

import matplotlib.pyplot as plt
import scipy.io as sio
import numpy as np
import math
import sys, os
import pylab as pl
import mce_data

#dp depend on the resolution of LC. 6001 pts from 0 to 60000 take 0.05e-12.
def single_pr_light(col,row,lc_cooked,calib,rnti,out_path,rnpsat,dp):
	
	if not out_path:
		doFigure = False
	else:
		doFigure = True
	
	yCalib_j = lc_cooked["calibFB"]
	biasCalib = lc_cooked["calibBIAS"]
	
	rr = (biasCalib/yCalib_j-1)*calib["R_SH"]
	pp = yCalib_j*yCalib_j*rr
	
	try:
		
		pphist,ppbinedge = np.histogram(pp,np.linspace(50.0e-12,max(pp),int(max(pp)/dp)))
		idxl = np.argmax(pphist)
		idxu = np.argmax(pphist)+1
		psat_loweredge = ppbinedge[idxl]
		psat_upperedge = ppbinedge[idxu]
		
		if rnpsat<rr[min(range(len(pp)), key=lambda i: abs(pp[i]-psat_loweredge))]:
			Nflag = 0
			while not (rnpsat>rr[min(range(len(pp)), key=lambda i: abs(pp[i]-psat_loweredge))] and rnpsat<rr[min(range(len(pp)), key=lambda i: abs(pp[i]-psat_upperedge))]):
				idxl -= 1
				idxu -= 1
				psat_loweredge = ppbinedge[idxl]
				psat_upperedge = ppbinedge[idxu]
				Nflag += 1
				if Nflag>2000:
					break
		
		if rnpsat>rr[min(range(len(pp)), key=lambda i: abs(pp[i]-psat_upperedge))]:
			Nflag = 0
			while not (rnpsat>rr[min(range(len(pp)), key=lambda i: abs(pp[i]-psat_loweredge))] and rnpsat<rr[min(range(len(pp)), key=lambda i: abs(pp[i]-psat_upperedge))]):
				idxl += 1
				idxu += 1
				psat_loweredge = ppbinedge[idxl]
				psat_upperedge = ppbinedge[idxu]
				Nflag += 1
				if Nflag>2000:
					break
		
		ra = [rr[ii] for ii in range(len(pp)) if pp[ii]<=psat_upperedge and pp[ii]>=psat_loweredge],[ipp for ipp in pp if ipp<=psat_upperedge and ipp>=psat_loweredge]
		psat = np.interp(rnpsat,[rr[ii] for ii in range(len(pp)) if pp[ii]<psat_upperedge and pp[ii]>psat_loweredge],[ipp for ipp in pp if ipp<psat_upperedge and ipp>psat_loweredge])
		#psat = np.mean([ipp for ipp in pp if ipp<psat_upperedge and ipp>psat_loweredge])
		rrhist,rrbinedge = np.histogram(rr,np.linspace(0.0,max(rr),int(max(rr)/0.2e-3)))
		rnti_loweredge = rrbinedge[np.argmax(rrhist)]
		rnti_upperedge = rrbinedge[np.argmax(rrhist)+1]
		rnti = np.mean([irr for irr in rr if irr<rnti_upperedge and irr>rnti_loweredge])

	except Exception:
		rnti = -1
		psat = -1
	#===================================#
	# Plot
	#===================================#
	if doFigure:
		pl.clf()
		pl.suptitle('Row %02d'%row + ' Col %02d'%col)
		pl.xlabel('R [mOhms]', fontsize=15)
		pl.ylabel('P [pW]', fontsize=15)
		pl.plot(rr*1.00e3,pp*1.00e12)
		if (not np.isnan(rnti)) and (not np.isinf(rnti)) and rnti>0:
			pl.xlim(0,rnti*1e3*1.2)
			pl.ylim(0,10)
		fn = os.path.join(out_path,'single_pr_row%02d'%row + '_col%02d_yes.png'%col)
		pl.axvline(x=rnti*1.0e3, color='r', linestyle='--')
		pl.axhline(y=psat*1.0e12, color='r', linestyle='--')		
		pl.grid()
		pl.savefig(fn)
		
	prdata = {}
	prdata["resistance"] = rr
	prdata["power"] = pp
	prdata["rnti"] = rnti
	prdata["psat"] = psat 	


		
	return prdata

# return psat if input valid rnpsat.
# make plots if given output path.

def single_pr_dark(col,row,lc_cooked,calib,rnti,out_path = [],rnpsat = -1):
	
	if not out_path:
		doFigure = False
	else:
		doFigure = True
	
	yCalib_j = lc_cooked["calibFB"]
	biasCalib = lc_cooked["calibBIAS"]
	
	rr = (biasCalib/yCalib_j-1)*calib["R_SH"]
	pp = yCalib_j*yCalib_j*rr
	
	try:
		
		pphist,ppbinedge = np.histogram(pp,np.linspace(0.0,max(pp),int(max(pp)/0.05e-12)))
		psat_loweredge = ppbinedge[np.argmax(pphist)]
		psat_upperedge = ppbinedge[np.argmax(pphist)+1]
		rr_psat = [rr[ii] for ii in range(len(pp)) if pp[ii]<psat_upperedge and pp[ii]>psat_loweredge]
		pp_psat = [ipp for ipp in pp if ipp<psat_upperedge and ipp>psat_loweredge]
		if rnpsat > 0:
			psat = np.interp(rnpsat, rr_psat, pp_psat)
		else:
			psat = -1

		#psat = np.mean([ipp for ipp in pp if ipp<psat_upperedge and ipp>psat_loweredge])

		rrhist,rrbinedge = np.histogram(rr,np.linspace(0.0,max(rr),int(max(rr)/0.2e-3)))
		rnti_loweredge = rrbinedge[np.argmax(rrhist)]
		rnti_upperedge = rrbinedge[np.argmax(rrhist)+1]
		rnti = np.mean([irr for irr in rr if irr<rnti_upperedge and irr>rnti_loweredge])

	except Exception:
		rr_psat = [-1]
		pp_psat = [-1]
		rnti = -1
		psat = -1
	#===================================#
	# Plot
	#===================================#
	if doFigure:
		pl.clf()
		pl.suptitle('Row %02d'%row + ' Col %02d'%col)
		pl.xlabel('R [mOhms]', fontsize=15)
		pl.ylabel('P [pW]', fontsize=15)
		pl.plot(rr*1.00e3,pp*1.00e12)
		if (not np.isnan(rnti)) and (not np.isinf(rnti)) and rnti>0:
			pl.xlim(0,rnti*1e3*1.2)
			pl.ylim(0,10)
		fn = os.path.join(out_path,'single_pr_row%02d'%row + '_col%02d_yes.png'%col)
		pl.axvline(x=rnti*1.0e3, color='r', linestyle='--')
		pl.axhline(y=psat*1.0e12, color='r', linestyle='--')		
		pl.grid()
		pl.savefig(fn)
		
	prdata = {}
	prdata["resistance"] = rr
	prdata["power"] = pp
	prdata["rnti"] = rnti
	prdata["rnpsat_low"] = rr_psat[0]
	prdata["rnpsat_high"] = rr_psat[-1]
	prdata["psat"] = psat


		
	return prdata
