import numpy as np
import mce_data
import sys,os
#sys.path.insert(0, "/home/cheng/analysis/ploting")
import easyplots as eps
import matplotlib.pyplot as plt
import pylab as pl

# input raw bias/fb in ADU, only for 1 det
# output biascalib, fbcalib with shift, rnti
def get_LCs(bias, fb, row, col, calib, fitrange = [], out_path = [], flip = 1):

	if not fitrange:
		fitrange = get_default_fitrange()

	fbcalib_ = -1*flip*fb*calib["FB_CAL"][0]
	biascalib_ = bias*calib["BIAS_CAL"][0]

	# Fit
	# fit range
	bias_fit = np.array([bias[i] for i in range((len(bias))) if (bias[i]>fitrange["rnti_low"] and bias[i]<fitrange["rnti_hgh"])])
	fb_fit   = np.array([fb[i] for i in range((len(fb))) if (bias[i]>fitrange["rnti_low"] and bias[i]<fitrange["rnti_hgh"])])
	# fit with poly1
	ceff  = np.polyfit(bias_fit*calib["BIAS_CAL"][0],-1*flip*fb_fit*calib["FB_CAL"][0],1)
	ffunc = np.poly1d(ceff)
	# resistance
	RR = calib["R_SH"]*(1.00/ceff[0]-1.00) #Ohms
	if RR<0:
		RR = float('nan')
	# fit sc
	bias_sc  = np.array([bias[i] for i in range((len(bias))) if (bias[i]>fitrange["sc_low"] and bias[i]<fitrange["sc_hgh"])])
	fb_sc    = np.array([fb[i] for i in range((len(fb))) if (bias[i]>fitrange["sc_low"] and bias[i]<fitrange["sc_hgh"])])
	# fit sc with polfb1
	ceff_sc  = np.polyfit(bias_sc*calib["BIAS_CAL"][0],-1*flip*fb_sc*calib["FB_CAL"][0],1)
	ffunc_sc = np.poly1d(ceff_sc)

	shift_h  = ffunc(0)
	fbcalib = fbcalib_ - shift_h
	biascalib = biascalib_
	rntifit  = ffunc(biascalib) - shift_h

	if out_path:
		fig,ax = eps.presetting(7.4,6,lx="Ib [uA]",ly="Ites [uA]")
		pl.suptitle('Row %02d'%row + ' Col %02d'%col)
		pl.plot(biascalib*1e6,fbcalib*1e6, biascalib*1e6, rntifit*1e6)
		pl.xlim(min(biascalib)*1e6,max(biascalib)*1e6)

		pl.text(0.6, 0.85, 'R = %.2f Ohms'%RR, fontsize=15, transform=ax.transAxes)
		pl.text(0.6, 0.75, 'SC slope = %.2f'%ceff_sc[0], fontsize=15, transform=ax.transAxes)

		fn = os.path.join(out_path,'single_iv_row%02d'%row + '_col%02d_yes.png'%col)
		eps.possetting(fig, ffn = fn, ifleg = False, ifgrid = True, ifshow = False)

	return biascalib,fbcalib,RR



def get_PRs(biascalib, fbcalib, row, col, calib, rnpsat, out_path = [], flip = 1):

	if not out_path:
		doFigure = False
	else:
		doFigure = True

	doFigure = True

	rr = (biascalib/fbcalib-1)*calib["R_SH"]
	pp = fbcalib*fbcalib*rr

	try:
		pphist,ppbinedge = np.histogram(pp,np.linspace(0.0,max(pp),int(max(pp)/0.05e-12)))
		#psat_loweredge = ppbinedge[np.argmax(pphist)]
		psat_loweredge = 1e-13
		#psat_upperedge = ppbinedge[np.argmax(pphist)+1]
		psat_upperedge = 6e-12
		rr_psat = [rr[ii] for ii in range(len(pp)) if pp[ii]<psat_upperedge and pp[ii]>psat_loweredge]
		pp_psat = np.array([ipp for ipp in pp if ipp<psat_upperedge and ipp>psat_loweredge])

		rr_ = abs(rr - rnpsat)
		ind_rr_ = np.argsort(rr_)
		ii = 0

		while True:
			ind = ind_rr_[ii]
			if (pp[ind]<psat_upperedge and pp[ind]>psat_loweredge):
				if rr[ind] < rnpsat:
					ind2 = ind-1
				elif rr[ind] > rnpsat:
					ind2 = ind+1
				else:
					ind2 = ind
				psat = pp[ind] + (pp[ind2]-pp[ind])*(rnpsat-rr[ind])/(rr[ind2]-rr[ind])
				break
			else:
				ii += 1

		#psat = np.mean([ipp for ipp in pp if ipp<psat_upperedge and ipp>psat_loweredge])

		rrhist,rrbinedge = np.histogram(rr,np.linspace(0.0,max(rr),int(max(rr)/0.2e-3)))
		rnti_loweredge = rrbinedge[np.argmax(rrhist)]
		rnti_upperedge = rrbinedge[np.argmax(rrhist)+1]
		rnti = np.mean([irr for irr in rr if irr<rnti_upperedge and irr>rnti_loweredge])


		positive_mask_p=np.where(pp>0)
		pp=pp[positive_mask_p]
		rr=rr[positive_mask_p]
		positive_mask_r=np.where(rr>0)
		pp=pp[positive_mask_r]
		rr=rr[positive_mask_r]


		idx_sort=np.argsort(rr)

		rr_sort=rr[idx_sort]
		pp_sort=pp[idx_sort]


		ind_new=np.argmax(pp_sort)
		R_max=rr_sort[ind_new]
		ind_psat=np.where(rr_sort>=0.8*R_max)[0][0]
		psat_new=pp_sort[ind_psat]
		r_psat=0.8*R_max
		# print('Rmax=', R_max)
		# print('ind_new=', ind_new)
		# print('rr=', rr)

	except Exception:
		rr_psat = [-1]
		pp_psat = [-1]
		rnti = -1
		psat = -1
		psat_new=-1
		r_psat=-1

	#===================================#
	# Plot
	#===================================#


	if doFigure:
		pl.clf()
		pl.suptitle('Row %02d'%row + ' Col %02d'%col)
		pl.xlabel('R [mOhms]', fontsize=15)
		pl.ylabel('P [pW]', fontsize=15)
		pl.plot(rr*1.00e3,pp*1.00e12)
		# if (not np.isnan(rnti)) and (not np.isinf(rnti)) and rnti>0:
		# 	pl.xlim(0,rnti*1e3*1.2)
		# 	pl.ylim(0,10)
		#fn = os.path.join(out_path,'single_pr_row%02d'%row + '_col%02d_yes.png'%col)
		pl.axvline(x=r_psat*1.0e3, color='r', linestyle='--')
		pl.axhline(y=psat_new*1.0e12, color='r', linestyle='--')
		pl.grid()
		#pl.savefig(fn)
		pl.close()



	return rr, pp, rnti, psat_new, r_psat



# return a default fitrange
# based on BA1 setup
# in ADU
def get_default_fitrange():
	fitrange = {}
	fitrange["rnal_low"] = 9000
	fitrange["rnal_hgh"] = 10000
	fitrange["rnti_low"] = 1900
	fitrange["rnti_hgh"] = 2000
	fitrange["sc_low"] = 0
	fitrange["sc_hgh"] = 50
	return fitrange
