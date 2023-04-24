import pickle
import numpy as np
'''
dccol : 1-8
dcpad : 1-10
mcecol: 0,1
mcerow: 0-32
'''

with open('Lx_mapping.pckl','rb') as fp:
        mce2det_dict = pickle.load(fp)
        det2mce_dict = pickle.load(fp)

def det2mce(iM, detcol,detrow,detpol):
        try:
                mcecol,mrow = det2mce_dict[detcol,detrow,detpol]
        except Exception as e:
                print("(%d,%d,%d) is not in use."%(detcol,detrow,detpol))
                mcecol,mrow = -1, -1
        return mcecol,mrow

def mce2det(mcecol,mcerow):
        iM = 0
        try:
                detcol,detrow,detpol = mce2det_dict[mcecol,mcerow]
                if(detcol == 'X'):
                        iM,detcol,detrow,detpol = -1,-1, -1,-1                        
        except Exception as e:
                print("(%s,%s) is not in use."%(mcecol,mcerow))
                iM,detcol,detrow,detpol = -1,-1, -1,-1

        return iM,detcol,detrow,detpol



