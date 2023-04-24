import os,sys
import numpy as np
import cPickle as pickle
import Lx_ModuleMapping as minfo
import triangle_mapping_Lx as tri
import histogram as hist
import shutil

module = 'L7'
date = '20230117_'
run = '%s_BA_%s'%(date,module)
fn = run+'.pkl'

infile = '/n/holylfs04/LABS/kovac_lab/www/bk_internal_websites/bicep_array/analysis_logbook/20230117_BA213_OE/plots/' + fn

d = pickle.load(open(infile,'r'))
