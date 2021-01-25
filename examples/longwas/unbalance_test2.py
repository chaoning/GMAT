import numpy as np
import pandas as pd
from gmat.gmatrix import agmat
from gmat.longwas.unbalance import unbalance_varcom
import logging

logging.basicConfig(level=logging.INFO)

bed_file = '../data/mouse_long/plink'
kin_lst = agmat(bed_file, inv=True, small_val=0.001, out_fmt='id_id_val')
data_file = '../data/mouse_long/phe.unbalance.txt'
id = 'ID'
tpoint = 'weak'
trait = 'trait'
kin_inv_file = '../data/mouse_long/plink.agiv2'
tfix = 'Sex'
prefix_outfile = '../data/mouse_long/unbalance_varcom'
#res_var = unbalance_varcom(data_file, id, tpoint, trait, kin_inv_file, tfix=None, prefix_outfile=prefix_outfile)
#print(res_var)


kin_file = '../data/mouse_long/plink.agrm2'
var_com = pd.read_csv("../data/mouse_long/unbalance_varcom.var", sep='\s+', header=0)

"""
from gmat.longwas.unbalance import unbalance_longwas_fixed
prefix_outfile = '../data/mouse_long/unbalance_longwas_fixed'
res_fixed = unbalance_longwas_fixed(data_file, id, tpoint, trait, bed_file, kin_file, var_com, snp_lst=None, tfix=None,
                             prefix_outfile=prefix_outfile)
"""

from gmat.longwas.unbalance import unbalance_condition_longwas_fixed
prefix_outfile = '../data/mouse_long/unbalance_condition_longwas_fixed3'
res_fixed = unbalance_condition_longwas_fixed(data_file, id, tpoint, trait, bed_file, kin_file, var_com, condition_snp="UNC30664376", snp_lst=None, tfix=None,
                             prefix_outfile=prefix_outfile)

"""
from gmat.longwas.unbalance import unbalance_longwas_fixed_permutation
res_fixed_perm = unbalance_longwas_fixed_permutation(data_file, id, tpoint, trait, bed_file, kin_file, var_com,
                                    tfix=None, prefix_outfile=prefix_outfile)

"""
"""
from gmat.longwas.unbalance import unbalance_longwas_trans

prefix_outfile = '../data/mouse_long/unbalance_longwas_trans'
res_trans = unbalance_longwas_trans(data_file, id, tpoint, trait, bed_file, kin_file, var_com, tfix=None,
                                    prefix_outfile=prefix_outfile)
"""
"""
from gmat.longwas.unbalance import unbalance_condition_longwas_trans

prefix_outfile = '../data/mouse_long/unbalance_condition_longwas_trans'
res_trans = unbalance_condition_longwas_trans(data_file, id, tpoint, trait, bed_file, kin_file, var_com, condition_snp="JAX00177214", tfix=None,
                                    prefix_outfile=prefix_outfile)
"""

"""
from gmat.longwas.unbalance import unbalance_longwas_trans_permutation
res_trans_perm = unbalance_longwas_trans_permutation(data_file, id, tpoint, trait, bed_file, kin_file, var_com,
                                    tfix=None, prefix_outfile=prefix_outfile)
"""

