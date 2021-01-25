"""
Analyze the mouse data. Partition the phenotypic variance into addtive + dominance + additve by additive + additive by dominance
+ dominance by dominance + residual

"""
import logging
logging.basicConfig(level=logging.INFO)

# prepare the phenotypic vector, design matrixed for fixed effects and random effects
from gmat.uvlmm import design_matrix_wemai_multi_gmat

pheno_file = '../data/mouse/pheno'
bed_file = '../data/mouse/plink'
y, xmat, zmat = design_matrix_wemai_multi_gmat(pheno_file, bed_file)

# Calculate the genomic relationship matrix
from gmat.gmatrix.gmatrix import agmat, dgmat_as
a = agmat(bed_file, inv=True)
b = dgmat_as(bed_file, inv=True)

zmat_lst = [zmat, zmat]
gmat_lst = [a[0], b[0]]
gmat_inv_lst = [a[1], b[1]]


from gmat.uvlmm.varcom import em_mme, em_vmat
em_vmat(y, xmat, zmat_lst, gmat_lst, init=None, maxiter=10, cc_par=1e-08)
em_mme(y, xmat, zmat_lst, gmat_inv_lst, init=None, maxiter=10, cc_par=1e-08)

from gmat.uvlmm.varcom import em_mme2, em_vmat2
em_mme2(y, xmat, zmat_lst, gmat_inv_lst, init=None, maxiter=10, cc_par=1e-08)
em_vmat2(y, xmat, zmat_lst, gmat_lst, init=None, maxiter=10, cc_par=1e-08)

from gmat.uvlmm.varcom import em_mme22, em_vmat22
em_mme22(y, xmat, zmat_lst, gmat_inv_lst, init=None, maxiter=10, cc_par=1e-08)
em_vmat22(y, xmat, zmat_lst, gmat_lst, init=None, maxiter=10, cc_par=1e-08)


from gmat.uvlmm.varcom import pxem_mme, pxem_vmat
# gmat_lst = [a[0], b[0], a[0] * a[0], a[0]*b[0], b[0]*b[0]]
pxem_mme(y, xmat, zmat_lst, gmat_inv_lst, init=None, maxiter=10, cc_par=1e-08)
pxem_vmat(y, xmat, zmat_lst, gmat_lst, init=None, maxiter=10, cc_par=1e-08)

from gmat.uvlmm.varcom import pxem_mme2, pxem_vmat2
pxem_mme2(y, xmat, zmat_lst, gmat_inv_lst, init=None, maxiter=100, cc_par=1e-08)
pxem_vmat2(y, xmat, zmat_lst, gmat_lst, init=None, maxiter=100, cc_par=1e-08)

from gmat.uvlmm.varcom import ai_vmat, ai_mme
ai_vmat(y, xmat, zmat_lst, gmat_lst, init=None, maxiter=100, cc_par=1e-08)
ai_mme(y, xmat, zmat_lst, gmat_inv_lst, init=None, maxiter=100, cc_par=1e-08)

from gmat.uvlmm.varcom import fisher_vmat, newton_vmat
newton_vmat(y, xmat, zmat_lst, gmat_lst, init=None, maxiter=100, cc_par=1e-08)
fisher_vmat(y, xmat, zmat_lst, gmat_lst, init=None, maxiter=100, cc_par=1e-08)



from gmat.uvlmm.varcom import weight_emai_vmat, weight_emai_mme
weight_emai_vmat(y, xmat, zmat_lst, gmat_lst, init=None, maxiter=100, cc_par=1e-08)
weight_emai_mme(y, xmat, zmat_lst, gmat_inv_lst, init=None, maxiter=100, cc_par=1e-08)

from gmat.uvlmm.varcom import weight_emnewton_vmat
weight_emnewton_vmat(y, xmat, zmat_lst, gmat_lst, init=None, maxiter=100, cc_par=1e-08)

from gmat.uvlmm.varcom import weight_emfisher_vmat
weight_emfisher_vmat(y, xmat, zmat_lst, gmat_lst, init=None, maxiter=100, cc_par=1e-08)


from gmat.uvlmm.varcom import weight_pxemai_vmat, weight_pxemai_mme
zmat_lst = [zmat, zmat, zmat, zmat, zmat]
gmat_lst = [a[0], b[0], a[0]*a[0], a[0]*b[0], b[0]*b[0]]
weight_pxemai_vmat(y, xmat, zmat_lst, gmat_lst, init=None, maxiter=100, cc_par=1e-08, step=1)
weight_pxemai_mme(y, xmat, zmat_lst, gmat_inv_lst, init=None, maxiter=100, cc_par=1e-08)





# Estimate the variances
from gmat.uvlmm.uvlmm_varcom import _wemai_multi_gmat
gmat_lst = [a[0], b[0]]
res = _wemai_multi_gmat(y, xmat, zmat, gmat_lst, init=None, maxiter=200, cc_par=1.0e-8, cc_gra=1.0e-6)



"""
Analyze the yeast data. Partition the phenotypic variance into addtive + dominance + additve by additive + additive by dominance
+ dominance by dominance + individual-specific residual + residual

"""
import numpy as np
import logging
logging.basicConfig(level=logging.INFO)

# prepare the phenotypic vector, design matrixed for fixed effects and random effects
from gmat.uvlmm import design_matrix_wemai_multi_gmat

pheno_file = '../data/yeast/CobaltChloride'
bed_file = '../data/yeast/CobaltChloride'
y, xmat, zmat = design_matrix_wemai_multi_gmat(pheno_file, bed_file)

# Calculate the genomic relationship matrix
from gmat.gmatrix.gmatrix import agmat
a = agmat(bed_file, inv=False)

# Estimate the variances
from gmat.uvlmm.uvlmm_varcom import wemai_multi_gmat
gmat_lst = [a[0], a[0]*a[0], np.eye(a[0].shape[0])]
res = wemai_multi_gmat(y, xmat, zmat, gmat_lst, init=None, maxiter=200, cc_par=1.0e-8, cc_gra=1.0e-6)



