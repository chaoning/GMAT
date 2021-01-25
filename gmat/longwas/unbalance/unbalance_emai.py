import math
import numpy as np
import pandas as pd
from scipy import linalg
from scipy import sparse
from scipy.sparse import csr_matrix, isspmatrix, hstack, vstack
import datetime
import logging
from gmat.common.is_positive import is_positive


from .common import *
from .pre_mat import *
from .iter_mat import *


def unbalance_emai(y, xmat, zmat, kin, init=None, maxiter=30, cc_par=1.0e-8, cc_gra=1.0e-6, em_weight_step=0.001):
    num_var_pos = 1
    for i in range(len(zmat)):
        num_var_pos += len(zmat[i])
    y_var = np.var(y)/num_var_pos
    num_record = y.shape[0]
    var_com = []
    eff_ind = [[0, xmat.shape[1]]]  # the index for all effects [start end]
    zmat_con_lst = []  # combined random matrix
    cov_dim = []  # the dim for covariance matrix
    vari = []
    varij = []
    varik = []
    i = 0
    for i in range(len(zmat)):
        temp = [eff_ind[i][-1]]
        zmat_con_lst.append(hstack(zmat[i]))
        cov_dim.append(len(zmat[i]))
        for j in range(len(zmat[i])):
            temp.append(temp[-1] + zmat[i][j].shape[1])
            for k in range(j + 1):
                vari.append(i+1)
                varij.append(j+1)
                varik.append(k+1)
                if j == k:
                    var_com.append(y_var)
                else:
                    var_com.append(0.0)
        eff_ind.append(temp)
    var_com.append(y_var)
    vari.append(i + 2)
    varij.append(1)
    varik.append(1)
    if init is None:
        var_com = np.array(var_com)
    else:
        if len(var_com) != len(init):
            logging.info('ERROR: The length of initial variances should be' + len(var_com))
            exit()
        else:
            var_com = np.array(init)
    var_com_update = np.array(var_com)
    logging.info('***prepare the MME**')
    zmat_con = hstack(zmat_con_lst)  # design matrix for random effects
    wmat = hstack([xmat, zmat_con])  # merged design matrix
    cmat_pure = np.dot(wmat.T, wmat)  # C matrix
    rhs_pure = wmat.T.dot(y)  # right hand
    # em weight vector
    if em_weight_step <= 0.0 or em_weight_step > 1.0:
        logging.info('ERROR: The em weight step should be between 0 (not include) and 1 (include)')
        exit()
    iter_count = 0
    cc_par_val = 1000.0
    cc_gra_val = 1000.0
    delta = 1000.0
    logging.info("initial variances: " + ' '.join(np.array(var_com, dtype=str)))
    covi_mat = pre_covi_mat(cov_dim, var_com)
    if covi_mat is None:
        logging.inf("ERROR: Initial variances is not positive define, please check!")
        exit()
    while iter_count < maxiter:
        iter_count += 1
        logging.info('***Start the iteration: ' + str(iter_count) + ' ***')
        logging.info("Prepare the coefficient matrix")
        cmat = (cmat_pure.multiply(1.0 / var_com[-1])).toarray()
        for i in range(len(cov_dim)):
            if isspmatrix(kin[i]):
                temp = sparse.kron(covi_mat[i], kin[i])
                temp = temp.toarray()
            else:
                temp = linalg.kron(covi_mat[i], kin[i])
            cmat[eff_ind[i + 1][0]:eff_ind[i + 1][-1], \
            eff_ind[i + 1][0]:eff_ind[i + 1][-1]] = \
                np.add(cmat[eff_ind[i + 1][0]:eff_ind[i + 1][-1], \
                       eff_ind[i + 1][0]:eff_ind[i + 1][-1]], temp)
        rhs_mat = np.divide(rhs_pure, var_com[-1])
        cmati = linalg.inv(cmat)
        eff = np.dot(cmati, rhs_mat)
        e = y - xmat.dot(eff[:eff_ind[0][1], :]) - zmat_con.dot(
            eff[eff_ind[0][1]:, :])
        # first-order derivative
        fd_mat = pre_fd_mat_x(cmati, kin, covi_mat, eff, eff_ind, e, cov_dim, zmat_con_lst, wmat, num_record, var_com)
        # AI matrix
        ai_mat = pre_ai_mat(cmati, covi_mat, eff, eff_ind, e, cov_dim, zmat_con_lst, wmat, var_com)
        # EM matrix
        em_mat = pre_em_mat(cov_dim, zmat_con_lst, num_record, var_com)
        # Increase em weight to guarantee variances positive
        gamma = -em_weight_step
        while gamma < 1.0:
            gamma = gamma + em_weight_step
            if gamma >= 1.0:
                gamma = 1.0
            wemai_mat = (1 - gamma) * ai_mat + gamma * em_mat
            delta = np.dot(linalg.inv(wemai_mat), fd_mat)
            var_com_update = var_com + delta
            covi_mat = pre_covi_mat(cov_dim, var_com_update)
            if covi_mat is not None:
                logging.info('EM weight value: ' + str(gamma))
                break
        logging.info('Updated variances: ' + ' '.join(np.array(var_com_update, dtype=str)))
        if covi_mat is None:
            logging.info("ERROR: Updated variances is not positive define!")
            exit()
        # Convergence criteria
        cc_par_val = np.sum(pow(delta, 2)) / np.sum(pow(var_com_update, 2))
        cc_par_val = np.sqrt(cc_par_val)
        cc_gra_val = np.sqrt(np.sum(pow(fd_mat, 2))) / len(var_com)
        var_com = var_com_update.copy()
        logging.info("Change in parameters, Norm of gradient vector: " + str(cc_par_val) + ', ' + str(cc_gra_val))
        if cc_par_val < cc_par and cc_gra_val < cc_gra:
            break
    if cc_par_val < cc_par and cc_gra_val < cc_gra:
        logging.info("Variances Converged")
    else:
        logging.info("Variances not Converged")
    var_pd = {'vari': vari,
              "varij": varij,
              "varik": varik,
              "var_val": var_com}
    var_pd = pd.DataFrame(var_pd, columns=['vari', "varij", "varik", "var_val"])
    return var_pd
