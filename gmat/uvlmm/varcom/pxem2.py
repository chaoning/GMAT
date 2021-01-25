import numpy as np
from scipy import linalg
from functools import reduce
import logging
from scipy.sparse import hstack


def pxem_mme2(y, xmat, zmat_lst, gmat_inv_lst, init=None, maxiter=100, cc_par=1.0e-8):
    """
    Estimate variance parameters with mme based on em
    :param y: Phenotypic vector
    :param xmat: the design matrix for fixed effect
    :param zmat_lst: A list of design matrix for random effects
    :param gmat_inv_lst: A list for the inversion of relationship matrix
    :param init: Default is None. A list for initial values of variance components
    :param maxiter: Default is 100. The maximum number of iteration times.
    :param cc_par: Convergence criteria for update vector.
    :return: The estimated variances.
    """
    num_var = len(gmat_inv_lst) + 1
    var_com = np.ones(num_var)
    if init is not None:
        var_com = np.array(init)
    var_com_new = np.array(var_com) * 100
    y = np.array(y).reshape(-1, 1)
    n = y.shape[0]
    xmat = np.array(xmat).reshape(n, -1)
    xmat_df = xmat.shape[1]
    df_index = [xmat_df]  # the index for the degree of freedom for different effects.
    for i in range(num_var - 1):
        df_index.append(df_index[-1] + gmat_inv_lst[i].shape[0])
    logging.info("Prepare the coefficient matrix without variance parameters")
    zmat_concat = hstack(zmat_lst)
    xzmat = hstack([xmat, zmat_concat])
    coef_null = xzmat.T.dot(xzmat)  # the coefficient matrix
    rhs_null = xzmat.T.dot(y)  # the right hand matrix
    h_mat = np.eye(n) - reduce(np.dot, [xmat, linalg.inv(np.dot(xmat.T, xmat)), xmat.T])
    whw_mat = xzmat.T.dot((xzmat.T.dot(h_mat)).T)
    logging.info("Start iteration")
    iter = 0
    cc_par_val = 10000.0
    while iter < maxiter:
        iter += 1
        logging.info('Round: {:d}'.format(iter))
        # coef matrix
        coef = coef_null.toarray() / var_com[-1]
        for k in range(num_var - 1):
            indexa = df_index[k]
            indexb = df_index[k + 1]
            coef[indexa:indexb, indexa:indexb] = coef[indexa:indexb, indexa:indexb] + gmat_inv_lst[k] / var_com[k]
        # right hand
        rhs_mat = rhs_null/var_com[-1]
        coef_inv = linalg.inv(coef)
        eff = np.dot(coef_inv, rhs_mat)
        smat = h_mat / var_com[-1]
        # error variance
        y_zu = y - zmat_concat.dot(eff[xmat_df:, :])
        var_com_new[-1] = np.sum(reduce(np.dot, [y_zu.T, h_mat, y_zu])) + \
                          np.sum(whw_mat * coef_inv)  # fast trace calculation
        var_com_new[-1] /= n - xmat_df
        # random variances
        symat = np.dot(smat, y)
        gamma_mat = np.zeros((num_var - 1, num_var - 1))
        gamma_vec = np.zeros(num_var - 1)
        for k in range(num_var - 1):
            indexa = df_index[k]
            indexb = df_index[k + 1]
            var_com_new[k] = np.sum(coef_inv[indexa:indexb, indexa:indexb]*gmat_inv_lst[k]) \
                             + np.sum(reduce(np.dot, [eff[indexa:indexb, :].T, gmat_inv_lst[k], eff[indexa:indexb, :]]))
            qk = gmat_inv_lst[k].shape[0]
            var_com_new[k] /= qk
            zu1 = zmat_lst[k].dot(eff[indexa:indexb, :])
            gamma1 = np.sum(np.dot(zu1.T, symat)) - \
                     np.sum(zmat_lst[k].T.dot(np.dot(smat, xmat)) * coef_inv[indexa:indexb, :xmat_df])
            gamma_vec[k] = gamma1
            for m in range(k+1):
                indexc = df_index[m]
                indexd = df_index[m + 1]
                zu2 = zmat_lst[m].dot(eff[indexc:indexd, :])
                gamma2 = np.sum(reduce(np.dot, [zu1.T, smat, zu2])) + \
                         np.sum(zmat_lst[k].T.dot((zmat_lst[m].T.dot(smat)).T) * coef_inv[indexa:indexb, indexc:indexd])
                gamma_mat[k, m] = gamma_mat[m, k] = gamma2
        # print(gamma_mat, gamma_vec)
        gamma = np.dot(linalg.inv(gamma_mat), gamma_vec)
        var_com_new[:-1] = var_com_new[:-1] * gamma * gamma
        # cc
        delta = np.array(var_com_new) - np.array(var_com)
        cc_par_val = np.sum(delta*delta)/np.sum(np.array(var_com_new)*np.array(var_com_new))
        cc_par_val = np.sqrt(cc_par_val)
        var_com = np.array(var_com_new)
        logging.info('Norm of update vector: {:e}'.format(cc_par_val))
        var_com_str = ' '.join(np.array(var_com, dtype=str))
        logging.info('Updated variances: {}'.format(var_com_str))
        if cc_par_val < cc_par:
            break
    if cc_par_val < cc_par:
        logging.info('Variances converged.')
    else:
        logging.info('Variances not converged.')
    return var_com


def pxem_vmat2(y, xmat, zmat_lst, gmat_lst, init=None, maxiter=100, cc_par=1.0e-8):
    """
    :param y:
    :param xmat:
    :param zmat_lst:
    :param gmat_lst:
    :param init:
    :param maxiter:
    :param cc_par:
    :return:
    """
    num_var = len(gmat_lst) + 1
    var_com = np.ones(num_var)
    if init is not None:
        var_com = np.array(init)
    var_com_new = np.array(var_com) * 100
    y = np.array(y).reshape(-1, 1)
    n = y.shape[0]
    xmat = np.array(xmat).reshape(n, -1)
    zgzmat_lst = []
    for val in range(num_var - 1):
        zgzmat_lst.append(zmat_lst[val].dot((zmat_lst[val].dot(gmat_lst[val])).T))
    h_mat = np.eye(n) - reduce(np.dot, [xmat, linalg.inv(np.dot(xmat.T, xmat)), xmat.T])
    logging.info('Initial variances: ' + ' '.join(list(np.array(var_com, dtype=str))))
    logging.info("#####Start the iteration#####\n\n")
    iter = 0
    cc_par_val = 1000.0
    while iter < maxiter:
        iter += 1
        logging.info('###Round: ' + str(iter) + '###')
        # V, the inverse of V and P
        vmat = np.diag([var_com[-1]] * n)
        for val in range(len(zgzmat_lst)):
            vmat += zgzmat_lst[val] * var_com[val]
        zgzmat = vmat - np.diag([var_com[-1]] * n)
        vmat = linalg.inv(vmat)
        vxmat = np.dot(vmat, xmat)
        xvxmat = np.dot(xmat.T, vxmat)
        xvxmat = linalg.inv(xvxmat)
        vxxvxxv_mat = reduce(np.dot, [vxmat, xvxmat, vxmat.T])
        pmat = vmat - vxxvxxv_mat
        pymat = np.dot(pmat, y)
        # Updated variances
        trace_term = np.trace(h_mat) * var_com[-1] - np.sum(h_mat*pmat) * var_com[-1] * var_com[-1]
        y_zu = y - np.dot(zgzmat, pymat)
        var_com_new[-1] = trace_term + np.sum(reduce(np.dot, [y_zu.T, h_mat, y_zu]))
        var_com_new[-1] /= n - xmat.shape[1]
        # random
        smat = h_mat / var_com[-1]
        symat = np.dot(smat, y)
        sxmat = np.dot(smat, xmat)
        gamma_mat = np.zeros((len(gmat_lst), len(gmat_lst)))
        gamma_vec = np.zeros(len(gmat_lst))
        for k in range(num_var - 1):
            var_com_new[k] = -np.sum(zgzmat_lst[k] * pmat) + np.sum(reduce(np.dot, [pymat.T, zgzmat_lst[k], pymat]))
            var_com_new[k] = var_com[k] + var_com[k] * var_com[k] / gmat_lst[k].shape[0] * var_com_new[k]
            zu1 = np.dot(zgzmat_lst[k], pymat) * var_com[k]
            gamma1 = np.sum(np.dot(zu1.T, symat)) + \
                     np.sum(zmat_lst[k].T.dot(sxmat) * reduce(np.dot, [(zmat_lst[k].dot(gmat_lst[k])).T, vxmat,
                                                                        xvxmat])) * var_com[k]
            gamma_vec[k] = gamma1
            for m in range(k+1):
                zu2 = np.dot(zgzmat_lst[m], pymat) * var_com[m]
                zjgizi = zmat_lst[m].dot((zmat_lst[k].dot(gmat_lst[k])).T) * var_com[k]
                zjgjzi = zmat_lst[m].dot((zmat_lst[k].dot(gmat_lst[m])).T) * var_com[m]
                if k == m:
                    gamma2_trace = np.sum(smat * zjgizi.T) - np.sum(np.dot(smat, zjgizi) * np.dot(zjgjzi.T, pmat))
                else:
                    gamma2_trace = - np.sum(np.dot(smat, zjgizi) * np.dot(zjgjzi.T, pmat))
                gamma2 = np.sum(reduce(np.dot, [zu1.T, smat, zu2])) + gamma2_trace
                gamma_mat[k, m] = gamma_mat[m, k] = gamma2
        # print(gamma_mat, gamma_vec)
        gamma = np.dot(linalg.inv(gamma_mat), gamma_vec)
        var_com_new[:-1] = var_com_new[:-1] * gamma * gamma
        delta = var_com_new - var_com
        cc_par_val = np.sum(delta * delta) / np.sum(var_com_new * var_com_new)
        cc_par_val = np.sqrt(cc_par_val)
        var_com = np.array(var_com_new)
        logging.info('Norm of update vector: ' + str(cc_par_val))
        logging.info('Updated variances: ' + ' '.join(list(np.array(var_com_new, dtype=str))))
        if cc_par_val < cc_par:
            break
    if cc_par_val < cc_par:
        logging.info('Variances converged.')
    else:
        logging.info('Variances not converged.')
    return var_com
