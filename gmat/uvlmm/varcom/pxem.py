import numpy as np
from scipy import linalg
from functools import reduce
import logging
from scipy.sparse import hstack


def pxem_mme(y, xmat, zmat_lst, gmat_inv_lst, init=None, maxiter=100, cc_par=1.0e-8):
    """
    Estimate variance parameters with mme based on pxem
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
    logging.info('Initial variances: ' + ' '.join(list(np.array(var_com, dtype=str))))
    logging.info("Start iteration")
    iter = 0
    cc_par_val = 10000.0
    while iter < maxiter:
        iter += 1
        logging.info('Round: {:d}'.format(iter))
        # coef matrix
        coef = coef_null.toarray() / var_com[-1]
        for k in range(num_var-1):
            indexa = df_index[k]
            indexb = df_index[k+1]
            coef[indexa:indexb, indexa:indexb] = coef[indexa:indexb, indexa:indexb] + gmat_inv_lst[k]/var_com[k]
        # right hand
        rhs_mat = rhs_null/var_com[-1]
        coef_inv = linalg.inv(coef)
        eff = np.dot(coef_inv, rhs_mat)
        # error variance
        e_hat = y - xzmat.dot(eff)
        var_com_new[-1] = np.sum(np.dot(e_hat.T, e_hat)) + \
                          np.sum((xzmat.T.dot(xzmat)).multiply(coef_inv))  # fast trace calculation
        var_com_new[-1] /= y.shape[0]
        # random variances
        y_xb = y - np.dot(xmat, eff[:xmat_df, :])
        gamma_mat = np.zeros((num_var - 1, num_var - 1))
        gamma_vec = np.zeros(num_var - 1)
        for k in range(num_var-1):
            indexa = df_index[k]
            indexb = df_index[k + 1]
            var_com_new[k] = np.sum(coef_inv[indexa:indexb, indexa:indexb]*gmat_inv_lst[k]) \
                             + np.sum(reduce(np.dot, [eff[indexa:indexb, :].T, gmat_inv_lst[k], eff[indexa:indexb, :]]))
            qk = gmat_inv_lst[k].shape[0]
            var_com_new[k] /= qk
            zu1 = zmat_lst[k].dot(eff[indexa:indexb, :])
            gamma1 = np.sum(np.dot(zu1.T, y_xb)) - \
                     np.sum(zmat_lst[k].T.dot(xmat) * coef_inv[indexa:indexb, :xmat_df])
            gamma_vec[k] = gamma1/var_com[-1]
            for m in range(k+1):
                indexc = df_index[m]
                indexd = df_index[m + 1]
                zu2 = zmat_lst[m].dot(eff[indexc:indexd, :])
                gamma2 = np.sum(np.dot(zu1.T, zu2)) + \
                         np.sum((zmat_lst[k].T.dot(zmat_lst[m])).multiply(coef_inv[indexa:indexb, indexc:indexd]))
                gamma_mat[k, m] = gamma2/var_com[-1]
        gamma_mat = gamma_mat + np.tril(gamma_mat, k=-1).T
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


def pxem_vmat(y, xmat, zmat_lst, gmat_lst, init=None, maxiter=100, cc_par=1.0e-8):
    """
    Estimate variance parameters with vmat based on pxem
    :param y: Phenotypic vector
    :param xmat: the design matrix for fixed effect
    :param zmat_lst: A list of design matrix for random effects
    :param gmat_lst: A list for relationship matrix
    :param init: Default is None. A list for initial values of variance components
    :param maxiter: Default is 100. The maximum number of iteration times.
    :param cc_par: Convergence criteria for update vector.
    :return: The estimated variances.
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
        vmat = linalg.inv(vmat)
        vxmat = np.dot(vmat, xmat)
        xvxmat = np.dot(xmat.T, vxmat)
        xvxmat = linalg.inv(xvxmat)
        pmat = vmat - reduce(np.dot, [vxmat, xvxmat, vxmat.T])
        pymat = np.dot(pmat, y)
        # Updated variances
        var_com_new[-1] = var_com[-1] - var_com[-1] * var_com[-1] * np.trace(pmat)/n + \
                          var_com[-1] * var_com[-1] * np.sum(np.dot(pymat.T, pymat))/n
        # random
        y_xb = y - reduce(np.dot, [xmat, xvxmat, vxmat.T, y])
        gamma_mat = np.zeros((num_var - 1, num_var - 1))
        gamma_vec = np.zeros(num_var - 1)
        for k in range(num_var - 1):
            var_com_new[k] = -np.sum(zgzmat_lst[k] * pmat) + np.sum(reduce(np.dot, [pymat.T, zgzmat_lst[k], pymat]))
            var_com_new[k] = var_com[k] + var_com[k] * var_com[k]/gmat_lst[k].shape[0] * var_com_new[k]
            zu1 = np.dot(zgzmat_lst[k], pymat) * var_com[k]
            gamma1 = np.sum(np.dot(zu1.T, y_xb)) + \
                     np.sum(zmat_lst[k].T.dot(xmat) * reduce(np.dot, [(zmat_lst[k].dot(gmat_lst[k])).T, vxmat, xvxmat])) * var_com[k]
            gamma_vec[k] = gamma1/var_com[-1]
            for m in range(k+1):
                zu2 = np.dot(zgzmat_lst[m], pymat) * var_com[m]
                zjgizi = zmat_lst[m].dot((zmat_lst[k].dot(gmat_lst[k])).T) * var_com[k]
                zjgjzi = zmat_lst[m].dot((zmat_lst[k].dot(gmat_lst[m])).T) * var_com[m]
                if k == m:
                    gamma2_trace = np.trace(zjgizi) - np.sum(zjgizi * np.dot(zjgjzi.T, pmat))
                else:
                    gamma2_trace = -np.sum(zjgizi * np.dot(zjgjzi.T, pmat))
                gamma2 = np.sum(np.dot(zu1.T, zu2)) + gamma2_trace
                gamma_mat[k, m] = gamma2 / var_com[-1]
        gamma_mat = gamma_mat + np.tril(gamma_mat, k=-1).T
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

