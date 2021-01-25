import numpy as np
from scipy import linalg
from functools import reduce
import logging
from scipy.sparse import hstack


def ai_vmat(y, xmat, zmat_lst, gmat_lst, init=None, maxiter=100, cc_par=1.0e-8):
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
    logging.info("#####Prepare#####")
    num_var = len(gmat_lst) + 1
    var_com = np.ones(num_var)
    if init is not None:
        var_com = np.array(init)
    y = np.array(y).reshape(-1, 1)
    n = y.shape[0]
    xmat = np.array(xmat).reshape(n, -1)
    zgzmat_lst = []
    for val in range(len(gmat_lst)):
        zgzmat_lst.append(zmat_lst[val].dot((zmat_lst[val].dot(gmat_lst[val])).T))
    logging.info('Initial variances: ' + ' '.join(list(np.array(var_com, dtype=str))))
    logging.info("#####Start the iteration#####\n\n")
    iter = 0
    cc_par_val = 1000.0
    while iter < maxiter:
        iter += 1
        logging.info('Round: {:d}'.format(iter))
        # V, the inverse of V and P
        vmat = np.diag([var_com[-1]] * n)
        for val in range(len(zgzmat_lst)):
            vmat += zgzmat_lst[val] * var_com[val]
        vmat = np.linalg.inv(vmat)
        vxmat = np.dot(vmat, xmat)
        xvxmat = np.dot(xmat.T, vxmat)
        xvxmat = linalg.inv(xvxmat)
        pmat = vmat - reduce(np.dot, [vxmat, xvxmat, vxmat.T])
        pymat = np.dot(pmat, y)
        # first partial derivatives and working variates
        fd_mat = np.zeros(num_var)
        wv_mat = []
        fd_mat[-1] = np.trace(pmat) - np.sum(np.dot(pymat.T, pymat))
        fd_mat[-1] = -0.5 * fd_mat[-1]
        for i in range(num_var - 1):
            fd_mat[i] = np.sum(pmat*zgzmat_lst[i]) - np.sum(reduce(np.dot, [pymat.T, zgzmat_lst[i], pymat]))
            fd_mat[i] = -0.5 * fd_mat[i]
            wv_mat.append(np.dot(zgzmat_lst[i], pymat))
        wv_mat.append(pymat)
        wv_mat = np.concatenate(wv_mat, axis=1)
        ai_mat = 0.5 * reduce(np.dot, [wv_mat.T, pmat, wv_mat])
        delta = np.dot(linalg.inv(ai_mat), fd_mat)
        var_com_new = var_com + delta
        cc_par_val = np.sum(delta * delta) / np.sum(np.array(var_com_new) * np.array(var_com_new))
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


def ai_mme(y, xmat, zmat_lst, gmat_inv_lst, init=None, maxiter=100, cc_par=1.0e-8):
    """
    Estimate variance parameters with mme based em
    :param y:
    :param xmat:
    :param gmat_inv_lst:
    :param init:
    :param maxiter:
    :param cc_par:
    :param cc_gra:
    :return:
    """
    num_var = len(gmat_inv_lst) + 1
    var_com = np.ones(num_var)
    if init is not None:
        var_com = np.array(init)
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
        for k in range(num_var - 1):
            indexa = df_index[k]
            indexb = df_index[k + 1]
            coef[indexa:indexb, indexa:indexb] = coef[indexa:indexb, indexa:indexb] + gmat_inv_lst[k] / var_com[k]
        # right hand
        rhs_mat = rhs_null / var_com[-1]
        coef_inv = linalg.inv(coef)
        eff = np.dot(coef_inv, rhs_mat)
        e_hat = y - xzmat.dot(eff)
        # first partial derivatives and working variates
        fd_mat = np.zeros(num_var)
        wv_mat = []
        for i in range(num_var - 1):
            indexa = df_index[i]
            indexb = df_index[i + 1]
            fd_mat[i] = np.sum(coef_inv[indexa:indexb, indexa:indexb] * gmat_inv_lst[i]) + \
                        np.sum(reduce(np.dot, [eff[indexa:indexb, :].T, gmat_inv_lst[i], eff[indexa:indexb, :]]))
            fd_mat[i] = -0.5 * gmat_inv_lst[i].shape[0]/var_com[i] + 0.5 * fd_mat[i] / (var_com[i] * var_com[i])
            wv_mati = zmat_lst[i].dot(eff[indexa:indexb, :]) / var_com[i]
            wv_mat.append(wv_mati)
        fd_mat[-1] = np.sum((xzmat.T.dot(xzmat)).multiply(coef_inv)) + np.sum(np.dot(e_hat.T, e_hat))
        fd_mat[-1] = -0.5*n/var_com[-1] + 0.5 * fd_mat[-1] / (var_com[-1] * var_com[-1])
        wv_mati = e_hat / var_com[-1]
        wv_mat.append(wv_mati)
        wv_mat = np.concatenate(wv_mat, axis=1)
        wvxz_mat = xzmat.T.dot(wv_mat)/var_com[-1]
        ai_mat = np.dot(wv_mat.T, wv_mat) / var_com[-1] - reduce(np.dot, [wvxz_mat.T, coef_inv, wvxz_mat])
        ai_mat = 0.5 * ai_mat
        delta = np.dot(linalg.inv(ai_mat), fd_mat)
        var_com_new = var_com + delta
        cc_par_val = np.sum(delta * delta) / np.sum(np.array(var_com_new) * np.array(var_com_new))
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
