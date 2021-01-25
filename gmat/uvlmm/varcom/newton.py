import numpy as np
from scipy import linalg
from functools import reduce
import logging
from scipy.sparse import hstack


def newton_vmat(y, xmat, zmat_lst, gmat_lst, init=None, maxiter=100, cc_par=1.0e-8):
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
        pzgzmat_lst = []
        fd_mat[-1] = np.trace(pmat) - np.sum(np.dot(pymat.T, pymat))
        fd_mat[-1] = -0.5 * fd_mat[-1]
        for i in range(num_var - 1):
            fd_mat[i] = np.sum(pmat*zgzmat_lst[i]) - np.sum(reduce(np.dot, [pymat.T, zgzmat_lst[i], pymat]))
            fd_mat[i] = -0.5 * fd_mat[i]
            wv_mat.append(np.dot(zgzmat_lst[i], pymat))
            pzgzmat_lst.append(np.dot(pmat, zgzmat_lst[i]))
        wv_mat.append(pymat)
        pzgzmat_lst.append(pmat)
        wv_mat = np.concatenate(wv_mat, axis=1)
        nr_mat = reduce(np.dot, [wv_mat.T, pmat, wv_mat])
        for i in range(len(pzgzmat_lst)):
            for j in range(i + 1):
                nr_mat[i, j] = 0.5*np.sum(pzgzmat_lst[i] * pzgzmat_lst[j].T) - nr_mat[i, j]
                nr_mat[j, i] = nr_mat[i, j]
        print(nr_mat)
        delta = np.dot(linalg.inv(nr_mat), fd_mat)
        var_com_new = var_com - delta
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
