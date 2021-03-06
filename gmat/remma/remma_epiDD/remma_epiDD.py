import numpy as np
from scipy import linalg
import pandas as pd
from functools import reduce
import gc
import time
import logging
from tqdm import tqdm
import sys
from scipy.stats import chi2

from gmat.process_plink.process_plink import read_plink, impute_geno
from gmat.uvlmm.design_matrix import design_matrix_wemai_multi_gmat


def _remma_epiDD(y, xmat, zmat, gmat_lst, var_com, bed_file, snp_lst_0=None, p_cut=0.0001, out_file='epiDD'):
    """
    dominance by dominance epistasis test by random SNP-BLUP model.
    :param y: phenotypic vector
    :param xmat: Designed matrix for fixed effect
    :param zmat: csr sparse matrix. Designed matrix for random effect.
    :param gmat_lst: A list for relationship matrix
    :param var_com: Estimated variances
    :param bed_file: the prefix for plink binary file
    :param snp_lst_0: the first SNP list for the SNP pairs. the min value is 0 and the max value is num_snp-2. The
    default value is None, which means list [0, num_snp-1)
    :param p_cut: put cut value. default value is 0.0001.
    :param out_file: output file. default value is 'remma_epiDD'.
    :return:
    """
    logging.info("Calculate the phenotypic covariance matrix and inversion")
    y = np.array(y).reshape(-1, 1)
    n = y.shape[0]
    xmat = np.array(xmat).reshape(n, -1)
    vmat = np.diag([var_com[-1]] * n)
    for val in range(len(gmat_lst)):
        vmat += zmat.dot((zmat.dot(gmat_lst[val])).T) * var_com[val]
    # del gmat_lst
    # gc.collect()
    vmat_inv = linalg.inv(vmat)
    logging.info("Calculate P matrix")
    vxmat = np.dot(vmat_inv, xmat)
    xvxmat = np.dot(xmat.T, vxmat)
    xvxmat = linalg.inv(xvxmat)
    pmat = reduce(np.dot, [vxmat, xvxmat, vxmat.T])
    pmat = vmat_inv - pmat
    pymat = zmat.T.dot(np.dot(pmat, y))
    # pvpmat = reduce(np.dot, [pmat, vmat, pmat])  # pvp = p
    pvpmat = zmat.T.dot((zmat.T.dot(pmat)).T)
    del vmat, vmat_inv, pmat
    gc.collect()
    logging.info("Read the SNP")
    np.savetxt(out_file, ['snp_0 snp_1 eff chi p_val'], fmt='%s')
    snp_mat = read_plink(bed_file)
    num_id, num_snp = snp_mat.shape
    if np.any(np.isnan(snp_mat)):
        logging.warning('Missing genotypes are imputed with random genotypes.')
        snp_mat = impute_geno(snp_mat)
    freq = np.sum(snp_mat, axis=0) / (2 * num_id)
    freq.shape = (1, -1)
    scale_vec = 2 * freq * (1 - freq)
    scale = np.sum(scale_vec * (1 - scale_vec))
    logging.info('The scaled factor is: {:.3f}'.format(scale))
    snp_mat[snp_mat > 1.5] = 0.0  # 2?????????0, ??????0???1???0??????
    snp_mat = snp_mat - scale_vec
    logging.info('Test')
    if snp_lst_0 is None:
        snp_lst_0 = range(num_snp - 1)
    else:
        if max(snp_lst_0) >= num_snp - 1 or min(snp_lst_0) < 0:
            logging.error('snp_lst_0 is out of range!')
            sys.exit()
    clock_t0 = time.perf_counter()
    cpu_t0 = time.process_time()
    for i in tqdm(snp_lst_0):
        epi_mat = snp_mat[:, i:(i+1)] * snp_mat[:, (i+1):]
        eff_vec = np.dot(epi_mat.T, pymat)
        var_vec = np.sum(epi_mat * np.dot(pvpmat, epi_mat), axis=0)
        var_vec = var_vec.reshape(len(var_vec), -1)
        chi_vec = eff_vec * eff_vec / var_vec
        p_vec = chi2.sf(chi_vec, 1)
        res = pd.DataFrame(
            {0: np.array([i]*(num_snp-i-1)), 1: np.arange((i+1), num_snp), 2: eff_vec[:, -1], 3: chi_vec[:, -1],
             4: p_vec[:, -1]})
        res = res[res[4] < p_cut]
        res.to_csv(out_file, sep=' ', header=False, index=False, mode='a')
    clock_t1 = time.perf_counter()
    cpu_t1 = time.process_time()
    logging.info("Running time: Clock time, {:.5f} sec; CPU time, {:.5f} sec.".format(clock_t1 - clock_t0, cpu_t1 - cpu_t0))
    return 0


def remma_epiDD(pheno_file, bed_file, gmat_lst, var_com, snp_lst_0=None, p_cut=1.0e-5, out_file='epiDD'):
    """
    Dominance by dominance epistasis test by random SNP-BLUP model.
    :param pheno_file: phenotypic file. The fist two columns are family id, individual id which are same as plink *.fam
    file. The third column is always ones for population mean. The last column is phenotypic values. The ohter covariate
    can be added between columns for population mean and phenotypic values.
    :param bed_file: the prefix for binary file
    :param gmat_lst: a list of genomic relationship matrixes.
    :param var_com: Estimated variances
    :param snp_lst_0: the first SNP list for the SNP pairs. the min value is 0 and the max value is num_snp-2. The
    default value is None, which means list [0, num_snp-1)
    :param p_cut: put cut value. default value is 1.0e-5.
    :param out_file: output file. default value is 'epiDD'.
    :return: 0
    """
    y, xmat, zmat = design_matrix_wemai_multi_gmat(pheno_file, bed_file)
    res = _remma_epiDD(y, xmat, zmat, gmat_lst, var_com, bed_file, snp_lst_0=snp_lst_0, p_cut=p_cut, out_file=out_file)
    return res


def _remma_epiDD_parallel(y, xmat, zmat, gmat_lst, var_com, bed_file, parallel, p_cut=1.0e-4, out_file='epiDD_parallel'):
    """
    Parallel version. Dominance by dominance epistasis test by random SNP-BLUP model.
    :param y: phenotypic vector
    :param xmat: Designed matrix for fixed effect
    :param zmat: csr sparse matrix. Designed matrix for random effect.
    :param gmat_lst: A list for relationship matrix
    :param var_com: Estimated variances
    :param bed_file: the prefix for plink binary file
    :param parallel: A list containing two integers. The first integer is the number of parts to parallel. The second
    integer is the part to run. For example, parallel = [3, 1], parallel = [3, 2] and parallel = [3, 3] mean to divide
    total number of tests into three parts and run parallelly.
    :param p_cut: put cut value. default value is 0.0001.
    :param out_file: The prefix for output file. default value is 'remma_epiDD_parallel'.
    :return: 0
    """
    logging.info("Parallel: " + str(parallel[0]) + ', ' + str(parallel[1]))
    bim_df = pd.read_csv(bed_file + '.bim', header=None)
    num_snp = bim_df.shape[0]
    num_snp_part = int(num_snp/(2*parallel[0]))
    snp_pos_0 = (parallel[1]-1) * num_snp_part
    snp_pos_1 = parallel[1] * num_snp_part
    snp_pos_2 = (2*parallel[0] - parallel[1]) * num_snp_part
    snp_pos_3 = (2*parallel[0] - parallel[1] + 1) * num_snp_part
    if parallel[1] == 1:
        snp_pos_3 = num_snp - 1
    logging.info('SNP position point: ' +
                 ','.join(list(np.array([snp_pos_0, snp_pos_1, snp_pos_2, snp_pos_3], dtype=str))))
    snp_list_0 = list(range(snp_pos_0, snp_pos_1)) + list(range(snp_pos_2, snp_pos_3))
    res = _remma_epiDD(y, xmat, zmat, gmat_lst, var_com, bed_file, snp_lst_0=snp_list_0, p_cut=p_cut,
                    out_file=out_file + '.' + str(parallel[1]))
    return res


def remma_epiDD_parallel(pheno_file, bed_file, gmat_lst, var_com, parallel, p_cut=1.0e-5, out_file='epiDD_parallel'):
    """
    Parallel version. Dominance by dominance epistasis test by random SNP-BLUP model.
    :param pheno_file: phenotypic file. The fist two columns are family id, individual id which are same as plink *.fam
    file. The third column is always ones for population mean. The last column is phenotypic values. The ohter covariate
    can be added between columns for population mean and phenotypic values.
    :param bed_file: the prefix for binary file
    :param gmat_lst: a list of genomic relationship matrixes.
    :param var_com: Estimated variances
    :param parallel: A list containing two integers. The first integer is the number of parts to parallel. The second
    integer is the part to run. For example, parallel = [3, 1], parallel = [3, 2] and parallel = [3, 3] mean to divide
    total number of tests into three parts and run parallelly.
    :param p_cut: put cut value. default value is 1.0e-5.
    :param out_file: output file. default value is 'epiDD_parallel'.
    :return: 0
    """
    y, xmat, zmat = design_matrix_wemai_multi_gmat(pheno_file, bed_file)
    res = _remma_epiDD_parallel(y, xmat, zmat, gmat_lst, var_com, bed_file, parallel, p_cut=p_cut, out_file=out_file)
    return res
