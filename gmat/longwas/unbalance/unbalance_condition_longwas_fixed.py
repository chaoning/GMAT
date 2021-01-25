import numpy as np
import pandas as pd
from scipy import linalg
from scipy import sparse
from scipy.sparse import csr_matrix, hstack
from scipy.stats import chi2
import datetime
import gc
import logging
from tqdm import tqdm
from patsy import dmatrix

from gmat.process_plink.process_plink import read_plink, impute_geno
from pysnptools.snpreader import Bed
from .common import *


def unbalance_condition_longwas_fixed(data_file, id, tpoint, trait, bed_file, kin_file, var_com, condition_snp, snp_lst=None,
            tfix=None, fix=None, forder=3, aorder=3, porder=3, na_method='omit',
                             prefix_outfile='unbalance_condition_longwas_fixed'):
    """
    the longitudinal GWAS for the unbalanced data treating the SNP as the time varied fixed effect.
    :param data_file: the data file. The first row is the variate names whose first initial position is alphabetical.
    For the class variates, the first letter must be capital; for the covariates (continuous variates), the first letter
    must be lowercase.
    :param id: A class variate name which indicates the individual id column in the data file.
    :param tpoint: A covariate names which indicates the time point column in the data file.
    :param trait: A variate name which indicates the analyzed trait column in the data file.
    :param bed_file: the prefix for the plink binary file.
    :param kin_file: the file for genomic relationship matrix. This file can be produced by
    gmat.gmatrix.agmat function using agmat(bed_file, inv=True, small_val=0.001, out_fmt='id_id_val')
    :param var_com: the estimated variance parameters by the gmat.longwas.unbalance.unbalance_varcom function.
    :param condition_snp: conditional snp
    :param snp_lst: the snp list to test. Default is None.
    :param tfix: A class variate name for the time varied fixed effect. Default value is None. Only one time varied
    fixed effect can be included in the current version.
    :param fix: Expression for the time independent fixed effect. Default value is None. An example:
    fix = "Sex + age + Season".
    :param forder: the order of Legendre polynomials for the time varied fixed effect. The default value is 3.
    :param aorder: the order of Legendre polynomials for the additive genetic effect. The default value is 3.
    :param porder: the order of Legendre polynomials for the permanent environment effect. The default value is 3.
    :param na_method: The method to deal with missing values. The default value is 'omit'. 'omit' method will delete the
    row with missing values. 'include' method will fill the missing values with the adjacent values.
    :param prefix_outfile: the prefix for the output file. Default is 'unbalance_longwas_fixed'.
    :return: A pandas data frame for the test result.
    """
    logging.info('################################')
    logging.info('###Prepare the related matrix###')
    logging.info('################################')
    if var_com.shape[0] != aorder*(aorder+1)/2 + aorder + 1 + porder*(porder+1)/2 + porder + 1 + 1:
        logging.info('ERROR: Variances do not match the data, please check')
        exit()
    logging.info('***Read the data file***')
    logging.info('Data file: ' + data_file)
    data_df = pd.read_csv(data_file, sep='\s+', header=0)
    logging.info('NA method: ' + na_method)
    if na_method == 'omit':
        data_df = data_df.dropna()
    elif na_method == 'include':
        data_df = data_df.fillna(method='ffill')
        data_df = data_df.fillna(method='bfill')
    else:
        logging.info('na_method does not exist: ' + na_method)
        exit()
    col_names = data_df.columns
    logging.info('The column names of data file: ' + ' '.join(list(col_names)))
    logging.info('Note: Variates beginning with a capital letter is converted into factors.')
    class_vec = []
    for val in col_names:
        if not val[0].isalpha():
            logging.info("The first character of columns names must be alphabet!")
            exit()
        if val[0] == val.capitalize()[0]:
            class_vec.append(val)
            data_df[val] = data_df[val].astype('str')
        else:
            try:
                data_df[val] = data_df[val].astype('float')
            except Exception as e:
                logging.info(e)
                logging.info(val + " may contain string, please check!")
                exit()
    logging.info('Individual column: ' + id)
    if id not in col_names:
        logging.info(id + ' is not in the data file, please check!')
        exit()
    if id not in class_vec:
        logging.info('The initial letter of {} should be capital'.format(id))
        exit()
    id_order = []
    id_arr = list(data_df[id])
    id_order.append(id_arr[0])
    for i in range(1, len(id_arr)):
        if id_arr[i] != id_arr[i - 1]:
            id_order.append(id_arr[i])
    id_in_data = set(data_df[id])
    if len(id_in_data) - len(id_order) != 0:
        logging.info('The data is not sored by individual ID!')
        exit()
    logging.info('Time points column: ' + tpoint)
    if tpoint not in col_names:
        logging.info(tpoint + ' is not in the data file, please check!')
        exit()
    if tpoint in class_vec:
        logging.info('The initial letter of {} should be lowercase'.format(tpoint))
        exit()
    logging.info('Trait column: ' + trait)
    if trait not in col_names:
        logging.info(trait + ' is not in the data file, please check!')
        exit()
    if trait in class_vec:
        logging.info('The initial letter of {} should be lowercase'.format(trait))
        exit()
    logging.info('Code factor variables of the data file: ' + ' '.join(list(class_vec)))
    code_val = {}
    code_dct = dct_2D()
    for val in class_vec:
        code_val[val] = 0
        temp = []
        for i in range(data_df.shape[0]):
            if data_df[val][i] not in code_dct[val]:
                code_val[val] += 1
                code_dct[val][data_df[val][i]] = str(code_val[val])
            temp.append(code_dct[val][data_df[val][i]])
        data_df[val] = np.array(temp)
    for val in class_vec:
        data_df[val] = data_df[val].astype('int')
    logging.info('***Build the design matrix for fixed effect***')
    logging.info('Time dependent fixed effect: ' + str(tfix))
    leg_fix = leg(data_df[tpoint], forder)
    if tfix == None:
        xmat_t = np.concatenate(leg_fix, axis=1)
        xmat_t = csr_matrix(xmat_t)
    else:
        if tfix not in class_vec:
            logging.info(tfix + ' is not the class variate')
            exit()
        row = np.array(range(data_df.shape[0]))
        col = np.array(data_df[tfix]) - 1
        val = np.array([1.0] * data_df.shape[0])
        tfix_mat = csr_matrix((val, (row, col)))
        xmat_t = []
        for i in range(len(leg_fix)):
            xmat_t.append(tfix_mat.multiply(leg_fix[i]))
        xmat_t = hstack(xmat_t)
        del row, col, val
        gc.collect()
    logging.info('Time independent fix effect: ' + str(fix))
    xmat_nt = None
    if fix == None:
        xmat_nt = None
    else:
        try:
            fix_exp = ''
            vec = fix.split('+')
            for i in vec:
                val = i.strip()
                if val in class_vec:
                    fix_exp += 'C(' + val + ')'
                else:
                    fix_exp += val
            xmat_nt = dmatrix(fix_exp, data_df)
            logging.info('The expression for fixed effect: ' + fix_exp)
        except Exception as e:
            logging.info(e + ': Check the fix effect expression.')
            exit()
        xmat_nt = csr_matrix(xmat_nt[:, 1:])
    xmat = hstack([xmat_t, xmat_nt])
    xmat = xmat.toarray()
    max_id = max(data_df[id]) + 1
    tmin = min(data_df[tpoint])
    tmax = max(data_df[tpoint])
    leg_lst = []  # legendre polynomials for time dependent fixed SNP effects, save for each individuals
    for i in range(1, max_id):
        leg_lst.append(leg_mt(data_df[data_df[id] == i][tpoint], tmax, tmin, forder))
    tpoint_vec = sorted(set(data_df[tpoint]))
    leg_tpoint_mat = leg_mt(np.array(tpoint_vec), tmax, tmin, forder)
    leg_tpoint_accum = np.sum(leg_tpoint_mat, axis=0)
    logging.info('***Read the kinship matrix***')
    logging.info('Kinship file: ' + kin_file)
    with open(kin_file) as fin:
        row = []
        col = []
        kin = []
        id_in_kin = {}
        for line in fin:
            arr = line.split()
            id_in_kin[arr[0]] = 1
            id_in_kin[arr[1]] = 1
            if arr[0] not in code_dct[id]:
                logging.info(arr[0] + ' is not in the kinship inversion file!')
                exit()
            if arr[1] not in code_dct[id]:
                logging.info(arr[1], 'is not in the kinship inversion file!')
                exit()
            row.append(int(code_dct[id][arr[0]]))
            col.append(int(code_dct[id][arr[1]]))
            kin.append(float(arr[2]))
    id_not_in_kin = list(set(code_dct[id].keys()) - set(id_in_kin.keys()))
    if len(id_not_in_kin) != 0:
        logging.info('The ID: {} in the data file is not in the kinship file!'.format(' '.join(id_not_in_kin)))
        exit()
    kin = csr_matrix((np.array(kin), (np.array(row) - 1, np.array(col) - 1))).toarray()
    kin = np.add(kin, kin.T)
    kin[np.diag_indices_from(kin)] = 0.5 * np.diag(kin)
    del row, col
    gc.collect()
    logging.info('***Build the dedign matrix for random effect***')
    logging.info('Legendre order for additive effects: ' + str(aorder))
    leg_add = leg(data_df[tpoint], aorder)
    row = np.array(range(data_df.shape[0]))
    col = np.array(data_df[id]) - 1
    val = np.array([1.0] * data_df.shape[0])
    add_mat = csr_matrix((val, (row, col)), shape=(data_df.shape[0], kin.shape[0]))
    zmat_add = []
    for i in range(len(leg_add)):
        zmat_add.append(add_mat.multiply(leg_add[i]))
    logging.info('Legendre order for permanent environmental effect: ' + str(porder))
    leg_per = leg(data_df[tpoint], porder)
    per_mat = csr_matrix((val, (row, col)))
    zmat_per = []
    for i in range(len(leg_per)):
        zmat_per.append((per_mat.multiply(leg_per[i])))
    del row, col, val
    gc.collect()
    zmat = [zmat_add, zmat_per]
    y = data_df[trait].values.reshape(data_df.shape[0], 1)
    # kin_inv = [kin_inv, sparse.eye(max(data_df[id]), format="csr")]
    logging.info('***Prepare the merged Z matrix***')
    eff_ind = [[0, xmat.shape[1]]]  # the index for all effects [start end]
    zmat_con_lst = []  # combined random matrix
    for i in range(len(zmat)):
        temp = [eff_ind[i][-1]]
        zmat_con_lst.append(hstack(zmat[i]))
        for j in range(len(zmat[i])):
            temp.append(temp[-1] + zmat[i][j].shape[1])
        eff_ind.append(temp)
    logging.info('***Calculate the phenotypic (co)variance***')
    add_cov = var_com.loc[var_com.loc[:, 'vari']==1, :]
    row = np.array(add_cov['varij']) - 1
    col = np.array(add_cov['varik']) - 1
    val = add_cov['var_val']
    add_cov = csr_matrix((val, (row, col))).toarray()
    add_cov = add_cov + np.tril(add_cov, k=-1).T
    per_cov = var_com.loc[var_com.loc[:, 'vari']==2, :]
    row = np.array(per_cov['varij']) - 1
    col = np.array(per_cov['varik']) - 1
    val = per_cov['var_val']
    per_cov = csr_matrix((val, (row, col))).toarray()
    per_cov = per_cov + np.tril(per_cov, k=-1).T
    res_var = np.array(var_com['var_val'])[-1]
    vmat = zmat_con_lst[0].dot((zmat_con_lst[0].dot(np.kron(add_cov, kin))).T)
    one_id = sparse.eye(zmat_con_lst[1].shape[1]/per_cov.shape[0])
    vmat = vmat + zmat_con_lst[1].dot((zmat_con_lst[1].dot(sparse.kron(per_cov, one_id))).T)
    vmat_diag = np.diag(vmat) + res_var
    np.fill_diagonal(vmat, vmat_diag)
    vmat = linalg.inv(vmat)
    logging.info('***Read the snp data***')
    # snp_mat = read_plink(bed_file)
    snp_on_disk = Bed(bed_file, count_A1=False)
    num_id = snp_on_disk.iid_count
    num_snp = snp_on_disk.sid_count
    logging.info("There are {:d} individuals and {:d} SNPs.".format(num_id, num_snp))
    fam_df = pd.read_csv(bed_file + '.fam', sep='\s+', header=None)
    id_geno = list(np.array(fam_df.iloc[:, 1], dtype=str))
    id_order_index = []
    for i in id_order:
        id_order_index.append(id_geno.index(i))
    if snp_lst is None:
        snp_lst = range(num_snp)
    snp_lst = list(snp_lst)
    if min(snp_lst) < 0 or max(snp_lst) >= num_snp:
        logging.info('The value in the snp list should be >= {} and < {}', 0, num_snp)
        exit()
    snp_mat = snp_on_disk[:, snp_lst].read().val
    if np.any(np.isnan(snp_mat)):
        logging.info('Missing genotypes are imputed with random genotypes.')
        snp_mat = impute_geno(snp_mat)
    snp_mat = snp_mat[id_order_index, :]
    condition_snp_index = snp_on_disk.sid_to_index([condition_snp])[0]
    condition_snp_val = snp_on_disk[:, condition_snp_index].read().val
    condition_snp_val = condition_snp_val[id_order_index, 0]
    # snp_mat = snp_mat[:, snp_lst]
    logging.info('#####################################################################')
    logging.info('###Start the fixed regression longitudinal GWAS for unbalance data###')
    logging.info('#####################################################################')
    chi_df = leg_lst[0].shape[1]
    eff_con_vec = []
    chi_con_vec = []
    p_con_vec = []
    p_min_con_vec = []
    p_accum_con_vec = []
    eff_vec = []
    chi_vec = []
    p_vec = []
    p_min_vec = []
    p_accum_vec = []
    snp_condition = list(map(lambda x, y: x*y, leg_lst, list(condition_snp_val)))
    snp_condition = np.concatenate(snp_condition, axis=0)
    xmat = np.concatenate((xmat, snp_condition), axis=1)
    for i in tqdm(range(snp_mat.shape[1])):
        snp_fix = list(map(lambda x, y: x*y, leg_lst, list(snp_mat[:, i])))
        snp_fix = np.concatenate(snp_fix, axis=0)
        snp_fix = np.concatenate((xmat, snp_fix), axis=1)
        xv = np.dot(snp_fix.T, vmat)
        xvx = np.dot(xv, snp_fix)
        snp_index = list(range(-chi_df, 0))
        con_snp_index = list(range(-2*chi_df, -chi_df))
        try:
            xvx = np.linalg.inv(xvx)
        except Exception as e:
            del e
            con_snp_index = snp_index[:]
            snp_fix = snp_fix[:, :-chi_df]
            xv = np.dot(snp_fix.T, vmat)
            xvx = np.dot(xv, snp_fix)
            xvx = np.linalg.inv(xvx)
        xvy = np.dot(xv, y)
        b = np.dot(xvx, xvy)
        eff = b[snp_index, -1]
        eff_var = xvx[snp_index, :]
        eff_var = eff_var[:, snp_index]
        chi_val = np.sum(np.dot(np.dot(eff.T, np.linalg.inv(eff_var)), eff))
        p_val = chi2.sf(chi_val, chi_df)
        eff_vec.append(eff)
        chi_vec.append(chi_val)
        p_vec.append(p_val)
        p_tpoint_vec = []
        for k in range(leg_tpoint_mat.shape[0]):
            eff_tpoint = np.sum(np.dot(leg_tpoint_mat[k, :], eff))
            eff_var_tpoint = np.sum(np.dot(leg_tpoint_mat[k, :], np.dot(eff_var, leg_tpoint_mat[k, :])))
            chi_tpoint = eff_tpoint*eff_tpoint/eff_var_tpoint
            p_tpoint = chi2.sf(chi_tpoint, 1)
            p_tpoint_vec.append(p_tpoint)
        p_min_vec.append(min(p_tpoint_vec))
        eff_accum = np.sum(np.dot(leg_tpoint_accum, eff))
        eff_var_accum = np.sum(np.dot(leg_tpoint_accum, np.dot(eff_var, leg_tpoint_accum)))
        chi_accum = eff_accum*eff_accum/eff_var_accum
        p_accum = chi2.sf(chi_accum, 1)
        p_accum_vec.append(p_accum)
        # conditional SNP
        eff = b[con_snp_index, -1]
        eff_var = xvx[con_snp_index, :]
        eff_var = eff_var[:, con_snp_index]
        chi_val = np.sum(np.dot(np.dot(eff.T, np.linalg.inv(eff_var)), eff))
        p_val = chi2.sf(chi_val, chi_df)
        eff_con_vec.append(eff)
        chi_con_vec.append(chi_val)
        p_con_vec.append(p_val)
        p_tpoint_vec = []
        for k in range(leg_tpoint_mat.shape[0]):
            eff_tpoint = np.sum(np.dot(leg_tpoint_mat[k, :], eff))
            eff_var_tpoint = np.sum(np.dot(leg_tpoint_mat[k, :], np.dot(eff_var, leg_tpoint_mat[k, :])))
            chi_tpoint = eff_tpoint * eff_tpoint / eff_var_tpoint
            p_tpoint = chi2.sf(chi_tpoint, 1)
            p_tpoint_vec.append(p_tpoint)
        p_min_con_vec.append(min(p_tpoint_vec))
        eff_accum = np.sum(np.dot(leg_tpoint_accum, eff))
        eff_var_accum = np.sum(np.dot(leg_tpoint_accum, np.dot(eff_var, leg_tpoint_accum)))
        chi_accum = eff_accum * eff_accum / eff_var_accum
        p_accum = chi2.sf(chi_accum, 1)
        p_accum_con_vec.append(p_accum)
    logging.info('Finish association analysis')
    logging.info('***Output***')
    snp_info_file = bed_file + '.bim'
    snp_info = pd.read_csv(snp_info_file, sep='\s+', header=None)
    res_df = snp_info.iloc[snp_lst, [0, 1, 3, 4, 5]]
    res_df.columns = ['chro', 'snp_ID', 'pos', 'allele1', 'allele2']
    res_df.loc[:, 'order'] = snp_lst
    res_df = res_df.iloc[:, [5, 0, 1, 2, 3, 4]]
    # conditional SNP
    eff_con_vec = np.array(eff_con_vec)
    for i in range(eff_con_vec.shape[1]):
        col_ind = 'con_eff' + str(i)
        res_df[col_ind] = eff_con_vec[:, i]
    res_df['con_chi_val'] = chi_con_vec
    res_df['con_p_val'] = p_con_vec
    res_df['con_p_min'] = p_min_con_vec
    res_df['con_p_accum'] = p_accum_con_vec
    # SNP
    eff_vec = np.array(eff_vec)
    for i in range(eff_vec.shape[1]):
        col_ind = 'eff' + str(i)
        res_df[col_ind] = eff_vec[:, i]
    res_df['chi_val'] = chi_vec
    res_df['p_val'] = p_vec
    res_df['p_min'] = p_min_vec
    res_df['p_accum'] = p_accum_vec
    out_file = prefix_outfile + '.res'
    res_df.to_csv(out_file, sep=' ', index=False)
    return res_df
