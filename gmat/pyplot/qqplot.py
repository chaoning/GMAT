import numpy as np
import logging
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def qqplot(parr):
    """
    :param parr:
    :return:
    """
    p_exp = np.arange(0, 1, 1/len(parr)) + 1/(len(parr)*2)
    p_obs = np.sort(parr)
    logging.info("The median of p values is: {}".format(np.median(p_obs)))
    plt.figure(figsize=(4, 4), dpi=600)
    plt.rc('font', family='Times New Roman')
    plt.scatter(-np.log10(p_exp), -np.log10(p_obs), s=1.0, c='red')
    plt.savefig('QQ.png')



