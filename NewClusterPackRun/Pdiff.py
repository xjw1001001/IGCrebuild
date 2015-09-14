import numpy as np
from scipy.linalg import expm
from CodonGeneconFunc import *
import matplotlib.pyplot as plt


def repack_Geneconv_mat(test):
    if test.Model == 'HKY':
        state_size = 4
        mat_basic = test.get_HKYBasic()
        mat_geneconv = np.zeros((state_size**2, state_size**2))
        mat_rcr = test.get_HKYGeneconv()[1]  # rcr stands for row column rate
    elif test.Model == 'MG94':
        state_size = 61
        mat_basic = test.get_MG94Basic()
        mat_geneconv = np.zeros((state_size ** 2, state_size ** 2))
        mat_rcr = test.get_MG94Geneconv_and_MG94()[1]  # rcr stands for row column rate
    row = [ state_size * i[0] + i[1] for i in mat_rcr['row']]
    col = [ state_size * i[0] + i[1] for i in mat_rcr['col']]
    mat_geneconv[row, col] = mat_rcr['rate']

    # now add in diagonal entries
    mat_basic = mat_basic - np.diag(mat_basic.sum(axis = 1))
    mat_geneconv = mat_geneconv - np.diag(mat_geneconv.sum(axis = 1))
    return mat_basic, mat_geneconv

def get_Pmat(mat_basic, mat_geneconv, basic_t):
    Pmat_basic = [np.matrix(expm(mat_basic * t)) for t in basic_t]
    Pmat_geneconv = [np.matrix(expm(mat_geneconv * t)) for t in basic_t]
    return Pmat_basic, Pmat_geneconv

def plot_pdiff(test,basic_t, p = None):
    if p == None:
        p = test.prior_distribution
    mat_basic, mat_geneconv = repack_Geneconv_mat(test)
    if test.Model == 'HKY':
        dup_prior = np.zeros((16))
        identical_state = np.zeros((16))
        for i in range(4):
            dup_prior[i * 4 + i] = test.prior_distribution[i]
            identical_state[i * 4 + i] = 1.0
    elif test.Model == 'MG94':
        dup_prior = np.zeros((61**2))
        identical_state = np.zeros((61**2))
        for i in range(61):
            dup_prior[i * 61 + i] = test.prior_distribution[i]
            identical_state[i * 61 + i] = 1.0

    #Pmat_basic, Pmat_geneconv = get_Pmat(mat_basic, mat_geneconv, basic_t)
    if test.tau.size > 1:
        tau_sum = test.tau.sum()
    else:
        tau_sum = 2 * test.tau
    tau_vec = np.array([tau_sum * test.omega if (not pair[0] == pair[1] and isNonsynonymous(pair[0], pair[1], test.codon_table))
                        else 0.0 if pair[0] == pair[1]
                        else tau_sum for i, pair in enumerate(product(test.codon_nonstop, repeat = 2))])
    basic_psame = []
    geneconv_psame = []
    mut_odds = []
    geneconv_odds = []
    for t in basic_t:
        Pmat_basic = np.matrix(expm(mat_basic * t))
        Pmat_geneconv = np.matrix(expm(mat_geneconv * t))
        basic_psame.append(np.matrix(test.prior_distribution) * np.power(Pmat_basic, 2))
        Pdist = np.matrix(dup_prior) * Pmat_geneconv
        geneconv_psame.append(Pdist) # used to store Pdist not exactly psame
        geneconv_chance = (tau_vec * Pdist.T)[0,0]
        geneconv_odds.append(geneconv_chance)
        mut_odds.append( -(mat_geneconv.diagonal() * Pdist.T)[0, 0] - geneconv_chance)
        
    #basic_psame = [np.matrix(test.prior_distribution) * np.power(mat, 2) for mat in Pmat_basic]
    #geneconv_psame = [np.matrix(dup_prior) * mat for mat in Pmat_geneconv]
    
    basic_pdiff = [1 - t.sum() for t in basic_psame]
    geneconv_pdiff = [1.0 - (identical_state * t.T)[0,0] for t in geneconv_psame]
    
    return basic_pdiff, geneconv_pdiff, mut_odds, geneconv_odds

def read_txt(summary_path, paralog, Model, Force = False, clock = False, Dir = False, gBGC = False, Spe = 'Yeast'):
    if not Force:
        prefix_summary = Model + '_'
    else:
        prefix_summary = 'Force_' + Model + '_'

    if Dir:
        prefix_summary = 'Dir_' + prefix_summary

    if gBGC:
        prefix_summary = 'gBGC_' + prefix_summary
        

    if clock:
        suffix_summary = '_clock_summary.txt'
    else:
        suffix_summary = '_nonclock_summary.txt'

    summary_file = summary_path + prefix_summary + '_'.join(paralog) + suffix_summary
    summary = np.genfromtxt(summary_file)

    if Model == 'HKY':
        num_para = 6 # pi_a pi_c pi_g pi_t kappa tau
    elif Model == 'MG94':
        num_para = 7 # pi_a pi_c pi_g pi_t kappa omega tau
    else:
        num_para = 0
        
    if Dir:
        num_para += 1
    if gBGC:
        num_para += 1

    if Spe == 'Yeast':
        num_blen = 12 # add tree branch lengths
    elif Spe == 'Primate':
        num_blen = 8
    # only works for Yeast tree now,  need to change in future
    x_process = [summary[2] + summary[4], summary[2] / (summary[2] + summary[4]), summary[3] / (summary[3] + summary[5])]
    x_process.extend(summary[6:(num_para + 2)])

    if gBGC:
        x_process[:-1] = np.log(x_process[:-1])
    else:
        x_process = np.log(x_process)
    x_rates = np.log(summary[(num_para + 2):(num_para + 2 + num_blen)])
    x = np.concatenate((x_process, x_rates))
    return x
                    

    

