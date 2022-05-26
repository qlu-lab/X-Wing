#!/usr/bin/env python

"""
Gibbs sampler for PANTHER (Polygenic prediction using Bayesian ANnoTation-dependent HorsEshoe Regression)

"""


import scipy as sp
from scipy import linalg 
from scipy import random
import numpy as np


# Full Bayesian Gibbs sampler when the current iteration returns effects from target population
def gibbs_fb_target(phi, snp_dict, beta_mrg, frq_dict, anno_matrix, idx_dict, n, ld_blk, blk_size, n_iter, n_burnin, thin, pop, curr_pop, chrom, out_dir, out_name, seed):
    print('### Run MCMC to obtain SNP posterior effects for %s population ###' % curr_pop)
    print('--- Begin MCMC ')
    if seed != None:
        random.seed(seed)

    # derived stats
    n_pst = (n_iter-n_burnin)/thin
    n_pop = len(pop)
    ind_curr_pop = pop.index(curr_pop)
    p_tot = len(snp_dict['SNP'])

    p = {}
    n_blk = {}
    het = {}
    for pp in range(n_pop):
        p[pp] = len(beta_mrg[pp])
        n_blk[pp] = len(ld_blk[pp])
        het[pp] = sp.sqrt(2.0*frq_dict[pp]*(1.0-frq_dict[pp]))


    K = {}
    n_anno = {}
    for pp in range(n_pop):
        K[pp] = len(anno_matrix[pp]) # K is the number of total annotation included for pp-th population
        n_anno[pp] = len(np.unique(anno_matrix[pp][0]))  # Number of category of k-th annotation for pp-th population

    n_grp = sp.zeros((p_tot,1))
    for jj in range(p_tot):
        for pp in range(n_pop):
            if jj in idx_dict[pp]:
                n_grp[jj] += 1

    # initialization
    beta = {}
    sigma = {}
    for pp in range(n_pop):
        beta[pp] = sp.zeros((p[pp],1))
        sigma[pp] = 1.0


    lamda = {}
    t = {}
    for pp in range(n_pop):
        lamda[pp] = {}
        t[pp] = {}
        for qq in range(n_anno[pp]):
            lamda[pp][qq] = np.ones(K[pp]) # lambda[x][y]  x ranges from 1:n_anno, y rangers from 1:K(groups);
            t[pp][qq] = np.ones(K[pp])

    psi = sp.ones((p_tot,1))
    h_prod = {}
    for pp in range(n_pop):
        h_prod[pp] = sp.ones((p[pp],1))

    if phi == None:
        phi = 1.0; phi_updt = True
    else:
        phi_updt = False

    beta_est = {}
    beta_sq_est = {}
    sigma_est = {}
    for pp in range(n_pop):
        beta_est[pp] = sp.zeros((p[pp],1))
        beta_sq_est[pp] = sp.zeros((p[pp],1))
        sigma_est[pp] = 0.0

    psi_est = sp.zeros((p_tot,1))
    phi_est = 0.0

    # MCMC
    for itr in range(1,n_iter+1):
        if itr % 100 == 0:
            print('... iter-' + str(itr) + ' ...')
        for pp in range(n_pop):
            mm = 0; quad = 0.0
            psi_pp = psi[idx_dict[pp]]
            h_prod_pp = h_prod[pp]
            for kk in range(n_blk[pp]):
                if blk_size[pp][kk] == 0:
                    continue
                else:
                    idx_blk = range(mm,mm+blk_size[pp][kk])
                    diag_tmp = np.array([max(1, x) for x in (1.0/(psi_pp[idx_blk]*phi*h_prod_pp[idx_blk]).T[0])])
                    dinvt = ld_blk[pp][kk] + sp.diag(diag_tmp)
                    dinvt_chol = linalg.cholesky(dinvt)
                    beta_tmp = linalg.solve_triangular(dinvt_chol, beta_mrg[pp][idx_blk], trans='T') \
                                + sp.sqrt(sigma[pp]/n[pp])*random.randn(len(idx_blk),1)
                    beta[pp][idx_blk] = linalg.solve_triangular(dinvt_chol, beta_tmp, trans='N')
                    quad += sp.dot(sp.dot(beta[pp][idx_blk].T, dinvt), beta[pp][idx_blk])
                    mm += blk_size[pp][kk]
            err = max(n[pp]/2.0*(1.0-2.0*sum(beta[pp]*beta_mrg[pp])+quad), n[pp]/2.0*sum(beta[pp]**2/psi_pp/phi/h_prod_pp))
            sigma[pp] = 1.0/random.gamma((n[pp]+p[pp])/2.0, 1.0/err)
        c = 1.0/random.gamma(1.0, 1.0/(1 + 1.0/psi))
        xx = sp.zeros((p_tot,1))
        for pp in range(n_pop):
            xx[idx_dict[pp]] += n[pp]*beta[pp]**2/(2* phi * sigma[pp] * h_prod[pp])
        for jj in range(p_tot):
            while True:
                try:
                    psi[jj] = 1.0/random.gamma(0.5*n_grp[jj] + 0.5, 1.0/(1.0/c[jj] + xx[jj])) # Add *h_prod[jj]
                except:
                    continue
                else:
                    break
        if phi_updt == True:
            v = 1.0/random.gamma(1.0, 1.0/(1 + 1.0/phi))
            zz = sp.zeros((1,1))
            p_all = 0
            for pp in range(n_pop):
                zz += n[pp]*sum(beta[pp]**2/(psi[idx_dict[pp]]*h_prod[pp]))/(2.0 * sigma[pp])
                p_all += p[pp]
            phi = 1.0/random.gamma((p_all+1.0)/2.0, 1.0/(1.0/v + zz))
        psi[psi*phi[0] > 1] = 1.0/phi[0]
        if (itr > n_burnin) and (itr % thin == 0):
            for pp in range(n_pop):
                beta_est[pp] = beta_est[pp] + beta[pp]/n_pst
                beta_sq_est[pp] = beta_sq_est[pp] + beta[pp]**2/n_pst
                sigma_est[pp] = sigma_est[pp] + sigma[pp]/n_pst
            psi_est = psi_est + psi/n_pst
            phi_est = phi_est + phi/n_pst


    # convert standardized beta to per-allele beta
    for pp in range(n_pop):
        beta_est[pp] /= het[pp]
        beta_sq_est[pp] /= het[pp]**2

    # write posterior effect sizes
    pp = ind_curr_pop
    # for pp in range(n_pop):
    if phi_updt == True:
        eff_file = out_dir + '/' + '%s_%s_pst_eff_phiauto_chr%d.txt' % (out_name, pop[pp], chrom)
    else:
        eff_file = out_dir + '/' + '%s_%s_pst_eff__phi%1.0e_chr%d.txt' % (out_name, pop[pp], phi, chrom)

    snp_pp = [snp_dict['SNP'][ii] for ii in idx_dict[pp]]
    bp_pp = [snp_dict['BP'][ii] for ii in idx_dict[pp]]
    a1_pp = [snp_dict['A1'][ii] for ii in idx_dict[pp]]
    a2_pp = [snp_dict['A2'][ii] for ii in idx_dict[pp]]

    with open(eff_file, 'w') as ff:
        for snp, bp, a1, a2, beta in zip(snp_pp, bp_pp, a1_pp, a2_pp, beta_est[pp]):
            ff.write('%d\t%s\t%d\t%s\t%s\t%.6e\n' % (chrom, snp, bp, a1, a2, beta))

    # print estimated phi
    if phi_updt == True:
        print('--- Estimated global shrinkage parameter: %1.2e ' % phi_est )

    print('--- Finish MCMC! ')


# Full Bayesian Gibbs sampler when the current iteration returns effects from non-target population
def gibbs_fb_other(phi, snp_dict, beta_mrg, frq_dict, anno_matrix, idx_dict, n, ld_blk, blk_size, n_iter, n_burnin, thin, pop, curr_pop, chrom, out_dir, out_name, seed):
    print('### Run MCMC to obtain SNP posterior effects for %s population ###' % curr_pop)
    print('--- Begin MCMC ')
    # seed
    if seed != None:
        random.seed(seed)

    # derived stats
    n_pst = (n_iter-n_burnin)/thin
    n_pop = len(pop)
    ind_curr_pop = pop.index(curr_pop)
    p_tot = len(snp_dict['SNP'])

    p = {}
    n_blk = {}
    het = {}
    for pp in range(n_pop):
        p[pp] = len(beta_mrg[pp])
        n_blk[pp] = len(ld_blk[pp])
        het[pp] = sp.sqrt(2.0*frq_dict[pp]*(1.0-frq_dict[pp]))


    K = {}
    n_anno = {}
    for pp in range(n_pop):
        K[pp] = len(anno_matrix[pp]) # K is the number of total annotation included for pp-th population
        n_anno[pp] = len(np.unique(anno_matrix[pp][0]))  # Number of category of k-th annotation for pp-th population

    n_grp = sp.zeros((p_tot,1))
    for jj in range(p_tot):
        for pp in range(n_pop):
            if jj in idx_dict[pp]:
                n_grp[jj] += 1

    # initialization
    beta = {}
    sigma = {}
    for pp in range(n_pop):
        beta[pp] = sp.zeros((p[pp],1))
        sigma[pp] = 1.0


    lamda = {}
    t = {}
    for pp in range(n_pop):
        lamda[pp] = {}
        t[pp] = {}
        for qq in range(n_anno[pp]):
            lamda[pp][qq] = np.ones(K[pp]) # lambda[x][y]  x ranges from 1:n_anno, y rangers from 1:K(groups);
            t[pp][qq] = np.ones(K[pp])

    psi = sp.ones((p_tot,1))
    h_prod = {}
    for pp in range(n_pop):
        h_prod[pp] = sp.ones((p[pp],1))

    if phi == None:
        phi = 1.0; phi_updt = True
    else:
        phi_updt = False

    beta_est = {}
    beta_sq_est = {}
    sigma_est = {}
    for pp in range(n_pop):
        beta_est[pp] = sp.zeros((p[pp],1))
        beta_sq_est[pp] = sp.zeros((p[pp],1))
        sigma_est[pp] = 0.0

    psi_est = sp.zeros((p_tot,1))
    phi_est = 0.0

    lamda_est = {}
    for pp in range(n_pop):
        lamda_est[pp] = {}
        for qq in range(n_anno[pp]):
            lamda_est[pp][qq] = np.zeros(K[pp])

    # MCMC
    for itr in range(1,n_iter+1):
        if itr % 100 == 0:
            print('... iter-' + str(itr) + ' ...')
        for pp in range(n_pop):
            mm = 0; quad = 0.0
            psi_pp = psi[idx_dict[pp]]
            h_prod_pp = h_prod[pp]
            for kk in range(n_blk[pp]):
                if blk_size[pp][kk] == 0:
                    continue
                else:
                    idx_blk = range(mm,mm+blk_size[pp][kk])
                    diag_tmp = np.array([max(1, x) for x in (1.0/(psi_pp[idx_blk]*phi*h_prod_pp[idx_blk]).T[0])])
                    dinvt = ld_blk[pp][kk] + sp.diag(diag_tmp)
                    dinvt_chol = linalg.cholesky(dinvt)
                    beta_tmp = linalg.solve_triangular(dinvt_chol, beta_mrg[pp][idx_blk], trans='T') \
                                + sp.sqrt(sigma[pp]/n[pp])*random.randn(len(idx_blk),1)
                    beta[pp][idx_blk] = linalg.solve_triangular(dinvt_chol, beta_tmp, trans='N')
                    quad += sp.dot(sp.dot(beta[pp][idx_blk].T, dinvt), beta[pp][idx_blk])
                    mm += blk_size[pp][kk]
            err = max(n[pp]/2.0*(1.0-2.0*sum(beta[pp]*beta_mrg[pp])+quad), n[pp]/2.0*sum(beta[pp]**2/psi_pp/phi/h_prod_pp))
            sigma[pp] = 1.0/random.gamma((n[pp]+p[pp])/2.0, 1.0/err)
        c = 1.0/random.gamma(1.0, 1.0/(1 + 1.0/psi))
        xx = sp.zeros((p_tot,1))
        for pp in range(n_pop):
            xx[idx_dict[pp]] += n[pp]*beta[pp]**2/(2* phi * sigma[pp] * h_prod[pp])
        for jj in range(p_tot):
            while True:
                try:
                    psi[jj] = 1.0/random.gamma(0.5*n_grp[jj] + 0.5, 1.0/(1.0/c[jj] + xx[jj])) # Add *h_prod[jj]
                except:
                    continue
                else:
                    break
        if phi_updt == True:
            v = 1.0/random.gamma(1.0, 1.0/(1 + 1.0/phi))
            zz = sp.zeros((1,1))
            p_all = 0
            for pp in range(n_pop):
                zz += n[pp]*sum(beta[pp]**2/(psi[idx_dict[pp]]*h_prod[pp]))/(2.0 * sigma[pp])
                p_all += p[pp]
            phi = 1.0/random.gamma((p_all+1.0)/2.0, 1.0/(1.0/v + zz))
        for pp in range(n_pop):
            for qq in range(n_anno[pp]):
                t[pp][qq] = 1.0/random.gamma(1.0, 1.0/(1 +  1.0/lamda[pp][qq]))
        ## lambda for all SNPs at current annotation. anno_matrix[kk] is the annotation group for all SNPs for current annotation
        for pp in range(n_pop):
                for kk in range(K[pp]):# K is the number of total annotation included for pp-th population; x rangers from 1:n_anno
                    snp_h = np.array([lamda[pp][x][kk] for x in anno_matrix[pp][kk]]).reshape(p[pp], 1) # lamda for all SNPs in the current annotation
                    for qq in range(n_anno[pp]): # qq is the number of category in each annotation
                        lamda_tmp = lamda[pp][qq][kk].copy() # Record lamda before update
                        yy = (1.0 * n[pp] * sum( (beta[pp]**2/(psi[idx_dict[pp]] * h_prod[pp] / snp_h))[anno_matrix[pp][kk] == qq] )/(2*sigma[pp]*phi))[0]
                        y_all = np.count_nonzero(anno_matrix[pp][kk] == qq)
                        lamda[pp][qq][kk] = 1.0/random.gamma( (y_all + 1.0)/2.0, 1.0/(yy + 1.0/t[pp][qq][kk] ))
                        h_prod[pp][anno_matrix[pp][kk] == qq] = h_prod[pp][anno_matrix[pp][kk] == qq]/lamda_tmp*lamda[pp][qq][kk]

        # posterior effects
        if (itr > n_burnin) and (itr % thin == 0):
            for pp in range(n_pop):
                beta_est[pp] = beta_est[pp] + beta[pp]/n_pst
                beta_sq_est[pp] = beta_sq_est[pp] + beta[pp]**2/n_pst
                sigma_est[pp] = sigma_est[pp] + sigma[pp]/n_pst
                for qq in range(n_anno[pp]):
                    lamda_est[pp][qq] = lamda_est[pp][qq] + lamda[pp][qq]/n_pst
            psi_est = psi_est + psi/n_pst
            phi_est = phi_est + phi/n_pst


    # convert standardized beta to per-allele beta
    for pp in range(n_pop):
        beta_est[pp] /= het[pp]
        beta_sq_est[pp] /= het[pp]**2

    # write posterior effect sizes
    pp = ind_curr_pop
    if phi_updt == True:
        eff_file = out_dir + '/' + '%s_%s_pst_eff_phiauto_chr%d.txt' % (out_name, pop[pp], chrom)
    else:
        eff_file = out_dir + '/' + '%s_%s_pst_eff_phi%1.0e_chr%d.txt' % (out_name, pop[pp], phi, chrom)

    snp_pp = [snp_dict['SNP'][ii] for ii in idx_dict[pp]]
    bp_pp = [snp_dict['BP'][ii] for ii in idx_dict[pp]]
    a1_pp = [snp_dict['A1'][ii] for ii in idx_dict[pp]]
    a2_pp = [snp_dict['A2'][ii] for ii in idx_dict[pp]]

    with open(eff_file, 'w') as ff:
        for snp, bp, a1, a2, beta in zip(snp_pp, bp_pp, a1_pp, a2_pp, beta_est[pp]):
            ff.write('%d\t%s\t%d\t%s\t%s\t%.6e\n' % (chrom, snp, bp, a1, a2, beta))

    # print estimated phi
    if phi_updt == True:
        print('--- Estimated global shrinkage parameter: %1.2e ' % phi_est )
    print('--- Finish MCMC! ')


# Tuning parameter Gibbs sampler when the current iteration returns effects from target population
def gibbs_tune_target(phi, snp_dict, beta_mrg, frq_dict, anno_matrix, idx_dict, n, ld_blk, blk_size, n_iter, n_burnin, thin, pop, curr_pop, chrom, out_dir, out_name, seed):
    print('### Run MCMC to obtain SNP posterior effects for %s population ###' % curr_pop)
    print('--- Begin MCMC ')
    if seed != None:
        random.seed(seed)

    # derived stats
    n_pst = (n_iter-n_burnin)/thin
    n_pop = len(pop)
    ind_curr_pop = pop.index(curr_pop)
    p_tot = len(snp_dict['SNP'])

    p = {}
    n_blk = {}
    het = {}
    for pp in range(n_pop):
        p[pp] = len(beta_mrg[pp])
        n_blk[pp] = len(ld_blk[pp])
        het[pp] = sp.sqrt(2.0*frq_dict[pp]*(1.0-frq_dict[pp]))


    K = {}
    n_anno = {}
    for pp in range(n_pop):
        K[pp] = len(anno_matrix[pp]) # K is the number of total annotation included for pp-th population
        n_anno[pp] = len(np.unique(anno_matrix[pp][0]))  # Number of category of k-th annotation for pp-th population

    n_grp = sp.zeros((p_tot,1))
    for jj in range(p_tot):
        for pp in range(n_pop):
            if jj in idx_dict[pp]:
                n_grp[jj] += 1

    # initialization
    beta = {}
    sigma = {}
    for pp in range(n_pop):
        beta[pp] = sp.zeros((p[pp],1))
        sigma[pp] = 1.0


    lamda = {}
    t = {}
    for pp in range(n_pop):
        lamda[pp] = {}
        t[pp] = {}
        for qq in range(n_anno[pp]):
            lamda[pp][qq] = np.ones(K[pp]) # lambda[x][y]  x ranges from 1:n_anno, y rangers from 1:K(groups);
            t[pp][qq] = np.ones(K[pp])

    psi = sp.ones((p_tot,1))
    h_prod = {}
    for pp in range(n_pop):
        h_prod[pp] = sp.ones((p[pp],1))

    if phi == None:
        phi = 1.0; phi_updt = True
    else:
        phi_updt = False

    beta_est = {}
    beta_sq_est = {}
    sigma_est = {}
    for pp in range(n_pop):
        beta_est[pp] = sp.zeros((p[pp],1))
        beta_sq_est[pp] = sp.zeros((p[pp],1))
        sigma_est[pp] = 0.0

    psi_est = sp.zeros((p_tot,1))
    phi_est = 0.0

    # MCMC
    for itr in range(1,n_iter+1):
        if itr % 100 == 0:
            print('... iter-' + str(itr) + ' ...')
        for pp in range(n_pop):
            mm = 0; quad = 0.0
            psi_pp = psi[idx_dict[pp]]
            h_prod_pp = h_prod[pp]
            for kk in range(n_blk[pp]):
                if blk_size[pp][kk] == 0:
                    continue
                else:
                    idx_blk = range(mm,mm+blk_size[pp][kk])
                    diag_tmp = np.array([max(1, x) for x in (1.0/(psi_pp[idx_blk]*phi*h_prod_pp[idx_blk]).T[0])])
                    dinvt = ld_blk[pp][kk] + sp.diag(diag_tmp)
                    dinvt_chol = linalg.cholesky(dinvt)
                    beta_tmp = linalg.solve_triangular(dinvt_chol, beta_mrg[pp][idx_blk], trans='T') \
                                + sp.sqrt(sigma[pp]/n[pp])*random.randn(len(idx_blk),1)
                    beta[pp][idx_blk] = linalg.solve_triangular(dinvt_chol, beta_tmp, trans='N')
                    quad += sp.dot(sp.dot(beta[pp][idx_blk].T, dinvt), beta[pp][idx_blk])
                    mm += blk_size[pp][kk]
            err = max(n[pp]/2.0*(1.0-2.0*sum(beta[pp]*beta_mrg[pp])+quad), n[pp]/2.0*sum(beta[pp]**2/psi_pp/phi/h_prod_pp))
            sigma[pp] = 1.0/random.gamma((n[pp]+p[pp])/2.0, 1.0/err)
        c = 1.0/random.gamma(1.0, 1.0/(1 + 1.0/psi))
        xx = sp.zeros((p_tot,1))
        for pp in range(n_pop):
            xx[idx_dict[pp]] += n[pp]*beta[pp]**2/(2* phi * sigma[pp] * h_prod[pp])
        for jj in range(p_tot):
            while True:
                try:
                    psi[jj] = 1.0/random.gamma(0.5*n_grp[jj] + 0.5, 1.0/(1.0/c[jj] + xx[jj])) # Add *h_prod[jj]
                except:
                    continue
                else:
                    break
        if phi_updt == True:
            v = 1.0/random.gamma(1.0, 1.0/(1 + 1.0/phi))
            zz = sp.zeros((1,1))
            p_all = 0
            for pp in range(n_pop):
                zz += n[pp]*sum(beta[pp]**2/(psi[idx_dict[pp]]*h_prod[pp]))/(2.0 * sigma[pp])
                p_all += p[pp]
            phi = 1.0/random.gamma((p_all+1.0)/2.0, 1.0/(1.0/v + zz))
        # print(phi)
        # print(h_prod)
        psi[psi*phi > 1] = 1.0/phi
        if (itr > n_burnin) and (itr % thin == 0):
            for pp in range(n_pop):
                beta_est[pp] = beta_est[pp] + beta[pp]/n_pst
                beta_sq_est[pp] = beta_sq_est[pp] + beta[pp]**2/n_pst
                sigma_est[pp] = sigma_est[pp] + sigma[pp]/n_pst
            psi_est = psi_est + psi/n_pst
            phi_est = phi_est + phi/n_pst


    # convert standardized beta to per-allele beta
    for pp in range(n_pop):
        beta_est[pp] /= het[pp]
        beta_sq_est[pp] /= het[pp]**2

    # write posterior effect sizes
    pp = ind_curr_pop
    # for pp in range(n_pop):
    if phi_updt == True:
        eff_file = out_dir + '/' + '%s_%s_pst_eff_phiauto_chr%d.txt' % (out_name, pop[pp], chrom)
    else:
        eff_file = out_dir + '/' + '%s_%s_pst_eff_phi%1.0e_chr%d.txt' % (out_name, pop[pp], phi, chrom)

    snp_pp = [snp_dict['SNP'][ii] for ii in idx_dict[pp]]
    bp_pp = [snp_dict['BP'][ii] for ii in idx_dict[pp]]
    a1_pp = [snp_dict['A1'][ii] for ii in idx_dict[pp]]
    a2_pp = [snp_dict['A2'][ii] for ii in idx_dict[pp]]

    with open(eff_file, 'w') as ff:
        for snp, bp, a1, a2, beta in zip(snp_pp, bp_pp, a1_pp, a2_pp, beta_est[pp]):
            ff.write('%d\t%s\t%d\t%s\t%s\t%.6e\n' % (chrom, snp, bp, a1, a2, beta))

    # print estimated phi
    if phi_updt == True:
        print('--- Estimated global shrinkage parameter: %1.2e ' % phi_est )

    print('--- Finish MCMC! ')


# Tuning parameter Gibbs sampler when the current iteration returns effects from non-target population
def gibbs_tune_other(phi, snp_dict, beta_mrg, frq_dict, anno_matrix, idx_dict, n, ld_blk, blk_size, n_iter, n_burnin, thin, pop, curr_pop, chrom, out_dir, out_name, seed):
    print('### Run MCMC to obtain SNP posterior effects for %s population ###' % curr_pop)
    print('--- Begin MCMC ')
    # seed
    if seed != None:
        random.seed(seed)

    # derived stats
    n_pst = (n_iter-n_burnin)/thin
    n_pop = len(pop)
    ind_curr_pop = pop.index(curr_pop)
    p_tot = len(snp_dict['SNP'])

    p = {}
    n_blk = {}
    het = {}
    for pp in range(n_pop):
        p[pp] = len(beta_mrg[pp])
        n_blk[pp] = len(ld_blk[pp])
        het[pp] = sp.sqrt(2.0*frq_dict[pp]*(1.0-frq_dict[pp]))


    K = {}
    n_anno = {}
    for pp in range(n_pop):
        K[pp] = len(anno_matrix[pp]) # K is the number of total annotation included for pp-th population
        n_anno[pp] = len(np.unique(anno_matrix[pp][0]))  # Number of category of k-th annotation for pp-th population

    n_grp = sp.zeros((p_tot,1))
    for jj in range(p_tot):
        for pp in range(n_pop):
            if jj in idx_dict[pp]:
                n_grp[jj] += 1

    # initialization
    beta = {}
    sigma = {}
    for pp in range(n_pop):
        beta[pp] = sp.zeros((p[pp],1))
        sigma[pp] = 1.0


    lamda = {}
    t = {}
    for pp in range(n_pop):
        lamda[pp] = {}
        t[pp] = {}
        for qq in range(n_anno[pp]):
            lamda[pp][qq] = np.ones(K[pp]) # lambda[x][y]  x ranges from 1:n_anno, y rangers from 1:K(groups);
            t[pp][qq] = np.ones(K[pp])

    psi = sp.ones((p_tot,1))
    h_prod = {}
    for pp in range(n_pop):
        h_prod[pp] = sp.ones((p[pp],1))

    if phi == None:
        phi = 1.0; phi_updt = True
    else:
        phi_updt = False

    beta_est = {}
    beta_sq_est = {}
    sigma_est = {}
    for pp in range(n_pop):
        beta_est[pp] = sp.zeros((p[pp],1))
        beta_sq_est[pp] = sp.zeros((p[pp],1))
        sigma_est[pp] = 0.0

    psi_est = sp.zeros((p_tot,1))
    phi_est = 0.0

    lamda_est = {}
    for pp in range(n_pop):
        lamda_est[pp] = {}
        for qq in range(n_anno[pp]):
            lamda_est[pp][qq] = np.zeros(K[pp])

    # MCMC
    for itr in range(1,n_iter+1):
        if itr % 100 == 0:
            print('... iter-' + str(itr) + ' ...')
        for pp in range(n_pop):
            mm = 0; quad = 0.0
            psi_pp = psi[idx_dict[pp]]
            h_prod_pp = h_prod[pp]
            for kk in range(n_blk[pp]):
                if blk_size[pp][kk] == 0:
                    continue
                else:
                    idx_blk = range(mm,mm+blk_size[pp][kk])
                    diag_tmp = np.array([max(1, x) for x in (1.0/(psi_pp[idx_blk]*phi*h_prod_pp[idx_blk]).T[0])])
                    dinvt = ld_blk[pp][kk] + sp.diag(diag_tmp)
                    dinvt_chol = linalg.cholesky(dinvt)
                    beta_tmp = linalg.solve_triangular(dinvt_chol, beta_mrg[pp][idx_blk], trans='T') \
                                + sp.sqrt(sigma[pp]/n[pp])*random.randn(len(idx_blk),1)
                    beta[pp][idx_blk] = linalg.solve_triangular(dinvt_chol, beta_tmp, trans='N')
                    quad += sp.dot(sp.dot(beta[pp][idx_blk].T, dinvt), beta[pp][idx_blk])
                    mm += blk_size[pp][kk]
            err = max(n[pp]/2.0*(1.0-2.0*sum(beta[pp]*beta_mrg[pp])+quad), n[pp]/2.0*sum(beta[pp]**2/psi_pp/phi/h_prod_pp))
            sigma[pp] = 1.0/random.gamma((n[pp]+p[pp])/2.0, 1.0/err)
        c = 1.0/random.gamma(1.0, 1.0/(1 + 1.0/psi))
        xx = sp.zeros((p_tot,1))
        for pp in range(n_pop):
            xx[idx_dict[pp]] += n[pp]*beta[pp]**2/(2* phi * sigma[pp] * h_prod[pp])
        for jj in range(p_tot):
            while True:
                try:
                    psi[jj] = 1.0/random.gamma(0.5*n_grp[jj] + 0.5, 1.0/(1.0/c[jj] + xx[jj])) # Add *h_prod[jj]
                except:
                    continue
                else:
                    break
        if phi_updt == True:
            v = 1.0/random.gamma(1.0, 1.0/(1 + 1.0/phi))
            zz = sp.zeros((1,1))
            p_all = 0
            for pp in range(n_pop):
                zz += n[pp]*sum(beta[pp]**2/(psi[idx_dict[pp]]*h_prod[pp]))/(2.0 * sigma[pp])
                p_all += p[pp]
            phi = 1.0/random.gamma((p_all+1.0)/2.0, 1.0/(1.0/v + zz))
        for pp in range(n_pop):
            for qq in range(n_anno[pp]):
                t[pp][qq] = 1.0/random.gamma(1.0, 1.0/(1 +  1.0/lamda[pp][qq]))
        ## lambda for all SNPs at current annotation. anno_matrix[kk] is the annotation group for all SNPs for current annotation
        for pp in range(n_pop):
            if pp == ind_curr_pop:
                for kk in range(K[pp]): # K is the number of total annotation; x rangers from 1:n_anno
                    snp_h = np.array([lamda[pp][x][kk] for x in anno_matrix[pp][kk]]).reshape(p[pp], 1) # lamda for all SNPs in the current annotation
                    for qq in range(n_anno[pp]): # qq is the number of category in each annotation
                        lamda_tmp = lamda[pp][qq][kk] # Record lamda before update
                        yy = (1.0 * n[pp] * sum( (beta[pp]**2/(psi[idx_dict[pp]] * h_prod[pp] / snp_h))[anno_matrix[pp][kk] == qq] )/(2*sigma[pp]*phi))[0]
                        y_all = np.count_nonzero(anno_matrix[pp][kk] == qq)
                        lamda[pp][qq][kk] = 1.0/random.gamma( (y_all + 1.0)/2.0, 1.0/(yy + 1.0/t[pp][qq][kk] ))
                        h_prod[pp][anno_matrix[pp][kk] == qq] = h_prod[pp][anno_matrix[pp][kk] == qq]/lamda_tmp*lamda[pp][qq][kk]
            else:
                h_prod[pp] = sp.ones((p[pp],1))

        # posterior effects
        if (itr > n_burnin) and (itr % thin == 0):
            for pp in range(n_pop):
                beta_est[pp] = beta_est[pp] + beta[pp]/n_pst
                beta_sq_est[pp] = beta_sq_est[pp] + beta[pp]**2/n_pst
                sigma_est[pp] = sigma_est[pp] + sigma[pp]/n_pst
                for qq in range(n_anno[pp]):
                    lamda_est[pp][qq] = lamda_est[pp][qq] + lamda[pp][qq]/n_pst
            psi_est = psi_est + psi/n_pst
            phi_est = phi_est + phi/n_pst


    # convert standardized beta to per-allele beta
    for pp in range(n_pop):
        beta_est[pp] /= het[pp]
        beta_sq_est[pp] /= het[pp]**2

    # write posterior effect sizes
    pp = ind_curr_pop
    if phi_updt == True:
        eff_file = out_dir + '/' + '%s_%s_pst_eff_phiauto_chr%d.txt' % (out_name, pop[pp], chrom)
    else:
        eff_file = out_dir + '/' + '%s_%s_pst_eff_phi%1.0e_chr%d.txt' % (out_name, pop[pp], phi, chrom)

    snp_pp = [snp_dict['SNP'][ii] for ii in idx_dict[pp]]
    bp_pp = [snp_dict['BP'][ii] for ii in idx_dict[pp]]
    a1_pp = [snp_dict['A1'][ii] for ii in idx_dict[pp]]
    a2_pp = [snp_dict['A2'][ii] for ii in idx_dict[pp]]

    with open(eff_file, 'w') as ff:
        for snp, bp, a1, a2, beta in zip(snp_pp, bp_pp, a1_pp, a2_pp, beta_est[pp]):
            ff.write('%d\t%s\t%d\t%s\t%s\t%.6e\n' % (chrom, snp, bp, a1, a2, beta))

    # print estimated phi
    if phi_updt == True:
        print('--- Estimated global shrinkage parameter: %1.2e ' % phi_est )
    print('--- Finish MCMC! ')

