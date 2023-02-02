#!/usr/bin/env python

"""
PANTHER(Polygenic prediction using Bayesian ANnoTation-dependent HorsEshoe Regression)

 - The detailed usage of PANTHER can be found in https://github.com/qlu-lab/X-Wing/wiki/1.-Improve-Cross%E2%80%90Population-polygenic-prediction-using-X%E2%80%90Wing
 - A read example of PANTHER can be found in https://github.com/qlu-lab/X-Wing/wiki/2.-Real-example-for-X%E2%80%90Wing

Usage:
python PANTHER.py \
--ref_dir PATH_TO_REFERENCE \
--bim_prefix PATH_TO_TARGET \
--sumstats PATH_TO_SUMSTATS \
--n_gwas GWAS_SAMPLE_SIZE \
--anno_file PATH_TO_ANNOTATION_MATRIX \
--pop POPULATION \
--target TARGET_POPULATION \
--out_name NAME_TO_OUTFILE \
--out_dir PATH_TO_OUTFILE \
--phi PHI \
--n_iter MCMC_ITERATIONS \
--n_burnin MCMC_BURIN \
--thin MCMC_THIN \
--chrom CHR \
--pst_pop OUT_PST_POPULATION \
--seed SEED

where the inputs in order are

 - PATH_TO_REFERENCE (required): Full path to the directory that contains the pre-computed LD matrix and SNP information.
 - PATH_TO_TARGET (required): Full path and the prefix of the bim file for target dataset, which is to provide the SNP list present in target dataset.
 - PATH_TO_SUMSTATS (required): Full paths to multiple GWAS summary statistics, separated by comma. The summary statistics files are in the format of:
        CHR	  SNP	       BP	    A1	A2	BETA	P
        1	  rs12238997   693731	G	A	-0.2768	0.7819
        1	  rs58276399   731718	C	T	-0.4127	0.6797
        1	  rs61770163   732032	C	A	-0.0565	0.9549
        ...
 - GWAS_SAMPLE_SIZE (required): GWAS sample size, in the same order of the GWAS summary statistics files, separated by comma.
 - PATH_TO_ANNOTATION_MATRIX (required): Full path to the multiple annotation matrix, in the same order of the GWAS summary statistics files, separated by comma. The annotation matrix files must have the following format (including the header line and order of the columns):
        CHR SNP A1 A2 Anno
        22 rs2186521 T C 0
        22 rs4911642 C T 0
        22 rs7287144 G A 0
        ...
 - POPULATION (required): Population of the GWAS sample, in the same order of the GWAS summary statistics files, separated by comma. AFR, AMR, EAS, EUR and SAS are allowed.
 - NAME_TO_OUTFILE (required): Full directory to the output posterior mean effect size estimates.
 - PATH_TO_OUTFILE (required): Prefix of the output posterior mean effect size estimates.
 - PHI (optional): Global shrinkage parameter. It you don't specify it, PANTHER will learn it from the data by fully Bayesian approach. Or you do a grid search (e.g., phi=1e-6, 1e-4, 1e-2, 1) to find the phi that maximized the predictive performance in validation set.
 - MCMC_ITERATIONS (optional): Number of MCMC sampling iterations. Default is 1,000 * the number of input GWAS summary statistics.
 - MCMC_BURIN (optional): Number of MCMC burnin iterations. Default is 500 * the number of input GWAS summary statistics.
 - MCMC_THIN (optional): Markov chain thinning factor. Default is 5.
 - CHR (optional): Chromosome fitted by the model, separated by comma, e.g., --chrom=1,3. We highly recommend parallel computation for the 22 autosomes. Default is iterating through 22 autosomes (may be time-consuming).
 - SEED (optional): Non-negative integer as seeds for random number generator.
"""

import os
import sys
import getopt

import munge_data
import gibbs_sampler
import numpy as np
import copy

TopHEAD = "*********************************************************************\n"
TopHEAD += "* Cross-population Weighting (X-Wing) \n"
TopHEAD += "* Version 1.0.0 \n"
TopHEAD += "* Step2: Polygenic prediction using Bayesian annotation-dependent horseshoe regression (PANTHER) \n"
TopHEAD += "* (C) Jiacheng Miao and Hanmin Guo \n"
TopHEAD += "* University of Wisconsin-Madison and Tsinghua University \n"
TopHEAD += "*  https://github.com/qlu-lab/X-Wing \n"
TopHEAD += "* GNU General Public License v3\n"
TopHEAD += "*********************************************************************\n"


def munge_param():
    long_opts_list = ['ref_dir=', 'bim_prefix=', 'anno_file=', 'sumstats=', 'phi=', 'n_gwas=', 'pop=', 'target_pop=',
                      'n_iter=', 'n_burnin=', 'thin=', 'out_dir=', 'out_name=', 'chrom=', 'pst_pop=','seed=', 'help']

    param_dict = {'ref_dir': None, 'bim_prefix': None, 'anno_file': None, 'sumstats': None,  'phi': None, 'n_gwas': None, 'pop': None, 'target_pop': None,
                  'n_iter': None, 'n_burnin': None, 'thin': 5, 'out_dir': None, 'out_name': None, 'chrom': range(1,23), 'pst_pop': None, 'seed': None}

    print('\n')

    if len(sys.argv) > 1:
        try:
            opts, args = getopt.getopt(sys.argv[1:], "h", long_opts_list)
        except:
            print('* Option not recognized.')
            print('* Use --help for usage information.\n')
            sys.exit(2)

        for opt, arg in opts:
            if opt == "-h" or opt == "--help":
                print(__doc__)
                sys.exit(0)
            elif opt == "--ref_dir": param_dict['ref_dir'] = arg
            elif opt == "--bim_prefix": param_dict['bim_prefix'] = arg
            elif opt == "--anno_file": param_dict['anno_file'] = arg.split(',')
            elif opt == "--sumstats": param_dict['sumstats'] = arg.split(',')
            elif opt == "--phi": param_dict['phi'] = float(arg)
            elif opt == "--n_gwas": param_dict['n_gwas'] = list(map(int,arg.split(',')))
            elif opt == "--pop": param_dict['pop'] = arg.split(',')
            elif opt == "--target_pop": param_dict['target_pop'] = arg
            elif opt == "--n_iter": param_dict['n_iter'] = int(arg)
            elif opt == "--n_burnin": param_dict['n_burnin'] = int(arg)
            elif opt == "--thin": param_dict['thin'] = int(arg)
            elif opt == "--out_dir": param_dict['out_dir'] = arg
            elif opt == "--out_name": param_dict['out_name'] = arg
            elif opt == "--chrom": param_dict['chrom'] = arg.split(',')
            elif opt == "--pst_pop": param_dict['pst_pop'] = arg
            elif opt == "--seed": param_dict['seed'] = int(arg)
    else:
        print(__doc__)
        sys.exit(0)


    if param_dict['ref_dir'] == None:
        print('* Must specify the directory to the reference panel using --ref_dir\n')
        sys.exit(2)
    elif param_dict['bim_prefix'] == None:
        print('* Must specify the directory and prefix of the bim file for the target data using --bim_prefix\n')
        sys.exit(2)
    elif param_dict['sumstats'] == None:
        print('* Must provide at least one summary statistics file using --sumstats\n')
        sys.exit(2)
    elif param_dict['n_gwas'] == None:
        print('* Must provide the sample size of the GWAS using --n_gwas\n')
        sys.exit(2)
    elif param_dict['pop'] == None:
        print('* Must specify the population of the GWAS sample using --pop\n')
        sys.exit(2)
    elif param_dict['out_dir'] == None:
        print('* Must specify the output directory using --out_dir\n')
        sys.exit(2)
    elif param_dict['out_name'] == None:
        print('* Must specify the prefix of the output file using --out_name\n')
        sys.exit(2)
    elif (len(param_dict['sumstats']) != len(param_dict['n_gwas']) or 
          len(param_dict['sumstats']) != len(param_dict['pop'])):
        print('* Length of sumstats, n_gwas and pop does not match\n')
        sys.exit(2)

    n_pop = len(param_dict['pop'])
    if param_dict['n_iter'] == None or param_dict['n_burnin'] == None:
        param_dict['n_iter'] = n_pop*1000
        param_dict['n_burnin'] = n_pop*500
    print('Options in effect:')
    print('python PANTHER.py')
    for key in param_dict:
        print('--%s=%s \\' % (key, param_dict[key]))

    print('\n')
    return param_dict


def main():
    header = TopHEAD    
    print(header)
    param_dict = munge_param()
    n_pop = len(param_dict['pop'])

    for chrom in param_dict['chrom']:
        print('### process chromosome %d ###' % int(chrom))
        n_pop = len(param_dict['pop'])
        if os.path.isfile(param_dict['ref_dir'] + '/snpinfo_mult_1kg_hm3'):
            ref = '1kg'
            ref_dict = munge_data.munge_ref(param_dict['ref_dir'] + '/snpinfo_mult_1kg_hm3', int(chrom), ref)
        elif os.path.isfile(param_dict['ref_dir'] + '/snpinfo_mult_ukbb_hm3'):
            ref = 'ukbb'
            ref_dict = munge_data.munge_ref(param_dict['ref_dir'] + '/snpinfo_mult_ukbb_hm3', int(chrom), ref)

        # Read the reference panel
        vld_dict = munge_data.munge_bim(param_dict['bim_prefix'], int(chrom))

        # Read the Anootation
        anno_dict = {}
        for pp in range(n_pop):
            anno_dict[pp] = munge_data.munge_anno(param_dict['anno_file'][pp], int(chrom), param_dict['pop'][pp])


        # Read the summary statistics and overlap with reference panel, annotation matrix, and validation set
        sst_dict = {}
        for pp in range(n_pop):
            sst_dict[pp] = munge_data.munge_sumstats(ref_dict, vld_dict, anno_dict[pp], param_dict['sumstats'][pp], param_dict['pop'][pp], param_dict['n_gwas'][pp])

        ld_blk = {}
        blk_size = {}
        for pp in range(n_pop):
            ld_blk[pp], blk_size[pp] = munge_data.munge_ldblk(param_dict['ref_dir'], sst_dict[pp], param_dict['pop'][pp], int(chrom), ref)

        snp_dict, beta_dict, frq_dict, idx_dict = munge_data.align_ldblk(ref_dict, vld_dict, sst_dict, n_pop, int(chrom))

        # Read the annotation matrix
        anno_matrix = {}
        for pp in range(n_pop):
            anno_matrix[pp] = munge_data.munge_anno_matrix(param_dict['anno_file'][pp], anno_dict[pp], snp_dict, int(chrom), idx_dict[pp])

        
        if param_dict['phi'] is None: # Full Bayesian
            for i in range(len(param_dict['pop'])):
                if (param_dict['pop'][i] == param_dict['target_pop'] and param_dict['pop'][i] == param_dict['pst_pop']): # Return posterior effects for Target population
                    gibbs_sampler.gibbs_fb_target(param_dict['phi'], snp_dict, beta_dict, frq_dict, anno_matrix, idx_dict, param_dict['n_gwas'], ld_blk, blk_size,
                    param_dict['n_iter'], param_dict['n_burnin'], param_dict['thin'], param_dict['pop'], param_dict['pop'][i], int(chrom),
                    param_dict['out_dir'], param_dict['out_name'], param_dict['seed'])
                else:  # Return posterior effects for non-Target population
                    if (param_dict['pop'][i] == param_dict['pst_pop']):
                        anno_matrix_tmp = copy.deepcopy(anno_matrix)
                        for j in range(len(param_dict['pop'])):
                            if j != i:
                                anno_matrix_tmp[j][0] = np.zeros(len(anno_matrix[j][0]), dtype=int)
                        gibbs_sampler.gibbs_fb_other(param_dict['phi'], snp_dict, beta_dict, frq_dict, anno_matrix_tmp, idx_dict, param_dict['n_gwas'], ld_blk, blk_size,
                        param_dict['n_iter'], param_dict['n_burnin'], param_dict['thin'], param_dict['pop'], param_dict['pop'][i], int(chrom),
                        param_dict['out_dir'], param_dict['out_name'], param_dict['seed'])
        else: # Tuning parameter
            for i in range(len(param_dict['pop'])): 
                if (param_dict['pop'][i] == param_dict['target_pop'] and param_dict['pop'][i] == param_dict['pst_pop']): # Return posterior effects for Target population
                    gibbs_sampler.gibbs_tune_target(param_dict['phi'], snp_dict, beta_dict, frq_dict, anno_matrix, idx_dict, param_dict['n_gwas'], ld_blk, blk_size,
                    param_dict['n_iter'], param_dict['n_burnin'], param_dict['thin'], param_dict['pop'], param_dict['pop'][i], int(chrom),
                    param_dict['out_dir'], param_dict['out_name'], param_dict['seed'])
                else:  # Return posterior effects for non-Target population
                    if (param_dict['pop'][i] == param_dict['pst_pop']):
                        anno_matrix_tmp = copy.deepcopy(anno_matrix)
                        for j in range(len(param_dict['pop'])):
                            if j != i:
                                anno_matrix_tmp[j][0] = np.zeros(len(anno_matrix[j][0]), dtype=int)
                        gibbs_sampler.gibbs_tune_other(param_dict['phi'], snp_dict, beta_dict, frq_dict, anno_matrix_tmp, idx_dict, param_dict['n_gwas'], ld_blk, blk_size,
                        param_dict['n_iter'], param_dict['n_burnin'], param_dict['thin'], param_dict['pop'], param_dict['pop'][i], int(chrom),
                        param_dict['out_dir'], param_dict['out_name'], param_dict['seed'])


        print('\n')
        print('### Finish PANTHER! ###')


if __name__ == '__main__':
    main()


