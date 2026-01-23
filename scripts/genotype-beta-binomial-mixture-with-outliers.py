# adapted from https://pyro.ai/examples/gmm.html

import argparse
import logging

import numpy as np
import allel
import matplotlib
import matplotlib.pyplot as plt
from collections import defaultdict
import tqdm
import pandas as pd
from sklearn.metrics import confusion_matrix

import torch
import pyro
import pyro.distributions as dist
from pyro import poutine
from pyro.infer.autoguide import AutoDelta
from pyro.optim import ClippedAdam
from pyro.infer import SVI, TraceEnum_ELBO, config_enumerate

logger = logging.getLogger('genotype-beta-binomial-mixture')
logging.basicConfig(level=logging.INFO)

def read_data(kanpig_vcf):
    callset = allel.read_vcf(kanpig_vcf, fields='*')

    # T indexes all records, t indexes subset with non-missing AD below
    chrom_T = callset['variants/CHROM'][:]
    pos_T = callset['variants/POS'][:]
    id_T = callset['variants/ID'][:]
    ref_T = callset['variants/REF'][:]
    alt_T = callset['variants/ALT'][:, 0]
    ref_count_T = callset['calldata/AD'][:, 0, 0]
    alt_count_T = callset['calldata/AD'][:, 0, 1]
    gt_Tp = callset['calldata/GT'][:, 0, :]
    gq_T = callset['calldata/GQ'][:, 0]
    sq_T = callset['calldata/SQ'][:, 0]
    sample = callset['samples'][0]

    return chrom_T, pos_T, id_T, ref_T, alt_T, ref_count_T, alt_count_T, gt_Tp, gq_T, sq_T, sample

@config_enumerate
def model(hyperparameters, total_count_t, alt_count_t=None, z_t=None, y_t=None):
    alpha_pi_k, alpha_pi_outlier, beta_pi_outlier, alpha_mu_k, beta_mu_k, lambda_nu_k = hyperparameters

    pi_k = pyro.sample('pi_k', dist.Dirichlet(alpha_pi_k))
    pi_outlier = pyro.sample('pi_outlier', dist.Beta(alpha_pi_outlier, beta_pi_outlier))

    with pyro.plate('components', alpha_pi_k.shape[0]):
        mu_k = pyro.sample('mu_k', dist.Beta(alpha_mu_k, beta_mu_k))
        nu_k = pyro.sample('nu_k', dist.Exponential(lambda_nu_k))
        alpha_k = mu_k * nu_k
        beta_k = (1. - mu_k) * nu_k

    with pyro.plate('loci', total_count_t.shape[0]):
        z_t = pyro.sample('z_t', dist.Categorical(pi_k), obs=z_t)
        y_t = pyro.sample('y_t', dist.Bernoulli(pi_outlier), obs=y_t)
        alt_count_t = pyro.sample('alt_count_t', dist.BetaBinomial(
            y_t.float() + (1 - y_t) * alpha_k[z_t], 
            y_t.float() + (1 - y_t) * beta_k[z_t], total_count_t), obs=alt_count_t)
    
    return pi_k, pi_outlier, mu_k, nu_k, z_t, y_t, alt_count_t

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--kanpig-vcf',
                        type=str,
                        required=True,
                        help='Input single-sample Kanpig VCF containing AD values.')

    parser.add_argument('--alpha-weights',
                        type=float,
                        nargs='*',
                        default=[1., 1., 1.],
                        help='Dirichlet concentration hyperparameter for component weights.')

    parser.add_argument('--alpha-outlier',
                        type=float,
                        nargs=1,
                        default=1.,
                        help='Beta alpha hyperparameter for outlier probability.')

    parser.add_argument('--beta-outlier',
                        type=float,
                        nargs=1,
                        default=100.,
                        help='Beta beta hyperparameter for outlier probability.')

    parser.add_argument('--alpha-means',
                        type=float,
                        nargs='*',
                        default=[1., 10., 10.],
                        help='Beta alpha hyperparameter for component means.')

    parser.add_argument('--beta-means',
                        type=float,
                        nargs='*',
                        default=[10., 10., 1.],
                        help='Beta beta hyperparameter for component means.')

    parser.add_argument('--lambda-precisions',
                        type=float,
                        nargs='*',
                        default=[0.01, 0.01, 0.01],
                        help='Exponential lambda hyperparameter for component precisions.')

    parser.add_argument('--n-init',
                        type=int,
                        default=100,
                        help='Number of initialization seeds.')

    parser.add_argument('--n-iter',
                        type=int,
                        default=500,
                        help='Number of iterations.')

    parser.add_argument('--learning-rate',
                        type=float,
                        default=0.01,
                        help='Learning rate for clipped Adam optimizer.')

    parser.add_argument('--n-count-bins',
                        type=int,
                        default=50,
                        help='Number of count bins for plotting histograms.')

    parser.add_argument('--max-qual',
                        type=int,
                        default=400,
                        help='Maximum GQ/SQ (can mitigate errors due to unsafe logs).')

    parser.add_argument('--output-prefix',
                        type=str,
                        required=True,
                        help='Prefix for output files.')
    
    parser.add_argument('--print-graphs',
                        type=int,
                        default=0,
                        help='Saves plots to disk.')

    args = parser.parse_args()

    # set hyperparameters
    assert len(set([len(args.alpha_weights), len(args.alpha_means), len(args.beta_means), len(args.lambda_precisions)])) == 1 # same number of components
    alpha_pi_k = torch.tensor(args.alpha_weights, dtype=torch.float32)
    alpha_pi_outlier = torch.tensor(args.alpha_outlier, dtype=torch.float32)
    beta_pi_outlier = torch.tensor(args.beta_outlier, dtype=torch.float32)
    alpha_mu_k = torch.tensor(args.alpha_means, dtype=torch.float32)
    beta_mu_k = torch.tensor(args.beta_means, dtype=torch.float32)
    lambda_nu_k = torch.tensor(args.lambda_precisions, dtype=torch.float32)
    n_components = alpha_pi_k.shape[0]
    hyperparameters = alpha_pi_k, alpha_pi_outlier, beta_pi_outlier, alpha_mu_k, beta_mu_k, lambda_nu_k

    # read data
    logger.info('Reading data...')
    chrom_T, pos_T, id_T, ref_T, alt_T, ref_count_T, alt_count_T, gt_Tp, gq_T, sq_T, sample = read_data(args.kanpig_vcf)
    missing_T = (ref_count_T == -1) | (alt_count_T == -1)
    ref_count_t = torch.tensor(ref_count_T[~missing_T])
    alt_count_t = torch.tensor(alt_count_T[~missing_T])
    total_count_t = ref_count_t + alt_count_t

    # initialize and set up SVI
    logger.info('Initializing...')
    optim = ClippedAdam({'lr': args.learning_rate})
    elbo = TraceEnum_ELBO(max_plate_nesting=1)

    def initialize(seed):
        global global_guide, svi
        pyro.set_rng_seed(seed)
        pyro.clear_param_store()
        global_guide = AutoDelta(
            poutine.block(model, expose=['pi_k', 'pi_outlier', 'mu_k', 'nu_k'])
        )
        svi = SVI(model, global_guide, optim, loss=elbo)
        return svi.loss(model, global_guide, hyperparameters, total_count_t, alt_count_t)

    loss, seed = min((initialize(seed), seed) for seed in tqdm.tqdm(range(args.n_init)))
    initialize(seed)
    logger.info(f'Initialized: seed = {seed}, loss = {loss}')

    gradient_norms = defaultdict(list)
    for name, value in pyro.get_param_store().named_parameters():
        value.register_hook(lambda g, name=name: gradient_norms[name].append(g.norm().item()))

    # perform inference
    logger.info('Performing inference...')
    losses = []
    for i in tqdm.tqdm(range(args.n_iter)):
        loss = svi.step(hyperparameters, total_count_t, alt_count_t)
        losses.append(loss)

    # print MAP estimates
    logger.info('MAP parameter estimates:')
    map_estimates = global_guide()
    for param in ['pi_k', 'pi_outlier', 'mu_k', 'nu_k']:
        param_map = map_estimates[param].data.numpy()
        logger.info(f'{param} = {param_map}')

    # calculate GQ and SQ
    data_ti = torch.tensor([[d, a, k, y] for d, a in zip(total_count_t, alt_count_t) for k in range(n_components) for y in [0, 1]])
    map_model = poutine.condition(model, data=global_guide())
    map_model_tr = poutine.trace(map_model).get_trace(hyperparameters, data_ti[:, 0], data_ti[:, 1], data_ti[:, 2], data_ti[:, 3].float())
    map_model_tr.compute_log_prob()
    loglike_tky = (map_model_tr.nodes['z_t']['log_prob'] + map_model_tr.nodes['y_t']['log_prob'] + map_model_tr.nodes['alt_count_t']['log_prob']).data.numpy().reshape(-1, n_components, 2)
    loglike_tk = loglike_tky.sum(axis=-1)
    gt_t = np.argmax(loglike_tk, axis=1)
    eps = 10.**(-args.max_qual / 10.)
    gq_t = 10. * np.diff(np.sort(np.log10(eps + np.exp(loglike_tk)), axis=1), n=1, axis=1)[:, -1]
    sq_t = -10. * np.log10(eps + np.exp(loglike_tk[:, 0] - torch.logsumexp(torch.tensor(loglike_tk), dim=1).data.numpy()))
    oq_t = 10 * np.diff(np.log10(eps + np.exp(loglike_tky.sum(axis=1))), n=1, axis=1)[:, -1]

    # output variant records to TSV; cat with the input VCF header to create a VCF suitable for use with bcftools annotate
    # in addition to sites with missing AD, we further exclude sites with zero depth to retain the original Kanpig flat priors/posteriors and GQ = 0
    logger.info('Writing VCF records for annotation...')
    output_df = pd.DataFrame()
    output_T = ~missing_T & (ref_count_T + alt_count_T > 0)
    output_t = ref_count_t + alt_count_t > 0
    output_df['CHROM'] = chrom_T[output_T]
    output_df['POS'] = pos_T[output_T]
    output_df['ID'] = id_T[output_T]
    output_df['REF'] = ref_T[output_T]
    output_df['ALT'] = alt_T[output_T]
    output_df['QUAL'] = '.'
    output_df['FILTER'] = '.'
    output_df['INFO'] = '.'
    output_df['FORMAT'] = 'GT:GQ:SQ:AD'
    output_df.loc[gt_t[output_t] == 0, 'GT'] = '0|0'
    output_df.loc[gt_t[output_t] == 2, 'GT'] = '1|1'
    output_df.loc[(gt_t[output_t] == 1) & (gt_Tp[output_T, 0] == 0) & (gt_Tp[output_T, 1] == 1), 'GT'] = '0|1'
    output_df.loc[(gt_t[output_t] == 1) & (gt_Tp[output_T, 0] == 1) & (gt_Tp[output_T, 1] == 0), 'GT'] = '1|0'
    output_df.loc[(gt_t[output_t] == 1) & (gt_Tp[output_T, 0] == gt_Tp[output_T, 1]), 'GT'] = '0|1' # arbitrarily choose phase if old GT is not het
    output_df.loc[(gt_t[output_t] != 2) & (gt_Tp[output_T, 0] != -1) & (gt_Tp[output_T, 1] == -1), 'GT'] = '0' # if old GT is haploid: new 0/0 and 0/1 -> 0
    output_df.loc[(gt_t[output_t] == 2) & (gt_Tp[output_T, 0] != -1) & (gt_Tp[output_T, 1] == -1), 'GT'] = '1' # if old GT is haploid: new 1/1 -> 1
    output_df['GQ'] = gq_t[output_t].round().astype(int)
    output_df['SQ'] = sq_t[output_t].round().astype(int)
    output_df['OQ'] = oq_t[output_t].round().astype(int)
    output_df['AD'] = [f'{r},{a}' for r, a in zip(ref_count_t[output_t], alt_count_t[output_t])]
    output_df[sample] = output_df['GT'].str.cat(output_df[['GQ', 'SQ', 'AD']].astype(str), sep=':')
    output_df.to_csv(f'{args.output_prefix}.annot.tsv', sep='\t', index=False, header=False,
                     columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample])

    # output TSV highlighting changes for debugging
    delta_df = pd.DataFrame()
    delta_df['CHROM'] = chrom_T[output_T]
    delta_df['POS'] = pos_T[output_T]
    delta_df['ID'] = id_T[output_T]
    delta_df['REF'] = ref_T[output_T]
    delta_df['ALT'] = alt_T[output_T]
    delta_df['GT_OLD_0'] = gt_Tp[output_T, 0]
    delta_df['GT_OLD_1'] = gt_Tp[output_T, 1]
    delta_df['GT_OLD'] = pd.Series(gt_Tp[output_T, 0].astype(str)).str.cat(gt_Tp[output_T, 1].astype(str), sep='|')
    delta_df.loc[(gt_Tp[output_T, 0] == -1) & (gt_Tp[output_T, 1] == -1), 'GT_OLD'] = './.'
    delta_df.loc[(gt_Tp[output_T, 0] != -1) & (gt_Tp[output_T, 1] == -1), 'GT_OLD'] = \
        gt_Tp[output_T][(gt_Tp[output_T, 0] != -1) & (gt_Tp[output_T, 1] == -1), 0].astype(str) # haploid
    delta_df['GT'] = output_df['GT']
    delta_df['GT_DELTA'] = delta_df['GT'] != delta_df['GT_OLD']
    delta_df['GQ_OLD'] = gq_T[output_T]
    delta_df['GQ'] = output_df['GQ']
    delta_df['SQ_OLD'] = sq_T[output_T]
    delta_df['SQ'] = output_df['SQ']
    delta_df['AD'] = output_df['AD']
    delta_df.to_csv(f'{args.output_prefix}.delta.tsv', sep='\t', index=False, header=False,
                    columns=['CHROM', 'POS', 'ID', 'REF', 'ALT',
                             'GT_OLD', 'GT', 'GT_DELTA', 'GQ_OLD', 'GQ', 'SQ_OLD', 'SQ', 'AD'])

    # print GT confusion matrix
    logger.info('GT confusion matrix:')
    logger.info('\n' + str(pd.crosstab(delta_df['GT_OLD'], delta_df['GT'], rownames=['GT_OLD'], colnames=['GT'])))

    # create all plots
    if args.print_graphs == 1:
        logger.info('Creating plots...')

        # plot loss
        plt.plot(losses)
        plt.xlabel('iteration')
        plt.ylabel('loss')
        plt.yscale('log')
        plt.savefig(f'{args.output_prefix}.loss.png', bbox_inches='tight')
        plt.close()

        # plot gradient norms
        for name, grad_norms in gradient_norms.items():
            plt.plot(grad_norms, label=name)
        plt.xlabel('iteration')
        plt.ylabel('gradient norm')
        plt.yscale('log')
        plt.legend(loc='best')
        plt.savefig(f'{args.output_prefix}.gradients.png', bbox_inches='tight')
        plt.close()

        # plot 2D histogram of counts
        count_bins = range(args.n_count_bins)
        plt.hist2d(ref_count_t, alt_count_t, bins=count_bins, norm=matplotlib.colors.LogNorm())
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('number', rotation=-90)
        plt.xlabel('ref count')
        plt.ylabel('alt count')
        plt.title('observed')
        plt.savefig(f'{args.output_prefix}.observed_counts_2d.png', bbox_inches='tight')
        plt.close()

        # simulate counts from posterior predictive and plot 2D histogram
        alt_count_t_sim = pyro.infer.Predictive(model, guide=global_guide, num_samples=1)(hyperparameters, total_count_t)['alt_count_t'][0]
        plt.hist2d(total_count_t - alt_count_t_sim, alt_count_t_sim, bins=count_bins, norm=matplotlib.colors.LogNorm())
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('number', rotation=-90)
        plt.xlabel('ref count')
        plt.ylabel('alt count')
        plt.title('simulated from posterior predictive')
        plt.savefig(f'{args.output_prefix}.simulated_counts_2d.png', bbox_inches='tight')
        plt.close()

        # plot 1D histograms of ref counts
        plt.hist(ref_count_t, bins=count_bins, histtype='step', log=True, label='observed')
        plt.hist(total_count_t - alt_count_t_sim, bins=count_bins, histtype='step', log=True, label='simulated from posterior predictive')
        plt.xlabel('ref count')
        plt.ylabel('number')
        plt.legend(loc='best')
        plt.savefig(f'{args.output_prefix}.counts_1d_ref.png', bbox_inches='tight')
        plt.close()

        # plot 1D histograms of alt counts
        plt.hist(alt_count_t, bins=count_bins, histtype='step', log=True, label='observed')
        plt.hist(alt_count_t_sim, bins=count_bins, histtype='step', log=True, label='simulated from posterior predictive')
        plt.xlabel('alt count')
        plt.ylabel('number')
        plt.legend(loc='best')
        plt.savefig(f'{args.output_prefix}.counts_1d_alt.png', bbox_inches='tight')
        plt.close()

        # plot new GQ
        plt.tricontourf(ref_count_t, alt_count_t, gq_t, levels=range(101))
        plt.xlim([0, args.n_count_bins])
        plt.ylim([0, args.n_count_bins])
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('GQ', rotation=-90)
        cbar.set_ticks(np.arange(0, 101, 20))
        plt.xlabel('ref count')
        plt.ylabel('alt count')
        plt.savefig(f'{args.output_prefix}.GQ.png', bbox_inches='tight')
        plt.close()

        # plot new SQ
        plt.tricontourf(ref_count_t, alt_count_t, sq_t, levels=range(101))
        plt.xlim([0, args.n_count_bins])
        plt.ylim([0, args.n_count_bins])
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('SQ', rotation=-90)
        cbar.set_ticks(np.arange(0, 101, 20))
        plt.xlabel('ref count')
        plt.ylabel('alt count')
        plt.savefig(f'{args.output_prefix}.SQ.png', bbox_inches='tight')
        plt.close()

        # scatterplot of old GQ vs. new GQ
        gq_max = np.max([delta_df['GQ_OLD'], delta_df['GQ']])
        plt.plot([0, gq_max], [0, gq_max], c='grey')
        plt.scatter(delta_df['GQ_OLD'], delta_df['GQ'], s=1)
        plt.xlabel('GQ_OLD')
        plt.ylabel('GQ')
        plt.savefig(f'{args.output_prefix}.GQcorr.png', bbox_inches='tight')
        plt.close()

        # scatterplot of old SQ vs. new SQ
        sq_max = np.max([delta_df['SQ_OLD'], delta_df['SQ']])
        plt.plot([0, sq_max], [0, sq_max], c='grey')
        plt.scatter(delta_df['SQ_OLD'], delta_df['SQ'], s=1)
        plt.xlabel('SQ_OLD')
        plt.ylabel('SQ')
        plt.savefig(f'{args.output_prefix}.SQcorr.png', bbox_inches='tight')
        plt.close()

    logger.info('Done.')

if __name__ == '__main__':
    main()
