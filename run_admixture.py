#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May  5 15:38:22 2017
Pairwise Fst, F2, F3, F3, ABBA-BABA, shared doubletons, admixture plots
requires PLINK, ADMIXTURE, vcftools, allel
@author: scott
"""

import allel
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import numpy as np
sns.set_style('white')
sns.set_style('ticks')


def plot_pca_coords(coords, model, pc1, pc2, ax, samples):
    """
    """
    sns.despine(ax=ax, offset=5)
    x = coords[:, pc1]
    y = coords[:, pc2]
    pop_colours = {
        'Haiti': '#FF0000',
        'Kenya': '#008000',
        'Mali': '#00FFFF',
        'PNG': '#90EE90'}
    sample_population = samples.population.values
    populations = samples.population.unique()
    for pop in populations:
        flt = (sample_population == pop)
        ax.plot(x[flt], y[flt], marker='o', linestyle=' ',
                color=pop_colours[pop],
                label=pop, markersize=6, mec='k', mew=.5)
    ax.set_xlabel('PC%s (%.1f%%)' % (pc1+1,
                  model.explained_variance_ratio_[pc1]*100))
    ax.set_ylabel('PC%s (%.1f%%)' % (pc2+1,
                  model.explained_variance_ratio_[pc2]*100))


def fig_pca(coords, model, title, sample_population, pc1, pc2):
    """
    """
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 1, 1)
    plot_pca_coords(coords, model, pc1, pc2, ax, sample_population)
    ax.legend(bbox_to_anchor=(1, 1), loc='upper left')
    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    fig.savefig("PCA_ax{}_ax{}.pdf".format(pc1, pc2), bbox_inches='tight')


def pca_fx(gtthin, df_samples, pcall, nchr):
    """
    """
    gnu = gtthin[:]  # uncompress
    # PCA
    coords, model = allel.pca(gnu, n_components=10, scaler='patterson')
    title = "PCA Chr:{}, var:{}".format(nchr, gnu.shape[0])
    if pcall:
        i = 0
        while i < 9:
            fig_pca(coords, model, title,
                    df_samples, i, i+1)
            i += 2
    else:
        fig_pca(coords, model, title, df_samples,
                0, 1)
    return(coords)


def admixture_prog_fx(plink):
    """
    """
    for k in range(5):
        command = "admixture --cv " + plink + " " + k + " | tee log" + k\
                  + ".out"
        proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        proc.wait()

    # PLOT

    return(None)


def jSFS_fx(jsfs, pop1, pop2):
    """
    """
    fig, ax = plt.subplots(figsize=(6, 6))
    allel.stats.plot_joint_sfs(jsfs, ax=ax)
    ax.set_ylabel('Alternate allele count, {}'.format(pop1))
    ax.set_xlabel('Alternate allele count, {}'.format(pop2))
    return(None)


def pairFST_fx(ac_subpopsdict, ac_subpops, genotypes_seg, pairlist, jsfs_r):
    """
    """
    for pair in pairlist:
        pop1 = pair.split("_")[0]
        pop2 = pair.split("_")[1]
        ac_seg = ac_subpops.compress(genotypes_seg[pair])
        # jsfs
        if jsfs_r:
            jsfs = allel.stats.joint_sfs(ac_seg[pop1][:, 1],
                                         ac_seg[pop2][:, 1])
            jSFS_fx(jsfs, pop1, pop2)
        # pairwise FST, std error
        fst, fst_se, _, _ = allel.stats.blockwise_hudson_fst(ac_seg[pop1],
                                                             ac_seg[pop2],
                                                             blen=2000)
        print("%s, Hudson's Fst: %.3f +/- %.3f" % (pair, fst, fst_se))
        # pairwise FST per chrom
        for nchr in ac_subpopsdict.keys():
            is_seg = ac_subpopsdict[nchr][pair].is_segregating()[:]
            ac_seg = ac_subpopsdict[nchr].compress(is_seg)
            fst, fst_se, _, _ = allel.stats.blockwise_hudson_fst(ac_seg[pop1],
                                                                 ac_seg[pop2],
                                                                 blen=2000)
            print("%s, %s, Hudson's Fst: %.3f +/- %.3f" % (nchr, pair, fst,
                                                           fst_se))
    return(None)


def Fstatistics_fx(chrdict, meta, subpops):
    """
    """
    return(None)


def doubletons_fx(ac_subpopscat, subpops):
    """
    double: dict(list)
        PNG: ([False, False, True ...
    """
    # count doubletons
    d2 = {}
    dshare = {}
    for pop in subpops.keys():
        d2[pop] = ac_subpopscat[pop].count_doubleton(allele=1)
    for pop in d2.keys():
        if "_" in pop:
            if "Africa" in pop:
                dshare[pop] = d2[pop.split("_")[0]] + d2["Mali"]\
                    + d2["Kenya"] - d2[pop]
            else:
                dshare[pop] = d2[pop.split("_")[0]] + d2[pop.split("_")[1]]\
                    - d2[pop]
        else:
            pass
    return(d2, dshare)
