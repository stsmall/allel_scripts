#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 16:30:34 2017
Assorted functions for using scikit allel
@author: scott
"""
import allel
import h5py
import pandas as pd
import pyfasta
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import bisect
import hmmlearn
from accessmasknp import access_mask_fx


def jackknife(x):
    """
    """

    vals = np.empty(x.shape, dtype=float)
    x = np.ma.asarray(x)
    x.mask = np.zeros(x.shape, dtype=bool)
    for i in range(x.size):
        x.mask[i] = True
        vals[i] = np.mean(x)
        x.mask[i] = False
    n = x.size
    try:
        sv = ((n - 1) / n) * np.sum((vals - vals.mean()) ** 2)
    except ZeroDivisionError:
        se = 0.0000000001
    se = np.sqrt(sv)
    m = np.mean(vals)
    z = m / se
    # here the forumula is actually m - 0 / se, 0 is expected value
    return(m, se, z, vals)


def roh(subpops, chrpos, chrdict, nchr, chrmaskdict):
    """Runs of Homozygosity for each individual sample
    phet_roh: prob of obs a het within a RoH ~mutation rate
    phet_nonroh: >1 prob of obs het outside of RoH, ~nucleotide div
    """
    m = access_mask_fx(chrmaskdict, nchr)
    for pop in subpops:  # Haiti
        for indx in pop:  # 1, 2, 3, 4, ... 8
            gn = chrdict[nchr].genotypes[:, indx]
            het_mask = gn.is_het()
            gt = gn.compress(het_mask)
            call_mask = gn.is_called()  # all sites called
            posmask = call_mask * het_mask
            pos = chrdict[nchr].positions[posmask]
            df, prop = allel.stats.roh.roh_mhmm(gt, pos, phet_roh=0.001,
                                                phet_nonroh=(0.0025, 0.01),
                                                is_accessible=m)
    return(None)


def makeh5fromvcf(vcfin, altnum, hf5):
    """
    """
    h5out = "{}.h5".format(vcfin.split(".")[:-2])
    if hf5:
        pass
    else:
        fieldsfromvcf = ['samples', 'calldata/GQ', 'variants/ALT',
                         'variants/REF', 'variants/QUAL', 'variants/CHROM',
                         'variants/POS', 'variants/AF', 'variants/AB',
                         'variants/MQM', 'variants/DP', 'calldata/DP',
                         'calldata/AD', 'calldata/GT']
        allel.vcf_to_hdf5(vcfin, h5out, fields=fieldsfromvcf,
                          types={'calldata/GQ': 'float32'}, alt_number=2)
    callset = h5py.File(h5out, mode='r')
    return(callset)


def loadgenome_extradata_fx(fasta_handle, gff3_handle, meta):
    """
    """
    genome = pyfasta.Fasta(fasta_handle, key_fn=lambda key: key.split()[0])
    gff3 = allel.FeatureTable.from_gff3(gff3_handle)
    meta = pd.read_csv(meta, delimiter=" ")
    return(genome, gff3, meta)


def plotvars(chrm, callset, window_size=10000, title=None, saved=False):
    """
    """
#    chrm = chrm.decode("utf-8")
    chrom = callset['variants/CHROM']
    chrom_mask = np.where(chrom[:] == chrm)
    pos = callset['variants/POS']
    p = pos[:][chrom_mask]
    varpos = allel.SortedIndex(p)
    # setup windows
    bins = np.arange(0, varpos.max(), window_size)
    # use window midpoints as x coordinate
    x = (bins[1:] + bins[:-1])/2
    # compute variant density in each window
    h, _ = np.histogram(varpos, bins=bins)
    y = h / window_size
    # plot
    fig, ax = plt.subplots(figsize=(12, 3))
    sns.despine(ax=ax, offset=10)
    ax.plot(x, y)
    ax.set_xlabel('Chromosome position (bp)')
    ax.set_ylabel('Variant density (bp$^{-1}$)')
    if title:
        ax.set_title(title)
    else:
        ax.set_title(chrm.decode("utf-8"))
    if saved:
        fig.savefig("{}.vars.pdf".format(chrm.decode("utf-8")),
                    bbox_inches='tight')


def misspos(chrm, callset, pc, window_size=10000, title=None, saved=False):
    """
    """
#    chrm = chrm.decode("utf-8")
    chrom = callset['variants/CHROM']
    chrom_mask = np.where(chrom[:] == chrm)
    pos = callset['variants/POS']
    p = pos[:][chrom_mask]
    varpos = allel.SortedIndex(p)
    bins = np.arange(0, varpos.max(), window_size)
    # use window midpoints as x coordinate
    x = bins
    miss_site = pc[:][chrom_mask]
    yy = []
    for i, j in enumerate(x):
        try:
            left = bisect.bisect_left(varpos, j)
            right = bisect.bisect_left(varpos, x[i+1]) - 1
            yy.append(np.mean(miss_site[left:right]))
        except:
            yy.append(0)
    y = np.array(yy)
    # plot
    fig, ax = plt.subplots(figsize=(12, 3))
    sns.despine(ax=ax, offset=10)
    ax.plot(x, y)
    ax.set_xlabel('Chromosome position (bp)')
    ax.set_ylabel('Missing Count')
    if title:
        ax.set_title(title)
    else:
        ax.set_title(chrm.decode("utf-8"))
    if saved:
        fig.savefig("{}.miss.pdf".format(chrm.decode("utf-8")),
                    bbox_inches='tight')


def plotstats(pc, title, saved=False):
    """
    """
    fig, ax = plt.subplots(figsize=(12, 4))
    sns.despine(ax=ax, offset=10)
    left = np.arange(len(pc))
    ax.bar(left, pc)
    ax.set_xlim(0, len(pc))
    ax.set_xlabel('Sample index')
    ax.set_ylabel('Percent calls')
    ax.set_title(title)
    if saved:
        fig.savefig("{}.pdf".format(title), bbox_inches='tight')