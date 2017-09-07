#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
chrom = callset['variants/CHROM']
# get specific chromosomes
chrom[np.where(chrom[:] == "chr")]
chrom[:]  # loads into numpy array from hf5 object
pos = callset['variants/POS']
# fast with Dask
gt = allel.GenotypeDaskArray(callset['calldata/GT'])
ac = gt.count_alleles(max_allele=2).compute
# gt = allel.GenotypeArray(callset['calldata/GT'])
# gt[variants, samples], gt[1:3, :] 2-4th variant all samples
# gt[:, 0:2], all variants, 1st and 2nd samples

@author: scott
"""

import allel
import h5py
import pandas as pd
import pyfasta
import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import bisect

from collections import defaultdict
# made functions
from allel_class import Chr
import run_admixture as admx
import run_diversity as div
from vcf2plink import vcf2plink_fx

parser = argparse.ArgumentParser()
parser.add_argument('-fa', '--fasta', help='path to fasta genome')
parser.add_argument('-gff3', '--ggf3feature', help="path to gff3")
parser.add_arguement('-v', "--vcf", help="path to vcf")
parser.add_argument('--h5file', action="store_true", help="h5 exists")
parser.add_arguement('-m', "--meta", required=True, help="path to meta data")
parser.add_arguement('-s', "--samples", type=int, help="number samples")
parser.add_arguement('-a', "--altnum", type=int, default=1,
                     help="number alt alleles")
args = parser.parse_args()


def makeh5fromvcf(vcfin, altnum, hf5):
    """
    """
    if hf5:
        pass
    else:
        h5out = "{}.h5".format(vcfin.split(".")[:-2])
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
    meta = pd.read_csv(args.meta, delimiter=" ")
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


if __name__ == '__main__':
    # load data args
    fasta_handle = args.fasta
    gff3_handle = args.gff3
    meta = pd.read_csv(args.meta, delimiter=" ")
    vcfin = args.vcf
    samp = args.samples
    alleles = args.altnum
    chrom = args.chromlist
    # start funcs
    callset = makeh5fromvcf(vcfin, alleles)
    genome, gff3, meta, chrlist = loadgenome_extradata_fx(fasta_handle,
                                                          gff3_handle, meta,
                                                          chrom)
    # some simple stats
    print("number of samples: {}".format(len(callset['samples'])))
    chrlist = np.unique(chrom[:])
    print("list of loaded chromosomes: {}".format(chrlist))
    # variants over chrom for each chromosome
    for c in chrlist:
        plotvars(c, callset)
    # genotype object as array
    gtd = allel.GenotypeDaskArray(callset['calldata/GT'])  # must add compute()
    gt = allel.GenotypeArray(callset['calldata/GT'])
    # stats
    het = gtd.is_het().compute()
    ac = gtd.count_alleles().compute()
    c_het = gtd.count_het(axis=1).compute()  # per site
    c_miss = gtd.count_missing(axis=1).compute()  # per site
    # singletons = np.where(hc == 1)
    pc_missing = gtd.count_missing(axis=0)[:].compute()  # per sample
    pc_het = gtd.count_het(axis=0)[:].compute()  # per sample
    n_variants = chrom[:].shape[0]
    dep = callset['calldata/DP']
    dp = np.mean(dep[:, :], axis=0)
    # function plots
    plotstats(pc_het/n_variants, 'Heterozygous')
    plotstats(pc_missing/n_variants, 'Missing')
    plotstats(dp, 'Depth')
    for c in chrlist:
        misspos(c, callset, pc_missing, window_size=10000, saved=True)

    # PCA
    thinmiss = {}
    chrpos = {}
    c = []
    p = []
    for nchr in chrlist:
        Chr1 = chrdict[nchr]
        Chr1.missing(Chr1.genotypes, Chr1.positions, 0)
        Chr1.rmvsingletons(Chr1.gtnomiss, Chr1.posnomiss)
        Chr1.ldfilter(Chr1.gtNOsingletons, Chr1.posNOsingletons, 1000, 500, .2, 1)
        thinmiss[nchr] = Chr1.ldpos
        c.append(Chr1.genotypes)
        p.append(Chr1.corr)
        chrpos[nchr] = Chr1.positions
        admx.pca_fx(Chr1.corr, meta, False, nchr)
    # concat
    chrcat_pca = chrcat_fx(p)
    chrcat_gt = chrcat_fx(c)
    # all chr PCA
    coords = admx.pca_fx(chrcat_fx(p), meta, False, "All")
    # ADMIXTURE
    vcf2plink_fx(thinmiss, vcf)
    admx.admixture_prog_fx('thinnedplink.bed')
    # FST, Fstatistics, doubletons, jsfs
    poplist = ["PNG", "Haiti", "Mali", "Kenya"]
    pairlist = ['PNG_Haiti', 'PNG_Mali', 'PNG_Kenya', 'Haiti_Kenya',
                'Haiti_Mali', 'Kenya_Mali', 'Haiti_Africa', 'PNG_Africa']
    png = meta[meta.population == "PNG"].index.tolist()
    mali = meta[meta.population == "Mali"].index.tolist()
    haiti = meta[meta.population == "Haiti"].index.tolist()
    kenya = meta[meta.population == "Kenya"].index.tolist()
    subpops = {
            'PNG': png,
            'Haiti': haiti,
            'Kenya': kenya,
            'Mali': mali,
            'Africa': mali + kenya,
            'PNG_Haiti': png + haiti,
            'PNG_Mali': png + mali,
            'PNG_Kenya': png + kenya,
            'Haiti_Kenya': haiti + kenya,
            'Haiti_Mali': haiti + mali,
            'Kenya_Mali': kenya + mali,
            'Haiti_Africa': haiti + mali + kenya,
            'PNG_Africa': png + mali + kenya
            }
    # allele counts for all subpops
    ac_subpopsdict = {}
    for nchr in chrlist:
        ac_subchr = chrdict[nchr].genotypes.count_alleles_subpops(
                subpops, max_allele=1)
        ac_subpopsdict[nchr] = ac_subchr
    # concatenate all chr for each pop
    ac_subpopscat = {}
    for pop in subpops.keys():
        ac_sub = []
        for nchr in chrlist:
            ac_sub.append(ac_subpopsdict[nchr][pop])
        ac_subpopscat[pop] = chrcat_fx(ac_sub)
    # structure stats
    d2, dshare = admx.doubletons_fx(ac_subpopscat, subpops)
    df_FST = admx.pairFST_fx(ac_subpopsdict, ac_subpopscat, chrcat_gt, subpops, poplist,
                    False, chrpos)
    fstats = admx.Fstatistics_fx(ac_subpopsdict, subpops)
    # diversity stats
    div.sfs_fx(ac_subpopscat)

    # load feature data
#    gff3.shape
#    np.unique(gff3.type)
#    is_CDS = gff3[(gff3.type == b'CDS')]
#    is_CDS.shape

#        x = []
#    for key in posDict:
#        x.append(list(posDict[key]))
#    allel.SortedMultiIndex(chromlist,[list(posDict)])
#    # stats
#    vtblCounter = 0
#    for name in vtblDict.keys():
#        print("{}:{}".format(name, len(vtblDict[name])))
#        vtblCounter += len(vtblDict[name])
#    print("total:{}".format(vtblCounter))
#    # stats
#    genosCounter = 0
#    for name in genosDict.keys():
#        print("{}:{}".format(name, len(genosDict[name])))
#        vtblCounter += len(genosDict[name])
#    print("total:{}".format(genosCounter))
    # start admixture
