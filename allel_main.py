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
import argparse
import numpy as np
# functions
import stat_fxs as sf

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
parser.add_argument('--h5', action="store_true", help="h5 exists")
parser.add_arguement('-m', "--meta", required=True, help="path to meta data")
parser.add_arguement('-s', "--samples", type=int, help="number samples")
parser.add_arguement('-a', "--altnum", type=int, default=1,
                     help="number alt alleles")
args = parser.parse_args()


if __name__ == '__main__':
    # load data args
    fasta_handle = args.fasta
    gff3_handle = args.gff3
    meta = args.meta
    vcfin = args.vcf
    samp = args.samples
    alleles = args.altnum
    # start funcs
    callset = sf.makeh5fromvcf(vcfin, alleles, args.h5)
    genome, gff3, meta = sf.loadgenome_extradata_fx(fasta_handle, gff3_handle,
                                                    meta)
    chrlist = np.unique(callset['variants/CHROM'][:])
    chrom = callset['variants/CHROM']
    # some simple stats
    print("number of samples: {}".format(len(callset['samples'])))
    print("list of loaded chromosomes: {}".format(chrlist))
    # variants over chrom for each chromosome
    for c in chrlist:
        sf.plotvars(c, callset)

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
    sf.plotstats(pc_het/n_variants, 'Heterozygous')
    sf.plotstats(pc_missing/n_variants, 'Missing')
    sf.plotstats(dp, 'Depth')
    for c in chrlist:
        sf.misspos(c, callset, pc_missing, window_size=10000, saved=True)

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
