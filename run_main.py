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
import pandas as pd
import pyfasta
import seaborn as sns
from collections import defaultdict
# functions
import astat
from allel_class import Chr
import apca as apca
import aplot as aplot
import autil as autil

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
    # callset = h5py.File(h5out, mode='r')
    return(None)


def loadgenome_extradata_fx(fasta_handle, gff3_handle, meta):
    """
    """
    genome = pyfasta.Fasta(fasta_handle, key_fn=lambda key: key.split()[0])
    gff3 = allel.FeatureTable.from_gff3(gff3_handle)
    meta = pd.read_csv(meta, delimiter=",")
    return(genome, gff3, meta)


if __name__ == '__main__':
    # load data args
    fasta_handle = args.fasta
    gff3_handle = args.gff3
    meta = args.meta
    vcfin = args.vcf
    samp = args.samples
    alleles = args.altnum
    # start funcs
    makeh5fromvcf(vcfin, alleles, args.h5)
    genome, gff3, meta = loadgenome_extradata_fx(fasta_handle, gff3_handle,
                                                 meta)

#    meta = "kirfol.meta.txt.csv"
#    meta = pd.read_csv(meta, delimiter=",")
#    meta.Population = meta.ChromForm
    # load chr class, calls, pop, chrm
    var = Chr('All', 'KirFol.2L.flt.h5')

    # define sub pops
    popdict = {}
    if var.pop is not "All":
        meta = meta.ix[meta.Population.isin(var.pop)]
    for pop in meta.Population.unique():
        popdict[pop] = meta[meta.Population == pop].index.tolist()
    # chromosome lists
    chrlist = np.unique(var.chrm[:])

    # plot colors
    pop2color = {}
    palette = sns.color_palette(n_colors=len(popdict.keys()))
    for pop in popdict.keys():
        pop2color[pop] = palette.pop()
#    popdict['all'] = list(range(len(meta.Population)))

    # number of samples and list of chromosomes
    n_samples = len(var.calls['samples'])
    print("number of samples: {}".format(n_samples))
    print("list of loaded chromosomes: {}".format(chrlist))

    # Sample GT stats
    gtd = allel.GenotypeDaskArray(var.calls['calldata/GT'])
    pc_missing = gtd.count_missing(axis=0)[:].compute()  # per sample
    miss = gtd.count_missing(axis=1)[:].compute()
    pc_het = gtd.count_het(axis=0)[:].compute()  # per sample
    n_variants = var.chrm[:].shape[0]
    dep = var.calls['calldata/DP']
    dp = np.mean(dep[:, :], axis=0)
    aplot.plotstats(pc_het/n_variants, 'Heterozygous', pop2color)
    aplot.plotstats(pc_missing/n_variants, 'Missing', pop2color)
    aplot.plotstats(dp, 'Depth', pop2color)

    # Chromosome GT Stats
    for c in chrlist:
        aplot.plotvars(c, var.calls)
    for c in chrlist:
        astat.misspos(c, var.calls, miss, n_samples, window_size=10000,
                      saved=False)

    # Population Stats
    chrstatdict = defaultdict(dict)
    for c in chrlist:
        var.geno(c, meta)
        print("\nStats for Chromosome {}\n".format(c))
        # allele count object
        ac_subpops = var.gt.count_alleles_subpops(popdict, max_allele=2)
        # population stats: SNP, singleton, doubleton, #het inds, #homref, alt
        for pop in popdict.keys():
            seg = ac_subpops[pop].count_segregating()
            sing = ac_subpops[pop].count_singleton()
            doub = ac_subpops[pop].count_doubleton()
            print("{} SNPs, singleton, doubleton: {} {} {}".format(pop, seg,
                                                                   sing, doub))
            gt_subpop = var.gt.take(popdict[pop], axis=1)
            het = gt_subpop.count_het()
            ref = gt_subpop.count_hom_alt()
            alt = gt_subpop.count_hom_ref()
            print("{} hets, homalt: {} {}".format(pop, het, alt))
            chrstatdict[c][pop] = (seg, sing, doub)
        priv = astat.privalleles(ac_subpops, popdict)
        print(priv)
        diff = astat.fixdiff(ac_subpops, popdict)
        print(diff)
    thinpos = {}
    pcadict = {}
    gnudict = {}
    for c in chrlist:
        var.geno(c, meta)
        # PCA and LD thin
        # LD thin has default mac of 2
        var.miss(var.gt, var.pos, 0)  # use only sites without missing data
        gn, thinp = autil.ldthin(var.gt_m, var.pos_m, "random")
        gnudict[c] = gn
        thinpos[c] = thinp
    for c in gnudict.keys():
        # PCA
        gn = gnudict[c]
        coords, model = apca.pca_fx(gn, meta, c, pop2color, False, var.pop)
        pcadict[c] = (coords, model)



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
