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
# functions
import astat
from allel_class import Chr
import apca as apca
import autil as autil
import adiff as ad
#import adiv as av

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

    meta = "kirfol.meta.txt.csv"
    meta = pd.read_csv(meta, delimiter=",")
    var = Chr('All', 'KirFol.2L.flt.h5')
    popdict = autil.subpops(var, meta)
    chrlist = np.unique(var.chrm[:])
    pop2color = autil.popcols(popdict)

    # number of samples and list of chromosomes
    n_samples = len(var.calls['samples\n'])
    print("number of samples: {}".format(n_samples))
    print("list of loaded chromosomes: {}".format(chrlist))

    # Sample GT stats
    miss = astat.gtstats(var.calls, pop2color, var.chrm.shape[0])
    # Chr stats
    astat.chrstats(chrlist, var.calls, miss, n_samples)
    # Population Stats
    diff, chrstatdict = astat.popstats(chrlist, meta, popdict, var)

    # PCA and LD thin
    thinpos = {}
    pcadict = {}
    gnudict = {}
    for c in chrlist:
        # LD thin
        var.geno(c, meta)
        var.miss(var.gt, var.pos, 0)  # use only sites without missing data
        gn, thinp = autil.ldthin(var.gt_m, var.pos_m, "thin", iters=5)
        gnudict[c] = gn
        thinpos[c] = thinp
    for c in gnudict.keys():
        # PCA
        gn = gnudict[c]
        coords, model = apca.pca_fx(gn, meta, c, pop2color, False, var.pop,
                                    bykary=True)
        pcadict[c] = (coords, model)

    # ADMIXTURE input files
    autil.vcf2plink(thinpos)  # make input files

    # FST, Fstatistics, doubletons, jsfs
    fstdict = {}
    dxydict = {}
    fstatdict = {}
    doubdict = {}
    for c in chrlist:
        var.geno(c, meta)
        print("\nStats for Chromosome {}\n".format(c))
        # allele count object
        ac_subpops = var.gt.count_alleles_subpops(popdict, max_allele=2)
        df_FST = ad.pairFST(c, ac_subpops, var, popdict)
        fstdict[c] = df_FST
        df_dxy = ad.Dxy(c, ac_subpops, var, popdict)
        dxydict[c] = df_dxy

        fstats = ad.Fstats(c, ac_subpops)
        d2, dshare = ad.doubletons(c, ac_subpops)

#    # Diversity statistics
#    for c in chrlist:
#        var.geno(c, meta)
#        print("\nStats for Chromosome {}\n".format(c))
#        # allele count object
#        ac_subpops = var.gt.count_alleles_subpops(popdict, max_allele=2)
#        av.sfs_fx(c, ac_subpops)

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
