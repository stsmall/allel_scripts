#!/usr/bin/env python3
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
gt.compress(mask)
gt.take(actual_indexes)
gt.subset()

@author: stsmall
"""
from __future__ import division
from __future__ import print_function
from imp import reload
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
import afst as afst
import adxy as adxy
import adiff as ad
from af234 import pF2
import asfs as asfs
import adiv as av
import ald as ald

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

    meta = "/home/scott/Documents/Wb/Wb_sWGA/data_files/allel/WbAllpops.47.info"
    meta = pd.read_csv(meta, delimiter=",")
    var = Chr('All', 'WbAllpops.impute.allchr.47.h5')
    popdict = autil.subpops(var, meta, bypop=True, bykary=False)
    pop2color = autil.popcols(popdict)
    chrlist = np.unique(var.chrm[:])
    # chrlist = del chrlist['Wb_ChrX_0']
    chrlen = {}
    with open("chr_info", 'r') as c:
        for line in c:
            x = line.strip().split()
            chrlen[x[0]] = int(x[1])

    # number of samples and list of chromosomes
    n_samples = len(var.calls['samples'])
    print("number of samples: {}".format(n_samples))
    print("list of loaded chromosomes: {}".format(chrlist))

    # Sample GT stats
    miss = astat.gtstats(var.calls, pop2color, var.chrm.shape[0])
    # Chr stats
    astat.chrstats(chrlist, var.calls, miss, n_samples)
    # Population Stats
    chrstatdict, diff, priv = astat.popstats(chrlist, meta, popdict, var)

    # PCA and LD thin
    thinpos = {}
    pcadict = {}
    gnudict = {}
    for c in chrlist:
        # LD thin
        print(c)
        var.geno(c, meta)
        var.miss(var.gt, var.pos, 0)  # use only sites without missing data
        gn, thinp = autil.ldthin(var.gt, var.pos, "thin", iters=5)
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

    # TODO: cut partial windows

    # FST, Fstatistics, doubletons, jsfs
    fstdict = {}
    dxydict = {}
    f2dict = {}
    doubdict = {}
    jsfsdict = {}
    for c in chrlist:
        var.geno(c, meta)
        # var.miss(var.gt, var.pos, .20)
        # var.mac(var.gt, var.pos, 1)
        print("\nStats for Chromosome {}\n".format(c))
        # allele count object
        ac_subpops = var.gt.count_alleles_subpops(popdict, max_allele=2)
        # TODO: error with Wb
        # df_FST = afst.pairFST(c, chrlen[c], ac_subpops, var, popdict,
        #                      plot=True)
        #fstdict[c] = df_FST
        df_dxy = adxy.pairDxy(c, chrlen[c], ac_subpops, var.pos, plot=True)
        dxydict[c] = df_dxy
#        jsfsdict = ad.jsfs(ac_subpops, save=False)
#        dshare = ad.shared_doubletons(ac_subpops)
#        doubdict[c] = dshare
# TODO: fix f2
#        f2 = pF2(ac_subpops)
#        f2dict[c] = f2

    # Diversity statistics
    sfsdict = {}
    pidict = {}
    tajddict = {}
    thetadict = {}
    for c in chrlist:
        var.geno(c, meta)
        print("\nStats for Chromosome {}\n".format(c))
        # var.miss(var.gt, var.pos, .20)
        # var.mac(var.gt, var.pos, 1)
        # allele count object
        ac_subpops = var.gt.count_alleles_subpops(popdict, max_allele=1)
        #sfs = asfs.sfs_plot(c, ac_subpops)
        #sfsdict[c] = sfs
        pi = av.pi(c, chrlen[c], ac_subpops, var.pos, plot=True)
        pidict[c] = pi
        d = av.tajd(c, chrlen[c], ac_subpops, var.pos, plot=True)
        tajddict[c] = d
        t = av.theta(c, chrlen[c], ac_subpops, var.pos, plot=True)
        thetadict[c] = t
    # diversity box plot: sumdict() autil then boxplot in aplot
    # LD decay plot
    lddict = {}
    for c in chrlist:
        var.geno(c, meta)
        print("\nStats for Chromosome {}\n".format(c))
        var.miss(var.gt, var.pos, .20)
        # allele count object
        ac_subpops = var.gt.count_alleles_subpops(popdict, max_allele=1)
        lddict[c] = ald.ld_decay(c, chrlen[c], ac_subpops, popdict,
                                 pop2color, var)




x = []
for p in tajddict.keys():
    x.append((tajddict[p]["Haiti"][2][0]))
m = np.concatenate(x).ravel()
n = m[~np.isnan(m)]
b,bins,patches = plt.hist(n, 50, density=True)
