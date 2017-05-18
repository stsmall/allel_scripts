#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 12:26:31 2017

@author: scott
"""

import allel
import numpy as np
import pyfasta
import argparse
import pandas as pd
import h5py
from collections import defaultdict
# made functions
from allel_class import Chr
from admixture import pca_fx
from vcf2plink import vcf2plink_fx
from admixture import admixture_prog_fx
from admixture import doubletons_fx
from admixture import pairFST_fx
from admixture import Fstatistics_fx
from diversity import sfs_fx

parser = argparse.ArgumentParser()
parser.add_argument('-fa', '--fasta',
                    help='path to fasta genome', required=True)
parser.add_argument('-gff3', '--ggf3feature', help="path to gff3")
parser.add_arguement('-v', "--vcf", help="path to vcf")
parser.add_argument('-h5', '--h5file', help="path to h5")
parser.add_arguement('-m', "--meta", help="path to meta data")
parser.add_arguement('-s', "--samples", type=int, help="number samples")
parser.add_arguement('-a', "--alleles", type=int, default=2,
                     help="number alleles")
args = parser.parse_args()


def loadgenome_extradata_fx(fasta_handle, gff3_handle, meta):
    """
    """
    genome = pyfasta.Fasta(fasta_handle, key_fn=lambda key: key.split()[0])
    gff3 = allel.FeatureTable.from_gff3(gff3_handle)
    return(genome, gff3)


def loadgenome_fx(callset_handle, chrlist):
    """
    """
    callset = h5py.File(callset_handle, mode='r')
    chrdict = {}
    for nchr in chrlist:
        Chr1 = Chr(nchr, callset)
        Chr1.loadpos(False)
        chrdict[nchr] = Chr1
    return(chrdict)


def chrcat_fx(l):
    """
    l: list of genotype or other allel objects by chromosome
    """
    chrcat = l[0]
    i = 1
    while i < len(l):
        chrcat = chrcat.concatenate(l[i])
        i += 1
    return(chrcat)

if __name__ == '__main__':
    fasta_handle = args.fasta
    gff3_handle = args.gff3
    meta = pd.read_csv(args.meta, delimiter=" ")
    vcf = args.vcf
    samp = args.samples
    alleles = args.alleles
    genome, gff3 = loadgenome_extradata_fx(fasta_handle, gff3_handle, meta)
    callset_handle = args.h5file

    chrlist = ["Wb_Chr1_0", "Wb_Chr1_1", "Wb_Chr2_0", "Wb_Chr2_1", "Wb_Chr2_2",
               "Wb_Chr2_3", "Wb_Chr3_0", "Wb_Chr3_1", "Wb_Chr4_0", "Wb_Chr4_1",
               "Wb_Chr4_2", "Wb_ChrX_0", "Wb_ChrX_1", "Wb_ChrX_2", "Wb_ChrX_3"]
    chrdict = loadgenome_fx(callset_handle, chrlist)
    # PCA
    thinmiss = {}
    c = []
    p = []
    for nchr in chrlist:
        Chr1 = chrdict[nchr]
        Chr1.missing(Chr1.genotypes, Chr1.positions, 0)
        Chr1.rmvsingletons(Chr1.gtnomiss, Chr1.posnomiss)
        Chr1.ldfilter(Chr1.gtNOsingletons, Chr1.posNOsingletons, 50, 5, .1, 1)
        thinmiss[nchr] = Chr1.ldpos
        c.append(Chr1.genotypes)
        p.append(Chr1.corr)
        pca_fx(Chr1.corr, meta, False, nchr)
    # concat
    chrcat_pca = chrcat_fx(p)
    chrcat_gt = chrcat_fx(c)
    # all chr PCA
    coords = pca_fx(chrcat_fx(p), meta, False, "All")
    # ADMIXTURE
    vcf2plink_fx(thinmiss, vcf)
    admixture_prog_fx('thinnedplink.bed')
    # FST, Fstatistics, doubletons, jsfs
    poplist = ["PNG", "Haiti", "Mali", "Kenya"]
    pairlist = ['PNG_Haiti', 'PNG_Mali', 'PNG_Kenya', 'Haiti_Kenya',
                'Haiti_Mali', 'Kenya_Mali', 'Haiti_Africa', 'Africa_PNG']
    png = meta[meta.population == "PNG"].index.tolist()
    mali = meta[meta.population == "Mali"].index.tolist()
    haiti = meta[meta.population == "Haiti"].index.tolist()
    kenya = meta[meta.population == "Kenya"].index.tolist()
    subpops = {
            'PNG': png,
            'Haiti': haiti,
            'Kenya': kenya,
            'Mali': mali,
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
    d2, dshare = doubletons_fx(ac_subpopscat, subpops)
    pairFST_fx(ac_subpopsdict, ac_subpopscat, pairlist, False)
    fstats = Fstatistics_fx(ac_subpopsdict)
    # diversity stats
    sfs_fx(ac_subpopscat)

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
