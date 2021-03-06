#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Sep 26 18:58:35 2017 @author: stsmall

dxy: pairwise distance / bases; low values possible introgression
    dxy as number of sequence diff between any 2 sequences, x and y, in two
    taxa, X and Y (divided by the number of sites), then dxy is the average
    distance between all sequences in the two species. # above assumes no
    variation in neutral mutation rate, low mutation rate can be mistaken for
    recent introgression

dmin: min(dxy), requires haplotypes
    minimum distance among all pairing of haplotypes in the 2 species. Pvalue
    by coalescent with no migration or from other parts of the genome average.
    # above assumes no variation in neutral mutation rate, low mutation rate
    can be mistaken for recent introgression

dout: (dXO + dYO) / 2
    average distance between species X and the Outgroup and Species Y and the
    Outgroup.

RND (relative node depth): dxy / dout
     Robust to low mutation rates like HKY test if neutrality. # not sensitive
     to low-frequency migrants.
    calculate Dxy in windows on haplotypeArray between Species X and Y
    calculate Dxy between Species X and Outgroup
    calculate Dxy between Species Y and Outgroup

Gmin: dmin / dxy
    low power, useful when migration prob and relative migration time is low
    This is likely due to the fact that as migrant lineages rise in frequency,
    dXY also gets lower. As a migrant haplotype approaches fixation, the ratio
    of dmin to dXY approaches 1

RNDmin: dmin / dout
    Similarly, like both dmin and Gmin, RNDmin should be sensitive to even rare
    migrant haplotypes. In addition, we expect RNDmin to be powerful even when
    migrants are high in frequency

D2 (bananas paper)
    distance between A & C when A,B are sister minus the distance between A & C
    when B,C are sister. Likely in ETE3 to identify trees/coordinates meeting
    the sister species requirement, then use sequence to calculate distance use
    phylip DNAdist or something like that; (d_AC | A,B) - (d_AC | B,C)

F_D: (Martin 2015)

R_D: (Racimo 2016)

Utwenty: (Racimo 2016, Jagoda 2018)

QninetyFive: (Racimo 2016, Jagoda 2018)

Dp_intro: (Racimo 2016)

Dp_combo: (Racimo 2016)

r2_intro: (Racimo 2016)

r2_combo: (Racimo 2016)

## determine significance by coalescent simulation w/ and w/out migration,
    also use genome average, X, in anopheles assuming that region is not
    introgressed.

"""

import allel
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import pandas as pd
from allel_class import Chr
import autil as autil
from itertools import combinations
sns.set_style('white')
sns.set_style('ticks')

parser = argparse.ArgumentParser()
parser.add_argument('-v', "--vcfFile", help="path to vcf")
parser.add_argument('--h5', action="store_true", help="h5 exists")
parser.add_argument('-m', "--meta", required=True, help="path to meta data")
parser.add_argument('-o', "--out", type=str, required=True)
parser.add_argument('-t', "--target", type=str, required=True)
parser.add_argument('-b', "--bait", type=str, required=True)
args = parser.parse_args()


def makeh5fromvcf(vcfin, altnum, hf5):
    """
    """
    h5out = "{}.h5".format(vcfin)
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


def Utwenty(pos, gt, out, target, bait, clen, tfreq=0.20, winSize=10000, winStep=None):
    """Calculates the U20 statistic from Racimo 2016
        U20(0.01, 0.20, 1.0)
        U50(0.01, 0.50, 1.0)
        U0(0.01, 0.0, 1.0)
    """

    # remove sites:
    # 1) <0.01 in outGroup
    gtO = gt.take(out, axis=1)
    ac = gtO.count_alleles()
    acf = ac.to_frequencies()
    outMask = acf[:, 1] <= 0.01
    pos = pos[outMask]
    gt = gt.compress(outMask, axis=0)
    # 2) >= .20 in target
    gtT = gt.take(target, axis=1)
    ac = gtT.count_alleles()
    acf = ac.to_frequencies()
    targetMask = acf[:, 1] >= tfreq
    pos = pos[targetMask]
    gt = gt.compress(targetMask, axis=0)
    # 3) == 1.0 freq in bait
    gtB = gt.take(bait, axis=1)
    ac = gtB.count_alleles()
    acf = ac.to_frequencies()
    baitMask = acf[:, 1] == 1.0
    pos = pos[baitMask]
    gt = gt.compress(baitMask, axis=0)
    # count sites remaining in sliding windows
    U20 = allel.windowed_count(pos, winSize, 1, clen, winStep)
    return(U20)


def QninetyFive(pos, gt, out, target, bait, clen, ofreq=0.01, winSize=10000, winStep=None):
    """Calculates the Q95 test statistic from Racimo 2016
        Q95(0.01, 1.0)
        Q95(0.10, 1.0)
    """

    # get 95% quantile of allele count for target
    gtT = gt.take(target, axis=1)
    ac = gtT.count_alleles()
    acf = ac.to_frequencies()
    acq95 = np.quantile(acf[:, 1], 0.95)
    # remove sites:
    # 1) NOT 95% quantile counts in target
    q95mask = acf[:, 1] >= acq95
    pos = pos[q95mask]
    gt = gt.compress(q95mask, axis=0)
    # 2) NOT 1% freq in outGroup
    gtO = gt.take(out, axis=1)
    ac = gtO.count_alleles()
    acf = ac.to_frequencies()
    prctmask = acf[:, 1] <= ofreq
    pos = pos[prctmask]
    gt = gt.compress(prctmask, axis=0)
    # 3) NOT 100% freq in bait
    gtB = gt.take(bait, axis=1)
    ac = gtB.count_alleles()
    acf = ac.to_frequencies()
    prctmask = acf[:, 1] == 1.0
    pos = pos[prctmask]
    gt = gt.compress(prctmask, axis=0)
    # count sites remaining in sliding windows
    Q95 = allel.windowed_count(pos, winSize, 1, clen, winStep)
    return(Q95)


def F_d(pos, gt, out, target, bait, ofreq=0.01, winSize=10000, winStep=None):
    """Fd from Martin 2015
    """
    return(None)


def Dout():
    """
    dout: (dXO + dYO) / 2
    average distance between species X and the Outgroup and Species Y and the
    Outgroup.
    """


def Dxy():
    """
    dxy: pairwise distance / bases; low values possible introgression
    dxy as number of sequence diff between any 2 sequences, x and y, in two
    taxa, X and Y (divided by the number of sites), then dxy is the average
    distance between all sequences in the two species. # above assumes no
    variation in neutral mutation rate, low mutation rate can be mistaken for
    recent introgression
    """
    # Dxy
    dxydict = {}
    for c in chrlen.keys():
        var = Chr('All', '{}.FSG.SNP.recode.h5'.format(c))
        var.geno(c, meta)
        # var.miss(var.gt, var.pos, .20)
        # var.mac(var.gt, var.pos, 1)
        print("\nStats for Chromosome {}\n".format(c))
        # allele count object
        ac_subpops = var.gt.count_alleles_subpops(popdict, max_allele=2)
        df_dxy = adxy.pairDxy(c, chrlen[c], ac_subpops, var.pos, plot=True)
        dxydict[c] = df_dxy

def Dmin():
    """
    REQUIRES PHASES, HAPLOTYPEARRAY
    dmin: min(dxy), requires haplotypes
    minimum distance among all pairing of haplotypes in the 2 species. Pvalue
    by coalescent with no migration or from other parts of the genome average.
    # above assumes no variation in neutral mutation rate, low mutation rate
    can be mistaken for recent introgression
    """


def RND(dxydict, sp1, out):
    """
    RND (relative node depth): dxy / dout
     Robust to low mutation rates like HKY test if neutrality. # not sensitive
     to low-frequency migrants.
    calculate Dxy in windows on haplotypeArray between Species X and Y
    calculate Dxy between Species X and Outgroup
    calculate Dxy between Species Y and Outgroup
    """
    out = "Riv"
    sp1 = ["Van", "Par", "FunMoz", "Long", "Like"]
    RNDdict = {}
    chrdict = {}
    for c in chrlen.keys():
        for i, j in combinations(sp1, 2):
            try:
                dxyP = dxydict[c]["{}-{}".format(i, j)][2][0]
            except KeyError:
                dxyP = dxydict[c]["{}-{}".format(j, i)][2][0]
            dxy1out = dxydict[c]["{}-{}".format(i, out)][2][0]
            dxy2out = dxydict[c]["{}-{}".format(j, out)][2][0]
            dxyout = (dxy1out + dxy2out) / 2.0
            chrdict["{}-{}".format(i, j)] = dxyP / dxyout
    RNDdict[c] = chrdict
        # window = RNDdict[c][][2][1]
    return(RNDdict)


def RNDmin():
    """
    RNDmin: dmin / dout
    Similarly, like both dmin and Gmin, RNDmin should be sensitive to even rare
    migrant haplotypes. In addition, we expect RNDmin to be powerful even when
    migrants are high in frequency
    """
#    Dmin/Dout


def adaptPlot(Q, chrom, name, p, save=True):
    """
    """
    fig, ax = plt.subplots(figsize=(10, 4))
    sns.despine(ax=ax, offset=5)
    title = "{}-{}".format(name, chrom)
    nx = Q[1]
    x = [(np.sum(i)-1)/2 for i in nx]  # need midpoints
    y = Q[0]
    ax.plot(x, y, lw=.5, label=p)
    ax.set_ylabel(name)
    ax.set_xlabel('Chromosome {} position (bp)'.format(chrom))
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    if save:
        fig.savefig("{}.pdf".format(title), bbox_inches='tight')
    return(None)


if __name__ == "__main__":
    makeh5fromvcf(args.vcfFile, 1, args.h5)
    meta = "AnopSG.40.info"
    #meta = "AnopSG.50.info"
    meta = args.meta
    meta = pd.read_csv(meta, delimiter=",")
    var = Chr('All',"AnfunSG.liftover.biSNP.40.NoMiss.AiTests.recode.h5")
    #var = Chr('All',"../AnfunSG.liftover.biSNP.50.NoMiss.AiTests.recode.h5")
    var = Chr('All', "{}.h5".format(args.vcfFile))
    popdict = autil.subpops(var, meta, bypop=True, bykary=False)
    pop2color = autil.popcols(popdict)
    chrlist = np.unique(var.chrm[:])
    chrlen = {"3L": 46588628,
              "3R": 43534138,
              "2L": 44382598,
              "2R": 54617495,
              "X": 20138207}
    pops = list(popdict.values())
    outgrp, target, bait = args.out, args.target, args.bait
    # outgrp, target, bait = "Fun", "Van", "Par"
    # outgrp, target, bait = "Par", "Van", "Fun"
    # outgrp, target, bait = "Par", "Fun", "Like"
    # outgrp, target, bait = "Par", "Like", "Fun"
    # outgrp, target, bait = "Van", "Long", "Par"
    # outgrp, target, bait = "Par", "Long", "Van"
    # outgrp, target, bait = "Fun", "Par", "Like"
    name = "{}_{}_{}".format(outgrp, target, bait)
    for c in chrlist:
        var.geno(c, meta)
        Q95 = QninetyFive(var.pos, var.gt, popdict[outgrp], popdict[target], popdict[bait], chrlen[c])
        adaptPlot(Q95, c, name, 'Q95')
        U20 = Utwenty(var.pos, var.gt, popdict[outgrp], popdict[target], popdict[bait], chrlen[c])
        adaptPlot(U20, c, name, 'U20')
        vq = np.stack([Q95[1][:, 0], Q95[1][:, 1], Q95[0]], axis=1)
        vu = np.stack([U20[1][:, 0], U20[1][:, 1], U20[0]], axis=1)
        v = np.hstack((vq, vu))
        np.savetxt("{}.{}.stats.out".format(name, c), v, '%5.2f')
        Q9599 = np.quantile(Q95[0], 0.99)
        U2099 = np.quantile(U20[0], 0.99)
        jointQuants = (Q95[0] >= Q9599) * (U20[0] >= U2099)
        Aiout = Q95[1][jointQuants]
        n = [[name, c]] * len(Aiout)
        outBed = np.hstack((n, Aiout))
        f = open("{}.{}.windows.out".format(name, c), 'w')
        for i in outBed:
            f.write("{}\n".format("\t".join(i)))
        f.close()
        #np.savetxt("{}.{}.windows.out".format(name, c), outBed, delimiter="\t")
