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
import numpy as np
import pandas as pd
from jackknife import jackknife
from sfs import jsfs_fx
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
    fig.savefig("{}.pdf".format(title), bbox_inches='tight')


def hist_var(model):
    """
    """
    fig, ax = plt.subplots()
    sns.despine(ax=ax)
    x = np.arange(10)
    y = model.explained_variance_ratio_ * 100
    ax.bar(x+.6, y, width=.8)
    ax.set_xticks(x+1)
    ax.set_xlim(0, 11)
    ax.set_xlabel('component')
    ax.set_ylabel('% variance explained')
    fig.savefig("PCA_histvar.pdf", bbox_inches='tight')


def pca_fx(gtthin, df_samples, pcall, nchr):
    """
    """
    gnu = gtthin[:]  # uncompress
    # PCA
    coords, model = allel.pca(gnu, n_components=10, scaler='standard')
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
    hist_var(model)
    return(coords, model)


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


def df_fst(ac1, ac2, pos, chrom, blen, pair):
    """
    """
    try:
        fst_hudson, se_hudson, vb_hudson, _ =\
            allel.stats.blockwise_hudson_fst(ac1, ac2, blen=blen)
        # fst_hudson, se_hudson, vb_hudson, _ =\
        # allel.stats.windowed_hudson_fst(pos, ac1, ac2, size=blen)
        # use the per-block average Fst as the Y coordinate
        y = vb_hudson
        # use the block centres as the X coordinate
        x = allel.stats.moving_statistic(pos,
                                         statistic=lambda v: (v[0] + v[-1])/2,
                                         size=blen)
        # chrom, pair
        c = [chrom] * len(x)
        p = [pair] * len(x)
    except ZeroDivisionError:
        x = [0]
        y = [0]
        c = [0]
        p = [0]
    return(x, y, c, p)


def df_dxy(ac1, ac2, pos, chrom, blen, pair):
    """
    """
    try:
        dxy = allel.stats.diversity.mean_pairwise_difference_between(ac1, ac2)
        dxy_blocks = allel.stats.moving_statistic(dxy, statistic=np.mean,
                                                  size=blen)
        y = dxy_blocks
        # use the block centres as the X coordinate
        x = allel.stats.moving_statistic(pos,
                                         statistic=lambda v: (v[0] + v[-1])/2,
                                         size=blen)
        # chrom, pair
        c = [chrom] * len(x)
        p = [pair] * len(x)
    except ZeroDivisionError:
        x = [0]
        y = [0]
        c = [0]
        p = [0]
    return(x, y, c, p)


def pairFST_fx(ac1, ac2, pop1, pop2, dip, prct, ac_subpopsdict, subpops, pair,
               chrpos, blenW, genotype, wc, pop1_idx, pop2_idx):
    """
    """
    pairslist = []
    chromlist = []
    poslist = []
    fstlist = []
    # weir and cockerham
    if wc:
        a, b, c = allel.stats.weir_cockerham_fst(genotype,
                                                 subpops=[pop1_idx, pop2_idx],
                                                 max_allele=1)
        snp_fst_wc = (a / (a + b + c))[:, 0]  # dist FST
        snp_fst_wc = snp_fst_wc[~np.isnan(snp_fst_wc)]
        fst_wc, se_wc, vb_wc, _ =\
            allel.stats.blockwise_weir_cockerham_fst(genotype,
                                                     subpops=[pop1_idx,
                                                              pop2_idx],
                                                     blen=2000, max_allele=1)
        print('%s-%s : %.04f +/- %.04f (Weir & Cockerham)' % (pop1, pop2,
                                                              fst_wc, se_wc))

    # hudson FST
    fst_hudson, se_hudson, vb_hudson, _ =\
        allel.stats.blockwise_hudson_fst(ac1, ac2, blen=blenW)
    print('{}-{} : {:.4f} +/- {:.4f} (Hudson)'.format(pop1, pop2, fst_hudson,
          se_hudson))
    # pairwise FST per chrom
    for nchr in ac_subpopsdict.keys():
        acu = allel.AlleleCountsArray(ac_subpopsdict[nchr][pop1][:] +
                                      ac_subpopsdict[nchr][pop2][:])
        miss = (acu[:, 0] + acu[:, 1]) > (dip * prct)
        seg = acu.is_segregating()
        flt = miss * seg
        blenWc = int(round(sum(flt)/5))
        ac1c = allel.AlleleCountsArray(ac_subpopsdict[nchr][pop1].compress(flt,
                                       axis=0)[:, :2])
        ac2c = allel.AlleleCountsArray(ac_subpopsdict[nchr][pop2].compress(flt,
                                       axis=0)[:, :2])
        pos = chrpos[nchr][flt]
        # hudson
        fst_hudson, se_hudson, vb_hudson, _ =\
            allel.stats.blockwise_hudson_fst(ac1c, ac2c, blen=blenWc)
        print("{}, {}, Hudson's Fst: {:.4f} +/- {:.4f}".format(nchr, pair,
              fst_hudson, se_hudson))
        # df for plotting across chrom
        x, y, c, p = df_fst(ac1c, ac2c, pos, nchr, 1000, pair)  # blen iswindow
        pairslist.extend(p)
        chromlist.extend(c)
        poslist.extend(x)
        fstlist.extend(y)

    return(pairslist, chromlist, poslist, fstlist)


def pairDXY_fx(ac1, ac2, pop1, pop2, dip, prct, ac_subpopsdict, subpops, pair,
               chrpos, blenW, genotype, wc):
    """
    """
    pairslist = []
    chromlist = []
    poslist = []
    dxylist = []
    # DXY
    f2 = allel.stats.admixture.patterson_f2(ac1, ac2)
    f2_blocks = allel.stats.moving_statistic(f2, statistic=np.mean, size=blenW)
    f2_m, f2_se, _, _ = jackknife(f2_blocks)
    print('{}-{} : {:.4f} +/- {:.4f} (f2)'.format(pop1, pop2, f2_m, f2_se))
    # DXY
    dxy = allel.stats.diversity.mean_pairwise_difference_between(ac1, ac2)
    dxy_blocks = allel.stats.moving_statistic(dxy, statistic=np.mean,
                                              size=blenW)
    dxy_m, dxy_se, _, _ = jackknife(dxy_blocks)
    print('{}-{} : {:.4f} +/- {:.4f} (Dxy)'.format(pop1, pop2, dxy_m, dxy_se))
    # pairwise Dxy per chrom
    for nchr in ac_subpopsdict.keys():
        acu = allel.AlleleCountsArray(ac_subpopsdict[nchr][pop1][:] +
                                      ac_subpopsdict[nchr][pop2][:])
        miss = (acu[:, 0] + acu[:, 1]) > (dip * prct)
        flt = miss
        blenWc = int(round(sum(flt)/5))
        ac1c = allel.AlleleCountsArray(ac_subpopsdict[nchr][pop1].compress(flt,
                                       axis=0)[:, :2])
        ac2c = allel.AlleleCountsArray(ac_subpopsdict[nchr][pop2].compress(flt,
                                       axis=0)[:, :2])
        pos = chrpos[nchr][flt]
        # f2
        f2 = allel.stats.admixture.patterson_f2(ac1c, ac2c)
        f2_blocks = allel.stats.moving_statistic(f2, statistic=np.mean,
                                                 size=blenWc)
        f2_m, f2_se, _, _ = jackknife(f2_blocks)
        print("{}, {}, f2: {:.4f} +/- {:.4f}".format(nchr, pair, f2_m, f2_se))
        # dxy
        dxy = allel.stats.diversity.mean_pairwise_difference_between(ac1c,
                                                                     ac2c)
        dxy_blocks = allel.stats.moving_statistic(dxy, statistic=np.mean,
                                                  size=blenWc)
        dxy_m, dxy_se, _, _ = jackknife(dxy_blocks)
        print("{}, {}, Dxy: {:.4f} +/- {:.4f}".format(nchr, pair, dxy_m,
              dxy_se))

        # df for plotting across chrom
        x, y, c, p = df_dxy(ac1c, ac2c, pos, nchr, 1000, pair)  # bleniswindow
        pairslist.extend(p)
        chromlist.extend(c)
        poslist.extend(x)
        dxylist.extend(y)
    return(pairslist, chromlist, poslist, dxylist)


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
    print("{}".format(d2))
    print("{}".format(dshare))
    return(d2, dshare)


def Fstatistics_fx(ac_subpopscat, subpops):
    """
    """
    poplist = ["Haiti", "Mali", "Kenya", "PNG"]
    pop1 = poplist[0]
    pop2 = poplist[1]
    pop3 = poplist[2]
    pop4 = poplist[3]
    pop1_idx = subpops[pop1]
    pop2_idx = subpops[pop2]
    pop3_idx = subpops[pop3]
    pop4_idx = subpops[pop4]
    acu = allel.AlleleCountsArray(ac_subpopscat[pop1][:] +
                                  ac_subpopscat[pop2][:] +
                                  ac_subpopscat[pop3][:] +
                                  ac_subpopscat[pop4][:])
    dip = (len(pop1_idx) + len(pop2_idx) + len(pop3_idx) + len(pop4_idx)) * 2.0
    # prc missing
    prct = .80
    miss = (acu[:, 0] + acu[:, 1]) > (dip * prct)
    seg = acu.is_segregating()
    flt = miss * seg
    blen = int(round(sum(flt)/5))
    aca = allel.AlleleCountsArray(ac_subpopscat[pop1].compress(flt,
                                  axis=0)[:, :2])
    acb = allel.AlleleCountsArray(ac_subpopscat[pop2].compress(flt,
                                  axis=0)[:, :2])
    acc = allel.AlleleCountsArray(ac_subpopscat[pop3].compress(flt,
                                  axis=0)[:, :2])
    acd = allel.AlleleCountsArray(ac_subpopscat[pop4].compress(flt,
                                  axis=0)[:, :2])
    # acc = test, aca = pop1, acb = pop2, acd = outgroup
    # Patternson F3 admixture
    f3, f3_se, f3zscore, _, _ = allel.stats.blockwise_patterson_f3(aca, acb,
                                                                   acc, blen)
    print("{}: {}, {}, {}".format([pop1, pop2, pop3], f3, f3_se, f3zscore))
    f3, f3_se, f3zscore, _, _ = allel.stats.blockwise_patterson_f3(acb, acc,
                                                                   aca, blen)
    print("{}: {}, {}, {}".format([pop2, pop3, pop1], f3, f3_se, f3zscore))
    f3, f3_se, f3zscore, _, _ = allel.stats.blockwise_patterson_f3(acc, aca,
                                                                   acb, blen)
    print("{}: {}, {}, {}".format([pop3, pop1, pop2], f3, f3_se, f3zscore))
    # Patterson F4
    # Mali is admixted of Kenya
    f4, f4_se, f4zscore, _, _ = allel.stats.blockwise_patterson_d(aca, acb,
                                                                  acc, acd,
                                                                  blen)
    print("{}: {}, {}, {}".format([pop1, pop2, pop3, pop4], f4, f4_se,
          f4zscore))
    # Kenya is admixed of Mali
    f4, f4_se, f4zscore, _, _ = allel.stats.blockwise_patterson_d(aca, acc,
                                                                  acb, acd,
                                                                  blen)
    print("{}: {}, {}, {}".format([pop1, pop3, pop2, pop4], f4, f4_se,
          f4zscore))
    f4, f4_se, f4zscore, _, _ = allel.stats.blockwise_patterson_d(acb, acc,
                                                                  aca, acd,
                                                                  blen)
    print("{}: {}, {}, {}".format([pop2, pop3, pop1, pop4], f4, f4_se,
          f4zscore))

    return(None)


def allelecount_pair(ac_subpopsdict, ac_subpopscat, chrcat_gt, subpops,
                     pairlist, jsfs_r, chrpos):
    """
    """
    pairslist_fst = []
    chromlist_fst = []
    poslist_fst = []
    fstlist = []
    pairslist_dxy = []
    chromlist_dxy = []
    poslist_dxy = []
    dxylist = []
    for pair in pairlist:
        pop1 = pair.split("_")[0]
        pop2 = pair.split("_")[1]
        pop1_idx = subpops[pop1]
        pop2_idx = subpops[pop2]
        acu = allel.AlleleCountsArray(ac_subpopscat[pop1][:] +
                                      ac_subpopscat[pop2][:])
        dip = (len(pop2_idx) + len(pop1_idx)) * 2.0
        # prc missing
        prct = .80
        miss = (acu[:, 0] + acu[:, 1]) > (dip * prct)
        seg = acu.is_segregating()
        flt_s = miss * seg
        flt_m = miss

        # FST
        ac1fst = allel.AlleleCountsArray(ac_subpopscat[pop1].compress(flt_s,
                                         axis=0)[:, :2])
        ac2fst = allel.AlleleCountsArray(ac_subpopscat[pop2].compress(flt_s,
                                         axis=0)[:, :2])
        print("{}".format(sum(flt_s)))
        genotype = chrcat_gt.compress(flt_s, axis=0)
        blenW = int(round(sum(flt_s)/10))
        pf, cf, xf, fst = pairFST_fx(ac1fst, ac2fst, pop1, pop2, dip, prct,
                                     ac_subpopsdict, subpops, pair, chrpos,
                                     blenW, genotype, False, pop1_idx,
                                     pop2_idx)
        # DXY, F2
        ac1 = allel.AlleleCountsArray(ac_subpopscat[pop1].compress(flt_m,
                                      axis=0)[:, :2])
        ac2 = allel.AlleleCountsArray(ac_subpopscat[pop2].compress(flt_m,
                                      axis=0)[:, :2])
        print("{}".format(sum(flt_m)))
        genotype = chrcat_gt.compress(flt_m, axis=0)
        blenW = int(round(sum(flt_m)/10))
        pdx, cdx, xdx, dxy = pairDXY_fx(ac1, ac2, pop1, pop2, dip, prct,
                                        ac_subpopsdict, subpops, pair, chrpos,
                                        blenW, genotype, False)

        # jsfs
        if jsfs_r:
            jsfs = allel.stats.joint_sfs(ac1[:, 1], ac2[:, 1])
            jsfsf = allel.stats.joint_sfs_folded(ac1[:, :2], ac2[:, :2])
            jsfs_fx(jsfs, jsfsf, pop1, pop2)

        pairslist_fst.extend(pf)
        chromlist_fst.extend(cf)
        poslist_fst.extend(xf)
        fstlist.extend(fst)
        pairslist_dxy.extend(pdx)
        chromlist_dxy.extend(cdx)
        poslist_dxy.extend(xdx)
        dxylist.extend(dxy)
    fstchrom_plot = pd.DataFrame({"pair": pairslist_fst,
                                  "chrom": chromlist_fst,
                                  "pos": poslist_fst,
                                  "fst": fstlist,
                                  })
    fstchrom_plot.to_csv("WbFSTchrom.csv")
    fstchrom_plot = pd.DataFrame({"pair": pairslist_dxy,
                                  "chrom": chromlist_dxy,
                                  "pos": poslist_dxy,
                                  "fst": dxylist,
                                  })
    fstchrom_plot.to_csv("WbDXYchrom.csv")

    return(None)
