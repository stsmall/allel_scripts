#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 15:40:01 2017

@author: scott
"""

import allel
from collections import defaultdict
import numpy as np
import seaborn as sns


def catdict(chrdict):
    """
    """
    diversity = defaultdict(list)
    for c in chrdict.keys():
        for p in chrdict[c].keys():
            diversity[p].append(chrdict[c][p][2][0])
    divsum = {}
    for p in diversity.keys():
        x = np.concatenate(diversity[p])
        divsum[p] = x[~np.isnan(x)]
    return(divsum)


def popcols(popdict):
    """Set plot colors for each subpop
    """
    pop2color = {}
    palette = sns.color_palette(n_colors=len(popdict.keys()))
    for pop in popdict.keys():
        pop2color[pop] = palette.pop()
#    if "Kiribina" in popdict.keys() and "Folonzo" in popdict.keys():
#        pop2color['Kiribina'] = '#FF0000'
#        pop2color['Folonzo'] = '#008000'
    return(pop2color)


def subpops(var, meta, bypop=False, bykary=False):
    """Define subpops
    """
    popdict = {}
    if bykary and not bypop:
        if var.pop is not "All":
            meta = meta.ix[meta.Population.isin(var.pop)]
        else:
            pass  # meta = meta
        for kar in meta.ChromForm.unique():
            popdict[kar] = meta[meta.ChromForm == kar].index.tolist()
    elif bypop and not bykary:
        if var.pop is not "All":
            meta = meta.ix[meta.Population.isin(var.pop)]
        for pop in meta.Population.unique():
            popdict[pop] = meta[meta.Population == pop].index.tolist()
    elif bypop and bykary:
        popdict = defaultdict(dict)
        if var.pop is not "All":
            meta = meta.ix[meta.Population.isin(var.pop)]
        else:
            pass
        for kar in meta.ChromForm.unique():
            kpop = meta.ix[meta.ChromForm == kar]
            for pop in kpop.Population.unique():
                popdict[kar][pop] = kpop[kpop.Population == pop].index.tolist()
    else:
        popdict["BurkinaFaso"] = meta.index.tolist()
    return(popdict)


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


def ldthin(geno, positions, method, rsize=50000, mac=1, size=100, step=20,
           thresh=.1, iters=1):
    """.take if coord, .compress if mask
    """
    ac = geno.count_alleles()
    # only biallelic and mac or 2
    pca_selection = (ac.max_allele() == 1) & (ac[:, :2].min(axis=1) > mac)

    print("Available site for PCA: {}".format(np.count_nonzero(pca_selection)))
    if 'random' in method:
        indices = np.nonzero(pca_selection)[0]
        indices_ds = np.random.choice(indices, size=rsize, replace=False)
        indices_ds.sort()
        genotypes_pca = geno.take(indices_ds, axis=0)
        # sites with missing data can return error
        gn = genotypes_pca.to_n_alt()[:]
        pos = indices_ds
        print("{} Random SNPs selected for PCA".format(gn.shape[0]))
    else:
        genotypes_pca = geno.compress(pca_selection, axis=0)
        gn = genotypes_pca.to_n_alt()
        pos = positions[pca_selection]
        for i in range(iters):
            loc_unlinked = allel.locate_unlinked(gn, size=size, step=step,
                                                 threshold=thresh)
            n = np.count_nonzero(loc_unlinked)
            n_remove = gn.shape[0] - n
            print("iteration {} retaining {} removing {} variants".format(i+1,
                  n, n_remove))
            gn = gn.compress(loc_unlinked, axis=0)
            pos = pos[loc_unlinked]
    return(gn, allel.SortedIndex(pos))


def vcf2plink(thindict):
    """Takes a thinned file from allel, thins the original vcf to match, then
    converts the vcf to plink. Returns plink formatted vcf to run in ADMIXTURE
    and Treemix
    Requires VCFTOOLS and PLINK1.9
    """
    with open("thin.pos", 'w') as tpos:
        for nchr in thindict.keys():
            for pos in thindict[nchr]:
                tpos.write("{}\t{}\n".format(nchr.decode("utf-8", "strict"),
                                             pos))
#    if vcf:
#        command = "vcftools --vcf " + vcf + " --positions thin.pos --recode --out \
#        thinnedsites"
#        proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
#        proc.wait()
#        # call plink
#        command = "plink --make-bed thinnedsites.recode.vcf --allow-extra-chrom \
#        --out thinnedplink"
#        proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
#        proc.wait()
    return(None)
