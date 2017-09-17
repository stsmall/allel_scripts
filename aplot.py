#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 16:30:34 2017
Assorted functions for using scikit allel
@author: scott
"""
import allel
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from allel.util import asarray_ndim, check_dim0_aligned, ensure_dim1_aligned
from itertools import combinations


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


def plotstats(pc, title, pop2colors, saved=False):
    """
    """
    fig, ax = plt.subplots(figsize=(12, 4))
    sns.despine(ax=ax, offset=10)
    left = np.arange(len(pc))
    colors = [pop2colors[p] for p in pop2colors.keys()]
    ax.bar(left, pc, color=colors)
    ax.set_xlim(0, len(pc))
    ax.set_xlabel('Sample index')
    ax.set_ylabel('Percent calls')
    ax.set_title(title)
#    handles, labels = ax.get_legend_handles_labels()
#    ax.legend(handles, labels, title="Population", bbox_to_anchor=(1, 1),
#              loc='upper left')
    handles = ["mtl.patches.Patch(color={}, label={})".format(pop2colors[p], p)
               for p in pop2colors.keys()]
    ax.legend(handles=handles, labels=list(pop2colors.keys()),
              title='Population', bbox_to_anchor=(1, 1), loc='upper left')
    if saved:
        fig.savefig("{}.pdf".format(title), bbox_inches='tight')


def plotmiss(x, y, title, chrm, saved):
    """
    """
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
