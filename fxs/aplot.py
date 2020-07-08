#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 16:30:34 2017
Assorted functions for using scikit allel
@author: scott
"""
import allel
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
mpl.rcParams['pdf.fonttype'] = 42


def plotvars2(chrm, ac_subpops, pos, pop2color, window_size=100000, title=None, saved=True):
    """
    """
    fig, ax = plt.subplots(figsize=(10, 4))
    sns.despine(ax=ax, offset=5)
    poplist = list(ac_subpops.keys())
    bins = np.arange(0, pos.max(), window_size)
    for pops in poplist:
        print(pops)
        acu = ac_subpops[pops]
        flt = acu.is_segregating()
        varpos = pos[flt]
        # use window midpoints as x coordinate
        x = (bins[1:] + bins[:-1])/2
        # compute variant density in each window
        h, _ = np.histogram(varpos, bins=bins)
        y = h / window_size
        ax.plot(x, y, lw=.5, label=pops, color=pop2color[pops])
    ax.set_xlabel('Chromosome position (bp)')
    ax.set_ylabel('Variant density (bp$^{-1}$)')
    ##### plot legend
    handle = []
    for p in poplist:
        handle.append(mlines.Line2D([], [], color=pop2color[p], label=p))
    ax.legend(handles=handle, title='Population', bbox_to_anchor=(1, 1), loc='upper left')    
    #####
    if title:
        ax.set_title(title)
    else:
        ax.set_title(chrm)
    if saved:
        fig.savefig("{}.vars.pdf".format(chrm), bbox_inches='tight')
    return(None)


def plotvars(chrm, callset, window_size=100000, title=None, saved=True):
    """
    """
    try:
        chrm = chrm.decode("utf-8")
    except AttributeError:
        chrm = chrm
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
        ax.set_title(chrm)
    if saved:
        fig.savefig("{}.vars.pdf".format(chrm), bbox_inches='tight')


def plotstats(pc, title, pop2color, saved=False):
    """
    """
    fig, ax = plt.subplots(figsize=(12, 4))
    sns.despine(ax=ax, offset=10)
    left = np.arange(len(pc))
    poplist = list(pop2color.keys())
    colors = [pop2color[p] for p in poplist]
    ax.bar(left, pc, color=colors)
    ax.set_xlim(0, len(pc))
    ax.set_xlabel('Sample index')
    ax.set_ylabel('Percent calls')
    ax.set_title(title)
    handle = []
    for p in poplist:
        handle.append(mlines.Line2D([], [], color=pop2color[p], label=p))
    ax.legend(handles=handle, title='Population', bbox_to_anchor=(1, 1), loc='upper left')    
    if saved:
        fig.savefig("{}.pdf".format(title), bbox_inches='tight')


def plotmiss(x, y, title, pop2color, chrm, saved):
    """
    """
    # plot
    try:
        chrm = chrm.decode("utf-8")
    except AttributeError:
        chrm = chrm
    poplist = list(pop2color.keys())
    fig, ax = plt.subplots(figsize=(12, 3))
    sns.despine(ax=ax, offset=10)
    colors = [pop2color[p] for p in poplist]
    ax.plot(x, y, color=colors)
    handle = []
    for p in poplist:
        handle.append(mlines.Line2D([], [], color=pop2color[p], label=p))
    ax.legend(handles=handle, title='Population', bbox_to_anchor=(1, 1), loc='upper left')    
    ax.set_xlabel('Chromosome position (bp)')
    ax.set_ylabel('Missing Count')
    if title:
        ax.set_title(title)
    else:
        ax.set_title(chrm)
    if saved:
        fig.savefig("{}.miss.pdf".format(chrm),
                    bbox_inches='tight')


def divboxplot(diversity, pop2color):
    """
    """
    poplist = list(pop2color.keys())
    fig, ax = plt.subplots(figsize=(4, 4))
    sns.despine(ax=ax, offset=5)
    lw = 1.5
    box = ax.boxplot(
        x=[diversity[pop] for pop in poplist],
        labels=poplist,  patch_artist=True,
        medianprops={"color": "k", "linewidth": lw},
        whiskerprops={"color": "k"},
        capprops={"color": "k"},
        showfliers=False,
        flierprops={"c": "k", "markersize": 2})
    ax.set_ylabel(r'$\pi$', rotation=0, fontsize=16)
    for patch, color in zip(box['boxes'], [pop2color[pop] for pop in poplist]):
        patch.set_facecolor(color)
        patch.set_linewidth(lw)
    fig.savefig("piboxplot.pdf", bbox_inches="tight")
