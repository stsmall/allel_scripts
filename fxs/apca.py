#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 15:38:44 2017

@author: scott
"""

import allel
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
sns.set_style('white')
sns.set_style('ticks')


def plot_pca_coords(coords, model, pc1, pc2, ax, samples, pop2color):
    """
    """
    sns.despine(ax=ax, offset=5)
    x = coords[:, pc1]
    y = coords[:, pc2]
    for pop in np.unique(samples):
        flt = (samples == pop)
        ax.plot(x[flt], y[flt], marker='o', linestyle=' ',
                color=pop2color[pop], label=pop, markersize=6, mec='k',
                mew=.5)
    ax.set_xlabel('PC%s (%.1f%%)' % (pc1+1,
                  model.explained_variance_ratio_[pc1]*100))
    ax.set_ylabel('PC%s (%.1f%%)' % (pc2+1,
                  model.explained_variance_ratio_[pc2]*100))


def fig_pca(coords, model, title, samples, pop2color, pc1, pc2):
    """
    """
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 2, 1)
    plot_pca_coords(coords, model, pc1, pc2, ax, samples, pop2color)
    ax.legend(bbox_to_anchor=(1, 1), loc='upper left')
    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    fig.savefig("{}.pdf".format(title), bbox_inches='tight')


def hist_var(model, title):
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
    fig.savefig("{}_histvar.pdf".format(title), bbox_inches='tight')


def pca_fx(gnu, meta, nchr, pop2color, pcAll, population, bykary=False):
    """
    gnu: genotype object transformed with .to_n_alt()
    pcAll:
    """
    # PCA
    coords, model = allel.pca(gnu, n_components=10, scaler='patterson')
#    corrds, model = allel.randomized_pca(gnu, n_components=10,
#                                         scaler='patterson')
    title = "PCA Chr:{}, var:{}".format(nchr, gnu.shape[0])
    if population is 'All':
        samples = meta.Population.values
    else:
        samples = meta.Population[meta.Population.isin(population)].values
    if bykary:
        s = meta.Population[meta.Population.isin(samples)].index.tolist()
        samples = meta.ChromForm[s].values
        pop2color['Kiribina'] = '#FF0000'
        pop2color['Folonzo'] = '#008000'
    if pcAll:
        i = 0
        while i < 9:
            fig_pca(coords, model, title, samples, pop2color, i, i+1)
            i += 2
    else:
        fig_pca(coords, model, title, samples, pop2color, 0, 1)
    hist_var(model, title)
    return(coords, model)
