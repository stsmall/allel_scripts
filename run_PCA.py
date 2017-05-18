#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 18:01:30 2017
Plots PCA for scikit-allel
@author: scott
"""

import matplotlib.pyplot as plt
import seaborn as sns
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
    fig.savefig("PCA_ax{}_ax{}.pdf".format(pc1, pc2), bbox_inches='tight')
