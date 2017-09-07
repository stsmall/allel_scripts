#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 16:09:58 2017

@author: scott
"""
import allel
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
sns.set_style('white')
sns.set_style('ticks')


def sfs_fx(ac_subpopscat):
    """
    note: should filter on segregating if only using subset of pops
    note: only biallelic if >1 allele
    is_biallelic_01 = ac_seg['all'].is_biallelic_01()[:]
    ac1 = ac_seg['BFM'].compress(is_biallelic_01, axis=0)[:, :2]
    ac2 = ac_seg['AOM'].compress(is_biallelic_01, axis=0)[:, :2]
    """
    ac1 = ac_subpopscat["PNG"]
    ac2 = ac_subpopscat["Haiti"]
    ac3 = ac_subpopscat["Mali"]
    ac4 = ac_subpopscat["Kenya"]
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.despine(ax=ax, offset=10)
    sfs1 = allel.stats.sfs_folded_scaled(ac1)
    allel.stats.plot_sfs_folded_scaled(sfs1, ax=ax, label='PNG',
                                       n=ac1.sum(axis=1).max())
    sfs2 = allel.stats.sfs_folded_scaled(ac2)
    allel.stats.plot_sfs_folded_scaled(sfs2, ax=ax, label='Haiti',
                                       n=ac2.sum(axis=1).max())
    sfs3 = allel.stats.sfs_folded_scaled(ac3)
    allel.stats.plot_sfs_folded_scaled(sfs3, ax=ax, label='Mali',
                                       n=ac3.sum(axis=1).max())
    sfs4 = allel.stats.sfs_folded_scaled(ac4)
    allel.stats.plot_sfs_folded_scaled(sfs4, ax=ax, label='Kenya',
                                       n=ac4.sum(axis=1).max())
    ax.legend()
    ax.set_title('Scaled folded site frequency spectra')
    ax.set_xlabel('minor allele frequency')
    fig.savefig("ScaledSFS.pdf", bbox_inches='tight')
    return(None)


def jsfs_fx(jsfs, jsfsf, pop1, pop2):
    """
    """
    fig, ax = plt.subplots(figsize=(6, 6))
    allel.stats.plot_joint_sfs(jsfs, ax=ax)
    ax.set_ylabel('Alternate allele count, {}'.format(pop1))
    ax.set_xlabel('Alternate allele count, {}'.format(pop2))
    fig.savefig("jSFS.{}-{}.pdf".format(pop1, pop2), bbox_inches='tight')
#    fig, ax = plt.subplots(figsize=(6, 6))
#    allel.stats.plot_joint_sfs_folded(jsfsf, ax=ax)
#    ax.set_ylabel('Alternate allele count, {}'.format(pop1))
#    ax.set_xlabel('Alternate allele count, {}'.format(pop2))
#    fig.savefig("jSFS-fold.{}-{}.pdf".format(pop1, pop2), bbox_inches='tight')
    return(None)
