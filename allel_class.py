#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 14:24:53 2017
Class for scikit-allel
Example:
Chr1_0 = Chr(callset, "Chr1_0")
Chr1_0.loadpos()
    Chr1_0.positions  # positions
    Chr1_0.vtbl  # variant table
    Chr1_0.genotypes  # genotypes
    Chr1_0.allelecounts  # allelecounts
Chr1_0.popsubset(["PNG", "Haiti", "Kenya", "Mali"])
    Chr1_0.subsample["PNG"]  # meta just in this pop
    Chr1_0.subsample_genotype["PNG"]  # genotypes just in this pop
Chr1_0.ldfilter(Chr1_0.genotypes, size=100, step=50, thresh=.1, n_iter=1)
    Chr1_0.ldthin
    Chr1_0.corr
    Chr1_0.gtNOsingletons

@author: ssmall2@nd.edu
"""
import numpy as np
import allel
import h5py


class Chr(object):
    def __init__(self, pop, calls):
        """
        name: chromosome to query
        calls: load data
            callset_fn = FOO.h5
            callset = h5py.File(callset_fn, mode='r')
            Chr1 = Chr("Wb_Chr1_0", callset)
        """
        self.calls = h5py.File(calls, mode='r')  # requires premade h5
        self.pop = pop
        self.chrm = self.calls['variants/CHROM']

    def geno(self, name, meta):
        """
        pos: positions in the chromosome
        gt: genotype array
        sm: attached meta data
        """
        # get chrom
        chrom = self.chrm
        chrom_mask = np.where(chrom[:] == name)
        p = self.calls['variants/POS']
        pos = p[:][chrom_mask]
        self.pos = allel.SortedIndex(pos)
        # get genotype data
        geno = allel.GenotypeArray(self.calls['calldata/GT'])
        gt = geno[:][chrom_mask]
        # get population data
        if self.pop == 'All':
            self.gt = gt
            self.sm = meta
        else:
            samples = meta.ix[meta.Population.isin(self.pop)].index.tolist()
            self.gt = gt[:, samples]
            self.sm = meta.ix[samples]

    def miss(self, genotypes, positions, prctmiss):
        """
        """
#        missarray = np.ones(genotypes.n_variants, dtype=bool)
#        is_missing = genotypes.is_missing()
#        for i, miss in enumerate(is_missing):
#            if sum(miss)/float(genotypes.n_samples) > prctmiss:
#                # prctmiss = 0 for non-missing sites
#                missarray[i] = False
        misscount = genotypes.count_missing(axis=1)
        missarray = misscount <= (genotypes.n_samples * prctmiss)
        # remove missing
        self.gt_m = genotypes.compress(missarray, axis=0)
        # adjust positions
        self.pos_m = positions[missarray]

    def seg(self, genotypes, positions):
        """
        """
        var_seg = genotypes.count_alleles().is_segregating()
        self.pos_s = positions[var_seg]
        self.gt_s = genotypes.compress(var_seg)

    def mac(self, genotypes, positions, mac):
        """
        mac: minor allele count
        """
        # identify singletons
        allelecounts = genotypes.count_alleles()
        fltsingle = (allelecounts.max_allele() == 1) &\
            (allelecounts[:, :2].min(axis=1) > mac)  # biallelic only
        # remove singletons
        self.gt_c = genotypes.compress(fltsingle, axis=0)
        # adjust positions
        self.pos_c = positions[fltsingle]
