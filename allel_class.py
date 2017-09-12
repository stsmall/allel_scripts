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


class Chr(object):
    def __init__(self, name, pop, calls):
        """
        name: chromosome to query
        calls: load data
            callset_fn = FOO.h5
            callset = h5py.File(callset_fn, mode='r')
            Chr1 = Chr("Wb_Chr1_0", callset)
        """
        self.name = name
        self.calls = calls
        self.pop = pop

    def loadpos(self, meta):
        """
        pos: positions in the chromosome
        gt: genotype array
        sm: attached meta data
        """
        # get chrom
        chrom = self.calls['variants/CHROM']
        chrom_mask = np.where(chrom[:] == self.name)
        # get chrom variant positions
        p = self.calls['variants/POS']
        pos = p[:][chrom_mask]
        self.pos = allel.SortedIndex(pos)
        # get genotype data
        geno = allel.GenotypeArray(self.calls['calldata/GT'])
        gt = geno[:][chrom_mask]
        # get population data
        if self.pop == 'ALL':
            self.gt = gt
            self.sm = meta
        elif type(p) is list:
            samples = []
            for p in self.pop:
                samples.append(np.where(meta.Population == p)[0])
            self.gt = gt[:, samples]
            self.sm = meta[samples]
        else:
            samples = np.where(meta.Population == self.pop)[0]
            self.gt = gt[:, samples]
            self.sm = meta[samples]

    def missing(self, genotypes, positions, prctmiss):
        """
        """
        missarray = np.ones(genotypes.n_variants, dtype=bool)
        is_missing = genotypes.is_missing()
        for i, miss in enumerate(is_missing):
            if sum(miss)/float(genotypes.n_samples) > prctmiss:
                missarray[i] = False
        # remove missing
        self.gt_nm = genotypes.compress(missarray, axis=0)
        # adjust positions
        self.pos_nm = positions[missarray]

    def isseg(self, genotypes, positions):
        """
        """
        var_seg = genotypes.count_alleles().is_segregating()
        self.pos_s = positions[var_seg]
        self.gt_s = genotypes.compress(var_seg)

    def rmvsingletons(self, genotypes, positions, mac):
        """
        mac: minor allele count
        """
        # identify singletons
        allelecounts = genotypes.count_alleles()
        fltsingle = (allelecounts.max_allele() == 1) &\
            (allelecounts[:, :2].min(axis=1) > mac)  # biallelic only
        # remove singletons
        self.gt_mac = genotypes.compress(fltsingle, axis=0)
        # adjust positions
        self.pos_mac = positions[fltsingle]

    def ldfilter(self, genotypes, positions, size, step, thresh,
                 n_iter):
        """
        genotype: genotype array; self.genotypes
        """
        # here it is now 2d, so no longer genotypeClass
        self.corr = genotypes.to_n_alt()
        self.ldpos = positions
        for i in range(n_iter):
            loc_unlinked = allel.locate_unlinked(self.corr, size=size,
                                                 step=step, threshold=thresh)
            n = np.count_nonzero(loc_unlinked)
            n_remove = self.corr.shape[0] - n
            print('iteration', i+1, 'retaining', n, 'removing',
                  n_remove, 'variants')
            self.corr = self.corr.compress(loc_unlinked, axis=0)
            self.ldpos = self.ldpos[loc_unlinked]
