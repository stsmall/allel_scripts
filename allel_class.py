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
    def __init__(self, name, calls):
        """
        name: chromosome to query
        calls: load data
            callset_fn = FOO
            callset = h5py.File(callset_fn, mode='r')
            Chr1 = Chr("Wb_Chr1_0", callset)
        """
        self.name = name
        self.calls = calls

    def loadpos(self, gt):
        """
        positions: positions in the chromosome
        vtbl: variant table with info data
        genotypes: genotype data
        """
        self.positions = allel.SortedIndex(self.calls[self.name]
                                           ["variants"]["POS"])
        self.vtbl = allel.VariantChunkedTable(self.calls[self.name]["variants"], names=["POS", "REF", 'ALT'])
        if not gt:
            self.haplotypes = allel.HaplotypeChunkedArray(self.calls[self.name]["calldata"])
            self.genotypes = self.haplotypes.to_genotypes(2)
        else:
            self.genotypes = allel.GenotypeChunkedArray(self.calls[self.name]["calldata"]["genotype"])
#        self.exonicPositions, _ = self.positions.locate_intersection_ranges(
#                features[features.seqid == self.name.encode('ascii')].start,
#                features[features.seqid == self.name.encode('ascii')].end)
#        # pull out exonic genotypes and genotype qualities
#        self.genotypes_exonic = allel.GenotypeArray(self.calls[self.name]
#                                                    ["calldata/genotype"]
#                                                    [self.exonicPositions, :])
#        self.vtbl_exonic = self.vtbl[:][self.exonicPositions]

    def missing(self, genotypes, positions, prctmiss):
        """
        """
        missarray = np.ones(genotypes.n_variants, dtype=bool)
        is_missing = genotypes.is_missing()[:]
        for i, miss in enumerate(is_missing):
            if sum(miss)/float(genotypes.n_samples) > prctmiss:
                missarray[i] = False
        # remove missing
        self.gtnomiss = genotypes.compress(missarray, axis=0)
        # adjust positions
        self.posnomiss = positions[missarray]

    def isseg(self, genotypes, positions):
        """
        """
        var_seg = genotypes[:, :].count_alleles().is_segregating()
        self.posisseg = positions[var_seg]
        self.gtisseg = genotypes.compress(var_seg)

    def rmvsingletons(self, genotypes, positions):
        """
        """
        # identify singletons
        allelecounts = genotypes.count_alleles()[:]
        fltsingle = (allelecounts.max_allele() == 1) &\
            (allelecounts[:, :2].min(axis=1) > 1)
        # remove singletons
        self.gtNOsingletons = genotypes.compress(fltsingle, axis=0)
        # adjust positions
        self.posNOsingletons = positions[fltsingle]

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
