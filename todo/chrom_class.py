#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 15:49:00 2019
@author: Scott T. Small

This module demonstrates documentation as specified by the `NumPy
Documentation HOWTO`_. Docstrings may extend over multiple lines. Sections
are created with a section header followed by an underline of equal length.

Example
-------
Examples can be given using either the ``Example`` or ``Examples``
sections. Sections support any reStructuredText formatting, including
literal blocks::

    $ python example_numpy.py


Section breaks are created with two blank lines. Section breaks are also
implicitly created anytime a new section starts. Section bodies *may* be
indented:

Notes
-----
    This is an example of an indented section. It's like any other section,
    but the body is indented to help it stand out from surrounding text.

If a section is indented, then a section break is created by
resuming unindented text.

Attributes
----------
module_level_variable1 : int
    Module level variables may be documented in either the ``Attributes``
    section of the module docstring, or in an inline docstring immediately
    following the variable.

    Either form is acceptable, but the two should not be mixed. Choose
    one convention to document module level variables and be consistent
    with it.


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
import allel
import numpy as np

class Chr(object):
    """Docstring


    Attributes
    ----------
    pop: str
        asdfasdf
    calls: str
        asdfklkasdlf



    """
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
        misscount = genotypes.count_missing(axis=1)
        # .20 will remove sites with > 20% missing data
        missarray = misscount <= (genotypes.n_samples * prctmiss)
        # remove missing
        self.gt = genotypes.compress(missarray, axis=0)
        # adjust positions
        self.pos = positions[missarray]

    def seg(self, genotypes, positions):
        """
        """
        var_seg = genotypes.count_alleles().is_segregating()
        self.pos = positions[var_seg]
        self.gt = genotypes.compress(var_seg)

    def mac(self, genotypes, positions, mac):
        """
        mac: minor allele count
        """
        # identify singletons
        allelecounts = genotypes.count_alleles()
        fltsingle = (allelecounts.max_allele() == 1) &\
            (allelecounts[:, :2].min(axis=1) > mac)  # biallelic only
        # remove singletons
        self.gt = genotypes.compress(fltsingle, axis=0)
        # adjust positions
        self.pos = positions[fltsingle]