#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 17:08:45 2018

@author: scott
"""
import matplotlib.pyplot as plt
import allel
import numpy as np
import h5py
chromlist = ["Wb_Chr1_0", "Wb_Chr2_0", "Wb_Chr2_1", "Wb_Chr2_2", "Wb_Chr3_1",
             "Wb_Chr4_0", "Wb_Chr4_1", "Wb_Chr4_2"]
chromlist = ["Wb_ChrX_0"]
window_size = 100
seldict = {}


def plth12(chromlist):
    """
    """
    for c in chromlist:
        # callset = h5py.File("PNG.phased.autosomal.recode.{}.h5".format(c), mode='r')
        callset = h5py.File("PNG.phased.X.recode.{}.h5".format(c), mode='r')
        samples = callset['samples'][:]
        sample_name = [sid.decode() for sid in samples.tolist()]
        g = allel.GenotypeChunkedArray(callset["calldata/GT"])
        h = g.to_haplotypes()
        pos = allel.SortedIndex(callset["variants/POS"][:])
        acc = h.count_alleles()[:, 1]
        # H12
        h12 = allel.moving_garud_h(h, window_size)[1]  # set window size
        h12_pos = []
        p = 0
        end = window_size
        i = 0
        while i < len(h12):
            stop = pos[end]
            while pos[p] < stop:
                h12_pos.append(h12[i])
                p += 1
            i += 1
            end += window_size
        while len(h12_pos) < len(pos):
            h12_pos.append(h12[-1])
        plt.plot(pos, h12_pos)
        plt.xlabel("{} genomic position".format(c))
        plt.ylabel("H12")
        plt.savefig("PNG.{}.H12.pdf".format(c))
        plt.clf()


def pltPi(chromlist):
    """
    """
    for c in chromlist:
        callset = h5py.File("PNG.phased.autosomal.recode.{}.h5".format(c), mode='r')
        # callset = h5py.File("PNG.phased.X.recode.{}.h5".format(c), mode='r')
        samples = callset['samples'][:]
        sample_name = [sid.decode() for sid in samples.tolist()]
        g = allel.GenotypeChunkedArray(callset["calldata/GT"])
        pos = allel.SortedIndex(callset["variants/POS"][:])
        acc = g.count_alleles()
        pi_windowed = allel.windowed_diversity(pos, acc, size=10)
        plt.plot(pos, h12_pos)
        plt.xlabel("{} genomic position".format(c))
        plt.ylabel("H12")
        plt.savefig("PNG.{}.H12.pdf".format(c))
        plt.clf()


plth12(chromlist)
