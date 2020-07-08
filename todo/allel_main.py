#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 15:10:39 2019
@author: Scott T. Small

The main code for running modules I have written associated with genomic
analysis using scikit-allel.

Example
-------
This script assumes that your vcf has been converted into zarr format.
Scikit-allel example for vcf_to_zarr() function found,
    `here <http://alimanfoo.github.io/2018/04/09/selecting-variants.html>`_.

    $ python allel_main.py -z FOO.zarr -m metaData.txt --help


Notes zarr to vcf
-----
    import allel
    import zarr
    vcf_path = "/dir/foo.vcf.gz
    zarr_path =/dir/foo.zarr
    allel.vcf_to_zarr(vcf_path, zarr_path, group='X', fields='*', log=sys.stdout, overwrite=True)
    # to view
    callset = zarr.open_group(zarr_path, mode='r')
    callset.tree(expand=True)


Fx notes
----------
    # load data
    callset = zarr.open_group(args.zarr, mode='r')
    callset.tree(expand=True)  # zarr tree
    # pos
    pos = allel.SortedIndex(callset["group/variants/POS"])
    loc_region = pos.locate_range(Start, Stop)  # genomic coords to array coords
    # genotype array
    gt_zarr = callset["{}/calldata/GT".format(args.group)]
    gt_region = allel.GenotypeArray(gt_zarr[loc_region])  # subset
    gt = allel.GenotypeArrya(gt_zarr)  # all data
    gt_variant_selection = gt.compress(loc_variant_selection, axis=0) # to select from gt array
    # Dask uses less memory but need to use compute()
    gt_dask = allel.GenotypeDaskArray(gt_zarr)
    gt_variant_selection = gt_dask.compress(loc_variant_selection, axis=0).compute()


    ac = gt.count_alleles(max_allele=2).compute
    gt = allel.GenotypeArray(callset['calldata/GT'])
    gt[variants, samples], gt[1:3, :] 2-4th variant all samples
    gt[:, 0:2], all variants, 1st and 2nd samples
    gt.compress(mask)
    gt.take(actual_indexes)
    gt.subset()

Attributes
----------
allel_class : class
    contains main class Chr. Chr calculates basic statistics about the data
    including missing sites, missing individuals,
popgen_class : class
    all functions that use population genetic data
plotting_class : class
    plotting functions for popgen and others
simulate_class : class
    simulations using msprime for when bootstrapping

"""
# import modules
import sys
import allel
import zarr
import numcodecs
import msprime
from bokeh.plotting import figure, show, ColumnDataSource, output_file
from bokeh.layouts import gridplot
import matplotlib as mpl
import numpy as np
import pandas as pd
import pyfasta
import argparse

# running Python3
if sys.version_info<(3,0,0):
  sys.stderr.write("You need python 3.0.0 or later to run this script\n")
  exit(1)

# check versions
print("allel {}".format(allel.__version__))
print("zarr {}".format(zarr.__version__))
print("numcodecs {}".format(numcodecs.__version__))
print("numpy {}".format(np.__version__))
print("pandas {}".format(pd.__version__))
print("msprime {}".format(msprime.__version__))
print("matplotlib {}".format(mpl.__version__))

# import custom modules
from allel_class import Chr



# options and args
parser = argparse.ArgumentParser()
parser.add_argument('-fa', "--fasta", help='path to fasta genome file masked')
parser.add_argument('-gff3', "--ggf3feature", help="path to gff3")
parser.add_argument('-zar', "--zarr", help="zarr file from vcf")
parser.add_argument('-g', "--group", help="chromosome to analyze")
parser.add_arguement('-m', "--meta", help="path to meta data")
parser.add_arguement('-s', "--samples", type=int, help="number samples")
parser.add_arguement('-a', "--altnum", type=int, default=1,
                     help="number alt alleles")
args = parser.parse_args()


if __name__ == "__main__":
    # from IPython.display import HTML
    callset = zarr.open_group(args.zarr, mode='r')
    # callset.tree(expand=True)
    gt_zarr = callset["{}/calldata/GT".format(args.group)]
    # gt_zarr.info
    pos = allel.SortedIndex(callset["{}/variants/POS".format(args.group)])
    gt = allel.GenotypeArray(gt_zarr)
    meta = pandas.read_csv(args.meta, sep='\t', usecols=['sample', 'pop', 'super_pop'])
    # meta.head()







def plot_fx(x, y, ax, title, y_label):
    ax.set_title(title)
    ax.set_ylabel(y_label)
    ax.plot(x, y)
    ax.margins(x=0, y=0)

fig, ax = plt.subplots()
plot_fx(x, y1, ax, "label")














