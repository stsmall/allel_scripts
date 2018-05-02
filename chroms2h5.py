#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 16:08:00 2018

@author: scott
"""
import allel
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-v', "--vcfFile", type=str, required=True,
                    help="vcf file")
parser.add_argument('-c', "--chromlist", nargs='+', required=True,
                    help="list of chromosomes")
args = parser.parse_args()


def makeh5fromvcf(vcfin, chromlist):
    """
    """
    for c in chromlist:
        h5out = "{}{}.h5".format(vcfin.rstrip("vcf"), c)
        fieldsfromvcf = ['samples', 'calldata/GQ', 'variants/ALT',
                         'variants/REF', 'variants/QUAL', 'variants/CHROM',
                         'variants/POS', 'variants/AF', 'variants/AB',
                         'variants/MQM', 'variants/DP', 'calldata/DP',
                         'calldata/AD', 'calldata/GT']
        allel.vcf_to_hdf5(vcfin, h5out, fields=fieldsfromvcf,
                          types={'calldata/GQ': 'float32'}, region=c)
    return(None)


if __name__ == "__main__":
    makeh5fromvcf(args.vcfFile, args.chromlist)
