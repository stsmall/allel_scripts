#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 13:54:20 2017

@author: scott
"""
import subprocess


def vcf2plink_fx(thindict, vcf):
    """Takes a thinned file from allel, thins the original vcf to match, then
    converts the vcf to plink. Returns plink formatted vcf to run in ADMIXTURE
    and Treemix
    """
    with open("thin.pos", 'w') as tpos:
        for nchr in thindict.keys():
            for pos in thindict[nchr]:
                tpos.write("{}\t{}\n".format(nchr, pos))
    command = "vcftools --vcf " + vcf + " --positions thin.pos --recode --out \
    thinnedsites"
    proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    proc.wait()
    # call plink
    command = "plink --make-bed thinnedsites.recode.vcf --allow-extra-chrom \
    --out thinnedplink"
    proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    proc.wait()
    return(None)
