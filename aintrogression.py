#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 18:58:35 2017 @author: stsmall

dxy: pairwise distance / bases; low values possible introgression
    dxy as number of sequence diff between any 2 sequences, x and y, in two
    taxa, X and Y (divided by the number of sites), then dxy is the average
    distance between all sequences in the two species. # above assumes no
    variation in neutral mutation rate, low mutation rate can be mistaken for
    recent introgression

dmin: min(dxy), requires haplotypes
    minimum distance among all pairing of haplotypes in the 2 species. Pvalue
    by coalescent with no migration or from other parts of the genome average.
    # above assumes no variation in neutral mutation rate, low mutation rate
    can be mistaken for recent introgression

dout: (dXO + dYO) / 2
    average distance between species X and the Outgroup and Species Y and the
    Outgroup.

RND (relative node depth): dxy / dout
     Robust to low mutation rates like HKY test if neutrality. # not sensitive
     to low-frequency migrants.

Gmin: dmin / dxy
    low power, useful when migration prob and relative migration time is low
    This is likely due to the fact that as migrant lineages rise in frequency,
    dXY also gets lower. As a migrant haplotype approaches fixation, the ratio
    of dmin to dXY approaches 1

RNDmin: dmin / dout
    Similarly, like both dmin and Gmin, RNDmin should be sensitive to even rare
    migrant haplotypes. In addition, we expect RNDmin to be powerful even when
    migrants are high in frequency

## determine significance by coalescent simulation w/ and w/out migration,
    also use genome average, X, in anopheles assuming that region is not
    introgressed.

1) call all SNPs on a single genome
2) call SNPs on specific reference genome
    a) carry over all coordinates for a new vcf
        i) Maffilter --vcf
        ii) mvftools  # not a traditional vcf so some issues with analysis
        iii) Emrich script
"""
allel.windowed_df
pairwise diff
Dxy
Dxy_min
Dxy_mean
pariwiseDistance
pairwiseDistancerankamongsamples
Dxy_vector
Snn