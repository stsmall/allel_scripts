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

F_D (Martin 2015)

R_D (Racimo 2016)

U20 (Racimo 2016, Jagoda 2018)

Q95 (Racimo 2016, Jagoda 2018)

Dp_intro (Racimo 2016)

Dp_combo (Racimo 2016)

r2_intro (Racimo 2016)

r2_combo (Racimo 2016)

"""

Dxy
Dxy_min
Dxy_mean
Dxy_vector
