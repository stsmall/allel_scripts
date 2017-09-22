#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 16:56:47 2017

@author: scott
"""


def Fstatistics_fx(ac_subpopscat, subpops):
    """
    """
    poplist = ["Haiti", "Mali", "Kenya", "PNG"]
    pop1 = poplist[0]
    pop2 = poplist[1]
    pop3 = poplist[2]
    pop4 = poplist[3]
    pop1_idx = subpops[pop1]
    pop2_idx = subpops[pop2]
    pop3_idx = subpops[pop3]
    pop4_idx = subpops[pop4]
    acu = allel.AlleleCountsArray(ac_subpopscat[pop1][:] +
                                  ac_subpopscat[pop2][:] +
                                  ac_subpopscat[pop3][:] +
                                  ac_subpopscat[pop4][:])
    dip = (len(pop1_idx) + len(pop2_idx) + len(pop3_idx) + len(pop4_idx)) * 2.0
    # prc missing
    prct = .80
    miss = (acu[:, 0] + acu[:, 1]) > (dip * prct)
    seg = acu.is_segregating()
    flt = miss * seg
    blen = int(round(sum(flt)/5))
    aca = allel.AlleleCountsArray(ac_subpopscat[pop1].compress(flt,
                                  axis=0)[:, :2])
    acb = allel.AlleleCountsArray(ac_subpopscat[pop2].compress(flt,
                                  axis=0)[:, :2])
    acc = allel.AlleleCountsArray(ac_subpopscat[pop3].compress(flt,
                                  axis=0)[:, :2])
    acd = allel.AlleleCountsArray(ac_subpopscat[pop4].compress(flt,
                                  axis=0)[:, :2])
    # acc = test, aca = pop1, acb = pop2, acd = outgroup
    f2 = allel.stats.admixture.patterson_f2(ac1, ac2)
    f2_blocks = allel.stats.moving_statistic(f2, statistic=np.mean, size=blenW)
    f2_m, f2_se, _, _ = jackknife(f2_blocks)
            # f2
        f2 = allel.stats.admixture.patterson_f2(ac1c, ac2c)
        f2_blocks = allel.stats.moving_statistic(f2, statistic=np.mean,
                                                 size=blenWc)
        f2_m, f2_se, _, _ = jackknife(f2_blocks)
        print("{}, {}, f2: {:.4f} +/- {:.4f}".format(nchr, pair, f2_m, f2_se))
    print('{}-{} : {:.4f} +/- {:.4f} (f2)'.format(pop1, pop2, f2_m, f2_se))
    # Patternson F3 admixture
    f3, f3_se, f3zscore, _, _ = allel.stats.blockwise_patterson_f3(aca, acb,
                                                                   acc, blen)
    print("{}: {}, {}, {}".format([pop1, pop2, pop3], f3, f3_se, f3zscore))
    f3, f3_se, f3zscore, _, _ = allel.stats.blockwise_patterson_f3(acb, acc,
                                                                   aca, blen)
    print("{}: {}, {}, {}".format([pop2, pop3, pop1], f3, f3_se, f3zscore))
    f3, f3_se, f3zscore, _, _ = allel.stats.blockwise_patterson_f3(acc, aca,
                                                                   acb, blen)
    print("{}: {}, {}, {}".format([pop3, pop1, pop2], f3, f3_se, f3zscore))
    # Patterson F4
    # Mali is admixted of Kenya
    f4, f4_se, f4zscore, _, _ = allel.stats.blockwise_patterson_d(aca, acb,
                                                                  acc, acd,
                                                                  blen)
    print("{}: {}, {}, {}".format([pop1, pop2, pop3, pop4], f4, f4_se,
          f4zscore))
    # Kenya is admixed of Mali
    f4, f4_se, f4zscore, _, _ = allel.stats.blockwise_patterson_d(aca, acc,
                                                                  acb, acd,
                                                                  blen)
    print("{}: {}, {}, {}".format([pop1, pop3, pop2, pop4], f4, f4_se,
          f4zscore))
    f4, f4_se, f4zscore, _, _ = allel.stats.blockwise_patterson_d(acb, acc,
                                                                  aca, acd,
                                                                  blen)
    print("{}: {}, {}, {}".format([pop2, pop3, pop1, pop4], f4, f4_se,
          f4zscore))

    return(None)
