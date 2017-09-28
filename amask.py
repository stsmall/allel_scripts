#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 16:28:17 2017

@author: scott
"""

import numpy as np


def access_mask(chrmaskdict, nchr):
    """
    """
    chrlen = chrmaskdict[nchr]
    m = np.ones(chrlen, dtype=bool)
    if len(chrmaskdict[nchr]) != 1:
        for mask in range(1, len(chrmaskdict[nchr]) + 1):
            start, stop = mask
            m[start:stop] = False
    return(m)
