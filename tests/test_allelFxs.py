#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 15:51:27 2019
@author: Scott T. Small

?can github run tests and check coverage upon push?

unit testing for functions

# tests code against versions of python and modules
conda install -c conda-forge tox

# tests assertions beyond listed
conda install -c conda-forge hypothesis

# how much of application is tested by unit tests
conda install -c anaconda coverage

# mock values useful for multi-level testing where func is dependent on another
# MagicMock
conda install mock

# decorator with pytest to avoid loading heavy examples
@pytest.fixture(scope='module')
pytest.raises(RuntimeError)
pytest.warns(RuntimeWarning)
    warnings.warn("")

**Note**

"""
import filecmp
import shutil
import gzip
from applymask2vcf import mask2dict as mask2dict
from applymask2vcf import applymask as applymask

def test_mask2dict():
    """Test of masking function
    """
    mask_file = "applymask.test.in.txt.gz"
    mask_dict = mask2dict(mask_file)
    x = {"500":["KirFol1"],
         "505":["KirFol1"],
         "11931":["KirFol1"],
         "12897":["KirFol1"],
         "13218":["KirFol1"],
         "17642":["KirFol1"]}
    assert(x == mask_dict)


def test_applymask():
    """Test of applying mask
    """
    mask_dict = {"500":["KirFol1"],
         "505":["KirFol1"],
         "11931":["KirFol1"],
         "12897":["KirFol1"],
         "13218":["KirFol1"],
         "17642":["KirFol1"]}
    applymask("1sample.vcf.gz", mask_dict, "masked.vcf", False)
    with gzip.open('masked.vcf.gz', 'rb') as f_in:
        with open('masked.vcf', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    assert(filecmp.cmp("masked.vcf", "masked.test.vcf") == True)

if __name__ == "__main__":
    test_mask2dict()
    test_applymask()
