import argparse
import allel
import numpy as np
import xpclr
import h5py
import sh
import pandas as pd


# FUNCTIONS
def load_hdf5_data(hdf5_fn, chrom, s1, s2):
    callset = h5py.File(hdf5_fn, mode='r')
    samples = callset['samples'][:]
    sample_name = [sid.decode() for sid in samples.tolist()]
    idx1 = np.array([sample_name.index(sid) for sid in s1])
    idx2 = np.array([sample_name.index(sid) for sid in s2])
    g = allel.GenotypeChunkedArray(callset["calldata/GT"])
    pos = allel.SortedIndex(callset["variants/POS"][:])
    return g.take(idx1, axis=1), g.take(idx2, axis=1), pos


def tabulate_results(chrom, model_li, null_li, selectionc, counts, windows,
                     edges):

    lidf = pd.DataFrame(np.vstack((model_li, null_li, selectionc, counts)).T,
                        columns=["modelL", "nullL", "sel_coef", "nSNPs"])

    # these are the nominal windows
    winf = pd.DataFrame(windows, columns=["start", "stop"])

    # these are the "real" windows. Gives a guide to how close we are.
    realf = pd.DataFrame(edges, columns=["pos_start", "pos_stop"])

    out = pd.concat([winf, realf, lidf], axis=1)

    out["xpclr"] = 2 * (out.modelL - out.nullL)
    out["xpclr_norm"] = (out.xpclr - np.nanmean(out.xpclr))/np.nanstd(out.xpclr)

    out.insert(0, "chrom", np.repeat(chrom, len(out)))

    string_id = ["{0}_{1:08d}_{2:08d}".format(r.chrom, r.pos_start, r.pos_stop)
                 for i, r in out.iterrows()]
    out.insert(0, "id", string_id)

    return out


# Argument parsing
psr = argparse.ArgumentParser(
    description='Tool to calculate XP-CLR as per Chen, Patterson, Reich 2010')

# files:
psr.add_argument('--out', "-O", required=True, help='output file')

# data:
psr.add_argument('--hdf5', "-I", required=False, help='input hdf5 filestem')

psr.add_argument('--samplesA', '-Sa', action='store', default=None,
                 dest='samplesA', required=False,
                 help='Which samples comprise population A. Comma separated')

psr.add_argument('--samplesB', '-Sb', action='store', default=None,
                 dest='samplesB', required=False,
                 help='Which samples comprise population B. Comma separated')

psr.add_argument('--chr', '-C', required=True, dest='chrom', action='store',
                 help='Which contig to use')

# parameters
psr.add_argument('--ld', '-L', dest='ldcutoff', action='store',
                 default=0.95, help='LD cutoff to apply for weighting')

psr.add_argument('--phased', '-P', dest='phased', action='store_true',
                 help='whether data is phased for more precise r2 calculation')

psr.add_argument('--maxsnps', '-M', dest='maxsnps', action='store',
                 default=200, help='max SNPs in a window')
psr.add_argument('--minsnps', '-N', dest='minsnps', action='store',
                 default=10, help='min SNPs in a window')

psr.add_argument('--rrate', '-R', dest='rrate', action='store',
                 default=1e-8, help='recombination rate per base')

psr.add_argument('--size', dest='size', action='store',
                 default=20000, help='window size', type=int)
psr.add_argument('--start', dest='start', action='store', type=int,
                 default=1, help='start base position for windows')
psr.add_argument('--stop', dest='stop', action='store', default=None, type=int,
                 help='stop base position for windows')
psr.add_argument('--step', dest='step', action='store', default=20000, type=int,
                 help='windows step for slide')

args = psr.parse_args()


fn = args.out
sh.touch(fn)    # so that we know at the start whether we have write access
print("XP_CLR/scikit-allel: {0}".format(xpclr.__version__))
chromosome = args.chrom
samples1 = [sample_id.strip() for sample_id in args.samplesA.split(",")]
samples2 = [sample_id.strip() for sample_id in args.samplesB.split(",")]

# if mode is "hdf5"
g1, g2, positions = load_hdf5_data(args.hdf5.strip(), chromosome,
                                   samples1, samples2)

ac1 = g1.count_alleles()
ac2 = g2.count_alleles()

print("TOTAL:     {0} SNPs are in the provided input files".format(g1.shape[0]))
try:
    multiallelic = (ac1[:, 2:].sum(1) > 0) | (ac2[:, 2:].sum(1) > 0)
    print("EXCLUDING: {0} SNPs as multiallelic ".format(multiallelic.sum()))
except IndexError:
    multiallelic = np.zeros(len(ac1), dtype=bool)
    print("EXCLUDING: 0 SNPs as multiallelic ")

# all missing in either
missing = (np.array(ac1.sum(1)) == 0) | (np.array(ac2.sum(1)) == 0)
print("EXCLUDING: {0} SNPs as missing in all samples in a population"
      .format(np.sum(missing & ~multiallelic)))

# drop if fixed in AC2,
fixed_p2 = np.array(ac2.is_non_segregating())
print("EXCLUDING: {0} SNPs as invariant in population 2"
      .format(np.sum(fixed_p2 & ~missing & ~multiallelic)))

# now compress all!
include = (~multiallelic & ~fixed_p2 & ~missing)
print("TOTAL:     {0} SNPs included in the analysis".format(include.sum()))

g1 = g1.compress(include, axis=0)
g2 = g2.compress(include, axis=0)
positions = np.compress(include, positions, axis=0)
print("..done dropping above SNPs from analysis.")

# check everything is as expected.
assert g1.shape[0] == g2.shape[0] == positions.shape[0]

# determine windows
if args.stop is None:
    args.stop = positions[-1]
spacing = np.arange(args.start, args.stop, args.step)
scan_windows = np.vstack([spacing, spacing - 1 + args.size]).T

# main function
modelL, nullL, selcoef, nsnps, snpedges = \
    xpclr.xpclr_scan(g1, g2, positions, scan_windows, geneticd=None,
                     ldcutoff=args.ldcutoff, phased=args.phased,
                     maxsnps=args.maxsnps, minsnps=args.minsnps,
                     rrate=args.rrate)

df = tabulate_results(chromosome, modelL, nullL, selcoef,
                      nsnps, scan_windows, snpedges)

df.to_csv(fn, sep="\t", index=False)