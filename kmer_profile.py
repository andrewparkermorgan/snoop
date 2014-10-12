#! /usr/bin/env python

## --- kmer_profile.py --- ##
##      Date: 13 August 2014
##      Udpated: 29 September 2014
##      Purpose: given a list of msBWTs, generate a (possibly-downsampled) k-mer kmer_profile

import os
import sys
import argparse
import numpy as np

import MUSCython.MultiStringBWTCython as ms

from common import *
import snoop.util
import snoop.kmers
import dna

## parse command-line arguments
parser = argparse.ArgumentParser(description = "Utility for generating a k-mer profile from a set of msBWTs.")
parser.add_argument(	"-M","--msbwt", nargs = "+", type = readable_dir,
			help = "one or more msBWT paths" )
parser.add_argument(	"-k","--kmer", type = int,
			default = 30,
			help = "value of k for which to enumerate k-mers [default: %(default)s]" )
parser.add_argument(	"-x","--minhits", type = int,
			default = 1,
			help = "only emit k-mers with (strictly) greater than this many occurrences in at least one msBWT [default: %(default)s]" )
parser.add_argument(	"-u","--unique", action = "store_true",
			help = "return only k-mers which are unique to at least one dataset (subject to abundance threshold)" )
parser.add_argument(	"-f","--sampling_rate", type = float,
			default = 1.0,
			help = "fraction by which to downsample k-mer counts [EXPERIMETNAL] [default: %(default)f]" )
parser.add_argument(	"-n","--niter", type = int,
			default = 100,
			help = "number of replicates to use for computing standard errors of divergence metric [EXPERIMETNAL] [default: %(default)f]" )
parser.add_argument(	"-N","--maxk", type = float,
			default = float("Inf"),
			help = "max number of k-mers to count; useful only for testing [default: %(default)f]" )
parser.add_argument(	"-s","--seed", type = int,
			default = None,
			help = "seed for RNG which controls downsampling; set it manually for reproducible results [default: %(default)s]" )
parser.add_argument(	"-o", "--output", type = writeable_or_stdout_handle,
			default = None,
			help = "file to which to write the (possibly-downsampled, abundance-thresholded) k-mer profile (use '-' for stdout) [default: %(default)s]" )
parser.add_argument(	"-v","--verbose", action = "store_true",
			default = False,
			help = "print status trace to stderr" )
args = parser.parse_args()


print "Using the following msBWTs:"
print args.msbwt

msbwt = util.load_bwts(args.msbwt)

## manually set seed for RNG if required
if args.seed is not None:
	np.random.seed(args.seed)

profiler = kmers.KmerProfiler(msbwt, args.kmer, args.minhits, args.sampling_rate, args.maxk)

total_jsd = 0
for (kmer, counts, jsd, cum_jsd) in profiler.profile():
	outfields = [kmer]
	outfields.extend(counts)
	outfields.append(jsd)
	total_jsd = cum_jsd

	args.output.write( "\t".join(outfields) + "\n" )

print "Total Jensen-Shannon divergence:", total_jsd
# np.savetxt(sys.stderr, jsd)
