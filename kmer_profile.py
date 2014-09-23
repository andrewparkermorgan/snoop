#! /usr/bin/env python

## --- kmer_profile.py --- ##
##      Date: 13 August 2014
##      Udpated: NA
##      Purpose: given a list of msBWTs, generate a (possibly-downsampled) k-mer kmer_profile

import os
import sys
import argparse
import numpy as np

import dna
from common import *

from MUSCython.MultiStringBWTCython import * # msBWT

## Matt's generator for the union of k-mers present in a list of msBWT objects
def kmerIterator(bwtList, k, threshold):
	'''
	@param bwtList - the list of BWTs (the actual objects) to search
	@param k - the length of the k-mers
	@param threshold - a k-mer must have more than the threshold to be returned
	@return - a generator for the k-mers
	'''
	vcLen = 6
	currentK = 0
	kmerStarts = np.zeros(dtype = "<u8", shape = (k+1, len(bwtList)))
	kmerEnds = np.zeros(dtype = "<u8", shape = (k+1, len(bwtList)))
	for x in range(0, len(bwtList)):
		kmerEnds[0, x] = bwtList[x].getTotalSize()
	
	status = np.zeros(dtype = "<u1", shape = (k+1, ))
	
	while currentK > 0 or status[0] < vcLen:
		if currentK == k:
			for x in range(0, len(bwtList)):
				if kmerStarts[currentK, x] != kmerEnds[currentK, x]:
					seq = bwtList[x].recoverString(kmerStarts[currentK, x])[0:k]
			yield (seq, kmerEnds[currentK]-kmerStarts[currentK])
			currentK -= 1
		elif status[currentK] == 0 or status[currentK] == 4:
			status[currentK] += 1
		elif status[currentK] == 6:
			currentK -= 1
		else:
			# we need to write the correct value into the kmerStarts and kmerEnds at currentK+1
			foundVals = False
			for x in range(0, len(bwtList)):
				fmStart = bwtList[x].getFullFMAtIndex(kmerStarts[currentK, x])
				fmEnd = bwtList[x].getFullFMAtIndex(kmerEnds[currentK, x])
				kmerStarts[currentK+1, x] = fmStart[status[currentK]]
				kmerEnds[currentK+1, x] = fmEnd[status[currentK]]
				
				if fmEnd[status[currentK]] - fmStart[status[currentK]] > threshold:
					foundVals = True
			
			# update status[currentK]
			status[currentK] += 1
			
			if foundVals:
				# increment currentK if we found something
				currentK += 1
				status[currentK] = 0

def profile_kmers(msbwt_paths, k, threshold = 1, sampling_rate = 1.0, n = 1, N = float("inf"), unique_only = False, out_stream = None, delim = "\t"):

	## initialize list of msBWTs and their sizes
	msbwt = []
	sizes = np.zeros(len(msbwt_paths), dtype = "u8")
	weights = np.zeros(len(msbwt_paths), dtype = "f8")

	for i in range(0, len(msbwt_paths)):
		msbwt.append( loadBWT(msbwt_paths[i]) )
		sizes[i] = msbwt[i].getTotalSize()
	weights = sizes/float(sum(sizes))

	total_jsd = np.zeros(n)
	j = 0

	## iterate using the k-mer generator
	for (kmer, counts) in kmerIterator(msbwt, k, threshold):

		## increment counter and exit early, if requested
		j += 1
		if j > N:
			break

		## extract string-formatted counts from numpy array
		str_counts = [ str(counts[i]) for i in range(0, len(msbwt)) ]
		is_unique = any(counts < threshold)

		if not args.unique:

			flag = np.ones(n)
			if not (sampling_rate == 1.0):
				flag = np.array(np.random.ranf(n) < sampling_rate, dtype = int)
			
			## compute element-wise Jensen-Shannon divergence (JSD)
			probs = [ (float(max(counts[i], threshold))/sizes[i]) for i in range(0, len(msbwt)) ]
			mean_prob = np.average(probs, weights = weights)
			logged = [ np.log(probs[i]/mean_prob)*probs[i] for i in range(0, len(msbwt)) ]
			jsd = np.average(logged, weights = weights)
			total_jsd += jsd * flag

			## write to output stream, if provided
			if out_stream is not None:
				out_stream.write( delim.join([ kmer, str_counts, str(jsd) ]) + "\n" )

		else:
			if is_unique:
				## also check for uniqueness of reverse-complement, since that's probably more biologically-meaningful
				rev_counts = np.array([ msbwt[i].countOccurrencesOfSeq(str(dna.revcomp(kmer))) for i in range(0, len(msbwt)) ])
				total_counts = counts + rev_counts
				if any(total_counts < threshold):
					out_stream.write( delim.join( [kmer].extend([ str(int(x)) for x in total_counts ]) ) )

	if args.unique:
		return 0
	else:
		return total_jsd

## parse command-line arguments
parser = argparse.ArgumentParser(description = "Utility for generating a k-mer profile from a set of msBWTs.")
parser.add_argument(	"-M","--msbwt", nargs = "*", type = readable_dir,
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
			help = "number of replicates to use for computign standard errors of divergence metric [EXPERIMETNAL] [default: %(default)f]" )
parser.add_argument(	"-N","--maxk", type = float,
			default = float("Inf"),
			help = "max number of k-mers to count; useful only for testing [default: %(default)f]" )
parser.add_argument(	"-s","--seed", type = int,
			default = None,
			help = "seed for RNG which controls downsampling; set it manually for reproducible results [default: %(default)s]" )
parser.add_argument(	"-v","--verbose", action = "store_true",
			default = False,
			help = "print a status trace to stderr" )
args = parser.parse_args()


print "Using the following msBWTs:"
print args.msbwt

## manually set seed for RNG if required
if args.seed is not None:
	np.random.seed(args.seed)

jsd = profile_kmers(args.msbwt, args.kmer, args.minhits, args.sampling_rate, args.niter, args.maxk, args.unique, None)
print "Total Jensen-Shannon divergence:", np.mean(jsd), np.std(jsd)
# np.savetxt(sys.stderr, jsd)
