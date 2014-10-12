#! /usr/bin/env python

## --- snoop/util.py --- ##
##	Date: 12 Oct 2014
##	Updated: NA
##	Purpose: Utility functions for retrieving and manipulating reads from msBWTs

import os
import sys
import numpy as np

import MUSCython.MultiStringBWTCython as ms

from common import *

## wrapper to load BWTs residing on a given list of paths, and return them as a list of objects
## returns None if no BWTs were loaded successfully
def load_bwts(bwt_dirs):

	msbwt = []
	for ff in bwt_dirs:
		if not readable_dir(ff):
			continue
		try:
			msbwt.append( ms.loadBWT(ff) )
		except Exception as e:
			sys.stderr.write("Couldn't load BWT at <{}>\n".format(ff))
			sys.stderr.write(str(e))

	if len(msbwt):
		return msbwt
	else:
		return None

## resurrect original read from the BWTed version: eg 'AC$GT' > 'GTAC'
## returns tuple of (resurrected read, dollar-sign index)
def unrotate_read(read):

	read = str(read)
	dollar = read.find("$")
	return ( read[ (dollar+1): ] + read[ 0:dollar ], dollar )

## given a bwt and a query string, return the (resurrected) hits
## returns tuple of (resurrected read, dollar-sign index)
def get_reads(bwt, query):

	reads = []
	dollars = []

	try:
		(s, e) = bwt.findIndicesOfStr(str(query).upper())
		for i in range(s,e):
			(r, d) = unrotate_read(bwt.recoverString(i))
			reads.append(r)
			dollars.append(d)
	except Exception as e:

		return (reads, dollars)

		print e
		return None

## given lists of (resurrected) reads and dollar-indices, 'pseudoalign' them by centering them on the original query
def pseudoalign(reads, dollars, spacer = "-"):

	nreads = len(reads)
	maxlen = max([ len(r) for r in reads ])

	stack = [ y for (x,y) in sorted(zip(dollars, range(0, nreads))) ]
	aln = [ [spacer]*(maxlen + max(dollars)-min(dollars)) for i in range(0, nreads) ]
	for ii in range(0, len(stack)):
		aln[stack[ii]][:] = spacer*(dollars[stack[ii]]-min(dollars)) + reads[stack[ii]] + "-"*(k-dollars[stack[ii]]-2)

	return aln

## step through pseudoalignment of reads and find sites which have variant alleles above some specified frequency
def call_variant_sites(aln, maf, ignore = ["-","N"]):

	nseq = len(aln)
	if nseq:
		width = len(aln[0])
		rez = [" "]*width
		maf_obs = []
		for i in range(0, width):
			seq = [ aln[j][i] for j in range(0, nseq) ]
			trimmed = [ x for x in seq if x not in ignore ]
			counts = collections.Counter(trimmed)
			total = float(sum(counts.values()))
			if len(trimmed):
				(nuc, maj) = counts.most_common()[0]
				maf_obs.append(1-float(maj)/total)
				if ((1.0 - maj/total) > maf):
					rez[i] = "*"
		# print maf_obs
		return "".join(rez)
	else:
		return None

## compute some kind of 'consistency score' which approximates the number of haplotypes represented in a pseudoalignment
## TODO: come up with some theoretical justification for this
def consistency_score(aln, maf, alphabet = "ACGT", eps = 0.00001):

	nseq = len(aln)
	if nseq:

		width = len(aln[0])
		freq = np.empty(shape = (width, len(alphabet)), type = float)
		freq.fill(eps)
		weights = np.zeros(shape = width, type = int)

		for i in range(0, width):

			seq = [ aln[j][i] for j in range(0, nseq) ]
			trimmed = [ x for x in seq if x in alphabet ]
			counts = collections.Counter(trimmed)
			total = sum(counts.values())

			k = 0
			for (nuc, freq) in counts.most_common():
				prob = float(freq/total)
				if prob > maf:
					freq[i][k] = prob
					k += 1
			weights[i] = total

		normfreq = freq/np.linalg.norm(freq, 1, 1).reshape(freq.shape[0], 1) ## column-wise frequencies
		colentropy = -1 * np.sum(normfreq * np.log2(normfreq), 1) ## column-wise Shannon information
		nalleles = np.power(len(alphabet), colentropy) ## 'expected number of alleles', column-wise
		adjweights = weights/np.sum(weights) * (nalleles - 1) ## adjusted weights: columns with less info count less; columns w/1 allle are downweighted

		return np.average(nalleles, weights = adjweights)

	else:
		return None
