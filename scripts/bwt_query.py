#! /usr/bin/env python

## --- snoop/bwt_query.py --- ##
##	Date: 22 Sep 2014
##	Updated: 14 Oct 2014
##	Purpose: Query a BWT at command-line; dump hits in a nicely-aligned format (like csbio website does) and estimate number of haplotypes present

## 	TODO
##	 * also query on reverse complement
##	 * figure out why alignments aren't stacking in proper order

import os
import sys
import argparse
import glob

import MUSCython.MultiStringBWTCython as ms

from snoop import util, io, dna, haplotype

parser = argparse.ArgumentParser(description = "Simple utility for dumping reads from msBWT in an alignment-like format")
parser.add_argument(	"-M", "--msbwt", type = str,
			required = True,
			help= "globbing expression to find directories containing msBWT components" )
parser.add_argument(	"-s", "--spacer", type = str,
			default = "-",
			help = "spacer character for printed alignments [default:%(default)s]" )
parser.add_argument(	"-f", "--maf", type = float,
			default = 0.1,
			help = "minimum 'minor allele frequency' to accept site as variant [default:%(default)f]" )
parser.add_argument(	"-a", "--alpha", type = float, nargs = "+",
			default = [1.0],
			help = "hyperparameter for prior on allele counts: alpha [default:%(default)f]" )
parser.add_argument(	"-b", "--beta", type = float, nargs = "+",
			default = [1.0],
			help = "hyperparameter for prior on allele counts: beta [default:%(default)f]" )
parser.add_argument(	"-r", "--revcomp", action = "store_true",
			default = False,
			help = "also search for this query's reverse complement [default: False]" )
parser.add_argument(	"queries", metavar = "K", type = str, nargs = "+",
			help = "queries; whitespace-separated if more than one" )
args = parser.parse_args()

bwt_dirs = glob.glob(args.msbwt)
msbwt = util.load_bwts(bwt_dirs)

for q in args.queries:
	for i in range(0, len(bwt_dirs)):
		print "Searching in: " + bwt_dirs[i]
		q = q.upper()
		reads = util.get_reads(msbwt[i], q, args.revcomp)
		reads.pseudoalign(k = len(q), spacer = args.spacer[0])
		if reads.alignment is not None:
			for j in range(0, len(reads.alignment)):
				print "".join(reads.alignment[j][:])
			varcodes = reads.call_variant_sites(args.maf)
			print varcodes
			print "Simple haplotype count:", reads.simply_count_haplotypes(args.maf)
			print "Consistency score:", reads.consistency_score(args.maf)
			print "Priors on haplotype frequencies:", zip(args.alpha, args.beta)
			print "Regularized number of haplotypes:", reads.count_haplotypes(args.alpha, args.beta)
			print ""
			print "Observed haplotypes:", reads.extract_haplotypes()[0]
			print varcodes
			for (h, n) in reads.haplotypes.iteritems():
				print str(h) + " " + str(n)

		else:
			print "no hits found for query: " + q
		print ""
