#! /usr/bin/env python

## --- bwt_query.py --- ##
##	Date: 22 Sep 2014
##	Updated: 12 Oct 2014
##	Purpose: Query a BWT at command-line, and dump hits in a nicely-aligned format (like csbio website does)

## 	TODO
##	 * also query on reverse complement
##	 * figure out why alignments aren't stacking in proper order

import os
import sys
import argparse
import glob
import collections

import MUSCython.MultiStringBWTCython as ms

from common import *
import snoop.util
import dna

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
parser.add_argument(	"queries", metavar = "K", type = str, nargs = "+",
			help = "queries; whitespace-separated if more than one" )
args = parser.parse_args()

bwt_dirs = glob.glob(args.msbwt)
msbwt = util.load_bwts(bwt_dirs)

for q in args.queries:
	for i in range(0, len(bwt_dirs)):
		print "Searching in: " + bwt_dirs[i]
		q = q.upper()
		(hits, dollars) = util.get_reads(msbwt[i], q)
		aln = util.pseudoalign(hits, dollars, args.spacer[0])
		if aln is not None:
			for j in range(0, len(aln)):
				print "".join(aln[j][:])
			varcodes = util.call_variant_sites(aln, args.maf, [ args.spacer[0], "N" ])
			print varcodes
		else:
			print "no hits found for query: " + q
		print ""
