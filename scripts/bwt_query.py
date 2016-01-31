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
parser.add_argument(	"-M", "--msbwt", type = io.readable_dir, nargs = "+",
			required = True,
			help= "path to one or more directories containing msBWT components" )
parser.add_argument(	"-s", "--spacer", type = str,
			default = "?",
			help = "spacer character for printed alignments, needs to be valid gap character for alignment software [default:%(default)s]" )
parser.add_argument(	"-f", "--maf", type = int,
			default = 1,
			help = "minimum 'minor allele frequency' (an integer number of reads) to accept site as variant [default:%(default)d]" )
parser.add_argument(	"-x", "--coverage", type = int,
			default = 5,
			help = "minimum number of reads covering a site to declare it callable [default:%(default)d]" )
parser.add_argument(	"-r", "--revcomp", action = "store_true",
			default = False,
			help = "also search for this query's reverse complement [default: False]" )
parser.add_argument(	"queries", metavar = "K", type = str, nargs = "+",
			help = "queries; whitespace-separated if more than one" )
args = parser.parse_args()

msbwt = util.load_bwts(args.msbwt)

def match_line(x):
	out = ""
	for i in x:
		if i:
			out += "*"
		else:
			out += " "
	return out

for q in args.queries:
	for i in range(0, len(args.msbwt)):
		print "Searching in: " + args.msbwt[i]
		q = q.upper()
		reads = util.get_reads(msbwt[i], q, args.revcomp)
		#reads.pseudoalign(k = len(q), spacer = args.spacer[0])
		reads.align()
		if reads.alignment is not None:
			rez = reads.find_variant_sites(args.maf, args.coverage, args.spacer)
			hap_count = reads.extract_haplotypes(rez[1])
			for r in reads.alignment.iterSeqs():
				print r
			print match_line(rez[0])
			print ""
			print "Counting variant positions with min coverage = {} and min MAF = {}.".format(args.coverage, args.maf)
			print "Number of variant positions: {} (of {} callable positions)".format(len(rez[1]), rez[2] - len(q))
			print hap_count
			#for j in range(0, len(reads.alignment)):
			#	print "".join(reads.alignment[j][:])
			#varcodes = reads.call_variant_sites(args.maf)
			# print varcodes
			# print "Simple haplotype count:", reads.simply_count_haplotypes(args.maf)
			#print "Consistency score:", reads.consistency_score(args.maf)
			#print "Priors on haplotype frequencies:", zip(args.alpha, args.beta)
			#print "Regularized number of haplotypes:", reads.count_haplotypes(args.alpha, args.beta)
			#print ""
			#print "Observed haplotypes:", reads.extract_haplotypes()[0]
			#print varcodes
			#for (h, n) in reads.haplotypes.iteritems():
			#	print str(h) + " " + str(n)

		else:
			print "no hits found for query: " + q
		print ""
