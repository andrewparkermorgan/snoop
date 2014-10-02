#! /usr/bin/env python

## --- bwt_query.py --- ##
##	Date: 22 Sep 2014
##	Updated: NA
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


def align_hits(query, bwt, spacer = "-"):

	try:
		(s, e) = bwt.findIndicesOfStr(query)
		reads = []
		dollars = []
		k = 0
		for i in range(s,e):
			r = bwt.recoverString(i)
			d = r.find("$")
			reads.append( r[(d+1):] + r[0:d] )
			dollars.append(d)
			k = max(k, len(r))

		stack = [ y for (x,y) in sorted(zip(dollars, range(0, e-s))) ]
		aln = [ [spacer]*(k + max(dollars)-min(dollars)) for i in range(0, e-s) ]
		for ii in range(0, len(stack)):
			aln[stack[ii]][:] = spacer*(dollars[stack[ii]]-min(dollars)) + reads[stack[ii]] + "-"*(k-dollars[stack[ii]]-2)

		return aln

	except Exception as e:
		print e
		return None


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

bwt_dirs = glob.glob(args.msbwt)
msbwt = []
for ff in bwt_dirs:
	if not readable_dir(ff):
		continue
	try:
		msbwt.append( ms.loadBWT(ff) )
	except Exception as e:
		print "Couldn't load BWT at <{}>".format(ff)
		print e

for q in args.queries:
	for i in range(0, len(bwt_dirs)):
		print "Searching in: " + bwt_dirs[i]
		q = q.upper()
		aln = align_hits(q, msbwt[i], args.spacer[0])
		if aln is not None:
			for j in range(0, len(aln)):
				print "".join(aln[j][:])
			varcodes = call_variant_sites(aln, args.maf, [ args.spacer[0], "N" ])
			print varcodes
		else:
			print "no hits found for query: " + q
		print ""
