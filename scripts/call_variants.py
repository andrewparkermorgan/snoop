#! /usr/bin/env python
"""
call_variants.py
Given pre-computed counts of k-mers from a reference sequence, use local assembly to identify variants.
"""

from __future__ import print_function

import os
import sys
import argparse
import numpy as np
import logging

import pybedtools as pbt

from snoop import util, io, assembly, swalign, rle
from snoop import kmers as ktools

parser = argparse.ArgumentParser(description = "Perform variant-calling by local de novo assembly, starting from reference k-mer counts.")
parser.add_argument(	"-M","--msbwt", nargs = "+", type = io.readable_dir,
 						required = True,
						help = "one or more BWTs with which to perform variant-calling" )
parser.add_argument(	"-w","--window", nargs = 2, type = int,
						default = [2,5],
						help = "minimum and maximum number of consecutive reference misses to trigger assembly [default: %(default)d]" )
parser.add_argument(	"-x","--coverage", type = int,
						required = False,
						help = "expected k-mer coverage of target sequence [default: estimated from data]" )
parser.add_argument(	"--counts", nargs = "?", type = argparse.FileType("r"),
						default = sys.stdin,
						help = "BED-like file of k-mer counts and coordinates [default: stdin]")
args = parser.parse_args()

## start error trace
logging.basicConfig(level = logging.DEBUG)
logging.StreamHandler(stream = sys.stderr)
logger = logging.getLogger()

## read in k-mer counts so I can be a bit lazy
## TODO: come up with a way to do this on a stream or in chunks
kmers = [ r for r in pbt.BedTool(args.counts) ]
counts = np.array([ int(r.score) for r in kmers])
if len(kmers) < 2:
	sys.exit("Not enough k-mer counts provided.")
ksize = len(kmers[0].name)
kstep = kmers[1].start - kmers[0].start
logger.info("Loaded counts for {} k-mers from file: <{}>".format(len(counts), args.counts.name))
logger.info("The k-mer size = {} and offset = {}.".format(ksize, kstep))
logger.info("Will attempt to assemble targets in size range [{}, {}].".format(args.window[0]*kstep, args.window[1]*kstep))

## load up BWTs
bwts = util.load_bwts(args.msbwt)
logger.info("Loaded {} msBWTs: {}".format(len(bwts), args.msbwt))

COVERAGE_QUANTILE = 30
if not args.coverage:
	coverage = np.percentile(counts, COVERAGE_QUANTILE)
else:
	coverage = args.coverage
min_coverage, max_coverage = int(0.3*coverage), int(3.0*coverage)
logger.info("Average coverage computed as {:02f} quantile from data: {}".format(float(COVERAGE_QUANTILE)/100, coverage))
logger.info("Will attempt to assemble targets with flanking coverage in range [{}, {}]".format(min_coverage, max_coverage))

## identify runs of zero-counts, which indicate a variant
zeros = rle.rle( np.array(counts == 0, dtype = np.int32) )
bounds = zeros.filter_lengths(lambda x: x >= min(args.window) and x <= max(args.window))

# write VCF header
SAMPLES = [ os.path.basename(x.rstrip("/")) for x in args.msbwt ]
print("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", *SAMPLES, sep = "\t")

## iterate on runs
for s,e in bounds:
	## we care about the zeros only
	if not zeros._original[s]:
		pass
	else:
		## can only work with runs having flanking information
		if s > 1 and e < (len(counts)-1):

			## check flanking coverage
			flank_left = int(kmers[s-1].score)
			flank_right = int(kmers[e].score)
			if min(flank_left, flank_right) < min_coverage or max(flank_left, flank_right) > max_coverage:
				continue

			## stitch together reference sequence
			ref = ""
			for r in kmers[(s-1):(e+1)]:
				#print(str(r).strip())
				ref += str(r.name)[:kstep]
			ref += str(kmers[e].name)[kstep:]

			## step through BWTs trying to assemble between flanking k-mers
			for bb, b in enumerate(bwts):
				alt, n = assembly.build_bridge(b, str(kmers[s-1].name), str(kmers[e].name))
				# check for successful assembly
				if not len(alt):
					# unsuccessful assembly
					logger.info("Couldn't bridge range {}:{}-{}".format(kmers[s-1].chrom, kmers[s-1].start+1, kmers[e].end))
				else:
					# assembly successful; align to reference
					ref_aln, alt_aln, ops = swalign.align(ref, alt[0])

					# check coverage of alt alleles
					alt_counts = []
					for k in ktools.kmerize(alt[0], ksize, kstep):
						if len(k) == ksize:
							alt_counts.append( util.count_reads(b, k) )
					support = np.median(alt_counts)
					#print(ref_aln)
					#print(ops)
					#print(alt_aln)
					# iterate on variant sites in alignment
					for offset, site_ref, site_alt in swalign.reconcile(ref_aln, alt_aln):
						# generate VCF entry
						CHROM = kmers[s-1].chrom
						POS = kmers[s-1].start+offset+1 #NB: must convert to 1-based
						ID = "."
						REF = site_ref
						ALT = site_alt
						QUAL = "."
						FILTER = "."
						INFO = "DPL={};DPR={};DP={}".format(flank_left, flank_right, support)
						FORMAT = "GT"
						print(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, "1/1", sep = "\t")
