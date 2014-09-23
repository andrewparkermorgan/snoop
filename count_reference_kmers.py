#! /nas02/apps/python-2.7.1/bin/python
## /usr/bin/env python

## --- count_reference_kmers.py --- ##
##	Date: 17 June 2014
##	Udpated: 8 July 2014 (allow parallel search of multiple BWTs via globbing on -M option)
##	Purpose: given a list of genomic intervals (in bed format) and a reference sequence, generate k-mers and count occurrences of each in a msBWT

import os
import sys
import argparse
import subprocess
import tempfile
import glob
import numpy

import dna
from common import *

import pybedtools
import pyfasta
from MUSCython.MultiStringBWTCython import * # msBWT

parser = argparse.ArgumentParser(description = "Utility for generating and counting occurences of k-mers from a reference sequence in a msBWT.")
parser.add_argument(	"-r","--regions", type = argparse.FileType("rU"),
			required = True,
			help = "list of regions (in bed format) for which to make k-mers" )
parser.add_argument(	"-g","--genome", type = readable_file,
			default = "$GENOMES/mm10/mm10.fa",
			help = "reference genome fasta file from which to obtain k-mer sequences" )
parser.add_argument(	"-k","--kmer", type = int,
			default = 30,
			help = "k-mer size in bp (ie. length of probe)" )
parser.add_argument(	"-s","--step", type = int,
			help = "step size for making k-mers; if kmer_1 = (1, k), kmer_2 = (1+s, k+s)" )
parser.add_argument(	"-M","--msbwt",
			default = "./",
			help = "globbing expression to find directories containing msBWT components" )
parser.add_argument(	"-S","--summarize", action = "store_true",
			default = False,
			help = "aggregate counts over the intervals in 'regions' file, instead of reporting k-mer-level counts" )
parser.add_argument(	"-x","--minhits", type = int,
			default = 3,
			help = "in summary mode, trim k-mers with (strictly) fewer than this many hits" )
parser.add_argument(	"-X","--maxhits", type = int,
			default = 10000,
			help = "in summary mode, trim k-mers with (strictly) greater than this many hits" )
parser.add_argument(	"-v","--verbose", action = "store_true",
			default = False,
			help = "print a status trace to stderr" )
args = parser.parse_args()


## connect to msBWT(s)
msbwt = []
bwt_dirs = glob.glob(args.msbwt)
for ff in bwt_dirs:
	
	if not readable_dir(ff):
		continue
	try:
		msbwt.append( loadBWT(ff) )
	except Exception as e:
		print "Couldn't load BWT at <{}>".format(ff)
		print e

## fail if couldn't connect to at least one msBWT
if len(msbwt) == 0:
	sys.exit("No valid BWTs were provided; quitting.")
else:
	sys.stderr.write("Using the following BWTs:\n{}\n".format(str(bwt_dirs)))
	if args.summarize:
		sys.stderr.write("Working in summary mode.\n")

## get sequences of k-mers via pyfasta
fasta = pyfasta.Fasta(args.genome, key_fn = lambda key: key.split()[0])

## load regions
regions = pybedtools.BedTool(args.regions)
if args.step is None:
	args.step = args.kmer

if args.verbose:
	sys.stderr.write("Starting query loop.\n")

## loop on regions
for r in regions:

	## print status
	if args.verbose:
		sys.stderr.write( "\t... {}: {} - {}\n".format(r.chrom, r.start, r.end) )

	## write this region to tempfile to allow 'bedtools makewindows' to find it
	(fh, fn) = tempfile.mkstemp(prefix = ".bed")
	with os.fdopen(fh, "w") as rfile:
		rfile.write(str(r))

	## retrieve sequence
	## for perverse reasons, the following will not work: bedtools makewindows requires the -w argument to precede -s
	# kmers = regions.window_maker(regions, w = args.kmer, s = args.step)

	## ... so do this ugly hack instead
	p = subprocess.Popen(	["windowMaker", "-b", fn, "-w", str(args.kmer), "-s", str(args.step), "-i", "srcwinnum"],
				stdout = subprocess.PIPE,
				stderr = subprocess.PIPE,
				bufsize = 1 )

	# NB: since this this next BedTool is reading on a pipe, it can be consumed ONCE and only once
	kmers = pybedtools.BedTool(p.stdout)
	## ... end hack

	## initialize counts for this region
	width = r.end - r.start
	agg = []
	zeros = 0

	for kmer in kmers:
		
		## skip kmers which are clipped by region boundaries; this produces boundary artefact but I will ignore it
		if (kmer.end - kmer.start) != args.kmer:
			continue
		# try:
		seq = dna.unmask( fasta[kmer.chrom][ kmer.start:kmer.end ] )
		## loop over all open BWTs and query them with this sequence
		# initialize counters for this kmer
		count_f = 0
		count_r = 0
		for i in range(0, len(msbwt)):
			count_f += msbwt[i].countOccurrencesOfSeq( str(seq) )
			count_r += msbwt[i].countOccurrencesOfSeq( str(dna.revcomp(seq)) )
			
		agg.append( count_f + count_r )

		## record this kmer as zero iff it returned zero across all BWTs
		if (count_f + count_r) == 0:
			zeros += 1

		## fine-grained mode: print results for each k-mer
		if not args.summarize:
			print kmer.chrom, kmer.start, kmer.end, seq, count_f, count_r, (count_f + count_r)

		# except Exception as e:
		# 	# print "Sequence was: {}; searched failed.".format(seq)
		# 	print type(e)
		# 	pass

	## coarse-grained mode: print only the summary over the region 
	if args.summarize:

		med, q, qqq, mu = ("NA", "NA", "NA", "NA")
		agg_trimmed = filter(lambda x: (x > args.minhits) and (x <= args.maxhits), agg)
		
		if len(agg_trimmed):
			med = numpy.median(agg_trimmed)
			q = numpy.percentile(agg_trimmed, 25.0)
			qqq = numpy.percentile(agg_trimmed, 75.0)
			mu = numpy.mean(agg_trimmed)

		print r.chrom, r.start, r.end, zeros, med, mu, q, qqq

	## clean up temp file
	try:
		os.remove(fn)
	except:
		pass

if args.verbose:
	sys.stderr.write("Done.\n")
