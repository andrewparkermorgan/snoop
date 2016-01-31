#! /nas02/apps/python-2.7.1/bin/python
## /usr/bin/env python

## --- snoop/count_reference_kmers.py --- ##
##	Date: 17 June 2014
##	Udpated: 27 Dec 2014 (allow bootstrap resampling of k-mer counts for getting CIs on sequence divergence)
##	Purpose: given a list of genomic intervals (in bed format) and a reference sequence, generate k-mers and count occurrences of each in a msBWT

import os
import sys
import argparse
import subprocess
import tempfile
import glob
import numpy as np

from snoop import util, io, dna

import pybedtools
import pyfasta
from MUSCython.MultiStringBWTCython import * # msBWT

parser = argparse.ArgumentParser(description = "Utility for generating and counting occurences of k-mers from a reference sequence in a msBWT.")
parser.add_argument(	"-r","--regions", nargs = "?", type = argparse.FileType("rU"),
			default = sys.stdin,
			help = "list of regions (in bed format) for which to make k-mers" )
parser.add_argument(	"-g","--genome", type = io.readable_file,
			default = "$GENOMES/mm10/mm10.fa",
			help = "reference genome fasta file from which to obtain k-mer sequences [default: %(default)s]" )
parser.add_argument(	"-k","--kmer", type = int,
			default = 30,
			help = "k-mer size in bp (ie. length of probe) [default: %(default)s]" )
parser.add_argument(	"-s","--step", type = int,
			help = "step size for making k-mers; if kmer_1 = (1, k), kmer_2 = (1+s, k+s) [default: same as k-mer size]" )
parser.add_argument(	"-M","--msbwt", nargs = "+",
			default = "./",
			help = "one or more directories containing msBWT components" )
parser.add_argument(	"-n","--normalize", type = io.readable_dir,
			default = None,
			help = "msBWT against which to normalize abundances of each queried k-mer" )
parser.add_argument(	"-R","--revcomp", action = "store_true",
			default = None,
			help = "normalize against both k-mer and its reverse complement (if reference msBWT is not from an assembled genome)" )
parser.add_argument(	"-a", "--alleles", action = "store_true",
			default = False,
			help = "compute 'allele scores': try to guess number of alleles from reads" )
parser.add_argument(	"--alpha", type = float, nargs = "+",
			default = 2.0,
			help = "hyperparameter for prior on allele counts: alpha [default:%(default)f]" )
parser.add_argument(	"--beta", type = float, nargs = "+",
			default = 2.0,
			help = "hyperparameter for prior on allele counts: beta [default:%(default)f]" )
parser.add_argument(	"-m", "--maf", type = float,
			default = 0.01,
			help = "minimum minor-allele frequency (MAF) at which to consider variant sites when counting alleles [default: %(default)f]" )
parser.add_argument(	"-S","--summarize", action = "store_true",
			default = False,
			help = "aggregate counts over the intervals in 'regions' file, instead of reporting k-mer-level counts" )
parser.add_argument(	"-x","--minhits", type = int,
			default = 2,
			help = "in summary mode, trim k-mers with (strictly) fewer than this many hits [default: %(default)s]" )
parser.add_argument(	"-X","--maxhits", type = int,
			default = 10000,
			help = "in summary mode, trim k-mers with (strictly) greater than this many hits [default: %(default)s]" )
parser.add_argument(	"-v","--verbose", action = "store_true",
			default = False,
			help = "print a status trace to stderr" )
parser.add_argument(	"-B", "--bootstrap", type = int,
			default = 0,
			help = "in summary mode, generate bootstrap CIs for count of zeros using this many replicates [default: %(default)s]" )
args = parser.parse_args()

## connect to msBWT(s)
bwt_dirs = args.msbwt
msbwt = util.load_bwts(bwt_dirs)

## fail if couldn't connect to at least one msBWT
if msbwt is None:
	sys.exit("No valid BWTs were provided; quitting.")
else:
	sys.stderr.write("Using the following BWTs:\n{}\n".format(str(bwt_dirs)))
	if args.summarize:
		sys.stderr.write("Working in summary mode.\n")

## try to connect to reference msBWT for normalizing
ref_bwt = None
if args.normalize is not None:
	try:
		ref_bwt = loadBWT(args.normalize)
	except Exception as e:
		print "Couldn't load reference BWT at <{}>".format(args.normalize)
		print e

if ref_bwt is not None:
	sys.stderr.write("Using the following BWT as reference:\n{}\n".format(args.normalize))

## get sequences of k-mers via pyfasta
fasta = pyfasta.Fasta(args.genome, key_fn = lambda key: key.split()[0])

## load regions
regions = pybedtools.BedTool(args.regions)
if args.step is None:
	args.step = args.kmer

if args.verbose:
	sys.stderr.write("Starting query loop.\n")

## simple function for bootstrap resampling of a np array
def resample(x):
	n = len(x)
	i = np.floor( np.random.rand(n)*len(x) ).astype(int)
	return x[i]

## apply some function to a bootstrap sample
def bootstrap(x, fn, **kwargs):
	return fn(resample(x, **kwargs))

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
	counts = []
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
		count_sum = 0
		count_ref = -1
		norm_count = 0
		for i in range(0, len(msbwt)):
			count_f += msbwt[i].countOccurrencesOfSeq( str(seq) )
			count_r += msbwt[i].countOccurrencesOfSeq( str(dna.revcomp(seq)) )
		count_sum = count_f + count_r
		counts.append(count_sum)

		## count occurrences in reference msBWT
		if args.normalize is not None:
			count_ref = ref_bwt.countOccurrencesOfSeq( str(seq) )
			if args.revcomp:
				count_ref += ref_bwt.countOccurrencesOfSeq( str(dna.revcomp(seq)) )

		## update aggregate counts; relevant for summary mode only
		if args.normalize is not None:
			if count_ref > 0:
				norm_count = float(count_sum)/count_ref
				agg.append(norm_count)
			else:
				# don't add this k-mer to aggregate count for this interval if it returned zero in reference BWT
				norm_count = "NA"
		else:
			agg.append(count_sum)

		## record this kmer as zero iff it returned zero across all BWTs
		if (count_sum) == 0:
			zeros += 1

		## compute allele scores, if requested
		## only makes sense in context of a single target msBWT, and for reads with reasonable number (1 < x < 5000?) of hits
		## NB: this is expensive, because it jumps all around FM-index to reconstruct reads
		(allele_score, nhaps, asymm) = ("NA","NA","NA")
		if args.alleles and len(msbwt) == 1 and count_sum > args.minhits and count_sum < args.maxhits:
			reads = util.get_reads(msbwt[0], str(seq), revcomp = True)
			reads.pseudoalign(k = args.kmer)
			if args.verbose:
				print reads.pretty_alignment()
				print reads.call_variant_sites(args.maf)
			(allele_score, haps) = reads.extract_haplotypes(args.maf)
			#(nhaps, asymm) = reads.consistency_score(args.maf)

		## fine-grained mode: print results for each k-mer
		if not args.summarize:
			print kmer.chrom, kmer.start, kmer.end, seq, count_f, count_r, count_sum, count_ref, norm_count, allele_score

		# except Exception as e:
		# 	# print "Sequence was: {}; searched failed.".format(seq)
		# 	print type(e)
		# 	pass

	## coarse-grained mode: print only the summary over the region
	if args.summarize:

		med, q, qqq, mu = ("NA", "NA", "NA", "NA")
		agg_trimmed = filter(lambda x: (x > args.minhits) and (x <= args.maxhits), agg)

		if len(agg_trimmed):
			med = np.median(agg_trimmed)
			if args.bootstrap == 0:
				## no bootstrapping; report order stats on k-mer counts
				q = np.percentile(agg_trimmed, 25.0)
				qqq = np.percentile(agg_trimmed, 75.0)
			else:
				## bootstrapping requested; report order stats on number of zeros
				x = [ bootstrap(np.array(counts), lambda z: sum(z == 0)) for i in range(0, args.bootstrap) ]
				q = np.percentile(x, 2.5)
				qqq = np.percentile(x, 97.5)
			mu = np.mean(agg_trimmed)

		print r.chrom, r.start, r.end, zeros, med, mu, q, qqq

	## clean up temp file
	try:
		os.remove(fn)
	except:
		pass

if args.verbose:
	sys.stderr.write("Done.\n")
