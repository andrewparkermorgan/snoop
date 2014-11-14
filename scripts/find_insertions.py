#! /nas02/apps/python-2.7.1/bin/python
## /usr/bin/env python

## --- snoop/find_insertions.py --- ##
##	Date: 31 Oct 2014
##	Udpated: NA
##	Purpose: Step through reference sequence and find insertion sites in a msBWT.
##		Procedure has two steps: (1) find variant sites from pre-computed k-mer counts by HMM; (2) backtrack and assemble boundaries.
##		Input is assumed to be in bed format: chr start end name score; with 'name' = k-mer sequence.

import os
import sys
import argparse
import subprocess
import tempfile
import glob
import collections

import numpy as np
from sklearn.hmm import GaussianHMM

from snoop import util, io, dna, kmers

import pybedtools
import pyfasta
from MUSCython.MultiStringBWTCython import * # msBWT

parser = argparse.ArgumentParser(description = "Tool for finding suspected TE insertions by using msBWT to assemble inward from boundaries of non-reference sites.")
parser.add_argument(	"-i","--infile", type = argparse.FileType("rU"),
			required = True,
			help = "pre-computed k-mer counts in bed format (chr start end k-mer count)" )
parser.add_argument(	"-M","--msbwt",
			default = "./",
			help = "globbing expression to find directories containing msBWT components" )
parser.add_argument(	"-c", "--complexity", type = float,
			default = 0.7,
			help = "exclude k-mers with sequence complexity (~= compression ratio) below this threshold [default: %(default)f]" )
parser.add_argument(	"-X", "--maxhits", type = int,
			default = 2000,
			help = "exclude k-mers with counts above this threshold [default: %(default)d]" )
parser.add_argument(	"-m", "--maf", type = int,
			default = 1,
			help = "minimum minor-allele frequency (MAF, as absolute count) at which to accept a haplotype [default: %(default)d]" )
parser.add_argument(	"--mu", type = float, nargs = "*",
			default = [0, 3.4],
			help = "prior means for mixture components [default: %(default)f]" )
parser.add_argument(	"--sigmasq", type = float, nargs = "*",
			default = [0.001, 1.0],
			help = "prior variances for mixture components [default: %(default)f]" )
parser.add_argument(	"-v","--verbose", action = "store_true",
			default = False,
			help = "print a status trace to stderr" )
args = parser.parse_args()

## given a msBWT and left and right boundary k-mers of a variant, get consensus sequence(s) just "inside"
def assemble_inward(bwt, left, right):

	k = len(left)

	queries = [left, right]
	breakpoints = []
	offsets = [ 0,1 ]

	for i in range(0, len(queries)):

		reads = util.get_reads(bwt, queries[i])
		reads.pseudoalign(k = k)

		if reads.alignment is not None:
			(score, haps) = reads.extract_haplotypes(args.maf)

			bp = {}
			start = reads.pivot + k
			end = len(reads.alignment[0])
			if offsets[i]:
				start = 0
				end = reads.pivot

			#print reads.pretty_alignment()
			for seq,count in haps.iteritems():
				if count > args.maf:
					bp.update({ seq[start:end] : count })
					# print seq[start:end], count

			breakpoints.append(bp)

	return breakpoints

def lazy_mean(x):
	return float(sum(x))/len(x)

## connect to msBWT(s)
bwt_dirs = glob.glob(args.msbwt)
msbwt = util.load_bwts(bwt_dirs)

## fail if couldn't connect to at least one msBWT
if msbwt is None:
	sys.exit("No valid BWTs were provided; quitting.")
else:
	sys.stderr.write("Using the following BWTs:\n{}\n".format(str(bwt_dirs)))

## load input file
kmer_iter = pybedtools.BedTool(args.infile)

## loop on regions and read counts
if args.verbose:
	sys.stderr.write("Reading precomputed k-mer counts...\n")

counts = []
kmer_stash = []
kmer_size = 0
i = 0
for kmer in kmer_iter:

	seq = kmer.name
	if not kmer_size:
		kmer_size = len(seq)
	cplx = dna.complexity(seq)
	count = int(kmer.score)

	if cplx > args.complexity and count < args.maxhits:
		trimmed_count = count
		if count < args.maf:
			trimmed_count = 0
		counts.append(trimmed_count)
		kmer_stash.append(kmer)
	i += 1

if not len(counts):
	sys.exit("No k-mer counts remain after filtering; check thresholds and try again.")

## fit HMM to counts
if len(args.mu) != len(args.sigmasq):
	sys.exit("Vectors of prior means and variances must be same length.")

counts = np.reshape(np.log1p(np.array(counts, dtype = "int")), (-1,1))
hmm = GaussianHMM( len(args.mu) )
hmm.fit([counts])

if args.verbose:
	sys.stderr.write("Fitting HMM to k-mer counts, assuming {} hidden states...\n".format(len(args.mu)))
	sys.stderr.write("means:\n" + str(hmm.means_) + "\n")
	sys.stderr.write("covariances:\n" + str(hmm.covars_) + "\n")
	sys.stderr.write("\n")
	sys.stderr.write("Processing possible variant sites...\n")
	sys.stderr.write("\trejecting haplotypes with read count < {}\n".format(args.maf))
	sys.stderr.write("\taccepting as TE/ME any haplotype with max count > {}\n".format(args.maxhits))

## find positions of transitions
states =  hmm.predict(counts)
breaks = np.where(np.diff(states))[0]
break_coords = []
break_kmers = []
break_counts = []
for j in range(0, breaks.shape[0]/2):
	i_start = breaks[2*j]
	i_end = breaks[2*j+1]+1
	break_coords.append( (kmer_stash[i_start].chrom, kmer_stash[i_start].start, kmer_stash[i_end].start) )
	break_kmers.append( (kmer_stash[i_start].name, kmer_stash[i_end].name) )
	break_counts.append( (kmer_stash[i_start].score, kmer_stash[i_end].score) )

for i in range(0, len(break_coords)):
	hap = assemble_inward(msbwt[0], break_kmers[i][0], break_kmers[i][1])
	if args.verbose:
		outline = [ break_coords[i][0], break_coords[i][1], break_counts[i][0], break_coords[i][2], break_counts[i][1] ]
		sys.stderr.write("\t".join([ str(x) for x in outline ]) + "\n")
		sys.stderr.write("\t" + str(hap) + "\n\n")
	for h in hap:
		hap_avg = {}
		for (seq, count) in h.iteritems():
			if count > args.maf:
				ungapped_seq = dna.ungap(seq)
				k_counts = []
				for k in kmers.kmerize(dna.ungap(seq), kmer_size):
					if dna.complexity(k) > args.complexity:
						k_counts.append( util.count_reads(msbwt[0], k) )
				if len(k_counts):
					hap_avg.update({ ungapped_seq: max(k_counts) })

		flag = "*"
		if len(hap_avg.keys()) > 1:
			sys.stderr.write("Warning: apparently there is >1 haplotype at this variant site.\n")
			flag = "+"
		for (hh, hc) in hap_avg.iteritems():
			if hc > args.maxhits:
				print break_coords[i][0], break_coords[i][1], break_coords[i][2], hh, hc, flag
				sys.stdout.flush() # force write; this lets me peek at output in almost-real time

if args.verbose:
	sys.stderr.write("Done.\n")
