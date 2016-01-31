#! /nas02/apps/python-2.7.1/bin/python
## /usr/bin/env python

## --- snoop/count_reference_kmers.py --- ##
##	Date: 16 April 2015
##	Udpated: NA
##	Purpose: count k-mers from sequences in a fasta file against given msBWTs

import os
import sys
import argparse

from Bio import SeqIO

from snoop import io, util, kmers

parser = argparse.ArgumentParser(description = "Utility for generating and counting occurences of k-mers from a reference sequence in a msBWT.")
parser.add_argument(	"-f","--fasta", type = argparse.FileType("rU"),
						help = "fasta file containing sequences to search with" )
parser.add_argument(	"-M", "--msbwt", type = io.readable_dir, nargs = "+",
						help = "one or more msBWTs in which to count k-mers" )
parser.add_argument(	"-k", "--kmer", type = int,
						default = 0,
						help = "k-mer size (set to 0 to use all of each sequence) [default: %(default)d]")
parser.add_argument(	"--normalize", action = "store_true",
						help = "normalize counts against total size of msBWTs [default: %(default)d]")
args = parser.parse_args()

bwts = util.BwtSet(args.msbwt)
fa = SeqIO.parse(args.fasta, format = "fasta")

sys.stderr.write( "Using the following msBWTs:\n{}".format(str(bwts)) )
sys.stderr.write( "Reporting counts as parts per billion.\n" )
if args.kmer > 0:
	sys.stderr.write( "Breaking sequences into k-mers of length {} for searches.\n".format(args.kmer) )
else:
	sys.stderr.write( "Searching with provided sequences as-is.\n" )

for seq in fa:

	if args.kmer > 0:
		k = args.kmer
	else:
		k = len(seq.seq)
	for subseq in kmers.kmerize(seq.seq, k):
		for bwtname,count in bwts.count(subseq, args.normalize).iteritems():
			print seq.name, subseq, bwtname, count
