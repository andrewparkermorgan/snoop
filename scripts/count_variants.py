#! /usr/bin/env python

## --- snoop/count_reference_kmers.py --- ##
##	Date: 9 Feb 2014
##	Udpated: NA
##	Purpose: use k-mers from reference sequence to pull out reads, then identify and count haplotypes in those reads outside the fixed k-mers

import os
import sys
import tempfile
import subprocess
import argparse

import MUSCython.MultiStringBWTCython as ms

import pyfasta
import pybedtools
from snoop import util, io, dna

parser = argparse.ArgumentParser(description = "Utility for generating and counting occurences of k-mers from a reference sequence in a msBWT.")
parser.add_argument(	"-r","--regions", type = argparse.FileType("rU"),
			required = True,
			help = "list of regions (in bed format) for which to make k-mers" )
parser.add_argument(	"-M","--msbwt", type = io.readable_dir, nargs = "+",
			required = True,
			help = "path to one or more directories containing msBWT components" )
parser.add_argument(	"-g","--genome", type = io.readable_file,
			default = "$GENOMES/mm10/mm10.fa",
			help = "reference genome fasta file from which to obtain k-mer sequences [default: %(default)s]" )
parser.add_argument(	"-k","--kmer", type = int,
			default = 30,
			help = "k-mer size in bp (ie. length of probe) [default: %(default)s]" )
parser.add_argument(	"-s","--step", type = int,
			help = "step size for making k-mers; if kmer_1 = (1, k), kmer_2 = (1+s, k+s) [default: same as k-mer size]" )
parser.add_argument(	"-f", "--maf", type = int,
			default = 1,
			help = "minimum 'minor allele frequency' (an integer number of reads) to accept site as variant [default:%(default)d]" )
parser.add_argument(	"-x", "--coverage", type = int,
			default = 5,
			help = "minimum number of reads covering a site to declare it callable [default:%(default)d]" )
parser.add_argument(	"-X", "--max-coverage", type = int,
			default = 1500,
			help = "maximum number of reads from a k-mer to attempt haplotype enumeration [default:%(default)d]" )
args = parser.parse_args()

## load msBWTs
msbwt = util.load_bwts(args.msbwt)

## fail if couldn't connect to at least one msBWT
if msbwt is None:
	sys.exit("No valid BWTs were provided; quitting.")
else:
	sys.stderr.write("Using the following BWTs:\n{}\n".format(str(args.msbwt)))
	sys.stderr.write("Coverage limits are [{}, {}] and MAF threshold is {}.".format(args.coverage, args.max_coverage, args.maf))

## get sequences of k-mers via pyfasta
fasta = pyfasta.Fasta(args.genome, key_fn = lambda key: key.split()[0])

## load regions
regions = pybedtools.BedTool(args.regions)
if args.step is None:
	args.step = args.kmer

## loop on regions
for r in regions:

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

	for kmer in kmers:

		## skip kmers which are clipped by region boundaries; this produces boundary artefact but I will ignore it
		if (kmer.end - kmer.start) != args.kmer:
			continue
		# try:
		seq = dna.unmask( fasta[kmer.chrom][ kmer.start:kmer.end ] )

		## loop over all open BWTs and query them with this sequence
		for i in range(0, len(msbwt)):

			bwtname = list( os.path.split(args.msbwt[i].rstrip("/")) ).pop()

			count_f = msbwt[i].countOccurrencesOfSeq( str(seq) )
			count_r = msbwt[i].countOccurrencesOfSeq( str(dna.revcomp(seq)) )

			if count_f + count_r <= args.max_coverage and count_f + count_r > 0:

				reads = util.get_reads(msbwt[i], seq, True)
				(nvar, varpos, nsites) = reads.find_variant_sites(args.maf, args.coverage, "?")
				haps = reads.extract_haplotypes(varpos, maf = args.maf)
				if nsites:
					print kmer.chrom, kmer.start, kmer.end, seq, dna.complexity(seq), float(len(varpos))/nsites
				#for h in haps:
				#	print kmer.chrom, kmer.start, kmer.end, seq, dna.complexity(seq), bwtname, h, haps[h]
