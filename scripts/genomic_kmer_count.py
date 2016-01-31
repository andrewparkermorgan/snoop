#! /usr/bin/env python

## --- snoop/genomic_kmer_count.py --- ##
##	Date: 28 June 2015
##	Udpated: NA
##	Purpose: count k-mers from reference sequence in a msBWT, aggregating in windows; compute bootstrapped estimate of sequence divergence

import os
import sys
import argparse
import subprocess
import tempfile

import pybedtools as pbt
from snoop import io

parser = argparse.ArgumentParser(description = "Count k-mers and average in windows (probably big) across a genome.")
parser.add_argument(	"-M", "--msbwt", type = io.readable_dir,
						default = "./",
						help = "path to directory containing msBWT components" )
parser.add_argument(	"-g","--genome", type = io.readable_file,
						default = "$GENOMES/mm10/mm10.fa",
						help = "reference genome fasta file from which to obtain k-mer sequences [default: %(default)s]" )
parser.add_argument(	"-k","--kmer", type = int,
						default = 31,
						help = "k-mer size in bp (ie. length of probe) [default: %(default)s]" )
parser.add_argument(	"-c", "--chrom", type = str, nargs = "*",
						help = "chromosome(s) to process; default is to do all" )
parser.add_argument(	"-w", "--window", type = int,
						default = 31000,
						help = "window size for calculations (bp)" )
parser.add_argument(	"-s", "--step", type = int,
						default = 15500,
						help = "step size between adjacent windows (bp)" )
parser.add_argument(	"-B", "--bootstrap", type = int,
						default = 0,
						help = "generate bootstrap CIs for count of zeros using this many replicates [default: %(default)s]" )
args = parser.parse_args()

## check input
if args.window % args.kmer:
	sys.exit("Window size should be an even multiple of k-mer size.")

## make temporary bed file with 1 entry per chromosome
(fh, fn) = tempfile.mkstemp()
with open(args.genome + ".fai", "r") as fai:
	for ll in fai:
		pieces = ll.strip().split()
		if pieces[0] in args.chrom or args.chrom is None:
			os.write(fh, "\t".join( [ pieces[0], str(3000000), pieces[1] ] ) + "\n")
os.close(fh)

## make windows on the chromosome
(winfh, winfn) = tempfile.mkstemp()
p = subprocess.Popen(	["windowMaker", "-b", fn, "-w", str(args.window), "-s", str(args.step), "-i", "srcwinnum"],
						stdout = subprocess.PIPE,
						stderr = subprocess.PIPE,
						bufsize = 0 )
for line in p.stdout:
	os.write(winfh, line)
os.close(winfh)

## launch command
bwt_id = os.path.basename(args.msbwt.strip("/"))
cmd = "count_reference_kmers.py -r {} -M {} -g {} -k {} --summarize --bootstrap {}"
cmd = cmd.format(winfn, args.msbwt, args.genome, args.kmer, args.bootstrap)
sys.stderr.write("-- Running command: '{}' --\n".format(cmd))
rez = subprocess.Popen(		cmd, shell = True,
							stdout = subprocess.PIPE,
							stderr = subprocess.PIPE,
							bufsize = 0 )

## slurp output, adding sample ID at the end of each line
for line in rez.stdout:
	print line.strip(), bwt_id
