#! /nas02/apps/python-2.7.1/bin/python
## /usr/bin/env python

## --- snoop/assemble_me_this.py --- ##
##	Date:		13 Nov 2014
##	Udpated:	NA
##	Purpose:	Do "greedy assembly" from a msBWT: first build de Bruijn graph, then traverse it by taking heaviest path at every branch point.

import os
import sys
import argparse
import re

from snoop import io, dna, kmers, assembly

parser = argparse.ArgumentParser(description = "Tool for doing simplistic/naive graph-based assembly using a msBWT, starting from a specified seed.")
parser.add_argument(	"-M","--msbwt", type = io.readable_dir,
			default = "./",
			help = "directory containing msBWT components" )
parser.add_argument(	"-s","--seed",
			required = True,
			help = "initiation seed -- assembly starts here" )
parser.add_argument(	"-e","--end_seeds", nargs = "+",
			help = "termination seed -- assembly stops if we see one of these" )
parser.add_argument(	"--backward", action = "store_false",
			help = "specify this option to assemble backward (ie. leftward) from seed; [default: forward (ie. rightward)]" )
parser.add_argument(	"-n","--maxnodes", type = int,
			default = 1000,
			help = "maximum size of assembly graph, in number of nodes [default: %(default)d]" )
parser.add_argument(	"-l","--maxlength", type = int,
			default = None,
			help = "maximum length of final sequence to assemble, in bp [default: %(default)d]" )
parser.add_argument(	"-x", "--minweight", type = int,
			default = 5,
			help = "minimum weight of path to keep assembling [default: %(default)d]" )
parser.add_argument(	"-X", "--maxweight", type = int,
			default = 3000,
			help = "maximum weight of path to keep assembling [default: %(default)d]" )
parser.add_argument(	"-k", "--kmer", type = int,
			default = 70,
			help = "k-mer size to use in building assembly graph [default: %(default)d]" )
parser.add_argument(	"--count_kmer", type = int,
			help = "k-mer size to use in calculating path weights [default: same as --kmer]" )
parser.add_argument(	"-o","--out",
			help = "prefix for output file: if specified, assembly graph will be dumped as pickle" )
parser.add_argument(	"--memmap", action = "store_true",
			help = "use memmapped files to initialize BWT, instead of loading into memory [default: False]" )
args = parser.parse_args()

## check command-line arguments
if not args.count_kmer:
	args.count_kmer = args.kmer
if len(args.seed) < args.kmer:
	sys.exit("Starting seed must be at least as long as specified k-mer size.")
if args.kmer < 1 or args.count_kmer < 1:
	sys.exit("Size of k-mers for both graph construction and path-weight calculation must be strictly positive.")
save_graph = False
if args.out:
	save_graph = True

bwtname = re.sub(r"/+$","", args.msbwt).split("/").pop()

## check assembly seeds
seed = dna.ungap(args.seed)
end_seeds = []
if args.end_seeds:
	for s in args.end_seeds:
		end_seeds.append( dna.ungap(s) )

## do assembly
seq = assembly.greedy_assemble( args.msbwt, seed = seed, end_seeds = end_seeds, k = args.kmer, count_k = args.count_kmer, direction = args.backward,
				min_weight = args.minweight, max_weight = args.maxweight,
				max_nodes = args.maxnodes, max_length = args.maxlength,
				save = save_graph, outprefix = args.out, memmap = args.memmap )

## dump result to stdout
print ">{}:{}-{}".format(bwtname, 1, len(seq))
for piece in kmers.kmerize(seq, 50, 50):
	print piece