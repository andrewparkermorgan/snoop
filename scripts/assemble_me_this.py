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

from snoop import util, io, dna, kmers, assembly

parser = argparse.ArgumentParser(description = "Tool for doing simplistic/naive graph-based assembly using a msBWT, starting from a specified seed.",
				parents = [ assembly.assembly_args() ] )
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
args = parser.parse_args()

## check command-line arguments
if not args.count_kmer:
	args.count_kmer = args.kmer
if len(args.seed) < args.kmer:
	sys.exit("Starting seed must be at least as long as specified k-mer size.")
if args.kmer < 1 or args.count_kmer < 1:
	sys.exit("Size of k-mers for both graph construction and path-weight calculation must be strictly positive.")
save_graph = False
if args.prefix:
	save_graph = True

bwtname = re.sub(r"/+$","", args.msbwt).split("/").pop()

## make sure seed sequences have proper form
seed = dna.ungap(args.seed.upper())
sys.stderr.write("Seed sequence is: {}\n".format(seed))
end_seeds = []
if args.end_seeds:
	for s in args.end_seeds:
		end_seeds.append( dna.ungap(s) )

## examine the seed, taking the first k-mer which returns a result
msbwts = util.load_bwts([args.msbwt])
seed_found = False
for seed_k in kmers.kmerize(seed, args.kmer, 1):
	x = util.count_reads(msbwts[0], seed_k)
	if x >= args.minweight:
		seed_found = True
		seed = seed_k
		print seed_k
		break

if not seed_found:
	sys.exit("Seed sequence doesn't contain any k-mers which meet the specified abundance threshold.")

## do assembly
seq = assembly.greedy_assemble( args.msbwt, seed = seed, end_seeds = end_seeds, k = args.kmer, count_k = args.count_kmer, direction = args.backward,
				min_weight = args.minweight, max_weight = args.maxweight,
				max_nodes = args.maxnodes, max_length = args.maxlength,
				save = save_graph, outprefix = args.prefix, memmap = args.memmap )

## dump result to stdout
print ">{}:{}-{}".format(bwtname, 1, len(seq))
for piece in kmers.kmerize(seq, args.kmer, args.kmer):
	print piece
