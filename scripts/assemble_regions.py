#! /nas02/apps/python-2.7.1/bin/python
## /usr/bin/env python

## --- snoop/assemble_regions.py --- ##
##	Date:		5 Jan 2015
##	Udpated:	NA
##	Purpose:	Do greedy assembly (via msBWT) of regions specified in a bed file

import os
import sys
import argparse
import collections
import re
import copy

import pybedtools
import pyfasta

from snoop import util, io, dna, kmers, assembly

parser = argparse.ArgumentParser(description = "Perform targeted assembly using msBWT, starting from seeds taken from specified regions of a reference genome.",
				parents = [ assembly.assembly_args() ] )
parser.add_argument(	"-M","--msbwt", type = io.readable_dir, nargs = "+",
			default = "./",
			help = "directory containing msBWT components" )
parser.add_argument(	"-r","--regions", type = argparse.FileType("rU"),
			required = True,
			help = "list of regions (in bed format) for which to take seed sequences" )
parser.add_argument(	"-g","--genome", type = io.readable_file,
			default = "$GENOMES/mm10/mm10.fa",
			help = "reference genome fasta file from which to obtain seed sequences [default: %(default)s]" )
parser.add_argument(	"-l","--minlength", type = int,
			default = 500,
			help = "maximum length of final sequence to assemble, in bp [default: no limit]" )
parser.add_argument(	"-o","--out", type = io.writeable_dir,
			default = "./",
			help = "where to dump fastas" )
args = parser.parse_args()

## check command-line arguments
if not args.count_kmer:
	args.count_kmer = args.kmer

## load seed regions and reference sequence
regions = pybedtools.BedTool(args.regions)
fasta = pyfasta.Fasta(args.genome, key_fn = lambda key: key.split()[0])

## load msBWTs
bwtnames = [ re.sub(r"/+$","", x).split("/").pop() for x in args.msbwt ]
msbwts = util.load_bwts(args.msbwt)
sys.stderr.write("Using msBWTS: {}\n".format(str(bwtnames)))

## data structure to hold assembled sequences
finished = []
contigs = collections.OrderedDict()
rez = open(os.path.join(args.out, "result.bed"), "w")

## utilty function to dump assembly results for a given region
def dump_results(outfile, fastadir, contigs):
	for r, d in contigs.iteritems():
		with open(os.path.join(fastadir, "{}.fa".format(r)), "w") as outfasta:
			for s,seq in d["seqs"].iteritems():
				outfile.write( "{}:{}\t{}\t{}\n".format( "\t".join(str(d["region"]).strip().split()[0:5]), s, len(seq), "*" ) )
				outfasta.write(">{}\n".format(s))
				for km in kmers.kmerize(seq, 50, 50):
					outfasta.write("{}\n".format(km))

## loop on regions
## regions may be grouped by the 'name' field in bed file, in which case only one region is assembled per group
j = 0
last_group = None
for r in regions:

	## make fake group name, if region hasn't been assigned one
	j += 1
	if r.name is None:
		r.name = "seq{}".format(j)

	## check if we need to dump results
	if r.name != last_group and j > 1:
		dump_results(rez, args.out, contigs)
		contigs = collections.OrderedDict()
		finished += last_group
	last_group = r.name

	## skip this region if we've already assembled successfully in another member of its group
	if r.name in finished:
		continue

	## make sure seed sequences have proper form
	width = r.end - r.start
	if width < args.kmer:
		continue
	seed = str(dna.unmask( fasta[r.chrom][ r.start:r.end ] ))

	sys.stderr.write("--- {}: {} - {} ({})---\n".format(r.chrom, r.start, r.end, r.name))

	## loop on BWTs
	all_found = True
	these_contigs = collections.defaultdict(str)
	for i in range(0, len(bwtnames)):

		sys.stderr.write("\t{} ...\n".format(bwtnames[i]))

		## first check that seed will work in this BWT
		## allow step-down of seed size to get started, but hold k-mer size for assembly constant
		seed_found = False
		k_start = args.kmer + 1
		while not seed_found:
			k_start -= 1
			## if seed is less than 21nt long, it's probably hopeless
			if k_start < 21:
				break
			for seed_k in kmers.kmerize(seed, k_start, 1):
				x = util.count_reads(msbwts[i], seed_k)
				if x >= args.minweight:
					seed_found = True
					seed = seed_k
					sys.stderr.write("\t\tseed accepted ({} bp): {}\n".format(len(seed_k), seed_k))
					break

		## seed not found in this sample: break loop
		if not seed_found:
			all_found = False
			break

		## do assembly
		seq = assembly.greedy_assemble( args.msbwt[i], seed = seed, k = args.kmer, count_k = args.count_kmer, direction = True,
					min_weight = args.minweight, max_weight = args.maxweight,
					max_nodes = args.maxnodes, max_length = args.maxlength,
					save = False, outprefix = args.out, memmap = args.memmap )

		if len(seq) < args.minlength:
			sys.stderr.write("\t\tsequence was only {} bp long; rejecting.\n".format(len(seq)))
			all_found = False
			break
		else:
			sys.stderr.write("\t\t>{}\n".format(bwtnames[i]))
			for km in kmers.kmerize(seq, 50, 50):
				sys.stderr.write("\t\t{}\n".format(km))
			these_contigs[ bwtnames[i] ] = seq

	if all_found:
		contigs[r.name] = { "region": r, "seqs": these_contigs }
		#print(contigs)

rez.close()
