#! /usr/bin/env python

import os
import sys
import itertools
import collections
import argparse

from cogent import LoadSeqs, DNA
from snoop import rle, io


## parse command-line arguments
parser = argparse.ArgumentParser(description = "Given an alignment of two sequences, take k-mers representing alternate alleles at variant sites, and use them to 'genotype' msBWTs.")
parser.add_argument(	"-a","--alignment", type = io.readable_file,
						help = "alignment of exactly two sequences in any reasonable format" )
parser.add_argument(	"-k","--kmer", type = int,
						default = 30,
						help = "k-mer size in bp (ie. length of probe) [default: %(default)s]" )
parser.add_argument(	"-r","--refs", nargs = 2,
						default = None,
						help = "reference samples to define haplotypes [default: %(default)s]" )

args = parser.parse_args()

## get aligned base (non-gap) nearest a given index in alignment, looking forward (stride > 0) or backward (stride < 0)
def nearest_nongap(seq, aln_pos, stride = 1):
	
	if stride == 0:
		raise ValueError("Stride of 0 will give infinite loop.")

	posmap = seq.gapMaps()[1]
	i = aln_pos
	while i not in posmap and i < len(seq) and i >= 0:
		i += stride
	if i in posmap:
		return posmap[i]
	else:
		raise KeyError("Can't find sequence matching this aligned position after striding to position {}".format(i))


def longest_common_substring(ss):

	if len(ss) > 1 and len(ss[0]) > 0:
		substr = ""
		for i in range(len(ss[0])):
			for j in range(len(ss[0])-i+1):
				if j > len(substr) and all(ss[0][i:i+j] in x for x in ss):
					substr = ss[0][i:i+j]
		return substr
	else:
		return None


## for integers <x>, find the smallest set of intervals of width <w> which covers them all
## (with thanks to Leonard McMillan)
def greedy_set_cover(x, w = 30, buff = 5, xmax = 1e9):

	## set bounds of range containing candidate intervals
	x = set(x)
	lb = max(0, min(x)-buff)
	ub = min(xmax, max(x)+buff)

	## make candidate intervals; initialize empty set of accepted ones
	candidates = [ range(i,i+w) for i in range(lb, ub) ]
	accepted = []

	## keep going while there are elements left to be covered
	while len(x) and len(candidates):

		## list of lists of positions in <x> which each candidate hits
		hits = []
		## how many such hits for each candidate
		nhits = []
		for j in range(0, len(candidates)):
			hits.append( [y in x for y in candidates[j]] )
			nhits.append( sum(hits[ len(hits)-1 ]) )

		## take candidate which hit the most (ignore ties)
		best = nhits.index( max(nhits) )
		
		## remove all elements in <x> covered by best candidate
		x = x - set(candidates[best])

		## add this candidate to the list of accepted intervals
		this_best = candidates.pop(best)
		accepted.append( (min(this_best), max(this_best)) )

	return (accepted, list(x))


def window_variant(aln, start, end, width = 0):

	offset = int((width - (end-start))/2)
	i_start = max(start - offset, 0)
	i_end = min(end + offset, len(aln))+1
	print "\n----------\nprocessing a variant site {}-{}...".format(start, end)
	print "expanding it by {}bp to {}-{}".format(width-(end-start), i_start, i_end)
	print aln[ i_start:i_end ]
	seqs = collections.defaultdict(list)
	for s in aln.Names:
		seq_start = nearest_nongap(aln.getGappedSeq(s), i_start, -1)
		seq_end = nearest_nongap(aln.getGappedSeq(s), i_end, +1)
		extra = width - (seq_end - seq_start)
		print "{}: extra with {}".format(s, extra)
		if extra > 0:
			new_offset = int(extra/2)
		else:
			new_offset = 0
		seq_start = max(seq_start - new_offset, 0)
		seq_end = min(seq_end + new_offset, len(aln.getSeq(s))-1)
		kmer = str(aln.getSeq(s)[ seq_start:(seq_end+1) ])
		seqs[kmer].append(s)

	
	shortest = min([ len(s) for s in seqs.keys() ])
	new_seqs = collections.defaultdict(list)
	for s in seqs.keys():
		new_seqs[s[ 0:shortest] ] = seqs[s]
	
	print new_seqs.keys()
	lcs = longest_common_substring(new_seqs.keys())
	print "LCS has length {}".format(len(lcs))

	if len(lcs) < shortest:
		return new_seqs


def invert_dict(d):

	dinv = collections.defaultdict(str)
	for k,v in d.iteritems():
		for x in v:
			dinv[x] = k

	return dinv


## load alignment
aln = LoadSeqs(args.alignment, aligned = True, moltype = DNA)

## initialize ref samples
if args.refs:
	ref_samples = args.refs
else:
	ref_samples = aln.Names

## get start and end of variable regions
sites = rle.rle.from_ones(aln.variablePositions(), len(aln))
print sites.reconstitute()
ivls = greedy_set_cover([ x for x,y in sites.filter_values() ], w = args.kmer, xmax = len(aln))

## check that specified reference sequences are in alignment
if args.refs:
	if not all([ r in aln.Names for r in args.refs ]):
		sys.exit("Not all specified reference sequences are present in the alignment.  Those present are:\n{}".format("\n".join(aln.Names)))

## iterate on variant sites and collect k-mers representing each allele
for i in range(0, len(ivls[0])):

	start = ivls[0][i][0]
	end = ivls[0][i][1]
	
	try:

		kmers = window_variant(aln, start, end, args.kmer)
		if kmers:
			kmers_inv = invert_dict(kmers)
			alleles = set([ kmers_inv[r] for r in ref_samples ])
	
			if len(alleles) > 1:
				print kmers_inv

	except Exception as e:
		# don't worry about edge effects for now
		print e
