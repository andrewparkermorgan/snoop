#! /usr/bin/env python

import os
import sys
import itertools
import collections
import argparse

import numpy as np
from sklearn.hmm import GaussianHMM
from sklearn.utils import check_random_state

from cogent import LoadSeqs, DNA
from snoop import util, rle, io, dna


## parse command-line arguments
parser = argparse.ArgumentParser(description = "Given an alignment of two sequences, take k-mers representing alternate alleles at variant sites, and use them to 'genotype' msBWTs.")
parser.add_argument(	"-a","--alignment", type = io.readable_file,
						help = "alignment of exactly two sequences in any reasonable format" )
parser.add_argument(	"-M","--msbwt", type = io.readable_dir, nargs = "+",
						help = "path(s) to one or more msBWTs" )
parser.add_argument(	"-k","--kmer", type = int,
						default = 30,
						help = "k-mer size in bp (ie. length of probe) [default: %(default)s]" )
parser.add_argument(	"-r","--refs", nargs = 2,
						default = None,
						help = "reference samples to define haplotypes [default: %(default)s]" )
parser.add_argument(	"-X","--coverage", type = float,
						default = 25.0,
						help = "prior mean coverage [default: %(default)f]" )
parser.add_argument(	"-x","--X_max", type = float,
						default = 3.0,
						help = "reject sites with coverage > X_MAX*COVERAGE [default: %(default)f]" )
parser.add_argument(	"--dryrun", type = int,
						default = 0,
						help = "do test run on first x bp of alignment; 0 to do full run [default: %(default)d]" )
parser.add_argument(	"--verbose", action = "store_true",
						help = "print counts to stderr [default: %(default)s]" )

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


def window_variant(aln, start, end, width = 0):

	offset = int((width - (end-start))/2)
	i_start = max(start - offset, 0)
	i_end = min(end + offset, len(aln))+1
	#print "\n----------\nprocessing a variant site {}-{}...".format(start, end)
	#print "expanding it by {}bp to {}-{}".format(width-(end-start), i_start, i_end)
	#print aln[ i_start:i_end ]
	seqs = collections.defaultdict(list)
	for s in aln.Names:
		seq_start = nearest_nongap(aln.getGappedSeq(s), i_start, -1)
		seq_end = nearest_nongap(aln.getGappedSeq(s), i_end, +1)
		extra = width - (seq_end - seq_start)
		#print "{}: extra with {}".format(s, extra)
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

	#print new_seqs.keys()
	lcs = longest_common_substring(new_seqs.keys())
	#print "LCS has length {}".format(len(lcs))

	if len(lcs) < shortest:
		return new_seqs



def invert_dict(d):

	dinv = collections.defaultdict(str)
	for k,v in d.iteritems():
		for x in v:
			dinv[x] = k

	return dinv

## load BWTs
bwts = util.BwtSet(args.msbwt[:1]) # allow only 1 msBWT for now
sys.stderr.write(str(bwts) + "\n")

## load alignment
aln = LoadSeqs(args.alignment, aligned = True, moltype = DNA)
if args.dryrun:
	aln = aln[ 0:min(args.dryrun, len(aln)) ]

## get start and end of variable regions
sys.stderr.write("finding variable sites in alignment...\n")
sites = rle.rle.from_ones(aln.variablePositions(), len(aln))
sys.stderr.write("doing greedy set cover...\n")
ivls = sites.greedy_set_cover(w = args.kmer, step = int(args.kmer/2), xmax = len(aln))

## check that specified reference sequences are in alignment
if args.refs:
	if not all([ r in aln.Names for r in args.refs ]):
		sys.exit("Not all specified reference sequences are present in the alignment.  Those present are:\n{}".format("\n".join(aln.Names)))

## initialize count arrays
if args.refs:
	ref_samples = args.refs
else:
	ref_samples = aln.Names

counts = np.ndarray([ len(ivls), len(ref_samples) ])
counts.fill(0.0)

## iterate on variant sites and collect k-mers representing each allele
for i in range(0, len(ivls)):

	start, end = ivls[i]

	try:

		kmers = window_variant(aln, start, end, args.kmer)
		if kmers:
			kmers_inv = invert_dict(kmers)
			alleles = set([ kmers_inv[r] for r in ref_samples ])
			if len(alleles) > 1:
				for j in range(0, len(ref_samples)):
					ref = ref_samples[j]
					k = kmers_inv[ref]
					for b,x in bwts.count(k, False).iteritems():
						if args.verbose:
							sys.stderr.write("{}\n".format("\t".join([ str(xx) for xx in [start, ref, k, b, x] ])))
						counts[i,j] = x

	except Exception as e:
		# just keep on going
		sys.stderr.write(str(e) + "\n")

	if not args.verbose and not i % 100:
		sys.stderr.write("{} of {} sites...\n".format(i, len(ivls)))

## initialize hmm parameters
rs = check_random_state(None) # fix RNG seed? maybe?

means = np.array([	[  0.0, 0.0 ],
					[ np.log1p(args.coverage), 0.0 ],
					[ 0.0, np.log1p(args.coverage) ],
					[ np.log1p(args.coverage/2), np.log1p(args.coverage/2) ],
					[ np.log1p(args.coverage), np.log1p(args.coverage) ] ])
cv = 1.0
covars = np.array([ [ 0.01, 0.01 ],
					[ cv, 0.01 ],
					[ 0.01, cv ],
					[ cv/2, cv/2 ],
					[ cv, cv ] ])
hidden = [ "private" ] + ref_samples + [ "heterozygous","pseudohet" ]

hmm = GaussianHMM(n_components = len(means), random_state = rs)
hmm._set_means(means)
hmm._set_covars(covars)

## filter sites; compute observation sequence as log(1+count)
keep = np.logical_and((counts.max(1) < args.X_max*args.coverage), (counts.sum(1) > -1.0))
counts = counts[ keep,: ]
obs = np.log1p(counts)
starts = np.array([ start for start,end in ivls ]).reshape( (len(ivls), 1) )
starts = starts[ keep,: ]

## run hmm
states =  hmm.decode(obs)

## print result to stdout
for i in range(0, counts.shape[0]):
	print starts[i,0], obs[i,0], obs[i,1], hidden[ states[1][i] ]
