#! /usr/bin/env python

## --- snoop/util.py --- ##
##	Date: 12 Oct 2014
##	Updated: 31 Oct 2014
##	Purpose: Utility functions for retrieving and manipulating reads from msBWTs

import os
import sys
import logging
import collections
import itertools
import operator
import fnmatch
import copy
import re
import numpy as np
import scipy.stats as sp

import MUSCython.MultiStringBWTCython as ms

#from cogent import LoadSeqs, DNA
#from cogent.evolve.pairwise_distance import TN93Pair
#from cogent.app import muscle_v38 as muscle
#from cogent.core.alignment import SequenceCollection, Alignment
#from cogent.parse.fasta import MinimalFastaParser

#from Bio.Cluster import cluster

from snoop import io, dna, haplotype, rle

ORIENTATION_FWD = 0
ORIENTATION_REV = 1

## stuff related to logging
LOG_STREAM = sys.stderr

def _init_logger(stream = LOG_STREAM):
	'''
	This code taken from Shunping's Lapels for initializing a logger
	'''
	logger = logging.getLogger('root')
	logger.setLevel(logging.DEBUG)
	ch = logging.StreamHandler(stream)
	ch.setLevel(logging.INFO)
	formatter = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s", "%Y-%m-%d %H:%M:%S")
	ch.setFormatter(formatter)
	logger.addHandler(ch)
	return logger

## container for doing stuff on multiple msBWTs
class BwtSet:

	def __init__(self, bwt_paths):

		self._bwts = []
		try:
			self._bwts = self._load_bwts(bwt_paths)
			self.sizes = { b: self._bwts[b].getSymbolCount(0) for b in self._bwts.keys() }
			self.names = self._bwts.keys()
		except Exception as e:
			raise ValueError("Coudln't load one or more msBWTs because: {}".format(str(e)))

	def __len__(self):
		return len(self._bwts)

	def __getitem__(self, key):
		return self._bwts[key]

	def __repr__(self):
		rez = "BwtSet with {} msBWTs:\n".format(len(self._bwts))
		for b in self._bwts.keys():
			rez += "\t{}: {} ({:.2f}M reads)\n".format(b, str(self._bwts[b]), self.sizes[b]/1e6)
		return rez

	def _load_bwts(self, bwt_dirs):

		msbwt = {}
		for ff in bwt_dirs:
			if not io.readable_dir(ff):
				continue
			name = os.path.basename(ff.rstrip("/"))
			msbwt.update( { name: ms.loadBWT(ff) } )

		if len(msbwt):
			return msbwt
		else:
			return None

	def count(self, seq, normalize = False):
		counts = collections.defaultdict(float)
		for b in self.names:
			if not normalize:
				sz = 1
			else:
				sz = self.sizes[b]/1e9
			counts[b] += count_reads(self._bwts[b], seq)/sz
		return counts

## abstract container for a single read as stored in a msBWT
class Read:

	def __init__(self, raw_seq = None, dollar = None, orientation = ORIENTATION_FWD):

		self._raw_seq = None
		self.offset = None
		self.orientation = orientation
		self.seq = None
		self.has_offset = False

		if raw_seq is not None:

			## check if we need to resurrect the read
			if "$" in raw_seq and dollar is None:
				self._raw_seq = raw_seq
				(seq, dollar) = _unrotate_read(raw_seq)
				self.seq = seq
				self.offset = dollar
				self.has_offset = True
			else:
				self.seq = raw_seq
				self.offset = dollar
				if dollar > len(raw_seq) - 1:
					raise IndexError
				else:
					self.has_offset = True

	def __len__(self):
		if self.seq is not None:
			return len(self.seq)
		else:
			return None

	def __repr__(self):
		return self.__str__()

	def __str__(self):

		ori = { str(ORIENTATION_FWD): "forward", str(ORIENTATION_REV): "reverse" }

		pretty = "Read object in {} orientation ".format(ori[self.orientation])
		if self.seq is None:
			pretty += "with no sequence\n"
		else:
			pretty += "of length {} bp".format(len(self.seq))
			if self.has_offset:
				pretty += " and offset {}".format(self.offset)
			pretty += "\n"
			pretty += "\tseq: {}\n".format(self.seq)
			if self._raw_seq is not None:
				pretty += "\traw: {}\n".format(self._raw_seq)

		return pretty

	def revcomp(self):
		if self.has_offset:
			self.offset = len(self.seq) - self.offset
		if self.seq is not None:
			self.seq = dna.revcomp(str(self.seq))


## abstract container for a set of Read objects
class ReadSet:

	## instantiate a ReadSet from an iterator of sequences
	def __init__(self, read_iter = [], alphabet = "ACGT", spacer = "-"):

		self.alphabet = alphabet
		self.alignment = None
		self.count = 0
		self.seq = []
		self.offset = []
		self.orientation = []
		self.score = None
		self.nhaplotypes = None
		self.haplotypes = None
		self._spacer = spacer
		self.pivot = None

		## iterate over input
		for obj in read_iter:
			if isinstance(obj, Read):
				## this is a read object; append it to set
				self.append(obj)
			elif len(obj) >= 2:
				## assume this is a tuple of (sequence, dollar-sign index, [orientation])
				self.seq.append(obj[0])
				self.offset.append(obj[1])
				if len(obj) == 3:
					self.orientation.append(obj[3])

		self.count = len(self.seq)


	## permit deep copy
	#def __deepcopy__(self):
	#	pass

	def _add_reads(self, reads):
		if reads.count > 0:
			self.seq.extend( reads.seq )
			self.offset.extend( reads.offset )
			self.orientation.extend( reads.orientation )
			self.count = len(self.seq)

	## allow concatenation of ReadSets using normal '+' operator
	def __add__(self, reads):
		newself = copy.deepcopy(self)
		newself._add_reads(reads)
		return newself

	## all these methods allow ReadSets to act as sequence-type objects
	def __len__(self):
		return self.count

	def __getitem__(self, key):
		if key < len(self.seq) and key >= 0:
			return Read(self.seq[key], self.offset[key], self.orientation[key])
		else:
			raise IndexError

	def __setitem__(self, key, read):
		if key < len(self.seq) and key >= 0:
			self.seq[key] = read.seq
			self.offset[key] = read.offset

	def _iter_reads(self):
		for i in range(0, len(self)):
			this_read = self.seq[i]
			if self.orientation[i] == ORIENTATION_REV:
				this_read = dna.revcomp(self.seq[i])
			yield this_read

	def items(self):
		for i in range(0, len(self.seq)):
			yield Read(self.seq[i], self.offset[i], self.orientation[i])

	def append(self, read):
		self.seq.append(read.seq)
		self.offset.append(read.offset)
		self.orientation.append(read.orientation)
		self.count += 1

	def extend(self, reads):
		self._add_reads(reads)

	def pop(self):
		self.count -= 1
		return Read(self.seq.pop(), self.offset.pop(), self.orientation.pop())

	def pseudoalign(self, k = 0, spacer = None):
		if spacer is not None and len(spacer):
			self._spacer = str(spacer)[0]
		if self.alignment is None:
			(self.alignment, self.pivot) = _pseudoalign(self.seq, self.offset, self.orientation, k, self._spacer)
		return self.alignment

	def align(self, force = True):
		if self.alignment is None or force:
			self.alignment = _align_reads(list(self._iter_reads()))
		return self.alignment

	def consistency_score(self, maf = 0.0, eps = 0.00001):
		if self.alignment is None:
			raise AttributeError
		if self.score is None:
			self.score = _consistency_score(self.alignment, maf, self.alphabet, eps)
		return self.score

	def call_variant_sites(self, maf = 0.0):
		if self.alignment is None:
			raise AttributeError
		return _call_variant_sites(self.alignment, maf, self.alphabet)

	def count_haplotypes(self, alpha = 1, beta = 1):
		if self.alignment is None:
			raise AttributeError
		(n, hapdict) = _count_haplotypes(self.alignment, alpha = alpha, beta = beta, alphabet = self.alphabet)
		self.nhaplotypes = n
		self.haplotypes = hapdict
		return self.nhaplotypes

	def simply_count_haplotypes(self, maf, dist = 100):
		if self.alignment is None:
			raise AttributeError
		(n, hapdict) = _simple_count_haplotypes(self.alignment, maf = maf, dist = dist, alphabet = self.alphabet)
		self.nhaplotypes = n
		self.haplotypes = hapdict
		return self.nhaplotypes

	def pretty_alignment(self):
		pretty = ""
		if self.alignment is not None:
			for j in range(0, len(self.alignment)):
				pretty += "\n" + "".join(self.alignment[j][:])

		return pretty

	#def extract_haplotypes(self, maf = 1):
	#
	#	if self.alignment is not None:
	#		haps = haplotype.extract_haplotypes(self.alignment)
	#		threshold = maf
	#		if maf < 1.0:
	#			threshold = float(len(haps))*maf
	#		hapdict = {}
	#		nhaps = 0
	#		for (row, seqs, scores) in haps:
	#			hapdict.update( {row.seq: scores[0][1]} )
	#			nhaps += (scores[0][1] > threshold)
	#
	#		self.haplotypes = hapdict
	#		self.nhaplotypes = nhaps
	#	return (self.nhaplotypes, self.haplotypes)

	## find variable sites in a true MSA of reads
	## return value a tuple: (pretty-printed indicator string, list of variable ranges, number of callable sites)
	def find_variant_sites(self, maf = 1, min_coverage = 5, spacer = None):

		if self.alignment is None:
			self.align(force = True)
			#print self.alignment
		# since reads don't all start at same site, next part is required to discriminate between true gaps and padding
		self.alignment = _mask_padding_gaps(self.alignment, spacer)

		if spacer is None:
			spacer = self._spacer

		passing = []
		variable = []
		pos_freq = self.alignment.columnFreqs()
		for i in range(0, len(pos_freq)):
			col = pos_freq[i]
			col_filtered = _freqs_no_gaps(col, self.alphabet + "-")
			total = _dict_sum_if(col_filtered)
			if total >= min_coverage:
				passing.append(i)
			variable.append( int(len(filter(lambda x: x >= maf, col_filtered.values())) > 1) )

		## convert to run-length encoding in order to avoid over-counting: consecutive variants are assumed a single event
		## TODO: impose the consecutive-event rule only for gaps; allow consecutive substitutions to be independent
		runs = rle.rle(variable)
		variable_pos = runs.filter_values()
		return (variable, variable_pos, len(passing))

	## given a list of variable ranges, create hapltoypes and count their abundance
	def extract_haplotypes(self, varpos, kmer = None, offset = None, maf = 1, allow_missing = True, missing_char = "?"):

		if self.alignment is None:
			raise ValueError

		#if not offset:
		#	if not kmer:
		#		raise ArgumentError("Must specify a k-mer or an offset in order to split alignment into halves.")
		#	else:
		#		offset = min([ str(s).find(kmer) for s in self.alignment.iterSeqs() ])
		#		if offset < 0:
		#			raise ArgumentError("k-mer not found in alignment.")
		#		offset = offset + len(kmer)/2

		## split alignment
		#alns = [ self.alignment[ :offset ], self.alignment[ offset: ] ]

		#for a in alns:
		haps = []
		for seq in self.alignment.iterSeqs():
			s = str(seq)
			this_hap = "".join([ s[i] for i,j in varpos ])
			if not (missing_char in this_hap) or allow_missing:
				#haps[this_hap] += 1
				haps.append(this_hap)

		unique_haps = _uniquify_haps(haps)

		## tally haplotypes
		counts = collections.defaultdict(int)
		for h in haps:
			matched = []
			unmatched = 0
			for u in unique_haps:
				if _match_with_missing(h, u):
					matched.append(u)
			## add to haplotype totals
			for u in matched:
				## if ambiguous match, put it in unmmatched pile
				if len(matched) > 1:
					unmatched += 1
				else:
					counts[u] += 1

			## allocate unmatched reads in proportion to matched reads (EM-style)
			unambig_total = _dict_sum_if(counts)
			for u in counts.keys():
				counts[u] += (float(counts[u])/unambig_total)*unmatched

		return _dict_filter(counts, lambda x: counts[x] >= maf)


	## cluster reads by pairwise distance, dump small clusters, and return consensus and count in each cluster
	## TODO: how to decide number of clusters???
	#def cluster_reads(self, min_coverage = 5):
	#
	#	raise NotImplementedError
	#
	#	if self.alignment is None:
	#		self.align(force = True)
	#
	#	dist_calc = TN93Pair(DNA, alignment = self.alignment)
	#	dist_calc.run()
	#	dist_mat = dist_calc.getPairwiseDistances()
	#
	#	clus = cluster.treecluster(distancematrix = dist_mat)


## convert leading and trailing gaps in an MSA to a special indicator
def _mask_padding_gaps(aln, mask = "?"):

	## iterate on sequences in alignment
	seqs = []
	names = []
	for seq in aln.iterSeqs():
		s = str(seq)
		m1 = re.search(r"^\-+", s)
		m2 = re.search(r"\-+$", s)
		## mask gaps preceding real sequence (left-padding)
		if m1:
			first_base = m1.span()[1]
			s = mask*first_base + s[first_base:]
		## mask gaps following real sequence (right-padding)
		if m2:
			last_base = m2.start()
			s = s[:last_base] + mask*(m2.span()[1] - m2.span()[0])
		seqs.append(s)
		names.append(seq.Name)

	## convert masked sequences to new alignment object and return it
	mt = seq.MolType
	rez = LoadSeqs(data = zip(names, seqs), moltype = mt, aligned = True)
	return rez

## sum values in a dictionary, ignoring keys, possibly applying a filtering function
def _dict_sum_if(d, criteria = None):
	return sum( filter(criteria, d.values()) )

def _dict_filter(d, criteria = None):
	k = filter(criteria, d.keys())
	return { i: d[i] for i in k }

## get counts of non-gap bases from count dictionary
def _freqs_no_gaps(pos_freq, alphabet = "ACTG"):
	out_dict = {}
	for k,v in pos_freq.iteritems():
		if k in alphabet:
			out_dict.update({k: v})
	return out_dict

## merge strings liberally to get rid of missing values
def _uniquify_haps(seqs):

	unique_haps = set(seqs)
	flags = [ "?" in x for x in unique_haps ]
	len_previous = len(unique_haps) + 1

	## iteratively merge candidate hapltoypes
	while any(flags) and len_previous > len(unique_haps):
		len_previous = len(unique_haps)
		combos = itertools.combinations(unique_haps, 2)
		unique_haps = set()
		used = set()
		for a,b in combos:
			matched = _match_with_missing(a,b)
			if matched:
				ab = _merge_with_missing(a,b)
				unique_haps.add(ab)
				for x in [a,b]:
					if "?" in x:
						used.add(x)
			else:
				unique_haps.update([a,b])
		unique_haps = unique_haps - used
		flags = [ "?" in x for x in unique_haps ]

	## go back and check that, for every pair of consecutive variant sites, that pair was actually observed
	observed_haps = set()
	for u in unique_haps:
		observed = False
		if len(u) == 1:
			observed_haps.add(u)
		else:
			for i in range(0, len(u)-1):
				observed = observed or (u[ i:i+1 ] in [ x[ i:i+1 ] for x in set(seqs) ])
			if observed:
				observed_haps.add(u)

	return observed_haps

def _argmin(x):
	val, idx = min([ (val, idx) for (idx, val) in enumerate(x) ])
	return idx

def _match_with_missing(s1, s2, missing_chars = "?*.N"):
	if len(s1) != len(s2):
		raise ValueError
	score = 0
	for i in range(0, len(s1)):
		if (s1[i] not in missing_chars) and (s2[i] not in missing_chars):
			score += int(s1[i] != s2[i])
	return score == 0

def _merge_with_missing(s1, s2, missing_chars = "?*.N"):

	if len(s1) != len(s2):
		raise ValueError
	rez = ""
	for i in range(0, len(s1)):
		if s1[i] in missing_chars:
			if s2[i] in missing_chars:
				rez += "?"
			else:
				rez += s2[i]
		else:
			if (s2[i] in missing_chars) or (s2[i] == s1[i]):
				rez += s1[i]
			else:
				raise ValueError

	return rez

class PseudoalignedRow:

	def __init__(self, seq = "", alphabet = "ACGT"):
		self.seq = seq
		self._alphabet = alphabet
		self.nonmissing = sum([ x in alphabet for x in seq ])

	def __len__(self):
		return len(self.seq)

	def __eq__(self, other):

		if len(self) != len(other):
			raise ArgumentError
		else:
			flag = True
			for i in range(0, len(self)):
				if (self.seq[i] in self._alphabet and other.seq[i] in self._alphabet):
					flag = flag and (self.seq[i] == other.seq[i])
			return flag

	def __str__(self):
		return "".join(self.seq) + "\n"

	## merge two aligned rows, assuming 'other' < 'self' (in the sense that the first nonmissing char in self occurs before that in other)
	def merge(self, other):

		if not isinstance(other, PseudoalignedRow):
			raise ArgumentError
		if len(self) != len(other):
			raise ArgumentError

		#print other.seq
		#print self.seq
		pattern = r"[{}]".format(re.escape(self._alphabet))
		m = re.search(pattern, self.seq)
		if m:
			right = self.seq[ m.start(): ]
			left = other.seq[ :m.start() ]
			return PseudoalignedRow(left + right, self._alphabet)
		else:
			## nothing to merge with; return self unchanged
			return self

	## return indices of differences between two aligned rows
	def diffs(self, other):

		if not isinstance(other, PseudoalignedRow):
			raise ArgumentError
		if len(self) != len(other):
			raise ArgumentError

		diffs = []
		for i in range(0, len(self)):
			if (self.seq[i] in self._alphabet and other.seq[i] in self._alphabet):
				if self.seq[i] != other.seq[i]:
					diffs.append(i)

		if len(diffs):
			return diffs
		else:
			return None


## wrapper to load BWTs residing on a given list of paths, and return them as a list of objects
## returns None if no BWTs were loaded successfully
def load_bwts(bwt_dirs):

	msbwt = []
	for ff in bwt_dirs:
		if not io.readable_dir(ff):
			continue
		try:
			msbwt.append( ms.loadBWT(ff) )
		except Exception as e:
			sys.stderr.write("Couldn't load BWT at <{}>\n".format(ff))
			sys.stderr.write(str(e))

	if len(msbwt):
		return msbwt
	else:
		return None

## resurrect original read from the BWTed version: eg 'AC$GT' > 'GTAC'
## returns tuple of (resurrected read, dollar-sign index)
def _unrotate_read(raw_seq = None):

	if raw_seq is None:
		return None
	else:
		raw_seq = str(raw_seq)
		dollar = raw_seq.find("$")
		return ( raw_seq[ (dollar+1): ] + raw_seq[ 0:dollar ], dollar )

## given a bwt and a query string, count occurrences in both orientations
## returns the count (as integer?)
def count_reads(bwt, query, revcomp = True):

	count = 0
	query = str(query).upper()
	count += bwt.countOccurrencesOfSeq(query)
	if revcomp:
		count += bwt.countOccurrencesOfSeq( dna.revcomp(query) )
	return count

## given a bwt and a query string, return the (resurrected) hits
## returns a ReadSet object
def get_reads(bwt, query, revcomp = False):

	queries = [ str(query).upper() ]
	ori = [ ORIENTATION_FWD ]
	if revcomp:
		queries.append( dna.revcomp(str(query).upper()) )
		ori.append(ORIENTATION_REV)

	try:
		reads = ReadSet()

		for q in range(0, len(queries)):
			(s, e) = bwt.findIndicesOfStr(queries[q])
			for i in range(s,e):
				reads.append( Read(bwt.recoverString(i), orientation = ori[q]) )

		return reads

	except Exception as e:
		print e
		return None

## given sequences, mulitple-align them with MUSCLE
## NB: this will probably fail if not enough reads and/or not enough overlap between them
## NB: requires muscle >=v3.8 in $PATH
def _align_reads(reads):

	seqs = LoadSeqs(data = [ ("read{}".format(i), reads[i]) for i in range(0, len(reads)) ],
					moltype = DNA, aligned = False)
	rez = muscle.muscle_seqs(str(seqs))
	alignment = dict(MinimalFastaParser(rez["StdOut"].readlines()))
	new_alignment = Alignment(alignment, MolType = DNA)
	#print new_alignment
	return new_alignment

## given lists of (resurrected) reads and dollar-indices, 'pseudoalign' them by centering them on the original query
def _pseudoalign(reads, dollars, ori, k = 0, spacer = "-"):

	if not len(reads):
		return (None, -1)

	nreads = len(reads)
	lens = [ len(r) for r in reads ]
	maxlen = max(lens)

	#print dollars
	for i in range(0, nreads):
		if ori[i] == ORIENTATION_REV:
			dollars[i] = lens[i] - dollars[i] + k
			reads[i] = dna.revcomp(str(reads[i]))
		#dollars[i] = dollars[i] + maxlen - len(reads[i])

	lefts = [ l-d for (d,l) in zip(dollars, lens) ]
	rights = [ d-k for d in dollars ]
	pivot = max(lefts)

	## pad reads on either size of target k-mer so tat all have same padded width
	aln = [ spacer for i in range(0, nreads) ]
	for i in range(0, nreads):
		aln[i] = spacer*(max(lefts) - lefts[i]) + reads[i] + spacer*(maxlen - dollars[i])

	## sort rows in alignment, then return it
	aln.sort()
	aln.reverse()
	return (aln, pivot)

## compute weighted Shannon entropy of a 1d numpy array
## alpha,beta are lists of parameters of negative binomial prior on haplotype counts
def _shannon_entropy(col, alpha, beta):

	if not isinstance(col, (np.ndarray,)):
		col = np.array(col, shape = len(col), dtype = float)

	if len(alpha) != len(beta):
		return ArgumentError

	priors = []
	for i in range(0, len(alpha)):
		priors.append( sp.nbinom(alpha[i], beta[i]) )

	col = np.sort(np.around(col))
	if col.shape[0] == 1:
		return 0.0
	else:
		## apply prior scaling
		weights = np.zeros_like(col, dtype = float)
		for i in range(0, col.shape[0]):
			weights[i] = np.max([ p.pmf(col[i]) for p in priors ])
		#print col
		#print weights
		freq = col * weights
		freq = freq/np.sum(freq)
		#print freq
		H = -1 * freq * np.log2(freq)
		#print H
		return np.nansum(H)

## TODO: FIX THIS; it is merging haplotypes that should be distinct
## merge sequences in a pseudoalignment into pseudohaplotypes and return their counts
def _merge_rows(aln, alphabet = "ACGT"):

	nseq = len(aln)
	if nseq == 1:
		return aln[0]
	else:
		j = 1
		seqs = [ PseudoalignedRow("".join(aln[0])) ]
		counts = [1.0]
		for i in range(1, nseq):
			#print "checking row {}; there are {} sequences now".format(i,j)
			row = PseudoalignedRow("".join(aln[i]))
			matches = []
			for s in range(0, j):
				if row == seqs[s]:
					seqs[s] = row.merge(seqs[s])
					matches.append(s)
			if len(matches):
				#print "there is a match"
				if len(matches) == 1:
					for m in matches:
						counts[m] += 1.0/len(matches)
			#if len(matches):
				#match_weights = [ counts[x] for x in matches ]
				#addto = match_weights.index(max(match_weights))
				#counts[addto] += 1
			else:
				#print "inserting new sequence into cache"
				#print row
				#print seqs[s]
				j += 1
				seqs.append(row)
				counts.append(1.0)

	return (counts, [ s.seq for s in seqs])

## estimate number of haplotypes present in reads, with negative binomial prior
def _count_haplotypes(aln, alpha = [1], beta = [1], alphabet = "ACGT"):

	nseq = len(aln)
	if nseq == 1:
		return 1.0
	else:
		(counts, seqs) = _merge_rows(aln, alphabet)
		# print dict(zip(seqs, counts))
		counts = np.array(counts)
		return ( np.power(2, _shannon_entropy(counts, alpha, beta)), dict(zip(seqs, counts)) )

## do naive haplotype counting
def _simple_count_haplotypes(aln, maf = 1.0, dist = 100, alphabet = "ACGT"):

	maf = float(maf)
	nseq = len(aln)
	if nseq == 1:
		return (1, { "".join(aln[1]) : 1 })
	else:
		(counts, seqs) = _merge_rows(aln, alphabet)
		total = sum(counts)
		if maf < 1.0:
			maf = round(float(total)*maf)

		seqs.sort()
		seqs.reverse()
		newseqs = [ PseudoalignedRow(seqs[i]) for i in range(0, len(counts)) if counts[i] > maf ]
		newhaps = [ newseqs[0] ]
		newcounts = [1]
		firstdiff = -1
		j = 0
		while len(newseqs) > 1:
			d = newseqs[0].diffs(newseqs[1])
			if d is not None:
				if min(d) - firstdiff < dist - 1:
					j += 1
					newhaps.append(newseqs[1])
					newcounts.append(1)
				firstdiff = max(d)
			else:
				newcounts[j] += 1
			newseqs.pop(1)

		# seqdict = dict(zip(seqs, counts))
		# print { h.seq: seqdict[h.seq] for h in newhaps }
		return ( len(newhaps), dict(zip([ h.seq for h in newhaps ], newcounts)) )


## step through pseudoalignment of reads and find sites which have variant alleles above some specified frequency
def _call_variant_sites(aln, maf, alphabet = "ACGT"):

	maf = float(maf)
	nseq = len(aln)
	if nseq:
		width = len(aln[0])
		rez = [" "]*width
		maf_obs = []
		for i in range(0, width):
			seq = [ aln[j][i] for j in range(0, nseq) ]
			trimmed = [ x for x in seq if x in alphabet ]
			counts = collections.Counter(trimmed)
			total = float(sum(counts.values()))
			if len(trimmed):
				(nuc, maj) = counts.most_common()[0]
				if (maf > 1):
					maf_obs.append(total-maj)
					if (total - maj >= maf):
						rez[i] = "*"
				else:
					maf_obs.append(1.0 - maj/total)
					if ((1.0 - maj/total) > maf):
						rez[i] = "*"
		# print maf_obs
		return "".join(rez)
	else:
		return None

## compute some kind of 'consistency score' which approximates the number of haplotypes represented in a pseudoalignment
## TODO: come up with some theoretical justification for this
def _consistency_score(aln, maf, alphabet = "ACGT", eps = 0.001):

	nseq = len(aln)
	if nseq:

		width = len(aln[0])
		freq = np.empty(shape = (width, len(alphabet)), dtype = float)
		freq.fill(eps)
		weights = np.zeros(shape = width, dtype = float)

		for i in range(0, width):

			seq = [ aln[j][i] for j in range(0, nseq) ]
			trimmed = [ x for x in seq if x in alphabet ]
			counts = collections.Counter(trimmed)
			total = sum(counts.values())

			k = 0
			for (nuc, times) in counts.most_common():
				prob = float(times)/float(total)
				if (maf < 1.0):
					if prob > maf:
						freq[i][k] = prob
						k += 1
				else:
					if times > maf:
						freq[i][k] = prob
						k += 1
			weights[i] = total

		## drop alignment columns with no acceptable characters
		nzcols = (weights > eps)
		weights = weights[ nzcols ]
		freq = freq[ nzcols,: ]

		## compute composition score: some approximation to the number of alleles present???
		colsums = np.linalg.norm(freq, 1, 1)
		normfreq = freq/colsums.reshape(freq.shape[0], 1) ## column-wise frequencies
		colentropy = -1 * np.sum(normfreq * np.log2(normfreq), 1) ## column-wise Shannon information
		nalleles = np.power(2, colentropy) ## 'expected number of alleles', column-wise
		adjweights = weights * (nalleles-1) ## adjusted weights: columns with less info count less; columns w/1 allle are downweighted

		# print nzcols, adjweights

		## compute skew score: whether or not reads extend evenly from either side of query
		skew = np.argmax(weights)/float(weights.shape[0])*2 - 1

		return (np.average(nalleles, weights = adjweights), skew)

	else:
		return None
