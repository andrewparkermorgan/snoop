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
import copy
import re
import numpy as np
import scipy.stats as sp

import MUSCython.MultiStringBWTCython as ms

from snoop import *

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

	def items(self):
		for i in range(0, len(self.seq)):
			yield Read(self.seq[i], self.offset[i], self.orientation[key])

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

	def extract_haplotypes(self, maf = 1):

		if self.alignment is not None:
			haps = haplotype.extract_haplotypes(self.alignment)
			threshold = maf
			if maf < 1.0:
				threshold = float(len(haps))*maf
			hapdict = {}
			nhaps = 0
			for (row, seqs, scores) in haps:
				hapdict.update( {row.seq: scores[0][1]} )
				nhaps += (scores[0][1] > threshold)

			self.haplotypes = hapdict
			self.nhaplotypes = nhaps
		return (self.nhaplotypes, self.haplotypes)

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
