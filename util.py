#! /usr/bin/env python

## --- snoop/util.py --- ##
##	Date: 12 Oct 2014
##	Updated: NA
##	Purpose: Utility functions for retrieving and manipulating reads from msBWTs

import os
import sys
import numpy as np

import MUSCython.MultiStringBWTCython as ms

from common import *

## abstract container for a single read as stored in a msBWT
class Read:

	def __init__(self, raw_seq = None, dollar = None):

		self._raw_seq = None
		self.offset = None
		self.seq = None
		self.has_offset = False

		if raw_seq is not None:

			## check if we need to resurrect the read
			if "$" in raw_seq and dollar is None:
				self._raw_seq = raw_seq
				(dollar, seq) = _unrotate_read(raw_seq)
				self.seq = seq
				self.offset = dollar
				self.has_offset = True
			else:
				self.seq = raw_seq
				self.offset = dollar
				if dollar > len(seq) - 1:
					raise IndexError
				else:
					self.has_offset = True

## abstract container for a set of Read objects
class ReadSet:

	def __init__(self, read_iter = [], alphabet = "ACGT", spacer = "-"):

		self.alphabet = alphabet
		self.alignment = None
		self.count = 0
		self.seq = []
		self.offset = []
		self.score = None
		self._spacer = spacer

		for (seq, dollar) in read_iter:
			self.seq.append(seq)
			self.offset.append(dollar)

		self.count = len(self.seq)

	def __len__(self):
		return self.count

	def __getitem__(self, key):
		if key < len(self.seq) and key > 0:
			return Read(self.seq[key], self.offset[key])
		else:
			raise IndexError

	def __setitem__(self, key, read):
		if key < len(self.seq) and key > 0:
			self.seq[key] = read.seq
			self.offset[key] = read.offset

	def items(self):
		for i in range(0, len(self.seq)):
			yield Read(self.seq[i], self.offset[i])

	def append(self, read):
		self.seq.append(read.seq)
		self.offset.append(read.offset)

	def pop(self):
		return Read(self.seq.pop(), self.offset.pop())

	def pseudoalign(self, spacer = None):
		if spacer is not None and len(spacer):
			self._spacer = str(spacer)[0]
		if self.alignment is None:
			self.alignment = _pseudoalign(self.seq, self.offset, self._spacer)
		return self.alignment

	def consistency_score(self, maf = 0.0, eps = 0.00001):
		if self.score is None:
			self.score = _consistency_score(self.alignment, maf, self.alphabet, eps)
		return self.score

	def call_variant_sites(self, maf = 0.0):
		if self.alignment is None:
			self.alignment = _pseudoalign(self.seq, self.offset, self._spacer)
		return _call_variant_sites(self.alignment, maf, self.alphabet)

## wrapper to load BWTs residing on a given list of paths, and return them as a list of objects
## returns None if no BWTs were loaded successfully
def load_bwts(bwt_dirs):

	msbwt = []
	for ff in bwt_dirs:
		if not readable_dir(ff):
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

## given a bwt and a query string, return the (resurrected) hits
## returns a ReadSet object
def get_reads(bwt, query):

	reads = ReadSet()

	try:
		(s, e) = bwt.findIndicesOfStr(str(query).upper())
		for i in range(s,e):
			read = Read(bwt.recoverString(i))
			reads.append(read)
		return reads

	except Exception as e:
		print e
		return None

## given lists of (resurrected) reads and dollar-indices, 'pseudoalign' them by centering them on the original query
def _pseudoalign(reads, dollars, spacer = "-"):

	nreads = len(reads)
	maxlen = max([ len(r) for r in reads ])

	stack = [ y for (x,y) in sorted(zip(dollars, range(0, nreads))) ]
	aln = [ [spacer]*(maxlen + max(dollars)-min(dollars)) for i in range(0, nreads) ]
	for ii in range(0, len(stack)):
		aln[stack[ii]][:] = spacer*(dollars[stack[ii]]-min(dollars)) + reads[stack[ii]] + "-"*(k-dollars[stack[ii]]-2)

	return aln

## step through pseudoalignment of reads and find sites which have variant alleles above some specified frequency
def _call_variant_sites(aln, maf, ignore = ["-","N"]):

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
				maf_obs.append(1-float(maj)/total)
				if ((1.0 - maj/total) > maf):
					rez[i] = "*"
		# print maf_obs
		return "".join(rez)
	else:
		return None

## compute some kind of 'consistency score' which approximates the number of haplotypes represented in a pseudoalignment
## TODO: come up with some theoretical justification for this
def _consistency_score(aln, maf, alphabet = "ACGT", eps = 0.00001):

	nseq = len(aln)
	if nseq:

		width = len(aln[0])
		freq = np.empty(shape = (width, len(alphabet)), type = float)
		freq.fill(eps)
		weights = np.zeros(shape = width, type = int)

		for i in range(0, width):

			seq = [ aln[j][i] for j in range(0, nseq) ]
			trimmed = [ x for x in seq if x in alphabet ]
			counts = collections.Counter(trimmed)
			total = sum(counts.values())

			k = 0
			for (nuc, freq) in counts.most_common():
				prob = float(freq/total)
				if prob > maf:
					freq[i][k] = prob
					k += 1
			weights[i] = total

		normfreq = freq/np.linalg.norm(freq, 1, 1).reshape(freq.shape[0], 1) ## column-wise frequencies
		colentropy = -1 * np.sum(normfreq * np.log2(normfreq), 1) ## column-wise Shannon information
		nalleles = np.power(len(alphabet), colentropy) ## 'expected number of alleles', column-wise
		adjweights = weights/np.sum(weights) * (nalleles - 1) ## adjusted weights: columns with less info count less; columns w/1 allle are downweighted

		return np.average(nalleles, weights = adjweights)

	else:
		return None
