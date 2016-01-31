#! /usr/bin/env python

## --- snoop/haplotype.py --- ##
##	Date: 12 Oct 2014
##	Updated: 31 Oct 2014
##	Purpose: Functions for extracting haplotypes from queries against a msBWT

import os
import sys
import collections
import itertools
import copy
import re
import numpy as np
import scipy.stats as sp

import MUSCython.MultiStringBWTCython as ms

from snoop import io, dna

class PseudoalignedRow:

	def __init__(self, seq = "", alphabet = "ACGT"):
		self.seq = seq
		self._alphabet = alphabet
		self.nonmissing = sum([ x in alphabet for x in seq ])

	def __len__(self):
		return len(self.seq)

	def __eq__(self, other):

		if not self._check_compatible(other):
			raise ArugmentError
		else:
			flag = True
			for i in range(0, len(self)):
				if (self.seq[i] in self._alphabet and other.seq[i] in self._alphabet):
					flag = flag and (self.seq[i] == other.seq[i])
			return flag

	def __str__(self):
		return "".join(self.seq) + "\n"

	## check that this row and another object are 'compatible': both of right class, both of same length
	def _check_compatible(self, other):
		flag  = True
		if not isinstance(other, PseudoalignedRow):
			flag = False
		if len(self) != len(other):
			flag = False
		return flag

	## merge two aligned rows, assuming 'other' < 'self' (in the sense that the first nonmissing char in self occurs before that in other)
	def merge(self, other):

		if not self._check_compatible(other):
			raise ArugmentError

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

		if not self._check_compatible(other):
			raise ArugmentError

		diffs = []
		for i in range(0, len(self)):
			if (self.seq[i] in self._alphabet and other.seq[i] in self._alphabet):
				if self.seq[i] != other.seq[i]:
					diffs.append(i)

		if len(diffs):
			return diffs
		else:
			return None

	## JMH 15 Oct 2014
	## score this row against another
	def score(self, other, counter):

		if not self._check_compatible(other):
			raise ArugmentError

		score = 0
		for x in xrange(0, len(self)):
			if ((self.seq[x] in self._alphabet) and
				(other.seq[x] in self._alphabet) and
				(self.seq[x] != other.seq[x])):
				if counter[x] == 0:
					raise Exception
				score += counter[x]
		return score


## JMH 15 Oct 2014
## generate consensus sequence from a list of pseudoaligned sequences
def consensus(seqList, alphabet = "ACGT"):

	if len(seqList) == 0:
		return '', []

	N = len(seqList[0])
	counter = np.zeros(dtype='<u8', shape=(N, ))
	posCount = [dict() for i in xrange(N)]

	for seq in seqList:
		for i, c in enumerate(seq):
			posCount[i][c] = posCount[i].get(c, 0) + 1
	result = ""

	for i in xrange(N):
		pairs = [(posCount[i][c], c) for c in posCount[i].iterkeys() if c in alphabet ]
		count, base = (N, "-") if (len(pairs) == 0) else max(pairs)
		for c, b in pairs:
			if b != base:
				counter[i] += c
		result += base

	return ( PseudoalignedRow(result), counter )

## JMH 15 Oct 2014
## extract haplotypes present in a pseudoalignment
# def extract_haplotypes(aln):
#
# 	nseq = len(aln)
# 	if not nseq:
# 		return None
#
# 	width = len(aln[0])
#
# 	finishedHaps = []
#
# 	previousConsensus = PseudoalignedRow("A"*width)
# 	currentConsensus, currentScorer = consensus(aln)
# 	currSeqs = aln
#
# 	while len(currSeqs) > 0 and previousConsensus.diffs(currentConsensus):
#
# 		nextSeqs = []
# 		consensusSeqs = []
#
# 		#we will fill in consensus Seqs downstream
# 		finishedHaps.append((currentConsensus, consensusSeqs, []))
#
# 		#first get all exact matching reads
# 		for seq in currSeqs:
# 			if not currentConsensus.diffs(PseudoalignedRow(seq)):
# 				consensusSeqs.append(seq)
# 			else:
# 				nextSeqs.append(seq)
#
# 		finishedHaps[-1][2].append((0, len(consensusSeqs)))
#
# 		#update these things
# 		previousConsensus = currentConsensus
# 		currSeqs = nextSeqs
# 		currentConsensus, currentScorer = consensus(currSeqs)
#
# 		#check if the next consensus is identical
# 		acceptedScore = 1
# 		while len(currSeqs) > 0 and not previousConsensus.diffs(currentConsensus):
# 			#print 'triggered', acceptedScore
# 			nextNextSeqs = []
# 			minScore = 0xFFFFFFFFFFFFFFFF
# 			for seq in nextSeqs:
# 				calcScore = currentConsensus.score(PseudoalignedRow(seq), currentScorer)
# 				if calcScore < minScore and calcScore > acceptedScore:
# 					minScore = calcScore
#
# 				if calcScore <= acceptedScore:
# 					consensusSeqs.append(seq)
# 				else:
# 					nextNextSeqs.append(seq)
# 			finishedHaps[-1][2].append((acceptedScore, len(nextSeqs)-len(nextNextSeqs)))
#
# 			nextSeqs = nextNextSeqs
# 			currSeqs = nextSeqs
# 			currentConsensus, currentScorer = consensus(currSeqs)
#
# 			#acceptedScore += 1
# 			acceptedScore = minScore
#
# 	for seq in currSeqs:
# 		consensusSeqs.append(seq)
#
# 	return finishedHaps
#
