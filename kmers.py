#! /usr/bin/env python

## --- snoop/kmers.py --- ##
##	Date: 12 Oct 2014
##	Updated: NA
##	Purpose: Utility functions for enumeration of k-mers in BWTs, plus associated counting and stats

import os
import sys
import numpy as np

## wrapper class for Matt's generator for the union of k-mers present in a list of msBWT objects
class KmerIterator:

	def __init__(self, bwts, k, threshold = 0):
		self._bwts = bwts
		self.k = k
		self.threshold = threshold

	def emit(self):

		nbwts = len(self._bwts)
		vcLen = 6
		currentK = 0
		kmerStarts = np.zeros(dtype = "<u8", shape = (self.k+1, nbwts))
		kmerEnds = np.zeros(dtype = "<u8", shape = (self.k+1, nbwts))
		for x in range(0, nbwts):
			kmerEnds[0, x] = self._bwts[x].getTotalSize()

		status = np.zeros(dtype = "<u1", shape = (self.k+1, ))

		while currentK > 0 or status[0] < vcLen:
			if currentK == self.k:
				for x in range(0, nbwts):
					if kmerStarts[currentK, x] != kmerEnds[currentK, x]:
						seq = self._bwts[x].recoverString(kmerStarts[currentK, x])[0:self.k]
				yield (seq, kmerEnds[currentK]-kmerStarts[currentK])
				currentK -= 1
			elif status[currentK] == 0 or status[currentK] == 4:
				status[currentK] += 1
			elif status[currentK] == 6:
				currentK -= 1
			else:
				# we need to write the correct value into the kmerStarts and kmerEnds at currentK+1
				foundVals = False
				for x in range(0, nbwts):
					fmStart = self._bwts[x].getFullFMAtIndex(kmerStarts[currentK, x])
					fmEnd = self._bwts[x].getFullFMAtIndex(kmerEnds[currentK, x])
					kmerStarts[currentK+1, x] = fmStart[status[currentK]]
					kmerEnds[currentK+1, x] = fmEnd[status[currentK]]

					if fmEnd[status[currentK]] - fmStart[status[currentK]] > self.threshold:
						foundVals = True

				# update status[currentK]
				status[currentK] += 1

				if foundVals:
					# increment currentK if we found something
					currentK += 1
					status[currentK] = 0

## wrapper class to enmerate k-mers and provide some statistics on their relative frequency
class KmerProfiler:

	def __init__(self, bwts, k, threshold = 1, sampling_rate = 1.0, maxiter = float("Inf")):

		self._bwts = bwts
		self._maxiter = maxiter
		self.k = k
		self.sampling_rate = sampling_rate
		self._threshold = threshold
		self._kmers = KmerIterator(self._bwts, self.k, self._threshold)

		self.sizes = np.zeros(len(self._bwts), dtype = "u8")
		self.weights = np.zeros(len(self._bwts), dtype = "f8")
		for i in range(0, len(self._bwts)):
			self.sizes[i] = self._bwts[i].getTotalSize()
		self.weights = self.sizes/float(sum(self.sizes))

		self._total_jsd = 0
		self._j = 0


	def reset(self):

		self._kmers = KmerIterator(self._bwts, self.k, self._threshold)
		self._j = 0
		self._total_jsd = 0

	## iterate using the k-mer generator
	## return tuple of (kmer, [list of counts], element JSD, cumulative JSD)
	def profile(self):

		for (kmer, counts) in self._kmers.emit():

			self._j += 1
			if self._j > self._maxiter:
				break

			(kmer, counts) = self._kmers.__next()__
			## downsampling by randomly including or excluding this k-mer for each sapmle... why??
			# flag = np.ones(n)
			flag = 1
			# if not (sampling_rate == 1.0):
				# flag = np.array(np.random.ranf(n) < sampling_rate, dtype = int)

			## perform downsampling
			if np.random.ranf() > self.sampling_rate:
				# print "skipping kmer '{}'".format(kmer)
				continue

			## compute element-wise Jensen-Shannon divergence (JSD)
			probs = [ (float(max(counts[i], threshold))/sizes[i]) for i in range(0, len(self._bwts)) ]
			mean_prob = np.average(probs, weights = self.weights)
			logged = [ np.log(probs[i]/mean_prob)*probs[i] for i in range(0, len(self._msbwt)) ]
			jsd = np.average(logged, weights = self.weights)
			self._total_jsd += jsd * flag

			yield (kmer, counts, jsd, total_jsd)
