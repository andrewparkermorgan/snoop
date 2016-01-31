#! /usr/bin/env python

## --- snoop/rle.py --- ##
##	Date: 16 April 2015
##	Updated: NA
##	Purpose: Simple run-length encoding scheme

import sys
import itertools
import operator

class rle:

	## see <https://docs.python.org/3/library/itertools.html#itertools.accumulate>
	def _accumulate(self, x, func = operator.add):
		it = iter(x)
		total = next(it)
		yield total
		for element in it:
			total = func(total, element)
			yield total

	def __init__(self, x = None):
		self._original = None
		self._iterator = None
		self._pairs = []
		self.values = []
		self.lengths = []
		self._breaks = []
		self.starts = []
		self.ends = []
		self.bounds = []
		if x is not None:
			self._original = x
			self._iterator = itertools.groupby(x)
			self._pairs = [ (k, len(list(g))) for k,g in self._iterator ]
			self.values = zip(*self._pairs)[0]
			self.lengths = zip(*self._pairs)[1]
			self._breaks = list(self._accumulate(self.lengths))
			self.starts = [0] + [ i for i in self._breaks[:-1] ]
			self.ends = list(self._breaks)
			self.bounds = zip(self.starts, self.ends)

	## special case: create a binary rle from positions of 1s
	@classmethod
	def from_ones(cls, x, maxlen):
		if max(x) >= maxlen:
			raise ArgumentError("Position of item exceeds specified sequence length.")
		y = [0]*maxlen
		for i in x:
			y[i] = 1
		return rle(y)

	def __len__(self):
		return len(self._pairs)

	def __str__(self):
		return self._pairs.__str__()

	def reconstitute(self):
		return [ str(x)*y for x,y in self._pairs ]

	## filter runs based on their value
	def filter_values(self, fn = lambda x: x):

		bounds = []
		for i in range(0, len(self)):
			if fn(self.values[i]):
				bounds.append( self.bounds[i] )

		return bounds

	## filter runs based on their value
	def filter_lengths(self, fn = lambda x: x > 0):

		bounds = []
		for i in range(0, len(self)):
			if fn(self.lengths[i]):
				bounds.append( self.bounds[i] )

		return bounds

	## get max element in list, but return 'middle' in case of ties
	def _argmax_middle(self, x):

			m = max(x)
			xx = [ i for i in range(0,len(x)) if x[i] == m ]
			ii = int(len(xx)/2)
			return xx[ii]

	## for integers <x>, find the smallest set of intervals of width <w> (step size <s> between them) which covers them all
	## (with thanks to Leonard McMillan)
	def greedy_set_cover(self, w = 30, step = 1, buff = 5, xmax = 1e9):

		## get boundaries of filtered runs
		x = [ s for s,e in self.filter_values() ]

		## set bounds of range containing candidate intervals
		x = set(x)
		lb = max(0, min(x)-buff)
		ub = min(xmax, max(x)+buff)

		## make candidate intervals; initialize empty set of accepted ones
		#sys.stderr.write(str(range(lb, ub, step)) + "\n")
		candidates = [ range(i,i+w) for i in range(lb, ub, step) if any (z in x for z in range(i,i+w)) ]
		accepted = []
		#sys.stderr.write("Input is {}\n".format(str(x)))
		#sys.stderr.write(str(candidates) + "\n")
		#sys.stderr.write("\t{} candidate sets\n".format(len(candidates)))

		## keep going while there are elements left to be covered
		while len(x) and len(candidates):

			#sys.stderr.write("\t{} points left to cover with {} candidates\n".format(len(x), len(candidates)))

			## list of lists of positions in <x> which each candidate hits
			hits = []
			## how many such hits for each candidate
			nhits = []
			for j in range(0, len(candidates)):
				x_hits = x.intersection(candidates[j])
				hits.append( x_hits )
				nhits.append(len(x_hits))

			## take candidate which hit the most (ignore ties)
			best = self._argmax_middle(nhits)
			
			## remove all elements in <x> covered by best candidate
			x = x - set(candidates[best])

			## add this candidate to the list of accepted intervals
			this_best = candidates.pop(best)
			accepted.append( (min(this_best), max(this_best)) )

		return sorted(accepted, key = operator.itemgetter(0))