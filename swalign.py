#! /usr/bin/env python

from __future__ import print_function

import os
import sys

import numpy as np

## constants for alignment operations
OP_I = 0
OP_D = 1
OP_M = 2

def align(ref, query, match = 10, mismatch = 0, gap_open = -2, gap_extend = -1):
	"""Perform global alignment of two sequences of similar length."""

	def score(one, two):
		if one == two:
			return match
		else:
			return mismatch

	m,n = len(query), len(ref)
	dp = np.zeros((m+1,n+1), dtype = np.int32)
	tb = np.zeros((m+1,n+1), dtype = np.int32)

	tb[:,0] = OP_I
	tb[0,:] = OP_D

	## fill DP matrix
	for j in xrange(1, n+1):
		for i in xrange(1, m+1):

			# match/mismatch
			m_score = dp[ i-1,j-1 ] + score( query[i-1], ref[j-1] )

			# insertion (gap in reference)
			i_score = 0
			if not dp[ i-1,j ] == 0:
				if tb[ i-1,j ] == OP_I:
					i_score = dp[ i-1,j ] + gap_extend
				else:
					i_score = dp[ i-1,j ] + gap_open

			# deletion (gap in query)
			d_score = 0
			if not dp[ i,j-1 ] == 0:
				if tb[ i,j-1 ] == OP_D:
					d_score = dp[ i,j-1 ] + gap_extend
				else:
					d_score = dp[ i,j-1 ] + gap_open

			# update scoring and traceback matrices
			scores = np.array([ i_score, d_score, m_score ])
			this_cell = np.max(scores)
			this_op = np.argmax(scores)
			dp[ i,j ] = this_cell
			tb[ i,j ] = this_op

	## do traceback
	row = m
	col = n
	query_aln, ref_aln, ops = "","",""
	while True:
		# quit if at upper left corner
		if row == 0 and col == 0:
			break
		score = dp[row,col]
		op = tb[row,col]
		if op == OP_M:
			#print "match"
			query_aln += query[row-1]
			ref_aln += ref[col-1]
			row -= 1
			col -= 1
			if ref[col] == query[row]:
				ops += "."
			else:
				ops += "*"
		elif op == OP_I:
			#print "ins"
			query_aln += query[row-1]
			ref_aln += "-"
			row -= 1
			ops += "i"
		elif op == OP_D:
			#print "del"
			query_aln += "-"
			ref_aln += ref[col-1]
			col -= 1
			ops += "d"

	return (ref_aln[::-1], query_aln[::-1], ops[::-1])

def reconcile(ref, alt):
	"""Given pairwise alignment, identify variant sites and make VCF-style POS/REF/ALT"""
	i = 0
	all_vars = []
	while i < len(ref):
		if ref[i] == alt[i]:
			# ref and alt match
			# are there active variants?
			if len(all_vars):
				# yes; finish it and emit
				yield all_vars.pop()
		else:
			# ref and alt don't match
			# are there active variants?
			if not len(all_vars):
				# no; start one
				if ref[i] == "-" or alt[i] == "-":
					# indel; backtrack one base to get REF allele
					all_vars.append( [ i-1, ref[i-1], alt[i-1] ] )
				else:
					all_vars.append( [i, "", ""] )
			# now operate on the active variant
			curr_var = all_vars[-1:][0]
			if ref[i] == "-":
				# insertion; extend ALT only
				curr_var[2] += alt[i]
			elif alt[i] == "-":
				# deletion; extend REF only
				curr_var[1] += ref[i]
			else:
				# [SM]NV; extend both alleles
				curr_var[1] += ref[i]
				curr_var[2] += alt[i]

		# walk to next base
		i += 1

	# done; emit any active variants
	if len(all_vars):
		yield all_vars.pop()


if __name__ == "__main__":

	ref, alt, stars = align("AAACTGGTTGTCC","AAACTGGCCTTGTAC")
	print(ref)
	print(alt)
	for v in reconcile(ref, alt):
		print(v)
