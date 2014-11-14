#! /usr/bin/env python

##	--- snoop/dna.py --- ##
##	Date: 18 June 2014
##	Updated: 14 Oct 2014
##	Purpose: utility functions for DNA sequences

import re
import zlib

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

## undo lowercase repeat-masking on DNA sequence
def unmask(dna):
	dnastr = "".join(dna)
	return dnastr.upper()

## mask lowercase repeats to Ns (or whatever)
def mask(dna, repeatchar = "N"):
	dnastr = "".join(dna)
	return re.sub(r"[a-z]", "N", dnastr)

def _as_seq_object(dna, alphabet = IUPAC.ambiguous_dna):

	if not isinstance(dna, Seq):
		dna = Seq(dna, alphabet)

	return dna

## reverse complement, including on-the-fly conversion to Bio.Seq
def revcomp(dna, alphabet = IUPAC.unambiguous_dna):

	dna = _as_seq_object(dna, alphabet)
	return str(dna.reverse_complement())

## check for ambiguous characters
def has_ambiguities(dna, alphabet = IUPAC.ambiguous_dna):

	dna = _as_seq_object(dna, alphabet)
	flag = any([ nt in dna for nt in alphabet.letters[4:] ])
	return flag

## strip gaps etc.
def ungap(dna, alphabet = IUPAC.unambiguous_dna):

	dna = _as_seq_object(dna, alphabet)
	ungapped = "".join([ nt for nt in dna if nt in alphabet.letters ])
	return ungapped

## quick-and-dirty estimation of sequence complexity, using compression ratio with Lempel-Ziv as proxy
def complexity(dna, alphabet = IUPAC.unambiguous_dna):

	dna = _as_seq_object(dna, alphabet)
	cplx = float(len(zlib.compress(str(dna))))/len(str(dna))
	return cplx
