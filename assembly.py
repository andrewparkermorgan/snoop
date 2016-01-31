#! /usr/bin/env python

## --- snoop/kmers.py --- ##
##	Date:		13 Nov 2014
##	Updated:	NA
##	Purpose:	Rudimentary de novo assembly tools for msBWTs; to be used with caution and maybe a grain of salt
##	Credits:	Most code taken verbatim from JMH's BWTPather.py

import os
import sys
import pickle
import collections
import argparse
import numpy as np

from snoop import util, dna
from MUS import TranscriptBuilder # actually does the assembly
from MUSCython import MultiStringBWTCython as MultiStringBWT

ASSEMBLE_FORWARD = True
ASSEMBLE_BACKWARD = False
ASSEMBLER_DEFAULTS = {
	'countK':		70,			# k-mer size to use in determining path weights
	'direction':		ASSEMBLE_FORWARD,	# assemble forward (right) or backward (left) from seed?
	'drawDollarTerminals':	False,			# follow paths ending in '$' to see where lots of reads end
	'isMerged':		False,			# are we using a merged msBWT? (JMH says this case is probably broken)
	'kmerSize':		70,			# k-mer size to use for graph construction
	'numNodes':		1000,			# maximum number of nodes to generate (NB: not same as final sequence length)
	'overloadThreshold':	3000,			# stop if path weight exceeds this threshold
	'pathThreshold':	5,			# stop if heaviest path weight doesn't meet this threshold
	'trackPairs':		False,			# keep track of read pairs
	'trackReads':		False,			# keep track of read IDs in each node (expensive)
	'useMemmap':		False }			# if False, load msBWT into memory first; up-front cost for runtime but downstream savings

def _build_graph(msbwt_path, seed, end_seeds = [], k = 70, count_k = 70, direction = ASSEMBLE_FORWARD,
		min_weight = 5, max_weight = 1.5e3, max_nodes = None,
		save = True, outprefix = "out", memmap = False,
		settings = {}, **kwargs):

	## make internal, fully-populated version of settings dictionary
	## override defaults with passed-in values
	_settings = collections.defaultdict()
	for e,v in ASSEMBLER_DEFAULTS.iteritems():
		if not e in settings:
			_settings[e] = v
		else:
			_settings[e] = settings[e]
	_settings["kmerSize"] = k
	_settings["pathThreshold"] = min_weight
	_settings["overloadThreshold"] = max_weight
	_settings["useMemmap"] = memmap
	_settings["countK"] = count_k
	_settings["direction"] = direction
	if max_nodes:
		_settings["numNodes"] = max_nodes

	sys.stderr.write(str(_settings) + "\n")

	## set up start and end seeds, trimming them to specifid k-mer size if needed
	seed = seed[0:k]
	for i in xrange(0, len(end_seeds)):
		end_seeds[i] = end_seeds[i][0:k]

	## set up output files
	fn = "{}.k{}.count{}.min{}.n{}".format(outprefix, k, _settings["countK"], min_weight, _settings["numNodes"])
	pickle_fn = fn + ".p"
	dot_fn = fn + ".dot.txt"

	## call the assembler
	logger = util._init_logger()
	tb = TranscriptBuilder.Assembler(msbwt_path, _settings, logger)
	tb.extendSeed(seed, end_seeds)

	## capture results
	(nodes, edges) = tb.getGraphResults()

	## dump to disk, if desired
	if save:
		pickle.dump((nodes, edges), open(pickle_fn, "wb"))

	return (nodes, edges)


def greedy_assemble(*args, **kwargs):

	print kwargs

	## set k-mer size
	k = 70
	if "k" in kwargs:
		k = kwargs["k"]
	kminus = k - 1

	## set maximum LENGTH of sequence to assemble
	maxlen = None
	if "max_length" in kwargs:
		maxlen = kwargs["max_length"]

	## construct assembly graph
	(nodes, edges) = _build_graph(*args, **kwargs)

	## initialize edge dictionary
	edgeDict = {}
	for x in xrange(0, len(nodes)):
		edgeDict[x] = []

	## fill in edge dictionary first
	for i, e in enumerate(edges):
		fromID = e[0]
		edgeDict[fromID].append(i)

	## traverse graph and assemble
	currentNode = 0
	seq = nodes[0][2]
	tailSeq = seq[ -kminus: ]
	visitedNodes = {0:True}
	while True:
		chosenEdge = -1
		chosenWeight = 0
		possibleEdges = edgeDict[currentNode]

		for eIndex in possibleEdges:
			toID = edges[eIndex][1]
			ew = edges[eIndex][2]

			if ew > chosenWeight:
				chosenEdge = eIndex
				chosenWeight = ew
			else:
				#previous choice had heavier weight
				pass

		if chosenEdge == -1 or visitedNodes.has_key(edges[chosenEdge][1]):
			#either there is no edge to choose or we already visited the chosen node at the end of the edge
			break
		else:

			#edge chosen to new node
			e = edges[chosenEdge]
			fromID = e[0]
			toID = e[1]

			print 'Chosen node = '+str(toID)

			n = nodes[toID]
			nodeSeq = n[2]

			if nodeSeq[0:kminus] != tailSeq:
				print tailSeq
				print nodeSeq[0:kminus]
				raise Exception("Different head and tail")
			else:
				visitedNodes[toID] = True
				currentNode = toID
				seq += nodeSeq[kminus:]
				tailSeq = nodeSeq[ -kminus: ]
				## check if we have extended sequence further than requested
				if maxlen:
					if len(seq) >= maxlen:
						seq = seq[0:maxlen]
						break

	return seq

## function to return command-line parser for parameters common to assembly tasks
def assembly_args():

	parser = argparse.ArgumentParser(add_help = False)
	parser.add_argument(	"-n","--maxnodes", type = int,
				default = 1000,
				help = "maximum size of assembly graph, in number of nodes [default: %(default)d]" )
	parser.add_argument(	"-L","--maxlength", type = int,
				default = None,
				help = "maximum length of final sequence to assemble, in bp [default: no limit]" )
	parser.add_argument(	"-x", "--minweight", type = int,
				default = 5,
				help = "minimum weight of path to keep assembling [default: %(default)d]" )
	parser.add_argument(	"-X", "--maxweight", type = int,
				default = 3000,
				help = "maximum weight of path to keep assembling [default: %(default)d]" )
	parser.add_argument(	"-k", "--kmer", type = int,
				default = 70,
				help = "k-mer size to use in building assembly graph [default: %(default)d]" )
	parser.add_argument(	"--count_kmer", type = int,
				help = "k-mer size to use in calculating path weights [default: same as --kmer]" )
	parser.add_argument(	"--prefix",
				help = "prefix for output file: if specified, assembly graph will be dumped as pickle" )
	parser.add_argument(	"--memmap", action = "store_true",
				help = "use memmapped files to initialize BWT, instead of loading into memory [default: False]" )

	return parser


def build_bridge(msbwt, seedKmer, targetKmer, tMin = 1, branchLim = 10, maxBranchLen = 250):
	"""
	Assemble the short "bridge" between two sequences expected to occur nearby on the template.
	@param msbwt - the MSBWT to use for searchs
	@param seedKmer - a k-mer to seed our bridging
	@param targetKmer - the target we are trying to bridge to
	@param tMin - the minimum k-count needed to consider the path
	@param branchLim - the maximum number of branches we will test
	@param maxBranchLen - the maximum length of a branch before giving up
	@return (ret, numBranched)
		ret - a list of bridges discovered; for most cases this is a list of length one
		numBranched - the number of branches we explored; if numBranched >= branchLim, this function was not 100% exhaustive
	"""

	#initialize to our input kmer
	ret = []
	possBridges = [ dna.unmask(dna.ungap(seedKmer)) ]
	targetKmer = dna.unmask(dna.ungap(targetKmer))
	kmerLen = len(seedKmer)

	#set up some easy values
	validChars = "ACGT"
	counts = np.zeros(dtype='<u8', shape=(len(validChars), ))
	numBranched = 0

	#print (seedKmer, targetKmer)

	#while we have things to explore, and we haven't explored too many, and we don't have a ridiculous number of possibilities
	while len(possBridges) > 0 and numBranched < branchLim:
		#get the bridge, the kmer, and the reverse kmer
		currBridge = possBridges.pop()
		numBranched += 1

		currKmer = currBridge[len(currBridge)-kmerLen:]
		revKmer = MultiStringBWT.reverseComplement(currKmer)

		#try to extend it on out
		while len(currBridge) < maxBranchLen:
			#get the counts for each possible extension
			for i, c in enumerate(validChars):
				counts[i] = msbwt.countOccurrencesOfSeq(currKmer+c)+msbwt.countOccurrencesOfSeq(dna.revcomp(c)+revKmer)

			#get the highest one
			maxPos = np.argmax(counts)
			maxSym = validChars[maxPos]

			#make sure the highest is high enough for us to consider it
			if counts[maxPos] >= tMin:
				if len(possBridges) < branchLim:
					#go through all the other possible extensions
					for i, c in enumerate(validChars):
						if i != maxPos and counts[i] >= tMin:
							#add the ones we aren't exploring right now if they're high enough
							possBridges.append(currBridge+c)

				#make sure the highest isn't too high
				#this extension meets our requirement so shift over to loop back around
				currBridge += maxSym
				currKmer = currKmer[1:]+maxSym
				revKmer = dna.revcomp(maxSym)+revKmer[0:len(revKmer)-1]
			else:
				#our BEST doesn't pass the threshold on this path, stop following
				break

			if currKmer.startswith(targetKmer):
				ret.append(currBridge)

	#return all our possibilities
	return (ret, numBranched)
