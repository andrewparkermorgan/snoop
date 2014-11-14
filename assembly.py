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
import numpy as np

from snoop import util
from MUS import TranscriptBuilder # actually does the assembly

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
		end_seeds[i] = end_seeds[i][0:kmerSize]

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