""" msBWT Tools<br>
&nbsp;&nbsp;&nbsp;&nbsp;Online tools for querying Multi-sring BWTs
"""

import markup
from markup import oneliner as element
import cgAdmin
import os, sys, time
import datetime
import glob
import locale
import numpy as np
locale.setlocale(locale.LC_ALL, 'en_US')

pysamPath = '/csbiosw/Ubuntu-10.04-x86_64/python_packages/lib/python/pysam-0.5-py2.6-linux-x86_64.egg/'
sys.path.append(pysamPath)

#import MultiStringBWT
#import pyximport
#import numpy as np
#pyximport.install(setup_args=dict(include_dirs=np.get_include()))
from MUSCython import MultiStringBWTCython as MultiStringBWT

MSBWTdir = "/playpen"

def indexPage(form):
	""" Main query page """
	panel = markup.page()
	panel.div(style="padding:50px 50px;")
	'''
	for bwtDir in glob.glob('/csbiodata/CEGS3x3BWT/*'):
		panel.h3(bwtDir)
		msbwt = MultiStringBWT.loadBWT(bwtDir)
		break
	'''
	panel.h3("Select Dataset:")
	available = sorted(glob.glob("%s/*/*msbwt.npy" % MSBWTdir))
	panel.div(style="padding:0px 0px 40px 120px;")
	panel.form(action="", method="POST", enctype="multipart/form-data")
	for i, dataset in enumerate(available):
		end = dataset.rfind('/')
		start = dataset.rfind('/', 0, end-1) + 1
		shorten = dataset[start:end]
		panel.input(type="checkbox", name="dataset", value=shorten)
		panel.add(shorten)
		panel.br()
	panel.div.close()

	panel.label("Search Pattern:")
	panel.input(type="text", name="pattern", size="100")
	panel.input(type="hidden", name="target", value="msAllele.Search")
	panel.input(type="submit", name="submit", value="Submit")
	panel.form.close()
	panel.div.close()
	return panel

WatsonComp = { 'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', '$':'$' }

def revComp(pat):
	return ''.join([WatsonComp[c] for c in reversed(pat)])
'''
def conSeq(seqList):
	N = len(seqList[0])
	posCount = [dict() for i in xrange(N)]
	for seq in seqList:
		for i, c in enumerate(seq):
			posCount[i][c] = posCount[i].get(c, 0) + 1
	result = ''
	for i in xrange(N):
		pairs = [(posCount[i][c], c) for c in posCount[i].iterkeys() if (c != '.' and c != '$')]
		count, base = (N, '.') if (len(pairs) == 0) else max(pairs)
		result += base.upper()
	return result
'''
def readCmp(read1, read2):
	i = len(read1)//2
	offset = 1
	while (read1[i] == read2[i]):
		i += offset
		offset = -(offset + 1)
		if ((i < 0) or (i >= len(read1))):
			return 0
	return ord(read1[i]) - ord(read2[i])

def SearchResponse(form):
	panel = markup.page()

	panel.script(type="text/javascript")
	panel.add("""
		function getSelectedText() {
			var hidden, submit;
			var selectedText=(window.getSelection ? window.getSelection() : document.getSelection ? document.getSelection() : document.selection.createRange().text);
			if (selectedText == "") {
				alert("You must select a subsequence");
				return false;
			} else {
				document.forms["SearchSelected"]["pattern"].value = selectedText;
			}
		}
	""")
	panel.script.close()

	panel.div(style="padding:50px 50px;")
	datasets = form.getvalue("dataset")
	if (datasets == None):
		panel.h3("ERROR: No datasets selected.")
		panel.div(align="center", style="padding: 30px 30px;")
		panel.input(type="button", value="New Search", onClick='self.location="./?run=msAllele"')
		panel.div.close()
		panel.div.close()
		return panel
	
	if isinstance(datasets, str):
		datasets = [datasets]

	pattern = form.getvalue("pattern")
	if (pattern == None):
		panel.h3("ERROR: No search pattern specified")
		panel.div(align="center", style="padding: 30px 30px;")
		panel.input(type="button", value="New Search", onClick='self.location="./?run=msAllele"')
		panel.div.close()
		panel.div.close()
		return panel
	pattern = pattern.upper()

	for dataset in datasets:
		panel.h3(dataset)
		bwtDirName = "%s/%s" % (MSBWTdir, dataset)
		filestat = os.stat(bwtDirName+"/comp_msbwt.npy")
		filesize = locale.format("%d", filestat.st_size, grouping=True)
		bwt = MultiStringBWT.loadBWT(bwtDirName)
		stringCount = locale.format("%d", bwt.getSymbolCount(0), grouping=True)
		baseCount = locale.format("%d", bwt.getTotalSize(), grouping=True)
		bitsPerBase = (8.0*filestat.st_size)/bwt.getTotalSize()
		panel.strong("%s: %s strings with %s bases and index size of %s bytes (%3.2f bits per base)<br />" % (dataset, stringCount, baseCount, filesize, bitsPerBase))
		panel.strong("Target: %s<br />" % (pattern))
		
		lo1, hi1 = bwt.findIndicesOfStr(pattern)
		lo2, hi2 = bwt.findIndicesOfStr(revComp(pattern))
		count = hi1 - lo1 + hi2 - lo2
		if (count > 10000):
			panel.add("Found %d times (%d forward, %d reverse-complemented)<br /><br />" % (count, hi1-lo1, hi2-lo2))
			panel.span("Too much data!", style="font-size: 180%;")
		elif count > 0:
			
			'''
			l = len(pattern)
			bufferLen = 101
			fixedSize = 2*bufferLen-l
			readlist = []
			
			for i in xrange(lo1, hi1):
				#pass
				suffix = bwt.recoverString(i)
				suffLen = len(suffix)
				end = suffix.find('$')
				beforePattern = suffLen-end-1
				read = ('.'*(bufferLen-l-beforePattern)+
						suffix[end+1:].lower()+
						suffix[:l]+
						suffix[l:end+1].lower())
				read += '.'*(fixedSize-len(read))
				readlist.append(read)
			
			for i in xrange(lo2, hi2):
				suffix = revComp(bwt.recoverString(i))
				suffLen = len(suffix)
				end = suffix.find('$')
				beforePattern = suffLen-end-l
				read = ('.'*(bufferLen-l-beforePattern)+
						suffix[end:-l].lower()+
						suffix[-l:]+
						suffix[:end].lower())
				read += '.'*(fixedSize-len(read))
				readlist.append(read)
			'''
			panel.add("Found %d times (%d forward, %d reverse-complemented)<br /><br />" % (count, hi1-lo1, hi2-lo2))
			panel.div(style="font-size:10px; font-family: monospace;")
			#margin = len(suffix)-l
			l = len(pattern)
			margin = 101-l

			haps = extractHaplotypes(bwt, pattern)
			if len(haps) > 0:
				consensusMain = (sorted(haps, key=lambda x: x[2][0][1], reverse=True))[0][0]

			panel.table(border='1')
			panel.tr()
			panel.th('Consensus')
			panel.th('Exact matches')
			panel.tr.close()
			extrasList = []
			for consensus, readlist, counts in sorted(haps, key=lambda x: x[2][0][1], reverse=True):
				#panel.strong('%s<span style="color: green;">%s</span>%s<br />' % (consensus[:margin].upper(), consensus[margin:margin+l].upper(), consensus[margin+l:].upper()))
				if counts[0][1] > 0:
					panel.tr()
					panel.td()
					panel.strong()
					output = ""
					for i, base in enumerate(consensus):
						if i == margin:
							output += '<span style="color: green;">'
						elif i == margin+l:
							output += '</span>'

						if(base != '$') and (base != '.') and (consensus[i] != '.') and (base.upper() != consensusMain[i].upper()):
							output += '<span style="background-color:yellow;">%s</span>' % base.upper()
						else:
							output += base.upper()
					panel.add(output)
					panel.strong.close()
					panel.td.close()
					panel.td(str(counts[0][1]))
					panel.tr.close()
				
				for read in readlist[counts[0][1]:]:
					extrasList.append(read)
			
			if len(extrasList) > 0:
				consensus, dummyVar = conSeq(extrasList)
				panel.tr()
				panel.td()
				panel.strong()
				output = ""
				for i, base in enumerate(consensus):
					if i == margin:
						output += '<span style="color: green;">'
					elif i == margin+l:
						output += '</span>'
					
					if(base != '$') and (base != '.') and (consensus[i] != '.') and (base.upper() != consensusMain[i].upper()):
						output += '<span style="background-color:yellow;">%s</span>' % base.upper()
					else:
						output += base.upper()
				panel.add(output)
				panel.strong.close()
				panel.td.close()
				panel.td('0')
				panel.tr.close()
			panel.table.close()
			
			for consensus, readlist, counts in sorted(haps, key=lambda x: x[2][0][1], reverse=True):
				#consensus = conSeq(readlist)
				#panel.add(consensus)
				#readlist.sort(cmp=readCmp)
				if counts[0][1] == 0:
					continue
				
				read = "."*margin + "*"*l + '.'*margin
				panel.add(read)
				panel.br()
				for read in readlist[0:counts[0][1]]:
					color = "red" if (read.find('$') > read.find(pattern)) else "blue"
					output = ""
					for i, base in enumerate(read):
						if (i == margin):
							output += '<span style="color: %s;">' % color
						elif (i == margin+l):
							output += '</span>'
						if (base != '$') and (base != '.') and (consensus[i] != '.') and (base.upper() != consensus[i].upper()):
							output += '<span style="background-color:yellow;">%s</span>' % base
						else:
							output += base
					output += '<br />'
					panel.add(output)
				panel.strong('%s<span style="color: green;">%s</span>%s<br />' % (consensus[:margin], consensus[margin:margin+l], consensus[margin+l:]))
				panel.br()
				panel.br()
		
			if len(extrasList) > 0:
				consensus, dummyVar = conSeq(extrasList)
				panel.add(consensus)
				extrasList.sort(cmp=readCmp)
				read = "."*margin + "*"*l + '.'*margin
				panel.add(read)
				panel.br()
				for read in extrasList:
					color = "red" if (read.find('$') > read.find(pattern)) else "blue"
					output = ""
					for i, base in enumerate(read):
						if (i == margin):
							output += '<span style="color: %s;">' % color
						elif (i == margin+l):
							output += '</span>'
						if (base != '$') and (base != '.') and (consensus[i] != '.') and (base.upper() != consensus[i].upper()):
							output += '<span style="background-color:yellow;">%s</span>' % base
						else:
							output += base
					output += '<br />'
					panel.add(output)
				panel.strong('%s<span style="color: green;">%s</span>%s<br />' % (consensus[:margin], consensus[margin:margin+l], consensus[margin+l:]))
				panel.br()
		
			panel.div.close()
		else:
			panel.add("Pattern not found<br /><br />")
	panel.form(action="", name="SearchSelected", method="POST", enctype="multipart/form-data", onsubmit='return getSelectedText()')
	panel.div(align="center", style="padding: 30px 30px;")
	panel.input(type="submit", name="submit", value="Search Selected")
	panel.input(type="button", value="New Search", onClick='self.location="./?run=msAllele"')
	for dataset in datasets:
		panel.input(type="hidden", name="dataset", value=dataset)
	panel.input(type="hidden", name="pattern", value=pattern)
	panel.input(type="hidden", name="target", value="msAllele.Search")
	panel.div.close()
	panel.form.close()
	panel.div.close()
	return panel

def extractHaplotypes(bwt, kmer):
	forwardIndices = bwt.findIndicesOfStr(kmer)
	revComp = MultiStringBWT.reverseComplement(kmer)
	reverseIndices = bwt.findIndicesOfStr(revComp)
	
	readLen = 101
	patternLen = len(kmer)
	totalBuffLen = 2*readLen-patternLen
	
	modifiedSeqs = []
	for i in xrange(forwardIndices[0], forwardIndices[1]):
		readSeq = bwt.recoverString(i)
		dollarPos = readSeq.find('$')
		
		#calcualte how many tailing '.' we need first, then construct the string from that info
		afterPattern = readLen-dollarPos-1
		modSeq = ('.'*(readLen-patternLen-afterPattern)+
				  readSeq[dollarPos+1:].lower()+
				  readSeq[0:patternLen]+
				  readSeq[patternLen:dollarPos+1].lower()+
				  '.'*(afterPattern))
		modifiedSeqs.append(modSeq)
	
	for i in xrange(reverseIndices[0], reverseIndices[1]):
		revCompSeq = bwt.recoverString(i)
		readSeq = MultiStringBWT.reverseComplement(revCompSeq)
		dollarPos = readSeq.find('$')
		
		#beforePattern = readLen-dollarPos
		afterPattern = readLen-dollarPos-patternLen
		modSeq = ('.'*(readLen-patternLen-afterPattern)+
				  readSeq[dollarPos:-patternLen].lower()+
				  readSeq[-patternLen:]+
				  readSeq[0:dollarPos].lower()+
				  '.'*(afterPattern))
		modifiedSeqs.append(modSeq)

	finishedHaps = []

	previousConsensus = 'A'*totalBuffLen
	currentConsensus, currentScorer = conSeq(modifiedSeqs)
	currSeqs = modifiedSeqs

	while len(currSeqs) > 0 and compareShiftedSeqs(previousConsensus, currentConsensus) > 0:
		nextSeqs = []
		consensusSeqs = []
		
		#we will fill in consensus Seqs downstream
		finishedHaps.append((currentConsensus, consensusSeqs, []))
		
		#first get all exact matching reads
		for seq in currSeqs:
			if compareShiftedSeqs(seq, currentConsensus) == 0:
				consensusSeqs.append(seq)
			else:
				nextSeqs.append(seq)
		
		finishedHaps[-1][2].append((0, len(consensusSeqs)))
		
		#update these things
		previousConsensus = currentConsensus
		currSeqs = nextSeqs
		currentConsensus, currentScorer = conSeq(currSeqs)
		
		#check if the next consensus is identical
		acceptedScore = 1
		while len(currSeqs) > 0 and compareShiftedSeqs(currentConsensus, previousConsensus) == 0:
			#print 'triggered', acceptedScore
			nextNextSeqs = []
			minScore = 0xFFFFFFFFFFFFFFFF
			for seq in nextSeqs:
				calcScore = scoreShiftedSeqs(seq, currentConsensus, currentScorer)
				if calcScore < minScore and calcScore > acceptedScore:
					minScore = calcScore
				
				if calcScore <= acceptedScore:
					consensusSeqs.append(seq)
				else:
					nextNextSeqs.append(seq)
			finishedHaps[-1][2].append((acceptedScore, len(nextSeqs)-len(nextNextSeqs)))
			
			nextSeqs = nextNextSeqs
			currSeqs = nextSeqs
			currentConsensus, currentScorer = conSeq(currSeqs)
			
			#acceptedScore += 1
			acceptedScore = minScore

	for seq in currSeqs:
		consensusSeqs.append(seq)

	return finishedHaps

def compareShiftedSeqs(seq1, seq2):
	'''
		compares two "shifted" sequences which have been offset such that they share a pattern in the middle GUARANTEED
		'''
	careSymbols = 'acgtACGT'
	l = len(seq1)
	diff = 0
	for x in xrange(0, l):
		if ((seq1[x] in careSymbols) and
			(seq2[x] in careSymbols) and
			(seq1[x] != seq2[x])):
			diff += 1
	return diff

def scoreShiftedSeqs(seq1, seq2, counter):
	careSymbols = 'acgtACGT'
	l = len(seq1)
	score = 0
	for x in xrange(0, l):
		if ((seq1[x] in careSymbols) and
			(seq2[x] in careSymbols) and
			(seq1[x] != seq2[x])):
			if counter[x] == 0:
				raise Exception
			score += counter[x]
	return score

def combineShiftedSeqs(seq1, seq2):
	careSymbols = 'acgtACGT'
	l = len(seq1)
	ret = ''
	for x in xrange(0, l):
		if (seq1[x] in careSymbols):
			ret += seq1[x]
		elif (seq2[x] in careSymbols):
			ret += seq2[x]
		else:
			ret += '.'
	return ret

def conSeq(seqList):
	if len(seqList) == 0:
		return '', []
	N = len(seqList[0])
	counter = np.zeros(dtype='<u8', shape=(N, ))
	posCount = [dict() for i in xrange(N)]
	for seq in seqList:
		for i, c in enumerate(seq):
			posCount[i][c] = posCount[i].get(c, 0) + 1
	result = ''
	for i in xrange(N):
		pairs = [(posCount[i][c], c) for c in posCount[i].iterkeys() if (c != '.' and c != '$')]
		count, base = (N, '.') if (len(pairs) == 0) else max(pairs)
		for c, b in pairs:
			if b != base:
				counter[i] += c
		result += base
	
	return (result, counter)

if __name__ == '__main__':
	print sys.argv[0]
