#! /usr/bin/env python

'''
Created on Nov 12, 2014
@author: Matt
'''

import logging
import numpy as np
import os
import pickle
import sys

from MUS import TranscriptBuilder

def initLogger():
    '''
    This code taken from Shunping's Lapels for initializing a logger
    '''
    global logger
    logger = logging.getLogger('root')
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s", "%Y-%m-%d %H:%M:%S")
    ch.setFormatter(formatter)
    logger.addHandler(ch)

def testTranscriptBuilder():
    #the bwt to load
    bwtDir = '/netscr/apm/msbwt/WSBEiJ'
    
    #the base extension for your output
    baseFN = '/nescr/apm/chr2/r2d2/work'
    
    #the seed you care about
    seed = 'GTTAATGTAGCTTAATAACAAAGCAAAGCACTGAAAATGCTTAGATGGATAATTGTATCCCATAAACACA'
    
    #if you have a particular k-mer you want to stop at, add it to this list
    endSeeds = []
    
    #Set these as well, make sure kmerSize = the length of your kmer
    kmerSize = 70
    
    #all counts will be based on k-mers of this size
    countK = 70
    
    #minimum threshold for a "path" to be consider real, 20 may be too high/low depending on the coverage
    threshold = 20
    
    #the maximum number of nodes to generate
    numNodes = 1000
    
    '''
    lots of info store here
    countK - see previous
    direction - True if you want to builder forward from your k-mer (backwards might not actually work right now anyways)
    drawDollarTerminals - if True, it will follow '$' paths to show large terminals for reads
    isMerged - only if a dataset is a merged dataset, also likely broken
    kmerSize - see previous
    numNodes - see previous
    overloadThreshold - if you get a count higher than this, the node will terminate
    pathThreshold - see previous
    trackPairs - if True, tracks read pairs if the information is available
    trackReads - if True, tracks where reads are in the nodes (causes high overhead)
    useMemmap - if False, the BWT is loaded into memory prior to extension, causes overhead at the start but leads to faster searches laster
    '''
    settingsDict = {'countK':countK,
                    'direction':True,
                    'drawDollarTerminals':False,
                    'isMerged':False,
                    'kmerSize':kmerSize,
                    'numNodes':numNodes,
                    'overloadThreshold':3000,
                    'pathThreshold':threshold,
                    'trackPairs':False,
                    'trackReads':False,
                    'useMemmap':False}
    
    seed = seed[0:kmerSize]
    for x in xrange(0, len(endSeeds)):
        endSeeds[x] = endSeeds[x][0:kmerSize]
    
    #this is where the output will be stored
    pickleFN = baseFN+str(kmerSize)+'mer_thresh'+str(threshold)+'_'+str(numNodes)+'n_'+str(countK)+'c.traversal.p'
    outputDotFile = baseFN+str(kmerSize)+'mer_thresh'+str(threshold)+'_'+str(numNodes)+'n_'+str(countK)+'c.dot.txt'
    
    #call the assembler
    tb = TranscriptBuilder.Assembler(bwtDir, settingsDict, logger)
    tb.extendSeed(seed, endSeeds)
    
    #retData is nodes, retEdges is edges
    (retData, retEdges) = tb.getGraphResults()
    
    #I don't recommend changing these for your purposes
    drawPairs = False
    drawInversions = False
    colorScheme = 'coverage'
    scales = [1]
    
    #saves the data to a pickle file and a dotty file
    saveTranscriptBuilderResults(retData, retEdges, pickleFN, outputDotFile, colorScheme, drawPairs, drawInversions, scales)
    
    #outFastaFN = 'this is a path to a fasta FN to write'
    #label = 'this is the label for the path in the fasta file'
    #createFastaFromTranscriptBuilder(pickleFN, outFastaFN, label)
    
def createFastaFromTranscriptBuilder(pickleFN, outFastaFN, seqLabel=""):
    '''
    This starts from node 0 (which contains the kmer seed) and greedily extends it on the most significant path
    '''
    #load the pickle
    (nodeData, edgeData) = pickle.load(open(pickleFN, 'rb'))
        
    edgeDict = {}
    for x in xrange(0, len(nodeData)):
        edgeDict[x] = []
    
    #build an edge dictionary first
    for i, e in enumerate(edgeData):
        fromID = e[0]
        edgeDict[fromID].append(i)
    
    currentNode = 0
    seq = nodeData[0][2]
    tailSeq = seq[-69:]
    visitedNodes = {0:True}
    
    while True:
        chosenEdge = -1
        chosenWeight = 0
        possibleEdges = edgeDict[currentNode]
        
        for eIndex in possibleEdges:
            toID = edgeData[eIndex][1]
            ew = edgeData[eIndex][2]
            
            if ew > chosenWeight:
                chosenEdge = eIndex
                chosenWeight = ew
            else:
                #previous choice had heavier weight
                pass
        
        if chosenEdge == -1 or visitedNodes.has_key(edgeData[chosenEdge][1]):
            #either there is no edge to choose or we already visited the chosen node at the end of the edge
            break
        else:
            
            #edge chosen to new node
            e = edgeData[chosenEdge]
            fromID = e[0]
            toID = e[1]
            
            print 'Chosen node = '+str(toID)
            
            n = nodeData[toID]
            nodeSeq = n[2]
            
            if nodeSeq[0:69] != tailSeq:
                print tailSeq
                print nodeSeq[0:69]
                raise Exception("Different head and tail")
            else:
                visitedNodes[toID] = True
                currentNode = toID
                seq += nodeSeq[69:]
                tailSeq = nodeSeq[-69:]
    
    print 'Writing sequence of length '+str(len(seq))+'...'
    
    fp = open(outFastaFN, 'w+')
    fp.write('>'+seqLabel+'\n')
    while len(seq) > 0:
        fp.write(seq[0:40]+'\n')
        seq = seq[40:]
    fp.close()
    
    print 'Fasta written.'

def saveTranscriptBuilderResults(nodeData, edgeData, pickleFN, dotFN, colorScheme, drawPairs, drawInversions, scales):
    #simple pickle save
    pickle.dump((nodeData, edgeData), open(pickleFN, 'wb'))
    
    #everything after is DOT parsing
    try:
        os.remove(dotFN)
    except:
        pass
    
    fp = open(dotFN, 'w+')
    fp.write("digraph g {\ngraph [\nrankdir=\"UD\"\n];\nnode [\nfontsize=16\nshape=record\n];\nedge [\n];\n")
    
    lowBlueThresh = 30
    highRedThresh = 500
    
    for k, blockTup in enumerate(nodeData):
        (nodeID, execID, seq, pileups, endReason, minDist, inversionSet, readSet, pairedNodes, sourceDistro) = blockTup
        fp.write("node"+str(nodeID)+" [\n")
        
        if colorScheme == 'merged':
            print sourceDistro
            for i in xrange(0, max(sourceDistro)+1):
                sourceDistro[i] = int(sourceDistro.get(i, 0)*scales[i])
            print sourceDistro
        
        #joinArray = [str(nodeID)+'('+str(execID)+')']
        joinArray = [str(nodeID)]
        if len(pileups) == 0:
            joinArray.append('0')
            joinArray.append('--')
        else:
            joinArray.append(str(len(pileups)))
            joinArray.append('{0:.1f}'.format(np.median(pileups)))
            #joinArray.append(str(int(np.median(pileups))))
        joinArray.append(str(minDist))
        
        joinArray.append(str(endReason))
        if drawInversions:
            joinArray.append(str(list(inversionSet)))
        joinArray.append(str(len(readSet)))
        
        if drawPairs:
            joinArray.append(str(pairedNodes))
        if colorScheme == 'merged':
            joinArray.append(str(sourceDistro))
        
        fp.write("label=\""+"|".join(joinArray)+"\"\n")
        fp.write("style=filled\n")
        
        #determine the color
        if colorScheme == 'coverage':
            if endReason != None and endReason[0:8] == 'END_SEED':
                fp.write('fillcolor=\"#55ff55\"\n')
            elif endReason != None and endReason == 'OVERLOAD':
                fp.write('fillcolor=\"#ffa500\"\n')
            elif len(pileups) > 0 and np.median(pileups) <= lowBlueThresh:
                fp.write('fillcolor=\"#ccccff\"\n')
            elif len(pileups) > 0 and np.median(pileups) >= highRedThresh:
                fp.write('fillcolor=\"#ffcccc\"\n')
        elif colorScheme == 'merged':
            tc = sum(sourceDistro.values())
            '''
            colorArrays = [(0, 255, 0),
                           (255, 0, 0),
                           (0, 0, 255)]
            '''
            
            colorArrays = [(0, 254, 0),
                           (0, 0, 254),
                           (127, 127, 127)]
            
            colorSubs = []
            for ctup in colorArrays:
                colorSubs.append([255-v for v in ctup])
            
            color = [255, 255, 255]
            for i in xrange(0, 3):
                ratio = float(sourceDistro.get(i, 0.0))/tc
                for j in xrange(0, 3):
                    color[j] -= int(ratio*colorSubs[i][j])
                
                #print i, colorSubs[i], ratio, color
            
            whiteText = True
            colorStr = '#'
            for i, c in enumerate(color):
                if c < 16:
                    colorStr += '0'
                colorStr += (hex(c))[2:]
                
                if c > 180 and (i == 0 or i == 1):
                    whiteText = False
            
            print color, colorStr, tc, sourceDistro
            fp.write('fillcolor=\"'+colorStr+'\"\n')
            if whiteText:
                fp.write('fontcolor=\"#ffffff\"\n')
            
        fp.write("];\n")
        
    for edge in edgeData:
        fp.write("node"+str(edge[0])+" -> node"+str(edge[1])+" [label=\""+str(edge[3])+"\", style="+edge[4]+"];\n")
    
    fp.write("}")
    fp.close()

if __name__ == '__main__':
    	initLogger()
	testTranscriptBuilder()
    
