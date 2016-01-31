#! /usr/bin/env python

## --- snoop/fa2bwt.py --- ##
##	Date:		14 Nov 2014
##	Udpated:	NA
##	Purpose:	Generate a msBWT from one or more fasta files, assuming sequences have non-uniform length

import os
import sys
import argparse

from snoop import util, io
from MUSCython import MultimergeCython as mm

parser = argparse.ArgumentParser(description = "Tool for doing simplistic/naive graph-based assembly using a msBWT, starting from a specified seed.")
parser.add_argument(	"-M","--msbwt", type = io.writeable_dir,
			required = True,
			help = "new directory in which to build msBWT components; best to provide it as absolute path" )
parser.add_argument(	"-t", "--threads", type = int,
			default = 1,
			help = "number of threads to use during msBWT construction [default: %(default)d]" )
parser.add_argument(	"fastas", metavar = "FASTA", nargs = "+",
			help = "path(s) to fasta files of sequences" )
args = parser.parse_args()

sys.stderr.write("\n=== msBWT construction from fasta ===\n")
sys.stderr.write("\tinput files: " + str(args.fastas) + "\n")
sys.stderr.write("\toutput directory: <" + args.msbwt + ">\n\n")

mm.createMSBWTFromFasta(args.fastas, args.msbwt, args.threads, 0, util._init_logger())