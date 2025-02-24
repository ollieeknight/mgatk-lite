#!/usr/bin/python

###############################################
# Filter reads / clip .bam files
# Don't print() anything!!!! writing to STDOUT
##############################################

import sys
import pysam

# Get command line arguments
bamfile, logfile, mtchr, proper_pair, NHmax, NMmax = sys.argv[1:7]

# Open the input BAM file
bam = pysam.AlignmentFile(bamfile, "rb")

# Modify the header to account for CRA v2 nonsense
new_header = str(bam.header).replace("@HD\tSO:coordinate", "@HD\tVN:1.5\tSO:coordinate")
out = pysam.AlignmentFile("-", "wb", text=new_header)

# Initialize counters
keepCount = 0
filtCount = 0

def filterReadTags(intags):
	'''
	Checks for aligner-specific read tags and filters
	'''
	for tg in intags:
		if (('NH' == tg[0] and int(tg[1]) > int(NHmax)) or
			(('NM' == tg[0] or 'nM' == tg[0]) and int(tg[1]) > int(NMmax))):
			return False
	return True

def pairing(read):
	'''
	Check if read is paired, properly paired, etc.
	'''
	return proper_pair != "True" or read.is_proper_pair

def processRead(read):
	'''
	Process each read and decide whether to keep or filter it
	'''
	global keepCount, filtCount
	if filterReadTags(read.tags) and read.reference_name == mtchr and pairing(read):
		keepCount += 1
		out.write(read)
	else:
		filtCount += 1

# Process each read in the BAM file
for read in bam:
	processRead(read)

# Write the log file with counts of kept and filtered reads
with open(logfile, 'w') as outfile:
	outfile.write(f"Kept {keepCount}\nRemoved {filtCount}\n")
