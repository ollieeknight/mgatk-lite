#!/usr/bin/python

###############################################
# Filter reads / clip .bam files
# Don't print() anything!!!! writing to STDOUT
##############################################

import sys
import pysam
import re

# Get command line arguments
bamfile, logfile, mtchr, proper_pair, NHmax, NMmax = sys.argv[1:7]

# Open the input BAM file
bam = pysam.AlignmentFile(bamfile, "rb")

# Modify the header to account for CRA v2 nonsense
header_string = str(bam.header)
if re.search(r"^@HD\s+SO:coordinate", header_string, re.MULTILINE):
    new_header = re.sub(r"^@HD\s+SO:coordinate", "@HD\tVN:1.5\tSO:coordinate", header_string, flags=re.MULTILINE)
else:
    new_header = header_string
out = pysam.AlignmentFile("-", "wb", text=new_header)

# Initialize counters
keepCount = 0
filtCount = 0

def filterReadTags(intags):
    '''
    Filter reads based on tag values
    '''
    nh = 0
    nm = 0
    for tag, value in intags:
        if tag == "NH":
            nh = value
        if tag == "NM":
            nm = value
    return nh <= int(NHmax) and nm <= int(NMmax)

def pairing(read):
    '''
    Check if reads are properly paired
    '''
    if proper_pair == "True":
        return read.is_proper_pair
    else:
        return True

def processRead(read):
    '''
    Process each read in the BAM file
    '''
    global keepCount, filtCount
    if filterReadTags(read.tags) and pairing(read):
        keepCount += 1
        out.write(read)
    else:
        filtCount += 1

# Process each read in the BAM file
for read in bam:
    processRead(read)

# Write the log file with counts of kept and filtered reads
with open(logfile, 'w') as outfile:
    outfile.write(f"Kept {keepCount} reads\n")
    outfile.write(f"Filtered {filtCount} reads\n")

# Close BAM files
bam.close()
out.close()