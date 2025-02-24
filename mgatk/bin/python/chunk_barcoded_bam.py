#!/usr/bin/python

import sys
import pysam
import os

# Get command line arguments
bamfile, outfolder, barcodeTag, bcfile, mtchr, umitag = sys.argv[1:7]

# Extract base names
base = os.path.basename(bamfile)
basename = os.path.basename(os.path.splitext(bcfile)[0])

def getBarcode(read, tag_get):
	'''
	Parse out the barcode per-read
	'''
	try:
		return read.get_tag(tag_get)
	except Exception:
		return "AA"

# Read in the barcodes
with open(bcfile) as barcode_file_handle:
	bc = set(line.strip() for line in barcode_file_handle)  # from list to set, so query becomes O(1)

# Open BAM file for reading and create output BAM file
bam = pysam.AlignmentFile(bamfile, "rb")
outname = os.path.join(outfolder, f"{basename}.bam")
out = pysam.AlignmentFile(outname, "wb", template=bam)

# Generate a list of faux UMIs
bases = "ACGT"
fauxdon = [a + b + c + d for a in bases for b in bases for c in bases for d in bases]

# Filter reads that match the set of possible barcodes for this sample
try:
	for read in bam.fetch(mtchr, multiple_iterators=False):
		barcode_id = getBarcode(read, barcodeTag)
		
		if barcode_id in bc:
			# Check for true UMI
			umi_id = getBarcode(read, umitag) if umitag != "XX" else ""
			
			# Create a fake UMI
			if barcode_id[-1].isnumeric():
				split_two = barcode_id.split("-")
				faux_umi = split_two[0] + umi_id + fauxdon[int(split_two[1]) - 1]
			else:
				faux_umi = barcode_id + umi_id
			
			# Add the fake UMI to the read tags
			read.tags = read.tags + [("MU", faux_umi)]
			out.write(read)

except OSError:
	print('Finished parsing bam')

# Close BAM files and create index for the output BAM
bam.close()
out.close()
pysam.index(outname)
