#!/usr/bin/python

###################################################
# Summarizes the total number of reads per position / strand
###################################################

import sys
import pysam
import numpy as np

# Parse command line arguments
bam_file = sys.argv[1]
barcodes_file = sys.argv[2]
out_pre = sys.argv[3]
max_bp = int(sys.argv[4])
base_qual = float(sys.argv[5])
fasta_file = sys.argv[6]
alignment_quality = float(sys.argv[7])
barcode_tag = sys.argv[8]

# Import barcodes
with open(barcodes_file) as barcode_file_handle:
	bcs = [x.strip() for x in barcode_file_handle.readlines()]

bam_input = pysam.AlignmentFile(bam_file, "rb")
dna_letters = ['A', 'C', 'G', 'T']

# Initialize coverage array
ca = np.zeros((len(bcs), max_bp, 4, 2), dtype=int)

def getBarcode(intags):
	'''
	Parse out the barcode per-read
	'''
	for tg in intags:
		if barcode_tag == tg[0]:
			return tg[1]
	return "NA"

# Process each read in the BAM file
for read in bam_input:
	s_idx = 1 if read.is_reverse else 0
	seq = read.seq
	quality = read.query_qualities
	align_qual_read = read.mapping_quality
	cell_barcode = getBarcode(read.tags)
	
	if cell_barcode != "NA":
		c_idx = bcs.index(cell_barcode)
		for q_idx, p_idx in read.get_aligned_pairs(True):
			if q_idx is not None and p_idx is not None and align_qual_read > alignment_quality:
				if quality[q_idx] > base_qual and seq[q_idx] in dna_letters:
					l_idx = dna_letters.index(seq[q_idx])
					ca[c_idx, p_idx, l_idx, s_idx] += 1

def writeSparseMatrixLetter(letter, letter_idx):
	'''
	Write the slice of the matrix associated with a specific letter
	'''
	out_file_fn = f"{out_pre}.{letter}.txt"
	with open(out_file_fn, "w") as file_handle_fn:
		for cell_idx, cell_name in enumerate(bcs):
			fw_vec = ca[cell_idx, :, letter_idx, 0].ravel()
			rev_vec = ca[cell_idx, :, letter_idx, 1].ravel()
			for i in range(max_bp):
				if fw_vec[i] > 0 or rev_vec[i] > 0:
					file_handle_fn.write(f"{i+1},{cell_name},{fw_vec[i]},{rev_vec[i]}\n")

# Write sparse matrices for each DNA letter
for idx, letter in enumerate(dna_letters):
	writeSparseMatrixLetter(letter, idx)

# Export per-base coverage and depth
out_file_depth = out_pre.replace("/temp/sparse_matrices/", "/logs/depth/") + "_depth.txt"
out_file_coverage = f"{out_pre}.coverage.txt"

with open(out_file_coverage, "w") as file_handle_cov, open(out_file_depth, "w") as file_handle_depth:
	for cell_idx, cell_name in enumerate(bcs):
		cov_vec = np.sum(ca[cell_idx, :, :, :], axis=(1, 2)).tolist()
		depth = round(sum(cov_vec) / len(cov_vec), 2)
		for i in range(max_bp):
			if cov_vec[i] > 0:
				file_handle_cov.write(f"{i+1},{cell_name},{cov_vec[i]}\n")
		file_handle_depth.write(f"{cell_name}\t{depth}\n")
