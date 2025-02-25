import itertools
import time
import shutil
import re
import os
import sys
import subprocess
import pysam
import filecmp
import math
from typing import List, Tuple, Dict

def string_hamming_distance(str1: str, str2: str) -> int:
	"""
	Fast hamming distance over 2 strings known to be of same length.
	In information theory, the Hamming distance between two strings of equal 
	length is the number of positions at which the corresponding symbols 
	are different.
	eg "karolin" and "kathrin" is 3.
	"""
	return sum(x != y for x, y in zip(str1, str2))

def rev_comp(seq: str) -> str:
	"""
	Fast Reverse Complement
	"""  
	tbl = str.maketrans('ATCGN', 'TAGCN')
	return seq.translate(tbl)[::-1]

def gettime() -> str: 
	"""
	Matches `date` in Linux
	"""
	return time.strftime("%a %b %d %X %Z %Y: ")

def findIdx(list1: List, list2: List) -> List[int]:
	"""
	Return the indices of list1 in list2
	"""
	return [i for i, x in enumerate(list1) if x in list2]

def check_R_dependencies(required_packages: List[str]) -> None:
	"""
	Determines whether or not R packages are properly installed
	"""
	R_path = shutil.which("R")
	if not R_path:
		sys.exit("ERROR: R is not installed or not in PATH.")
	
	installed_packages = os.popen(f"{R_path} -e \"installed.packages()\" | awk '{{print $1}}' | sort | uniq").read().strip().split("\n")
	missing_packages = set(required_packages) - set(installed_packages)
	if missing_packages:
		sys.exit(f"ERROR: cannot find the following R package(s): {', '.join(missing_packages)}\n"
				 "Install them in your R console and then try rerunning mgatk (but there may be other missing dependencies).")

def check_software_dependencies(tools: List[str]) -> None:
	"""
	Determines whether or not required software tools are properly installed
	"""
	for tool in tools:
		if not shutil.which(tool):
			sys.exit(f"ERROR: cannot find {tool} in environment; add it to user PATH environment")

def check_pip_dependencies(required_packages: List[str]) -> None:
	"""
	Determines whether or not pip packages are properly installed.
	Supports version checking with syntax package==version
	"""
	import pkg_resources
	missing_packages = []
	
	for package in required_packages:
		if "==" in package:
			# Handle versioned packages
			name, version = package.split("==")
			try:
				installed_version = pkg_resources.get_distribution(name).version
				if installed_version != version:
					missing_packages.append(f"{name} (requires version {version}, found {installed_version})")
			except pkg_resources.DistributionNotFound:
				missing_packages.append(package)
		else:
			# Handle packages without version requirements
			try:
				__import__(package)
			except ImportError:
				missing_packages.append(package)
	
	if missing_packages:
		sys.exit(f"ERROR: cannot find the following pip package(s): {', '.join(missing_packages)}\n"
				 "Install them using pip and then try rerunning the script.")

def parse_fasta(filename: str) -> Dict[str, str]:
	"""
	Imports specified .fasta file
	"""
	sequences = {}
	with open(filename) as f:
		for line in f:
			if line.startswith('>'):
				name = line[1:].strip()
				sequences[name] = ''
			else:
				sequences[name] += line.strip()
	return sequences

def verify_bai(bamfile: str) -> None:
	"""
	Function that indexes bam file from input if missing
	"""
	bai_file = bamfile + ".bai"
	if not os.path.exists(bai_file):
		pysam.index(bamfile)

def split_chunk_file(one_barcode_file: str, script_dir: str, input: str, bcbd: str, barcode_tag: str, mito_chr: str, umi_barcode: str) -> None:
	chunk_bam_py = os.path.join(script_dir, "bin/python/chunk_barcoded_bam.py")
	pycall = ['python', chunk_bam_py, input, bcbd, barcode_tag, one_barcode_file, mito_chr, umi_barcode]
	subprocess.run(pycall, check=True)

def verify_sample_mitobam(bam: str, mito_chr: str, mito_length: int) -> bool:
	idxs = pysam.idxstats(bam).split("\n")
	nReads = 0
	bam_length = 0
	for i in idxs:
		if i.split("\t")[0] == mito_chr:
			bam_length = int(i.split("\t")[1])
			nReads = int(i.split("\t")[2])
	if mito_length == -9:
		mito_length = bam_length
	return bam_length == mito_length and nReads > 0

def handle_fasta_inference(mito_genome: str, supported_genomes: List[str], script_dir: str, of: str, write_files: bool = True) -> Tuple[str, str, int]:
	"""
	Determines what's going on with the mitochondrial genome
	based on user input / existing data
	"""
	if mito_genome in supported_genomes:
		fastaf = os.path.join(script_dir, "bin/anno/fasta", f"{mito_genome}.fasta")
	elif os.path.exists(mito_genome):
		fastaf = mito_genome
	else:
		sys.exit(f'ERROR: Could not find file {mito_genome}; QUITTING')

	fasta = parse_fasta(fastaf)
	if len(fasta) != 1:
		sys.exit('ERROR: .fasta file has multiple chromosomes; supply file with only 1; QUITTING')

	mito_genome, mito_seq = list(fasta.items())[0]
	mito_length = len(mito_seq)

	if write_files:
		make_folder(os.path.join(of, "fasta"))
		make_folder(os.path.join(of, "final"))

	newfastaf = os.path.join(of, "fasta", f"{mito_genome}.fasta")
	writeFA = not os.path.exists(newfastaf) or not filecmp.cmp(fastaf, newfastaf, shallow=False)

	if writeFA and write_files:
		shutil.copyfile(fastaf, newfastaf)
		fastaf = newfastaf
		pysam.faidx(fastaf)
		with open(os.path.join(of, "final", f"{mito_genome}_refAllele.txt"), 'w') as f:
			for b, base in enumerate(mito_seq, 1):
				f.write(f"{b}\t{base}\n")
	return fastaf, mito_genome, mito_length

def make_folder(folder: str) -> None:
	"""
	Function to only make a given folder if it does not already exist
	"""
	os.makedirs(folder, exist_ok=True)

def file_len(fname: str) -> int:
	with open(fname) as f:
		return sum(1 for _ in f)

def split_barcodes_file(barcode_file: str, nsamples: int, output: str) -> List[str]:
	"""
	Function to split barcode file into smaller files if needed
	"""
	n_samples_observed = file_len(barcode_file)
	if n_samples_observed < nsamples or nsamples == 0:
		return [barcode_file]

	total_files = math.ceil(n_samples_observed / nsamples)
	lines_per_file = nsamples
	full_output_folder = os.path.join(output, "temp", "barcode_files")
	make_folder(full_output_folder)

	barcodes_files = []
	with open(barcode_file) as bigfile:
		for lineno, line in enumerate(bigfile):
			if lineno % lines_per_file == 0:
				if lineno > 0:
					smallfile.close()
				small_filename = os.path.join(full_output_folder, f"barcodes.{len(barcodes_files) + 1}.txt")
				smallfile = open(small_filename, "w")
				barcodes_files.append(small_filename)
			smallfile.write(line)
		if smallfile:
			smallfile.close()
	return barcodes_files