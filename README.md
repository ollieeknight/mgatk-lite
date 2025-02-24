# mgatk lite: Mitochondrial Genome Analysis Toolkit

**mgatk** is a comprehensive toolkit designed for the analysis of mitochondrial genomes, originally developed by [Caleb Lareau](https://github.com/caleblareau). The original provides a suite of tools to process and analyze mitochondrial DNA (mtDNA) from single-cell sequencing data, particularly from 10x Genomics platforms.  

*This* version, however, has now removed most flexible functionality, with the purpose of optimizing and streamlining the `tenx` mode to work with 10x Genomics scATAC/ASAP-seq data.

Some of the options may not work as in our lab we generally run the program as:

```sh
mgatk tenx \
    -i ${project_outs}/${library}/outs/possorted_bam.bam \
    -n output \
    -o ${project_outs}/${library}/mgatk \
    -c 1 \
    -bt CB \
    -b ${project_outs}/${library}/outs/filtered_peak_bc_matrix/barcodes.tsv \
    --skip-R \
    --remove-snakemake
```

## Features

- **Single-cell mtDNA Analysis**: Process and analyze mtDNA from single-cell ATAC-seq
- **Variant Calling**: Identify and quantify mtDNA variants across single cells.
- **Quality Control**: Perform extensive quality control on mtDNA sequencing data.

## Installation

To install **mgatk**, use the following command:

```sh
conda create -y -n mgatk openjdk r-base r-data.table r-matrix bioconductor-genomicranges bioconductor-summarizedexperiment 
pip install pulp==2.7.0 matplotlib
git clone https://github.com/ollieeknight/mgatk
cd mgatk
pip install .
```

## Test usage

```sh
git clone https://github.com/ollieeknight/mgatk
cd mgatk/tests/

mgatk tenx -i barcode/test_barcode.bam -n bc1 -o bc1dmem -bt CB -b barcode/test_barcodes.txt -c 2

2024-12-03 11:08:47,969 - INFO - mgatk version 1.0.0
2024-12-03 11:08:48,307 - INFO - Configuration for analysis:
2024-12-03 11:08:48,307 - INFO - Mode                      : tenx
2024-12-03 11:08:48,307 - INFO - Input                     : barcode/test_barcode.bam
2024-12-03 11:08:48,307 - INFO - Output                    : bc1dmem
2024-12-03 11:08:48,307 - INFO - Name                      : bc1
2024-12-03 11:08:48,307 - INFO - Mitochondrial genome      : rCRS
2024-12-03 11:08:48,307 - INFO - Number of cores           : 2
2024-12-03 11:08:48,307 - INFO - Barcode tag               : CB
2024-12-03 11:08:48,308 - INFO - Barcodes                  : barcode/test_barcodes.txt
2024-12-03 11:08:48,308 - INFO - Min barcode reads         : 1000
2024-12-03 11:08:48,308 - INFO - NHmax                     : 1
2024-12-03 11:08:48,308 - INFO - NMmax                     : 4
2024-12-03 11:08:48,308 - INFO - Remove duplicates         : False
2024-12-03 11:08:48,308 - INFO - UMI barcode               : XX
2024-12-03 11:08:48,308 - INFO - Handle overlap            : False
2024-12-03 11:08:48,308 - INFO - Low coverage threshold    : 10
2024-12-03 11:08:48,308 - INFO - Max java memory           : 8000m
2024-12-03 11:08:48,308 - INFO - Proper pairs              : False
2024-12-03 11:08:48,308 - INFO - Base quality              : 0
2024-12-03 11:08:48,308 - INFO - Alignment quality         : 0
2024-12-03 11:08:48,308 - INFO - Emit base qualities       : False
2024-12-03 11:08:48,308 - INFO - Number of samples         : 0
2024-12-03 11:08:48,308 - INFO - Keep samples              : ALL
2024-12-03 11:08:48,308 - INFO - Ignore samples            : NONE
2024-12-03 11:08:48,308 - INFO - Keep temp files           : False
2024-12-03 11:08:48,308 - INFO - Keep QC .bam files        : False
2024-12-03 11:08:48,308 - INFO - Skip R                    : False
2024-12-03 11:08:48,308 - INFO - Snake stdout              : False
2024-12-03 11:08:48,308 - INFO - Remove .snakemake         : False
2024-12-03 11:08:48,330 - INFO - Using fasta file: bc1dmem/fasta/chrM.fasta
2024-12-03 11:08:50,618 - INFO - Executing Snakemake command: snakemake --snakefile /home/olive/bin/mgatk/mgatk/bin/snake/Snakefile.tenx --cores 2 --config cfp="bc1dmem/.internal/parseltongue/snake.scatter.yaml" > bc1dmem/logs/snakemake.log 2>&1
2024-12-03 11:09:22,163 - INFO - Successfully created final output files. Analysis complete.
```

### Options

- `--input, -i`: Input .bam file from 10x single cell ATAC library.
- `--output, -o`: Output directory for analysis. Default is `mgatk`.
- `--name, -n`: Prefix for project name. Default is `mgatk`.
- `--mito-genome, -g`: Mitochondrial genome configuration. Default is `rCRS`.
- `--ncores, -c`: Number of cores to run the main job in parallel. Default is 1.
- `--barcode-tag, -bt`: Read tag to separate single cells.
- `--barcodes, -b`: Path to a file containing known barcodes.
- `--min-barcode-reads, -mb`: Minimum number of mitochondrial reads for a barcode to be genotyped. Default is 1000.
- `--NHmax`: Maximum number of read alignments allowed. Default is 1.
- `--NMmax`: Maximum number of paired mismatches allowed. Default is 4.
- `--remove-duplicates, -rd`: Remove duplicate reads.
- `--umi-barcode, -ub`: Read tag to specify the UMI tag when removing duplicates.
- `--handle-overlap, -ho`: Only count each base in the overlap region between a pair of reads once.
- `--low-coverage-threshold, -lc`: Variant count for each cell will be ignored below this threshold. Default is 10.
- `--max-javamem, -jm`: Maximum memory for Java for running duplicate removal per core. Default is 8000m.
- `--proper-pairs, -pp`: Require reads to be properly paired.
- `--base-qual, -q`: Minimum base quality for inclusion in the genotype count. Default is 0.
- `--alignment-quality, -aq`: Minimum alignment quality to include read in genotype. Default is 0.
- `--emit-base-qualities, -eb`: Output mean base quality per alt allele as part of the final output.
- `--nsamples, -ns`: Number of samples/cells to be processed per iteration. Default is 0 (all).
- `--keep-samples, -k`: Comma-separated list of sample names to keep. Default is ALL.
- `--ignore-samples, -x`: Comma-separated list of sample names to ignore. Default is NONE.
- `--keep-temp-files, -z`: Keep all intermediate files.
- `--keep-qc-bams, -qc`: Keep the quality-controlled bams after processing.
- `--skip-R, -sr`: Generate plain-text only output.
- `--snake-stdout, -so`: Write Snakemake log to stdout.
- `--remove-snakemake, -rs`: Delete the .snakemake directory once successfully run.

For more details, visit the [mgatk GitHub repository](https://github.com/caleblareau/mgatk).

## Contact

The [**mgatk**](https://github.com/caleblareau/mgatk) package was developed and is maintained by [Caleb Lareau](https://www.mskcc.org/research/ski/labs/caleb-lareau), but this repo is maintained by [Oliver Knight](https://immunologie.charite.de/metas/person/person/address_detail/oliver_knight/).