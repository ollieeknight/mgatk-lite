import click
import os
import sys
import shutil
import glob
import logging
import math
import pysam
from pkg_resources import get_distribution
from subprocess import call, check_call
from .mgatkHelp import *
from ruamel.yaml import YAML
from ruamel.yaml.scalarstring import SingleQuotedScalarString as sqs
from multiprocessing import Pool
from itertools import repeat

# Setup logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@click.command()
@click.version_option()
@click.option('--input', '-i', default=".", required=True, help='Input .bam file from ASAP- or DOGMA-seq 10x single cell ATAC library')
@click.option('--output', '-o', default="mgatk", help='Output directory for analysis. Default = mgatk')
@click.option('--name', '-n', default="mgatk", help='Prefix for project name. Default = mgatk')
@click.option('--mito-genome', '-g', default="rCRS", required=True, help='Mitochondrial genome configuration. Default = rCRS.')
@click.option('--ncores', '-c', default=1, help='Number of cores to run the main job in parallel.')
@click.option('--barcode-tag', '-bt', default="X", help='Read tag (generally two letters) to separate single cells')
@click.option('--barcodes', '-b', default="", help='Path to a file containing known barcodes, must be .tsv, so gunzip if from a paired DOGMA-seq mapping')
@click.option('--min-barcode-reads', '-mb', default=1000, help='Minimum number of mitochondrial reads for a barcode to be genotyped. Default = 1000.')
@click.option('--NHmax', default=1, help='Maximum number of read alignments allowed as governed by the NH flag. Default = 1.')
@click.option('--NMmax', default=4, help='Maximum number of paired mismatches allowed represented by the NM/nM tags. Default = 4.')
@click.option('--remove-duplicates', '-rd', is_flag=True, help='Remove duplicate (presumably PCR) reads')
@click.option('--umi-barcode', '-ub', default="XX", help='Read tag (generally two letters) to specify the UMI tag when removing duplicates for genotyping.')
@click.option('--handle-overlap', '-ho', is_flag=True, help='Only count each base in the overlap region between a pair of reads once')
@click.option('--low-coverage-threshold', '-lc', default=10, help='Variant count for each cell will be ignored below this when calculating VMR')
@click.option('--max-javamem', '-jm', default="8000m", help='Maximum memory for java for running duplicate removal per core. Default = 8000m.')
@click.option('--proper-pairs', '-pp', is_flag=True, help='Require reads to be properly paired.')
@click.option('--base-qual', '-q', default=0, help='Minimum base quality for inclusion in the genotype count. Default = 0.')
@click.option('--alignment-quality', '-aq', default=0, help='Minimum alignment quality to include read in genotype. Default = 0.')
@click.option('--emit-base-qualities', '-eb', is_flag=True, help='Output mean base quality per alt allele as part of the final output.')
@click.option('--nsamples', '-ns', default=0, help='The number of samples / cells to be processed per iteration; Default = 0, all. Supply an integer')
@click.option('--keep-samples', '-k', default="ALL", help='Comma separated list of sample names to keep; ALL (special string) by default. Sample refers to basename of .bam file')
@click.option('--ignore-samples', '-x', default="NONE", help='Comma separated list of sample names to ignore; NONE (special string) by default. Sample refers to basename of .bam file')
@click.option('--keep-temp-files', '-z', is_flag=True, help='Add this flag to keep all intermediate files.')
@click.option('--keep-qc-bams', '-qc', is_flag=True, help='Add this flag to keep the quality-controlled bams after processing.')
@click.option('--skip-R', '-sr', is_flag=True, help='Generate plain-text only output. Otherwise, this generates a .rds object that can be immediately read into R for downstream analysis.')
@click.option('--snake-stdout', '-so', is_flag=True, help='Write snakemake log to stdout rather than a file. May be necessary for certain HPC environments.')
@click.option('--remove-snakemake', '-rs', is_flag=True, help='Delete the .snakemake directory once successfully run.')

def main(input, output, name, mito_genome, ncores, barcode_tag, barcodes, min_barcode_reads, nhmax, nmmax, remove_duplicates, umi_barcode, handle_overlap, low_coverage_threshold, max_javamem, proper_pairs, base_qual, alignment_quality, emit_base_qualities, nsamples, keep_samples, ignore_samples, keep_temp_files, keep_qc_bams, skip_r, snake_stdout, remove_snakemake):
    """
    mgatk: a mitochondrial genome analysis toolkit.
    See https://github.com/ollieeknight/mgatk for more details.
    """
    script_dir = os.path.dirname(os.path.realpath(__file__))
    cwd = os.getcwd()
    __version__ = get_distribution('mgatk').version
    logger.info("mgatk version %s" % __version__)

    # Verify dependencies
    logger.info("Verifying software dependencies...")
    check_software_dependencies(['R', 'java'])
    logger.info("Verifying R dependencies...")
    check_R_dependencies(["data.table", "SummarizedExperiment", "GenomicRanges", "Matrix"])
    logger.info("Verifying pip dependencies...")
    check_pip_dependencies(['matplotlib', 'pulp'])

    ncores = int(ncores)
    logger.info(f"Using {ncores} cores for parallel processing.")

    # Determine which genomes are available
    logger.info("Determining available mitochondrial genomes...")
    rawsg = glob.glob(script_dir + "/bin/anno/fasta/*.fasta")
    supported_genomes = [x.replace(script_dir + "/bin/anno/fasta/", "").replace(".fasta", "") for x in rawsg]
    logger.info(f"Supported genomes: {supported_genomes}")

    # Input argument is assumed to be a .bam file
    filename, file_extension = os.path.splitext(input)
    if file_extension != ".bam":
        logger.error(f'The input should be an individual .bam file.')
        sys.exit(1)
    if not os.path.exists(input):
        logger.error(f'No file found called "{input}"; please specify a valid .bam file.')
        sys.exit(1)
    if not os.path.exists(input + ".bai"):
        logger.info(f"Attempting to index: {input}")
        pysam.index(input)
        if not os.path.exists(input + ".bai"):
            logger.error('Cannot find index file for input .bam file; please ensure that the .bam file is indexed.')
            sys.exit(1)

    # Determine whether or not we have been supplied barcodes
    if os.path.exists(barcodes) and barcodes != "":
        barcode_known = True
        logger.info("Using provided barcode list.")
    else:
        logger.error('Must specify a known barcode list')
        sys.exit('Must specify a known barcode list with --barcodes')

    # Log the configuration details as a table
    config_details = [
        ("Input", input),
        ("Output", output),
        ("Name", name),
        ("Mitochondrial genome", mito_genome),
        ("Number of cores", ncores),
        ("Barcode tag", barcode_tag),
        ("Barcodes", barcodes),
        ("Min barcode reads", min_barcode_reads),
        ("NHmax", nhmax),
        ("NMmax", nmmax),
        ("Remove duplicates", remove_duplicates),
        ("UMI barcode", umi_barcode),
        ("Handle overlap", handle_overlap),
        ("Low coverage threshold", low_coverage_threshold),
        ("Max java memory", max_javamem),
        ("Proper pairs", proper_pairs),
        ("Base quality", base_qual),
        ("Alignment quality", alignment_quality),
        ("Emit base qualities", emit_base_qualities),
        ("Number of samples", nsamples),
        ("Keep samples", keep_samples),
        ("Ignore samples", ignore_samples),
        ("Keep temp files", keep_temp_files),
        ("Keep QC .bam files", keep_qc_bams),
        ("Skip R", skip_r),
        ("Snake stdout", snake_stdout),
        ("Remove .snakemake", remove_snakemake)
    ]

    logger.info("Configuration for analysis:")
    for key, value in config_details:
        logger.info(f"{key:26}: {value}")

    # Make temporary directory of inputs
    of = output
    tf = of + "/temp"
    bcbd = tf + "/barcoded_bams"  # bcdb = barcoded bam directory
    folders = [of, tf, bcbd, of + "/final"]
    logger.info("Creating output directories...")
    mkfolderout = [make_folder(x) for x in folders]

    # Handle fasta requirements
    logger.info("Handling fasta requirements...")
    fastaf, mito_chr, mito_length = handle_fasta_inference(mito_genome, supported_genomes, script_dir, of)
    logger.info(f"Using fasta file: {fastaf}")
    idxs = pysam.idxstats(input).split("\n")

    # Handle common mtDNA reference genome errors
    bam_length = 0
    for i in idxs:
        if i.split("\t")[0] == mito_chr:
            bam_length = int(i.split("\t")[1])

    if mito_length == bam_length:
        pass
    elif bam_length == 16569:
        fastaf, mito_chr, mito_length = handle_fasta_inference("rCRS", supported_genomes, script_dir, of)
    else:
        logger.info(f"User specified mitochondrial genome does NOT match .bam file; available references are GRCh37, GRCh38, NC_012920, rCRS, GRCm38")
        sys.exit(1)

    # Split barcodes file for parallel processing
    logger.info("Checking if barcode files and BAM files already exist...")
    barcode_files = split_barcodes_file(barcodes, math.ceil(file_len(barcodes) / ncores), output)
    samples = [os.path.basename(os.path.splitext(sample)[0]) for sample in barcode_files]
    samplebams = [of + "/temp/barcoded_bams/" + sample + ".bam" for sample in samples]

    # Check if all barcode files and BAM files already exist
    all_files_exist = all(os.path.exists(bf) for bf in barcode_files) and all(os.path.exists(sb) for sb in samplebams)

    if all_files_exist:
        logger.info("All barcode files and BAM files already exist. Skipping splitting and parallel processing.")
    else:
        logger.info("Splitting barcode file for parallel processing...")
        barcode_files = split_barcodes_file(barcodes, math.ceil(file_len(barcodes) / ncores), output)
        samples = [os.path.basename(os.path.splitext(sample)[0]) for sample in barcode_files]
        samplebams = [of + "/temp/barcoded_bams/" + sample + ".bam" for sample in samples]

        # Parallel processing of barcode files
        logger.info("Starting parallel processing of barcode files...")
        pool = Pool(processes=ncores)
        pool.starmap(split_chunk_file, zip(barcode_files, repeat(script_dir), repeat(input), repeat(bcbd), repeat(barcode_tag), repeat(mito_chr), repeat(umi_barcode)))
        pool.close()
        pool.join()
        logger.info("Finished parallel processing of barcode files.")

    # Create necessary directories for output and logs
    of = output
    tf = of + "/temp"
    logs = of + "/logs"
    folders = [
        logs, of + "/logs/filter", of + "/logs/depth", of + "/fasta", of + "/.internal",
        of + "/.internal/parseltongue", of + "/.internal/samples", of + "/final",
        tf, tf + "/ready_bam", tf + "/temp_bam", tf + "/sparse_matrices"]
    logger.info("Creating additional output directories...")
    mkfolderout = [make_folder(x) for x in folders]

    if remove_duplicates:
        make_folder(of + "/logs/rmdupslogs")

    # Create README files for internal directories
    if not os.path.exists(of + "/.internal/README"):
        with open(of + "/.internal/README", 'w') as outfile:
            outfile.write("This folder creates important (small) intermediate; don't modify it.\n\n")
    if not os.path.exists(of + "/.internal/parseltongue/README"):
        with open(of + "/.internal/parseltongue/README", 'w') as outfile:
            outfile.write("This folder creates intermediate output to be interpreted by Snakemake; don't modify it.\n\n")
    if not os.path.exists(of + "/.internal/samples/README"):
        with open(of + "/.internal/samples/README", 'w') as outfile:
            outfile.write("This folder creates samples to be interpreted by Snakemake; don't modify it.\n\n")

    # Write sample bam file paths to internal directory
    for i in range(len(samples)):
        with open(of + "/.internal/samples/" + samples[i] + ".bam.txt", 'w') as outfile:
            outfile.write(samplebams[i])

    # Create YAML configuration file for Snakemake
    logger.info("Creating Snakemake configuration file...")
    dict1 = {
        'input_directory': sqs(input), 'output_directory': sqs(output), 'script_dir': sqs(script_dir),
        'fasta_file': sqs(fastaf), 'mito_chr': sqs(mito_chr), 'mito_length': sqs(mito_length), 'name': sqs(name),
        'base_qual': sqs(base_qual), 'remove_duplicates': sqs(remove_duplicates), 'handle_overlap': sqs(handle_overlap),
        'low_coverage_threshold': sqs(low_coverage_threshold), 'barcode_tag': sqs(barcode_tag), 'umi_barcode': sqs(umi_barcode),
        'alignment_quality': sqs(alignment_quality), 'emit_base_qualities': sqs(emit_base_qualities),
        'proper_paired': sqs(proper_pairs), 'NHmax': sqs(nhmax), 'NMmax': sqs(nmmax), 'max_javamem': sqs(max_javamem)
    }

    y_s = of + "/.internal/parseltongue/snake.scatter.yaml"
    with open(y_s, 'w') as yaml_file:
        yaml = YAML()
        yaml.default_flow_style = False
        yaml.dump(dict1, yaml_file)
    logger.info("Snakemake configuration file created.")

    # Write configuration details to configuration.txt
    logger.info("Writing configuration details to file...")
    with open(f"{logs}/configuration.txt", 'a') as param_file:
        param_file.write("Configuration details:\n")
        for key, value in config_details:
            param_file.write(f"{key:25}: {value}\n")
    logger.info("Configuration details written to file.")

    import subprocess

    # Path to the lock file
    lock_file = os.path.join(script_dir, ".snakemake", "lock")

    # Check if the lock file exists
    if os.path.exists(lock_file):
        # Unlock the working directory
        snakecmd_unlock = f'snakemake --unlock --snakefile {script_dir}/bin/snake/Snakefile.tenx'
        logger.info(f"Unlocking Snakemake working directory: {snakecmd_unlock}")

        try:
            result = subprocess.run(snakecmd_unlock, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            logger.info("Successfully unlocked Snakemake working directory.")
        except subprocess.CalledProcessError as e:
            logger.error("Failed to unlock Snakemake working directory.")
            sys.exit(1)
    else:
        logger.info("No lock file found. Proceeding without unlocking.")

    # Execute Snakemake command
    snake_log = logs + "/snakemake.log"
    snake_log_out = "" if snake_stdout else f'> {snake_log} 2>&1'
    snakecmd_tenx = f'snakemake --snakefile {script_dir}/bin/snake/Snakefile.tenx --cores {ncores} --config cfp="{y_s}" {snake_log_out}'

    logger.info(f"Executing Snakemake command: {snakecmd_tenx}")

    try:
        result = subprocess.run(snakecmd_tenx, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        logger.info(f"Snakemake failed, see {snake_log} for details.")
        sys.exit(1)

    if not skip_r:
        logger.info("Running R script to generate .rds object...")
        Rcall = f"Rscript {script_dir}/bin/R/toRDS.R {output}/final {name}"
        os.system(Rcall)
        logger.info("R script completed.")

    if remove_snakemake:
        logger.info("Removing .snakemake directory as --remove-snakemake was specified.")
        shutil.rmtree(".snakemake")

    if keep_qc_bams:
        logger.info("Final bams retained since --keep-qc-bams was specified.")
        shutil.move(of + "/temp/ready_bam", of + "/qc_bam")

    if not keep_temp_files:
        logger.info("Removing temporary files...")
        for folder in [of + "/fasta", of + "/.internal", of + "/temp"]:
            shutil.rmtree(folder)
        logger.info("Temporary files removed.")
    else:
        logger.info("Temporary files not deleted since --keep-temp-files was specified.")

    logger.info("Successfully created final output files. Analysis complete.")