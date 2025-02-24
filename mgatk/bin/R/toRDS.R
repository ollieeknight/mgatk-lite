#!/usr/bin/env Rscript
suppressMessages(suppressWarnings(library(tools)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(SummarizedExperiment)))
suppressMessages(suppressWarnings(library(GenomicRanges)))
suppressMessages(suppressWarnings(library(data.table)))

options(warn=-1)

# Function to import data from a file
importDT <- function(file) {
  ext <- tools::file_ext(file)
  if (ext == "gz") {
    return(suppressMessages(data.table::fread(paste0("zcat < ", file), stringsAsFactors = TRUE)))
  } else if (ext %in% c("txt", "csv", "tsv")) {
    return(suppressMessages(data.table::fread(file, stringsAsFactors = TRUE)))
  } else {
    stop("Provide a valid file format (.gz, .txt, .csv, or .tsv)")
  }
}

# Function to import sparse matrices
importSMs <- function(file, samplesOrder, maxpos, maxsamples) {
  dt <- importDT(file)
  dt[[2]] <- factor(as.character(dt[[2]]), levels = samplesOrder)
  
  base_quals <- ncol(dt) == 6
  count_fw_idx <- if (base_quals) 3 else 3
  count_rev_idx <- if (base_quals) 5 else 4
  
  counts_fw <- Matrix::sparseMatrix(i = c(dt[[1]], maxpos), j = c(as.numeric(dt[[2]]), maxsamples), x = c(dt[[count_fw_idx]], 0))
  counts_rev <- Matrix::sparseMatrix(i = c(dt[[1]], maxpos), j = c(as.numeric(dt[[2]]), maxsamples), x = c(dt[[count_rev_idx]], 0))
  
  if (base_quals) {
    qual_fw <- Matrix::sparseMatrix(i = c(dt[[1]], maxpos), j = c(as.numeric(dt[[2]]), maxsamples), x = c(dt[[4]], 0))
    qual_rev <- Matrix::sparseMatrix(i = c(dt[[1]], maxpos), j = c(as.numeric(dt[[2]]), maxsamples), x = c(dt[[6]], 0))
    return(list("counts_fw" = counts_fw, "qual_fw" = qual_fw, "counts_rev" = counts_rev, "qual_rev" = qual_rev))
  } else {
    return(list("counts_fw" = counts_fw, "counts_rev" = counts_rev))
  }
}

# Function to import mitochondrial data
importMito.explicit <- function(Afile, Cfile, Gfile, Tfile, coverageFile, depthFile, referenceAlleleFile, mitoChr = "chrM") {
  variantFiles <- list(Afile, Cfile, Gfile, Tfile)
  metaFiles <- list(coverageFile, depthFile, referenceAlleleFile)
  
  lapply(c(variantFiles, metaFiles), function(file) stopifnot(length(file) == 1))
  
  cov <- importDT(coverageFile)
  ref <- importDT(referenceAlleleFile)
  maxpos <- max(ref[[1]])
  samplesOrder <- levels(cov[[2]])
  maxsamples <- length(samplesOrder)
  
  covmat <- Matrix::sparseMatrix(i = c(cov[[1]], maxpos), j = c(as.numeric(cov[[2]]), maxsamples), x = c(cov[[3]], 0))
  
  ACGT <- lapply(variantFiles, importSMs, samplesOrder, maxpos, maxsamples)
  names(ACGT) <- c("A", "C", "G", "T")
  
  depth <- data.frame(importDT(depthFile))
  sdf <- merge(data.frame(sample = samplesOrder), depth, by.x = "sample", by.y = "V1")
  rownames(sdf) <- samplesOrder
  colnames(sdf) <- c("sample", "depth")
  
  row_g_cov <- GenomicRanges::GRanges(seqnames = mitoChr, IRanges::IRanges(1:maxpos, width = 1))
  GenomicRanges::mcols(row_g_cov) <- data.frame(refAllele = toupper(ref[[2]][1:maxpos]))
  
  refallele <- data.frame(pos = 1:maxpos, ref = as.character(toupper(ref[[2]][1:maxpos])), stringsAsFactors = FALSE)
  signac_counts <- do.call(rbind, lapply(ACGT, function(x) rbind(x[["counts_fw"]], x[["counts_rev"]])))
  
  signac_object <- list("counts" = signac_counts, "depth" = sdf, "refallele" = refallele)
  
  assays_list <- list(
    "A_counts_fw" = ACGT[["A"]][["counts_fw"]], "A_counts_rev" = ACGT[["A"]][["counts_rev"]],
    "C_counts_fw" = ACGT[["C"]][["counts_fw"]], "C_counts_rev" = ACGT[["C"]][["counts_rev"]],
    "G_counts_fw" = ACGT[["G"]][["counts_fw"]], "G_counts_rev" = ACGT[["G"]][["counts_rev"]],
    "T_counts_fw" = ACGT[["T"]][["counts_fw"]], "T_counts_rev" = ACGT[["T"]][["counts_rev"]],
    "coverage" = covmat
  )
  
  if (length(ACGT[["A"]]) > 2) {
    assays_list <- c(assays_list, list(
      "A_qual_fw" = ACGT[["A"]][["qual_fw"]], "A_qual_rev" = ACGT[["A"]][["qual_rev"]],
      "C_qual_fw" = ACGT[["C"]][["qual_fw"]], "C_qual_rev" = ACGT[["C"]][["qual_rev"]],
      "G_qual_fw" = ACGT[["G"]][["qual_fw"]], "G_qual_rev" = ACGT[["G"]][["qual_rev"]],
      "T_qual_fw" = ACGT[["T"]][["qual_fw"]], "T_qual_rev" = ACGT[["T"]][["qual_rev"]]
    ))
  }
  
  SE <- SummarizedExperiment::SummarizedExperiment(assays = assays_list, colData = S4Vectors::DataFrame(sdf), rowData = row_g_cov)
  return(list(SE, signac_object))
}

# Function to parse the folder hierarchy
importMito <- function(folder, ...) {
  files <- list.files(folder, full.names = TRUE)
  
  checkGrep <- function(hit) {
    if (length(hit) != 1) {
      stop("Improper folder specification; file missing / extra file present. See documentation")
    } else {
      return(hit)
    }
  }
  
  Afile <- files[checkGrep(grep(".A.txt", files))]
  Cfile <- files[checkGrep(grep(".C.txt", files))]
  Gfile <- files[checkGrep(grep(".G.txt", files))]
  Tfile <- files[checkGrep(grep(".T.txt", files))]
  coverageFile <- files[checkGrep(grep(".coverage.txt", files))]
  depthFile <- files[checkGrep(grep(".depthTable.txt", files))]
  referenceAlleleFile <- files[checkGrep(grep("refAllele.txt", files))]
  
  sv <- strsplit(gsub("_refAllele.txt", "", basename(referenceAlleleFile)), split = "[.]")[[1]]
  mitoChr <- sv[length(sv)]
  
  SElist <- importMito.explicit(Afile, Cfile, Gfile, Tfile, coverageFile, depthFile, referenceAlleleFile, mitoChr, ...)
  return(SElist)
}

# Command line i/o
args <- commandArgs(trailingOnly = TRUE)
folder <- args[1]
name <- args[2]

SElist <- importMito(folder)
saveRDS(SElist[[1]], file = paste0(folder, "/", name, ".rds"))
saveRDS(SElist[[2]], file = paste0(folder, "/", name, ".signac.rds"))
