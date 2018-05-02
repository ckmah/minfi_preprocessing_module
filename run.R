suppressMessages(suppressWarnings(library("optparse")))
suppressMessages(suppressWarnings(library("minfi")))

# Parse input arguments
parser = OptionParser()

parser <- add_option(parser, c("-d", "--data"), help = "zip or gzip containing DNA methylation microarray data")
parser <- add_option(parser, c("-n", "--normalization"), help = "normalization method")
parser <- add_option(parser, c("-t", "--outputtype"), help = "beta, m-values, or MethylSet (.rds file) output value type [default beta]")

args = parse_args(parser)

if (args$outputtype == "MethylSet" && args$normalization != "None") {
  write("'None' normalization must be performed to obtain the MethylSet output type.", stdout())
  stop()
}


# Path to folder containing raw experiment .IDAT files. Will read Illumina format
# SampleSheet.csv to load data if found.
unzip(args$data)  # extract
data.paths <- unzip(args$data, list = TRUE)  # list extracted
data.folder <- paste(getwd(), strsplit(data.paths[[1]][1], "/")[[1]][1], sep = "/")

# Try Load sample sheet
if (any(grepl(".csv$", list.files(data.folder)))) {
  targets <- read.metharray.sheet(data.folder)

  # remove unannotated samples
  targets.rmdups <- targets[targets$Basename != "character(0)", ]
  experiment.rgset <- read.metharray.exp(targets = targets.rmdups)

  # Recursively find all .idat files in data.folder if sample sheet not found
} else {
  experiment.rgset <- read.metharray.exp(base = data.folder, recursive = TRUE)
}

# Convert to MethylSet
methyl.set <- preprocessRaw(experiment.rgset)

# QC plots
pdf("qcPlots.pdf")

# Plot median log2 meth vs unmeth intensities
qc <- getQC(methyl.set)
plotQC(qc)

# Plot Beta value distributions
phenoData <- pData(methyl.set)
if ("Sample_Group" %in% names(phenoData)) {
  densityPlot(methyl.set, sampGroups = phenoData$Sample_Group)
} else {
  densityPlot(methyl.set)
}
dev.off()

# Normalization function
performPreprocessing <- function(rg.set, preprocess.method = "") {
  if (preprocess.method == "preprocessFunnorm") {
    write("Perform background subtraction with dye-bias normalization and infer between-array technical variation.",
      stdout())
    m.set <- preprocessFunnorm(rg.set)
  } else if (preprocess.method == "preprocessQuantile") {
    write("Perform stratified quantile normalization.", stdout())
    m.set <- preprocessQuantile(rg.set)
  } else if (preprocess.method == "None") {
    write("Preprocessing without normalization.", stdout())
    m.set <- preprocessRaw(rg.set)
  } else {
    write("Preprocessing method not recognized. Available methods: noob, raw.",
      stdout())
    return()
  }

  # Convert to GenomicRatioSet as needed
  if (class(m.set) == "MethylSet" && args$outputtype != "MethylSet") {
    gm.set <- mapToGenome(m.set)
    gr.set <- ratioConvert(gm.set)
    return(gr.set)
  } else {
    return(m.set)
  }
}

output.set <- performPreprocessing(experiment.rgset, args$normalization)

write("Removing loci with SNPs...", stdout())
snps <- getSnpInfo(output.set)
output.set <- addSnpInfo(output.set)
output.set <- dropLociWithSnps(output.set, snps = c("SBE", "CpG"))

# Convert to output format
write("Saving output...", stdout())
if (args$outputtype == "beta") {
  output.data <- getBeta(output.set)
} else if (args$outputtype == "m-values") {
  output.data <- getM(output.set)
} else if (args$outputtype == "MethylSet") {
  output.data <- output.set
} else {
  write("Output file format not recognized. Available formats: 'beta', 'm-values', or 'MethylSet'",
    stdout())
  return()
}

# Write to file as appropriate
if (class(output.data) == "MethylSet") {
  saveRDS(output.data, "methyl_set.rds")
} else {
  write.table(output.data, file.path(getwd(), paste("methyl-", args$outputtype, ".txt", sep="")), quote = FALSE, sep = "\t")
}

write("Cleaning up intermediate files...", stdout())
unlink("Demo Data EPIC", recursive = TRUE)
write("Done.", stdout())
