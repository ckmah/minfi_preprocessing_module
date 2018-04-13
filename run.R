suppressMessages(suppressWarnings(library("optparse")))
suppressMessages(suppressWarnings(library("minfi")))

# Parse input arguments
parser = OptionParser()

parser <- add_option(parser, c("-d", "--data"), help = "zip or gzip containing DNA methylation microarray data")
parser <- add_option(parser, c("-n", "--normalization"), default = "None", help = "normalization method")
parser <- add_option(parser, c("-t", "--outputtype"), default = "beta", help = "beta or m-values output value type [default beta]")

args = parse_args(parser)

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

  if (class(m.set) == "MethylSet") {
    pdf("qcPlots.pdf")

    # Plot median log2 meth vs unmeth intensities
    qc <- getQC(m.set)
    plotQC(qc)

    # Plot Beta value distributions
    phenoData <- pData(m.set)
    if ("Sample_Group" %in% names(phenoData)) {
      densityPlot(m.set, sampGroups = phenoData$Sample_Group)
    } else {
      densityPlot(m.set)
    }
    dev.off()

    # Convert to GenomicRatioSet as needed
    gm.set <- mapToGenome(m.set)
    gr.set <- ratioConvert(gm.set)
  } else {
    gr.set <- m.set
  }

  return(gr.set)
}

gr.set <- performPreprocessing(experiment.rgset, args$normalization)

write("Removing loci with SNPs...", stdout())
snps <- getSnpInfo(gr.set)
gr.set <- addSnpInfo(gr.set)
gr.set <- dropLociWithSnps(gr.set, snps = c("SBE", "CpG"))

# Convert to output format
write("Saving output...", stdout())
if (args$outputtype == "beta") {
  output.data <- getBeta(gr.set)
} else if (args$outputtype == "m-values") {
  output.data <- getM(gr.set)
} else {
  write("Output file format not recognized. Available formats: 'beta','m-values'",
    stdout())
  return()
}

write.table(output.data, file.path(getwd(), "methyl.txt"), quote = FALSE, sep = "\t")

write("Cleaning up intermediate files...")
unlink("Demo Data EPIC", recursive = TRUE)
write("Done.", stdout())
