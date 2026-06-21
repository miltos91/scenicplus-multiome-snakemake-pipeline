library(argparse)

# ---- Parameters (passed in by Snakemake from config.yaml) ------------------
# Instead of reading each DAR table by hand, we loop over the population names
# listed in the config. input_format picks how the DARs arrive:
#   "csv" = Signac FindDAR output DARs_{population}.csv  (convert to BED)
#   "bed" = ready-made {population}.bed files            (just copy chrom/start/end)
parser <- ArgumentParser()
parser$add_argument("--input_dir",    required = TRUE)
parser$add_argument("--output_dir",   required = TRUE)
parser$add_argument("--populations",  required = TRUE,
                    help = "Comma-separated list of full population names")
parser$add_argument("--input_format", default = "csv")
args <- parser$parse_args()

options(scipen = 999)

input_dir    <- args$input_dir
outdir       <- args$output_dir
populations  <- trimws(strsplit(args$populations, ",")[[1]])
input_format <- tolower(args$input_format)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Loop over every population and write one BED file each
for (name in populations) {

  out_path <- file.path(outdir, paste0(name, ".bed"))

  # ---- Ready-made BED: just keep the first three columns (chrom/start/end) ----
  if (input_format == "bed") {
    src <- file.path(input_dir, paste0(name, ".bed"))
    if (!file.exists(src)) {
      message(paste("Skipping:", name, "- BED not found."))
      next
    }
    bed <- read.table(src, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    if (nrow(bed) == 0 || ncol(bed) < 3) {
      message(paste("Skipping:", name, "- BED empty or has < 3 columns."))
      next
    }
    write.table(
      bed[, 1:3],
      file      = out_path,
      sep       = "\t",
      quote     = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    next
  }

  # ---- CSV (Signac FindDAR output): row names are peaks like chr1-1000-2000 ----
  csv_path <- file.path(input_dir, paste0("DARs_", name, ".csv"))
  if (!file.exists(csv_path)) {
    message(paste("Skipping:", name, "- file not found."))
    next
  }
  df <- read.csv(csv_path, row.names = 1)

  if (nrow(df) == 0) {
    message(paste("Skipping:", name, "- No significant DARs found."))
    next
  }

  peaks <- rownames(df)

  coords <- do.call(rbind, strsplit(peaks, "-", fixed = TRUE))

  bed <- data.frame(
    chrom = coords[, 1],
    start = as.integer(coords[, 2]),
    end   = as.integer(coords[, 3]),
    stringsAsFactors = FALSE
  )

  write.table(
    bed,
    file      = out_path,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}
