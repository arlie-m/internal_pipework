3-Barcoding
================
Arlie McCarthy
2026-06-26

## Importing the sequence files

``` r
# Set the directory containing the sequence files
dir <- here::here("data", "renamed_sequences")

# Get a list of all files in the directory
files <- list.files(dir, pattern = "\\.ab1$", full.names = TRUE)

# Create a vector to store the file names with suffixes
file_names <- character(length(files))

# Process each file and create a SangerRead object
sanger_reads <- lapply(files, function(file) {
  if (grepl("_F\\.ab1$", file)) {
    read_type <- "Forward Read"
  } else {
    read_type <- "Reverse Read"
  }
  
  SangerRead(
    readFeature = read_type,
    readFileName = file,
    geneticCode = GENETIC_CODE,
    TrimmingMethod = "M2",
    M1TrimmingCutoff = NULL,
    M2CutoffQualityScore = 20,
    M2SlidingWindowSize = 10,
    baseNumPerRow = 100,
    heightPerRow = 200,
    signalRatioCutoff = 0.33,
    showTrimmed = TRUE
  )
})
```

# Save processed and trimmed sequences

``` r
lapply(sanger_reads, qualityBasePlot)
lapply(sanger_reads, writeFastaSR, outputDir = here::here("outputs", "renamed_sequences_fasta"))
```
