# salicylate_barcode_experiment
R scripts used for analysis of mixed, differentially barcoded mutant/wild-type communities

## Overall analysis procedure

- Create sample metadata spreadsheets (one for library association and one for each experiment)
- Collect/concatenate sequence files as necessary and distribute into analysis directories
- Remove reverse reads as Bartender requires only forward reads
- Calculate read numbers for each samples from sequence files
- Enter read numbers into respective sample metadata spreadsheets
- Separately commence bartender analysis for library association and then for each experiment, according to Bartender instructions
- Analyze library_seqs to create a Barcode/Phenotype dictionary (barcode_metadata.csv). This includes discarding all barcodes that have more than one genotype assigned to them
- Copy Barcode/Genotype dictionary to experimental sample analysis directories
- Anaylze experimental samples
