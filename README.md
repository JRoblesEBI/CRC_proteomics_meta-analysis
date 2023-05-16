# CRC_proteomics_meta-analysis
DATA

Batch analysis output files can be found in "Solid samples" and "Secreted samples" folders. Each folder contains:
  1. The protein expression matrix (ppb.iBAQ): "proteinGroups_ppb_final.txt"
  2. The bin-transformed protein abundances (ranked bins from 1 to 5): "proteinGroups_ppb_final_binning.txt"
  3. The annotation file containing sample name, dataset, condition and source: "Sample_classification.txt"

CODE

R script used for performing Cox regression analysis to the protein expression matrix can be found in "survival analysis" folder: "Cox_regression.R"
