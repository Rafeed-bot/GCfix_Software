#!/bin/bash

# R has to be installed
## Dependencies: "tidyverse", "data.table", "svglite", "argparse", "arrow"

maindir=$(realpath $PWD)
cd $maindir

find $maindir/example/NDRCoverage -type f -name "*-ndr_coverage.csv.gz" > "${maindir}/ndr_coverage_input.list"

geneset="${maindir}/example/geneset/GTEx_Blood_TPM-500_sites.csv"
label="example"

Rscript ./1_ndr_range.R \
    --maindir "${maindir}" \
    --covlist "${maindir}/ndr_coverage_input.list" \
    --geneset "$geneset" \
    --label "$label"

ndr="${maindir}/example/NDR/example-ndr.parquet"
    
Rscript ./2_ndr_plots.R \
    --maindir "${maindir}" \
    --ndr "$ndr" \
    --geneset "$geneset" \
    --label "$label"
