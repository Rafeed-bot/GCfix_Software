import pysam
import os
import subprocess
import pandas as pd
import numpy as np


input_bam = 'tagged_sample.bam' # bam file path tagged with "GC" tag by the software
bin_locations = pd.read_csv('../hg19__GRCh37/correction_factor_estimation_bin_locations.csv').values.tolist() # load the genomic regions you are interested in
total_read_no = 0
pybam = pysam.AlignmentFile(input_bam, "rb") # opening tagged bam
for bin_location in bin_locations:
    contig, start, end = bin_location[0], int(bin_location[1]), int(bin_location[2]) # obtaining single genomic region
    for read in pybam.fetch(contig, start, end): # fetching reads from the genomic region
        correction_factor = read.get_tag('GC') # getting GC correction tag value for the fetched read
        total_read_no += correction_factor # adding correction factor to the total read no instead of adding +1 
pybam.close() # closing bam
print(total_read_no)