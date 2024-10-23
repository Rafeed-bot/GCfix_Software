import pysam
import numpy as np
import pandas as pd


input_bam = 'sample.bam' # input bam file for GC correction
correction_csv = 'sample_correction_factors.csv' # estimated GC correction factors by the software
reference_genome_path = '../hg19__GRCh37/ref.fa' # input bam was aligned to this reference genome
bin_locations = pd.read_csv('../hg19__GRCh37/correction_factor_estimation_bin_locations.csv').values.tolist() # you are interested in reads from these genomic regions
start_length, end_length = 51, 400 # fragment length range
lag = 10 # required for getting GC content from reference genome for each read

# loading the correction factors into an array of size (end_length - start_length + 1) X 101 -> 101 represents GC content of 0%, 1%, ..., 100%
correction_factors = np.array( pd.read_csv(correction_csv).values.tolist() ) 
pybam = pysam.AlignmentFile(input_bam, "rb") # loading the input bam
reference_genome = pysam.FastaFile(reference_genome_path) # loading the ref genome
total_fragment_no = 0
for bin_location in bin_locations:
    contig, start, end = bin_location[0], int(bin_location[1]), int(bin_location[2]) # obtaining single genomic region
    ref_contig = contig[3:] # getting appropriate name for reference genome contig
    for read in pybam.fetch(contig, start, end): # # fetching reads from the genomic region
        # read filterings applied -> you can change these as per your requirement
        if (read.flag & 1024 == 0) and (read.flag & 2048 == 0) and (read.flag & 4 == 0) and (read.flag & 8 == 0):
            if read.mapping_quality >= 30 and (read.reference_name == read.next_reference_name):
                length = read.template_length
                if length>=start_length and length<=end_length:
                    # getting GC content of the read using reference genome
                    ref_start = read.reference_start + lag
                    ref_end = read.reference_start + length - lag
                    ref_seq = reference_genome.fetch(ref_contig, ref_start, ref_end)
                    if len(ref_seq)!=0:
                        GC_cnt = ref_seq.count('G') + ref_seq.count('C') + ref_seq.count('g') + ref_seq.count('c')
                        AT_cnt = ref_seq.count('A') + ref_seq.count('T') + ref_seq.count('a') + ref_seq.count('t')
                        total_cnt = GC_cnt + AT_cnt
                        GC_content = float( GC_cnt/total_cnt )
                        GC_content = round(round(GC_content, 2) * 100)
                        valid_percent = float( total_cnt/len(ref_seq) )
                        # we want reads to have at least 90% A, T, C, G -> a sign of good quality reads
                        if valid_percent>=0.9: 
                            # obtaining correction factor for the read according to corresponding GC content and fragment length
                            correction_factor = correction_factors[length-start_length, GC_content]
                            total_fragment_no += correction_factor # adding correction factor to the total read no instead of adding +1 
pybam.close() # closing bam
print(total_fragment_no)

                            