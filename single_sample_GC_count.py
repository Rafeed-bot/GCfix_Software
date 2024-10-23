import pysam
import subprocess
import pandas as pd
import numpy as np
import sys
from multiprocessing import Pool


input_bam = sys.argv[1] 
GC_count_npy = sys.argv[2] 
mapq = int(sys.argv[3])
start_len = int(sys.argv[4])
end_len = int(sys.argv[5])
CPU = int(sys.argv[6])
reference_genome_path = sys.argv[7]
bin_location_file = sys.argv[8]

bin_locations = pd.read_csv(bin_location_file).values.tolist()
lag = 10
final_GC_array = np.zeros((end_len-start_len+1, 101))


def count_GC(sub_bin_locations):
    global reference_genome_path, input_bam, start_len, end_len
    GC_array = np.zeros((end_len-start_len+1, 101))
    pybam = pysam.AlignmentFile(input_bam, "rb")
    command = 'samtools view -H {bamPath} | grep @SQ | grep chr'.format(bamPath=input_bam)
    p = subprocess.run(command, shell=True, capture_output=True, text=True)
    bam_flag = len(p.stdout)

    reference_genome = pysam.FastaFile(reference_genome_path)
    all_ref_contigs = reference_genome.references
    ref_flag = 0
    if 'chr' in all_ref_contigs[0]:
        ref_flag = 1

    for bin_location in sub_bin_locations:
        contig, start, end = bin_location[0], int(bin_location[1]), int(bin_location[2])
        ref_contig = contig
        if ref_flag==0: # checks whether ref contig has 'chr' or not
            ref_contig = ref_contig[3:]
        if bam_flag==0: # checks whether bam contig has 'chr' or not
            contig = contig[3:]
        for read in pybam.fetch(contig, start, end): # window
            if (read.flag & 1024 == 0) and (read.flag & 2048 == 0) and (read.flag & 4 == 0) and (read.flag & 8 == 0):
                if read.mapping_quality >= mapq and (read.reference_name == read.next_reference_name):
                    length = read.template_length
                    if length>=start_len and length<=end_len:
                        ref_start = read.reference_start + lag
                        ref_end = read.reference_start + length - lag
                        ref_seq = reference_genome.fetch(ref_contig, ref_start, ref_end)
                        if len(ref_seq)!=0:
                            GC_cnt = ref_seq.count('G') + ref_seq.count('C') + ref_seq.count('g') + ref_seq.count('c')
                            AT_cnt = ref_seq.count('A') + ref_seq.count('T') + ref_seq.count('a') + ref_seq.count('t')
                            total_cnt = GC_cnt + AT_cnt
                            valid_percent = float( total_cnt/len(ref_seq) )
                            if valid_percent>=0.9:
                                GC_content = float( GC_cnt/total_cnt )
                                GC_content = round(round(GC_content, 2) * 100)
                                GC_array[length-start_len][GC_content] += 1
    return GC_array


p = Pool(processes=CPU)
sub_bin_locations = np.array_split(bin_locations, CPU)
GC_array_list = p.map(count_GC, sub_bin_locations, 1)

p.close()
p.join()

for i, GC_array in enumerate(GC_array_list):
    final_GC_array = final_GC_array + GC_array
np.save(GC_count_npy, final_GC_array)