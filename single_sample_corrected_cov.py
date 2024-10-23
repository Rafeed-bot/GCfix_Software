import pysam
import numpy as np
import os, sys, subprocess, time, shutil
from multiprocessing import Pool
import pandas as pd


input_bam = sys.argv[1]
output_cov_csv = sys.argv[2]
correction_path = sys.argv[3]
mapq = int(sys.argv[4])
start_len = int(sys.argv[5])
end_len = int(sys.argv[6])
CPU = int(sys.argv[7])
reference_genome_path = sys.argv[8]
bin_location_file = sys.argv[9]

bin_locations = pd.read_csv(bin_location_file).values.tolist()
bin_no = len(bin_locations)
pro_bin_locations = []
for i in range(len(bin_locations)):
    pro_bin_locations.append([i, bin_locations[i][0], bin_locations[i][1], bin_locations[i][2]])
lag = 10
final_cov_array = np.zeros((bin_no))

command = f'samtools view -H {input_bam} | grep @SQ | grep chr'
p = subprocess.run(command, shell=True, capture_output=True, text=True)
flag = len(p.stdout)


def get_cov(sub_bin_locations):
    global reference_genome_path, input_bam, bin_no, correction_path, flag
    reference_genome = pysam.FastaFile(reference_genome_path)
    all_ref_contigs = reference_genome.references
    ref_flag = 0
    if 'chr' in all_ref_contigs[0]:
        ref_flag = 1
    correction_factors = np.array( pd.read_csv(correction_path).values.tolist() )

    cov_array = np.zeros((bin_no))
    with pysam.AlignmentFile(input_bam, "rb") as bam_in:
        for bin_location in sub_bin_locations:
            bin_index, contig, start, end = int(bin_location[0]), bin_location[1], int(bin_location[2]), int(bin_location[3])
            ref_contig = contig if ref_flag==1 else contig[3:]
            contig = contig[3:] if flag==0 else contig
            for read in bam_in.fetch(contig=contig, start=start, end=end): 
                if (read.flag & 1024 == 0) and (read.flag & 2048 == 0) and (read.flag & 4 == 0) and (read.flag & 8 == 0):
                    if read.mapping_quality>=mapq and (read.reference_name == read.next_reference_name):
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
                                    correction_factor = correction_factors[length-51, GC_content]
                                    cov_array[bin_index] += correction_factor
    return cov_array


p = Pool(processes=CPU)
sub_bin_locations = np.array_split(pro_bin_locations, CPU)
cov_array_list = p.map(get_cov, sub_bin_locations, 1)

p.close()
p.join()

for i, cov_array in enumerate(cov_array_list):
    final_cov_array = final_cov_array + cov_array
csv_list = []
for i in range(len(bin_locations)):
    contig, start, end = bin_locations[i][0], bin_locations[i][1], bin_locations[i][2]
    coverage = final_cov_array[i]
    csv_list.append([contig, start, end, coverage])
df = pd.DataFrame(csv_list)
df.to_csv(output_cov_csv, index=False, header=['contig', 'start', 'end', 'coverage'])