import pysam
import numpy as np
import os, sys, subprocess, time, shutil
from multiprocessing import Pool
import pandas as pd


input_bam = sys.argv[1]
output_bam = sys.argv[2]
correction_path = sys.argv[3]
mapq = int(sys.argv[4])
start_len = int(sys.argv[5])
end_len = int(sys.argv[6])
CPU = int(sys.argv[7])
reference_genome_path = sys.argv[8]
bin_location_file = sys.argv[9]
temp_folder = sys.argv[10]
if os.path.isdir(temp_folder)==False:
    os.mkdir(temp_folder)

bin_locations = pd.read_csv(bin_location_file).values.tolist()
bin_no = len(bin_locations)
lag = 10

command = f'samtools view -H {input_bam} | grep @SQ | grep chr'
p = subprocess.run(command, shell=True, capture_output=True, text=True)
flag = len(p.stdout)


def tag_bam(sub_bin_locations):
    global reference_genome_path, input_bam, bin_no, correction_path, flag
    reference_genome = pysam.FastaFile(reference_genome_path)
    all_ref_contigs = reference_genome.references
    ref_flag = 0
    if 'chr' in all_ref_contigs[0]:
        ref_flag = 1
    correction_factors = np.array( pd.read_csv(correction_path).values.tolist() )
    contig, start, end = sub_bin_locations[0][0], sub_bin_locations[0][1], sub_bin_locations[0][2]
    output_bam = f'{temp_folder}/split_{contig}_{start}_{end}.bam'

    with pysam.AlignmentFile(input_bam, "rb") as bam_in:
        with pysam.AlignmentFile(output_bam, "wb", header=bam_in.header) as bam_out:
            for bin_location in sub_bin_locations:
                contig, start, end = bin_location[0], int(bin_location[1]), int(bin_location[2])
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
                                        read.set_tag("GC", correction_factor)
                                        bam_out.write(read)
    return output_bam



p = Pool(processes=CPU)
sub_bin_locations = np.array_split(bin_locations, CPU)
_ = p.map(tag_bam, sub_bin_locations, 1)

p.close()
p.join()


subprocess.run(f'samtools merge -@ 32 {temp_folder}/merged.bam {temp_folder}/split_*.bam', shell=True)
subprocess.run(f'samtools sort -@ 32 -o {output_bam} {temp_folder}/merged.bam', shell=True)
subprocess.run(f'samtools index -@ 32 {output_bam}', shell=True)
shutil.rmtree(temp_folder)
