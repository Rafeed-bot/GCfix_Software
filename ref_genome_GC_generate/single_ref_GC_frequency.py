import pysam
import pandas as pd
import numpy as np
import sys
from multiprocessing import Pool
import time


fragment_length = int(sys.argv[1]) # one length from 51-400
lag = 10
mappable_regions_path = sys.argv[2]
ref_seq_path = sys.argv[3]
out_dir = sys.argv[4]
out_file = f'{out_dir}/{fragment_length}.npy'
CPU = int(sys.argv[5])

bin_locations = pd.read_csv(mappable_regions_path).values.tolist()
bin_no = len(bin_locations)
GC_window_size = fragment_length - 2*lag
final_GC_array = np.zeros((101))


def count_GC(sub_bin_locations):
    global ref_seq_path, fragment_length, GC_window_size
    GC_array = np.zeros((101))
    ref_seq = pysam.FastaFile(ref_seq_path)
    all_ref_contigs = ref_seq.references
    ref_flag = 0
    if 'chr' in all_ref_contigs[0]:
        ref_flag = 1

    for bin_location in sub_bin_locations:
        chrom, start, end = bin_location[0], int(bin_location[1]), int(bin_location[2])
        ref_contig = chrom
        if ref_flag==0:
            ref_contig = ref_contig[3:]
        fetched = ref_seq.fetch(ref_contig, start, end)
        fetched = np.array(list(fetched.upper()))
        fetched[np.isin(fetched, ['A','T','W'])] = 0
        fetched[np.isin(fetched, ['C','G','S'])] = 1
        rng = np.random.default_rng(start)
        fetched[np.isin(fetched, ['N','R','Y','K','M','B','D','H','V'])] = rng.integers(2, size=len(fetched[np.isin(fetched, ['N','R','Y','K','M','B','D','H','V'])])) 
        fetched = fetched.astype(float)
        n=len(fetched)
        window_sum = int(sum(fetched[lag : fragment_length-lag]))

        GC_content = float( window_sum/GC_window_size )
        GC_content = round(round(GC_content, 2) * 100)
        GC_array[GC_content] += 1
        for j in range(n-fragment_length):
            window_sum = int(window_sum - fetched[j+lag] + fetched[j+fragment_length-lag])
            GC_content = float( window_sum/GC_window_size )
            GC_content = round(round(GC_content, 2) * 100)
            GC_array[GC_content] += 1
    return GC_array

p = Pool(processes=CPU)
sub_bin_locations = np.array_split(bin_locations, CPU)
GC_array_list = p.map(count_GC, sub_bin_locations, 1)

p.close()
p.join()

for i, GC_array in enumerate(GC_array_list):
    final_GC_array = final_GC_array + GC_array
np.save(out_file, final_GC_array)