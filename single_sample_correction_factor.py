import numpy as np
import sys
import pandas as pd
import statsmodels.api as sm


sample_GC_file = sys.argv[1]
corr_csv = sys.argv[2]
start_len, end_len = int(sys.argv[3]), int(sys.argv[4])
ref_genome_GC_npy = sys.argv[5]

lower_cnt, upper_cnt = 10, 90
lag = 10


def interpolate_zero_elements(array):
    first_non_zero_index = np.argmax(array != 0)
    last_non_zero_index = len(array) - 1 - np.argmax(array[::-1] != 0)
    inner_array = array[first_non_zero_index:last_non_zero_index + 1]
    non_zero_indices = np.nonzero(inner_array)[0]
    zero_indices = np.where(inner_array == 0)[0]
    inner_array[zero_indices] = np.interp(zero_indices, non_zero_indices, inner_array[non_zero_indices])
    array[first_non_zero_index:last_non_zero_index + 1] = inner_array
    return array

def process_length_level(len_sample):
    processed_sample = np.zeros( len_sample.shape )
    for i in range(len_sample.shape[0]):
        if np.sum(len_sample[i, :])==0:
            processed_sample[i, :] = len_sample[i, :]
        else:
            processed_sample[i, :] = interpolate_zero_elements(len_sample[i, :])
    return processed_sample

def remove_outlier(arr):
    thr = 3
    high = np.percentile(arr, 100-thr)
    low = np.percentile(arr, thr)
    arr = np.where(arr>high, high, arr)
    arr = np.where(arr<low, low, arr)
    return arr

def construct_GC_bias(sample_lengths, ref_lengths, X):
    global lower_cnt, upper_cnt
    sample = np.sum(sample_lengths, axis=0)
    ref = np.sum(ref_lengths, axis=0)
    safe_ref = ref + 1
    GC_array = sample / safe_ref
    
    GC_array = GC_array[lower_cnt: upper_cnt+1]
    GC_array = remove_outlier(GC_array)
    GC_array = GC_array / np.mean(GC_array)
    
    lowess = sm.nonparametric.lowess(GC_array, X, frac=0.1)
    GC_array = lowess[:, 1]
    GC_array = GC_array[lag:-lag]
    GC_array = GC_array / np.mean(GC_array)
    return GC_array

def post_process_outlier(arr):
    min_thr = min(1/20, np.percentile(arr, 1)) # this thr can be 0 or less also
    positive_min = np.min(arr[arr>0])
    final_thr = max(min_thr, positive_min) # ensures +ve GC bias value
    arr[arr<final_thr] = final_thr


len_group_dic = {}
for length in range(start_len, end_len+1):
    start, end = length-2, length+2
    if start<start_len:
        start, end = start_len, start_len+4
    elif end>end_len:
        start, end = end_len-4, end_len
    len_group_dic[length] = [start, end]

ref_GC_array = np.load(ref_genome_GC_npy)
ref_len_GC_array = process_length_level(ref_GC_array)
sample_GC_array = np.load(sample_GC_file)
sample_len_GC_array = process_length_level(sample_GC_array)

final_GC_array = np.zeros((end_len-start_len+1, 61))
correction_factors = np.zeros((end_len-start_len+1, 101))
X = [float(i/100) for i in range(lower_cnt, upper_cnt+1)]
for length in range(start_len, end_len+1):
    start_len_ind = len_group_dic[length][0] - start_len
    end_len_ind = len_group_dic[length][1] - start_len
    if np.sum(sample_len_GC_array[length-start_len, :])==0:
        final_GC_array[length-start_len, :] = np.ones(final_GC_array.shape[1])
    else:
        final_GC_array[length-start_len, :] = construct_GC_bias(sample_len_GC_array[start_len_ind: end_len_ind+1, :], 
                            ref_len_GC_array[start_len_ind: end_len_ind+1, :], X)
post_process_outlier(final_GC_array)

correction_factors[:, lower_cnt+lag: upper_cnt-lag+1] = 1.0/final_GC_array
columns = [str(i) for i in range(0, 101)]
df = pd.DataFrame(correction_factors.tolist())
df.to_csv(corr_csv, index=False, header=columns)