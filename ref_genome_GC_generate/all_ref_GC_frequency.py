import subprocess, os, argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--ref', type=str, required=True, help='Reference genome path')
parser.add_argument('--output_folder', type=str, required=True, help='Folder where expected GC content frequencies for different fragment lengths will be stored')
parser.add_argument('--start_length', type=int, default=51, help='Smallest fragment length of interest')
parser.add_argument('--end_length', type=int, default=400, help='Largest fragment length of interest')
parser.add_argument('--CPU', type=int, default=32, help='Number of CPU to use for parallel processing')
parser.add_argument('--regions', type=str, required=True, help='Genome region list based on which expected GC content frequencies will be calculated')

args = parser.parse_args()
ref_genome_path = args.ref
output_folder = args.output_folder 
start_len, end_len = args.start_length, args.end_length
CPU = args.CPU
bin_location_file = args.regions

for i in range(start_len, end_len+1):
    length = i
    command = f'python single_ref_GC_frequency.py {length} {bin_location_file} {ref_genome_path} {output_folder} {CPU}'
    subprocess.run(command, shell=True)

final_GC_array = np.zeros((end_len-start_len+1, 101))
for i in range(start_len, end_len+1):
    final_GC_array[i-start_len, :] = np.load(f'{output_folder}/{i}.npy')
    os.remove(f'{output_folder}/{i}.npy')
np.save(f'{output_folder}/ref_genome_GC.npy', final_GC_array)
