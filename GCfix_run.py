import subprocess, os, time, argparse, shutil

parser = argparse.ArgumentParser()
parser.add_argument('--input_folder', type=str, required=True, help='Input folder path full of bam files to be corrected')
parser.add_argument('--output_folder', type=str, required=True, help='Output folder where the correction factors will be created')
parser.add_argument('--MAPQ', type=int, default=30, help='Mapping quality of reads used for GC correction')
parser.add_argument('--start_length', type=int, default=51, help='Smallest fragment length of interest')
parser.add_argument('--end_length', type=int, default=400, help='Largest fragment length of interest')
parser.add_argument('--CPU', type=int, default=32, help='Number of CPU to use for parallel processing')
parser.add_argument('--ref', type=str, default='hg19__GRCh37/ref.fa', help='Reference genome to use for GC correction factor estimation')
parser.add_argument('--ref_genome_GC', type=str, default='hg19__GRCh37/ref_genome_GC.npy', help='Expected GC count frequencies from reference genome (.npy file path)')
parser.add_argument('--GC_estimation_regions', type=str, default='hg19__GRCh37/correction_factor_estimation_bin_locations.csv', help='Genome region list based on which GC bias will be estimated')
parser.add_argument('--GC_tagging_flag', type=int, default=0, help='Whether to create tagged bam files with GC correction factor (1 for tagging)')
parser.add_argument('--GC_tagging_regions', type=str, default='hg19__GRCh37/GC_tagging_bin_locations.csv', help='Genome region list where the reads need to be tagged with their GC correction factors')
parser.add_argument('--GC_tagging_folder', type=str, default='Sample_Output/Tagged_Bams/', help='Output folder to create tagged bam files')
parser.add_argument('--coverage_flag', type=int, default=0, help='Whether to create corrected coverage profile or not')
parser.add_argument('--coverage_regions', type=str, default='hg19__GRCh37/correction_factor_estimation_bin_locations.csv', help='Genome regions/bins where the corrected coverage is desired')
parser.add_argument('--coverage_folder', type=str, default='Sample_Output/Coverage_Profiles/', help='Output folder where the corrected coverage profiles will be created for the input samples')

args = parser.parse_args()
input_folder, output_folder = args.input_folder, args.output_folder 
mapq = args.MAPQ
start_len, end_len = args.start_length, args.end_length
CPU = args.CPU
ref_genome_path = args.ref
ref_genome_GC_npy = args.ref_genome_GC
GC_estimation_region_file = args.GC_estimation_regions
GC_flag = args.GC_tagging_flag
GC_tagging_region_file = args.GC_tagging_regions
GC_tag_folder = args.GC_tagging_folder
cov_flag = args.coverage_flag
cov_region_file = args.coverage_regions
cov_folder = args.coverage_folder

for file_ in os.listdir(input_folder):
    if file_.endswith('.bam'):
        start = time.time()
        print(file_[:-4])

        input = f'{input_folder}/{file_}'
        output_npy = f'{output_folder}/{file_[:-4]}.npy'
        command = f'python single_sample_GC_count.py {input} {output_npy} {mapq} {start_len} {end_len} {CPU} {ref_genome_path} {GC_estimation_region_file}'
        subprocess.run(command, shell=True)

        input = output_npy
        output_csv = f'{output_folder}/{file_[:-4]}__correction_factors.csv'
        command = f'python single_sample_correction_factor.py {input} {output_csv} {start_len} {end_len} {ref_genome_GC_npy}'
        subprocess.run(command, shell=True)
        os.remove(input)

        if GC_flag==1:
            if os.path.isdir(GC_tag_folder)==False:
                os.mkdir(GC_tag_folder)
            input_bam = f'{input_folder}/{file_}'
            output_bam = f'{GC_tag_folder}/{file_}'
            correction_path = output_csv
            temp_folder = f'{GC_tag_folder}/temp'
            command = f'python single_sample_GC_tagging.py {input_bam} {output_bam} {correction_path} {mapq} {start_len} {end_len} {CPU} {ref_genome_path} {GC_tagging_region_file} {temp_folder}'
            subprocess.run(command, shell=True)

        if cov_flag==1:
            if os.path.isdir(cov_folder)==False:
                os.mkdir(cov_folder)
            input_bam = f'{input_folder}/{file_}'
            output_cov = f'{cov_folder}/{file_[:-4]}__corrected_coverage.csv'
            correction_path = output_csv
            command = f'python single_sample_corrected_cov.py {input_bam} {output_cov} {correction_path} {mapq} {start_len} {end_len} {CPU} {ref_genome_path} {cov_region_file}'
            subprocess.run(command, shell=True)

        end = time.time()
        print( f'{round((end - start)/60, 3)} minutes' )
        print()
