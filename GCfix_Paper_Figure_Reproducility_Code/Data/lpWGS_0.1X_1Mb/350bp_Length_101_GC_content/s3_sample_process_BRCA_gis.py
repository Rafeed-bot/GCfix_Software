import subprocess

dest_dir = '/rfs-storageservice/GIS/Projects/LOCCG/rafeed/Deep_Samples/BRCA_gis/'
directory = 'Segmentation_Project/deep_brca_gis_rmdup/'
bucket = 's3://fragment-share.store.genome.sg/' + directory

file_names = ['E10c', 'E2c', 'E6c', 'E7c', 'E8c', 'D14', 'D19', 'D23', 'D7', 'D9']
files = []
for file_name in file_names: 
    files.append(file_name + '.bam')
    files.append(file_name + '.bam.bai')

for file_ in files:
    print(file_)
    # filePath = directory + file_
    # command = 'aws s3api restore-object --bucket fragment-share.store.genome.sg --key {filePath} --restore-request \'{{\"Days\":10,\"GlacierJobParameters\":{{\"Tier\":\"Standard\"}}}}\' --profile peter'.format(filePath=filePath)
    # subprocess.run(command, shell=True)

    filePath = bucket + file_
    command = 'aws s3 cp {filePath} {dest_dir} --profile peter'.format(filePath=filePath, dest_dir=dest_dir)
    subprocess.run(command, shell=True)
