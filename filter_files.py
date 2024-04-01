import h5py
from os import listdir, mkdir
from tqdm import tqdm
from Bio.SeqIO import parse, write

def filter_fast5_by_seqids(fast5_input_folder, fast5_output_folder, selected_reads):

    target_files = listdir(fast5_input_folder)

    mkdir(fast5_output_folder)

    selected_reads = {'read_' + s for s in selected_reads}

    target_files = [t for t in target_files if '.fast5' in t]

    for f in tqdm(target_files):

        with h5py.File(f'{fast5_input_folder}/{f}', 'r') as input_f5_file:

            with h5py.File(f'{fast5_output_folder}/{f}', 'w') as output_f5_file:
                
                # copy main attributes
                output_f5_file['/'].attrs['file_type'] = input_f5_file['/'].attrs['file_type']
                output_f5_file['/'].attrs['file_version'] = input_f5_file['/'].attrs['file_version']

                for i in input_f5_file.items():
                     
                    if i[0] not in selected_reads:
                        continue

                    output_f5_file.create_group(f'/{i[0]}')
                    input_f5_file.copy(input_f5_file[f'/{i[0]}'], output_f5_file[f'/{i[0]}'], 'Group')



def filter_fastq_by_seqids(fastq_input_folder, fastq_output_folder, selected_reads):
    
    
    target_files = listdir(fastq_input_folder)

    mkdir(fastq_output_folder)

    target_files = [t for t in target_files if '.fq' in t or '.fastq' in t]

    print(target_files)

    for f in tqdm(target_files):

        input_fq_file = parse(f'{fastq_input_folder}/{f}', 'fastq')

        with open(f'{fastq_output_folder}/{f}', 'w') as output_fq_file:

            for rec in input_fq_file:
                if rec.description.split(' ')[0] not in selected_reads:
                    continue
                
                write(rec, output_fq_file, 'fastq')



    