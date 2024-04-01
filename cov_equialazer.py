
import pysam
import numpy as np
from argparse import ArgumentParser
import os
import matplotlib.pyplot as plt

from filter_files import filter_fast5_by_seqids, filter_fastq_by_seqids

parser = ArgumentParser()

parser.add_argument('--bam', type=str, help='Sorted BAM file', required=True)
parser.add_argument('--fast5_dir', type=str, help='Multi fast5 directory', required=True)
parser.add_argument('--fastq', type=str, help='fastq file (single file)', required=True)
parser.add_argument('--cov', type=float, help='required mean coverage', default=60)
parser.add_argument('-o', type=str, help='output dir', required=True)

args = parser.parse_args()


delta = 0.01
STEP = 0
COV_THRESHOLD = args.cov


try:
    os.mkdir(args.o)
except FileExistsError:
    raise FileExistsError('Output folder exists!')


print('Bam file processing...')
samfile = pysam.AlignmentFile(args.bam, 'rb')

read_data = {}

for read in samfile:

    if read.is_mapped == False:
        continue

    if read.reference_name not in read_data:
        read_data[read.reference_name] = []
        
    read_data[read.reference_name].append(
        (read.query_name, read.reference_start, read.reference_end)
    )

max_coords = {
    ref: read_data[ref][-1][-1] for ref in read_data
}


global_selected_reads = set()

for ref in read_data:

    print(f'\nContig {ref} procesing...')


    source_coverage = np.zeros((max_coords[ref]))
    for rec in read_data[ref]:
        read, start, end = rec
        source_coverage[start:end] += 1


    coverage = np.zeros((max_coords[ref]))
    selected_reads = set()

    i = 1

    prev_cov = -1000

    firstskip = 0
    secondskip = 0
    thirdskip = 0

    while np.mean(coverage) < COV_THRESHOLD:
        print(f'Iteration {i}')
        current_coord = -STEP - 1

        for rec in read_data[ref]:
            
            read, start, end = rec
            
            if start <= current_coord + STEP:
                firstskip += 1
                continue
            
            if read in selected_reads:
                secondskip += 1
                continue    
                
            if np.mean(coverage[start:end]) > COV_THRESHOLD*2:
                thirdskip += 1
                current_coord = start
                continue
            
            selected_reads.add(read)
            current_coord = start

            coverage[start:end] += 1
        
        i += 1
        
        next_cov = np.mean(coverage)
        
        print(f'Collected coverage is {round(next_cov,2)}X')
        
        if next_cov - prev_cov < delta:
            print('Complete!')
            break
            
        else:
            prev_cov = next_cov


    print(f'{len(selected_reads)} were selected for contig {ref}.')
    with open(f'{args.o}/selected_reads.{ref}.tmp', 'w') as f:
        for read in selected_reads:
            f.write(f'{read}\n')


    
    source_percantage = len(source_coverage[source_coverage>args.cov])*100/len(coverage)
    source_percantage = round(source_percantage, 2)
    
    percantage = len(coverage[coverage>args.cov])*100/len(coverage)
    percantage = round(percantage, 2)

    print(f'{percantage}% of contig {ref} positions have coverage > {args.cov}X')
    
    fig, axs = plt.subplots(2,1, figsize=(20,11))
    
    axs[0].plot(source_coverage, linewidth=0.5)
    axs[0].hlines(xmin=0, xmax=max_coords[ref], y=args.cov, color='red', linestyle='--')
    axs[0].set_xlabel('position in chromosome')
    axs[0].set_title(f'Input coverage : contig {ref}\n{source_percantage}% positions have coverage > {args.cov}X')
    axs[0].set_ylim(-10, max(coverage))
    
    axs[1].plot(coverage, linewidth=0.5)
    axs[1].hlines(xmin=0, xmax=max_coords[ref], y=args.cov, color='red', linestyle='--')
    axs[1].set_xlabel('position in chromosome')
    axs[1].set_title(f'Cov Equalizer : contig {ref}\n{percantage}% positions have coverage > {args.cov}X')
    axs[1].set_ylim(-10, max(coverage))

    plt.tight_layout()
    plt.savefig(f'{args.o}/cov_plot.{ref}.pdf', format='pdf')

    global_selected_reads = global_selected_reads.union(selected_reads)



print()
print('Filtering fastq files...')
filter_fastq_by_seqids(args.fastq, f'{args.o}/filtered_fastq/', global_selected_reads)


print()
print('Filtering fast5 files...')
filter_fast5_by_seqids(args.fast5_dir, f'{args.o}/filtered_fast5/', global_selected_reads)





