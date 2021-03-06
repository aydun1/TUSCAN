#!/usr/bin/env python

#Given a sequence, this program uses a Random Forest to predict the activity
#via a feature matrix, using a predefined master model created using a training set

# USAGE:
# Users can use the TUSCAN model in 3 ways:
# 1. Only supply a fasta file. The entire sequence will be analysed. Output will show: Chrom: N/A. Candidate site locations will be 0-based.
# 2. Supply a fasta file with known chrom/position information. The entire sequence will be analysed. Output will show supplied chromosome info, and position relative to supplied start/finish positions. 
# 3. Supply a fasta file, include chrom/position information for a region to extract. Include the -e flag with no arguments. The region of interest (start pos -> end pos) will be extracted and analysed. Output will reflect start position relative to supplied start/end position and chrom info.

# Users can also supply a thread number for multithreading using -t. If none supplied, default thread depending on available cores will be set. (-t [THREADS (int)])
# Users can also supply an output file name (recommended), otherwise output is written to TUSCAN_output.txt and overwritten upon each iteration (-o [FILENAME])
# Users must also specify whether regression or classification is required (-m [Regression/Classification])

# NOTE: Final output will be in arbitrary order due to multithreading. 
# Output is tab-separated and easy for users to manipulate using UNIX sort or other flavours. 

import argparse
import os
import re
import sys

from collections import namedtuple
from multiprocessing import Process, Queue, cpu_count
from sklearn.externals import joblib

NucleotideLocation = namedtuple('NucleotideLocation', ['nucleotide', 'location'])
DinucleotideLocation = namedtuple('DinucleotideLocation', ['dinucleotide', 'location'])

REGRESSION_NUCLEOTIDES_OF_INTEREST = (
    NucleotideLocation(nucleotide='T', location=4),
    NucleotideLocation(nucleotide='C', location=7),
    NucleotideLocation(nucleotide='C', location=10),
    NucleotideLocation(nucleotide='T', location=17),
    NucleotideLocation(nucleotide='C', location=20),
    NucleotideLocation(nucleotide='T', location=20),
    NucleotideLocation(nucleotide='G', location=21),
    NucleotideLocation(nucleotide='T', location=21),
    NucleotideLocation(nucleotide='G', location=22),
    NucleotideLocation(nucleotide='T', location=22),
    NucleotideLocation(nucleotide='C', location=24),
    NucleotideLocation(nucleotide='G', location=24),
    NucleotideLocation(nucleotide='T', location=24)
)

CLASSIFICATION_NUCLEOTIDES_OF_INTEREST = (
    NucleotideLocation(nucleotide='G', location=5),
    NucleotideLocation(nucleotide='T', location=11),
    NucleotideLocation(nucleotide='C', location=12),
    NucleotideLocation(nucleotide='A', location=16),
    NucleotideLocation(nucleotide='T', location=17),
    NucleotideLocation(nucleotide='C', location=20),
    NucleotideLocation(nucleotide='T', location=20),
    NucleotideLocation(nucleotide='T', location=22),
    NucleotideLocation(nucleotide='T', location=23),
    NucleotideLocation(nucleotide='C', location=24),
    NucleotideLocation(nucleotide='G', location=24),
    NucleotideLocation(nucleotide='T', location=24),
)

REGRESSION_DINUCLEOTIDES_OF_INTEREST = (
    DinucleotideLocation(dinucleotide='AC', location=1),
    DinucleotideLocation(dinucleotide='AC', location=2),
    DinucleotideLocation(dinucleotide='CA', location=3),
    DinucleotideLocation(dinucleotide='TT', location=4),
    DinucleotideLocation(dinucleotide='GA', location=5),
    DinucleotideLocation(dinucleotide='CT', location=6),
    DinucleotideLocation(dinucleotide='AC', location=8),
    DinucleotideLocation(dinucleotide='CC', location=8),
    DinucleotideLocation(dinucleotide='GA', location=8),
    DinucleotideLocation(dinucleotide='TT', location=9),
    DinucleotideLocation(dinucleotide='AT', location=10),
    DinucleotideLocation(dinucleotide='CG', location=11),
    DinucleotideLocation(dinucleotide='GA', location=12),
    DinucleotideLocation(dinucleotide='CC', location=14),
    DinucleotideLocation(dinucleotide='GA', location=15),
    DinucleotideLocation(dinucleotide='CC', location=16),
    DinucleotideLocation(dinucleotide='GG', location=16),
    DinucleotideLocation(dinucleotide='TT', location=16),
    DinucleotideLocation(dinucleotide='CT', location=17),
    DinucleotideLocation(dinucleotide='AA', location=18),
    DinucleotideLocation(dinucleotide='GG', location=19),
    DinucleotideLocation(dinucleotide='AT', location=20),
    DinucleotideLocation(dinucleotide='CC', location=20),
    DinucleotideLocation(dinucleotide='CG', location=20),
    DinucleotideLocation(dinucleotide='CT', location=20),
    DinucleotideLocation(dinucleotide='GG', location=20),
    DinucleotideLocation(dinucleotide='TA', location=21),
    DinucleotideLocation(dinucleotide='TG', location=21),
    DinucleotideLocation(dinucleotide='CC', location=22),
    DinucleotideLocation(dinucleotide='GA', location=22),
    DinucleotideLocation(dinucleotide='TA', location=22),
    DinucleotideLocation(dinucleotide='CG', location=23),
    DinucleotideLocation(dinucleotide='GA', location=23),
    DinucleotideLocation(dinucleotide='GG', location=23),
    DinucleotideLocation(dinucleotide='TG', location=23),
    DinucleotideLocation(dinucleotide='GA', location=24),
    DinucleotideLocation(dinucleotide='GT', location=24),
    DinucleotideLocation(dinucleotide='TC', location=24),
)

CLASSIFICATION_DINUCLEOTIDES_OF_INTEREST = (
    DinucleotideLocation(dinucleotide='CG', location=11),
    DinucleotideLocation(dinucleotide='GA', location=15),
    DinucleotideLocation(dinucleotide='TT', location=16),
    DinucleotideLocation(dinucleotide='CC', location=20),
    DinucleotideLocation(dinucleotide='TA', location=22),
    DinucleotideLocation(dinucleotide='CG', location=23),
    DinucleotideLocation(dinucleotide='TC', location=23),
    DinucleotideLocation(dinucleotide='TG', location=23),
    DinucleotideLocation(dinucleotide='CC', location=24),
    DinucleotideLocation(dinucleotide='GA', location=24),
    DinucleotideLocation(dinucleotide='GC', location=24),
    DinucleotideLocation(dinucleotide='GT', location=24),
    DinucleotideLocation(dinucleotide='TC', location=24),
)

CLASSIFICATION_DINUCLEOTIDES = [
    'AA',
    'AC',
    'AG',
    'AT',
    'CA',
    'CC',
    'CG',
    'CT',
    'GA',
    'GC',
    'GG',
    'GT',
    'TA',
    'TC',
    'TG',
    'TT',
]

REGRESSION_DINUCLEOTIDES = [
    'CA',
    'CT',
    'GC',
    'TC',
    'TG',
    'TT',
]


#determines gc content of given sequence
def gc(seq):
    return round((seq.count('C') + seq.count('G')) / len(seq) * 100, 2)


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement[base] for base in reversed(sequence))


#counts appearance of dinucleotides in sequence
def di_content(seq, dinucleotides_to_count, features, start_index):
    for idx, dinucleotide in enumerate(dinucleotides_to_count):
        count = start = 0
        while True:
            start = seq.find(dinucleotide, start) + 1
            if start:
                count += 1
            else:
                features[start_index+idx] = count
                break


#checks if specific PAM is present in sequence
def pam(seq, features, index):
    if seq[24:28] == 'TGGT':
        features[index] = 1


#checks if given position-specific nucleotides are present in sequence
def nucleotide(seq, nucleotides_of_interest, features, start_index):
    for idx, nucleotide_loc in enumerate(nucleotides_of_interest):
        if seq[nucleotide_loc.location-1] == nucleotide_loc.nucleotide:
            features[start_index+idx] = 1


#checks if given position-specific dinucleotides are present in sequence
def dinucleotide(seq, dinucleotides_of_interest, features, start_index):
    #-1 is since a sequence of length N has N-1 dinucleotides
    for idx, dinucleotides_loc in enumerate(dinucleotides_of_interest):
        location = dinucleotides_loc.location
        if seq[location-1:location+1] == dinucleotides_loc.dinucleotide:
            features[start_index+idx] = 1


#generates a feature vector from a given 30 nucleotide sequence
def get_features(seq):
    if is_regression:
        features = [0] * 63
        features[0] = gc(seq)
        features[1] = seq.count('A')
        features[2] = seq.count('C')
        features[3] = seq.count('G')
        features[4] = seq.count('T')
        di_content(seq, REGRESSION_DINUCLEOTIDES, features, 5)
        nucleotide(seq, REGRESSION_NUCLEOTIDES_OF_INTEREST, features, 11)
        dinucleotide(seq, REGRESSION_DINUCLEOTIDES_OF_INTEREST, features, 24)
        pam(seq, features, 62)
    else:
        features = [0] * 46
        features[0] = gc(seq)
        features[1] = seq.count('A')
        features[2] = seq.count('C')
        features[3] = seq.count('G')
        features[4] = seq.count('T')
        di_content(seq, CLASSIFICATION_DINUCLEOTIDES, features, 5)
        nucleotide(seq, CLASSIFICATION_NUCLEOTIDES_OF_INTEREST, features, 21)
        dinucleotide(seq, CLASSIFICATION_DINUCLEOTIDES_OF_INTEREST, features, 33)
    return features


def output_sequences(sequences, metadata, output_queue):
    scores = rf.predict([get_features(s) for s in sequences])
    for m in zip(metadata, scores):
        output_queue.put(m[0] + [m[1]])


def score_sequences(matches_queue, output_queue, is_reverse):
    strand = '-' if is_reverse else '+'
    sequences = []
    metadata = []
    while True:
        match_start, sequence = matches_queue.get()
        if match_start is None:
            if sequences:
                output_sequences(sequences, metadata, output_queue)
            output_queue.put(None)
            break
        if is_reverse:
            sequence_end_pos = end - match_start - 3 - 1
            sequence_start_pos = sequence_end_pos - 23 + 1 
        else:
            sequence_start_pos = start + match_start + 3 + 2
            sequence_end_pos = sequence_start_pos + 23 - 1

        sequences.append(sequence)
        metadata.append([chrom, sequence_start_pos, sequence_end_pos, strand, sequence[4:-3]])
        if len(sequences) >= 10000:
            output_sequences(sequences, metadata, output_queue)
            sequences = []
            metadata = []


def fill_queue(matches, matches_queue):
    for m in matches:
        matches_queue.put((m.start(), m.group(1)), block=True)
    for x in range(num_threads):
        matches_queue.put((None, None), block=True)


if __name__ == '__main__':
    try:
        num_cores = cpu_count()
    except NotImplementedError:
        num_cores = 4

    #Parse Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', required=False, help='Genome')
    parser.add_argument('-o', required=False, help='Output file')
    parser.add_argument('-s', required=False, help='Start Position')
    parser.add_argument('-f', required=False, help='End Position')
    parser.add_argument('-c', required=False, help='Chromosome')
    parser.add_argument('-m', required=True, help='Type of Model: (Regression/Classification)')
    parser.add_argument('-t', required=False, default=num_cores, help='Number of threads (default: 4)')
    parser.add_argument('-e', required=False, action='store_true', help='If you want to excise a region inside a supplied genome, include this flag with no arguments, exclude it to analyse the full genome')

    args = parser.parse_args()

    sequence = ''
    is_regression = False
    if args.m == 'Regression':
        is_regression = True
    elif args.m != 'Classification':
        print('Please specify Regression or Classification')
        sys.exit()

    chrom = args.c if args.c else 'N/A'
    genome = args.g
    output_file = args.o if args.o else 'TUSCAN_output.txt'
    num_threads = args.t
    extract = args.e
    #if a chromosome, genome, start and stop positions are given
    if extract:
        import pybedtools
        print('Extracting region')
        start = args.s
        end = args.f
        #turn region information into string and write to BED file
        region = '{}\t{}\t{}'.format(chrom, start, end)
        with open('customRegion.bed', 'w') as f:
            f.write(region)
        #extract region of interest from genome
        c = pybedtools.BedTool('customRegion.bed')
        d = c.sequence(fi = genome)
        with open(d.seqfn, 'r') as g:
            sequence = [line.rstrip().upper() for line in g if not line.startswith('>')][-1]
        os.remove('customRegion.bed')
    #else if a FASTA sequence has been nominated
    elif genome:
        #extract from file:
        print('Analysing entire supplied sequence.\nIf you wish to analyse a sub-region, supply start, end and chromosome flags and include the -e flag')
        with open(genome, 'r') as f:
            sequence = ''.join(line.rstrip().upper() for line in f if not line.startswith('>'))
        start = int(args.s) - 1 if args.s else 0
        end = int(args.f) - 1 if args.f else len(sequence)
    else:
        print('Please supply a sequence to be analysed')
        sys.exit()

    if not sequence:
        print('If using the extract -e flag, please supply ALL OF THESE and make sure they are correct: start, end AND chromosome flags')
        sys.exit() 

    #Find and store all sequences + location information
    guide_re = re.compile(r'(?=([ACTG]{25}GG[ACTG]{3}))')
    matches = guide_re.finditer(sequence)
    matches_rev = guide_re.finditer(reverse_complement(sequence))
    
    # Load the appropriate model
    file_path = os.path.dirname(os.path.realpath(__file__))
    file_name = 'rfModelregressor.joblib' if is_regression else 'rfModelclassifier.joblib'
    rf = joblib.load(os.path.join(file_path, file_name))

    output_queue = Queue()
    matches_queue = Queue(maxsize=num_threads*2)

    print('Analysing')
    with open(output_file, 'w') as f:
        layout = '{!s:5} {!s:10} {!s:10} {!s:8} {!s:34} {!s:15}\n'
        f.write(layout.format('Chrom', 'Start', 'End', 'Strand', 'Candidate_sgRNA', 'TUSCAN_Score'))
        for is_reverse, match_type in enumerate((matches, matches_rev)):
            matches_process = Process(target=fill_queue, args=(match_type, matches_queue))
            matches_process.start()
            processes = [Process(target=score_sequences, args=(matches_queue, output_queue, is_reverse)) for x in range(num_threads)]
            for p in processes:
                p.start()
            num_done = 0
            while num_done < num_threads:
                output = output_queue.get()
                if output is None:
                    num_done += 1
                else:
                    f.write(layout.format(*output))
            matches_process.join()
            for process in processes:
                process.join()
