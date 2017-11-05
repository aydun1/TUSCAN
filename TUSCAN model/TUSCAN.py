#!/usr/bin/env python

#Given a sequence, this program uses a Random Forest to predict the activity 
#via a feature matrix, using a predefined master model created using a training set

# Can use a bed file, fasta file, or plain list of sequences

import sys
import FeatureMatrix
import pickle
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.externals import joblib
import numpy
import pybedtools
import argparse
from collections import OrderedDict
import os

#Valid DNA letters
VALID = 'ACTG'
#collect directory
directory = os.path.dirname(os.path.realpath(__file__))

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def reverse_complement(sequence):
    rc_map = {'A':'T','C':'G','G':'C','T':'A'}
    return ''.join([rc_map[B] for B in sequence][::-1])


def preread_bed(input_name, genome):
    with open(input_name, 'r') as f:
        bed_line = f.readline()
    bed_type = bed_line.split()
    bed_columns = len(bed_type)
    c = pybedtools.BedTool(input_name)

    if bed_columns == 3:
        has_name = False
        bed_handle = c.sequence(fi=genome)
    elif bed_columns > 3:
        has_name = True
        bed_handle = c.sequence(name=True, fi=genome)
    else:
        eprint('Invalid bed file, must have at least 3 columns')
        sys.exit()

    return has_name, bed_handle


#get important features
def get_important_features(model_type):
    if model_type == 'Regression':
        features_file = os.path.join(directory, 'rf_features_regression.txt')
    elif model_type == 'Classification':
        features_file = os.path.join(directory, 'rf_features_classification.txt')
    else:
        eprint('Invalid model type, must be Classification or Regression')
        sys.exit()

    with open(features_file, 'r') as f:
        important_features = [b.strip('"') for line in f for b in line.split()]

    return important_features


#Store sequences for recall later
def read_bed_or_fa(input_name, input_type, has_name):
    if input_type == 'bed':
        file_name = d.seqfn
    else:
        file_name = input_name
    f = open(file_name, 'r')
    s = OrderedDict()
    z = 1
    for line in f:
        line = line.rstrip()
        if line[0] == '>':
            if has_name:
                k = line[1:]
            else:
                k = (z+1)/2
        else:
            line = line.upper()
            if not all(i in VALID for i in line):
                count = z
                if (input_type == 'bed'):
                    count = z/2
                eprint('The sequence at line {} contains an invalid nucleotide'.format(count))
            elif len(line) != 30:
                count = z
                if (input_type == 'bed'):
                    count = z/2
                eprint('The sequence at line {} is not 30 base pairs, it is {}'.format(count, len(line)))
            elif line[25:27] != 'GG':
                reverseLine = reverse_complement(line)
                if reverseLine[25:27] != 'GG':
                    count = z
                    if (input_type == 'bed'):
                        count = z/2
                    eprint('The sequence at line {} does not have a PAM motif'.format(count))
                else:
                    s[k] = {}
                    s[k]['seq'] = reverseLine
                    s[k]['dir'] = '-'
                    count = z
                    if (input_type == 'bed'):
                        count = z/2
                    eprint('The sequence at line {} was in the negative orientation'.format(count))
            else:
                s[k] = {}
                s[k]['seq'] = line
                s[k]['dir'] = '+'
        z += 1
    f.close()
    return s


def read_txt(input_name, input_type):
    s = OrderedDict()
    count = 1
    f = open(input_name, 'r')
    for line in f:
        line = line.rstrip()
        line = line.upper()
        if not all(i in VALID for i in line):
            eprint('The sequence at line {} contains an invalid nucleotide'.format(count))
        elif len(line) != 30:
            eprint('The sequence at line {} is not 30 base pairs, it is {}'.format(count, len(line)))
        elif line[25:27] != 'GG':
            reverseLine = reverse_complement(line)
            if reverseLine[25:27] != 'GG':
                eprint('The sequence at line {} does not have a PAM motif'.format(count))
            else:
                s[count] = {}
                s[count]['seq'] = reverseLine
                s[count]['dir'] = '-'
                eprint('The sequence at line {} was in the negative orientation'.format(count))
        else:        
            s[count] = {}    
            s[count]['seq'] = line
            s[count]['dir'] = '+'
            count += 1
    f.close()
    return s


#Read in sequence data to create feature matrix
#generate feature matrix
def create_feature_matrix(parsed_input_file, matrix_out):
    stdout_ = sys.stdout
    sys.stdout = open(matrix_out, 'w')
    FeatureMatrix.main([parsed_input_file])
    sys.stdout.close()
    sys.stdout = stdout_


#grabs list of features
def get_feature_list(important_features, matrix_out):
    with open(matrix_out, 'r') as f:
        features = f.readline()
    features = features.split()
    num_features = len(features)
    features = features[1:]

    #gets index of important features
    a = [features.index(i) for i in important_features]

    data = numpy.genfromtxt(matrix_out, dtype = 'f8', skip_header = 1, usecols = range(1, num_features))

    if len(data.shape) == 1:
        data = numpy.array([data])

    return data[:, a]


def run_regressor(feature_list):
    rfm = os.path.join(directory, 'rfModelregressor.joblib')
    with open(rfm, 'rb') as f:
        rf = joblib.load(f)
    scores = rf.predict(feature_list)
    return scores


def run_classifier(feature_list):
    cfm = os.path.join(directory, 'rfModelclassifier.joblib')
    with open(cfm, 'rb') as f:
        rf = joblib.load(f)
    scores = rf.predict(feature_list)
    return scores


def create_output(scores, parsed_input_file, output_file):
    layout = '{!s:50} {!s:31} {!s:15} {!s:3}\n'
    header = layout.format('ID', 'Sequence', 'Score', 'Dir')
    if output_file:
        with open(output_file, 'w') as f:
            f.write(header)
            for idx, a in enumerate(parsed_input_file):
                f.write(layout.format(a, parsed_input_file[a]['seq'], scores[idx], parsed_input_file[a]['dir']))
    else:
        print(header)
        for idx, a in enumerate(parsed_input_file):
            print(layout.format(a, parsed_input_file[a]['seq'], scores[idx], parsed_input_file[a]['dir']))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', required=False, help='Genome')
    parser.add_argument('-o', required=False, help='Output file')
    parser.add_argument('-i', required=True, help='Input file')
    parser.add_argument('-t', required=True, help='Type of input file (bed/txt/fa)')
    parser.add_argument('-m', required=True, help='Type of Model: (Regression/Classification)')

    args = parser.parse_args()

    #arguments to pass to matrix builder
    matrix_out = str(args.i[:-len(str(args.t))-1]) + '_matrix.txt'

    if (args.t == 'bed'):
        if not args.g:
            print('Genome required with bed file')
            sys.exit()
        has_name, d = preread_bed(args.i, args.g)
    else:
        has_name = True

    important_features = get_important_features(args.m)

    if (args.t == 'bed' or args.t == 'fa'):
        parsed_input_file = read_bed_or_fa(args.i, args.t, has_name)
    elif (args.t == 'txt'):
        parsed_input_file = read_txt(args.i, args.t)
    else:
        print('Usage [-t bed/fa/txt]')
        sys.exit()

    create_feature_matrix(parsed_input_file, matrix_out)

    feature_list = get_feature_list(important_features, matrix_out)

    if args.m == 'Regression':
        scores = run_regressor(feature_list)
    elif args.m == 'Classification':
        scores = run_classifier(feature_list)
    else:
        print('Usage: [-m Regression] or [-m Classification]')

    create_output(scores, parsed_input_file, args.o)


if __name__== '__main__':
    main()
