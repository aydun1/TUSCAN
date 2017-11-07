#!/usr/bin/env python

#This program creates a feature matrix given a fasta (.fa) and corresponding bed (.bed) file
#Written by Daniel Reti, November 2016

#USAGE: ./this.program.py file.fa [file.bed]
#If no bed file, then score is not in present in matrix
#If both fa and bed provided, then only sequences with both score and appear in .fa will be in matrix
#TODO: May need to swap around score and sign in bed split - see TODO below

from __future__ import print_function
import sys

#determines gc content of given sequence
def gc(seq):
    l = len(seq)
    n = 0
    for r in seq:
        if r == 'C' or r == 'G':
            n += 1
    return str(round(float(n)/l * 100, 2))
        
#determines number of a given base in a given sequence
def content(seq, base):
    n = 0
    for r in seq:
        if r == base:
            n += 1
    return str(n)

def di_content(seq, di):
    the_list = []
    for i in di:
        count = 0
        for j in range(len(seq)-1):
            if (seq[j:j+2] == i):
                count += 1
        the_list.append(str(count))
    return the_list

def pam(seq, di, bases):
    l = [0] * 16
    for i in di:
        multiply = bases.index(seq[24])
        offset = bases.index(seq[27])
        l[multiply * 4 + offset] = 1
    the_list = [str(n) for n in l]
    return the_list
        

#Nucleotide distribution across sequence
#Position 1 is 5' end, Position 21 is the N in the NGG PAM (3' end)
#Sequence structure: 5' 20merNGG 3'
#C7 means a C at position 7
def nucleotide(seq, bases):
    l = [0] * (len(seq)) * 4
    for i in range(len(seq)):
        count = 4 * i + bases.index(seq[i]) 
        l[count] = 1
    the_list = [str(n) for n in l]
    return the_list
    
#Dinucleotide distribution accross sequence
#CT7 means C at 7 and T at 8 
#structure outlined above 
def dinucleotide(seq, di):
    #-1 is since a sequence of length N has N-1 dinucleotides
    l = [0] * ((len(seq))-1) * 16
    for i in range((len(seq)) - 1):
        l[16 * i + di.index(seq[i:i+2])] = 1
    the_list = [str(n) for n in l]
    return the_list
        
#global variables

def main(args):
    bases = ['A', 'C', 'G', 'T']
    # This can be adjusted if the sequence length changes down the track
    seqLength = 30
    di = [a + b for a in bases for b in bases]
    d = args[0]

    #reads bed file (arg2 in cmd line)
    if (len(args) > 1):
        g = open(args[1], 'r')
        for line in g:
            #TODO: may need to swap around sign and score as required
            #if score last
            #(ch, start, end, name, sign, score) = line.split()
            #if sign last
            (ch, start, end, name, score, sign) = line.split()
            if name in d:
                d[name]['score'] = score
        g.close()

    #view feature matrix
    #headers

    h1 = ['Name', 'GC_', 'A', 'C', 'G', 'T']

    #with bed file (scores)
    if (len(args) > 1):
        h1.append('Activity')

    # range is from 1 up to sequence length (+1 accounts for shift from zero)
    h2 = ['{}{}'.format(j, i) for i in range(1, seqLength + 1) for j in bases]
    h3 = ['{}{}'.format(j, i) for i in range(1, seqLength) for j in di]
    h4 = ['{}'.format(di[i]) for i in range(len(di))]
    h5 = ['{}GG{}'.format(di[i][0], di[i][1]) for i in range(len(di))]

    print(' '.join(h1 + h2 + h3 + h4 + h5))

    #matrix data with bed file
    if (len(args) > 1):
        for a in d:
            if 'score' in d[a]:
                seq = d[a]['seq']
                a_list = [a, gc(seq), content(seq, 'A'), content(seq, 'C'), content(seq, 'G'), content(seq, 'T'), d[a]['score']] + \
                nucleotide(seq, bases) + \
                dinucleotide(seq, di) + \
                di_content(seq, di) + \
                pam(seq, di, bases)
                print(' '.join(a_list))
    else:
    #matrix data without bed file
        for a in d:
            seq = d[a]['seq']
            a_list = [a, gc(seq), content(seq, 'A'), content(seq, 'C'), content(seq, 'G'), content(seq, 'T')] + \
            nucleotide(seq, bases) + \
            dinucleotide(seq, di) + \
            di_content(seq, di) + \
            pam(seq, di, bases)
            print(' '.join(a_list))


if __name__ == '__main__':
    main(sys.argv)
