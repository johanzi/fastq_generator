#!/usr/bin/env python

# AUTHOR: Johan Zicola
# DATE: 2018-07-19

import argparse
import os
import random
import sys

'''
	USAGE: This script contains several functions to generate fasta 
        and fastq files. It can generate a random fasta file, generate 
        random fastq files (paired end or single end) of specified length
        and insert size (for paired end). The script can also generate
        single-end and paired-end fastq files based on a specified fasta
        reference file. Coverage and read size can be specified.
'''



#Functions

def generate_DNA(sequence_size):
    seq = []
    while sequence_size != 0:
        nucleotide = random.choice('ACGT')
        seq.append(nucleotide)
        sequence_size -= 1

    return ''.join(seq)


def generate_fasta(name_seq, sequence_size):
    print('>'+name_seq)
    print(generate_DNA(sequence_size))


def reverse_complement(sequence):
    my_dna = str(sequence.upper())
    new_str1 = str.replace(my_dna, 'A','t')
    new_str2 = str.replace(new_str1, 'T','a')
    new_str3 = str.replace(new_str2, 'C','g')
    new_str4 = str.replace(new_str3, 'G','c')
    new_str5 = new_str4.upper()
    new_str6 = new_str5[::-1]
    return new_str6


def generate_fastq(sequence, header="@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:1297 1:N:0:TGACCAAT"):
    # Head contains: Instrument name, run ID, flowcell ID, flowcell lane,
    # tile number within flowcell lane, x-coordinate, y-coordinate,
    # member of a pair (1 or 2), filtered or not (Y or N), control bits, index sequence
    #header = "@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:1297 1:N:0:TGACCAAT"
    # Assuming best quality Phred score for Illumina 1.8+
    quality = len(sequence) * 'I'
    print(header)
    print(sequence)
    print('+')
    print(quality)


def generate_random_fastq_SE(sequence_size, nb_seq):
    index = 1

    while nb_seq > 0:
        sequence = generate_DNA(sequence_size)
        header = "@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:"+str(index)+" 1:N:0:TGACCAAT"
        generate_fastq(sequence, header)
        nb_seq -=1
        index += 1

        
def generate_random_fastq_PE(sequence_size, nb_seq):

    nb_seq_initial = nb_seq
    
    index = 1
    while nb_seq > 0:
        sequence = generate_DNA(sequence_size)
        header = "@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:"+str(index)+"#0/1"
        generate_fastq(sequence, header)
        nb_seq -=1
        index += 1

    #Reset nb_seq
    nb_seq = nb_seq_initial
    
    index = 1
    while nb_seq > 0:
        sequence = generate_DNA(sequence_size)
        header = "@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:"+str(index)+"#0/2"
        generate_fastq(sequence, header)
        nb_seq -=1
        index += 1


def generate_mapped_fastq_SE(ref_fasta, sequence_size, coverage):
    #Retrieve sequences from fasta file
    fasta =  open(ref_fasta,"r")
    fasta = fasta.read()
    fasta = fasta.split("\n")

    list_seq = []

    for seq in fasta:
        if seq[:1] != ">":
            seq = seq.strip()
            list_seq.append(seq)

    #Remove empty line if there is
    list_seq = [_f for _f in list_seq if _f]
    


    #Create all possible kmers for each sequence of the fasta file
    
    kmers = []
    
    while coverage > 0:
        start = 0 + coverage
        for dna in list_seq:
            while start + sequence_size < len(dna):
                    kmer = dna[start:start+sequence_size]
                    kmers.append(kmer)
                    start += sequence_size
        coverage -= 1

        
    index = 1
    for kmer in kmers:
        header = "@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:"+str(index)+" 1:N:0:ATGT"
        generate_fastq(kmer, header)
        index += 1



# Function generates a set of paired-end reads. The first read maps to the
# forward strand and the second on the reverse strand. The output needs to 
# be processed to separate reads 1 and 2 (e.g. with tail and head functions)

def generate_mapped_fastq_PE(ref_fasta, sequence_size, insertion_size, coverage):
    #Retrieve sequences from fasta file
    fasta =  open(ref_fasta,"r")
    fasta = fasta.read()
    fasta = fasta.split("\n")

    list_seq = []

    for seq in fasta:
        if seq[:1] != ">":
            seq = seq.strip()
            list_seq.append(seq)

    #Remove empty line if there is
    list_seq = [_f for _f in list_seq if _f]
    
    # Initialize the 2 lists
    kmers1 = []
    kmers2 = []

    # Get read sequences starting from index coverage-1 so that the last coverage
    # the index 0 (coverage = 1 - 1). Reads containing read1 (kmer1), an insert, 
    # read2 (kmer2). The coverage will switch the reads of 1 bp along the sequence.
    while coverage > 0:
        start = -1 + coverage
        for dna in list_seq:
            while (start + sequence_size*2 + insertion_size) < len(dna):
                kmer1 = dna[start:start+sequence_size]
                kmers1.append(kmer1)
                kmer2 = dna[start+insertion_size+sequence_size:start+insertion_size+sequence_size*2]
                kmer2 = reverse_complement(kmer2)
                kmers2.append(kmer2)
                start += (sequence_size*2)+insertion_size
        coverage -= 1
    
    # To link the paired reads, they need to contain in their name the same
    # Y coordinate (based on observation 35517.end1.fq and 35517.end2.fq headers
    index = 1
    for kmer1 in kmers1:
        header = "@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:"+str(index)+"#0/1"
        generate_fastq(kmer1, header)
        index += 1

    index = 1
    for kmer2 in kmers2:
        header = "@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:"+str(index)+"#0/2"
        generate_fastq(kmer2, header)
        index += 1



def main():
    if len(sys.argv) == 1 or  sys.argv[1] == '-h':
        print('''
    fastq_generator
    
    Work with Python2.7
    
    Author: Johan Zicola (johan.zicola@gmail.com)
    
    Date: 2017-09-09
        
    usage: fastq_generator.py [-h] function

    functions available: generate_fasta, generate_random_fastq_SE, generate_random_fastq_PE, generate_mapped_fastq_SE, generate_mapped_fastq_PE. 

    For more information for each function, enter function name followed by "-h" (e.g.  python fastq_generator.py generate_fasta -h')''')
        sys.exit()
    
    if sys.argv[1] == "generate_fasta":
        #Display help if no argument of "-h"
        if len(sys.argv) == 2 or  sys.argv[2] == '-h':
            print('''
    generate_fasta creates a de novo fasta file of chosen
    name and size. The sequences are generated randomly.
    
    usage: generate_fasta name_seq sequence_size
    
    arguments:
    name_seq = name of the sequence (string)
    sequence_size = length of the sequence in bp (integer)''')
            sys.exit()
        else:
            name_seq = str(sys.argv[2])
        
        #Test if sequence_size argument is given and if integer
        try:
            sys.argv[3] = int(sys.argv[3])
        except IndexError:
            sys.exit("Error: sequence_size missing")
        except ValueError:
            sys.exit("Error: sequence_size should be an integer")
        
        #Call the function
        sequence_size = sys.argv[3]
        generate_fasta(name_seq, sequence_size)
    
    elif sys.argv[1] == "generate_random_fastq_SE":
        if len(sys.argv) == 2 or  sys.argv[2] == '-h':
            print('''
    generate_random_fastq_SE generates random sets of single-end reads
    
    usage: generate_random_fastq_SE sequence_size nb_seq
                     
    arguments:
    sequence_size = length of the sequence in bp (integer)       
    nb_seq = number of sequences required (integer)''') 
            sys.exit()
        else:
            sequence_size = sys.argv[2]
            try:
                sequence_size = int(sequence_size)
            except ValueError:
                sys.exit("Error: sequence_size should be an integer")  
        
        try: 
            sys.argv[3] = int(sys.argv[3])
        except IndexError:
            sys.exit("Error: nb_seq missing")
        except ValueError:
            sys.exit("Error: nb_seq should be an integer")        

        #Call the function
        nb_seq = sys.argv[3]
        generate_random_fastq_SE(sequence_size, nb_seq)
 
    elif sys.argv[1] == "generate_random_fastq_PE":
        if len(sys.argv) == 2 or  sys.argv[2] == '-h':
            print('''
    generate_random_fastq_PE generates random sets of paired-end reads
    
    usage: generate_random_fastq_SE sequence_size nb_seq
                     
    arguments:
    sequence_size = length of the sequence in bp (integer)       
    nb_seq = number of sequences required (integer)''') 
            sys.exit()
        else:
            sequence_size = sys.argv[2]
            try:
                sequence_size = int(sequence_size)
            except ValueError:
                sys.exit("Error: sequence_size should be an integer")  
        
        try: 
            sys.argv[3] = int(sys.argv[3])
        except IndexError:
            sys.exit("Error: nb_seq missing")
        except ValueError:
            sys.exit("Error: nb_seq should be an integer")        

        #Call the function
        nb_seq = sys.argv[3]
        generate_random_fastq_PE(sequence_size, nb_seq)
     
    elif sys.argv[1] == "generate_mapped_fastq_SE":
        if len(sys.argv) == 2 or  sys.argv[2] == '-h':
            print('''
    generate_mapped_fastq_SE generates sets of single-end reads from a reference fasta file
    
    usage: generate_mapped_fastq_SE ref_fasta sequence_size coverage
                     
    arguments:
    ref_fasta = reference file in fasta format (fasta file)
    sequence_size = length of the sequence in bp (integer)
    coverage = number of reads covering each bp (integer)''')
            sys.exit()
        else:
            ref_fasta = sys.argv[2]
        
        try: 
            sys.argv[3] = int(sys.argv[3])
        except IndexError:
            sys.exit("Error: sequence_size missing")
        except ValueError:
            sys.exit("Error: sequence_size should be an integer") 
    
        
        try: 
            sys.argv[4] = int(sys.argv[4])
        except IndexError:
            sys.exit("Error: coverage missing")
        except ValueError:
            sys.exit("Error: coverage should be an integer")
        
        #Call the function
        sequence_size = sys.argv[3]
        coverage = sys.argv[4]
        generate_mapped_fastq_SE(ref_fasta, sequence_size, coverage)
        
    elif sys.argv[1] == "generate_mapped_fastq_PE":
        if len(sys.argv) == 2 or  sys.argv[2] == '-h':
            print('''
    generate_mapped_fastq_PE generates sets of paired-end reads from a reference fasta file
    
    usage: generate_mapped_fastq_PE ref_fasta sequence_size  insertion_size coverage
     
    arguments:
    ref_fasta = reference file in fasta format (fasta file)
    sequence_size = length of the sequence in bp (integer)       
    insertion_size = distance between the paired reads in bp (integer)
    coverage = number of reads covering each bp (integer)''')
            sys.exit()
        else:
            ref_fasta = sys.argv[2]
        
        # test argument sequence_size
        try: 
            sys.argv[3] = int(sys.argv[3])
        except IndexError:
            sys.exit("Error: sequence_size missing")
        except ValueError:
            sys.exit("Error: sequence_size should be an integer") 
        
        # test argument insertion_size
        try: 
            sys.argv[4] = int(sys.argv[4])
        except IndexError:
            sys.exit("Error: insertion_size missing")
        except ValueError:
            sys.exit("Error: insertion_size should be an integer") 
        
        # test argument coverage
        try: 
            sys.argv[5] = int(sys.argv[5])
        except IndexError:
            sys.exit("Error: coverage missing")
        except ValueError:
            sys.exit("Error: coverage should be an integer")         
        
        #Call the function
        sequence_size = sys.argv[3]
        insertion_size = sys.argv[4]
        coverage = sys.argv[5]
        generate_mapped_fastq_PE(ref_fasta, sequence_size, insertion_size, coverage)    
        
        
        
        
    
if __name__ == "__main__":
    sys.exit(main())

    
    
