# fastq_generator
A suite of tools to generate random or mapped fastq files.
The script is compatible with python2.7 and python3.4.

# Table of content
- [fastq_generator](#fastq-generator)
- [Table of content](#table-of-content)
  * [Usage](#usage)
  * [Generate random fasta files](#generate-random-fasta-files)
  * [Generate random fastq files with single-end reads](#generate-random-fastq-files-with-single-end-reads)
  * [Generate random fastq files with paired-end reads](#generate-random-fastq-files-with-paired-end-reads)
  * [Generate fastq files SE mapping on a specified fasta reference](#generate-fastq-files-se-mapping-on-a-specified-fasta-reference)
  * [Generate fastq files PE mapping on a specified fasta reference](#generate-fastq-files-pe-mapping-on-a-specified-fasta-reference)
  * [Generate reads mapping to the reverse strand](#generate-reads-mapping-to-the-reverse-strand)
  * [Authors](#authors)
  * [License](#license)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>

## Usage
The different functions of fastq_generator are located within one python script.
One can call the script in python without argument or with '-h' to get the help message.

## Generate random fasta files

generate_fasta creates a de novo fasta file of chosen
name and size. The sequences are generated randomly.

usage: generate_fasta name_seq sequence_size

arguments:
name_seq = name of the sequence (string)
sequence_size = length of the sequence in bp (integer)

```
python fastq_generator.py generate_fasta Chr1 1000 > chr1.fa

cat chr1.fa

>Chr1
TGGGATCTCTAAGCACTTGTACTCGATGCATCTCCACGATTACCGTTGCCATAGCTGCGAAACGATATAAAGGGTTTGGCTGTGCATGGTGTTCTCTTCAACGCCTGTACCAAACTTTCTACAATTGTCGGTGACCTCTTCCGGAAGCTGAGAATCGCTAATGAGACAAACGGCTTACTCCCTGAGTAGTACGACGTCCAAGCGAAACTATCACCACGGACTTTTTGTGTATAACTGCGTTTATCTAAAACAAGGTCGTGAGGCCTGATTGCGCACTTACGTCTGGGCCAAGACCGTGTGGACGTCTCCAACCAACACTGCTCGATGCTGCCTTGAAGATCTCGCCATCTCGAGTCAGTAACATCCCAGTGGAACACTTGAGACAGGGAAATCAGGCGCGTCTTCACGACAGGACTCTCCCCAGAGGTCCTCACTGAAGAGGTGCTATCGGTCAACCGTTCGTGGGGCCCGTTAGGCGTCTGATCTACCCAGGTCATCTTGAGACTCCCGACACCCGCATTAATTGTACCGAAATTACAAGCCAACCAAGACATTCTCGCCAATTCGTATTTGCAGCAGGCGCAAGAACATATTAGTCAATTTGCGTTCCCGGGTATTTAAGCTCAAGCAATTCTGATGTATGTAGAGACCCGTTTAGGCGTGCGTCCACAGATGGTTCGTGAGGAACCGACCGCGGCTCGCTATGCAATCAACCTACAGGTAGTGACACGCCGCAGCGGTAAACAGTGTTGTTTCGCCTTATGTCCGACCGGGTATTGTGTAATCCTCATTGAAACCCAGAGGAGCGGTATCTTAACCCATGCACATCCCTCCCCGACCACTTTCGGGCCGACCTACGCCTCTCTTTTCGTTGACTGCAGAGGCCCACCGAGCGGTCCTGATTTCATTACTTGGATAACCCGAACCTCCGATAGACCTTAGTGCCCGGTCGAGTCAGCCGTCCGCATTTCTCGATATGTCGGAGCACCAGGTGAGTATC

```

## Generate random fastq files with single-end reads

generate_random_fastq_SE generates random sets of single-end reads
usage: generate_random_fastq_SE sequence_size nb_seq
	     
arguments:
sequence_size = length of the sequence in bp (integer)       
nb_seq = number of sequences required (integer)

```
python fastq_generator.py generate_random_fastq_SE 100 10 > random.fastq

head random.fastq

@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:1 1:N:0:TGACCAAT
AGGCCTATATAAGCTCTGAAGGACAGTATTGAGGCGAGCTTACCCCAAATCGGGTGATCCATGACGGAAACCCTTCGCGCTCGTTTGGCTTTAGGTCTTC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:2 1:N:0:TGACCAAT
TCATAACCCCTAGTCCCGTCTGCTGAATTAATAAGGCAAATTACTCGCTTGAGGACCGAACAAGGGATGCCGACTTTTGGGCCCAGAGCTTGTTTAACTT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

```

## Generate random fastq files with paired-end reads

generate_random_fastq_PE generates random sets of paired-end reads
usage: generate_random_fastq_SE sequence_size nb_seq
arguments:
sequence_size = length of the sequence in bp (integer)       
nb_seq = number of sequences required (integer)

*Note that the left reads (/1) and right reads (/2) are within the same file.

```
python fastq_generator.py generate_random_fastq_PE 100  10 > random.fastq

head random.fastq

@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:1#0/1
TAACACTTAGGGGGCTGCTACGCAGGGGCGCAAAGGCGGTGCCACTGATTTTGCTCGAGACATAATCACGCCATATCCACGGCAGCTTTACATGGATTAG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:2#0/1
GGCAGGGCGAGACATGGTAACGGAGACAGGTCTGCTGACGGTCTTTAACGAATATCACTCTTCACCGTACGGGTTTAGTCTATGAGTCCAAGACGCGTGA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:3#0/1
GGGTATCGTTCGCCACGGAACAGCTATGGGTTACAGCACCTCGGTATTCGCTTCTGCCAGTGGCGCTAATGTAGGAGAAAACTGATGCCGTGGATCCCAA

```

## Generate fastq files SE mapping on a specified fasta reference

generate_mapped_fastq_SE generates sets of single-end reads from a reference fasta file
usage: generate_mapped_fastq_SE ref_fasta sequence_size coverage
arguments:
ref_fasta = reference file in fasta format (fasta file)
sequence_size = length of the sequence in bp (integer)
coverage = number of reads covering each bp (integer)

```
python fastq_generator.py generate_mapped_fastq_SE chr1.fa 100 2 > chr1_SE.fastq

head chr1_SE.fastq

@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:1 1:N:0:ATGT
CGTTGCTCTGAATAGGCCGGGAGGACGACTAGGCCGATTCCACGAACTGCACAGCGCCTATAATTTAGAAGGCGGGCACTGGAGATCTAAGATGGTTTTG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:2 1:N:0:ATGT
TCTACTAAAAACCTTTCATCTTTCATGGTCTTTGCCCTAACATTATGGGCCGCGGCCGTTGCAGGAGTAGACCTGATTACTCCTCACTGTATAGTTATCA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:3 1:N:0:ATGT
GCCTTCTCTCAAACTGGCTTCGAATTACTTTAAGAATTTCTCTTATAAGCTACGCTGACCTACACCGCTTTCCAGGCGTGGTACGACAAAGCGGGGACCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:4 1:N:0:ATGT
CCAGAGCGCGAAAGAAATGCACTACTAGGAAGTTGGGCTCTACATTGTTCACAATCGAGAGGTTCTACCTGGCTCTGCACTGTGACTCAACAGCGGCGCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII


```

## Generate fastq files PE mapping on a specified fasta reference

generate_mapped_fastq_PE generates sets of paired-end reads from a reference fasta file
usage: generate_mapped_fastq_PE ref_fasta sequence_size  insertion_size coverage

arguments:
ref_fasta = reference file in fasta format (fasta file)
sequence_size = length of the sequence in bp (integer)       
insertion_size = distance between the paired reads in bp (integer)
coverage = number of reads covering each bp (integer)


```
python fastq_generator.py generate_mapped_fastq_PE chr1.fa 100  100 2 > chr1_PE.fastq

head chr1_PE.fastq

@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:1#0/1
GCGTTGCTCTGAATAGGCCGGGAGGACGACTAGGCCGATTCCACGAACTGCACAGCGCCTATAATTTAGAAGGCGGGCACTGGAGATCTAAGATGGTTTT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:2#0/1
GCCAGAGCGCGAAAGAAATGCACTACTAGGAAGTTGGGCTCTACATTGTTCACAATCGAGAGGTTCTACCTGGCTCTGCACTGTGACTCAACAGCGGCGC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:3#0/1
TGTGGGACTAAACATAGGTCGTTTGTTGCTCCCCGCGCAACTGATCAGGATTAGGCCGTGCGGCAGAAGCACCATCTTTTGCACGCCGCTTGTTTACCTT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:4#0/1
GGCGTTGCTCTGAATAGGCCGGGAGGACGACTAGGCCGATTCCACGAACTGCACAGCGCCTATAATTTAGAAGGCGGGCACTGGAGATCTAAGATGGTTT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@FAKE-SEQ:1:FAKE-FLOWCELL-ID:1:1:0:5#0/1
CGCCAGAGCGCGAAAGAAATGCACTACTAGGAAGTTGGGCTCTACATTGTTCACAATCGAGAGGTTCTACCTGGCTCTGCACTGTGACTCAACAGCGGCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

```

## Generate reads mapping to the reverse strand

```
cat chr1.fa | python reverse_complement.py > chr1_reverse_complement.fa

python fastq_generator.py generate_mapped_fastq_SE chr1_reverse_complement.fa 100 2 > chr1_SE_rc.fastq

# Merge fastq files with reads mapping to both + and - strands
cat chr1.fastq chr1_SE_rc.fastq > chr1_both_strands.fastq

```


## Authors

* **Johan Zicola** - *Initial work* - [johanzi](https://github.com/johanzi)


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details



