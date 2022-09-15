#Provide complementary strand for my_dna given sequence
# In standard input
# For instance
# cat sequence.fa | python reverse_complement.py -

import string
import os
import sys

def reverse_complement(my_dna):
#my_dna = "ACCTTGG"
    my_dna = str(my_dna)
    my_dna = my_dna.upper()
    new_str1 = str.replace(my_dna, 'A','t')
    new_str2 = str.replace(new_str1, 'T','a')
    new_str3 = str.replace(new_str2, 'C','g')
    new_str4 = str.replace(new_str3, 'G','c')
    new_str5 = new_str4.upper()
    new_str6 = new_str5[::-1]

    return new_str6


if not sys.stdin.isatty():
    for line in sys.stdin.readlines():
        if line[0] == ">":
            print(line.strip())
        else:
            rc = reverse_complement(line)
            print(rc.strip())

