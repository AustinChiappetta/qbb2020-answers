#!/usr/bin/env python3

import sys

def calcN50(ref, contigs):
     """
     Reports the N50 from a ref sequence.fai and contigs.fai
     Usage: calc50(ref_file.fai, contigs_file.fai)
     Expecting ref_file.fai to be one line long, contigs_file.fai can be multiple lines long"""

     ref_file = open(ref, 'r')
     contigs_file = open(contigs, 'r')

     #1. Calculate ref length/2
     for line in ref_file:
          fields = line.split()
          ref_length = int(fields[1])
          target_length = ref_length/2
          #print('target length =', target_length)
     ref_file.close()

     #2. Grab contigs lengths, sort
     lengths = []
     for line in contigs_file:
          fields = line.split()
          contig_length = int(fields[1])
          lengths += [contig_length]
     #print('this is the lengths list:', lengths)
     lengths.sort(reverse = True) #sort contig lengths largest to smallest
     #print('this is lengths list after sorting', lengths)
     contigs_file.close()

     #3. Find N50:
     size = 0
     for length in lengths:
          size += length
          if target_length <= size:
               N50 = length
               break
     print(N50)
    
if __name__ == "__main__":
     calcN50(sys.argv[1], sys.argv[2])
