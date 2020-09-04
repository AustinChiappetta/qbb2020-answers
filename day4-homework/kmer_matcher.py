#!/usr/bin/env python3

#The matcher should take three arguments:
#kmer_matcher.py <target.fa> <query.fa> <k>

#The script should find k-mer matches and for each write:
#target_sequence_name    target_start    query_start k-mer

import sys
from fasta_iterator_class import FASTAReader

def query_kmer_pos(seq, kmer_n):
    #1. read in sequence
    seq_reader = FASTAReader(open(seq))
    query_kmers = {}

    k = kmer_n
    
    #2. divide query into kmers
    #Make a dictionary (keys are the sequence of kmers,
    #values are lists for the start positions)
    
    for seq_id, sequence in seq_reader:
        for i in range(0, len(sequence) - k + 1):
            #if i > 20:
                #break
            kmer = sequence[i:i + k]
            query_kmers.setdefault(kmer, [])
            query_kmers[kmer] += [i]
        
    return(query_kmers)
        

#print(kmer_pos('/Users/cmdb/qbb2020-answers/day4-homework/droYak2_seq.fa', 11))
        

#3. target sequence kmers
    #make a new dictionary of dictionaries
    #outer keys are kmers
    #inner keys are target seq names
    #values are a list of starting positions for that 
    #kmer in that sq
    #{
    #'kmer1': {
                    #'target_seqid1' : [9, 80, 123],
                    #'target_seqid2' : [80, 314]
                        #}
     #'kmer2': {
                    #'target_seqid1' : [9, 45, 98],
                    #'target_seqid3' : [34, 56]
                        #}
    #}
    
def targets_kmer_pos(seqs, kmer_n):
    #read in target sequences
    seq_reader = FASTAReader(open(seqs))
    targets_kmers = {}

    k = kmer_n
    
    for seq_id, sequence in seq_reader:
        for i in range(0, len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            targets_kmers.setdefault(kmer, {})
            targets_kmers[kmer].setdefault(seq_id, [])
            targets_kmers[kmer][seq_id] += [i]

    return(targets_kmers)

#print(targets_kmer_pos('/Users/cmdb/qbb2020-answers/day4-homework/subset.fa', 11))

#Loop through each dictionary and do the following:
#4. Figure out which keys are shared between the two
    #*loop through keys in query
    #*check if key in target dict
        #if it is in the dict:
            #store/print seq_name, starts in target, 
            #starts in query, kmer_seq
            
def kmer_matcher(query_dict, targets_dict):
    for kmer_seq in query_dict.keys():
        if kmer_seq in targets_dict.keys():
            for targets in targets_dict[kmer_seq]:
                print(targets, targets_dict[kmer_seq][targets],
                     query_dict[kmer_seq], kmer_seq)
                
if __name__ == "__main__":
    
    #q_dict = query_kmer_pos('/Users/cmdb/qbb2020-answers/day4-homework/droYak2_seq.fa', 11)
    q_dict = query_kmer_pos(query, k)

    #t_dict = targets_kmer_pos('/Users/cmdb/qbb2020-answers/day4-homework/subset.fa', 11)
    t_dict = targets_kmer_pos(target, k)

    kmer_matcher(q_dict, t_dict)
        
        
