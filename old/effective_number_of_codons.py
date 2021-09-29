#!/usr/bin/env python3
# effective number of codons function
# based on:
# Estimating the “Effective Number of Codons”: The Wright Way of Determining Codon Homozygosity Leads to Superior Estimates
# Anders Fuglsang
# GENETICS February 1, 2006 vol. 172 no. 2 1301-1307; https://doi.org/10.1534/genetics.105.049643
# https://www.genetics.org/content/172/2/1301

import collections
import sys
import os

def count_codons(seq):
    codons = [ seq.upper()[i:i+3] for i in range(0,len(seq),3) ] # split seq by codons
    codons = list(filter(lambda x: len(x) == 3, codons)) # remove incomplete codons
    codons = list(filter(lambda x: x not in ['TAA','TAG','TGA'], codons)) # remove stops
    codons = list(filter(lambda x: get_aa(x) is not None, codons)) # remove codons with missing data
    return(dict(collections.Counter(codons)))

def get_aa(cod):
    codon2aa = {
        #1
        'ATG' : 'M', # Met
        'TGG' : 'W', # Trp

        'TTT' : 'F', 'TTC' : 'F', # Phe
        'TAT' : 'Y', 'TAC' : 'Y', # Tyr
        'CAT' : 'H', 'CAC' : 'H', # His
        'CAA' : 'Q', 'CAG' : 'Q', # Gln
        'AAT' : 'N', 'AAC' : 'N', # Asn
        'AAA' : 'K', 'AAG' : 'K', # Lys
        'GAT' : 'D', 'GAC' : 'D', # Asp
        'GAA' : 'E', 'GAG' : 'E', # Glu
        'TGT' : 'C', 'TGC' : 'C', # Cys

        #3
        'ATT' : 'I', 'ATC' : 'I', 'ATA' : 'I', # Ile

        #4
        'GTT' : 'V', 'GTC' : 'V', 'GTA' : 'V', 'GTG' : 'V', # Val
        'CCT' : 'P', 'CCC' : 'P', 'CCA' : 'P', 'CCG' : 'P', # Pro
        'ACT' : 'T', 'ACC' : 'T', 'ACA' : 'T', 'ACG' : 'T', # Thr
        'GCT' : 'A', 'GCC' : 'A', 'GCA' : 'A', 'GCG' : 'A', # Ala
        'GGT' : 'G', 'GGC' : 'G', 'GGA' : 'G', 'GGG' : 'G', # Gly

        #6
        'TTA' : 'L', 'TTG' : 'L', 'CTT' : 'L', 'CTC' : 'L', 'CTA' : 'L', 'CTG' : 'L', # Leu
        'TCT' : 'S', 'TCC' : 'S', 'TCA' : 'S', 'TCG' : 'S', 'AGT' : 'S', 'AGC' : 'S', # Ser
        'CGT' : 'R', 'CGC' : 'R', 'CGA' : 'R', 'CGG' : 'R', 'AGA' : 'R', 'AGG' : 'R' # Arg
    }
    if codon2aa.get(cod):
        return(codon2aa[cod])

def get_Faa(codons_in_aa, count):
    if len(codons_in_aa) != 0:
        total = sum([ count[x] for x,y in codons_in_aa ])
        Faa = sum([ (count[x]/total) ** 2 for x,y in codons_in_aa ])
        return(Faa)
    else:
        return(0)

def ENC(count):
    codons = count.keys()
    all_aa = [ get_aa(c) for c in codons ]
    uniq_aa = list(set(all_aa))
    all_Faa = []
    for aa in uniq_aa:
        codons_in_aa = list(filter(lambda x: x[1] == aa, list(zip(codons,all_aa))))
        all_Faa.append(get_Faa(codons_in_aa, count))
    Nc = sum([ 1/x for x in all_Faa ])
    return(Nc)

if __name__ == '__main__':
    # check if something has been piped in
    if not os.isatty(0):
        print("gene,ENC,length")
        head = ''
        for line in sys.stdin.readlines():
            if line[0] == ">":
                if len(head) != 0 and len(seq) != 0:
                        count = count_codons(seq)
                        total_codons = sum(list(count.values()))
                        ENC_est = ENC(count)
                        print(head[1:],ENC_est,total_codons,sep=',')
                head = line.rstrip()
                seq = ''
            else:
                seq += line.rstrip()
    else:
        input_file = sys.argv[1]
        print("gene,ENC,length")
        with open(input_file, 'r') as f:
            head = ''
            for line in f:
                if line[0] == ">":
                    if len(head) != 0 and len(seq) != 0:
                            count = count_codons(seq)
                            total_codons = sum(list(count.values()))
                            ENC_est = ENC(count)
                            print(head[1:],ENC_est,total_codons,sep=',')
                    head = line.rstrip()
                    seq = ''
                else:
                    seq += line.rstrip()

