#!/usr/bin/env python

#By Mike Schmitt and Brendan Kohrn
#Modified by Dana Nachmanson on November 29, 2016:  
#This version includes the lengths of the insertions and deletions
#Version 1.33

"""
This script gives position-specific mutation frequencies from a tagcounts file given as stdin.

The output is tab-delimited and specifies:
chromosome number, template base, nucleotide position, depth, mutations to T, C, G, A, insertions, deletions, N's


"""

from argparse import ArgumentParser
import sys
import re
import csv
import string

def MutPos(o, f, fOut):
    lines = f.readlines()

    chrom=[]
    pos=[]
    muts = []
    depths = []
    template =[]
    Tcount=[]
    Ccount=[]
    Gcount=[]
    Acount=[]
    insert = []
    inscount=[]
    delet = []
    delcount=[]
    Ncount=[]
    for i,line in enumerate( lines ):

          linebins = line.split()

    #convert sequence information to uppercase
          linebins[4] = linebins[4].replace('t','T')
          linebins[4] = linebins[4].replace('c','C')
          linebins[4] = linebins[4].replace('g','G')
          linebins[4] = linebins[4].replace('a','A')
          linebins[4] = linebins[4].replace('n','N')

    #remove start line, end line, and N entries, as well as 1st and last nucleotide of a read.
          linebins[4] = re.sub('\$','',linebins[4])
          linebins[4] = re.sub('\^.','',linebins[4])    
          #linebins[4] = linebins[4].replace('N','')      

    #create a dictionary of and remove insertions
          ins = {}
          newIns = map(int, re.findall(r'\+\d+', linebins[4]))
          insSeq = map(str, re.findall(r'\+\w+', linebins[4]))
          for j in xrange(len(newIns)):
              length = newIns[j]
              #insertion is indicated as ref base + insertion
              insertion = str(linebins[2]) + insSeq[j][-(length):]
              if insertion in ins.keys():
                  ins[insertion] += 1
              else:
                  ins[insertion] = 1 
              rmStr = r'\+' + str(length) + "."*length
              linebins[4] = re.sub(rmStr, '', linebins[4])
          #check if dictionary is empty
          if not ins:
              ins['-'] = 0
              

    #create a dictionary of and remove deletions
          dels = {}
          newDels = map(str, re.findall(r'-\d+', linebins[4]))
          delSeq = map(str, re.findall(r'-\w+', linebins[4]))
          for k in xrange(len(newDels)):
              length = int(newDels[k][1:])
              #deletion reference is indicated as ref base + insertion
              deletion = str(linebins[2]) + delSeq[k][-(length):]
              if deletion in dels.keys():
                  dels[deletion] += 1
              else:
                  dels[deletion] = 1
              rmStr = r'-' + str(length) + "."*length
              linebins[4] = re.sub(rmStr, '', linebins[4])
          linebins[4] = linebins[4].replace('*','')
          #check if dictionary is empty
          if not dels:
              dels['-'] = 0

    #count depth                                                                                                
          depth = int(linebins[3]) - linebins[4].count('N')
                                                        
    #skip lines that do not meet filtering criteria
          if    (
                depth < o.mindepth
                or
                ((float(max(linebins[4].count('T'),linebins[4].count('C'),linebins[4].count('G'),linebins[4].count('A'), (max(newIns.count(n) for n in list(set(newIns))) if newIns != [] else 0), (max(newDels.count(m) for m in list(set(newDels))) if newDels != [] else 0))) / float(depth)) > o.clonal_max)
                or
                ((float(max(linebins[4].count('T'),linebins[4].count('C'),linebins[4].count('G'),linebins[4].count('A'), (max(newIns.count(n) for n in list(set(newIns))) if newIns != [] else 0), (max(newDels.count(m) for m in list(set(newDels))) if newDels != [] else 0))) / float(depth)) < o.clonal_min)
                or    
                (max(float(linebins[4].count('T')),float(linebins[4].count('C')),float(linebins[4].count('G')),float(linebins[4].count('A')),(max(newIns.count(n) for n in list(set(newIns))) if newIns != [] else 0), (max(newDels.count(m) for m in list(set(newDels))) if newDels != [] else 0)) < o.num_muts)
                ):
                pass
                                                        
          else:

    #count position-specific mutation frequency
                                                        
                mut = linebins[4].count('T') + linebins[4].count('C') + linebins[4].count('G') + linebins[4].count('A') + len(newIns) + len(newDels)
                Tcount.append(linebins[4].count('T'))
                Ccount.append(linebins[4].count('C'))
                Gcount.append(linebins[4].count('G'))
                Acount.append(linebins[4].count('A'))
                Ncount.append(linebins[4].count('N'))
                inscount.append('|'.join(map(str,ins.values())))
                insert.append('|'.join(map(str,ins.keys())))
                delcount.append('|'.join(map(str,dels.values())))
                delet.append('|'.join(map(str,dels.keys())))
                chrom.append(linebins[0])
                pos.append(linebins[1])
                template.append(linebins[2])
                depths.append(depth)
                muts.append(mut)
                                                      
    script_output=zip(chrom, template, pos, depths, muts, Tcount, Ccount, Gcount, Acount, Ncount, inscount, insert, delcount, delet)
    
    csv_writer = csv.writer(fOut, delimiter='\t')
#    csv_writer.writerow('Chr' + '\t' + 'Ref Base' + '\t' + 'Coordinate' + '\t' + 
#                        'DCS Reads' + '\t' + '#Muts'  + '\t' + 'T'  + '\t' + 
#                        'C'  + '\t' + 'G' + '\t' + 'A'  + '\t' + 'Ns' + '\t' + 
#                        'InsCount' + '\t' + 'Ins' + '\t' + 'DelsCount' + '\t' + 
#                        'Dels' + '\n')
    csv_writer.writerow(['CHR','REF','POS', 'DCS','MUTS',
                         'T','C','G','A','Ns','INS', 'INS_ALT',
                         'DEL' ,'DEL_REF'])
    csv_writer.writerows(script_output)


def main():
    parser = ArgumentParser()
    parser.add_argument('-i', '--infile', action ='store', dest = 'inFile', help = 'An imput file. If None, defaults to stdin. [%(default)s]', default = None)
    parser.add_argument('-o', '--outfile', action = 'store', dest = 'outFile', help = 'A filename for the output file.  If None, outputs to stdout.  [%(default)s]', default = None)
    parser.add_argument("-d", "--depth", action="store", type=int, dest="mindepth", 
                      help="Minimum depth for counting mutations at a site [%(default)s]", default=20)
    parser.add_argument("-c", "--min_clonality", action="store", type=float, dest="clonal_min",
                      help="Cutoff of mutant reads for scoring a clonal mutation [%(default)s]", default=0)
    parser.add_argument("-C", "--max_clonality", action="store", type=float, dest="clonal_max",
                      help="Cutoff of mutant reads for scoring a clonal mutation [%(default)s]", default=1)
    parser.add_argument("-n", "--num_muts", action="store", type=int, dest="num_muts",
                      help="Minimum number of mutations for scoring a site [%(default)s]", default=0)
    o = parser.parse_args()
    if o.inFile != None:
        f = open(o.inFile, 'r')
    else:
        f = sys.stdin
    if o.outFile != None:
        fOut = open(o.outFile, 'wb')
    else:
        fOut = sys.stdout
    
    MutPos(o, f, fOut)
        


if __name__ == "__main__":
    main()
