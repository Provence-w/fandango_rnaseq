#!/usr/bin/env python
"""
This is a simple script to calculate an n-bp coverage graph based on input coordinates...
May also depend on flags (e.g. mapping quality, properly paired reads only)
"""
import sys
import pysam
import re


#Basic counting: does not work for split alignments
def process_locus_single(samfile,locus_chr, locus_start,locus_end):
  
  length = locus_end - locus_start  
  
  total_locus = samfile.count(locus_chr, locus_start, locus_end)
  
  counts = [0]*length
  for pileupcolumn in samfile.pileup(locus_chr, locus_start, locus_end):

    cur_rel_pos = pileupcolumn.pos - locus_start		
    if((cur_rel_pos>=0) and (cur_rel_pos<length)):
      
      ##Simplest: just count reads that overlap this region... :
      counts[cur_rel_pos]= pileupcolumn.n
  
  return [total_locus, counts]
  
 

#Counting that tries to take split alignments into account... an alternative take on it
def process_locus_split(samfile,locus_chr, locus_start,locus_end):
  
  loci = {}    
  length = locus_end - locus_start  
  
  total_locus = 0
  
  counts = [0]*length
  for alignedRead in samfile.fetch(locus_chr, locus_start, locus_end):
    
    #print alignedRead.cigar   
    #print alignedRead.alen   
    #print alignedRead.positions
    for pos in alignedRead.positions:
      cur_rel_pos = pos - locus_start
      if((cur_rel_pos>=0) and (cur_rel_pos<length)):
	  loci[alignedRead.qname] = 1
	  counts[cur_rel_pos] += 1  
    #sys.exit()
    #if(pileupread.alignment.seq[pileupread.qpos] in nts):
    #loci[pileupread.alignment.qname] = 1
    #counts[cur_rel_pos] += 1  

  total_locus = len(loci)
              
  return [total_locus, counts]


def create_coverage_graph(alignments, input_bed, length, minimum):
  samfile = pysam.Samfile( alignments, "rb" )
   
  feat_file=open(input_bed,"r")

  for line in feat_file:
    
    line=line.strip()   
    flist=line.split("\t")
    locus_chr = flist[0]
    locus_start = int(flist[1])
    locus_end = int(flist[2])
    locus_name = flist[3]
    if((locus_end - locus_start) != length):
      sys.stderr.write('Feature of incorrect size\n')
      continue
        
    #Choose which way to process depending on the input file... (maybe pass as parameter)?    
    #[total_locus, counts] = process_locus_single(samfile, locus_chr, locus_start, locus_end)        
    [total_locus, counts] = process_locus_split(samfile, locus_chr, locus_start, locus_end)
	
    #only print if there is a minimum number of reads in the locus...
    if(total_locus > minimum):     
      sys.stdout.write(locus_name)
      for value in counts:
	#Normalize...
	value = 1.0 * value / total_locus
	sys.stdout.write("\t%f" % value)
      sys.stdout.write("\n")

    

instructions = "e.g. python create_coverage_graph.py alignments.bam features.bed len [min=50]"

if __name__=="__main__":
    if len(sys.argv) == 4:
        create_coverage_graph(sys.argv[1],sys.argv[2], int(sys.argv[3]), 50)
    else:
	if len(sys.argv) == 5:
            create_coverage_graph(sys.argv[1],sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
	else:
            print instructions
