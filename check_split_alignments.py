#!/usr/bin/env python
"""
This is a simple script to calculate an n-bp coverage graph based on input coordinates...
May also depend on flags (e.g. mapping quality, properly paired reads only)
"""
import sys
import pysam
import re


#Counting that tries to take split alignments into account... an alternative take on it
def process_locus_split(samfile,locus_chr, locus_start,locus_end):
  
  #Bam cigar codes:
  # 0 -> match
  # 1 -> insert
  # 2 -> deletion
  # 3 -> skip (split alignment)
  # 4 -> soft clipping (not part of alignment)
  # 5 -> hard clipping (not part of alignment)
  # 6 -> padding
  # 7 -> sequence match (cigarx)
  # 8 -> sequence mismatch (cigarx)

  length = locus_end - locus_start  

  total_locus = 0
  unsplit = 0
  well_split = 0
  mis_split = 0  
 
  for alignedRead in samfile.fetch(locus_chr, locus_start, locus_end):
    
    try:
      #have a minimum flank of the read around either the known splice site, or the split in the read... (?)
      min_flank = 10

      is_potential_unsplit = False
      positions = alignedRead.positions
      cur_pos = 0
      for (event, length) in alignedRead.cigar:
	#print str(event)+"\t"+str(length)
	if(event == 0): 
	  #print "Match\t"+str(positions[cur_pos]+1)+"\t"+str(positions[cur_pos+length-1]+1)
	  #potentially un-spliced, but depends on events that happen after that..
	  if(((positions[cur_pos]+1)< locus_start-min_flank) and ((positions[cur_pos+length-1]+1)> locus_start+min_flank)): #require flanking(?)
	    is_potential_unsplit = True
	  cur_pos += length-1 # To keep it 0-based (check!)
	#split
	if(event == 3):
	   #cannot count if beginning of split is before intron or end of split after intron (may be just different exon being used)
	   #positions from bam always come 0-based, so need to add +1 to them...
	   #print "Split: "+locus_chr+":"+str(positions[cur_pos]+1)+"-"+str(positions[cur_pos+1]+1)+" ("+str(locus_start)+"-"+str(locus_end)+")"
	   
	   is_potential_unsplit = False
	   ##In this case, even if exon usage is different, we can tell there was a mis-split
	   if(((positions[cur_pos]+1)>locus_start) and ((positions[cur_pos]+1)<locus_end)):
	     total_locus += 1
	     mis_split += 1
	     break
	   ##In this case, even if exon usage is different, we can tell there was a mis-split
	   if(((positions[cur_pos+1]+1)>locus_start) and ((positions[cur_pos+1]+1)<locus_end)):
	     total_locus += 1
	     mis_split += 1
	     break
	   ##Normal split: check carefully if positions fit...
	   if(((positions[cur_pos]+1)==locus_start) and ((positions[cur_pos+1]+1)==locus_end)):
	     total_locus += 1
	     well_split += 1
	     break
	   ##Other situations of split we cannot assess as it may be just different exons being used...

      if(is_potential_unsplit):
      	total_locus += 1
      	unsplit += 1
    
    except IndexError:
      sys.stderr.write(alignedRead.cigar)
      sys.stderr.write(alignedRead.positions)
      sys.exit()
      
    
    #print "Totals: %f\t%f\t%f\t%f" % (total_locus, unsplit, well_split, mis_split)
    #sys.exit()
              
  return [total_locus, unsplit, well_split, mis_split]


def check_split_alignments(alignments, introns):

  samfile = pysam.Samfile( alignments, "rb" )
   
  #Arbitrary, pass as parameter...
  #minimum = 50
  #minimum = 10
  
  feat_file=open(introns,"r")

  #sys.stdout.write("Locus\tTotal\tUnsplit\tWell_Split\tMis_Split\n")

  for line in feat_file:
    
    line=line.strip()   
    flist=line.split("\t")
    locus_chr = flist[0]
    locus_start = int(flist[1])
    locus_end = int(flist[2])
    locus_name = flist[3]
    
    #split should be in locus_start-1 and locus_end + 1
    #Carefull with off-by-1 errors with BED (0-based) and others 1-based!!
    #[total_locus, unsplit, well_split, mis_split] = process_locus_split(samfile, locus_chr, locus_start-1, locus_end+1)
    [total_locus, unsplit, well_split, mis_split] = process_locus_split(samfile, locus_chr, locus_start-1, locus_end)
	
    #only print if there is a minimum number of reads in the locus... (?). No... make that afterwards
    #if(total_locus > minimum):     
    #sys.stdout.write(locus_name)
    sys.stdout.write("\t".join(flist))
    #unsplit = 1.0 * unsplit / total_locus
    #well_split = 1.0 * well_split / total_locus
    #mis_split = 1.0 * mis_split / total_locus
    sys.stdout.write("\t%d\t%d\t%d\t%d\n" % (total_locus, unsplit, well_split, mis_split))
    #sys.stdout.write("\t%f\t%f\n" % (well_split, mis_split))
    

instructions = "e.g. python check_split_alignments.py alignments.bam introns.bed"

if __name__=="__main__":
    if len(sys.argv) == 3:
        check_split_alignments(sys.argv[1],sys.argv[2])
    else:
        print instructions
