#!/usr/bin/env python
"""

Usage:
    extract_splice_sites.py referece.fasta intron_annotation intron_list splice_site strand

Reference is the reference genome in fasta format (indexed with samtools faidx)
E.g. genome/Drosophila_melanogaster.BDGP5.66.dna.toplevel.fa

Intron annotation contains a tabbed list of introns, with positions and strand
e.g. results/Drosophila_melanogaster.BDGP5.66.gtf.exons.tab.introns.tab.no_exon_overlap.tab.remove_dups.tab.Fly_Mut_B_s_5_tophat_maternal

The output (sent to STDOUT) consists of fasta with sequence of splice sites according to the given annotation, site and strand

"""
import pysam
import string
import sys

def scan_features(ref,intron_annot,site,chosen_strand):
                  

    fastaFile = pysam.Fastafile(ref)
    intron_file = open(intron_annot,"r")
    for nline in intron_file:
        nline = nline.strip()
        flist = nline.split("\t")
        chr = flist[0]  
        start = flist[1]
        end = flist[2]  
        name = flist[3]   
        strand = flist[4]

        #separate two fasta files...
        if(strand == '+' and chosen_strand == '+'):
	    	if(site == '5p'):
	        	print ">"+name+"\n"+fastaFile.fetch(chr, -10+int(start), 5+int(start))
		if(site == '3p'):
	        	print ">"+name+"\n"+fastaFile.fetch(chr, -8+int(end), 7+int(end))		
        if(strand == '-' and chosen_strand == '-'):    
		# Need to reverse complement...
		complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
	    	if(site == '3p'):
	        	print ">"+name+"\n"+fastaFile.fetch(chr, -8+int(start), 7+int(start)).translate(complements)[::-1]
		if(site == '5p'):
	        	print ">"+name+"\n"+fastaFile.fetch(chr, -7+int(end), 8+int(end)).translate(complements)[::-1]		

    intron_file.close()
    fastaFile.close()
    
if __name__ == "__main__":
  if len(sys.argv) == 5:
    scan_features(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
  else:
    print "extract_splice_sites.py reference.fasta intron_annotation [3p|5p] [+|-]"
