use warnings;
use strict;

#Requires a file with 'Chr\tStart\tEnd\tExon_id\tStrand\tTranscript_id\tGene_id
#Exons of each transcript must be on correct order in terms of position
#replaces : with - in IDs because some software (R?) cannot handle :

open(FILE,$ARGV[0]);
my $cur_transcript;
my $cur_gene;
my $cur_start;
my $cur_end;
my $cur_id;
my $cur_strand;
while(<FILE>){
  chomp;
  my ($chr,$start,$end,$exon_id,$strand,$transcript,$gene)=split(/\t/);
  $exon_id=~s/:/-/;
  if($cur_transcript && ($cur_transcript eq $transcript)){
    	#print $chr."\t".($cur_end+3)."\t".($start-3)."\t".$cur_id."_".$id."\n";
	if($strand eq $cur_strand){
		if($strand eq '+'){
			print $chr."\t".($cur_end+1)."\t".($start)."\t".$cur_id."_".$exon_id."\t".$cur_strand."\t".$cur_transcript."\t".$cur_gene."\n";
		} else {
			print $chr."\t".($end+1)."\t".($cur_start)."\t".$cur_id."_".$exon_id."\t".$cur_strand."\t".$cur_transcript."\t".$cur_gene."\n";
	
		}
	} else {
		print STDERR $cur_transcript." ".$exon_id." has a weird reverse exon\n";	
	}
    
  }
  $cur_gene = $gene;
  $cur_transcript = $transcript;
  $cur_start = $start;
  $cur_strand = $strand;
  $cur_end = $end;
  $cur_id=$exon_id;
  
}
close FILE;
