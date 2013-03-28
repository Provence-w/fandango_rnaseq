use warnings;
use strict;

if(!$ARGV[0] || !$ARGV[1] || !$ARGV[2]){
  print STDERR "Usage make_intron_exon_boundaries exon.tab intron.tab length [3p]\n";
  exit;
}

my $exon_file=$ARGV[0];
my $intron_file=$ARGV[1];
my $len=$ARGV[2];
my $side='5p';
if($ARGV[3]){
	$side='3p';
}

#Get exon length
my %exon_len;
open(FILE,$exon_file);
while(<FILE>){
  chomp;
  my ($ch,$st,$en,$id)=split(/\t/);
  $exon_len{$id}=$en-$st;
}
close FILE;

open(FILE,$intron_file);
while(<FILE>){
	chomp;
	my ($ch,$st,$en,$id,$str,$tr,$gn) = split(/\t/);
	if(!$ch || !$st || !$en || !$id){
	  print STDERR "Intron file must be: 'Chr\tStart\tEnd\tExon5p_Exon3p'\n";
	  exit;
	}
	$id=~s/-/:/g;
	my $intron_len = $en-$st;
	if($intron_len<$len){
	  print STDERR "Ignoring intron smaller than $len : $id \n";
	  next;
	}
	$id =~/([^_]+_\d+)_([^_]+_\d+)/;
	my $exon_5p = $1;
	my $exon_3p = $2;
	
	my $exon;
	if($side eq '5p'){
		$exon = $exon_5p;
	} else {
		$exon = $exon_3p;
	}

	#my ($exon_5p,$exon_3p)=split(/_/,$id);
	if(!defined($exon_len{$exon})){
	  print STDERR "Ignoring intron where could not check $side exon: $id \n";
	  next;
	}
	if($exon_len{$exon} < $len ){
	  print STDERR "Ignoring intron where $side exon was smaller than $len : $id \n";
	  next;
	}
	
	if($side eq '5p'){
		print $ch."\t".($st-$len)."\t".($st+$len)."\t".$id."\t".$str."\t".$tr."\t".$gn."\n";
	} else {
		print $ch."\t".($en-$len)."\t".($en+$len)."\t".$id."\t".$str."\t".$tr."\t".$gn."\n";		
	}
	#print $cols[0]."\t".$cols[1]."\t".$cols[2]."\n";
}
close FILE;
