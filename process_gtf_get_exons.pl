use warnings;
use strict;

#e.g. genome/Drosophila_melanogaster.BDGP5.66.gtf
my $gtf_file = $ARGV[0];

open(FILE,$gtf_file);
while(<FILE>){
	chomp;
	#3R      protein_coding  exon    380     1913    .       +       .        gene_id "FBgn0037213"; transcript_id "FBtr0078962"; exon_number "1"; gene_name "CG12581"; transcript_name "CG12581-RA"; seqedit "false";
	my ($chr,$type,$class,$start,$end,undef,$strand,undef,$attrs)=split(/\t/);
	if($class eq 'exon'){
		if($attrs =~/transcript_id\s+\"([^\"]+)\"/){
			my $transcript = $1;
			if($attrs =~/exon_number\s+\"(\d+)\"/){				
				my $exon_number = $1;
				$attrs =~/gene_id\s+\"([^\"]+)\"/;
				my $gene_id = $1;
				print $chr."\t".$start."\t".$end."\t".$transcript."_".$exon_number."\t".$strand."\t".$transcript."\t".$gene_id."\n";
			} else {
				print STDERR "WARNING NO EXON NUMBER\n";
			}
		} else {
			print STDERR "WARNING NO TRANSCRIPT ID\n";
		}
	}
}
close FILE;
