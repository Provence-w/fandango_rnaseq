use warnings;
use strict;

my $maternal_genes = $ARGV[0];
my $zygotic_genes = $ARGV[1];
my $base_graph = $ARGV[2];
my $column_with_gene = $ARGV[3];

my $conversion_file = $ARGV[4];

if(!$maternal_genes){
	print "Usage: split_general_zygotic_maternal.pl MATERNAL_GENE_LIST ZYGOTIC_GENE_LIST FILE COLUMN_WITH_GENE_ID (0-based) [CONVERSION_FILE]\n";
	exit;
}

my %convert;
if($conversion_file){
	open(FILE,$conversion_file);
	while(<FILE>){
		chomp;
		my ($symb,$gene)=split(/\t/);
		$convert{$symb}=$gene;
	}
	close FILE;
}

my %maternal;
#open(FILE,"expression_data/Gelbart.2010.10.13/maternal_genes_p4000.txt");
open(FILE,$maternal_genes);
while(<FILE>){
  chomp;
  $maternal{$_}=1;
}
close FILE;

my %zygotic;
#open(FILE,"expression_data/Gelbart.2010.10.13/zygotic_genes_l400_p1000.txt");
open(FILE,$zygotic_genes);
while(<FILE>){
  chomp;
  $zygotic{$_}=1;
}
close FILE;

open(FOZ,">".$base_graph."_zygotic");
open(FOM,">".$base_graph."_maternal");
open(FOO,">".$base_graph."_other");
open(FILE,$base_graph);
while(<FILE>){
  chomp;
  my @cols = split(/\s+/);
  my $gene = $cols[$column_with_gene];
  $gene =~ s/"//g; #Just in case... 

  if($conversion_file){
	$gene = $convert{$gene};
  }

  if($zygotic{$gene}){
    print FOZ $_."\n";
  } elsif($maternal{$gene}){
    print FOM $_."\n";
  } else {
    print FOO $_."\n";
  }
}
close FILE;
close FOZ;
close FOM;
close FOO;
