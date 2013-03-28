use warnings;
use strict;

my %info;

open(FILE,$ARGV[0]);
while(<FILE>){
	chomp;
	my @fields = split(/\t/);
	if(!defined($info{$fields[0]."_".$fields[1]."_".$fields[2]})){
		print join("\t",@fields)."\n";
	}
	$info{$fields[0]."_".$fields[1]."_".$fields[2]}  = 1;
}
close FILE;
