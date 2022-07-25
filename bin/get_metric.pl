#!/usr/bin/perl -w


use strict;


my $file=$ARGV[0];
my $contrast=$ARGV[1];

my %stats=();
open(INF, "$file");
while(<INF>){
	chomp $_;
	my @a=split("\t",$_);
	($a[0] eq '') && next;

	my $fc=$a[2];
	my $gene=$a[0];
	my $pval=$a[4];
	if($pval == 0 ){
		print $_."\n";
		$pval=10e-308;
	}

	my $stat;
	if ($fc>0){
		$stat = -log($pval);
	}else{
		$stat=(-1) * -log($pval);
	}	
	$stats{$gene}=$stat;
		
}
close(INF);

open(OUTF, ">${contrast}.rnk");
print OUTF "Gene_name\tStat\n";
foreach my $gene(sort_numerically_by_values(\%stats)){
	print OUTF $gene."\t".$stats{$gene}."\n";
}
close(OUTF);


sub sort_numerically_by_values {
	my $hs = shift;
	return sort {$hs->{$b} <=> $hs->{$a}} keys %$hs;
}


