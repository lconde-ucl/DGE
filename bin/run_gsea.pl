#!/usr/bin/perl -w


use strict;
use warnings;

no warnings qw/uninitialized/;

use Getopt::Long qw(GetOptionsFromArray :config pass_through no_auto_abbrev bundling);
use File::Basename;



my $gsea="/home/regmond/Scratch/software/gsea-3.0.jar";

my $rand_id=`od -N 4 -t uL -An /dev/urandom | tr -d " " | tr -d "\n"`;

my @args=@ARGV;

my $rnk='';
my $gmx='';	
my $min_set="";
my $max_set="";
my $perm="";
my $output="";
GetOptionsFromArray (
    \@args,
    "rnk=s" => \$rnk,
    "gmx=s"   => \$gmx,
    "out=s"   => \$output,
    "min_set=i"   => \$min_set,
    "max_set=i"   => \$max_set,
    "perm=i"   => \$perm,
    "<>"   => \&print_usage
) or die "\n";

if ($output eq '') {
	$output="gsea_results";
}		
if (!-e "$rnk") {
	usage_gsea("Can't find rnk file <$rnk>");
	return;
}	
if (!-e "$gmx") {
	usage_gsea("Can't find gmx (gene set) file <$gmx>");
	return;
}
		
if($min_set ne ''){
	if(!&isnum($min_set) || ($min_set < 0)){
		usage_gsea("<$min_set> is not a valid number. Please insert an integer > 0");
		return;
	}	
}else{
	$min_set=15;
}
if($max_set ne ''){
	if(!&isnum($max_set) || ($max_set < 0)){
		usage_gsea("<$max_set> is not a valid number. Please insert an integer > 0");
		return;
	}	
}else{
	$max_set=500;
}
if($perm ne ''){
	if(!&isnum($perm) || ($perm < 0)){
		usage_gsea("<$perm> is not a valid number. Please insert an integer > 0");
		return;
	}	
}else{
	$perm=1000;
}

my @suffix=(".txt",".gmt",".gmx",".rnk");
my $basename_rnk = basename($rnk, @suffix);
my $basename_gmx = basename($gmx, @suffix);
	

my $runstr="java -Xmx8g -cp $gsea xtools.gsea.GseaPreranked";
$runstr.=" -gmx $gmx";
$runstr.=" -collapse false -mode Max_probe -norm meandiv -nperm $perm";
$runstr.=" -rnk $rnk";
$runstr.=" -scoring_scheme weighted";
$runstr.=" -rpt_label gsea";
$runstr.=" -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp";
$runstr.=" -set_max $max_set -set_min $min_set -zip_report false";
$runstr.=" -out ./${output}_${basename_rnk}_${basename_gmx}_${rand_id}";
$runstr.=" -gui false";
`$runstr`;	
	
my @genesets=();
my $dir2 = "${output}_${basename_rnk}_${basename_gmx}_${rand_id}";
opendir(DIR, "./${dir2}") || die "MAL\n";
my @dirs= grep { /^gsea.GseaPreranked/ } readdir (DIR);
my $dir=$dirs[0];
opendir(DIR2, "${dir2}/${dir}") || die "MAL\n";;
my @files= grep { /.html$/ && !/^index.html$/ && !/^gsea_report_for_na/ && !/pos_snapshot/ && !/neg_snapshot/} readdir (DIR2);
foreach my $file(@files){
	$file=~s/\.html//;
	push(@genesets, $file);
}
close(DIR2);
close(DIR);
	
my $runstr2="java -Xmx8g -cp $gsea xtools.gsea.LeadingEdgeTool";
$runstr2.=" -dir ./${output}_${basename_rnk}_${basename_gmx}_${rand_id}/$dir";
$runstr2.=" -out ./${output}_${basename_rnk}_${basename_gmx}_${rand_id}";
$runstr2.=" -rpt_label gsea";
$runstr2.=" -extraPlots TRUE -gui FALSE";
$runstr2.=" -gsets ".join(",", @genesets);
`$runstr2`;


open(OUTF, ">gsea_table.txt");
print OUTF "RANK\tGENESET\tGENESET_ORIGINAL_SIZE\tGENESET_SIZE_IN_DATA\tLEADING_EDGE_GENES\tRATIO\tRATIO_ORIGINAL_SIZE\tFDR_p\tNES\n";

#- orioginal geneset sizes
my %originalsize=();
open(INF, "./$dir2/$dir/gene_set_sizes.xls");
while(<INF>){
	chomp $_;
	my @a=split("\t",$_);
	($a[0] eq 'NAME') && next;
	$originalsize{$a[0]}=$a[1];
}
close(INF);
opendir(D, "./$dir2/$dir");
my @files2= grep { /.xls$/ && /^gsea_report_for_na/} readdir (D);
foreach my $file(@files2){
	open(INF, "./$dir2/$dir/$file") || die "$! ./$dir2/$dir/$file\n";
	while(<INF>){
		chomp $_;
		my @a=split("\t",$_);
		($a[0] eq 'NAME') && next;
		
		(!-e "./$dir2/$dir/$a[0].xls") && next;
		
		my $count=0;
		open(INF2, "./$dir2/$dir/$a[0].xls");
		while(<INF2>){
			chomp $_;
			my @a=split("\t",$_);
			($a[0] eq 'NAME') && next;
			
			if($a[7] eq 'Yes'){
				$count++;
			}
		}
		close(INF2);
		my $ratio=sprintf("%.2f", $count/$a[3]);
		my $ratiooriginal=sprintf("%.2f", $count/$originalsize{$a[0]});
		print OUTF $rnk."\t".$a[0]."\t".$originalsize{$a[0]}."\t".$a[3]."\t".$count."\t".$ratio."\t".$ratiooriginal."\t".$a[7]."\t".$a[5]."\n";
	}
	close(INF);
}
closedir(D);
closedir(DIR);

#`PERM=$perm FILE=gsea_table_${basename_rnk}_${basename_gmx}.txt R CMD BATCH plot.R`;


sub isnum ($) {
    return 0 if $_[0] eq '';
    $_[0] ^ $_[0] ? 0 : 1
}


sub print_usage {

	my @a=@_;

	my $usage0="\t";
	my $usage1="\tProgram: Run_GSEA3.0";
	my $usage2 = "\tDescription: Runs gene set enrichment analysis (GSEA3.0) given a ranked list of genes (-rnk) and a gene set (-gmx).";
	my $usage3 = "\tContact: Lucia Conde <l.conde\@ucl.ac.uk>";
	my $usage4="\t";
	my $usage5="\tUsage:   run_gsea --rnk FILE.rnk --gmx FILE.gmt [options]";
	my $usage6="\t";
	my $usage7="\tRequired: --rnk FILE	Ranked list of genes. Needs a column with gene names and a column with the stat";
	my $usage8="\t          --gmx FILE		Gene sets in GMT format (transposed GMX) (https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29)";
	my $usage9="\t";
	my $usage10="\tOptions: --min_set NUM     Ignore gene sets that contain less than NUM genes [15]";
	my $usage11="\t         --max_set NUM	Ignore gene sets that contain more than NUM genes [500]";
	my $usage12="\t         --perm NUM	Number of permutations [1000]";
	my $usage13="\t         --out TEXT	Outputdir prefix [gsea_results]";
        my $usage14="\t";

	print "$usage0\n$usage1\n$usage2\n$usage3";
	print "$usage4\n$usage5\n$usage6\n$usage7\n";
	print "$usage8\n$usage9\n$usage10\n";
	print "$usage11\n$usage12\n$usage13\t$usage14\n";

	die "ERROR: Invalid argument '@a'. Please check above the usage of this script\n";

}

sub usage_gsea {
	my $error=shift;
	die qq(
	Program: Run_GSEA3.0
	Description: Runs gene set enrichment analysis (GSEA3.0) given a ranked list of genes (-rnk) and a gene set (-gmx)
	Contact: Lucia Conde <l.conde\@ucl.ac.uk>

	Usage:   run_gsea --rnk FILE.rnk --gmx FILE.gmt [options]

	Required: --rnk FILE	Ranked list of genes. Needs a column with gene names and a column with the stat
	          --gmx FILE		Gene sets in GMT format (transposed GMX) (https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29)
	Options: --min_set NUM	Ignore gene sets that contain less than NUM genes [15]
		 --perm NUM	Number of permutations [1000]
	         --max_set NUM	Ignore gene sets that contain more than NUM genes [500]
  		 --out TEXT	Outputdir prefix [gsea_results]

	ERROR: $error

	);
}
