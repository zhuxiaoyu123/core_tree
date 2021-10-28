#! usr/bin/perl
#written by xiaoyuzhu 2021-10-27
use Getopt::Long;
use Term::ANSIColor qw(:constants); #color
use POSIX qw(strftime); #time
use strict;

my ($proteins, $outdir, $tree, $help);
GetOptions
(
	"p|proteins=s" => \$proteins,          # string
	"o|outdir=s" => \$outdir,
	"t|tree=s" => \$tree,
	"h|help" => \$help                 # flag
);

################################################################################
# usage
############################# usage begin ######################################
my $usage= << "USAGE";

Author: xiaoyuzhu 2021-10-27
Github: https://github.com/zhuxiaoyu123/core_tree

Function: This script is used to construct phylogenetic tree based on multi core genes. 
The core genes dosen't mean they present in every species, instead they can be absent in
some species. You only need to provide fasta files which contain homologous genes from different 
genomes and ensure sequence names in the right format as described in "Requirement" part.

Dependency: MAFTT, TrimAl, fasttree or iqtree and seqkit. Please ensure they have been installed ahead.

Example: 
##concatenate sequences:       
perl $0 -p faa -o core_tree
##concatenate sequences and build tree:
perl $0 -p faa -o core_tree -t fasttree

Options:
	-p -proteins <dir>           required, the path to protein sequences file of each marker gene (files must end with .faa)
	-o -outdir <dir>             required, the output folder
	-t -tree <tree menthod>      optional, two options: "iqtree" or "fasttree"
	-h -help                     print help information

Requirement:
The sequence names must be composed of genome name and gene id with "__" as separator:
>MAG6666__gene001
>MAG7777__gene002
To meet this requirement, we suggest you add the genome names in front of each gene once you get the gene sequence of each genome,
which is a good habit and is benefit for your downstream analyses.
You can simply conduct the following command in your directory of faa files:
"for filename in *.faa; do base=\$(basename \$filename .faa); echo \$base; sed -i "s/^>/>\${base}__/g" \$filename; done"

Note: 
When you apply iqtree to build tree, it may inform you:
"OMP: Info #270: omp_set_nested routine deprecated, please use omp_set_max_active_levels instead."
This message can be ignored.

USAGE
############################## usage end #######################################

if ($help || !(defined $proteins) || !(defined $outdir)){
	print BRIGHT_YELLOW, $usage."\n", RESET;
	exit;
}

my $start_time=time;
print strftime("\nStart time is %Y-%m-%d %H:%M:%S\n", localtime(time));

#generate the genome array
chomp (my $genomes= `cat $proteins/*faa | grep \">\" | sed 's/>//g'| sed 's/__.*//g' | sort | uniq`);
my @genomes = split/\n/, $genomes;

#generate the gene array
chomp (my $genes = `cd $proteins && ls *faa | sed 's/.faa//g'`);
my @genes = split /\n/, $genes;


#identify multi-copy genes and retain only one sequence
system "mkdir $outdir";
print "starting to identify multi-copy genes and retain only one sequence for each ...\n\n";
foreach my $genes (@genes) {
	my @arr;
	open OUTPUT, '>>', "$outdir/$genes\_temp.faa" || die "Cannot creat the file: $outdir/$genes\_temp.faa\n";
	open RP, '<', "$proteins/$genes.faa" || die "Cannot open the file: $proteins/$genes.faa\n";
	$/ = ">";
	while (<RP>) {
		chomp;
		my $gene_list = (split/__/,(split/\n/,$_)[0])[0];  ##基因组名按两个下划线分隔
		if ((grep {$_ eq $gene_list} @arr) || !(defined $gene_list)){
			next;
		}else{
			print OUTPUT ">" . $_;
			push @arr, $gene_list;
		}	
	}
}
#add empty sequences names to single-copy files
system "touch $outdir/temp.faa"; 
foreach my $genes (@genes) {
	my %hash;
	open TEMP, '<', "$outdir/$genes\_temp.faa";
	$/ = ">";
	while (<TEMP>) {
		chomp;	
		my $gene_list = (split/__/,(split/\n/,$_)[0])[0]; 
		$hash{$gene_list} = $_;
	}
	open OUT, '>>', "$outdir/$genes\_single_copy.faa";
	foreach my $genome (@genomes) {
		if (exists $hash{$genome}) {      
			print OUT ">" . $hash{$genome};	
		}else{
			print OUT ">"."$genome\__ribosomal_protein_$genes\n";
		}
	}
}
system "rm -rf $outdir/*temp.faa";

##align and trim
print "starting to align sequences ...\n\n";
foreach my $genes (@genes) {
system "mafft --quiet --maxiterate 1000 --localpair  $outdir/$genes\_single_copy.faa > $outdir/$genes\_mafft.faa";
}
print "starting to trim sequences ...\n\n";
foreach my $genes (@genes) {
system "trimal -in $outdir/$genes\_mafft.faa -out $outdir/$genes\_temp1.faa -automated1 -keepseqs"; #-automated1 is suitable for ML tree
system "seqkit seq -w 0 $outdir/$genes\_temp1.faa >> $outdir/$genes\_trimal.faa && rm -rf $outdir/*_temp1.faa";
}

print "starting to concatenate sequences ...\n\n";
my %newhash;
foreach my $genes (@genes) {
	open FINAL, '<', "$outdir/$genes\_trimal.faa" || die "Cannot open the file: $outdir/$genes\_trimal.faa\n";
	$/ = ">";
	while (<FINAL>) {
		chomp;	
		my $gene_list = (split/__/,(split/\n/,$_)[0])[0]; 
		my $sequence = (split/\n/,$_)[1];
		$newhash{$gene_list}.= $sequence;
	}
}
open CONCAT, '>>', "$outdir/01.Concatenated.faa" || die "Cannot creat the file: 01.Concatenated.faa\n";
foreach (@genomes) {
	print CONCAT ">".$_."\n";
	print CONCAT $newhash{$_}."\n";

}

##Construct trees using Fastree or IQtree
if ($tree eq "iqtree") {
	print "starting to constuct tree with iqtree ...\n";
	system "cd $outdir && iqtree -s 01.Concatenated.faa -alrt 1000 -bb 1000 --quiet --prefix 02.RP_iqtree";
}elsif ($tree eq "fasttree") {
	print "\nstarting to constuct tree with fasttree ...\nFasttree log:\n";
	system "FastTreeMP -gamma $outdir/01.Concatenated.faa > $outdir/02.RP_fasttree.nwk";
}

##Record the program running time!
my $duration_time=time-$start_time;
print strftime("\n\nEnd time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Running time: $duration_time seconds\.\n";
print "Congratulations!!!\n\n";