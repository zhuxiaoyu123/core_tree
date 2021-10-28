# core_tree
Construct ML tree based on multi genes

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
