einsi --thread 6 1-OrderProtSeq/Parvoviridae.faa > 2-MSA_mafft/Parvoviridae_aligned.fa
trimal -in 2-MSA_mafft/Parvoviridae_aligned.fa -out 3-trimAl/Parvoviridae_trimal.fa -gt 0.8 -cons 5
iqtree -s 3-trimAl/Parvoviridae_trimal.fa --prefix 6-PhylogenyTree/Parvoviridae_tree/Parvoviridae --redo -T 6 --mem 48G --ufboot 1000 --boot-trees
