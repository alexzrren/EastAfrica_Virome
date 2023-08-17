einsi --thread 16 1-OrderProtSeq/Flaviviridae.faa > 2-MSA_mafft/Flaviviridae_aligned.fa
trimal -in 2-MSA_mafft/Flaviviridae_aligned.fa -out 3-trimAl/Flaviviridae_trimal.fa -gt 0.8 -cons 5
iqtree -s 3-trimAl/Flaviviridae_trimal.fa --prefix 6-PhylogenyTree/Flaviviridae_tree/Flaviviridae --redo -T 16 --mem 80G --ufboot 1000 --boot-trees
