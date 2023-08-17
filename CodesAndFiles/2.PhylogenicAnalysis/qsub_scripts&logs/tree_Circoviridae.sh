einsi --thread 6 1-OrderProtSeq/Circoviridae.faa > 2-MSA_mafft/Circoviridae_aligned.fa
trimal -in 2-MSA_mafft/Circoviridae_aligned.fa -out 3-trimAl/Circoviridae_trimal.fa -gt 0.8 -cons 5
iqtree -s 3-trimAl/Circoviridae_trimal.fa --prefix 6-PhylogenyTree/Circoviridae_tree/Circoviridae --redo -T 6 --mem 48G --ufboot 1000 --boot-trees
