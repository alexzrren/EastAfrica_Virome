einsi --thread 6 1-OrderProtSeq/Papillomaviridae.faa > 2-MSA_mafft/Papillomaviridae_aligned.fa
trimal -in 2-MSA_mafft/Papillomaviridae_aligned.fa -out 3-trimAl/Papillomaviridae_trimal.fa -gt 0.8 -cons 5
iqtree -s 3-trimAl/Papillomaviridae_trimal.fa --prefix 6-PhylogenyTree/Papillomaviridae_tree/Papillomaviridae --redo -T 6 --mem 48G --ufboot 1000 --boot-trees
