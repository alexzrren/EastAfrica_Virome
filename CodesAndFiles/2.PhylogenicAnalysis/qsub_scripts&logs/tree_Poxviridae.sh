einsi --thread 6 1-OrderProtSeq/Poxviridae.faa > 2-MSA_mafft/Poxviridae_aligned.fa
trimal -in 2-MSA_mafft/Poxviridae_aligned.fa -out 3-trimAl/Poxviridae_trimal.fa -gt 0.8 -cons 5
iqtree -s 3-trimAl/Poxviridae_trimal.fa --prefix 6-PhylogenyTree/Poxviridae_tree/Poxviridae --redo -T 6 --mem 48G --ufboot 1000 --boot-trees
