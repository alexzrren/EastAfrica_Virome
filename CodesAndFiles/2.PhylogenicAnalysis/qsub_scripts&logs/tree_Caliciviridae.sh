einsi --thread 6 1-OrderProtSeq/Caliciviridae.faa > 2-MSA_mafft/Caliciviridae_aligned.fa
trimal -in 2-MSA_mafft/Caliciviridae_aligned.fa -out 3-trimAl/Caliciviridae_trimal.fa -gt 0.8 -cons 5
iqtree -s 3-trimAl/Caliciviridae_trimal.fa --prefix 6-PhylogenyTree/Caliciviridae_tree/Caliciviridae --redo -T 6 --mem 48G --ufboot 1000 --boot-trees
