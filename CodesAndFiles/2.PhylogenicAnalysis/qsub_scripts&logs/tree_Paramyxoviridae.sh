einsi --thread 6 1-OrderProtSeq/Paramyxoviridae.faa > 2-MSA_mafft/Paramyxoviridae_aligned.fa
trimal -in 2-MSA_mafft/Paramyxoviridae_aligned.fa -out 3-trimAl/Paramyxoviridae_trimal.fa -gt 0.8 -cons 5
iqtree -s 3-trimAl/Paramyxoviridae_trimal.fa --prefix 6-PhylogenyTree/Paramyxoviridae_tree/Paramyxoviridae --redo -T 6 --mem 48G --ufboot 1000 --boot-trees
