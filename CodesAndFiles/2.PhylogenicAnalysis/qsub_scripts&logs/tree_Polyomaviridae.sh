einsi --thread 6 1-OrderProtSeq/Polyomaviridae.faa > 2-MSA_mafft/Polyomaviridae_aligned.fa
trimal -in 2-MSA_mafft/Polyomaviridae_aligned.fa -out 3-trimAl/Polyomaviridae_trimal.fa -gt 0.8 -cons 5
iqtree -s 3-trimAl/Polyomaviridae_trimal.fa --prefix 6-PhylogenyTree/Polyomaviridae_tree/Polyomaviridae --redo -T 6 --mem 48G --ufboot 1000 --boot-trees
