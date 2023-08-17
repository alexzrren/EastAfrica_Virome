einsi --thread 6 1-OrderProtSeq/Reoviridae.faa > 2-MSA_mafft/Reoviridae_aligned.fa
trimal -in 2-MSA_mafft/Reoviridae_aligned.fa -out 3-trimAl/Reoviridae_trimal.fa -gt 0.8 -cons 5
iqtree -s 3-trimAl/Reoviridae_trimal.fa --prefix 6-PhylogenyTree/Reoviridae_tree/Reoviridae --redo -T 6 --mem 48G --ufboot 1000 --boot-trees
