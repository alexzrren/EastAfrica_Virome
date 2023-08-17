einsi --thread 6 1-OrderProtSeq/Astroviridae.faa > 2-MSA_mafft/Astroviridae_aligned.fa
trimal -in 2-MSA_mafft/Astroviridae_aligned.fa -out 3-trimAl/Astroviridae_trimal.fa -gt 0.8 -cons 5
iqtree -s 3-trimAl/Astroviridae_trimal.fa --prefix 6-PhylogenyTree/Astroviridae_tree/Astroviridae --redo -T 6 --mem 48G --ufboot 1000 --boot-trees
