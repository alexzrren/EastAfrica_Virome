einsi --thread 6 1-OrderProtSeq/Herpesviridae.faa > 2-MSA_mafft/Herpesviridae_aligned.fa
trimal -in 2-MSA_mafft/Herpesviridae_aligned.fa -out 3-trimAl/Herpesviridae_trimal.fa -gt 0.8 -cons 5
iqtree -s 3-trimAl/Herpesviridae_trimal.fa --prefix 6-PhylogenyTree/Herpesviridae_tree/Herpesviridae --redo -T 6 --mem 48G --ufboot 1000 --boot-trees
