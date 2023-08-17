einsi --thread 6 1-OrderProtSeq/Coronaviridae.faa > 2-MSA_mafft/Coronaviridae_aligned.fa
trimal -in 2-MSA_mafft/Coronaviridae_aligned.fa -out 3-trimAl/Coronaviridae_trimal.fa -gt 0.8 -cons 5
iqtree -s 3-trimAl/Coronaviridae_trimal.fa --prefix 6-PhylogenyTree/Coronaviridae_tree/Coronaviridae --redo -T 6 --mem 48G --ufboot 1000 --boot-trees
