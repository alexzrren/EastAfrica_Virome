einsi --thread 6 1-OrderProtSeq/Adenoviridae.faa > 2-MSA_mafft/Adenoviridae_aligned.fa
trimal -in 2-MSA_mafft/Adenoviridae_aligned.fa -out 3-trimAl/Adenoviridae_trimal.fa -gt 0.8 -cons 5
iqtree -s 3-trimAl/Adenoviridae_trimal.fa --prefix 6-PhylogenyTree/Adenoviridae_tree/Adenoviridae --redo -T 6 --mem 48G --ufboot 1000 --boot-trees
