einsi --thread 6 1-OrderProtSeq/Picornaviridae.faa > 2-MSA_mafft/Picornaviridae_aligned.fa
trimal -in 2-MSA_mafft/Picornaviridae_aligned.fa -out 3-trimAl/Picornaviridae_trimal.fa -gt 0.8 -cons 5
iqtree -s 3-trimAl/Picornaviridae_trimal.fa --prefix 6-PhylogenyTree/Picornaviridae_tree/Picornaviridae --redo -T 6 --mem 48G --ufboot 1000 --boot-trees
