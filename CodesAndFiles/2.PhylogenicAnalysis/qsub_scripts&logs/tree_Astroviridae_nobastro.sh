einsi --thread 16 1-OrderProtSeq/Astroviridae_nobastro.faa > 2-MSA_mafft/Astroviridae_nobastro_aligned.fa
trimal -in 2-MSA_mafft/Astroviridae_nobastro_aligned.fa -out 3-trimAl/Astroviridae_nobastro_trimal.fa -gt 0.8 -cons 5
iqtree -s 3-trimAl/Astroviridae_nobastro_trimal.fa --prefix 6-PhylogenyTree/Astroviridae_nobastro_tree/Astroviridae_nobastro --redo -T 16 --mem 80G --ufboot 1000 --boot-trees
