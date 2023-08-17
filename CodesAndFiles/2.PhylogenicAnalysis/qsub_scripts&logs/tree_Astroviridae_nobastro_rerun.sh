einsi --thread 16 1-OrderProtSeq/Astroviridae_nobastro_rerun.faa > 2-MSA_mafft/Astroviridae_nobastro_rerun_aligned.fa
trimal -in 2-MSA_mafft/Astroviridae_nobastro_rerun_aligned.fa -out 3-trimAl/Astroviridae_nobastro_rerun_trimal.fa -gt 0.8 -cons 5
iqtree -s 3-trimAl/Astroviridae_nobastro_rerun_trimal.fa --prefix 6-PhylogenyTree/Astroviridae_nobastro_rerun_tree/Astroviridae_nobastro_rerun --redo -T 16 --mem 80G --ufboot 1000 --boot-trees
