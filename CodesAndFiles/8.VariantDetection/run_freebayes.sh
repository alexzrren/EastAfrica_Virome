cat requant_group.tsv | cut -f2,3 | sort -u | while read -r vsname seqid
do
echo "cat $PWD/2.rglist/$vsname.rglist | xargs samtools index -M
freebayes -f /jdfssz1/ST_HEALTH/P20Z10200N0206/renzirui/EA_redo/ViralSequences/EastAfrica_256viralstrains.improved.fna -L $PWD/2.rglist/$vsname.rglist -p 1 -r $seqid -0 > $PWD/3.variantcalling/$vsname.vcf" > qsub/$vsname.fb.job
qsub -clear -wd $PWD/qsub -P P20Z10200N0206 -binding linear:1 -l vf=8g,p=1 -q st.q qsub/$vsname.fb.job
done