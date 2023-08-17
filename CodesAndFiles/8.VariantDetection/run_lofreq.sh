for vsname in `cat requant_group.tsv | cut -f2 | sort -u`
do
echo $vsname
mkdir -p 8.lofreq_result/$vsname
cat requant_group.tsv | grep -w $vsname | while read -r sn vsname ref
do
echo $vsname $sn
echo "lofreq call -f /jdfssz1/ST_HEALTH/P20Z10200N0206/renzirui/EA_redo/ViralSequences/EastAfrica_256viralstrains.improved.fna -o $PWD/8.lofreq_result/$vsname/$sn.lofreq.vcf $PWD/1.mappedbam/$vsname/$sn.bam" > qsub/$sn.$vsname.lofreq.job
qsub -clear -wd $PWD/qsub -q st.q -P P20Z10200N0206 -l vf=2g,p=1 -binding linear:1 qsub/$sn.$vsname.lofreq.job
done
done