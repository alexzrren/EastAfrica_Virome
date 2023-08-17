index=../../ViralSequences/bowtie2_index/EastAfrica_256viralstrains.improved

sortedbam=$PWD/1.sortedbam
coverage=$PWD/2.coverage

mkdir -p qsub $sortedbam $coverage

cat ../../Metadata/EA_1282sample_paths.tsv | grep R2007002118 | while read -r sn sample fq1 fq2 assembly
do
echo "bowtie2 -x $index -1 $fq1 -2 $fq2 -p 8 --very-sensitive-local | samtools view -bS -F 4 | samtools sort -@ 3 -o $sortedbam/${sn}_sorted.bam
samtools index $sortedbam/${sn}_sorted.bam
samtools coverage -q 0 --ff UNMAP,SECONDARY,QCFAIL,DUP -H $sortedbam/${sn}_sorted.bam | sed 's/^/'"$sn"'\t/g' > $coverage/${sn}_coverge.txt" > qsub/${sn}_cvg.orig.sh
qsub -clear -wd $PWD/qsub -binding linear:4 -P P20Z10200N0206 -q st.q -l vf=1g,num_proc=8 qsub/${sn}_cvg.orig.sh

echo "seqkit stats $fq1 $fq2 > $PWD/0.fqstat/$sn.stats.txt" > $PWD/qsub/$sn.fqstat.sh
qsub -clear -wd $PWD/qsub -binding linear:1 -P P20Z10200N0206 -q st_short.q -l vf=1g,num_proc=1 $PWD/qsub/$sn.fqstat.sh
done