WDIR=/jdfssz1/ST_HEALTH/P20Z10200N0206/renzirui/Wuhan_EastAfrica/HostFstCalculation
NCPU=7
while read -r sampleid sname fq1 fq2 species  
do
echo "bowtie2 -x $reference_index -1 $fq1 -2 $fq2 -p 6 | samtools view -bS -F 4 | samtools sort -@ 1 -o $WDIR/1-mappedBAM/${sampleid}.bam " > $WDIR/qsub/bowtie_${sampleid}.sh
echo "samtools addreplacerg -r "ID:$sampleid" -r "SM:$sampleid" -o $WDIR/2-addRG/$species/${sampleid}.bam $WDIR/1-mappedBAM/${sampleid}.bam" >> $WDIR/qsub/bowtie_${sampleid}.sh

qsub -clear -wd /jdfssz1/ST_HEALTH/P20Z10200N0206/renzirui/Wuhan_EastAfrica/HostFstCalculation/qsub -q st.q -P P20Z10200N0206 -binding linear:7 -l vf=8G,num_proc=7 $WDIR/qsub/bowtie_${sampleid}.sh
done < $WDIR/sampleinfo_usedfeces7bats.tsv

for species in `cat sid2species.tsv | awk '{print $2}' | sort | uniq`
do
echo $species
ls $PWD/2-addRG/$species/* | grep -f <(grep $species sid2species.tsv | awk '{print $1}') > GroupList/${species}.rglist

echo "freebayes -f $WDIR/REF-BatMT/${species}.fna -L $WDIR/GroupList/${species}.rglist -p 1 > $WDIR/3-VCF_Ploidy1/${species}.vcf" > qsub/freebayes_p1_${species}.sh
qsub -clear -wd $WDIR/qsub -q st.q -P P20Z10200N0206 -binding linear:1 -l vf=64G,num_proc=1 qsub/freebayes_p1_${species}.sh
echo "freebayes -f $WDIR/REF-BatMT/${species}.fna -L $WDIR/GroupList/${species}.rglist -p 2 > $WDIR/3-VCF_Ploidy2/${species}.vcf" > qsub/freebayes_p2_${species}.sh
qsub -clear -wd $WDIR/qsub -q st.q -P P20Z10200N0206 -binding linear:1 -l vf=64G,num_proc=1 qsub/freebayes_p2_${species}.sh
done

for vspec in `ls bowtie_mapping`
do
    echo $vspec
#REMOVE MISSING LOC
    vcftools --vcf variant_calling/orig/${vspec}.vcf --max-missing 0.75 --max-alleles 2 --minQ 30 --minDP 1 --max-alleles 2 --remove-filtered-all --recode --recode-INFO-all --stdout | vcfsnps | bgzip -c > variant_calling/qc1/${vspec}_snp.vcf.gz

#REMOVE MISSING INDV
    vcftools --gzvcf variant_calling/qc1/${vspec}_snp.vcf.gz --missing-indv
    awk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv
    vcftools --gzvcf variant_calling/qc1/${vspec}_snp.vcf.gz --remove lowDP.indv --remove-filtered-all --recode --recode-INFO-all --stdout | vcfsnps  > variant_calling/qc2/${vspec}_snp.vcf

#VCF TO FASTA
    cat variant_calling/qc2/${vspec}_snp.vcf | /ldfssz1/ST_INFECTION/P17Z10200N0536_Echinococcus/USER/wangdaxi/Github/HPGAP_HPC/Tools/vcf-to-tab > popgenome.tab
    ./vcf_tab_to_fasta_alignment.pl -i popgenome.tab > snp_alignment/${vspec}.fasta
    echo $vspec done  
done
#CLEANING
rm popgenome.tab* -f
rm lowDP.indv -f
rm out.imiss -f
