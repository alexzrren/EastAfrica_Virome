cat filtered_requantresult.tsv | grep -f ge4samp.vsname.list | while read -r sn vsname;
do
sitenum=`zcat 9.lofreq_iSNV/$vsname/$sn.iSNV.vcf.gz | grep -v '#' | grep -v 'iSNV=REF$' | wc -l`
if [ $sitenum -gt 0 ]; then
mkdir -p 9.1.lofreq_iSNV_ALTonly/$vsname
zcat 9.lofreq_iSNV/$vsname/$sn.iSNV.vcf.gz | grep -v 'iSNV=REF$' | sed 's/;iSNV=ALT//g' | bgzip -c > 9.1.lofreq_iSNV_ALTonly/$vsname/$sn.ALTiSNV.vcf.gz
bcftools index 9.1.lofreq_iSNV_ALTonly/$vsname/$sn.ALTiSNV.vcf.gz

mkdir -p 9.2.lofreq_iSNV_annotated/$vsname

echo "/ldfssz1/ST_HEALTH/P20Z10200N0206/P20Z10200N0206_pathogendb/renzirui/miniconda3/bin/java -jar /jdfssz1/ST_HEALTH/P20Z10200N0206/renzirui/Softwares/snpEff/snpEff.jar eff -v EA_$vsname.prodigal $PWD/9.lofreq_iSNV/$vsname/$sn.iSNV.vcf.gz > $PWD/9.2.lofreq_iSNV_annotated/$vsname/$sn.ALTiSNV.vcf" > snpeff_qsub/$vsname.$sn.snpeff.job
qsub -clear -wd $PWD/snpeff_qsub -q st_short.q -P P20Z10200N0206 -binding linear:1 -l vf=1g,p=1 snpeff_qsub/$vsname.$sn.snpeff.job

fi
done