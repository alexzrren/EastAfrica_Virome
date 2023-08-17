cat requant_group.tsv  | awk '{print $2,$3}'  OFS='\t' | sort -u  | while read -r vsname seqid;
do
echo $vsname
resultdir=prodigal_annotation/EA_$vsname.prodigal
mkdir -p prodigal_annotation/EA_$vsname.prodigal
seqkit grep -p $seqid ../ViralSequences/EastAfrica_256viralstrains.improved.fna > $resultdir/sequences.fa 
prodigal -g 1 -f gff -o $resultdir/genes.gff -i $resultdir/sequences.fa -a $resultdir/protein.fa -d $resultdir/cds.fa -p meta
echo -e "\n#\nEA_$vsname.prodigal.genome : Viruses" >> /jdfssz1/ST_HEALTH/P20Z10200N0206/renzirui/Softwares/snpEff/snpEff.config
ln -s /jdfssz1/ST_HEALTH/P20Z10200N0206/renzirui/EA_redo/ViralStrain_PopGenetics/$resultdir  /jdfssz1/ST_HEALTH/P20Z10200N0206/renzirui/Softwares/snpEff/data
cd /jdfssz1/ST_HEALTH/P20Z10200N0206/renzirui/Softwares/snpEff/
echo "pwd cd to : $PWD"
echo "java -jar /jdfssz1/ST_HEALTH/P20Z10200N0206/renzirui/Softwares/snpEff/snpEff.jar build -gff3 -v EA_$vsname.prodigal -noCheckCDS -noCheckProtein" | qsub -clear -wd /jdfssz1/ST_HEALTH/P20Z10200N0206/renzirui/Softwares/snpEff -binding linear:1 -P P20Z10200N0206 -l vf=1g,p=1 -q st.q -N $vsname.snpEff.job
cd -
echo "pwd cd to : $PWD"
done