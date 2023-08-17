# read config info into this namespace
#localrules:all,QC_gzipFQ,QC_seqkit,megahit_seqkit,extract_id,CAT_extract_id
#测试1
#configfile: "config.yaml"

SAMPLES = {}
with open(config['params']['samples'], 'rb') as f:
    for line in f:
        parts = line.decode('utf-8').split()
        key = parts[0]
        SAMPLES[key] = [parts[1], parts[2], parts[3]]

rule all:
        input:
                expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_CatRNAVirus_checkv"],"quality_summary.tsv"),sample=SAMPLES),
                expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.CatRNAVirus_nt.blastn"),sample=SAMPLES),
                expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"{sample}_virus_sofa.txt"),sample=SAMPLES),
                expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["cox1_findhost"],"{sample}.centrifuge.KSout.Ssum"),sample=SAMPLES),
                expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_rdrp"],"{sample}.finish.combine"),sample=SAMPLES),
                #expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit"],"final.virfinder"),sample=SAMPLES),
                expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["closespecies"],"{sample}.contigspecies_info"),sample=SAMPLES),
                expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.blast"),sample=SAMPLES)



##########################################01.QC-fastp-soapnuke-prinseq-rmrRNA############################################
rule QC_fastp:
        input:
                raw_fq1 = lambda wildcards: SAMPLES[wildcards.sample][0],
                raw_fq2 = lambda wildcards: SAMPLES[wildcards.sample][1]
        output:
                fastp_fq1 = temp(os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "{sample}.fastp_1.fq.gz")),
                fastp_fq2 = temp(os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "{sample}.fastp_2.fq.gz")),
                fastp_html = os.path.join(config["output"]["relative"],"{sample}", config["output"]["fastp"], "{sample}.report.html"),
                fastp_json = os.path.join(config["output"]["relative"],"{sample}", config["output"]["fastp"], "{sample}.report.json"),
                raw_fq1_stat = protected(os.path.join(config["output"]["relative"], "{sample}", config["output"]["stats"], "raw_1_seqkit_stats.txt")),
                raw_fq2_stat = protected(os.path.join(config["output"]["relative"], "{sample}", config["output"]["stats"], "raw_2_seqkit_stats.txt")),
                fastp_fq1_stat = protected(os.path.join(config["output"]["relative"], "{sample}", config["output"]["stats"], "fastp_1_seqkit_stats.txt")),
                fastp_fq2_stat = protected(os.path.join(config["output"]["relative"], "{sample}", config["output"]["stats"], "fastp_2_seqkit_stats.txt"))
        params:
                fastp = config["softwares"]["fastp"],
                fastp_threads = config["fastp"]["threads"],
                f_adapter = config["fastp"]["f_adapter"],
                r_adapter = config["fastp"]["r_adapter"],
                quality = config["fastp"]["q"],
                un_qualified = config["fastp"]["u"],
                length_required = config["fastp"]["length_required"]

        log:
                os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "logs/{sample}.fastp.logs") 
        benchmark:
                os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "benchmarks/{sample}.fastp.benchmark.txt")
        shell:
                '''
                {params.fastp} -w {params.fastp_threads} \
                -i {input.raw_fq1} -I {input.raw_fq2} \
                --adapter_sequence {params.f_adapter} \
                --adapter_sequence_r2 {params.r_adapter} \
                --detect_adapter_for_pe \
                -q {params.quality} \
                --length_required {params.length_required} \
                -n 2 -y -c -p \
                --disable_trim_poly_g \
                -o {output.fastp_fq1} -O {output.fastp_fq2} \
                -h {output.fastp_html} -j {output.fastp_json}
                seqkit stats -j 2 -a  {input.raw_fq1} > {output.raw_fq1_stat}
                seqkit stats -j 2 -a  {input.raw_fq2} > {output.raw_fq2_stat}
                seqkit stats -j 2 -a  {output.fastp_fq1} > {output.fastp_fq1_stat}
                seqkit stats -j 2 -a  {output.fastp_fq2} > {output.fastp_fq2_stat}
                '''


rule QC_SOAPnuke_V2:
        input:
                fastp_fq1 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["fastp"],"{sample}.fastp_1.fq.gz"),
                fastp_fq2 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["fastp"],"{sample}.fastp_2.fq.gz")
        output:
                soapnuke_fq1 = temp(os.path.join(config["output"]["relative"],"{sample}",config["output"]["soapnuke"],"{sample}_1.fq.gz")),
                soapnuke_fq2 = temp(os.path.join(config["output"]["relative"],"{sample}",config["output"]["soapnuke"],"{sample}_2.fq.gz")),
                soapnuke_fq1_stat = os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"soapnuke_1_seqkit_stats.txt"),
                soapnuke_fq2_stat = os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"soapnuke_2_seqkit_stats.txt")
        log:
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["soapnuke"],"{sample}.logs")
        benchmark:
                os.path.join(config["output"]["relative"],"{sample}/benchmark/QC/{sample}.soapnuke.benchmark.txt")
        params:
                soapnuke_V2 = config["softwares"]["soapnuke_V2"],
                f_adapter = config["soapnuke_V2"]["f_adapter"],
                r_adapter = config["soapnuke_V2"]["r_adapter"],
                configfile = config["soapnuke_V2"]["configfile"],
                rg= os.path.join(config["output"]["relative"],"{sample}",config["output"]["soapnuke"]),
                cg="{sample}_1.fq.gz",
                dg="{sample}_2.fq.gz",
                c = config["soapnuke_V2"]["c"]
        shell:
                """
                {params.soapnuke_V2} filter -T 4 \
                -1 {input.fastp_fq1} \
                -2 {input.fastp_fq2} \
                -f {params.f_adapter} \
                -r {params.r_adapter} \
                -C {params.cg} -D {params.dg} \
                -c {params.c} \
                -l 20 \
                -q 0.2 \
                -n 0.02 \
                -4 50 \
                -o {params.rg} 
                seqkit stats -j 2 -a  {output.soapnuke_fq1} > {output.soapnuke_fq1_stat}
                seqkit stats -j 2 -a  {output.soapnuke_fq2} > {output.soapnuke_fq1_stat}
                """




rule QC_prinseq:
        input:
                soapnuke_fq1 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["soapnuke"],"{sample}_1.fq.gz"),
                soapnuke_fq2 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["soapnuke"],"{sample}_2.fq.gz")
        output:
                prinseq_fq1 = temp(os.path.join(config["output"]["relative"],"{sample}",config["output"]["prinseq"],"{sample}_good_out_R1.fastq.gz")),
                prinseq_fq2 = temp(os.path.join(config["output"]["relative"],"{sample}",config["output"]["prinseq"],"{sample}_good_out_R2.fastq.gz")),
                prinseq_fq1_stat = os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"prinseq_1_seqkit_stats.txt"),
                prinseq_fq2_stat = os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"prinseq_2_seqkit_stats.txt")
        log:
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["prinseq"],"{sample}.logs")
        benchmark:
                os.path.join(config["output"]["relative"],"{sample}/benchmark/QC/{sample}.prinseq.benchmark.txt")
        params:
                prinseq = config["softwares"]["prinseq"],
                workdir = os.path.join(config["output"]["relative"],"{sample}",config["output"]["prinseq"],"{sample}")
        shell:
                '''
                {params.prinseq} \
                -threads 4 \
                -fastq {input.soapnuke_fq1} -fastq2 {input.soapnuke_fq2} \
                -lc_entropy=0.5 \
                -lc_dust=0.5 \
                -out_gz \
                -out_good {output.prinseq_fq1} \
                -out_good2 {output.prinseq_fq2} \
                -out_single /dev/null -out_single2 /dev/null -out_bad /dev/null -out_bad2 /dev/null 
                seqkit stats -j 2 -a  {output.prinseq_fq1} > {output.prinseq_fq1_stat}
                seqkit stats -j 2 -a  {output.prinseq_fq1} > {output.prinseq_fq2_stat}
                '''



rule QC_SortMeRNA_V4:
        input:
                prinseq_fq1 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["prinseq"],"{sample}_good_out_R1.fastq.gz"),
                prinseq_fq2 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["prinseq"],"{sample}_good_out_R2.fastq.gz")
        output:
                sortmerna_fq1 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rmrRNA_fwd.fq.gz"),
                sortmerna_fq2 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rmrRNA_rev.fq.gz"),
                sortmerna_fq1_stat = os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"sortmerna_1_seqkit_stats.txt"),
                sortmerna_fq2_stat = os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"sortmerna_2_seqkit_stats.txt"),
                sortmerna_rRNA_fq1 = temp(os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rRNA_fwd.fq.gz")),
                sortmerna_rRNA_fq2 = temp(os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rRNA_rev.fq.gz")),
        log:
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"rmrRNA.log")
        benchmark:
                os.path.join(config["output"]["relative"],"{sample}/benchmark/QC/{sample}.rmrRNA1.benchmark.txt")
        params:
                sortmerna_V4 = config["softwares"]["sortmerna_V4"],
                rRNA = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rRNA"),
                rmrRNA = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rmrRNA"),
                workdir = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"]),
                rfam_5p8s = config["database"]["SortMeRNA"]["rfam_5p8s"],
                rfam_5s = config["database"]["SortMeRNA"]["rfam_5s"],
                silva_arc_16s = config["database"]["SortMeRNA"]["silva_arc_16s"],
                silva_arc_23s = config["database"]["SortMeRNA"]["silva_arc_23s"],
                silva_bac_16s = config["database"]["SortMeRNA"]["silva_bac_16s"],
                silva_bac_23s = config["database"]["SortMeRNA"]["silva_bac_23s"],
                silva_euk_18s = config["database"]["SortMeRNA"]["silva_euk_18s"],
                silva_euk_28s = config["database"]["SortMeRNA"]["silva_euk_28s"],
                indexdir = config["database"]["SortMeRNA"]["indexdir"]
        shell:
                """
                {params.sortmerna_V4} \
                --ref {params.rfam_5p8s} \
                --ref {params.rfam_5s} \
                --ref {params.silva_arc_16s} \
                --ref {params.silva_arc_23s} \
                --ref {params.silva_bac_16s} \
                --ref {params.silva_bac_23s} \
                --ref {params.silva_euk_18s} \
                --ref {params.silva_euk_28s} \
                --reads {input.prinseq_fq1} --reads {input.prinseq_fq2} \
                --workdir {params.workdir} --idx-dir {params.indexdir} \
                --fastx --aligned {params.rRNA} --other {params.rmrRNA}  \
                --no-best --num_alignments 1 --paired_out  --out2 --zip-out --threads 4
                seqkit stats -j 2 -a  {input.prinseq_fq1} > {output.sortmerna_fq1_stat}
                seqkit stats -j 2 -a  {input.prinseq_fq2} > {output.sortmerna_fq2_stat}
                """



rule QC_Check:
        input:
                sortmerna_fq1_stat = os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"sortmerna_1_seqkit_stats.txt"),
                sortmerna_fq2_stat = os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"sortmerna_2_seqkit_stats.txt")
        output:
                check_QC = os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"check_QC.txt"),
                rmrRNA_done = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"rmrRNA_done.txt"),
        shell:
                """
                seqkit=/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/bin/seqkit

                pid_1=$(cat {input.sortmerna_fq1_stat}|wc -l)
                pid_2=$(cat {input.sortmerna_fq2_stat}|wc -l)

                if [ -f {output.check_QC} ];
                then
                        rm -f {output.check_QC}
                fi

                touch {output.rmrRNA_done}
                #判断一端是否完整，如果完整，将1输出到output1，如果不完整就删除output1，输出到output0
                if [ $pid_1 -eq 2 ];
                then
                        echo -e "1" >> {output.rmrRNA_done}
                else
                        echo -e "1" >> {output.check_QC}
                        rm -f {output.rmrRNA_done}
                fi

                if [ $pid_2 -eq 2 ];
                then
                        echo -e "2" >> {output.rmrRNA_done}
                else
                        echo -e "2" >> {output.check_QC}
                        rm -f {output.rmrRNA_done}
                fi

                if  [[ $pid_1 -eq 2 ]] $$ [[ $pid_2 -eq 2 ]];
                then
                        touch {output.rmrRNA_done}
                else
                        rm -f {output.rmrRNA_done}
                fi
                touch {output.check_QC}
                echo end at:`date`
                """


##############################################02.assemble-megahit/trinity/metaspades##############################################
rule Assembly_megahit:
        input:
                sortmerna_fq1 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rmrRNA_fwd.fq.gz"),
                sortmerna_fq2 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rmrRNA_rev.fq.gz"),
                done = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"rmrRNA_done.txt")
        output:
                fasta = os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit"],"final.contigs.fa"),
                fasta_stat = os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"megahit_seqkit_stats.txt"),
                fasta_fx2tab = os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"megahit_seqkit_fx2tab.txt")
        benchmark:
                os.path.join(config["output"]["relative"],"{sample}/benchmark/assemble/{sample}.megahit.benchmark.txt")
        params:
                megahit = config["softwares"]["megahit"],
                workdir = os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit"]),
                interdir = os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit"],"intermediate_contigs")
        shell:
                '''
                {params.megahit} \
                -1 {input.sortmerna_fq1} -2 {input.sortmerna_fq2} -o {params.workdir} -f
                seqkit stats -j 2 -a  {output.fasta} > {output.fasta_stat}
                seqkit fx2tab -j 2 -l -n -i {output.fasta} > {output.fasta_fx2tab}
                rm -rf {params.interdir}
                '''



##############################################contigs注释##############################################
rule Annotation_CAT_RefseqNR:
        input:
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["blastx"],"{sample}.blast.fasta")
        output:
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.contig2classification.txt"),
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.ORF2LCA.txt"),
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.predicted_proteins.faa"),
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.contig.tax"),
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.contig.official.tax"),
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.ORF.tax"),
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.summary.txt"),
                temp(os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.alignment.diamond"))
        benchmark:
                os.path.join(config["output"]["relative"],"{sample}/benchmark/down/{sample}.CAT_megahit1.benchmark.txt")
        params:
                CAT = config["softwares"]["CAT"],
                prodigal = config["softwares"]["prodigal"],
                diamond = config["softwares"]["diamond"],
                workdir = os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}"),
                db = config["database"]["cat_db"],
        log:
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"CAT.log")
        threads:4
        shell:
                '''
                {params.CAT} contigs \
                -n 4 --force --block_size 5.0 --index_chunks 1 -c {input} \
                -o {params.workdir} --no_stars \
                -d {params.db} \
                -t /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_tax/ \
                --path_to_prodigal {params.prodigal} \
                --path_to_diamond {params.prodigal}

                {params.CAT} add_names -i {output[0]} -o {output[4]} \
                -t /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_tax/ --only_official
                {params.CAT} add_names -i {output[0]} -o {output[3]} \
                -t /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_tax/
                {params.CAT} add_names -i {output[1]} -o {output[5]} \
                -t /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_tax/ --only_official

                /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_pack/CAT summarise \
                -c {input[0]} -i {output[4]} -o {output[6]} 
                '''


rule Annotation_get_RNAVirus:
        input:
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["blastx"],"{sample}.blast.fasta"),
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.contig.tax"),
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.blast")
        output:
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.CatRNARdRpContigs.fasta"),
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.CatRNARdRpContigs.id"),
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.CatRNAVirusContigs.fasta"),
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.CatRNAVirusContigs.id"),
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.CatContigsProtein.blastp2rdrp.fasta"),
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.CatContigsProtein.blastp2rdrp.id"),
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.virus.contig.tax"),
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.virus.contig.info"),
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.CAT_virus_failed"),
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.CAT_RNA_virus_failed")
        shell:
                '''
                seqkit=/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/bin/seqkit
                #获取CAT的信息
                echo start CAT_extract_visual at:`date`
                grep 'Viruses (superkingdom)' {input[1]} |awk -v OFS='\t' -v FS='\t' '{{print $0}}'|sort|uniq > {output[6]}

                touch {output[6]}

                pid_1=$(cat {output[6]}|wc -l)
                if [ $pid_1 -eq 0 ];
                then
                        echo -e "virus" >> {output[8]}
                else
                        python3 /hwfssz5/ST_INFECTION/GlobalDatabase/share/wangyifei_handover/genome_type/CAT_info_simplify.py \
                        /hwfssz5/ST_INFECTION/GlobalDatabase/share/wangyifei_handover/genome_type/ICTV_genometype.index {output[6]} {output[7]}
                fi


                echo end CAT_extract_visual at:`date`
                #获取比对到RdRp的contigs
                echo start RNA_RdRp.fasta at:`date`
                awk '{{print $1}}' {input[2]}|awk -F _ '{{print $1"_"$2}}'|sort|uniq > {output[5]}
                $seqkit grep -f {output[5]} {input[0]} > {output[4]}

                #获取CAT注释为RNA病毒的contigs
                awk '$4 ~ /RNA/ {{print $1}}' {output[7]}|sort|uniq > {output[3]}
                touch {output[3]}

                pid_2=$(cat {output[3]}|wc -l)
                if [ $pid_2 -eq 0 ];
                then
                        echo -e "RNA_virus" >> {output[9]}
                else
                        $seqkit grep -f {output[3]} {input[0]} > {output[2]}
                fi

                #对CAT注释为RNA病毒的contigs和比对到RdRp的contigs取交集
                sort {output[5]} {output[3]} |uniq -d > {output[1]}
                #$seqkit grep -f {output[1]} {input[0]} > {output[0]}

                touch {output[2]}
                touch {output[0]}
                touch {output[1]}
                touch {output[5]}
                touch {output[4]}
                touch {output[8]}
                touch {output[9]}
                echo end RNA_RdRp.fasta at:`date`
                '''


rule RNAVirus_Checkvgenome:
        input:
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.CatRNAVirusContigs.fasta")
        output:
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_CatRNAVirus_checkv"],"quality_summary.tsv")  
        params:
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_CatRNAVirus_checkv"])
        shell:
                '''
                /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/envs.multi-user/checkv/bin/checkv end_to_end {input[0]} {params[0]} \
                -d /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/envs.multi-user/checkv/checkv-db-v0.6
                '''


######################################## BLASTP RdRp #############################################################


rule RNAVirus_Blastp_VirusRdRp:
        input:
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.predicted_proteins.faa")
        output:
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.blast")
        benchmark:
                os.path.join(config["output"]["relative"],"{sample}/benchmark/down/{sample}.rdrp_blastp.benchmark.txt")
        shell:
                '''
                /hwfssz5/ST_INFECTION/GlobalDatabase/user/fengqikai/software/diamond blastp \
                -d /hwfssz5/ST_INFECTION/GlobalDatabase/share/public_database/rdrp_AAseq/rdrp.mBioANDrefseq \
                -q {input} -k 1 -f 6 -e 1e-5 -o {output[0]}


                '''


##################################### BLASTN NT for False Positive Removal #########################################
rule RNAVirus_Blastn_nt:
        input:
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.CatRNAVirusContigs.fasta")
        output:
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.CatRNAVirus_nt.blastn"),
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.CatRNAVirus_nt.blastn2")
        benchmark:
                os.path.join(config["output"]["relative"],"{sample}/benchmark/down/{sample}.blastn_nt.benchmark.txt")
        shell:
                '''
                export BLASTDB=/ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/liqian6/blastdb_20201224/
                /hwfssz5/ST_INFECTION/GlobalDatabase/user/liqian6/tools/ncbi-blast-2.10.1+/bin/blastn \
                -db /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/liqian6/blastdb_20201224/nt  \
                -query {input[0]}  -out {output[0]} -num_threads 32 \
                -outfmt "6 qacc sacc qlen slen pident evalue bitscore mismatch staxids sscinames scomnames sblastnames sskingdoms" \
                -max_target_seqs 1
                /hwfssz5/ST_INFECTION/GlobalDatabase/user/liqian6/tools/ncbi-blast-2.10.1+/bin/blastn \
                -db /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/liqian6/blastdb_20201224/nt  \
                -query {input[0]}  -out {output[1]} -num_threads 32 \
                -outfmt "6 qacc sacc qlen slen pident evalue bitscore mismatch staxids sscinames scomnames sblastnames sskingdoms"
                '''