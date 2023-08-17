for hspec in `cat homosamples.txt | awk '{print $3}' | sort -u`
do
	#echo $hspec
	for gc in `cat homosamples.txt | grep $hspec | awk '{print $2}' | sort -u`
	do
		gcnum=`echo $gc | sed 's/G//g'`
		for gc2 in `cat homosamples.txt | grep $hspec | awk '{print $2}' | sort -u`
		do
		gcnum2=`echo $gc2 | sed 's/G//g'`
		if [ $gcnum -lt $gcnum2 ];then
			
			cat homosamples.txt | grep $hspec | grep -w $gc | awk '{print $1}' > keep1.indv
			cat homosamples.txt | grep $hspec | grep -w $gc2 | awk '{print $1}' > keep2.indv
			#vcftools --vcf VCF_afterQC_forHetCalc/$hspec.filtered.snp.vcf --weir-fst-pop keep1.indv --weir-fst-pop keep1.indv --out tmp #--fst-window-size 5000 --fst-window-step 1000 --stdout # | sed 1d | sed 's/^/'"$hspec"'\t'"$gc"'-'"$gc2"'\t/g'
			echo $hspec $gc $gc2 `wc -l keep1.indv | awk '{print $1}'` `wc -l keep2.indv | awk '{print $1}'` `vcftools --vcf VCF_afterQC_perBatSpecies/$hspec.vcf --weir-fst-pop keep1.indv --weir-fst-pop keep2.indv --out tmp 2>&1 | grep 'weighted Fst' | awk -F': ' '{print $2}'`
		fi		
		done
	done
done