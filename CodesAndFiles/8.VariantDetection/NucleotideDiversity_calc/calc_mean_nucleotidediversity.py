#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 23:43:07 2023

@author: zirui
"""

import pandas as pd
import BioinfoPy.FileIO as BioIO
import BioinfoPy.Utility as BioUtils
from tqdm import tqdm
from glob import glob
import statistics as stat


vsinfo = pd.read_excel('/Users/zirui/BGI_Projects/WIV_EastAfrica_Clean/Tables/2_VSinfo_v2.1.xlsx')
vsname2len = vsinfo[['ViralStrain_Name', 'Reference_Length (bp)']].set_index('ViralStrain_Name').to_dict()['Reference_Length (bp)']

requant = pd.read_excel('/Users/zirui/BGI_Projects/WIV_EastAfrica_Clean/Tables/3_Requantification_v1.3.xlsx', index_col=0).reset_index(drop=True)
requant_new = requant[['vsname', 'sn', 'meandepth', 'coverage', 'RPM', 'revised_gps']].sort_values('RPM', ascending=False).drop_duplicates(subset=['vsname', 'sn'], keep='first').copy()

sn2gps = requant[['sn', 'revised_gps']].drop_duplicates().set_index('sn').to_dict()['revised_gps']

paths = glob('/Users/zirui/BGI_Projects/WIV_EastAfrica_Clean/ViralStrain_PopGenome/consensus_SNPseq/*')

result = []
result_full = []
for path in paths:
    vsname = path.split('/')[-1].split('.')[0]
    if vsname == 'ParV-1A':
        continue
    fastadict = BioIO.fasta2dict(path)
    seqids = set(fastadict.keys()) & set(requant_new.loc[requant_new['vsname']==vsname, 'sn'])
    sampsites = set(requant_new.loc[requant_new['vsname']==vsname, 'revised_gps'])
    if len(seqids) < 4 or len(sampsites) < 2:
        continue
    seqids = list(seqids)
    pidentlist = []
    
    samegps_mismatch = []
    diffgps_mismatch = []

    for i in range(len(seqids)):
        for j in range(i+1, len(seqids)):
            seqid1 = seqids[i]
            seqid2 = seqids[j]
            same_gps = sn2gps[seqid1]==sn2gps[seqid2]
            match, mismatch, gap, length, pident = BioUtils.calc_aai(fastadict[seqid1], fastadict[seqid2])
            
            pidentlist.append(mismatch)
            result_full.append([vsname, seqid1, seqid2, mismatch/(length-gap), same_gps])
            if same_gps:
                samegps_mismatch.append(mismatch)
            else:
                diffgps_mismatch.append(mismatch)
            
    mean_pident = stat.mean(pidentlist)
    if len(samegps_mismatch) and len(diffgps_mismatch):
        estimate_fst = (stat.mean(diffgps_mismatch) - stat.mean(samegps_mismatch)) / stat.mean(diffgps_mismatch)
    else:
        continue
    vcfpath = '/Users/zirui/BGI_Projects/WIV_EastAfrica_Clean/ViralStrain_PopGenome/NucleotideDiversity/Popgenome_stat_cSNV/{vsname}_snp.vcf'.format(vsname=vsname)
    code, stdout, stderr = BioUtils.runcommand(r"cat {path} | sed 's/\t0:/\t0\/0:/g' | sed 's/\t1:/\t1\/1:/g' | vcftools --vcf - --TajimaD {length} --out calc_temp".format(path=vcfpath, length=vsname2len[vsname]))
    print(r"cat {path} | sed 's/\t0:/\t0\/0:/g' | sed 's/\t1:/\t1\/1:/g' | vcftools --vcf - --TajimaD {length} --out {vsname}".format(path=path, length=vsname2len[vsname], vsname=vsname))
    if code:
        print('Error exit with code' + str(code))
        print(stderr.decode())
    with open('calc_temp.Tajima.D') as fd:
        n_snp, tajimaD = fd.read().splitlines()[1].split()[-2:]
    
    result.append([vsname, len(seqids), len(sampsites), int(n_snp), mean_pident/vsname2len[vsname]*100, estimate_fst, float(tajimaD)])
    
resultdf = pd.DataFrame(result, columns=['vsname', 'num_individual', 'num_samplesite', 'num_snp', 'mean_diversity', 'estimated_Fst', 'tajimaD'])
#resultdf.to_excel('NucleotideDiversity_kept.xlsx', index=False)

fullresultdf = pd.DataFrame(result_full, columns=['vsname', 'sn1', 'sn2', 'pmismatch', 'same_gps'])
fullresultdf.to_excel('NucleotideDiversity_full.xlsx', index=False)

