#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 22 22:32:55 2023

@author: zirui
"""

import pandas as pd
from glob import glob
import BioinfoPy.FileIO as BioIO
import BioinfoPy.Utility as Bioutils


msapaths = glob('/Users/zirui/BGI_Projects/WIV_EastAfrica_Clean/ViralProtein/EA_PhyloMSA_trimal/*')

vsname2seqid = pd.read_excel('/Users/zirui/BGI_Projects/WIV_EastAfrica/Picante_test/Assembly_sequences.xlsx')
vs2seqid_dict = vsname2seqid.set_index('PhylogenicTree_SeqID').to_dict(orient='dict')['vsname']


result = []
for msapath in msapaths:
    fastadict = BioIO.fasta2dict(msapath)
    vfamily = msapath.split('/')[-1].split('_')[0]
    markeridlist = list(set(vs2seqid_dict.keys()) & set(fastadict.keys()))
    for markerid in markeridlist:
        markerseq = fastadict[markerid]
        msalen = len(markerseq)
        min_pmismatch = 9999999
        min_seqid = ''
        for seqid in fastadict.keys():
            if seqid == markerid:
                continue
            else:
                _, mismatch, gap, _, _ = Bioutils.calc_aai(markerseq, fastadict[seqid])
                if (msalen - gap) == 0:
                    continue
                pismatch = mismatch / (msalen - gap) * 100
                if pismatch < min_pmismatch:
                    min_pmismatch = pismatch
                    if seqid in markeridlist:
                        seqid = vs2seqid_dict[seqid]
                        
                    min_seqid = seqid
        result.append([vfamily, vs2seqid_dict[markerid], min_seqid, min_pmismatch, msalen])
        
        
resultdf = pd.DataFrame(result, columns=['vfamily', 'vsname', 'nearest_seqid', 'nearest_pmismatch', 'msalen'])
resultdf.to_excel('/Users/zirui/BGI_Projects/WIV_EastAfrica_Clean/ViralProtein/AAI_output.xlsx')



result = []
for msapath in msapaths:
    fastadict = BioIO.fasta2dict(msapath)
    vfamily = msapath.split('/')[-1].split('_')[0]
    markeridlist = list(set(vs2seqid_dict.keys()) & set(fastadict.keys()))
    for markerid in markeridlist:
        markerseq = fastadict[markerid]
        msalen = len(markerseq)
        min_pmismatch = 9999999
        min_seqid = ''
        for seqid in set(fastadict.keys()) - set(markeridlist):
            _, mismatch, gap, _, _ = Bioutils.calc_aai(markerseq, fastadict[seqid])
            if (msalen - gap) == 0:
                continue
            pismatch = mismatch / (msalen - gap) * 100
            if pismatch < min_pmismatch:
                min_pmismatch = pismatch
                if seqid in markeridlist:
                    seqid = vs2seqid_dict[seqid]
                    
                min_seqid = seqid
        result.append([vfamily, vs2seqid_dict[markerid], min_seqid, min_pmismatch, msalen])

resultdf = pd.DataFrame(result, columns=['vfamily', 'vsname', 'nearest_seqid', 'nearest_pmismatch', 'msalen'])
resultdf.to_excel('/Users/zirui/BGI_Projects/WIV_EastAfrica_Clean/ViralProtein/AAI_pubonly.xlsx')