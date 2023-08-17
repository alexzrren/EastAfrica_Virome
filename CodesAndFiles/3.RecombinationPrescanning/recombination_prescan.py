#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 21:26:00 2022

@author: zirui
"""
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import BioinfoPy.FileIO as BioIO
import BioinfoPy.Utility as Bioutils
import os
import sys
import numpy as np
from tqdm import tqdm
from recan.simgen import Simgen
#sns.set_theme(style='darkgrid')

def sseqid_substitute(seqid, vsnamedict, ledict):
    if seqid not in vsnamedict.keys():
        return '%-13s(len=%s)' % (seqid, lendict[seqid])
    else:
        vsname = vsnamedict[seqid]
        return '%-13s(len=%s)' % (vsname, lendict[seqid])
    
    
def assemblyseq_identify(seqid, assemblyids):
    if seqid in assemblyids:
        return 'Assembly'
    else:
        return 'Public'
    

#读入病毒株信息表，完成病毒株ID替换
vsdata = pd.read_excel('/Users/zirui/BGI_Projects/WIV_EastAfrica/Update_Figure&Table/TableS2.Summary of virus strains(add_DNAvirus)v7_dw.xlsx')
vsnamedict = vsdata[['Reference_Sequence', 'ViralStrain_Name']].set_index('Reference_Sequence').to_dict(orient='dict')['ViralStrain_Name']

#读入比对序列组装序列长度表
with open('/Users/zirui/BGI_Projects/WIV_EastAfrica/SlideWindow_BLASTN/subdb_lenstat.tsv') as lenstat:
    lendict = { key: value for key, value in map(lambda x: x.split(), lenstat.read().splitlines()) }

blastresult_path = '/Users/zirui/BGI_Projects/WIV_EastAfrica/SlideWindow_BLASTN/slidewindow_blast/'
blastresults = os.listdir(blastresult_path)
filenum = len(blastresults)
for i, resultfile in enumerate(blastresults):
    print('[%d/%d] processing...' %(i+1, filenum))
    #blastn基础数据读取，处理获取query序列ID以及window信息
    blastdf = BioIO.read_blast(blastresult_path+'/'+resultfile)
    assemblyseq = '_'.join(blastdf.iloc[0,0].split('_')[:-1])
    vsname = vsnamedict[assemblyseq]
    blastdf['windowstart'] = blastdf['qseqid'].apply(lambda x: int(x.split('_')[-1].split(':')[0])-1)
    blastdf['heatmap_rownames'] = blastdf['sseqid'].apply(lambda x: sseqid_substitute(x, vsnamedict, lendict))
    blastdf['SeqGroup'] = blastdf['sseqid'].apply(lambda x: assemblyseq_identify(x, list(vsnamedict.keys())))
    #除去每个组装序列BLASTN结果中的自身比对结果
    subject_seqs = list(set(blastdf['sseqid'].drop_duplicates()) - set([assemblyseq]))
    
    if not len(subject_seqs): #如果除去自身没有其他序列，不进行后续的分析
        continue

    start_locs = list(blastdf['windowstart'].drop_duplicates())
    heatmap_data = []
    bitscore_sumdict = dict()
    heatmap_mask = [] #用于seaborn或matplotlib热图绘制，遮蔽NA点
    
    # sseq: Subject SequenceID
    for sseq in subject_seqs:
        dfsubject = blastdf[blastdf['sseqid']==sseq]
        heatmap_rowname = dfsubject.reset_index().loc[0,'heatmap_rownames']
        subject_identity = [sseq, heatmap_rowname, dfsubject.reset_index().loc[0, 'SeqGroup']]       #记录每个window上的subject的BLASTN一致性
        subject_mask = []               #寻找NA点，记录Mask信息
        bitscore_sumdict[sseq] = 0      #计算每个Subject序列在所有window上的bitscore总和
        for windowstart in start_locs:
            windowdata = dfsubject[dfsubject['windowstart']==windowstart]
            
            if len(windowdata):
                subject_identity.append(windowdata.reset_index().loc[0,'pident'])
                bitscore_sumdict[sseq] += windowdata.reset_index().loc[0,'bitscore']
                subject_mask.append(False)
            else:
                subject_identity.append(None)
                subject_mask.append(True)
                
        heatmap_data.append(subject_identity)
        heatmap_mask.append(subject_mask)
        
    
    #将上述循环提取的信息转换为dataframe或numpy数据结构，用于R/seaborn绘图
    dfheatmap = pd.DataFrame(heatmap_data, columns=['sseqid','heatmap_rownames','SeqGroup']+start_locs)
    mask = np.matrix(heatmap_mask)
    dfheatmap.set_index('sseqid', drop=True, inplace=True)
    
    dfheatmap_posmax = dfheatmap.iloc[:,2:].dropna(axis=1,how='all').replace(np.nan, 0).idxmax()
    dfheatmap_keepseq = dfheatmap_posmax.drop_duplicates()
    draw_heatmap = True
    i = 0
    while len(list(dfheatmap_keepseq)) < 2:
        print('Iteration %d' % (i))
        i += 1
        dfheatmap = dfheatmap.drop(list(dfheatmap_keepseq)[0])
        if len(dfheatmap) < 2:
            print('No recombination event detected![1]')
            draw_heatmap = False
            break
        dfheatmap_posmax = dfheatmap.iloc[:,2:].dropna(axis=1,how='all').replace(np.nan, 0).idxmax()
        dfheatmap_keepseq = dfheatmap_posmax.drop_duplicates()

    
    #生成最大一致性Marker矩阵
    dfmarker = dfheatmap.copy()
    dfmarker.loc[:,:]='NA'
    dfmarker = dfmarker.drop(['heatmap_rownames', 'SeqGroup'], axis=1)
    
    for pos in dfheatmap_posmax.index:
        maxseq = dfheatmap_posmax[pos]
        dfmarker.loc[maxseq ,pos] = '*'
    dfmarker = dfmarker.replace('NA','')
    
    #将热图绘制数据导出为csv文件
    csv_outfile = '/Users/zirui/BGI_Projects/WIV_EastAfrica/SlideWindow_BLASTN/slidewindow_heatmap_data/'+vsname+'.csv'
    csv_outfile_mk = '/Users/zirui/BGI_Projects/WIV_EastAfrica/SlideWindow_BLASTN/slidewindow_heatmap_data/'+vsname+'_marker.csv'
    dfheatmap.to_csv(csv_outfile)
    dfmarker.to_csv(csv_outfile_mk)
    
    #调用R脚本绘图
    #print('/usr/local/bin/Rscript /Users/zirui/BGI_Projects/WIV_EastAfrica/EastAfrica_Plot/script/Recombination_HeatMap.R '+vsname)
    if draw_heatmap:
        (code, stdout, stderr) = Bioutils.runcommand('/usr/local/bin/Rscript /Users/zirui/BGI_Projects/WIV_EastAfrica/EastAfrica_Plot/script/Recombination_HeatMap.R '+vsname)
        if code:
            print('R exit with non-zero code %d' %code)
            print('Error when processing %s\n>ErrorInfo------------------'% vsname)
            print(stdout.decode())
            print(stderr.decode())

        dfheatmap_keepseq = dfheatmap_keepseq.append(pd.Series([assemblyseq]))
        dfheatmap_keepseq.to_csv('.seqids_csv.tmp', index=None, header=None)
        (code, stdout, stderr) = Bioutils.runcommand('/Users/zirui/Softwares/seqkit grep -f .seqids_csv.tmp Sequence_Original/%s_db.fa > Sequence_Selected/%s.fa' % (vsname, vsname))
        if code:
            print('seqkit grep exit with non-zero code %d' %code)
            print('Error when processing %s\n>ErrorInfo------------------'% vsname)
            print(stdout.decode())
            print(stderr.decode())
        
        os.remove('.seqids_csv.tmp')
    
        continue
    
        (code, stdout, stderr) = Bioutils.runcommand('/usr/local/bin/mafft --thread 8 --adjustdirection --auto Sequence_Selected/%s.fa > Sequence_Aligned_Clean/%s.fa' % (vsname, vsname))    
        if code:
            print('mafft exit with non-zero code %d' %code)
            print('Error when processing %s\n>ErrorInfo------------------'% vsname)
            print(stdout.decode())
            print(stderr.decode())
            
            
        if not os.path.exists('Sequence_Aligned_Clean/%s.fa' % (vsname)):
            print('Skip '+vsname)
            continue
        sim_obj = Simgen('Sequence_Aligned_Clean/%s.fa' % (vsname))
        sd = sim_obj.__dict__
        for i in range(len(sd['_records'])):
            if sd['_records'][i].id == assemblyseq:
                assemblyseqindex = i
                break
        
        skiplist = ['PoxV-1A']
        if vsname in skiplist:
            continue
        sim_obj.simgen(window=500, shift=50, pot_rec=assemblyseqindex)
        try:
            distdata_wide = sim_obj.get_data().T
        except:
            #distdata_wide = sim_obj.get_data()
            sys.exit(1)

        fig, ax = plt.subplots()

        for seq in distdata_wide.columns:
            if seq in vsnamedict.keys():
                seqvsname = vsnamedict[seq]
                plt.plot(distdata_wide.index, distdata_wide[seq], lw=0.5, label=seqvsname)
            else:
                plt.plot(distdata_wide.index, distdata_wide[seq], lw=0.5, label=seq)
            print(seq)
        
        plt.axis('on')
        legend = ax.legend(loc = 0, prop = {'size':5})
        frame = legend.get_frame()
        frame.set_alpha(1)
        frame.set_facecolor('none')
        plt.title(vsname)
        ax.spines['left'].set_color('black')
        ax.spines['top'].set_color('black')
        ax.spines['right'].set_color('black')
        ax.spines['bottom'].set_color('black')
        plt.xticks()
        plt.yticks()
        fig.savefig('/Users/zirui/BGI_Projects/WIV_EastAfrica/SlideWindow_BLASTN/Simplot/%s.pdf' % (vsname), transparent=True)
        fig.clf()
        
