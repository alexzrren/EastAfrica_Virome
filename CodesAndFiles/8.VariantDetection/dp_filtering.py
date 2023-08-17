import os
from tqdm import tqdm 
import pandas as pd
import sys


with open('requant_group.tsv') as fd:
    requant = fd.read().splitlines()
    
requant_vsnames = set([ group.split()[1] for group in requant ])
print(requant_vsnames)

def parsing_fbvcf(path):
    '''
    获取freebayes vcf的
    '''
    fbsites = []
    with open(path) as fd:
        for line in fd.readlines():
            if line.startswith('#'):
                continue
            else:
                ref, pos, *_ = line.split()
                fbsites.append(pos)
    return fbsites


def vcf2dataframe(vcfpath):
    '''
    读取一个vcf转换成dataframe 列与vcfheader的列一致
    '''
    for line in open(vcfpath).readlines():
        if line.startswith('#CHROM'):
            vcfdfheader = line.lstrip('#').split()
            break
    vcfdf = pd.read_table(vcfpath, comment='#', header=None, names=vcfdfheader)
    return vcfdf


def vcf2varsitesdict(vcfpath):
    '''
    读取vcf文件的变异
    '''
    varsitedict = dict()
    vcfdf = vcf2dataframe(vcfpath)
    for sample in list(vcfdf.columns)[9:]:
        varsitedict[sample] = list(vcfdf.loc[vcfdf[sample].apply(lambda x: x.startswith('1:')), 'POS'])
    return varsitedict
        

def compare_sites(altisnv_sites, fbvcfdict):
    '''
    比较isnv和
    '''
    output_dict = dict()
    for sample, sitelist in fbvcfdict.items():
        isec_sites = set(map(int, altisnv_sites)) & set(sitelist)
        if len(isec_sites) >= 3:
            output_dict[sample] = len(isec_sites)
    return output_dict


fbvcf_sitedict = dict()
fbvcf_sitepersamp_dict = dict()
for vsname in requant_vsnames:
    vcfpath = os.path.join('/jdfssz1/ST_HEALTH/P20Z10200N0206/renzirui/EA_redo/ViralStrain_PopGenetics', '3.variantcalling', vsname+'.vcf')
    try:
        fbvcf_sites = parsing_fbvcf(vcfpath)
        varsitedict = vcf2varsitesdict(vcfpath)
    except FileNotFoundError:
        fbvcf_sites = []
        varsitedict = dict()
    
    fbvcf_sitedict[vsname] = fbvcf_sites
    fbvcf_sitepersamp_dict[vsname] = varsitedict
print(fbvcf_sitepersamp_dict)
    

resultlist = []
for group in (requant):
    sn, vsname, refid = group.split()
    lofreq_vcfpath = os.path.join('/jdfssz1/ST_HEALTH/P20Z10200N0206/renzirui/EA_redo/ViralStrain_PopGenetics/8.lofreq_result', vsname, sn+'.lofreq.vcf')
    if not os.path.exists(lofreq_vcfpath):
        print(lofreq_vcfpath +': file not found!')
        continue
    output_isnvpath = os.path.join('/jdfssz1/ST_HEALTH/P20Z10200N0206/renzirui/EA_redo/ViralStrain_PopGenetics/9.lofreq_iSNV', vsname)
    if not os.path.exists(output_isnvpath):
        os.mkdir(output_isnvpath)
    output_content = ''

    
    dpsum, altisnv, refisnv, sitenum = 0,0,0,0
    var_direction = { 'A->G':0, 'G->A':0, 'T->C':0, 'C->T':0 }
    altisnv_sites = []
    for line in open(lofreq_vcfpath).readlines():
        if line.startswith('#'):
            output_content += line
        else:
            chrom, pos, siteid, ref, alt, qual, filter, info = line.strip().split()
            sitenum += 1
            info_dict = {  item.split('=')[0]:item.split('=')[1] for item in info.split(';') }
            dpsum += int(info_dict['DP'])
            info_dict['DP4'] = list(map(int, info_dict['DP4'].split(',')))
            if int(info_dict['DP'])>50 and ((info_dict['DP4'][0] + info_dict['DP4'][1]) / int(info_dict['DP']))>0.05 and ((info_dict['DP4'][2] + info_dict['DP4'][3]) / int(info_dict['DP']))>0.05:
                if (info_dict['DP4'][0] + info_dict['DP4'][1]) < (info_dict['DP4'][2] + info_dict['DP4'][3]):
                    isnv_type = 'REF'
                    refisnv += 1
                else:
                    isnv_type = 'ALT'
                    altisnv += 1
                    altisnv_sites.append(pos)
                    var_type = ref + '->' + alt
                    if var_type in var_direction.keys():
                        var_direction[var_type] += 1
                    
                info_dict['DP4'] = ','.join(list(map(str,info_dict['DP4'])))
                
                info = ';'.join([ k+'='+v for k,v in info_dict.items() ]) + ';iSNV='+isnv_type
                lineout = '\t'.join([chrom, pos, siteid, ref, alt, qual, filter, info])
                output_content += lineout+'\n'
    
    if sitenum == 0:
        continue
    compare_result = compare_sites(altisnv_sites, fbvcf_sitepersamp_dict[vsname])
    compare_result = dict(sorted(compare_result.items(), key=lambda x: x[1], reverse=True))
    if len(compare_result):
        print(vsname, sn, compare_result)
    resultlist.append([sn, refid, vsname, dpsum/(sitenum), refisnv+altisnv, refisnv, altisnv, len(set(fbvcf_sitedict[vsname]) & set(altisnv_sites)),len(fbvcf_sitedict[vsname])] + list(var_direction.values()) + [  var_direction['T->C']+var_direction['C->T'],  var_direction['A->G']+var_direction['G->A'] ])
    # with open(os.path.join(output_isnvpath, sn+'.iSNV.vcf'), 'w') as fdout:
        # fdout.write(output_content)


resultdf = pd.DataFrame(resultlist, columns=['sn', 'refid','vsname', 'meandp', 'iSNV', 'REF_iSNV', 'ALT_iSNV', 'intersect_sites', 'fbsites'] + list(var_direction.keys()) + ['T<->C', 'A<->G'])


basenum = pd.read_table('refseq_baseperc.tsv', header=None, names=['refid', 'length', r'%A', r'%T', r'%C', r'%G', r'%GC', r'%AT'])
#resultdf = resultdf.merge(basenum, on='refid', how='left')


resultdf.to_excel('output.xlsx', index=None)