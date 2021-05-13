import pandas as pd
import numpy as np
import csv
import subprocess
import os

#curl 192.168.1.10/api/report/sample-csv?gdna_data_id=

def combine_info():
    id_list = []
    title = ['Data_id','type','gene','chr','NM','exon','base','AA']
    for indexs in df2.index:
        [data_id,var_type,gene,gene_chr,nm,exon,base,aa,vaf_tumor,vaf_blood,cosmic,start,end,ref,alt,copy_number_tumor,copy_number_blood,virus,virus_vaf,sample] = df2.loc[indexs].values#参数读取
        if sample == '血液':
            sample = 'blood'
        else:
            sample = 'tumor'
        
        if isinstance(aa,str):
            if '*' in aa[-1]:
                aa = aa[:-1] + 'X'
        #----------------------title------------------------------
        if str(data_id) + '_' + sample not in id_list:
            id_list += [str(data_id) + '_' + sample]
            if sample == 'tumor' and (vaf_blood > -1 or copy_number_blood > -1):
                id_list += [str(data_id) + '_' + 'blood']

        #------------------------病毒数据保存---------------------
        if isinstance(virus,float) == False:
            if virus not in total['virus']:
                total['virus'][virus] = {'virus':virus}
            total['virus'][virus][str(data_id) + '_' + sample] = virus_vaf
        #----------------------过滤空行----------------------------
        if isinstance(var_type,float):
            continue
        #---------------------扩增数据保存-------------------------
        if var_type == 'cnv':
            if gene not in total['cnv']:
                total['cnv'][gene]= {'type':'cnv','gene':gene,'base':base}
            if '/' in aa:
                total['cnv'][gene][str(data_id) + '_tumor'] = df2.loc[indexs].values[15]
                total['cnv'][gene][str(data_id) + '_blood'] = df2.loc[indexs].values[16]
            else:
                total['cnv'][gene][str(data_id) + '_' + sample] = str(aa)
        #-----------------------突变，融合数据保存----------------------
        else:
            list1 = ['Data_id','type','gene','chr','NM','exon','base','AA','cosmic','start','end','ref','alt']
            list2 = [data_id,var_type,gene,gene_chr,nm,exon,base,aa,cosmic,start,end,ref,alt]
            info = dict(map(lambda x,y:[x,y],list1,list2))
            target = gene+'.*'+exon+'.*'+base+'.*'+aa
            if target not in total[var_type]:
                total[var_type][target] = info
            if vaf_tumor >= 0:
                if vaf_blood >= 0:
                    total[var_type][target][str(data_id) + '_tumor'] = df2.loc[indexs].values[8]
                    total[var_type][target][str(data_id) + '_blood'] = df2.loc[indexs].values[9]
                else:
                    total[var_type][target][str(data_id) + '_' + sample] = vaf_tumor
            else:
                total[var_type][target][str(data_id) + '_blood'] = vaf_blood
        
    id_list = sorted(id_list) #ID排序
    title += id_list + [ffpe_id,'cosmic','start','end','ref','alt']
    df = pd.DataFrame(columns=title)

    return df

def Snv():
    for i in total['snv']:
        if '?' in i.split('.*')[-1]:
            (tmp_error,tmp_data) = subprocess.getstatusoutput("grep {S1}.*{S2} {S3}".format(S1=i.split('.*')[0],S2=i.split('.*')[-2].split('>')[0],S3=snv_file))
        else:
            (tmp_error,tmp_data) = subprocess.getstatusoutput("grep {S1}.*{S2} {S3}".format(S1=i.split('.*')[0],S2=i.split('.*')[-1],S3=snv_file))
        lines = tmp_data.split('\n')
        if lines == ['']:
            continue
        for j in lines:
            if float(j.split('%')[1].split(':')[-1]) > float(j.split('%')[0].split(':')[-1])*2:
                total['snv'][i][ffpe_id] = j.split('%')[1].split(':')[-1]

def Indel():
    for i in total['indel']:
        (tmp_error,tmp_data) = subprocess.getstatusoutput("grep {S1}.*{S2}: {S3}".format(S1=i.split('.*')[0],S2=i.split('.*')[1],S3=indel_file))
        lines = tmp_data.split('\n')
        for j in lines:
            if j.split('\t')[1:5] == [str(total['indel'][i]['start']),str(total['indel'][i]['end']),total['indel'][i]['ref'],total['indel'][i]['alt']]:
                if float(j.split('%')[1].split(':')[-1]) > float(j.split('%')[0].split(':')[-1])*2:
                    total['indel'][i][ffpe_id] = j.split('%')[1].split(':')[-1]
    '''
    with open('tmp_indel.csv','w',newline='') as w:
        writer = csv.writer(w)
        writer.writerows(output_indel_tmp)
    os.system('Rscript vepSomatic.R tmp_indel.csv')

    with open('tmp_indel.utf8vep.csv','r') as indel_f:
        data = csv.reader(indel_f)
        for i in data:
            target = i[2]+'.*'+i[5]+'.*'+i[6]+'.*'+i[7]
            if target in total['indel']:
                total['indel'][target][sample_id][:8]
                for j in id_list:
                    if j in total['indel'][target]:
                        tmp.append(total['indel'][target][j])
                    else:
                        tmp.append('')
                    #if j not in total['indel'][i]:
                if i[8] == '':
                    tmp.append(i[9])
                else:
                    tmp.append(i[8])
                tmp += i[-5:]
                output_indel.append(tmp)
    with open('before_indel.utf8vep.csv','a') as indel_w:
        writer = csv.writer(indel_w)
        writer.writerows(output_indel_before)
    return output_indel
    '''

def Virus():
    with open(virus_file,'r') as f:
        data = f.readlines()
        for i in data:
            if i.split('\t')[0] in total['virus']:
                total['virus'][i.split('\t')[0]] = i.split('\t')[1][:-3]

def Cnv():    #写入CNV数据到新表 （完成）
    with open(cnv_file,'r') as f:
        data = f.read().splitlines()
        for i in data:
            if i.split(',')[0] in total['cnv']:
                total['cnv'][i.split(',')[0]][ffpe_id] = i.split(',')[1]

def get_cnv():
    for i in os.listdir('/x01_haplox/hapreports/0.上机信息汇总/'):
        if '_share.xlsx' in i:
            info_sheet = '/x01_haplox/hapreports/0.上机信息汇总/' + i
    df = pd.read_excel(info_sheet)
    base_list = []
    id_list = []
    cnv_list = [cnv_file]
    with open(before_file) as f:
        data = f.readlines()
        for i in data:
            id_list.append(i.split(',')[0])
        id_list = id_list[1:]
    for i in df['建库信息'].tolist():
        if isinstance(i,float) == False and '_' in i:
            if i.split('_')[-1] in id_list:
                base_list.append([df[df.建库信息 == i].index.tolist()[0],i])
    for i in base_list:
        if 'cfdna' in i[1]:
            get_dir = '/haplox/rawout/{S1}/{S2}/cnv/'.format(S1=i[0],S2=i[1])
        else:
            get_dir = '/haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/cnv/'.format(S1=i[0],S2=i[1])
        for j in os.listdir(get_dir):
            if 'genes188in680.csv' in j:
                before_cnv_file = get_dir + j
            elif '_rg_cnv_result.csv' in j:
                before_cnv_file = get_dir + j
        if before_cnv_file not in cnv_list:
            cnv_list.append(before_cnv_file)
    
    df = pd.read_csv(cnv_list[0])
    df.columns = list(map(lambda x: x + '_' + cnv_list[0].split('/')[-3].split('_')[-1],df.columns))

    for i in cnv_list[1:]:
        df1 = pd.read_csv(i)
        df1.columns = list(map(lambda x: x + '_' + i.split('/')[-3].split('_')[-1],df1.columns))
        df = df.join(df1, how='outer')
    df.to_csv('defense_cnv.csv',index=False)

def main():
    #---------------合并扩增信息------------------------
    get_cnv()
    #---------------合并守护信息------------------------
    df = combine_info()

    if len(total['cnv']) > 0:
        Cnv()

    if len(total['virus']) > 0:
        Virus()


    if len(total['snv']) > 0:
        Snv()

    if len(total['indel']) > 0:
        Indel()

    #--------------------output---------------------
    for i in total:
        for j in total[i]:
            df = df.append([total[i][j]], ignore_index=True)
    df.to_csv("defense_var.csv", index=False)


location = os.getcwd()
sample_id = location.split('/')[-1]
ffpe_id = location.split('_')[-1]
for i in os.listdir(location):
    if 'gdna' in i and 'bai' in i:
        gdna_id = i.split('_')[-2]
    if 'indel_MrBam.txt' in i:
        indel_file = location + '/' + i
    if 'snv_MrBam.txt' in i:
        snv_file = location + '/' + i
os.system('curl 192.168.1.10/api/report/sample-csv?gdna_data_id={S1} > sample-{S2}.csv'.format(S1=gdna_id,S2=ffpe_id))
before_file = 'sample-{S2}.csv'.format(S2=ffpe_id)
cnv_file = location + '/cnv/{S1}_rg_cnv.genes188in680.csv'.format(S1=sample_id)
virus_file = location + '/virus/virus_tumor_normal_result.txt'
fusion_file = location + '/fusionscan/{S1}_fusion.csv'.format(S1=sample_id)
#-------------参数初始化-----------------
total = {'snv':{},'indel':{},'cnv':{},'fusion':{},'virus':{}}

df2 = pd.read_csv(before_file,encoding='utf-8')
if len(df2.index) == 0 : #判断是否有连续样本
    main()
