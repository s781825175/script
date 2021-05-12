#coding:utf-8
#!/bin/bash
#Rscript /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/nobias_getMrBam_txt.R /haplox/rawout/181103_A00250_0061_AH7YN3DSXX/{S2} {S3}
import os
import argparse
import xlwt
import csv
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir", nargs='?', action="store" , type=str, help="annokb sample directory")
args = parser.parse_args()

def txt_xls(in_file,out_file):
    f=open(in_file,"rt")
    x=0
    y=0
    xls=xlwt.Workbook()
    sheet = xls.add_sheet('sheet1',cell_overwrite_ok=True)
    while True:
          line = f.readline()
          if not line:
             break
          for i in line.split('\t'):
             item=i.strip().decode('utf8')
             sheet.write(x,y,item)
             y+=1
          x+=1
          y=0
    f.close()
    xls.save(out_file)

def cnv(cnv_x):
    filename = cnv_x
    with open(filename) as c:
        reader = csv.reader(c)
        cnv = list(reader)
        cnv_num = '0'   
        for i in cnv:
            if i[1] != 'cnv':
                cnv_num = cnv_num + '+' + i[1]
        avg = (eval(cnv_num)/(len(cnv)-1))
    return avg

def Annokb(path,hap,dat,gdna,num):
    if hap != gdna:
        os.system('Rscript /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/getFFPE_cfDNA_snv_indel_v2.R \
        /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/ffpe_vs_pbl/genes464.bed \
        {S1}/{S3}.snv-nobias-GB18030-baseline.csv T \
        {S1}/{S3}.snv-nobias-GB18030-baseline-genes464.csv 0.1'.format(S1=path,S3=dat))

        os.system('Rscript /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/getFFPE_cfDNA_snv_indel_v2.R \
        /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/ffpe_vs_pbl/genes464.bed \
        {S1}/{S3}.indel-nobias-GB18030-baseline.csv T \
        {S1}/{S3}.indel-nobias-GB18030-baseline-genes464.csv 0.1'.format(S1=path,S3=dat))
    else:
        os.system('Rscript /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/getFFPE_cfDNA_snv_indel_v2.R \
        /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/ffpe_vs_pbl/genes464.bed \
        {S1}/{S3}_snv-GB18030-baseline.csv T \
        {S1}/{S3}.snv-nobias-GB18030-baseline-genes464.csv 0.1'.format(S1=path,S3=dat))

        os.system('Rscript /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/getFFPE_cfDNA_snv_indel_v2.R \
        /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/ffpe_vs_pbl/genes464.bed \
        {S1}/{S3}_indel-GB18030-baseline.csv T \
        {S1}/{S3}.indel-nobias-GB18030-baseline-genes464.csv 0.1'.format(S1=path,S3=dat))


    #-----------------------------------------------------------------------------------------      
    os.system('python /thinker/net/tools/Annokb.py \
    -f {S1}/{S3}.indel-nobias-GB18030-baseline-genes464.csv  \
    -t indel -o {S1}/Annokb_mrbam_{S3}.indel-nobias-GB18030-baseline-genes464.csv'.format(S1=path,S3=dat))
    if 'ffpe' in path:
        sample = 'ffpe'
    elif gdna == hap:
        sample = 'ffpe'
    else:
        sample = 'cf'
        


    os.system('python /thinker/net/tools/Annokb.py \
    -f {S1}/{S3}.snv-nobias-GB18030-baseline-genes464.csv  \
    -t snv -o {S1}/Annokb_mrbam_{S3}.snv-nobias-GB18030-baseline-genes464.csv'.format(S1=path,S3=dat))

    os.system('python /haplox/users/yangbo/futionbase.py -f {S1}/fusionscan/{S2}_fusion.json -b {S1}/{S2}_rg.bam'.format(S1=path,S2=hap))
    #-----------------------------------------------------------------------------------------          
                                                      
    os.system('perl /haplox/users/ZhouYQ/germline/bin/Gene2Disease.pl /haplox/users/ZhouYQ/germline/bin/Database/female_cancer.list \
    {S1}/germline/result/{S4}.germline.txt \
    {S1}/germline/result/{S4}.cancer.female.txt \
    {S1}/germline/result/{S4}.nocancer.female.txt \
    {S1}/germline/result'.format(S1=path,S4=gdna)) 
                                                    
    os.system('perl /haplox/users/ZhouYQ/germline/bin/Gene2Disease.pl /haplox/users/ZhouYQ/germline/bin/Database/male_cancer.list \
    {S1}/germline/result/{S4}.germline.txt \
    {S1}/germline/result/{S4}.cancer.male.txt \
    {S1}/germline/result/{S4}.nocancer.male.txt \
    {S1}/germline/result '.format(S1=path,S4=gdna))

    os.system('perl /haplox/users/liaowt/Script/germline/germline_trans_v2.pl \
    {S1}/germline/result//{S4}.cancer.female.txt  \
    {S1}/germline/result//{S4}_trans.cancer.female.txt'.format(S1=path,S4=gdna)) 

    os.system('perl /haplox/users/liaowt/Script/germline/germline_trans_v2.pl \
    {S1}/germline/result//{S4}.cancer.male.txt  \
    {S1}/germline/result//{S4}_trans.cancer.male.txt'.format(S1=path,S4=gdna))
    

    cnv_list = path + '/cnv/' + hap + '_rg_cnv.chrX_genes.csv'
    avg = cnv(cnv_list)
    os.system('curl haplab.haplox.net/api/report/chemotherapy/{S5} -F "qc=@{S1}/germline/result/{S4}.information.txt" '.format(S1=path,S2=hap,S4=gdna,S5=num))
    os.system('curl haplab.haplox.net/api/report/chemotherapy/{S5} -F "chem=@{S1}/germline/result/{S4}.chem_451.txt" '.format(S1=path,S2=hap,S4=gdna,S5=num))
    if avg >= 1.5:
        os.system("Rscript /haplox/tools/mhgvs/myvep/vepMerge.R {S1}/germline/result/{S4}_trans.cancer.female.txt {S1}/germline/result/sample.filter.hg19_multianno.txt".format(S1=path,S2=hap,S4=gdna,S5=num))
        os.system("Rscript /haplox/tools/mhgvs/myvep/vepMerge.R {S1}/germline/result/{S4}_trans.cancer.female.txt {S1}/germline/result/gdna.filter.hg19_multianno.txt".format(S1=path,S2=hap,S4=gdna,S5=num))
        os.system('curl haplab.haplox.net/api/report/chemotherapy/{S5} -F "germline=@{S1}/germline/result/{S4}_trans.cancer.female.utf8vep.txt" '.format(S1=path,S2=hap,S4=gdna,S5=num))
    if avg < 1.5:
     
        os.system("Rscript /haplox/tools/mhgvs/myvep/vepMerge.R {S1}/germline/result/{S4}_trans.cancer.male.txt {S1}/germline/result/sample.filter.hg19_multianno.txt".format(S1=path,S2=hap,S4=gdna,S5=num))
        os.system("Rscript /haplox/tools/mhgvs/myvep/vepMerge.R {S1}/germline/result/{S4}_trans.cancer.male.txt {S1}/germline/result/gdna.filter.hg19_multianno.txt".format(S1=path,S2=hap,S4=gdna,S5=num))
        os.system('curl haplab.haplox.net/api/report/chemotherapy/{S5} -F "germline=@{S1}/germline/result/{S4}_trans.cancer.male.utf8vep.txt" '.format(S1=path,S2=hap,S4=gdna,S5=num))

    return sample,avg

def cnv_updata(path,num,sample,tmp,avg):
    title = ["data_id","gene","result","tumor","cfdna","type","chr","nm","region"]
    amplification_black = ['BRCA2','MLH1','PMS2']
    loss_black = ['CSF1R','FGF1','NPM1','PDGFRB','FGFR4','CCND3','CDKN1A','VEGFA','IRS2','FLT4','DAXX']
    cnv_dir = path + '/' + tmp + '/'
    a=0
    ERBB2 = False
    for i in os.listdir(cnv_dir):
        if i.endswith("_rg_cnv_result.csv") and i.startswith("S"):
            with open(cnv_dir + i) as cnv_list:
                point_reader = csv.reader(cnv_list)
                for row in point_reader:
                    if a==0:
                        a += 1
                        df = pd.DataFrame ([title] , columns = title)
                        df.to_csv(cnv_dir + "result_cnv.csv",encoding="utf-8",index=False, header=False)
                        continue
                    if float(row[1]) >= 3 and row[0] not in amplification_black:
                        if row[0] == 'ERBB2':
                            ERBB2 = True
                        if sample == 'ffpe':
                            inside = [num,row[0],"扩增",row[1],'',"cnv",row[-3],row[-2],row[-1]]
                        else:
                            inside = [num,row[0],"扩增",'',row[1],"cnv",row[-3],row[-2],row[-1]]
                        df = pd.DataFrame ([inside] , columns = title)
                        df.to_csv(cnv_dir + "result_cnv.csv",encoding="utf-8",index=False, header=False, mode='a')
                        a += 1
                    if float(row[1]) <= 1.3 and row[0] not in loss_black:
                        if avg < 1.5 and row[-3] == 'chrX': 
                            continue
                        if sample == 'ffpe':
                            inside = [num,row[0],"缺失",row[1],'',"cnv",row[-3],row[-2],row[-1]]
                        else:
                            inside = [num,row[0],"缺失",'',row[1],"cnv",row[-3],row[-2],row[-1]]
                        df = pd.DataFrame ([inside] , columns = title)
                        df.to_csv(cnv_dir + "result_cnv.csv",encoding="utf-8",index=False, header=False, mode='a')
                        a += 1
    if ERBB2 == True and tmp == 'cnv':
        diffpanel(path,num,sample,avg)
    elif a > 1:
        os.system('curl haplab.haplox.net/api/report/csv?type=cnv -F "import_file=@{S1}result_cnv.csv" '.format(S1=cnv_dir))

def diffpanel(path,num,sample,avg):
    os.mkdir('{S1}/cnv_single'.format(S1=path))
    os.system('/haplox/users/huang/myGit/cnvkit/rawPython/bin/python /haplox/users/huang/myGit/cnvkit/cnv.py -s {S1}/{S2}_rg.bam -b /haplox/ref/bed/605panel_v2_primary_targets_gene.bed -r /haplox/ref/GATK/ucsc.hg19/ucsc.hg19.fasta -o {S1}/cnv_single'.format(S1=path,S2=path.split('/')[-1]))
    os.system('Rscript /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/getCnv.R /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/cnv/cnv.csv  {S1}/cnv_single/{S2}_rg_cnv_result.txt  {S1}/cnv_single/{S2}_rg_cnv_result.csv'.format(S1=path,S2=path.split('/')[-1]))
    os.system('cp -r {S1}/cnv_single {S1}/{S2}'.format(S1=path,S2=num))
    tmp = 'cnv_single'
    cnv_updata(path,num,sample,tmp,avg)

def main():
    path = os.getcwd()
    if args.dir:
        path = args.dir
    hap = path.split('/')[-1]
    dat = hap.split('_')[1] + '_' + hap.split('_')[-1]
    gdna_dir = path + '/germline/result/'
    for i in os.listdir(gdna_dir):
        if i.endswith('chem_451.txt'):
            gdna = i.split('.')[0]
    num = hap.split('_')[-1]
    sample,avg = Annokb(path,hap,dat,gdna,num)
    target = ['utf8vep.txt','_trans.cancer.male.txt','information.txt','chem_451.txt','Target_451.txt','_trans.cancer.female.txt']
    for i in os.listdir(gdna_dir):
        for j in target:
            if i.endswith(j):
                pat = gdna_dir + '/' + i              
                if j in target[1:4]:
                    txt_xls(pat,(pat[:-3]+'xls'))
                if j not in target[1:4]:
                    txt_xls(pat,(pat.replace('.','_')[:-4] + '.xls'))
    if sample == 'ffpe':  #把血液栏的丰度转移到组织栏
        for i in os.listdir(path):
            if 'Annokb_mrbam' in i:
                tmp = []
                Annokb_file =  open(path + '/' + i,'r')
                data = csv.reader(Annokb_file)
                for j in data:
                    if j[1] == 'gene':
                        tmp.append(j)
                        continue
                    j[7],j[8]=j[8],j[7]
                    tmp.append(j)
                Annokb_file.close()
                Annokb_file = open(path + '/2_' + i,'w')
                writer = csv.writer(Annokb_file)
                writer.writerows(tmp)
                Annokb_file.close()    
    os.system('Rscript /haplox/users/yangbo/605_last_pair_mutscan_fusion_virus_sheet.R {S1}/'.format(S1=path))
    tmp = 'cnv'
    cnv_updata(path,num,sample,tmp,avg)

if __name__ == "__main__":
    main()

