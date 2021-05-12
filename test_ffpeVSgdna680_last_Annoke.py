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

def grep_sample(fi,dire,sample_name):
    filter_data = []
    clean_data = []
    with open(dire + '/' +  fi,'r') as f:
        data = csv.reader(f)
        for i in data:
            if i[1] == 'gene':
                clean_data.append(i)
                filter_data.append(i)
                continue
            print("grep '{S1}.*{S2}' {S3}/*/*hg19*txt".format(S1=i[1],S2=i[6],S3=dire[::-1].split('/',1)[-1][::-1]))
            tmp_data = os.popen("grep '{S1}.*{S2}' {S3}/*/*hg19*txt".format(S1=i[1],S2=i[6],S3=dire[::-1].split('/',1)[-1][::-1]))
            tmp_data = tmp_data.readlines()
            sample = []
            order = []
            
            for j in tmp_data:
                for k in j.split(':'):
                    if '%' in k:
                        if sample_name in j.split(':')[0]:
                            sample.append(float(k[:-1]))
                        else:
                            order.append(float(k[:-1]))
            if order == []:
                clean_data.append(i)
                continue

            elif int(sample[1]) > 10 or len(tmp_data) < 5 or (int(sample[1]) > 5 and sum(order[1:][::2])/len(order[::2]) < 2 ):
                clean_data.append(i)
                continue
            else:
                filter_data.append(i)
                continue
    with open(dire + '/clean_{S1}'.format(S1=fi),'w') as w:
        writer = csv.writer(w)
        writer.writerows(clean_data)    

    with open(dire + '/filter_{S1}'.format(S1=fi),'w') as w:
        writer = csv.writer(w)
        writer.writerows(filter_data)    

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
        os.system('Rscript /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/nobias_getMrBam_txt_V3.R {S1} {S3}'.format(S1=path,S3=dat))
        os.system('Rscript /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/getFFPE_cfDNA_snv_indel_v3.R \
        /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/ffpe_vs_pbl/genes528in680.bed \
        {S1}/{S3}.snv-nobias-GB18030-baseline.csv T \
        {S1}/{S3}.snv-nobias-GB18030-baseline-genes528in680.csv 0.1'.format(S1=path,S3=dat))

        os.system('Rscript /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/getFFPE_cfDNA_snv_indel_v3.R \
        /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/ffpe_vs_pbl/genes528in680.bed \
        {S1}/{S3}.indel-nobias-GB18030-baseline.csv T \
        {S1}/{S3}.indel-nobias-GB18030-baseline-genes528in680.csv 0.1'.format(S1=path,S3=dat))
    else:
        os.system('Rscript /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/getFFPE_cfDNA_snv_indel_v2.R \
        /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/ffpe_vs_pbl/genes528in680.bed \
        {S1}/{S3}_snv-GB18030-baseline.csv T \
        {S1}/{S3}.snv-nobias-GB18030-baseline-genes528in680.csv 0.5'.format(S1=path,S3=dat))

        os.system('Rscript /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/getFFPE_cfDNA_snv_indel_v2.R \
        /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/ffpe_vs_pbl/genes528in680.bed \
        {S1}/{S3}_indel-GB18030-baseline.csv T \
        {S1}/{S3}.indel-nobias-GB18030-baseline-genes528in680.csv 0.5'.format(S1=path,S3=dat))


    #-----------------------------------------------------------------------------------------      
    os.system('python /thinker/net/tools/Annokb.py \
    -f {S1}/{S3}.indel-nobias-GB18030-baseline-genes528in680.csv  \
    -t indel -o {S1}/Annokb_mrbam_{S3}.indel-nobias-GB18030-baseline-genes528in680.csv'.format(S1=path,S3=dat))

    os.system('python /thinker/net/tools/Annokb.py \
    -f {S1}/{S3}.snv-nobias-GB18030-baseline-genes528in680.csv  \
    -t snv -o {S1}/Annokb_mrbam_{S3}.snv-nobias-GB18030-baseline-genes528in680.csv'.format(S1=path,S3=dat))
    
    if os.path.exists(path + '/' + 'Annokb_mrbam_{S3}.snv-nobias-GB18030-baseline-genes528in680.csv'.format(S3=dat)):
        grep_sample('Annokb_mrbam_{S3}.snv-nobias-GB18030-baseline-genes528in680.csv'.format(S3=dat),path,hap)

    os.system('python /haplox/users/yangbo/script/combine_delins.py {S1}/Annokb_mrbam_{S3}.indel-nobias-GB18030-baseline-genes528in680.csv {S1}/Annokb_mrbam_{S3}.snv-nobias-GB18030-baseline-genes528in680.csv'.format(S1=path,S3=dat))


    os.system('Rscript /haplox/tools/mhgvs/myvep/vepSomatic.R {S1}/Annokb_mrbam_{S3}.snv-nobias-GB18030-baseline-genes528in680.csv'.format(S1=path,S3=dat))

    os.system('Rscript /haplox/tools/mhgvs/myvep/vepSomatic.R {S1}/Annokb_mrbam_{S3}.indel-nobias-GB18030-baseline-genes528in680.csv'.format(S1=path,S3=dat))

    os.system('python /haplox/users/yangbo/futionbase.py -f {S1}/fusionscan/{S2}_fusion.json -b {S1}/{S2}_rg.bam'.format(S1=path,S2=hap))

    #-----------------------------------------------------------------------------------------          

def germline_updata(path,hap,gdna,num):
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

    os.system('perl /haplox/users/wenger/script/germline_trans_v2.pl \
    {S1}/germline/result//{S4}.cancer.female.txt  \
    {S1}/germline/result//{S4}_trans.cancer.female.txt'.format(S1=path,S4=gdna)) 

    os.system('perl /haplox/users/wenger/script/germline_trans_v2.pl \
    {S1}/germline/result//{S4}.cancer.male.txt  \
    {S1}/germline/result//{S4}_trans.cancer.male.txt'.format(S1=path,S4=gdna))

    os.system("Rscript /haplox/tools/mhgvs/myvep/vepMerge.R {S1}/germline/result/{S4}_trans.germline_680.cancer.txt {S1}/germline/result/gdna.filter.hg19_multianno.txt".format(S1=path,S2=hap,S4=gdna,S5=num))
    os.system("Rscript /haplox/tools/mhgvs/myvep/vepMerge.R {S1}/germline/result/{S4}_trans.germline_DDR.cancer.txt {S1}/germline/result/gdna.filter.hg19_multianno.txt".format(S1=path,S2=hap,S4=gdna,S5=num))
    
    cnv_list = path + '/cnv/' + hap + '_rg_cnv.chrX_genes.csv'
    avg = cnv(cnv_list)

    os.system('curl haplab.haplox.net/api/report/chemotherapy/{S5} -F "qc=@{S1}/germline/result/{S4}.information.txt" '.format(S1=path,S2=hap,S4=gdna,S5=num))
    os.system('curl haplab.haplox.net/api/report/chemotherapy/{S5} -F "germline=@{S1}/germline/result/{S4}_trans.germline_680.cancer.utf8vep.txt"'.format(S1=path,S2=hap,S4=gdna,S5=num))
    os.system('curl haplab.haplox.net/api/report/chemotherapy/{S5} -F "chem=@{S1}/germline/result/{S4}.chem_451.txt" '.format(S1=path,S2=hap,S4=gdna,S5=num))
#    os.system('curl haplab.haplox.net/api/report/chemotherapy/{S5} -F "germline=@{S1}/germline/result/{S4}_trans.germline_680.cancer.txt" '.format(S1=path,S2=hap,S4=gdna,S5=num))
    os.system('curl 192.168.1.10/api/report/hla-v?data_id={S5} -F "import_file=@{S1}/HLA/{S4}_sortbam-hla.csv" '.format(S1=path,S2=hap,S4=gdna,S5=num))
#    os.system('curl 192.168.1.10/api/report/hrr?data_id={S5} -F "import_file=@{S1}/germline/result/{S4}_trans.germline_DDR.cancer.txt" '.format(S1=path,S2=hap,S4=gdna,S5=num))
    os.system('curl haplab.haplox.net/api/report/hrr?data_id={S5} -F "import_file=@{S1}/germline/result/{S4}_trans.germline_DDR.cancer.utf8vep.txt"'.format(S1=path,S2=hap,S4=gdna,S5=num))
    return avg

def cnv_updata(path,cnv_file,num,sample,avg):
    title = ["data_id","gene","result","tumor","cfdna","type","chr","nm","region"]
    amplification_black = ['BRCA2','MLH1','PMS2']
    loss_black = ['CSF1R','FGF1','NPM1','PDGFRB','FGFR4','CCND3','CDKN1A','VEGFA','IRS2','FLT4','DAXX']
    with open(cnv_file) as cnv_list:
        point_reader = csv.reader(cnv_list)
        a = 0
        for row in point_reader:
            if a==0:
                a += 1
                df = pd.DataFrame ([title] , columns = title)
                df.to_csv(cnv_file[:-4] + '_gain.csv',encoding="utf-8",index=False, header=False)
                continue
            if float(row[1]) >= 3 and row[0] not in amplification_black:
                if sample == 'ffpe':
                    inside = [num,row[0],"扩增",row[1],'',"cnv",row[-3],row[-2],row[-1]]
                else:
                    inside = [num,row[0],"扩增",'',row[1],"cnv",row[-3],row[-2],row[-1]]
                df = pd.DataFrame ([inside] , columns = title)
                df.to_csv(cnv_file[:-4] + '_gain.csv',encoding="utf-8",index=False, header=False, mode='a')
                a += 1
            if float(row[1]) <= 1.3 and row[0] not in loss_black:
                if avg < 1.5 and row[-3] == 'chrX': 
                    continue
                if sample == 'ffpe':
                    inside = [num,row[0],"缺失",row[1],'',"cnv",row[-3],row[-2],row[-1]]
                else:
                    inside = [num,row[0],"缺失",'',row[1],"cnv",row[-3],row[-2],row[-1]]
                df = pd.DataFrame ([inside] , columns = title)
                df.to_csv(cnv_file[:-4] + '_gain.csv',encoding="utf-8",index=False, header=False, mode='a')
                a += 1
    if a > 1:
        os.system('curl haplab.haplox.net/api/report/csv?type=cnv -F "import_file=@{S1}" '.format(S1=cnv_file[:-4] + '_gain.csv'))

def fusion_updata(path,fusionfile,num,sample): # fusion上传到系统后台
    title = ['data_id','gene','result','tumor','cfdna','type']
    result_gain = []
    hot_gene = ['ALK','ROS1','RET','NRG1','NTRK1','PAX7','FGFR2','KMT2A','ETV6','NAB2','DDIT3','NUTM1','FUS','RARA','SS18','PRKACA','PAX8','PAX3','TMPRSS2','EWSR1','PDGFB','FGFR3','BRAF','MYC','NR4A3','NTRK2','TFE3','PRKACA']
    with open(fusionfile,'r') as f:
        data = csv.reader(f)
        switch = False
        for i in data:
            if i[1] == '4':
                switch = True
                tmp = []
                if i[0] != 'fusion':
                    tmp.append(i[0])
                elif len(result_gain) > 1:
                    continue
                else:
                    tmp.append(num)
                if i[2] in hot_gene:
                    tmp.append(i[2])
                else:
                    tmp.append(i[3])
                if i[0] != 'fusion':
                    print(i)
                    if i[12] == 'reversed':
                        tmp_str = '{S1}-exon{S3}-{S2}-exon{S4}融合'.format(S1=i[3],S2=i[2],S3=eval(i[11]),S4=eval(i[9]+'+1'))
                        tmp.append(tmp_str)
                    else:
                        tmp_str = '{S1}-exon{S3}-{S2}-exon{S4}融合'.format(S1=i[2],S2=i[3],S3=eval(i[9]),S4=eval(i[11]+'+1'))
                        tmp.append(tmp_str)
                else:
                    if 'not' in i[9] or 'not' in i[11]:
                        tmp.append('{S1}-exon{S3}-{S2}-exon{S4}融合'.format(S1=i[2],S2=i[3],S3=i[9],S4=i[11]))
                    else:
                        tmp.append('{S1}-exon{S3}-{S2}-exon{S4}融合'.format(S1=i[2],S2=i[3],S3=i[9].split('-')[1],S4=i[11].split('-')[1]))
                if sample == 'ffpe':
                    if i[2] in hot_gene:
                        tmp.append(i[4])
                    else:
                        tmp.append(i[5])
                    tmp.append('')
                else:
                    tmp.append('')
                    if i[2] in hot_gene:
                        tmp.append(i[4])
                    else:
                        tmp.append(i[5])
                tmp.append('fusion')
                result_gain.append(tmp)
    if switch == True:
        df = pd.DataFrame ([title] , columns = title)
        df.to_csv('{S1}/fusionscan/fusion_gain.csv'.format(S1=path),encoding="utf-8",index=False, header=False)
        for row in result_gain:
            df = pd.DataFrame ([row])
            df.to_csv('{S1}/fusionscan/fusion_gain.csv'.format(S1=path),encoding="utf-8",index=False, header=False, mode='a')
        os.system('curl haplab.haplox.net/api/report/csv?type=fusion -F "import_file=@{S1}/fusionscan/fusion_gain.csv"'.format(S1=path))

def main():
    path = os.getcwd()
    if args.dir:
        path = args.dir
    hap = path.split('/')[-1]
    dat = hap.split('_')[1] + '_' + hap.split('_')[-1]
    gdna_dir = path + '/germline/result/'
    for i in os.listdir(path + '/germline/result/'):
        if i.endswith('chem_451.txt'):
            gdna = i.split('.')[0]
    num = hap.split('_')[-1]

    if 'ffpe' in path:
        sample = 'ffpe'
    elif gdna == hap:
        sample = 'ffpe'
    else:
        sample = 'cf'
    
    Annokb(path,hap,dat,gdna,num)
    avg = germline_updata(path,hap,gdna,num)
    target = ['_trans.germline_DDR.cancer.utf8vep.txt','_trans.germline_680.cancer.utf8vep.txt','_trans.cancer.male.txt','information.txt','chem_451.txt','Target_451.txt','_trans.cancer.female.txt','_trans.germline_680.txt','_trans.germline_DDR.cancer.txt']
    for i in os.listdir(gdna_dir):
        for j in target:
            if i.endswith(j):
                pat = gdna_dir + '/' + i              
                if j in target[1:4]:
                    txt_xls(pat,(pat[:-3]+'xls'))
                if j not in target[1:4]:
                    txt_xls(pat,(pat.replace('.','_')[:-4] + '.xls'))

    cnvfile = '{S1}/cnv/{S2}_rg_cnv.genes188in680.csv'.format(S1=path,S2=hap)
    cnv_updata(path,cnvfile,num,sample,avg)
    fusionfile = '{S1}/fusionscan/{S2}_fusion.csv'.format(S1=path,S2=hap)
    fusion_updata(path,fusionfile,num,sample)
    os.system('Rscript /haplox/users/yangbo/680_last_pair_mutscan_fusion_virus_sheet.R {S1}/'.format(S1=path))

if __name__ == "__main__":
    main()


