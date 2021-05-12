#coding=utf-8
import csv
import io
import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--ffpedna", nargs='?', action="store" , type=str, help="ffpedna file name")
parser.add_argument("-c", "--cfdna", nargs='?', action="store" , type=str, help="cfdna file name")
parser.add_argument("-o", "--output", nargs='?', action="store" , type=str, help="output file name")
args = parser.parse_args()
file1 = args.ffpedna
file2 = args.cfdna
test = args.output
importfile = (args.output[::-1].split('/',1))[1][::-1] + '/importfile.csv' 
path = args.ffpedna[::-1].split('/',2)[-1][::-1]
num = path.split('_')[-1]

def fusion_updata(path,fusionfile,num): # fusion上传到系统后台
    title = ['data_id','gene','result','tumor','cfdna','type']
    result_gain = []
    output = []
    output_list = {}
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

                if i[0] == num:
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

    for i in result_gain:
        if i[2] not in output_list:
            output_list[i[2]] = i
        else:
            if i[3] == '':
                output_list[i[2]][4] = i[4]
            else:
                output_list[i[2]][3] = i[3]
    for i in output_list:
        output.append(output_list[i])
    
    if switch == True:
        df = pd.DataFrame ([title] , columns = title)
        df.to_csv('{S1}/fusionscan/fusion_gain.csv'.format(S1=path),encoding="utf-8",index=False, header=False)
        for row in output:
            df = pd.DataFrame ([row])
            df.to_csv('{S1}/fusionscan/fusion_gain.csv'.format(S1=path),encoding="utf-8",index=False, header=False, mode='a')
        #os.system('curl haplab.haplox.net/api/report/csv?type=fusion -F "import_file=@{S1}/fusionscan/fusion_gain.csv"'.format(S1=path))

tmp = []
tmp2 = []
tmp3 = []
with open(file1 , 'r',  encoding="gb2312") as ffpe:
    fusion1 = csv.reader(ffpe)
    for row in fusion1:
        if len(row) <13:
            if int(row[1]) > 1:
                tmp3.append(row)
            t_row = row
            t_row[0]='ffpedna'
            tmp2.append(t_row)
            continue
        if row[13].isdigit():
            if int(row[1]) > 1:
                tmp3.append(row)
            if int(row[13]) <1000:
                tmp.append(row)

with open(file2, 'r',  encoding="gb2312") as ffpe2:
    fusion2 = csv.reader(ffpe2)
    for row in fusion2:
        if len(row) <13:
            if int(row[1]) > 1:
                tmp3.append(row)
            t_row = row
            t_row[0]='cfdna'
            tmp2.append(t_row)
            continue
        if row[13].isdigit():
            if int(row[1]) > 1:
                tmp3.append(row)
            if int(row[13]) <1000:
                tmp.append(row)

with open(test, 'w',  newline='', encoding="gb2312") as out:
    csv_write = csv.writer(out,dialect='excel')
    t = ['hapnum','level','left_gene','right_gene','left_vaf','right_vaf','total','unique','left_exon','left_exonid','right_exon','right_exonid','strand','baseline_count','total_average','unique_average','health_baseline','health_total','health_unique','blood_baseline','blood_total','blood_unique']
    csv_write.writerow(t)         
    for i in tmp:
        csv_write.writerow(i)
    for i in tmp2:
        csv_write.writerow(i)

with open(importfile, 'a',  newline='', encoding="gb2312") as out2:
    csv_write = csv.writer(out2,dialect='excel')
    t = ['hapnum','level','left_gene','right_gene','left_vaf','right_vaf','total','unique','left_exon','left_exonid','right_exon','right_exonid','strand','baseline_count','total_average','unique_average','health_baseline','health_total','health_unique','blood_baseline','blood_total','blood_unique']
    csv_write.writerow(t)
    for i in tmp3:
        csv_write.writerow(i)

fusion_updata(path,test,num)
