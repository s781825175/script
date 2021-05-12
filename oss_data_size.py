import subprocess
import re
import csv
import os

def data_size(data):
    pattern = re.compile(r' (\d+) ')
    num = 0
    c=0
    with open(data,'r',encoding='utf-8') as f:
        data=csv.reader(f)
        for i in data:
            if len(i) == 1:
                continue
            if '_' not in i[0]:
                (status1, output1) = subprocess.getstatusoutput('ossutil ls oss://sz-hapseq/rawfq/MGI_merge/{S2}/{S3}'.format(S1=i[0][:4],S2=i[0],S3=i[1]))
                a = pattern.findall(output1)
                print('ossutil ls oss://sz-hapseq/rawfq/MGI_merge/{S2}/{S3}'.format(S1=i[0][:4],S2=i[0],S3=i[1]))
            else:
                (status1, output1) = subprocess.getstatusoutput('ossutil ls oss://sz-hapseq/rawfq/20{S1}/{S2}_clinic/{S3}'.format(S1=i[0][:4],S2=i[0],S3=i[1]))
                a = pattern.findall(output1)
                if a == []:
                    print('ossutil ls oss://sz-hapseq/rawfq/20{S1}/{S2}/{S3}'.format(S1=i[0][:4],S2=i[0],S3=i[1]))
                    (status1, output1) = subprocess.getstatusoutput('ossutil ls oss://sz-hapseq/rawfq/20{S1}/{S2}/{S3}'.format(S1=i[0][:4],S2=i[0],S3=i[1]))
                    a = pattern.findall(output1)
                else:
                    print('ossutil ls oss://sz-hapseq/rawfq/20{S1}/{S2}_clinic/{S3}'.format(S1=i[0][:4],S2=i[0],S3=i[1]))
            if a == []:
                print('error')
                continue
            num += int(a[0])
            num += int(a[1])
                #if 'gdna' in i[1]:
                   # continue
            c += 1        
    print(num,c)

def data_download(data):
    pattern = re.compile(r'oss:.*?.fastq.gz')
    with open(data,'r',encoding='utf-8') as f:
        data=csv.reader(f)
        num = 0
        for i in data:
            if len(i) == 1:
                print(i)
                continue
            num += 1
            (status1, output1) = subprocess.getstatusoutput('ossutil ls oss://sz-hapseq/rawfq/20{S1}/{S2}_clinic/{S3}'.format(S1=i[0][:4],S2=i[0],S3=i[1]))
            a = pattern.findall(output1)
            #print(a)
            for j in a:
                print('ossutil cp ' + j + ' .')
                os.system('ossutil cp ' + j + ' ./xu')
        print(num)

data = 'ck3.csv'
data_size(data)
#data_download(data)
