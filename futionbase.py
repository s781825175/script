#coding=utf-8
import pymysql
import io
import json
import argparse
import re
import csv
import os
import commands
import logging
import time


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", nargs='?', action="store" , type=str, help="json file name")
parser.add_argument("-d", "--directory", nargs='?', action="store" , type=str, help="directory and subdirectory ")
parser.add_argument("--delete", nargs='?', action="store" , type=str, help="delete hupnum")
parser.add_argument("-b", "--bam", nargs='?', action="store" , type=str, help="rg_bam") 
args = parser.parse_args()

def todb(data):     #传入数据
    plus = [['healthgdna_wesplus','wesgdnabl'],['ffpedna_wesplus','wesplusfb'],['healthgdna','gdnabl'],['healthcfdna','healthcfdb'],['cfdna','cfdnadb']]
    cursor = connect.cursor()
    sql = "INSERT INTO fusiondb (hapnum,time,total_count,unique_count,left_gene,left_chr,left_exon,left_exonid,left_strand,right_gene,right_chr,right_exon,right_exonid,right_strand,left_position,right_position,special) VALUES  ('%d', '%s', '%d','%d', '%s','%s','%s','%d','%s', '%s','%s','%s','%d','%s','%d','%d','%d' )"
    cursor.execute(sql % data)
    connect.commit()
    cursor.close()

    cursor = connect.cursor()
    
    for i in plus:
        if i[0] in jsonfile:
            sql = "INSERT INTO {S1} (hapnum,time,total_count,unique_count,left_gene,left_chr,left_exon,left_exonid,left_strand,right_gene,right_chr,right_exon,right_exonid,right_strand,left_position,right_position) VALUES  ('%d', '%s', '%d','%d', '%s','%s','%s','%d','%s', '%s','%s','%s','%d','%s','%d','%d')".format(S1=i[1])
            break
        else:
            sql = "INSERT INTO ffpednadb (hapnum,time,total_count,unique_count,left_gene,left_chr,left_exon,left_exonid,left_strand,right_gene,right_chr,right_exon,right_exonid,right_strand,left_position,right_position) VALUES  ('%d', '%s', '%d','%d', '%s','%s','%s','%d','%s', '%s','%s','%s','%d','%s','%d','%d')"

    cursor.execute(sql % data[:-1])
    connect.commit()
    cursor.close()
    
def getbaseline(data,tumor):   #从数据库中获取基线
    cursor = connect.cursor()
    sql = "SELECT count(*) FROM fusiondb WHERE left_gene = '%s' AND left_exon = '%s' AND left_exonid = '%d' AND right_gene = '%s' AND right_exon = '%s' AND right_exonid = '%d' "
    sql_avg = "SELECT AVG(unique_count), AVG(total_count) FROM fusiondb WHERE left_gene = '%s' AND left_exon = '%s' AND left_exonid = '%d' AND right_gene = '%s' AND right_exon = '%s' AND right_exonid = '%d' "
    cursor.execute(sql % data)
    count_tu = cursor.fetchall()     #条目数
    count = count_tu[0][0]

    cursor = connect.cursor()
    cursor.execute(sql_avg % data)
    avg_count = cursor.fetchall()
    for row in avg_count:
        unique_count = row[0]   #平均数
        total_count = row[1]
    cursor.close()

    cursor = connect.cursor()
    sql = "SELECT count(*) FROM healthcfdb WHERE left_gene = '%s' AND left_exon = '%s' AND left_exonid = '%d' AND right_gene = '%s' AND right_exon = '%s' AND right_exonid = '%d' "
    sql_avg = "SELECT AVG(unique_count), AVG(total_count) FROM healthcfdb WHERE left_gene = '%s' AND left_exon = '%s' AND left_exonid = '%d' AND right_gene = '%s' AND right_exon = '%s' AND right_exonid = '%d' "
    cursor.execute(sql % data)
    count_tu = cursor.fetchall()     #条目数
    count_hea = count_tu[0][0]

    cursor = connect.cursor()
    cursor.execute(sql_avg % data)
    avg_count = cursor.fetchall()
    for row in avg_count:
        unique_count_hea = row[0]   #平均数
        total_count_hea = row[1]
    cursor.close()

    cursor = connect.cursor()
    if tumor == "ffpe":
        sql = "SELECT count(*) FROM ffpednadb WHERE left_gene = '%s' AND left_exon = '%s' AND left_exonid = '%d' AND right_gene = '%s' AND right_exon = '%s' AND right_exonid = '%d' "
        sql_avg = "SELECT AVG(unique_count), AVG(total_count) FROM ffpednadb WHERE left_gene = '%s' AND left_exon = '%s' AND left_exonid = '%d' AND right_gene = '%s' AND right_exon = '%s' AND right_exonid = '%d' "
    elif tumor == "cfdna":
        sql = "SELECT count(*) FROM cfdnadb WHERE left_gene = '%s' AND left_exon = '%s' AND left_exonid = '%d' AND right_gene = '%s' AND right_exon = '%s' AND right_exonid = '%d' "
        sql_avg = "SELECT AVG(unique_count), AVG(total_count) FROM cfdnadb WHERE left_gene = '%s' AND left_exon = '%s' AND left_exonid = '%d' AND right_gene = '%s' AND right_exon = '%s' AND right_exonid = '%d' "
    else:
        cursor.close()
        return [count,total_count,unique_count,count_hea,total_count_hea,unique_count_hea]

    cursor.execute(sql % data)
    count_tu = cursor.fetchall()     #条目数
    count_tumor = count_tu[0][0]

    cursor = connect.cursor()
    cursor.execute(sql_avg % data)
    avg_count = cursor.fetchall()
    for row in avg_count:
        unique_count_tumor = row[0]   #平均数
        total_count_tumor = row[1]
    cursor.close()
    result = [count,total_count,unique_count,count_hea,total_count_hea,unique_count_hea,count_tumor,total_count_tumor,unique_count_tumor]
    for i in range(len(result)):
        if result[i] == None:
            result[i] = 0
    return result
    
def factera_baseline(data):   #从数据库中获取基线
    cursor = connect.cursor()
    sql = "SELECT count(*) FROM factera_fulldb WHERE Region1 = '%s' AND Region2 = '%s'"
    sql_avg = "SELECT AVG(Break_support1), AVG(Break_support2), AVG(Total_depth) FROM factera_fulldb WHERE Region1 = '%s' AND Region2 = '%s'"
    cursor.execute(sql % data)
    count_tu = cursor.fetchall()     #条目数
    count = count_tu[0][0]

    cursor = connect.cursor()
    cursor.execute(sql_avg % data)
    avg_count = cursor.fetchall()
    for row in avg_count:
        left_count = row[0]   #平均数
        right_count = row[1]
        total_count = row[2]
    cursor.close()

    result = [count,left_count,right_count,total_count]
    for i in range(len(result)):
        if result[i] == None:
            result[i] = 0
    return result

def getjson():   #json是传入的文件，提取json中的信息
    global jsonfile
    logging.info(jsonfile)
    tmp = []
    try:
        setting = json.load(f, strict=False)
    except ValueError as e:
        print(e)
        return tmp
    try:
        num = re.compile(r'_(\d{4,6})[_|.]')   #正则匹配reads位置，返回一个列表
        hapnum = int(num.findall(jsonfile.split('/')[-1])[0])
    except AttributeError as e:
        return tmp
    except ValueError as e:
        return tmp
    except IndexError as e:
        return tmp
    if hapnum > 332767:
        return tmp

    try:
        time = (setting['time'])[:10]
    except KeyError as e:
        print(e)
        return tmp
    try:
        fusion = dict.keys(setting['fusions'])
    except KeyError as e:
        print(e)
        return tmp

    unique = []
    left_gene = []
    left_chr = []
    left_exon = []
    left_exonid = []
    left_strand = []
    left_position = []
    right_gene = []
    right_chr = []
    right_exon = []
    right_exonid = []
    right_strand = []
    right_position = []
    total = []
    a=0
    for i in fusion:
        if setting['fusions'][i]['left']['gene_name'] == setting['fusions'][i]['right']['gene_name']: #去除左右相同基因名
            continue 
        left_gene.append(setting['fusions'][i]['left']['gene_name'])
        left_chr.append(setting['fusions'][i]['left']['gene_chr'])
        left_exon.append(setting['fusions'][i]['left']['exon_or_intron'])
        left_exonid.append(setting['fusions'][i]['left']['exon_or_intron_id'])
        left_strand.append(setting['fusions'][i]['left']['strand'])
        right_gene.append(setting['fusions'][i]['right']['gene_name'])
        right_chr.append(setting['fusions'][i]['right']['gene_chr'])
        right_exon.append(setting['fusions'][i]['right']['exon_or_intron'])
        right_exonid.append(setting['fusions'][i]['right']['exon_or_intron_id'])
        right_strand.append(setting['fusions'][i]['right']['strand'])
        unique.append(setting['fusions'][i]['unique'])
        pattern = re.search(r'total: (\d+)',i,re.M|re.I)
        total.append(int(pattern.group(1)))
        reposition = re.compile(r'[X|\d+]:(\d+)')   #正则匹配reads位置，返回一个列表
        result = reposition.findall(i) 
        left_position.append(int(result[0]))
        right_position.append(int(result[1]))
    tmp = [hapnum,time,total,unique,left_gene,left_chr,left_exon,left_exonid,left_strand,right_gene,right_chr,right_exon,right_exonid,right_strand,left_position,right_position]
    return tmp

def list_all_files(directory):    #遍历目录及其子目录获取json
    _files = []
    listdb = os.listdir(directory) #列出文件夹下所有的目录与文件
    for i in range(0,len(listdb)):
        path = os.path.join(directory,listdb[i])
        if os.path.isdir(path):
            _files.extend(list_all_files(path))
        if os.path.isfile(path):
            _files.append(path)
            if path.endswith('.json'):
                jsondb.append(path)
    return jsondb

def dup(data):    #检测数据库中是否有相同ID的信息
    cursor = connect.cursor()
    sql = "SELECT id FROM fusiondb WHERE hapnum = '%d' "
    cursor.execute(sql % data)
    if cursor.rowcount != 0:
        cursor.close()
        fdb = 0
    else:
        fdb = 1
    cursor.close()
    cursor = connect.cursor()
    sql = "SELECT id FROM hapfusion WHERE hapnum = '%d' "
    cursor.execute(sql % data)
    if cursor.rowcount != 0:
        cursor.close()
        hdb = 0
    else:
        hdb = 1
    cursor.close()
    return [fdb,hdb]

def deldate(num):  
    cursor = connect.cursor()
    sql = "DELETE FROM fusiondb WHERE hapnum = '%d' "
    cursor.execute(sql % num)
    connect.commit()
    cursor.close()
    print("delete success")

def tocosf(data):   #上传数据到cosf数据表中
    cursor = connect.cursor()
    sql = "INSERT INTO hapfusion (hapnum,time,total_count,unique_count,left_gene,left_chr,left_exon,left_exonid,left_strand,right_gene,right_chr,right_exon,right_exonid,right_strand,left_position,right_position) VALUES  ('%d', '%s', '%d','%d', '%s','%s','%s','%d','%s', '%s','%s','%s','%d','%s','%d','%d')"
    cursor.execute(sql % data)
    connect.commit()
    cursor.close()

def comparecosf():  #从数据表中获取COSF信息
    cursor = connect.cursor()
    sql = "select * from cosf "
    cosfdata=[]
    cursor.execute(sql)
    rows = cursor.fetchall()
    for row in rows:
        cosfdata.append([row[1],row[2]])
    cursor.close()
    return cosfdata

def compareoncokb():
    cursor = connect.cursor()
    sql = "select * from oncokb "
    oncokbdata=[]
    cursor.execute(sql)
    rows = cursor.fetchall()
    for row in rows:
        oncokbdata.append([row[1],row[2]])
    cursor.close()
    return oncokbdata

def get_depth(haptotal):
    if args.bam:
        (status1, output1) = commands.getstatusoutput('samtools depth -r '+ haptotal[5] + ':' + str(haptotal[14]) + '-' + str(haptotal[14]) + ' ' + args.bam)  
        (status2, output2) = commands.getstatusoutput('samtools depth -r '+ haptotal[5] + ':' + str(haptotal[15]) + '-' + str(haptotal[15]) + ' ' + args.bam)
    else:
        output1 = output2 = '1\t10000000'
    if 'chr' not in output1 or 'chr' not in output2:
        output1 = output2 = '1\t10000000'
    if output1:    #计算丰度
        depth1 = output1.split('\t')[-1]

        try:
            left_vaf = float(haptotal[3]*200) / int(depth1)
        except ZeroDivisionError as e:
            print("ZeroDivisionError", e)
            left_vaf = 1
        if left_vaf > 100 or left_vaf < 0.1:
            left_vaf = 0
    else:
        left_vaf = 0
    if output2:
        depth2 = output2.split('\t')[-1]
        try:
            right_vaf = float(haptotal[3]*200) / int(depth2)
        except ZeroDivisionError as e:
            print("ZeroDivisionError", e)
            right_vaf = 1
        if right_vaf > 100 or right_vaf < 0.1:
            right_vaf = 0
    else:
        right_vaf = 0
 
    return [round(left_vaf,2),round(right_vaf,2)]

def get_factera_data(data):
    cursor = connect.cursor()
    sql = "SELECT * FROM gene_hg19 WHERE chr = '%s' AND position1 < '%d' and position2 > '%d'"
    cursor.execute(sql % data)
    count_tu = cursor.fetchall()
    if count_tu == ():
        cursor = connect.cursor()
        sql = "select * from gene_hg19 where chr = '%s' and position1 > '%d' order by (position1-'%d') limit 1;"
        cursor.execute(sql % data)
        gene1 = cursor.fetchall()

        cursor = connect.cursor()
        sql = "select * from gene_hg19 where chr = '%s' and position2 < '%d' order by abs(position2-'%d') limit 1;"
        cursor.execute(sql % data)
        gene2 = cursor.fetchall()
        return [gene1,gene2]
    else:
        return count_tu

def get_exon(data,position):
    exon = []
    for row in data:
        gene = row.split('name')[1].split(';')[0][1:]
        nm = row.split('ID')[1].split(';')[0][1:]
        plus = row.split('\t')[6]
        exon_num = len(row.split('CDS'))-1
        if exon_num == 0:
            exon.append([gene,nm,plus,exon_num])
            continue
        num = 0
        for i in row.split('CDS'):
            if i.split('\t')[1] == 'UCSC':
                continue
            if position < int(i.split('\t')[1]):
                 if plus == '+':
                    exon.append([gene,nm,plus,num])
                 else:
                    exon.append([gene,nm,plus,exon_num-num])
                 break
            if position > int(i.split('\t')[2]):
                 num += 1
    return exon

def get_result(data,position):
    if () in data:
        return [1,1]
    if isinstance(data,list):
        gene1 = data[0][0][4].split('name')[1].split(';')[0][1:]
        gene2 = data[1][0][4].split('name')[1].split(';')[0][1:]
        return [gene1,gene2]
    else:
        infomation = []
        for i in data:
            if len(i) < 2:
                break
            infomation.append(i[4])
        exon = get_exon(infomation,position)
        return exon

def fusionscan_writein(total):
    start = time.clock()
    write_inlist = []
    for i in range(len(total[3])):
        haptotal = (total[0],total[1],total[2][i],total[3][i],total[4][i],total[5][i],total[6][i],total[7][i],total[8][i],total[9][i],total[10][i],total[11][i],total[12][i],total[13][i],total[14][i],total[15][i])
        #hapnum,time,total,unique,left_gene,left_chr,left_exon,left_exonid,left_strand,right_gene,right_chr,right_exon,right_exonid,right_strand,left_position,right_position
        getdata = (total[4][i],total[6][i],total[7][i],total[9][i],total[11][i],total[12][i])
        if "heal" in jsonfile:
            ut = getbaseline(getdata,"heathdna")
        elif "cf" in jsonfile:
            ut = getbaseline(getdata,"cfdna")
        else:
            ut = getbaseline(getdata,"ffpe")
        special = 0
        hot_gene = ['ALK','ROS1','RET','PRKACA']
        for j in hot_gene:
            if j in total[4][i].split('_')[0] or j in total[9][i].split('_')[0]:
                special += 1
        if [(total[4][i]).split('_')[0],(total[9][i]).split('_')[0]] in cosfdata or [(total[9][i]).split('_')[0],(total[4][i]).split('_')[0]] in cosfdata:
            special += 1
        if [(total[4][i]).split('_')[0],(total[9][i]).split('_')[0]] in oncokbdata or [(total[9][i]).split('_')[0],(total[4][i]).split('_')[0]] in oncokbdata:
            special += 1
        if booldb[1] == 1 and special > 0:
            tocosf(haptotal)
        if (total[4][i]).split('_')[0] != (total[9][i]).split('_')[0] and total[3][i] > 5:
            special +=1
        if ut[0] > 200:
            special -=1
        depth = [0,0]
        if special > 0:
            depth = get_depth(haptotal)
        left_depth = depth[0]
        right_depth = depth[1]
        
        tuptotal = (total[0],total[1],total[2][i],total[3][i],total[4][i],total[5][i],total[6][i],total[7][i],total[8][i],total[9][i],total[10][i],total[11][i],total[12][i],total[13][i],total[14][i],total[15][i],special)
        #if booldb[0] == 1: # (写入数据库)
            #todb(tuptotal)
            #logging.info("update data")

        setin = [total[0],special,(total[4][i]).split('_')[0],(total[9][i]).split('_')[0],str(left_depth),str(right_depth),total[2][i],total[3][i],total[6][i],total[7][i],total[11][i],total[12][i],total[8][i]] + ut
        #hapnum,special,left_gene,right_gene,'left_vaf','right_vaf',total,unique,total_count,total_average,unique_average,left_exon,left_exonid,right_exon,right_exonid,strand
        write_inlist.append(setin)
        
        del(tuptotal)
        del(getdata)
    for i in write_inlist:
        csv_write.writerow(i)
    elapsed = (time.clock() - start)
    print("fusionscan_writein Time used:",elapsed)

def fusion_writein(data):
    start = time.clock()
    write_inlist = []
    csv_inlist = [['id','level','Region1','Region2','left_vaf','right_vaf','left_support','right_support','left_exon','left_exonid','right_exon','right_exonid','strand','count','left_avg_support','right_avg_support','avg_depth']]
    path = (args.file[::-1].split('/',2)[-1])[::-1] + '/fusion'
    print(path)
    for i in os.listdir(path):
        if 'important' not in i:
            with open(path + '/' + i,'r') as out:
                a=0
                for line in out:
                    if a == 0:
                        a+=1
                        continue
                    line = line.split('\t')
                    write_inlist.append([data]+line[:-3]+[line[-2]])

    for i in write_inlist:
        special = 0
        hot_gene = ['ALK','ROS1','RET','NRG1','NTRK1','PAX7','FGFR2','KMT2A','ETV6','NAB2','DDIT3','NUTM1','FUS','RARA','SS18','PRKACA','PAX8','PAX3','TMPRSS2','EWSR1','PDGFB','FGFR3','BRAF','MYC','NR4A3','NTRK2','TFE3','PRKACA']
        for j in hot_gene:
            if j in i[2] or j in i[3]:
                special += 1
        if [i[2],i[3]] in cosfdata or [i[3],i[2]] in cosfdata:
            special += 1
        if [i[2],i[3]] in cosfdata or [i[3],i[2]] in oncokbdata:
            special += 1
        if i[2] != i[3] and (float(i[6])/float(i[17]) > 0.02 or float(i[7])/float(i[17]) > 0.02):
            special +=1
        getdata = (i[2],i[3])
        ut = factera_baseline(getdata)
        if ut[0] >= 200:
            special -=1
        setin = ['fusion',special,i[2],i[3],float(i[6])/float(i[17]),float(i[7])/float(i[17]),i[6],i[7],'intron','','intron','',''] + ut
        chr1,position1 = i[4].split(':')[0],int(i[4].split(':')[1])
        chr2,position2 = i[5].split(':')[0],int(i[5].split(':')[1])
        data1 = (chr1,position1,position1)
        data2 = (chr2,position2,position2)
        data1 = get_factera_data(data1)
        data2 = get_factera_data(data2)
        data1 = get_result(data1,position1)
        data2 = get_result(data2,position2)
        if data1 == []:
            data1 = [1,1]
        if data2 == []:
            data2 = [1,1]
        if isinstance(data1[0],unicode):
            setin[9] = data1[0] + '-' + data1[1]
        elif data1 == [1,1]:
            setin[9] = 'not found'
        else:
            for gene in data1:
                setin[9] = gene[0] + '-' + str(gene[-1])
                    
        if isinstance(data2[0],unicode):
            setin[11] = data2[0] + '-' + data2[1]
        elif data2 == [1,1]:
            setin[11] = 'not found'
        else:
            for gene in data2:
                setin[11] = gene[0] + '-' + str(gene[-1])

        csv_inlist.append(setin)
    for i in csv_inlist:
        csv_write.writerow(i)
    cursor = connect.cursor()
    sql = "SELECT id FROM factera_fulldb WHERE hapnum = '%d' "
    cursor.execute(sql % data)
    if cursor.rowcount != 0:
        cursor.close()
        fdb = 0
    else:
        fdb = 1

    if fdb == 1:
        for i in write_inlist:
            in_data = tuple(i)
            cursor = connect.cursor()
            sql = "INSERT INTO factera_fulldb (hapnum,Est_Type,Region1,Region2,Break1,Break2,Break_support1,Break_support2,Break_offset,Orientation,Order1,Order2,Break_depth,Proper_pair_support,Unmapped_support,Improper_pair_support,Paired_end_depth,Total_depth,Fusion_seq,Non_templated_seq) VALUES  ('%d', '%s', '%s','%s','%s', '%s', '%s', '%s','%s','%s', '%s','%s', '%s','%s','%s', '%s','%s', '%s','%s', '%s')"
            cursor.execute(sql % in_data)
            connect.commit()
            cursor.close()

    elapsed = (time.clock() - start)
    print("fusion_writein Time used:",elapsed)

def check_json(jsonfile):
    jsonfile1 = open(jsonfile)
    data = jsonfile1.read()
    jsonfile1.close()
    qual = '"qual":"(.*?)"\n'
    qual_split = re.compile(qual).findall(data)
    if '"' not in str(qual_split):
        return
    else:
        for i in qual_split:
            if '"' in i:
                data = data.split(i)[0] + i.replace('\"','\'') + data.split(i)[1]
        with open(jsonfile,'w') as w:
            w.write(data)

def main():
    start = time.clock()
    print(jsonfile)
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(filename)s %(levelname)s %(message)s',
                        datefmt='%a, %d %b %Y %H:%M:%S', filename=log, filemode='w')
    total=getjson()
    if total == 0:
        return 
    if len(total) > 0:
        global booldb,cosfdata,oncokbdata,csv_write
        booldb = dup(total[0])
        cosfdata = comparecosf()
        oncokbdata = compareoncokb()
        if total[3] == []:  #判断json中是否存在信息
            logging.info("fusion don't exist")
            return 0 
        csv_write = csv.writer(w,dialect='excel')
        if "cf" in jsonfile:
            title = ['hapnum','level','left_gene','right_gene','left_vaf','right_vaf','total','unique','left_exon','left_exonid','right_exon','right_exonid','strand','baseline_count','total_average','unique_average','health_baseline','health_total','health_unique','blood_baseline','blood_total','blood_unique']
        elif "cf" not in jsonfile:
            title = ['hapnum','level','left_gene','right_gene','left_vaf','right_vaf','total','unique','left_exon','left_exonid','right_exon','right_exonid','strand','baseline_count','total_average','unique_average','health_baseline','health_total','health_unique','tumor_baseline','tumor_total','tumor_unique']
        csv_write.writerow(title)
        fusionscan_writein(total)
        fusion_writein(total[0])
        f.close()
        w.close()
    elapsed = (time.clock() - start)
    print("main Time used:",elapsed)


start = time.clock()
connect = pymysql.connect(host="192.168.1.11",user="root",passwd="haplox2017",db="mutation",port=3306,charset="utf8") #host为服务器地址，db为数据库名，port为接口
log = "fusion.log"
if args.delete!=None:
    delnum = int(args.delete)
    if type(args.delete) != int:
        print('please input a int')
        exit(0)
    deldate(delnum)
if args.file!=None:
    jsonfile=args.file
    check_json(jsonfile)
    if jsonfile.endswith('.json') == 0:
        print('please input a json file')
        exit(0)
    if not os.path.exists(jsonfile):
        os.system("mv {S2} {S1}".format(S2 = jsonfile[::-1].split('/',1)[1][::-1] + "/*.json" ,S1 = jsonfile))
    f = io.open(jsonfile,'r',encoding='utf-8')
    w = io.open((jsonfile[:-5]+'.csv') , 'wb')
    main()
if args.directory!=None:
    jsondb=[]
    jsonfile_list = list_all_files(args.directory)
    if jsondb==[]:
        print("this directory is empty")
        exit(1)
    for jsonfile in jsonfile_list:
        try:
            f = io.open(jsonfile,'r',encoding='utf-8')
        except IOError as e:
            continue
        try:
            w = io.open((jsonfile[:-5]+'.csv') , 'wb')
        except IOError as e:
            print("IOError", e)
            continue
        else:
            main()
connect.close() 
elapsed = (time.clock() - start)
print("total Time used:",elapsed)



