import xlrd
from multiprocessing import Pool
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("excel", type=str, help="minimum number of breakpoint-spanning reads required for output")
args = parser.parse_args()

opts = vars(args)


def long_time_task(data):
    if data[2] != 'vs':
        #os.mkdir('/haplox/rawout/{S1}/{S2}/fupy_{S3}'.format(S1=data[0],S2=data[1],S3=data[1].split('_')[-1]))
        os.system('python3 /haplox/users/yangbo/factrea/py/pyfusion.py /haplox/rawout/{S1}/{S2}/{S2}_rg.bam -o /haplox/rawout/{S1}/{S2}/fupy_{S3}'.format(S1=data[0],S2=data[1],S3=data[1].split('_')[-1]))
        print('python3 /haplox/users/yangbo/factrea/py/pyfusion.py /haplox/rawout/{S1}/{S2}/{S2}_rg.bam -o /haplox/rawout/{S1}/{S2}/fupy_{S3}'.format(S1=data[0],S2=data[1],S3=data[1].split('_')[-1]))

if __name__=='__main__':
    xlsfile = opts['excel']
    book = xlrd.open_workbook(xlsfile)
    sheet0 = book.sheet_by_index(0)
    sheet_name = book.sheet_names()[0]
    sheet1 = book.sheet_by_name(sheet_name)
    nrows = sheet0.nrows

    p = Pool(4)

    for row in range(nrows):   #与上机信息表匹配
        data = sheet0.row_values(row)
        print(data)
        p.apply_async(long_time_task, args=(data,))

    p.close()
    p.join()



