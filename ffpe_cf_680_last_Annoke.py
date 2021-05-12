#coding=utf-8
import csv
import xlrd
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", nargs='?', action="store" , type=str, help="json file name")
parser.add_argument("-d", "--directory", nargs='?', action="store" , type=str, help="directory and subdirectory ")
args = parser.parse_args()

def make_shell(ffpe_data_dir,ffpe_data_id,cf_data_dir,cf_data_id,g_data_dir,g_data_id,data_id):
    shell = '''
#!/bin/bash
Rscript /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/getCnv_ff_cf.R /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/cnv/{S2}_rg_cnv.genes188in680.csv /haplox/rawout/{S3}/{S4}/cnv/{S4}_rg_cnv.genes188in680.csv /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/cnv/{S2}_rg_cnv_tumor_cfdna.genes188in680.csv 
#-----------------------------------------------------------------------------------------
Rscript /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/getFFPE_cfDNA_snv_indel_v2.R /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/ffpe_vs_pbl/genes528in680.bed /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/{S8}_{S9}.indel-nobias-GB18030-baseline.csv /haplox/rawout/{S3}/{S4}/{S10}_{S6}.indel-nobias-GB18030-baseline.csv /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/mrbam_{S8}_{S9}_indel_combine_genes528in680.csv 0.1 
Rscript /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/getFFPE_cfDNA_snv_indel_v2.R /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/ffpe_vs_pbl/genes528in680.bed /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/{S8}_{S9}.snv-nobias-GB18030-baseline.csv /haplox/rawout/{S3}/{S4}/{S10}_{S6}.snv-nobias-GB18030-baseline.csv /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/mrbam_{S8}_{S9}_snv_combine_genes528in680.csv 0.1 
python /thinker/net/tools/Annokb.py -f /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/mrbam_{S8}_{S9}_indel_combine_genes528in680.csv  -t indel -o /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/Annokb_mrbam_{S8}_{S9}_indel_combine_genes528in680.csv 
python /thinker/net/tools/Annokb.py -f /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/mrbam_{S8}_{S9}_snv_combine_genes528in680.csv  -t snv -o /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/Annokb_mrbam_{S8}_{S9}_snv_combine_genes528in680.csv 
#-----------------------------------------------------------------------------------------
python /haplox/users/yangbo/futionbase.py -f /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/fusionscan/{S2}_fusion.json -b /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/{S2}_rg.bam 
python /haplox/users/yangbo/futionbase.py -f /haplox/rawout/{S3}/{S4}/fusionscan/{S4}_fusion.json -b /haplox/rawout/{S3}/{S4}/{S4}_rg.bam 
python3 /haplox/users/yangbo/Annokb_fusion.py -f /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/fusionscan/{S2}_fusion.csv -c /haplox/rawout/{S3}/{S4}/fusionscan/{S4}_fusion.csv -o /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/fusionscan/{S2}_ffpe_vs_cfdna_fusion.csv 
#-----------------------------------------------------------------------------------------
perl /haplox/users/ZhouYQ/germline/bin/Gene2Disease.pl /haplox/users/ZhouYQ/germline/bin/Database/female_cancer.list /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result/{S5}.germline.txt /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result/{S5}.cancer.female.txt /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result/{S5}.nocancer.female.txt /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result 
perl /haplox/users/ZhouYQ/germline/bin/Gene2Disease.pl /haplox/users/ZhouYQ/germline/bin/Database/male_cancer.list /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result/{S5}.germline.txt /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result/{S5}.cancer.male.txt /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result/{S5}.nocancer.male.txt /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result 
perl /haplox/users/liaowt/Script/germline/germline_trans_v2.pl /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result/{S5}.cancer.female.txt /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result/{S5}_trans.cancer.female.txt 
perl /haplox/users/liaowt/Script/germline/germline_trans_v2.pl /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result/{S5}.cancer.male.txt /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result/{S5}_trans.cancer.male.txt 
python /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/zhaoys/txt_xls.py -i /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result/{S5}_trans.cancer.female.txt -o /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result/{S5}_trans.cancer.female.xls 
python /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/zhaoys/txt_xls.py -i /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result/{S5}_trans.cancer.male.txt -o /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result/{S5}_trans.cancer.male.xls 
python /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/zhaoys/txt_xls.py -i /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result/{S5}.information.txt -o /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result/{S5}.information.xls 
python /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/zhaoys/txt_xls.py -i /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result/{S5}.chem_451.txt -o /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result/{S5}.chem_451.xls 
python /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/zhaoys/txt_xls.py -i /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result/{S5}.Target_451.txt -o /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result/{S5}.Target_451.xls 
curl haplab.haplox.net/api/report/chemotherapy/{S9} -F "qc=@/haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result/{S5}.information.txt" 
curl haplab.haplox.net/api/report/chemotherapy/{S9} -F "chem=@/haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result/{S5}.chem_451.txt" 
curl haplab.haplox.net/api/report/chemotherapy/{S9} -F "germline=@/haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result/{S5}_trans.germline_680.cancer.txt" 
curl haplab.haplox.net/api/report/hla-v?data_id={S9} -F "import_file=@/haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/HLA/{S5}_sortbam-hla.csv" 
curl haplab.haplox.net/api/report/hrr?data_id={S9} -F "import_file=@/haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/germline/result/{S5}_trans.germline_DDR.cancer.txt" 
curl haplab.haplox.net/api/report/cfdna/{S9}?cfdna_data_id={S6}
curl haplab.haplox.net/api/report/gdna/{S9}?gdna_data_id={S7}
Rscript /haplox/users/liaowt/Script/somatic/result_combine/last_pair_mutscan_fusion_virus_sheet_680_combine.R /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/ /haplox/rawout/{S3}/{S4}/ 
'''.format(S1=ffpe_data_dir,S2=ffpe_data_id,S3=cf_data_dir,S4=cf_data_id,S5=g_data_id,S6=cf_data_id.split('_')[-1],S7=g_data_id.split('_')[-1],S8=ffpe_data_id.split('_')[1],S9=data_id,S10=cf_data_id.split('_')[1])
    with open('last_tumor_cfdna_gdna_huang_680.sh','w') as f:
        f.write(shell)

path = os.getcwd()
for i in os.listdir('/x01_haplox/hapreports/0.上机信息汇总'):
    if '_share.xlsx' in i:
        xlsfile = '/x01_haplox/hapreports/0.上机信息汇总/' + i
book = xlrd.open_workbook(xlsfile)
sheet0 = book.sheet_by_index(0)
sheet_name = book.sheet_names()[0]
sheet1 = book.sheet_by_name(sheet_name)
nrows = sheet0.nrows
name = {}
data_id = int(path.split('_')[-1])
ffpe_data_dir = path.split('/')[3]
ffpe_data_id = path.split('/')[-1]
tmp1 = tmp2 = data_id
for row in range(nrows):
    if sheet0.cell_value(row, 3) != '':
        if sheet0.cell_value(row, 3) not in name:
            name[sheet0.cell_value(row, 3)] = [sheet0.row_values(row)]
        else:
            name[sheet0.cell_value(row, 3)].append(sheet0.row_values(row))
    if path.split('/')[-1] == sheet0.cell_value(row, 1):
        sample = sheet0.cell_value(row, 3)
for i in name[sample]:
    if 'cfdna' in i[1]:
        if tmp1 > abs(int(i[1].split('_')[-1]) - data_id):
            cf_data_dir = i[0]
            cf_data_id = i[1]
            tmp1 = abs(int(i[1].split('_')[-1]) - data_id)
    if 'gdna' in i[1]:
        if tmp2 > abs(int(i[1].split('_')[-1]) - data_id):
            g_data_dir = i[0]
            g_data_id = i[1]
            tmp2 = abs(int(i[1].split('_')[-1]) - data_id)

if args.file:
    cf_data_id = args.file
    cf_data_dir = args.directory

make_shell(ffpe_data_dir,ffpe_data_id,cf_data_dir,cf_data_id,g_data_dir,g_data_id,str(data_id))

