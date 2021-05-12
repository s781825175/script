import os

path = os.getcwd()
sample = path.split('/')[-1]
print(path,sample)
os.mkdir(path + '/cnv_single')

os.system('/haplox/users/huang/myGit/cnvkit/rawPython/bin/python /haplox/users/huang/myGit/cnvkit/cnv.py -s {S1}/{S2}_rg.bam -b /haplox/ref/bed/605panel_v2_primary_targets_gene.bed -r /haplox/ref/GATK/ucsc.hg19/ucsc.hg19.fasta -o {S1}/cnv_single'.format(S1=path,S2=sample))

os.system('Rscript /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/getCnv.R /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/cnv/cnv.csv  {S1}/cnv_single/{S2}_rg_cnv_result.txt {S1}/cnv_single/{S2}_rg_cnv_result.csv'.format(S1=path,S2=sample))

#curl haplab.haplox.net/api/report/csv?type=cnv -F "import_file=@ {file} ")
