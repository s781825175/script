import os



for j in range(1,13):
    if j<10:
        j = '0'+str(j) 
    os.system('ossutil ls oss://sz-hapres/haplox/hapyun/2019{S1}/pair_IDT_HPCH_ -a | xargs -n 1 -P 15 ossutil rm -f'.format(S1=j))

for j in range(1,13):
    if j<10:
        j = '0'+str(j)
    os.system('ossutil ls oss://sz-hapres/haplox/hapyun/2019{S1}/pair_wes_IDT_no_facter_ -a | xargs -n 1 -P 15 ossutil rm -f'.format(S1=j))
