read -p "print input project name:" S1
ossutil ls oss://sz-hapres/haplox/hapyun/202007/ -s | grep ".*${S1}.*.fq.gz$" | xargs -n 1 -P 15 ossutil rm -f
ossutil ls oss://sz-hapres/haplox/hapyun/202007/ -s | grep ".*${S1}.*.bam$" | xargs -n 1 -P 15 ossutil rm -f
ossutil ls oss://sz-hapres/haplox/hapyun/202007/ -s | grep ".*${S1}.*.depth$" | xargs -n 1 -P 15 ossutil rm -f
ossutil ls oss://sz-hapres/haplox/hapyun/202007/ -s | grep ".*${S1}.*.pileup$" | xargs -n 1 -P 15 ossutil rm -f
ossutil ls oss://sz-hapres/haplox/hapyun/202007/ -s | grep ".*${S1}.*.bai$" | xargs -n 1 -P 15 ossutil rm -f
ossutil ls oss://sz-hapres/haplox/hapyun/202007/ -s | grep ".*${S1}.*.avinputs$" | xargs -n 1 -P 15 ossutil rm -f
ossutil ls oss://sz-hapres/haplox/hapyun/202007/ -s | grep ".*${S1}.*.bam_temp_len.png$" | xargs -n 1 -P 15 ossutil rm -f
ossutil ls oss://sz-hapres/haplox/hapyun/202007/ -s | grep ".*${S1}.*.rfq" | xargs -n 1 -P 15 ossutil rm -f
ossutil ls oss://sz-hapres/haplox/hapyun/202007/ -s | grep ".*${S1}.*.fastq.gz" | xargs -n 1 -P 15 ossutil rm -f


