#!/usr/bin/env python
import sys, os, pysam, re
import logging

if len(sys.argv) <= 2:
    print("Usage: ubam bamfile 1:1000-2000|bed_file <output_file>\n")
    exit(-1)

if not os.path.exists(sys.argv[1]):
    print("%s not exist\n" % sys.argv[1])
    exit(-1)

mode = re.match(r'(\w+):(\d+)-(\d+)',sys.argv[2])
if mode:
    chrom = mode.group(1)
    start = mode.group(2)
    end = mode.group(3)

    if len(sys.argv) == 3:
        suffix = os.path.basename(sys.argv[1])
        output_file = './' + sys.argv[2] + '_' + suffix
    else:
        output_file = sys.argv[3]

    bam = pysam.AlignmentFile(sys.argv[1],'rb')
    reads = bam.fetch('chr'+ chrom, max(1,int(start)-200), int(end) + 200)
    output_file = pysam.AlignmentFile(output_file, 'wb', template = bam)
    count = 0
    for read in reads:
        count += 1
        output_file.write(read)
    
    print("%d reads extracted\n " % count)
    output_file.close()
    bam.close()

else:
    exit(-1)
