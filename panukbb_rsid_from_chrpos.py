#!/usr/bin/python3
'''This python script 1000g plink .bim file for a dictionary to map rsids for TwoSampleMR
##Camilla's edit
'''

import gzip
import re
import sys
import argparse
import os
import math

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Script to find intersection')
    parser.add_argument('-b', '--bim',
                        help='path to 1000g .bim file',
                        required='True'
                        )
    parser.add_argument('-f', '--first',
                        help='path to summary statistics file',
                        required='True'
                        )
    parser.add_argument('-s', '--second',
                        help='path to summary statistics file',
                        required='True'
                        )
    parser.add_argument('-o', '--output',
                        help='output file prefix',
                        required='True'
                        )
    return parser.parse_args(args)

#retrieve command line arguments

args = check_arg(sys.argv[1:])
bim_file = args.bim
first_infile = args.first
second_infile = args.second
outfile = args.output


#build dictionary from microbiome file with
# c:pos as keys and needed data as list of values (including rsid)

bimdict = dict()
firstdict = dict()
#count = 0

for line in open(bim_file):
#for line in open('/home/camilla/all_phase3.rsid.maf001.dups2.removed.bim'):
    arr = line.strip().encode('utf-8').split()
    chr, rsid, posincm, pos, allele1, allele2 = arr[0:6]
    cpos = chr + b":" + pos
    bimdict[cpos] = rsid

for line in gzip.open(first_infile):
#for line in gzip.open("continuous-1160-both_sexes.tsv.bgz"):
    arr = line.strip().split()
    chr, pos, ref_other, alt_effect, af, beta, se, lnpval = arr[0:8]
    cpos = chr + b":" + pos
    if cpos in bimdict:
      snp_id = bimdict[cpos]
      firstdict[cpos] = [snp_id, chr, pos, ref_other, alt_effect, beta, se, lnpval]
    #count = count + 1
    #if count > 1 and lnpval != b"NA":
    #  firstp = math.e**(float(lnpval))
    
outfile = "/home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/intersection"
out = gzip.open(outfile + ".txt.gz","wb")
out.write(b"SNP\tCHR\tPOS\tfirst_other\tfirst_effect\tfirst_beta\tfirst_se\tfirst_lnp\tfirst_p\tsecond_effect\tsecond_other\tsecond_beta\tsecond_se\tsecond_p\n")

with gzip.open(second_infile) as f:
#with gzip.open('/home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/phecode-411.4-both_sexes.tsv.bgz') as f:
    for line in f:
        arr = line.strip().split()
        #use hq (high-quality) meta results
        chr, pos, ref_other, alt_effect, afcases, afcontrols, beta, se, lnpval = arr[0:9]
        cpos = chr + b":" + pos
        #retrieve microbiome data
        if cpos in firstdict:
            firstdata = firstdict[cpos]
            #retrieve microbiome pval for conditional test
            if firstdata[7] != b"NA" and firstdata[7] != b"pval_meta_hq":
              firstp = math.e**(float(firstdata[7]))
              if lnpval != b"NA" and lnpval != b"neglog10_pval_meta_hq": 
              #convert pan-ukb p-val from ln(pval) to pval
              #pval = math.e**(float(lnpval))
              #convert pan-ukb p-val from neglog10pval to pval???????
                pval = 10**(-(float(lnpval)))
                #only write SNPs with P< 0.0001 in at least 1 GWAS
                if firstp < 0.0001 or pval < 0.0001:
                #encode float as bytes
                  first_pval_bytestring = str(firstp).encode('utf-8')
                  pval_bytestring = str(pval).encode('utf-8')
                  out.write(b"\t".join(firstdata) + b"\t" + first_pval_bytestring + b"\t" + alt_effect + b"\t" + ref_other + b"\t" + beta + b"\t" + se + b"\t" + pval_bytestring + b"\n")

out.close()

