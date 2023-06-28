#!/usr/bin/python3
'''This python script takes one gwas summary statistics file from
mibiogen (microbiome GWAS) and one from pan-uk biobank
as input, takes intersection, recalculates pukbb pvals, and outputs
needed data for TwoSampleMR
##Camilla's edit
'''
import gzip
import re
import sys
import argparse
import os
import math

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Script to retrieve rsids')
    parser.add_argument('-m', '--micro',
                        help='path to microbiome GWAS summary statistics file',
                        required='True'
                        )
    parser.add_argument('-p', '--panukb',
                        help='path to pan-UKBiobank GWAS summary statistics file',
                        required='True'
                        )
    parser.add_argument('-o', '--output',
                        help='output file prefix',
                        required='True'
                        )
    return parser.parse_args(args)

#retrieve command line arguments

args = check_arg(sys.argv[1:])
micro_infile = args.micro
pan_infile = args.panukb
outfile = args.output


#build dictionary from microbiome file with
# c:pos as keys and needed data as list of values (including rsid)

microdict = dict()

for line in gzip.open(micro_infile):
#for line in gzip.open('/home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/33462485-GCST90017070-EFO_0007874-Build37.f.tsv.gz'):
    arr = line.strip().split()
    chr, pos, rsid, other_allele, effect_allele, beta, se, pval = arr[0:8]
    cpos = chr + b":" + pos
    microdict[cpos] = [rsid, chr, pos, effect_allele, other_allele, beta, se, pval]

#out = gzip.open('/home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/intersectionfile.txt.gz', 'wb')
out = gzip.open(outfile + ".txt.gz","wb")
out.write(b"SNP\tCHR\tPOS\tmicro_effect\tmicro_other\tmicro_beta\tmicro_se\tmicro_p\tpanukb_effect\tpanukb_other\tpanukb_beta\tpanukb_se\tpanukb_p\n")

#with gzip.open('/home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/phecode-411.4-both_sexes.tsv.bgz') as f:
with gzip.open(pan_infile) as f:
    for line in f:
        arr = line.strip().split()
        #use hq (high-quality) meta results
        chr, pos, ref_other, alt_effect, af, beta, se, lnpval = arr[0:8]
        cpos = chr + b":" + pos
        #retrieve microbiome data
        if cpos in microdict:
            microdata = microdict[cpos]
            #retrieve microbiome pval for conditional test
            microp = (float(microdata[7]))
            #retrieve rsid for conditional test
            rsid = microdata[0]
            if rsid != b"NA" and lnpval != b"NA": #need rsid for twosampleMR, some microbiome SNPs don't have rsids
                #convert pan-ukb p-val from ln(pval) to pval
                pval = math.e**(float(lnpval))
                #convert pan-ukb p-val from neglog10pval to pval
                #pval = 10**(-(float(lnpval)))
                #only write SNPs with P< 0.0001 in at least 1 GWAS
                if microp < 0.0001 or pval < 0.0001:
                    #encode float as bytes
                    pval_bytestring = str(pval).encode('utf-8')
                    out.write(b"\t".join(microdata) + b"\t" + alt_effect + b"\t" + ref_other + b"\t" + beta + b"\t" + se + b"\t" + pval_bytestring + b"\n")

out.close()

