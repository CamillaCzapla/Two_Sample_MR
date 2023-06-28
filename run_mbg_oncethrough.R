setwd('/home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics')
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
#suppressMessages(library(argparse))
suppressMessages(library(TwoSampleMR))
suppressMessages(library(ieugwasr))
'%&%' = function(a,b) paste (a,b,sep='')

##plotting - only change file names and phenotype names (and column names inside python code if needed)
output_file_path = "/home/camilla/TwoSampleMR/MR"
File.Name = "33462485-GCST90017115-EFO_0007874-Build37.f.tsv.gz"
File.Name2 = "continuous-1160-both_sexes.tsv.bgz"
#for file naming - can be personalized differently
sumstats2 = "continuous-1160-both_sexes"
phenotype1 = "Gut microbiota abundance (phylum Lentisphaerae id.2238)"
phenotype2 = "sleep duration"
#for naming intersection file - not as important when plotting one test at a time
sumstats = "mbfile"
my_command_1 = "python3 /home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/00_mibiogen_panukbb_4_twosampleMR_heart.py --micro /home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/" %&% File.Name %&%" --panukb /home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/" %&% File.Name2 %&% " --output /home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/" %&% sumstats %&% "_intersection"
#call python script with system() in order to have a smaller file for read_outcome_data()
#otherwise, a Cstack usage error
system(my_command_1)
filename_string = sumstats %&% "_intersection.txt.gz"
microbiome_exp_dat <- read_exposure_data(
  filename = filename_string,
  #clump = TRUE,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "micro_beta",
  se_col = "micro_se",
  effect_allele_col = "micro_effect",
  other_allele_col = "micro_other",
  pval_col = "micro_p",
  #min_pval = 1e-04,
  chr_col = "CHR",
  pos_col = "POS"
  )
                 
microbiome_exp_dat$exposure <- phenotype2
#filtering microbiome data to pval <= 1e-04 (pvalue in MiBioGen paper was 1e-05)
microbiome_exp_dat <- filter(microbiome_exp_dat, pval.exposure <= 1e-05)
#microbiome_exp_dat <-filter(microbiome_exp_dat, pval.exposure <= 5e-08)
chd_out_dat <- read_outcome_data(
       #snps = microbiome_exp_dat$SNP,
       filename = filename_string,
       #filename = "intersection.tsv",
       #snps = microbiome_exp_dat$SNP,
       sep = "\t",
       #phenotype_col = "Phenotype",
       snp_col = "SNP",
       beta_col = "panukb_beta",
       se_col = "panukb_se",
       #eaf_col = "eaf",
       effect_allele_col = "panukb_effect",
       other_allele_col = "panukb_other",
       pval_col = "panukb_p",
       #units_col = "units",
       #ncase_col = "ncase",
       chr_col = "CHR",
       pos_col = "POS"
       )
                                  
chd_out_dat$outcome <- phenotype1
##HARMONISE DATA (Intersection/CREATES NEW COMBINED DATA FRAME WITH EXPOSURE AND OUTCOME DATA)
dat_for <- harmonise_data(microbiome_exp_dat, chd_out_dat)

#clump
if (dim(dat_for)[1] != 0) {
  dat_for <- dplyr::rename(dat_for, rsid = SNP)
  dat_for <- dplyr::rename(dat_for, pval = pval.exposure)
  
  #hacked ld_clump_local function with --allow-extra-chr flag
  #skip if not necessary?
  ld_clump_local_hacked <- function(dat, clump_kb, clump_r2, clump_p, bfile, plink_bin, allow_extra_chr)
  {
    
    # Make textfile
    shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
    fn <- tempfile()
    write.table(data.frame(SNP=dat[["rsid"]], P=dat[["pval"]]), file=fn, row.names=F, col.names=T, quote=F)
    
    fun2 <- paste0(
      shQuote(plink_bin, type=shell),
      " --bfile ", shQuote(bfile, type=shell),
      " --keep /home/wheelerlab3/bioi397/homeworks/HW2/1000g_plink2/EUR_subset",
      " --clump ", shQuote(fn, type=shell), 
      " --clump-p1 ", clump_p, 
      " --clump-r2 ", clump_r2, 
      " --clump-kb ", clump_kb, 
      " --allow-extra-chr", #allow_extra_chr,
      " --out ", shQuote(fn, type=shell)
    )
    system(fun2)
    res <- read.table(paste(fn, ".clumped", sep=""), header=T)
    unlink(paste(fn, "*", sep=""))
    y <- subset(dat, !dat[["rsid"]] %in% res[["SNP"]])
    if(nrow(y) > 0)
    {
      message("Removing ", length(y[["rsid"]]), " of ", nrow(dat), " variants due to LD with other variants or absence from LD reference panel")
    }
    return(subset(dat, dat[["rsid"]] %in% res[["SNP"]]))
  }
  
  ##LD clumping locally
  ##Error HERE
  dat_for <- ld_clump_local_hacked(dat_for, clump_kb = 10000, clump_r2 = 0.001, clump_p = 1, bfile = "/home/camilla/all_phase3.rsid.maf001.dups2.removed", 
                                   plink_bin = "/usr/local/bin/plink")
  
  #change column names back
  dat_for <- dplyr::rename(dat_for, SNP = rsid)
  dat_for <- dplyr::rename(dat_for, pval.exposure = pval)
  
}
                                  
##Perform MR
res_for <- mr(dat_for)
p1 <- mr_scatter_plot(res_for, dat_for)
p1[[1]]
#plotting
res_for <- mr(dat_for, method_list = c("mr_weighted_median", "mr_ivw", "mr_simple_mode", "mr_weighted_mode"))
res_for <- mr(dat_for, method_list = c("mr_ivw"))
jpeg(file= output_file_path %&% "mbg11316" %&% "_vs_" %&% phenotype2 %&% "_scatterplot.png", width = 660, height = 300)
p1[[1]]
dev.off()
                                  

chd_exp_dat <- read_exposure_data(
  filename = filename_string,
  #clump = TRUE,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "panukb_beta",
  se_col = "panukb_se",
  #eaf_col = "eaf",
  effect_allele_col = "panukb_effect",
  other_allele_col = "panukb_other",
  pval_col = "panukb_p",
  #min_pval = 5e-08,
  chr_col = "CHR",
  pos_col = "POS"
)

chd_exp_dat$exposure <- phenotype2

chd_exp_dat <- filter(chd_exp_dat, pval.exposure <= 5e-08)

#chd_out_dat <- filter(chd_out_dat, pval.exposure < 5e-08)

microbiome_out_dat <- read_outcome_data(
  filename = filename_string, ##output file name from python script
  #clump = TRUE,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "micro_beta",
  se_col = "micro_se",
  effect_allele_col = "micro_effect",
  other_allele_col = "micro_other",
  pval_col = "micro_p",
  #min_pval = 1e-04,
  chr_col = "CHR",
  pos_col = "POS"
)

microbiome_out_dat$outcome <- from_dict
#microbiome_out_dat <- filter(microbiome_out_dat, pval.outcome < 1e-04)

##HARMONISE REVERSED DATA (CREATES NEW COMBINED DATA FRAME WITH EXPOSURE AND OUTCOME DATA)
dat_rev <- harmonise_data(chd_exp_dat, microbiome_out_dat)

if (dim(dat_rev)[1] != 0) {
  dat_rev <- dplyr::rename(dat_rev, rsid = SNP)
  dat_rev <- dplyr::rename(dat_rev, pval = pval.exposure)
  
  dat_rev <- ld_clump_local_hacked(dat_rev, clump_kb = 10000, clump_r2 = 0.001, clump_p = 1, bfile = "/home/camilla/all_phase3.rsid.maf001.dups2.removed", 
                                   plink_bin = "/usr/local/bin/plink")
  
  #change column names back
  dat_rev <- dplyr::rename(dat_rev, SNP = rsid)
  dat_rev <- dplyr::rename(dat_rev, pval.exposure = pval)
}

##Perform MR
res_rev <- mr(dat_rev)
p1 <- mr_scatter_plot(res_rev, dat_rev)
p1[[1]]
#plotting
res_rev <- mr(dat_rev, method_list = c("mr_weighted_median", "mr_ivw", "mr_simple_mode", "mr_weighted_mode"))
res_rev <- mr(dat_rev, method_list = c("mr_ivw"))
jpeg(file= output_file_path %&% "mbg11316" %&% "_vs_" %&% phenotype2 %&% "_scatterplot.png", width = 660, height = 300)
p1[[1]]
dev.off()

tiff(filename = "Rplottry.tiff",
     width = 480, height = 480, units = "px", pointsize = 12,
     #compression = c("none", "rle", "lzw", "jpeg", "zip", "lzw+p", "zip+p"),
     bg = "white", res = NA,
     type = c("cairo"))