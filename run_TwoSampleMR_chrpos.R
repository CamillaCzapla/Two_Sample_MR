##running TwoSampleMR for one phenotype of interest

setwd('/home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics')

suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
#suppressMessages(library(argparse))
suppressMessages(library(TwoSampleMR))
suppressMessages(library(ieugwasr))
'%&%' = function(a,b) paste (a,b,sep='')

#output_file_path = "/home/camilla/TwoSampleMR/MR"
output_file_path = "/home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/results_other/"
#File.Name = "33462485-GCST90017115-EFO_0007874-Build37.f.tsv.gz"
File.Name = "continuous-1160-both_sexes.tsv.bgz"
File.Name2 = "phecode-411.4-both_sexes.tsv.bgz"
sumstats1 = "continuous-1160-both_sexes"
sumstats2 = "phecode-411.4-both_sexes"
phenotype1 = "sleep duration"
phenotype2 = "coronary atherosclerosis"

my_command_1 = "python3 /home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/panukbb_rsid_from_chrpos.py --bim /home/camilla/all_phase3.rsid.maf001.dups2.removed.bim --first /home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/" %&% File.Name %&% " --second /home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/" %&% File.Name2 %&% " --output /home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/intersection"

#call python script with system() in order to have a smaller file with necessary info for read_outcome_data()
#otherwise, a Cstack usage error
system(my_command_1)

##reading in exposure data
first_exp_dat <- read_exposure_data(
  filename = "intersection.txt.gz",
  #clump = TRUE,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "first_beta",
  se_col = "first_se",
  effect_allele_col = "first_effect",
  other_allele_col = "first_other",
  pval_col = "first_p",
  chr_col = "CHR",
  pos_col = "POS"
)

#first_exp_dat <- dplyr::rename(first_exp_dat, chr = chr.exposure)
#first_exp_dat <- dplyr::rename(first_exp_dat, pos = SNP)

##renaming exposure column to phenotype name to identify the genera for each test
first_exp_dat$exposure <- phenotype1

first_exp_dat <- filter(first_exp_dat, pval.exposure <= 5e-08)

second_out_dat <- read_outcome_data(
  #snps = microbiome_exp_dat$SNP,
  filename = "intersection.txt.gz",
  #filename = "intersection.tsv",
  #snps = microbiome_exp_dat$SNP,
  sep = "\t",
  #phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "second_beta",
  se_col = "second_se",
  #eaf_col = "eaf",
  effect_allele_col = "second_effect",
  other_allele_col = "second_other",
  pval_col = "second_p",
  #units_col = "units",
  #ncase_col = "ncase",
  #ncontrol_col = "ncontrol",
  #samplesize_col = "samplesize",
  #gene_col = "gene",
  #id_col = "id",
  #min_pval = 5e-08,
  #log_pval = FALSE,
  chr_col = "CHR",
  pos_col = "POS"
)


#second_out_dat <- dplyr::rename(second_out_dat, chr = chr.outcome)
#second_out_dat <- dplyr::rename(second_out_dat, pos = SNP)
second_out_dat$outcome <- phenotype2

#harmonise_data_hacked <- function (exposure_dat, outcome_dat, action = 2) 
#{
#  stopifnot(all(action %in% 1:3))
#  #check_required_columns(exposure_dat, "exposure")
#  #check_required_columns(outcome_dat, "outcome")
#  res.tab <- merge(outcome_dat, exposure_dat, by = c("chr", "pos"))
#  #old one
#  #res.tab <- merge(outcome_dat, exposure_dat, by = "SNP")
#  ncombinations <- length(unique(res.tab$id.outcome))
#  if (length(action) == 1) {
#    action <- rep(action, ncombinations)
#  }
#  else if (length(action) != ncombinations) {
#    stop("Action argument must be of length 1 (where the same action will be used for all outcomes), or number of unique id.outcome values (where each outcome will use a different action value)")
#  }
#  res.tab <- harmonise_cleanup_variables(res.tab)
#  if (nrow(res.tab) == 0) {
#    return(res.tab)
#  }
#  d <- data.frame(id.outcome = unique(res.tab$id.outcome), 
#                  action = action)
#  res.tab <- merge(res.tab, d, by = "id.outcome")
#  combs <- subset(res.tab, !duplicated(paste(id.exposure, id.outcome)), 
#                  select = c(id.exposure, id.outcome))
#  fix.tab <- list()
#  mr_cols <- c("beta.exposure", "beta.outcome", "se.exposure", 
##               "se.outcome")
#  for (i in 1:nrow(combs)) {
#    x <- subset(res.tab, id.exposure == combs$id.exposure[i] & 
#                  id.outcome == combs$id.outcome[i])
#    message("Harmonising ", x$exposure[1], " (", x$id.exposure[1], 
#            ") and ", x$outcome[1], " (", x$id.outcome[1], ")")
#    x <- harmonise(x, 0.08, x$action[1])
#    attr(x, "log")[["candidate_variants"]] <- sum(exposure_dat$id.exposure == 
#                                                    x$id.exposure[1])
#    attr(x, "log")[["variants_absent_from_reference"]] <- sum(exposure_dat$id.exposure == 
#                                                                x$id.exposure[1]) - nrow(x)
#    x$mr_keep[apply(x[, mr_cols], 1, function(y) any(is.na(y)))] <- FALSE
#    attr(x, "log")[["total_variants"]] <- nrow(x)
#    attr(x, "log")[["total_variants_for_mr"]] <- sum(x$mr_keep)
#    attr(x, "log")[["proxy_variants"]] <- ifelse(is.null(x$proxy.outcome), 
#                                                 0, sum(x$proxy.outcome, na.rm = TRUE))
#    fix.tab[[i]] <- x
#  }
#  jlog <- plyr::rbind.fill(lapply(fix.tab, function(x) attr(x, 
#                                                            "log")))
#  fix.tab <- plyr::rbind.fill(fix.tab)
#  attr(fix.tab, "log") <- jlog
#  if (!"samplesize.outcome" %in% names(fix.tab)) {
#    fix.tab$samplesize.outcome <- NA
#  }
#  return(fix.tab)
#}

dat_for <- harmonise_data(first_exp_dat, second_out_dat)

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
##Sensitivity analyses
if (dim(dat_for)[1] != 0) {
  het_for<- mr_heterogeneity(dat_for)
  plt_for<- mr_pleiotropy_test(dat_for)
  sin_for<- mr_singlesnp(dat_for)
}

##saving files one at a time
if (dim(dat_for)[1] != 0) {
  fwrite(dat_for, file = output_file_path %&% sumstats1 %&% "_vs_" %&% sumstats2 %&% "_MR_data.csv", quote = FALSE, na = "NA")
}
if (dim(res_for)[1] != 0) {
  fwrite(res_for, file = output_file_path %&% sumstats1 %&% "_vs_" %&% sumstats2 %&% "_MR_results.csv", quote = FALSE, na = "NA")
}
exist_check_het_for <- exists("het_for")
if (exist_check_het_for == TRUE && dim(het_for)[1] != 0) {
  fwrite(het_for, file = output_file_path %&% sumstats1 %&% "_vs_" %&% sumstats2 %&% "_MR_heterogeneity.csv", quote = FALSE, na = "NA")
}
exist_check_plt_for <- exists("plt_for")
if (exist_check_plt_for == TRUE && dim(plt_for)[1] != 0) {
  fwrite(plt_for, file = output_file_path %&% sumstats1 %&% "_vs_" %&% sumstats2 %&% "_MR_pleiotropy.csv", quote = FALSE, na = "NA")
}
exist_check_sin_for <- exists("sin_for")
if (exist_check_sin_for == TRUE && dim(sin_for)[1] != 0) {
  fwrite(sin_for, file = output_file_path %&% sumstats1 %&% "_vs_" %&% sumstats2 %&% "_MR_singlesnp.csv", quote = FALSE, na = "NA")
}

##Now we do the reverse

second_exp_dat <- read_exposure_data(
  filename = "intersection.txt.gz",
  #clump = TRUE,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "second_beta",
  se_col = "second_se",
  effect_allele_col = "second_effect",
  other_allele_col = "second_other",
  pval_col = "second_p",
  chr_col = "CHR",
  pos_col = "POS"
)

second_exp_dat$exposure <- phenotype2

second_exp_dat <- filter(second_exp_dat, pval.exposure <= 5e-08)

first_out_dat <- read_outcome_data(
  filename = "intersection.txt.gz", ##output file name from python script
  #clump = TRUE,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "first_beta",
  se_col = "first_se",
  effect_allele_col = "first_effect",
  other_allele_col = "first_other",
  pval_col = "first_p",
  #min_pval = 1e-04,
  chr_col = "CHR",
  pos_col = "POS"
)

first_out_dat$outcome <- phenotype1

##HARMONISE REVERSED DATA (CREATES NEW COMBINED DATA FRAME WITH EXPOSURE AND OUTCOME DATA)
dat_rev <- harmonise_data(second_exp_dat, first_out_dat)

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

##Sensitivity analyses
if (dim(dat_rev)[1] != 0) {
  het_rev <- mr_heterogeneity(dat_rev)
  plt_rev <- mr_pleiotropy_test(dat_rev)
  sin_rev <- mr_singlesnp(dat_rev)
}

##saving files one at a time
if (dim(dat_rev)[1] != 0) {
  fwrite(dat_rev, file = output_file_path %&% sumstats2 %&% "_vs_" %&% sumstats1 %&% "_MR_data.csv", quote = FALSE, na = "NA")
}
if (dim(res_rev)[1] != 0) {
  fwrite(res_rev, file = output_file_path %&% sumstats2 %&% "_vs_" %&% sumstats1 %&% "_MR_results.csv", quote = FALSE, na = "NA")
}
exist_check_het_rev <- exists("het_rev")
if (exist_check_het_rev == TRUE && dim(het_rev)[1] != 0) {
  fwrite(het_rev, file = output_file_path %&% sumstats2 %&% "_vs_" %&% sumstats1 %&% "_MR_heterogeneity.csv", quote = FALSE, na = "NA")
}
exist_check_plt_rev <- exists("plt_rev")
if (exist_check_plt_rev == TRUE && dim(plt_rev)[1] != 0) {
  fwrite(plt_rev, file = output_file_path %&% sumstats2 %&% "_vs_" %&% sumstats1 %&% "_MR_pleiotropy.csv", quote = FALSE, na = "NA")
}
exist_check_sin_rev <- exists("sin_rev")
if (exist_check_sin_rev == TRUE && dim(sin_rev)[1] != 0) {
  fwrite(sin_rev, file = output_file_path %&% sumstats2 %&% "_vs_" %&% sumstats1 %&% "_MR_singlesnp.csv", quote = FALSE, na = "NA")
}

#plotting
res_for <- mr(dat_for, method_list = c("mr_weighted_median", "mr_ivw", "mr_simple_mode", "mr_weighted_mode"))
p1 <- mr_scatter_plot(res_for, dat_for)
png(file="/home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/results_other/" %&% sumstats1 %&% "_vs_" %&% sumstats2 %&% "_scatterplot.png",
    width=600, height=550)
p1[[1]]
dev.off()

res_rev <- mr(dat_rev, method_list = c("mr_weighted_median", "mr_ivw", "mr_simple_mode", "mr_weighted_mode"))
p2 <- mr_scatter_plot(res_rev, dat_rev)
png(file="/home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/results_other/" %&% sumstats2 %&% "_vs_" %&% sumstats1 %&% "_scatterplot.png",
    width=600, height=550)
p2[[1]]
dev.off()