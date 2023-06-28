# Two_Sample_MRBidirectional MR can help determine causality between risk factors and disease. When Microbiome genera abundancy phenotypes were the exposure, the SNPs are filtered by 1e-05. When cardiometabolic or sleep phenotypes were the exposure, the SNPs are filtered by 5e-08.

For ~200 genera of gut microbiota, TwoSampleMR can be run in a loop format with the shell script included. For other traits, I used a different script. Each R script calls on a python script to create an intersection file and edit the pval column. *PANUKBB stores pval as either lnpval or neglog10pval which must be converted for TwoSamepleMR

loop TwoSampleMR script: /home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/run_TwoSampleMR_intersection_then_clumping.R

python script for loop: /home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/00_mibiogen_panukbb_4_twosampleMR_heart.py

once through TwoSampleMR script: /home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/run_TwoSampleMR_chrpos.R

python script for once through: /home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/panukbb_rsid_from_chrpos.py

*python script must be edited for the amount of columns in the sum stats file

*R script variables must be edited for the names of the sumstats files and phenotype descriptions

link to MR project summary Google Doc: https://drive.google.com/drive/folders/1edSytZiMqrDo-g1mreTJ7r2kSfdOlU7X?usp=sharing
