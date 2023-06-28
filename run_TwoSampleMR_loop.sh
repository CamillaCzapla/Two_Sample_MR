##for loop for TwoSampleMR (splitting jobs)

for n in {1..12}
do nohup Rscript /home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/run_TwoSampleMR_intersection_then_clumping.R --mbfile $n > ${n}-run_Ischemic_Heart_Disease_EURclump_outcomeunfiltered_TwoSampleMR.out &
done