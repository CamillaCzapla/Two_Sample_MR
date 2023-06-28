
##running TwoSampleMR through all files in a directory
setwd('/home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics')

suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(argparse))
suppressMessages(library(TwoSampleMR))
suppressMessages(library(ieugwasr))
'%&%' = function(a,b) paste (a,b,sep='')

parser <- ArgumentParser()
parser$add_argument('--mbfile', help='microbiome text file number - expecting a number between 1 and 12')
args <- parser$parse_args()

#dictionary of study accession number in filename to microbiome phenotype
microbiome_phenotypes <- c(
  "GCST90016908" = "Gut microbiota abundance (class Actinobacteria id.419)",
  "GCST90016909" = "Gut microbiota abundance (class Alphaproteobacteria id.2379)",
  "GCST90016910" = "Gut microbiota abundance (class Bacilli id.1673)",
  "GCST90016911" = "Gut microbiota abundance (class Bacteroidia id.912)",
  "GCST90016912" = "Gut microbiota abundance (class Betaproteobacteria id.2867)",
  "GCST90016913" = "Gut microbiota abundance (class Clostridia id.1859)",
  "GCST90016914" = "Gut microbiota abundance (class Coriobacteriia id.809)",
  "GCST90016915" = "Gut microbiota abundance (class Deltaproteobacteria id.3087)",
  "GCST90016916" = "Gut microbiota abundance (class Erysipelotrichia id.2147)", 
  "GCST90016917" = "Gut microbiota abundance (class Gammaproteobacteria id.3303)",
  "GCST90016918" = "Gut microbiota abundance (class Lentisphaeria id.2250)",
  "GCST90016919" = "Gut microbiota abundance (class Melainabacteria id.1589)", 
  "GCST90016920" = "Gut microbiota abundance (class Methanobacteria id.119)",
  "GCST90016921" = "Gut microbiota abundance (class Mollicutes id.3920)",
  "GCST90016922" = "Gut microbiota abundance (class Negativicutes id.2164)",
  "GCST90016923" = "Gut microbiota abundance (class Verrucomicrobiae id.4029)",
  "GCST90016924" = "Gut microbiota abundance (family Acidaminococcaceae id.2166)",
  "GCST90016925" = "Gut microbiota abundance (family Actinomycetaceae id.421)",
  "GCST90016926" = "Gut microbiota abundance (family Alcaligenaceae id.2875)",
  "GCST90016927" = "Gut microbiota abundance (family Bacteroidaceae id.917)",
  "GCST90016928" = "Gut microbiota abundance (family Bacteroidales S24 7group id.11173)",
  "GCST90016929" = "Gut microbiota abundance (family Bifidobacteriaceae id.433)",
  "GCST90016930" = "Gut microbiota abundance (family Christensenellaceae id.1866)",
  "GCST90016931" = "Gut microbiota abundance (family Clostridiaceae1 id.1869)", 
  "GCST90016932" = "Gut microbiota abundance (family Clostridiales vadin BB60 group id.11286)",
  "GCST90016933" = "Gut microbiota abundance (family Coriobacteriaceae id.811)",
  "GCST90016934" = "Gut microbiota abundance (family Defluviitaleaceae id.1924)",
  "GCST90016935" = "Gut microbiota abundance (family Desulfovibrionaceae id.3169)",
  "GCST90016936" = "Gut microbiota abundance (family Enterobacteriaceae id.3469)",
  "GCST90016937" = "Gut microbiota abundance (family Erysipelotrichaceae id.2149)",
  "GCST90016938"= "Gut microbiota abundance (family Family XI id.1936)", 
  "GCST90016939" = "Gut microbiota abundance (family Family XIII id.1957)",
  "GCST90016940" = "Gut microbiota abundance (family Lachnospiraceae id.1987)",
  "GCST90016941" = "Gut microbiota abundance (family Lactobacillaceae id.1836)",
  "GCST90016942" = "Gut microbiota abundance (family Methanobacteriaceae id.121)",
  "GCST90016943" = "Gut microbiota abundance (family Oxalobacteraceae id.2966)",
  "GCST90016944" = "Gut microbiota abundance (family Pasteurellaceae id.3689)",
  "GCST90016945" = "Gut microbiota abundance (family Peptococcaceae id.2024)",
  "GCST90016946" = "Gut microbiota abundance (family Peptostreptococcaceae id.2042)",
  "GCST90016947" = "Gut microbiota abundance (family Porphyromonadaceae id.943)",
  "GCST90016948" = "Gut microbiota abundance (family Prevotellaceae id.960)",
  "GCST90016949" = "Gut microbiota abundance (family Rhodospirillaceae id.2717)",
  "GCST90016950" = "Gut microbiota abundance (family Rikenellaceae id.967)",
  "GCST90016951" = "Gut microbiota abundance (family Ruminococcaceae id.2050)",
  "GCST90016952" = "Gut microbiota abundance (family Streptococcaceae id.1850)",
  "GCST90017058" = "Gut microbiota abundance (genus Ruminococcaceae UCG010 id.11367)",
  "GCST90017059" = "Gut microbiota abundance (genus Ruminococcaceae UCG011 id.11368)",
  "GCST90017060" = "Gut microbiota abundance (genus Ruminococcaceae UCG013 id.11370)",
  "GCST90017061" = "Gut microbiota abundance (genus Ruminococcaceae UCG014 id.11371)",
  "GCST90017062" = "Gut microbiota abundance (genus Ruminococcus1 id.11373)",
  "GCST90017063" = "Gut microbiota abundance (genus Ruminococcus2 id.11374)",
  "GCST90017064" = "Gut microbiota abundance (genus Ruminococcus gauvreauii group id.11342)",
  "GCST90017065" = "Gut microbiota abundance (genus Ruminococcus gnavus group id.14376)",
  "GCST90017066" = "Gut microbiota abundance (genus Ruminococcus torques group id.14377)",
  "GCST90017067" = "Gut microbiota abundance (genus Sellimonas id.14369)",
  "GCST90017068" = "Gut microbiota abundance (genus Senegalimassilia id.11160)",
  "GCST90017069" = "Gut microbiota abundance (genus Slackia id.825)",
  "GCST90017070" = "Gut microbiota abundance (genus Streptococcus id.1853)",
  "GCST90017071" = "Gut microbiota abundance (genus Subdoligranulum id.2070)",
  "GCST90017072" = "Gut microbiota abundance (genus Sutterella id.2896)",
  "GCST90017073" = "Gut microbiota abundance (genus Terrisporobacter id.11348)",
  "GCST90017074" = "Gut microbiota abundance (genus Turicibacter id.2162)",
  "GCST90017075" = "Gut microbiota abundance (genus Tyzzerella3 id.11335)",
  "GCST90017076" = "Gut microbiota abundance (unknown genus id.1000000073)",
  "GCST90017077" = "Gut microbiota abundance (unknown genus id.1000001215)",
  "GCST90017078" = "Gut microbiota abundance (unknown genus id.1000005472)",
  "GCST90017079" = "Gut microbiota abundance (unknown genus id.1000005479)",
  "GCST90017080" = "Gut microbiota abundance (unknown genus id.1000006162)",
  "GCST90017081" = "Gut microbiota abundance (unknown genus id.1868)",
  "GCST90017082" = "Gut microbiota abundance (unknown genus id.2001)",
  "GCST90017083" = "Gut microbiota abundance (unknown genus id.2041)",
  "GCST90017084" = "Gut microbiota abundance (unknown genus id.2071)",
  "GCST90017085" = "Gut microbiota abundance (unknown genus id.2755)",
  "GCST90017086" = "Gut microbiota abundance (unknown genus id.826)",
  "GCST90017087" = "Gut microbiota abundance (unknown genus id.959)",
  "GCST90017033" = "Gut microbiota abundance (genus Methanobrevibacter id.123)",
  "GCST90017034" = "Gut microbiota abundance (genus Odoribacter id.952)",
  "GCST90017035" = "Gut microbiota abundance (genus Olsenella id.822)",
  "GCST90017036" = "Gut microbiota abundance (genus Oscillibacter id.2063)",
  "GCST90017037" = "Gut microbiota abundance (genus Oscillospira id.2064)",
  "GCST90017038" = "Gut microbiota abundance (genus Oxalobacter id.2978)",
  "GCST90017039" = "Gut microbiota abundance (genus Parabacteroides id.954)",
  "GCST90017040" = "Gut microbiota abundance (genus Paraprevotella id.962)",
  "GCST90017041" = "Gut microbiota abundance (genus Parasutterella id.2892)",
  "GCST90017042" = "Gut microbiota abundance (genus Peptococcus id.2037)",
  "GCST90017043" = "Gut microbiota abundance (genus Phascolarctobacterium id.2168)",
  "GCST90017044" = "Gut microbiota abundance (genus Prevotella7 id.11182)",
  "GCST90017045" = "Gut microbiota abundance (genus Prevotella9 id.11183)",
  "GCST90017046" = "Gut microbiota abundance (genus Rikenellaceae RC9 gut group id.11191)",
  "GCST90017047" = "Gut microbiota abundance (genus Romboutsia id.11347)",
  "GCST90017048" = "Gut microbiota abundance (genus Roseburia id.2012)",
  "GCST90017049" = "Gut microbiota abundance (genus Ruminiclostridium5 id.11355)",
  "GCST90017050" = "Gut microbiota abundance (genus Ruminiclostridium6 id.11356)",
  "GCST90017051" = "Gut microbiota abundance (genus Ruminiclostridium9 id.11357)",
  "GCST90017052" = "Gut microbiota abundance (genus Ruminococcaceae NK4A214 group id.11358)",
  "GCST90017053" = "Gut microbiota abundance (genus Ruminococcaceae UCG002 id.11360)",
  "GCST90017054" = "Gut microbiota abundance (genus Ruminococcaceae UCG003 id.11361)",
  "GCST90017055" = "Gut microbiota abundance (genus Ruminococcaceae UCG004 id.11362)",
  "GCST90017056" = "Gut microbiota abundance (genus Ruminococcaceae UCG005 id.11363)",
  "GCST90017057" = "Gut microbiota abundance (genus Ruminococcaceae UCG009 id.11366)",
  "GCST90017008" = "Gut microbiota abundance (genus Family XIII AD3011 group id.11293)",
  "GCST90017009" = "Gut microbiota abundance (genus Family XIII UCG001 id.11294)",
  "GCST90017010" = "Gut microbiota abundance (genus Flavonifractor id.2059)",
  "GCST90017011" = "Gut microbiota abundance (genus Fusicatenibacter id.11305)",
  "GCST90017012" = "Gut microbiota abundance (genus Gordonibacter id.821)",
  "GCST90017013" = "Gut microbiota abundance (genus Haemophilus id.3698)",
  "GCST90017014" = "Gut microbiota abundance (genus Holdemanella id.11393)",
  "GCST90017015" = "Gut microbiota abundance (genus Holdemania id.2157)",
  "GCST90017016" = "Gut microbiota abundance (genus Howardella id.2000)",
  "GCST90017017" = "Gut microbiota abundance (genus Hungatella id.11306)",
  "GCST90017018" = "Gut microbiota abundance (genus Intestinibacter id.11345)",
  "GCST90017019" = "Gut microbiota abundance (genus Intestinimonas id.2062)",
  "GCST90017020" = "Gut microbiota abundance (genus Lachnoclostridium id.11308)",
  "GCST90017021" = "Gut microbiota abundance (genus Lachnospiraceae FCS020 group id.11314)",
  "GCST90017022" = "Gut microbiota abundance (genus Lachnospiraceae NC2004 group id.11316)",
  "GCST90017023" = "Gut microbiota abundance (genus Lachnospiraceae ND3007 group id.11317)",
  "GCST90017024" = "Gut microbiota abundance (genus Lachnospiraceae NK4A136 group id.11319)",
  "GCST90017025" = "Gut microbiota abundance (genus Lachnospiraceae UCG001 id.11321)",
  "GCST90017026" = "Gut microbiota abundance (genus Lachnospiraceae UCG004 id.11324)",
  "GCST90017027" = "Gut microbiota abundance (genus Lachnospiraceae UCG008 id.11328)",
  "GCST90017028" = "Gut microbiota abundance (genus Lachnospiraceae UCG010 id.11330)",
  "GCST90017029" = "Gut microbiota abundance (genus Lachnospira id.2004)",
  "GCST90017030" = "Gut microbiota abundance (genus Lactobacillus id.1837)",
  "GCST90017031" = "Gut microbiota abundance (genus Lactococcus id.1851)",
  "GCST90017032" = "Gut microbiota abundance (genus Marvinbryantia id.2005)",
  "GCST90017127" = "Gut microbiota presence (genus Anaerostipes id.1991)",
  "GCST90017128" = "Gut microbiota presence (family Bacteroidales S24 7group id.11173)",
  "GCST90017129" = "Gut microbiota presence (unknown genus id.1000005479)",
  "GCST90017130" = "Gut microbiota alpha diversity (Inverse Simpson index)",
  "GCST90017131" = "Gut microbiota alpha diversity (Shannon index)",
  "GCST90017132" = "Gut microbiota alpha diversity (Simpson index)",
  "GCST90016964" = "Gut microbiota abundance (genus Alloprevotella id.961)",
  "GCST90016965" = "Gut microbiota abundance (genus Anaerofilum id.2053)",
  "GCST90016966" = "Gut microbiota abundance (genus Anaerostipes id.1991)",
  "GCST90016967" = "Gut microbiota abundance (genus Anaerotruncus id.2054)",
  "GCST90016968" = "Gut microbiota abundance (genus Bacteroides id.918)",
  "GCST90016969" = "Gut microbiota abundance (genus Barnesiella id.944)",
  "GCST90016970" = "Gut microbiota abundance (genus Bifidobacterium id.436)",
  "GCST90016971" = "Gut microbiota abundance (genus Bilophila id.3170)",
  "GCST90016972" = "Gut microbiota abundance (genus Blautia id.1992)",
  "GCST90016973" = "Gut microbiota abundance (genus Butyricicoccus id.2055)",
  "GCST90016974" = "Gut microbiota abundance (genus Butyricimonas id.945)",
  "GCST90016975" = "Gut microbiota abundance (genus Butyrivibrio id.1993)",
  "GCST90016976" = "Gut microbiota abundance (genus Candidatus Soleaferrea id.11350)",
  "GCST90016977" = "Gut microbiota abundance (genus Catenibacterium id.2153)",
  "GCST90016978" = "Gut microbiota abundance (genus Christensenellaceae R 7group id.11283)",
  "GCST90016979" = "Gut microbiota abundance (genus Clostridium innocuum group id.14397)",
  "GCST90016980" = "Gut microbiota abundance (genus Clostridium sensustricto1 id.1873)",
  "GCST90016981" = "Gut microbiota abundance (genus Collinsella id.815)",
  "GCST90016982" = "Gut microbiota abundance (genus Coprobacter id.949)",
  "GCST90016983" = "Gut microbiota abundance (genus Coprococcus1 id.11301)",
  "GCST90016984" = "Gut microbiota abundance (genus Coprococcus2 id.11302)",
  "GCST90016985" = "Gut microbiota abundance (genus Coprococcus3 id.11303)",
  "GCST90016986" = "Gut microbiota abundance (genus Defluviitaleaceae UCG011 id.11287)",
  "GCST90016987" = "Gut microbiota abundance (genus Desulfovibrio id.3173)",
  "GCST90016988" = "Gut microbiota abundance (genus Dialister id.2183)",
  "GCST90016989" = "Gut microbiota abundance (genus Dorea id.1997)",
  "GCST90016990" = "Gut microbiota abundance (genus Eggerthella id.819)",
  "GCST90016991" = "Gut microbiota abundance (genus Eisenbergiella id.11304)",
  "GCST90016992" = "Gut microbiota abundance (genus Enterorhabdus id.820)",
  "GCST90016993" = "Gut microbiota abundance (genus Erysipelatoclostridium id.11381)",
  "GCST90016994" = "Gut microbiota abundance (genus Erysipelotrichaceae UCG003 id.11384)",
  "GCST90016995" = "Gut microbiota abundance (genus Escherichia Shigella id.3504)",
  "GCST90016996" = "Gut microbiota abundance (genus Eubacterium brachy group id.11296)",
  "GCST90016997" = "Gut microbiota abundance (genus Eubacterium coprostanoligenes group id.11375)",
  "GCST90016998" = "Gut microbiota abundance (genus Eubacterium eligens group id.14372)",
  "GCST90016999" = "Gut microbiota abundance (genus Eubacterium fissicatena group id.14373)",
  "GCST90017000" = "Gut microbiota abundance (genus Eubacterium hallii group id.11338)",
  "GCST90017001" = "Gut microbiota abundance (genus Eubacterium nodatum group id.11297)",
  "GCST90017002" = "Gut microbiota abundance (genus Eubacterium oxidoreducens group id.11339)",
  "GCST90017003" = "Gut microbiota abundance (genus Eubacterium rectale group id.14374)",
  "GCST90017004" = "Gut microbiota abundance (genus Eubacterium ruminantium group id.11340)",
  "GCST90017005" = "Gut microbiota abundance (genus Eubacterium ventriosum group id.11341)",
  "GCST90017006" = "Gut microbiota abundance (genus Eubacterium xylanophilum group id.14375)",
  "GCST90017007" = "Gut microbiota abundance (genus Faecalibacterium id.2057)",
  "GCST90017088" = "Gut microbiota abundance (genus Veillonella id.2198)",
  "GCST90017089" = "Gut microbiota abundance (genus Victivallis id.2256)",
  "GCST90017090" = "Gut microbiota abundance (order Actinomycetales id.420)",
  "GCST90017091" = "Gut microbiota abundance (order Bacillales id.1674)",
  "GCST90017092" = "Gut microbiota abundance (order Bacteroidales id.913)",
  "GCST90017093" = "Gut microbiota abundance (order Bifidobacteriales id.432)",
  "GCST90017094" = "Gut microbiota abundance (order Burkholderiales id.2874)",
  "GCST90017095" = "Gut microbiota abundance (order Clostridiales id.1863)",
  "GCST90017096" = "Gut microbiota abundance (order Coriobacteriales id.810)",
  "GCST90017097" = "Gut microbiota abundance (order Desulfovibrionales id.3156)",
  "GCST90017098" = "Gut microbiota abundance (order Enterobacteriales id.3468)",
  "GCST90017099" = "Gut microbiota abundance (order Erysipelotrichales id.2148)",
  "GCST90017100" = "Gut microbiota abundance (order Gastranaerophilales id.1591)",
  "GCST90016953" = "Gut microbiota abundance (unknown family id.1000001214)",
  "GCST90016954" = "Gut microbiota abundance (unknown family id.1000005471)",
  "GCST90016955" = "Gut microbiota abundance (unknown family id.1000006161)",
  "GCST90016956" = "Gut microbiota abundance (family Veillonellaceae id.2172)",
  "GCST90016957" = "Gut microbiota abundance (family Verrucomicrobiaceae id.4036)",
  "GCST90016958" = "Gut microbiota abundance (family Victivallaceae id.2255)",
  "GCST90016959" = "Gut microbiota abundance (genus Actinomyces id.423)",
  "GCST90016960" = "Gut microbiota abundance (genus Adlercreutzia id.812)",
  "GCST90016961" = "Gut microbiota abundance (genus Akkermansia id.4037)",
  "GCST90016962" = "Gut microbiota abundance (genus Alistipes id.968)",
  "GCST90016963" = "Gut microbiota abundance (genus Allisonella id.2174)",
  "GCST90017101" = "Gut microbiota abundance (order Lactobacillales id.1800)",
  "GCST90017102" = "Gut microbiota abundance (order Methanobacteriales id.120)",
  "GCST90017103" = "Gut microbiota abundance (order Mollicutes RF9 id.11579)",
  "GCST90017104" = "Gut microbiota abundance (order NB1n id.3953)",
  "GCST90017105" = "Gut microbiota abundance (order Pasteurellales id.3688)",
  "GCST90017106" = "Gut microbiota abundance (order Rhodospirillales id.2667)",
  "GCST90017107" = "Gut microbiota abundance (order Selenomonadales id.2165)",
  "GCST90017108" = "Gut microbiota abundance (order Verrucomicrobiales id.4030)",
  "GCST90017109" = "Gut microbiota abundance (order Victivallales id.2254)",
  "GCST90017110" = "Gut microbiota abundance (phylum Actinobacteria id.400)",
  "GCST90017111" = "Gut microbiota abundance (phylum Bacteroidetes id.905)",
  "GCST90017112" = "Gut microbiota abundance (phylum Cyanobacteria id.1500)",
  "GCST90017113" = "Gut microbiota abundance (phylum Euryarchaeota id.55)",
  "GCST90017114" = "Gut microbiota abundance (phylum Firmicutes id.1672)",
  "GCST90017115" = "Gut microbiota abundance (phylum Lentisphaerae id.2238)",
  "GCST90017116" = "Gut microbiota abundance (phylum Proteobacteria id.2375)",
  "GCST90017117" = "Gut microbiota abundance (phylum Tenericutes id.3919)",
  "GCST90017118" = "Gut microbiota abundance (phylum Verrucomicrobia id.3982)",
  "GCST90017119" = "Gut microbiota presence (genus Turicibacter id.2162)",
  "GCST90017120" = "Gut microbiota presence (genus Ruminiclostridium6 id.11356)",
  "GCST90017121" = "Gut microbiota presence (genus Odoribacter id.952)",
  "GCST90017122" = "Gut microbiota presence (family Peptococcaceae id.2024)",
  "GCST90017123" = "Gut microbiota presence (genus Enterorhabdus id.820)",
  "GCST90017124" = "Gut microbiota presence (genus Coprococcus1 id.11301)",
  "GCST90017125" = "Gut microbiota presence (family Pseudomonadaceae id.3718)",
  "GCST90017126" = "Gut microbiota presence (genus Lachnospiraceae UCG 010 id.11330)"
)

#pulling name of microbiome file from one of the microbiome file lists (there are 12)
mbfile_list <- fread('/home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/MiBioGen_' %&% args$mbfile %&% '.txt',
                     header = FALSE) %>% pull(V1)

#for loop for each file
#I print file name in order to know which file any errors in the .out file pertain to
for (File.Name in mbfile_list){
  print(File.Name)
  ##need phenotype name from dictionary key
  uniquename <- substr(File.Name, 10, 21)
  from_dict <- microbiome_phenotypes[uniquename]
  
  ##reading in exposure data
  microbiome_exp_dat <- read_exposure_data(
    filename = File.Name,
    #clump = TRUE,
    sep = "\t",
    snp_col = "variant_id",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "p_value",
  )
  
  ##renaming exposure column to file name to identify the genera for each test
  microbiome_exp_dat$exposure <- from_dict
  
  #filtering microbiome exposure data to pval threshold (pvalue in MiBioGen paper was 1e-05)
  microbiome_exp_dat <- filter(microbiome_exp_dat, pval.exposure <= 1e-05)
  
  ##here I will use a python script which uses the rsids in the microbiome data to give PUKBB data rsids using chr/pos in both datasets
  ##all genome build 37 in PUKBB and MIBIOGEN files
  ##python script will also create a temp file with the intersection so I can read in a smaller filtered file into Rscript
  ##PUKBB files too big to read in - Cstack usage error
  
  sumstats = "phecode-411-both_sexes"
  my_command_1 = "python3 /home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/00_mibiogen_panukbb_4_twosampleMR.py --micro /home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/" %&% File.Name %&%" --panukb /home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/" %&% sumstats %&% ".tsv.bgz --output /home/camilla/TwoSampleMR/MiBioGen_Summary_Statistics/" %&% args$mbfile %&% "_intersection"
  system(my_command_1)
  
  ##name of temp file created by python scipt
  filename_string = args$mbfile %&% "_intersection.txt.gz"

  
  #read in outcome locally
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
    #ncontrol_col = "ncontrol",
    #samplesize_col = "samplesize",
    #gene_col = "gene",
    #id_col = "id",
    #min_pval = 5e-08,
    #log_pval = FALSE,
    #chr_col = "chr",
    #pos_col = "pos"
  )
  
  ##indicating outcome phenotype from PUKBB
  chd_out_dat$outcome <- "Ischemic Heart Disease"
  
  ##we DO NOT filter the outcome, only the exposure
  #chd_out_dat <- filter(chd_out_dat, pval.outcome < 5e-08)
  
  ##HARMONISE DATA (Intersection/CREATES NEW COMBINED DATA FRAME WITH EXPOSURE AND OUTCOME DATA)
  dat_for <- harmonise_data(microbiome_exp_dat, chd_out_dat)
  
  #hacked ld_clump_local function with --allow-extra-chr flag
  #We started clumping after the intersection is taken (using pval.exposure!)
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
  
  #if the dimensions of the data frame are not 0 then there are SNPs to work with and the analysis
  #can be run, otherwise error messages
  ##I change column names for the clumping function, and then I change them back.
  if (dim(dat_for)[1] != 0) {
  dat_for <- dplyr::rename(dat_for, rsid = SNP)
  dat_for <- dplyr::rename(dat_for, pval = pval.exposure)
  
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
  
  ##writing files
  if (dim(dat_for)[1] != 0) {
  fwrite(dat_for, file = sumstats %&% "_vs_" %&% from_dict %&% "_MR_data.csv", quote = FALSE, na = "NA")
  }
  if (dim(res_for)[1] != 0) {
  fwrite(res_for, file = sumstats %&% "_vs_" %&% from_dict %&% "_MR_results.csv", quote = FALSE, na = "NA")
  }
  exist_check_het_for <- exists("het_for")
  if (exist_check_het_for == TRUE && dim(het_for)[1] != 0) {
  fwrite(het_for, file = sumstats %&% "_vs_" %&% from_dict %&% "_MR_heterogeneity.csv", quote = FALSE, na = "NA")
  }
  exist_check_plt_for <- exists("plt_for")
  if (exist_check_plt_for == TRUE && dim(plt_for)[1] != 0) {
  fwrite(plt_for, file = sumstats %&% "_vs_" %&% from_dict %&% "_MR_pleiotropy.csv", quote = FALSE, na = "NA")
  }
  exist_check_sin_for <- exists("sin_for")
  if (exist_check_sin_for == TRUE && dim(sin_for)[1] != 0) {
  fwrite(sin_for, file = sumstats %&% "_vs_" %&% from_dict %&% "_MR_singlesnp.csv", quote = FALSE, na = "NA")
  }
  
  #call python script with system() in order to have a smaller file for read_outcome_data()
  #otherwise, a Cstack usage error
  system(my_command_1)
  
  ##reading in previous outcome as exposure data with temp file
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
  )
  
  chd_exp_dat$exposure <- "Ischemic Heart Disease"
  
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
  )
  
  microbiome_out_dat$outcome <- from_dict
  ##We DO NOT filter outcome pval
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
  
  ##Sensitivity analyses
  if (dim(dat_rev)[1] != 0) {
  het_rev <- mr_heterogeneity(dat_rev)
  plt_rev <- mr_pleiotropy_test(dat_rev)
  sin_rev <- mr_singlesnp(dat_rev)
  }
  
  ##saving files one at a time
  if (dim(dat_rev)[1] != 0) {
  fwrite(dat_rev, file = from_dict %&% "_vs_" sumstats %&% "_MR_data.csv", quote = FALSE, na = "NA")
  }
  if (dim(res_rev)[1] != 0) {
  fwrite(res_rev, file = from_dict %&% "_vs_" sumstats %&% "_MR_results.csv", quote = FALSE, na = "NA")
  }
  exist_check_het_rev <- exists("het_rev")
  if (exist_check_het_rev == TRUE && dim(het_rev)[1] != 0) {
  fwrite(het_rev, file = from_dict %&% "_vs_" sumstats %&% "_MR_heterogeneity.csv", quote = FALSE, na = "NA")
  }
  exist_check_plt_rev <- exists("plt_rev")
  if (exist_check_plt_rev == TRUE && dim(plt_rev)[1] != 0) {
  fwrite(plt_rev, file = from_dict %&% "_vs_" sumstats %&% "_MR_pleiotropy.csv", quote = FALSE, na = "NA")
  }
  exist_check_sin_rev <- exists("sin_rev")
  if (exist_check_sin_rev == TRUE && dim(sin_rev)[1] != 0) {
  fwrite(sin_rev, file = from_dict %&% "_vs_" sumstats %&% "_MR_singlesnp.csv", quote = FALSE, na = "NA")
  }
  
  ##filtering for bonferroni all five tests with pval <= 0.05
  ##in progress
  
  #clear variables for next loop iteration
  rm(dat_for, res_for, het_for, plt_for, sin_for, dat_rev, res_rev, het_rev, plt_rev, sin_rev)
  
}