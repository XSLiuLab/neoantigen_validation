# Prepare data for neoantigen prediction


library(data.table)
# Melanoma ----------------------------------------------------------------

maf.mela = fread("Neoantigen-Melanoma/neo_melanoma.maf.txt")
maf.mela[, Tumor_Sample_Barcode := paste0("sample", Tumor_Sample_Barcode)]
colnames(maf.mela) = c("Tumor_Sample_Barcode",
                       "Hugo_Symbol",
                       "Transcript_ID",
                       "Variant_Classification",
                       "Variant_Type",
                       "NCBI_Build",
                       "Chromosome",
                       "Start_position",
                       "End_position",
                       "Strand",
                       "Reference_Allele",
                       "Tumor_Seq_Allele1",
                       "Tumor_Seq_Allele2",
                       "cDNA_Change",
                       "Codon_Change",
                       "Protein_Change")
  

fwrite(maf.mela, file = "neo_mela.maf", sep = "\t")

hla.mela = fread("Neoantigen-Melanoma/hla_melanoma.txt")
colnames(hla.mela) = c("Tumor_Sample_Barcode", "A1", "A2", "B1", "B2")

correct = function(x) {
  ifelse(nchar(x) == 4, paste0("0", x), x)
}

hla.mela$A1 =  correct(hla.mela$A1)
hla.mela$A2 =  correct(hla.mela$A2)
hla.mela$B1 =  correct(hla.mela$B1)
hla.mela$B2 =  correct(hla.mela$B2)

hla.mela[, `:=`(Tumor_Sample_Barcode = paste0("sample", Tumor_Sample_Barcode),
                HLA = paste(
                  paste("HLA-A", A1, sep = "*"),
                  paste("HLA-A", A2, sep = "*"),
                  paste("HLA-B", B1, sep = "*"),
                  paste("HLA-B", B2, sep = "*"),
                  sep = ","
                ))]
hla.mela[, c("A1","A2", "B1", "B2"):=NULL]
hla.mela
fwrite(hla.mela, file = "neo_mela_hla.tsv", sep = "\t")


# GBM ---------------------------------------------------------------------

maf.gbm = fread("Neoantigen-glioblastoma/neo_gbm.maf.txt")
maf.gbm[, Tumor_Sample_Barcode := paste0("sample", Tumor_Sample_Barcode)]
fwrite(maf.gbm, file = "neo_gbm.maf", sep = "\t")

table(maf.gbm$Tumor_Sample_Barcode)

library(readxl)
library(dplyr)
hla.gbm = read_excel("Neoantigen-glioblastoma/41586_2018_792_MOESM5_ESM.xlsx", skip = 4, col_names = FALSE)
head(hla.gbm)
hla.gbm = hla.gbm[, c(1,6)]
colnames(hla.gbm) = c("Tumor_Sample_Barcode", "HLA")
hla.gbm = unique(hla.gbm)

hla.gbm %>% 
  dplyr::filter(HLA != "NA") %>% 
  dplyr::mutate(Tumor_Sample_Barcode = paste0("sample", Tumor_Sample_Barcode), 
                HLA = paste("HLA", paste0(substr(HLA, 1, 1), "*", substring(HLA, 2)), sep = "-")) %>% 
  group_by(Tumor_Sample_Barcode) %>% 
  summarise(HLA = paste(HLA, collapse = ",")) -> hla.gbm
readr::write_tsv(hla.gbm, path = "neo_gbm_hla.tsv")


maf.gbm = maf.gbm[paste0("sample", 1:8), on = "Tumor_Sample_Barcode"]
fwrite(maf.gbm, file = "neo_gbm.maf", sep = "\t")

#upload to HPC
#sync-upload -n neo_* -d '/public/home/liuxs/wangshx/pipeline_test'

