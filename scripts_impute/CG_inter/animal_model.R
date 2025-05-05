# Loading the necessary libraries
library(cmdstanr)
options(brms.backend = "cmdstanr")  # Set cmdstanr as the backend for brms
library(brms)
library(Matrix)


# Command line arguments for SLURM array task
args <- commandArgs(trailingOnly = TRUE)
st <- as.numeric(args[1])  # Start gene index
end <- as.numeric(args[2])  # End gene index

temp_dir = "/projects/ag-demeaux/labraha3/EpiDom_A_lyrata/Animal_model/tmp"
Sys.setenv(TMPDIR = temp_dir)

# Set the directory where the input data is located
setwd("/projects/ag-demeaux/labraha3/EpiDom_A_lyrata/Animal_model/inter_impute/inter_CG/input_files")

# Load the data without setting row names
df_transposed_inter <- read.csv("gBM_CG_for_animal_model_inter_04_12.csv", header = TRUE, row.names = 1, check.names = FALSE, sep= ",")


Amat_inter = read.csv("Amat_inter_DNA_methylation.csv", header = TRUE, check.names = FALSE)
Dmat_inter = readMM("Dmat_inter_DNA_methylation.txt")

a = Amat_inter[, 1]
Amat_inter = Amat_inter[,-1]
colnames(Amat_inter) = a
rownames(Amat_inter) = a
colnames(Dmat_inter) = colnames(Amat_inter)
rownames(Dmat_inter) = rownames(Amat_inter)

# Load random effects
batch_info = read.csv("batch_info_methyl_inter.csv", header = T, sep = ",")
dam_info = read.csv("dam_info_methyl_inter.csv", header = T)
tray_info = read.csv("tray_info_methyl_inter.csv", header = T)
samplingdate_info = read.csv("samplingdata_info_methyl_inter.csv", header = T)

colnames(batch_info)[1] = "animal"
colnames(dam_info)[1] = "animal"
colnames(tray_info)[1] = "animal"
colnames(samplingdate_info)[1] = "animal"

Amat = Amat_inter
Dmat = Dmat_inter

variancesCalc = function(myfile, st, end, batch_info, dam_info, tray_info, samplingdate_info, outdir, individual_files = FALSE) {
  out1 = list()
  counts = myfile[, st:end]
  for (iter in 1:length(colnames(counts))) {
    temp = as.data.frame(counts[, iter])
    gene = colnames(counts)[iter]
    output_gene_file = paste0(outdir, "/output_inter_CG_", gene, ".csv")
    if (individual_files && file.exists(output_gene_file)) {
      next
    }
    colnames(temp) = c("counts")
    temp$animal = as.character(rownames(counts))
    temp$dom = as.character(rownames(counts))
    temp$counts = as.numeric(temp$counts)
    temp$batch = batch_info$batch[match(temp$animal, batch_info$animal)]
    temp$batch = factor(temp$batch)
    temp$dam = dam_info$dam[match(temp$animal, dam_info$animal)]
    temp$dam = factor(temp$dam)
    temp$tray = tray_info$tray[match(temp$animal, tray_info$animal)]
    temp$tray = factor(temp$tray)
    temp$samplingdate = samplingdate_info$samplingdate[match(temp$animal, samplingdate_info$animal)]
    temp$samplingdate = factor(temp$samplingdate)
    print(head(temp))
    
    myformula = brmsformula(counts ~ 1 + (1 | gr(animal, cov = Amat)) + (1 | gr(dom, cov = Dmat)) + (1 | dam) + (1 | batch) + (1 | tray) + (1 | samplingdate))
    fit = brm(myformula, data = temp, family = skew_normal(), data2 = list(Amat = Amat, Dmat = Dmat), chains = 4, threads = threading(4), cores = 4, iter = 4000, warmup = 1000, thin = 2)
    
    Va = median(unlist(VarCorr(fit, summary = FALSE)$animal))
    Vd = median(unlist(VarCorr(fit, summary = FALSE)$dom))
    Vr = median(unlist(VarCorr(fit, summary = FALSE)$residual))
    Vm = median(unlist(VarCorr(fit, summary = FALSE)$dam))
    Vb = median(unlist(VarCorr(fit, summary = FALSE)$batch))  
    Vt = median(unlist(VarCorr(fit, summary = FALSE)$tray))
    Vs = median(unlist(VarCorr(fit, summary = FALSE)$samplingdate))
    ha = round(Va / (Va + Vm + Vr + Vd + Vb + Vt + Vs), digits = 2)
    hd = round(Vd / (Va + Vm + Vr + Vd + Vb + Vt + Vs), digits = 2)
    hr = round(Vr / (Va + Vm + Vr + Vd + Vb + Vt + Vs), digits = 2)
    hm = round(Vm / (Va + Vm + Vr + Vd + Vb + Vt + Vs), digits = 2)
    hb = round(Vb / (Va + Vm + Vr + Vd + Vb + Vt + Vs), digits = 2)
    ht = round(Vt / (Va + Vm + Vr + Vd + Vb + Vt + Vs), digits = 2)
    hs = round(Vs / (Va + Vm + Vr + Vd + Vb + Vt + Vs), digits = 2)
    outp = c(as.character(gene), Va, Vd, Vm, Vr, Vb, Vt, Vs, ha, hd, hm, hr, hb, ht, hs)
    out1[[iter]] = outp
    if (individual_files) {
      individual_output = data.frame(t(outp))
      write.table(individual_output, file = output_gene_file, sep = "\t", row.names = FALSE, col.names = TRUE)
    }
  }
  out1_df = do.call(rbind, out1)
  write.table(out1_df, file = paste0(outdir, "/output_", st, "_", end, ".csv"), sep = "\t", row.names = FALSE, col.names = TRUE)
  return(out1_df)
}

## Apply the function using the start and end gene indices
variancesCalc(df_transposed_inter, st, end, batch_info, dam_info, tray_info, samplingdate_info, outdir = "/projects/ag-demeaux/labraha3/EpiDom_A_lyrata/Animal_model/inter_impute/inter_CG/", individual_files = TRUE)