# Load necessary libraries
library(cmdstanr)
options(brms.backend = "cmdstanr")  # Use cmdstanr as backend
library(brms)
library(Matrix)

# Read command line arguments from SLURM script
args <- commandArgs(trailingOnly = TRUE)
st <- as.numeric(args[1])  # Start gene index
end <- as.numeric(args[2])  # End gene index
dataset_type <- args[3]  # "observed" or "permuted"

# Set working directories
temp_dir <- "/projects/ag-demeaux/labraha3/EpiDom_A_lyrata/Animal_model/tmp"
Sys.setenv(TMPDIR = temp_dir)

setwd("/projects/ag-demeaux/labraha3/EpiDom_A_lyrata/Animal_model/intra_impute/intra_CG/input_files")

# Load gene expression data
df_transposed_intra <- read.csv("gBM_CG_for_animal_model_intra_04_12.csv", header = TRUE, row.names = 1, check.names = FALSE, sep= ",")

# Load relatedness matrices
Amat_intra <- read.csv("Amat_intra_DNA_methylation.csv", header = TRUE, check.names = FALSE)
Dmat_intra <- readMM("Dmat_intra_DNA_methylation.txt")

# Process matrices
a <- Amat_intra[, 1]
Amat_intra <- Amat_intra[, -1]
colnames(Amat_intra) <- a
rownames(Amat_intra) <- a
colnames(Dmat_intra) <- colnames(Amat_intra)
rownames(Dmat_intra) <- rownames(Amat_intra)

# Load random effect variables
batch_info <- read.csv("batch_info_methyl_intra.csv", header = TRUE, sep = ",")
dam_info <- read.csv("dam_info_methyl_intra.csv", header = TRUE)
tray_info <- read.csv("tray_info_methyl_intra.csv", header = TRUE)
samplingdate_info <- read.csv("samplingdata_info_methyl_intra.csv", header = TRUE)

# Ensure column names are correctly assigned
colnames(batch_info)[1] <- "animal"
colnames(dam_info)[1] <- "animal"
colnames(tray_info)[1] <- "animal"
colnames(samplingdate_info)[1] <- "animal"

# Use correct matrices based on dataset type
if (dataset_type == "observed") {
    message("Processing observed dataset...")
    Amat <- Amat_intra
    Dmat <- Dmat_intra
    output_dir <- "/projects/ag-demeaux/labraha3/EpiDom_A_lyrata/Animal_model/intra_impute/intra_CG/output_files_random/"
} else if (dataset_type == "permuted") {
    message("Processing permuted dataset...")
    Amat <- diag(nrow(Amat_intra))  # Identity matrix for permuted data
    Dmat <- diag(nrow(Dmat_intra))  # Identity matrix for permuted data
    output_dir <- "/projects/ag-demeaux/labraha3/EpiDom_A_lyrata/Animal_model/intra_impute/intra_CG/permuted_output/"

    # Shuffle individual assignments for permutation test
    set.seed(123)  # Ensure reproducibility
    shuffled_animals <- sample(batch_info$animal)
    batch_info$animal <- shuffled_animals
    dam_info$animal <- shuffled_animals
    tray_info$animal <- shuffled_animals
    samplingdate_info$animal <- shuffled_animals
}

# Function to run variance partitioning
variancesCalc <- function(myfile, st, end, batch_info, dam_info, tray_info, samplingdate_info, outdir, Amat, Dmat, individual_files = FALSE) {
  out1 <- list()
  counts <- myfile[, st:end]

  for (iter in 1:length(colnames(counts))) {
    temp <- as.data.frame(counts[, iter])
    gene <- colnames(counts)[iter]
    output_gene_file <- paste0(outdir, "/output_intra_CG_", gene, ".csv")

    if (individual_files && file.exists(output_gene_file)) {
      next
    }

    colnames(temp) <- c("counts")
    temp$animal <- as.character(rownames(counts))
    temp$dom <- as.character(rownames(counts))
    temp$counts <- as.numeric(temp$counts)

    # Match random effects
    temp$batch <- factor(batch_info$batch[match(temp$animal, batch_info$animal)])
    temp$dam <- factor(dam_info$dam[match(temp$animal, dam_info$animal)])
    temp$tray <- factor(tray_info$tray[match(temp$animal, tray_info$animal)])
    temp$samplingdate <- factor(samplingdate_info$samplingdate[match(temp$animal, samplingdate_info$animal)])

    # Define the model formula
    myformula <- brmsformula(counts ~ 1 + (1 | gr(animal, cov = Amat)) + (1 | gr(dom, cov = Dmat)) + 
                            (1 | dam) + (1 | batch) + (1 | tray) + (1 | samplingdate))

    # Fit the Bayesian mixed model
    fit <- brm(myformula, data = temp, family = skew_normal(), 
               data2 = list(Amat = Amat, Dmat = Dmat), 
               chains = 4, threads = threading(4), cores = 4, 
               iter = 4000, warmup = 1000, thin = 2)

    # Extract variance components
    Va <- median(unlist(VarCorr(fit, summary = FALSE)$animal))
    Vd <- median(unlist(VarCorr(fit, summary = FALSE)$dom))
    Vr <- median(unlist(VarCorr(fit, summary = FALSE)$residual))
    Vm <- median(unlist(VarCorr(fit, summary = FALSE)$dam))
    Vb <- median(unlist(VarCorr(fit, summary = FALSE)$batch))
    Vt <- median(unlist(VarCorr(fit, summary = FALSE)$tray))
    Vs <- median(unlist(VarCorr(fit, summary = FALSE)$samplingdate))

    total_var <- Va + Vm + Vr + Vd + Vb + Vt + Vs
    ha <- round(Va / total_var, 2)
    hd <- round(Vd / total_var, 2)

    outp <- c(as.character(gene), Va, Vd, Vm, Vr, Vb, Vt, Vs, ha, hd)
    out1[[iter]] <- outp

    if (individual_files) {
      write.table(data.frame(t(outp)), file = output_gene_file, sep = "\t", row.names = FALSE, col.names = TRUE)
    }
  }

  out1_df <- do.call(rbind, out1)
  write.table(out1_df, file = paste0(outdir, "/output_", st, "_", end, ".csv"), sep = "\t", row.names = FALSE, col.names = TRUE)
  return(out1_df)
}

# Run variance partitioning
variancesCalc(df_transposed_intra, st, end, batch_info, dam_info, tray_info, samplingdate_info, 
              outdir = output_dir, Amat = Amat, Dmat = Dmat, individual_files = TRUE)

message("Finished processing ", dataset_type, " dataset for genes ", st, " to ", end)