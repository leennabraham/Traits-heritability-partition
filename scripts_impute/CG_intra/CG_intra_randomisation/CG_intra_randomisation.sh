#!/bin/bash
#SBATCH --job-name=brms_full        # Job name for full dataset run
#SBATCH --output=logs/output_%A_%a.log  # Standard output log for each array job
#SBATCH --error=logs/error_%A_%a.log    # Error log for each array job
#SBATCH --array=1-82                  # Array range (20,496 genes / 250 per job ≈ 82 jobs)
#SBATCH --ntasks=1                     # Run a single task per array job
#SBATCH --account=ag-demeaux
#SBATCH --cpus-per-task=8              # Number of CPU cores per task
#SBATCH --mem=64G                      # Memory per task
#SBATCH --time=72:00:00                # Adjusted time for full dataset

# Load necessary modules
module load lang/R/4.2.1-foss-2022a

# Define variables for full run
total_genes=20496                     # Running on all 20,496 genes
genes_per_job=250                      # Number of genes per job
start_gene=$(( (SLURM_ARRAY_TASK_ID - 1) * genes_per_job + 1 ))
end_gene=$(( SLURM_ARRAY_TASK_ID * genes_per_job ))

# Adjust the end_gene if it exceeds total_genes
if [ $end_gene -gt $total_genes ]; then
    end_gene=$total_genes
fi

echo "Processing genes from $start_gene to $end_gene for task ID $SLURM_ARRAY_TASK_ID"

# Define script path
R_SCRIPT="/projects/ag-demeaux/labraha3/EpiDom_A_lyrata/Animal_model/scripts_impute/CG_intra/CG_intra_randomisation/CG_intra_randomised.R"

# Run observed data analysis
echo "Running observed data analysis..."
Rscript $R_SCRIPT $start_gene $end_gene "observed"

# Run permuted data analysis
echo "Running permuted data analysis..."
Rscript $R_SCRIPT $start_gene $end_gene "permuted"