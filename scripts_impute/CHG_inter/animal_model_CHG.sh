#!/bin/bash
#SBATCH --job-name=brms_model        # Job name
#SBATCH --output=logs/output_%A_%a.log       # Standard output log for each array job
#SBATCH --error=logs/error_%A_%a.log         # Error log for each array job
#SBATCH --array=1-44              # Define array range (e.g., 1-100 for 100 tasks)
#SBATCH --ntasks=1                   # Run a single task per array job
#SBATCH --account=ag-demeaux
#SBATCH --cpus-per-task=8            # Number of CPU cores per task
#SBATCH --mem=64G                    # Memory per task
#SBATCH --time=78:00:00              # Time limit hrs:min:sec

# Load necessary modules

module load lang/R/4.2.1-foss-2022a


# Define variables for array job
total_genes=11000                     # Total number of genes
genes_per_job=250                     # Number of genes per job
start_gene=$(( (SLURM_ARRAY_TASK_ID - 1) * genes_per_job + 1 ))
end_gene=$(( SLURM_ARRAY_TASK_ID * genes_per_job ))

# Adjust the end_gene if it exceeds the total number of genes
if [ $end_gene -gt $total_genes ]; then
    end_gene=$total_genes
fi

echo "Processing genes from $start_gene to $end_gene for task ID $SLURM_ARRAY_TASK_ID"

# Run the R script with start and end gene as arguments
Rscript /projects/ag-demeaux/labraha3/EpiDom_A_lyrata/Animal_model/scripts_impute/CHG_inter/animal_model_CHG.R $start_gene $end_gene