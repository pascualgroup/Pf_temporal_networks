#!/bin/bash
#SBATCH --job-name=Res
#SBATCH --time=05:00:00
#SBATCH --output=slurm_output/Res_%A_%a.out
#SBATCH --error=slurm_output/Res_%A_%a.err
#SBATCH --array=4,5,6
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=32000
#SBATCH --partition=broadwl
##SBATCH --mail-type=END
##SBATCH --mail-user=pilosofs@uchicago.edu

# Load modules
module load R

# Initialize
exp='001' # The experiment ID
base_folder='/scratch/midway2/pilosofs/PLOS_Biol/' # Where all the source files are and where the experiment folders will be created

# Set up
cd $base_folder # Go to base folder

# Run R script
Rscript Results_wrapper_Midway_plosbiol.R 'S' '001'
Rscript Results_wrapper_Midway_plosbiol.R 'N' '001'
Rscript Results_wrapper_Midway_plosbiol.R 'G' '001'
