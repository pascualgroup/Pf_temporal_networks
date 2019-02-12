#!/bin/bash
#SBATCH --job-name=XXX
#SBATCH --time=XXX
#SBATCH --output=slurm_output/XXX_%A_%a.out
#SBATCH --error=slurm_output/XXX_%A_%a.err
#SBATCH --array=XXX
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=XXX
#SBATCH --partition=broadwl
##SBATCH --mail-type=END
##SBATCH --mail-user=pilosofs@uchicago.edu

# Load modules
module load R
module load gcc/6.1

# Initialize
PS='XXX' # the parameter space
scenario='XXX' # Can be: 'S', 'G' or 'N'
exp='XXX' # The experiment ID
layers='XXX'
cutoff_prob=XXX
base_folder='/scratch/midway2/pilosofs/PLOS_Biol/' # Where all the source files are and where the experiment folders will be created
modularity_exp='XXX'


# Set up
prog='get_data_Midway.R'
cd $base_folder # Go to base folder
if test $modularity_exp == '0';
then
  echo modularity_exp 0
  mkdir -p $base_folder'/Results/'$PS'_'$scenario
  cp -n Infomap_v01926 $base_folder'/Results/'$PS'_'$scenario
  cp -n $prog $base_folder'/Results/'$PS'_'$scenario
  cp -n functions.R $base_folder'/Results/'$PS'_'$scenario
  cd $base_folder'/Results/'$PS'_'$scenario
fi

if test $modularity_exp != '0';
then
  echo modularity_exp not 0
  mkdir -p $base_folder'/Results/'$PS'_'$scenario'_'$modularity_exp
  cp -n Infomap_v01926 $base_folder'/Results/'$PS'_'$scenario'_'$modularity_exp
  cp -n $prog $base_folder'/Results/'$PS'_'$scenario'_'$modularity_exp
  cp -n functions.R $base_folder'/Results/'$PS'_'$scenario'_'$modularity_exp
  cd $base_folder'/Results/'$PS'_'$scenario'_'$modularity_exp
fi
