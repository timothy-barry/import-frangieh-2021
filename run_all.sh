# Note: This script should be run on an HPC with at least 250GB of RAM.
module load R/R-4.1.2
source ~/.research_config
source 1_download_data.sh
Rscript 2_convert_to_odm.R
Rscript 3_preliminary_QC.R