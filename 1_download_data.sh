######################################################################
#
# Download raw data from Frangieh et al, 2021
#
# Files come from two locations: Broad's Single Cell Portal
# and supplementary tables from Nature Genetics. 
#
######################################################################

# get paths to local files
source ~/.research_config

# check if data directory already exists
if [ -d $LOCAL_FRANGIEH_2020_DATA_DIR ] 
then
  echo "Raw data directory already exists!"
  exit 1
else
  # make the data directory
  mkdir $LOCAL_FRANGIEH_2020_DATA_DIR
fi

# navigate to data directory
cd $LOCAL_FRANGIEH_2020_DATA_DIR

####### 1. Download files from single cell portal #########

# download the file with all the download links 
# NOTE: the URL contains a single-use authorization code that must be generated
# at https://singlecell.broadinstitute.org/single_cell/study/SCP1064/multi-modal-pooled-perturb-cite-seq-screens-in-patient-models-define-novel-mechanisms-of-cancer-immune-evasion#/
wget --no-check-certificate "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP1064&auth_code=oNv7WFN8&directory=all&context=study" -O cfg.txt

# download the files
echo "Downloading single cell portal files!"
curl -K cfg.txt

# remove the file cfg.txt
rm cfg.txt

# move the directory SCP1064 to single-cell-portal
mv SCP1064 single-cell-portal

# unzip files
echo "Unzipping single cell portal files!"
gunzip single-cell-portal/other/RNA_expression.csv.gz
gunzip single-cell-portal/other/raw_CITE_expression.csv.gz 

####### 1. Download supplementary tables from Nature Genetics #########

echo "Downloading supplementary tables!"
mkdir supp_tables
cd supp_tables
supp_tab_url="https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-021-00779-1/MediaObjects/41588_2021_779_MOESM3_ESM.xlsx"
wget supp_tab_url

echo "Done!"