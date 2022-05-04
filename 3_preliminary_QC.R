###############################################################################
#
# Preliminary QC for Frangieh et al (2020) data.
#
# The goal of the preliminary QC is to (a) retain only cells with exactly one
# perturbation and (b) split the data into three separate ODMs based on
# experimental condition, since the three experimental conditions were analyzed
# separately.
#
###############################################################################

### retrieve top-level data directory ###
frangieh_dir <- .get_config_path("LOCAL_FRANGIEH_2021_DATA_DIR")

# get paths to processed directories for each modality
exper_name <- "perturb-cite-seq"
processed_gene_dir <- sprintf(
  "%sprocessed/%s/gene",
  frangieh_dir, exper_name
)
processed_gRNA_dir <- sprintf(
  "%sprocessed/%s/gRNA",
  frangieh_dir, exper_name
)
processed_protein_dir <- sprintf(
  "%sprocessed/%s/protein",
  frangieh_dir, exper_name
)

# read gene ODM
gene_odm_fp <- sprintf("%s/gene_expression_matrix.odm", processed_gene_dir)
gene_metadata_fp <- sprintf("%s/gene_expression_metadata.rds", processed_gene_dir)
gene_odm <- ondisc::read_odm(gene_odm_fp, gene_metadata_fp)

# read gRNA ODM
gRNA_odm_fp <- sprintf("%s/gRNA_assignments_ungrouped.odm", processed_gRNA_dir)
gRNA_metadata_fp <- sprintf("%s/gRNA_assignments_ungrouped_metadata.rds", processed_gRNA_dir)
gRNA_odm <- ondisc::read_odm(gRNA_odm_fp, gRNA_metadata_fp)

# read protein ODM
protein_odm_fp <- sprintf("%s/protein_expression_matrix.odm", processed_protein_dir)
protein_metadata_fp <- sprintf("%s/protein_expression_metadata.rds", processed_protein_dir)
protein_odm <- ondisc::read_odm(protein_odm_fp, protein_metadata_fp)

# extract info on experimental condition and MOI from gene and gRNA ODMs
condition_MOI_info <- dplyr::left_join(
  gene_odm |>
    ondisc::get_cell_covariates() |>
    tibble::rownames_to_column(var = "cell_barcode") |>
    dplyr::select(cell_barcode, condition),
  gRNA_odm |>
    ondisc::get_cell_covariates() |>
    tibble::rownames_to_column(var = "cell_barcode") |>
    dplyr::select(cell_barcode, n_nonzero),
  by = "cell_barcode"
)

# extract the three conditions
conditions <- condition_MOI_info |>
  dplyr::pull(condition) |>
  unique()

# for each condition, subset the gene, gRNA, and protein ODMs
for (condition in conditions) {
  # extract cells in the given condition, with exactly one gRNA
  cells_to_keep <- condition_MOI_info |>
    dplyr::filter(condition == !!condition, n_nonzero == 1) |>
    dplyr::pull(cell_barcode)

  # rename conditions to eliminate the nonstandard gamma character, and change
  # to lowercase
  condition_name <- gsub("γ", "-gamma", condition) |> tolower()
  
  # subset and save gene ODM
  gene_metadata_fp <- sprintf("%s/gene_expression_metadata_%s.rds", processed_gene_dir, condition_name)
  gene_odm[, cells_to_keep] |>
    ondisc::save_odm(metadata_fp = gene_metadata_fp)

  # subset and save gRNA ODM
  gRNA_metadata_fp <- sprintf("%s/gRNA_assignments_ungrouped_metadata_%s.rds", processed_gRNA_dir, condition_name)
  gRNA_odm[, cells_to_keep] |>
    ondisc::save_odm(metadata_fp = gRNA_metadata_fp)

  # subset and save protein ODM
  protein_metadata_fp <- sprintf("%s/protein_expression_metadata_%s.rds", processed_protein_dir, condition_name)
  protein_odm[, cells_to_keep] |>
    ondisc::save_odm(metadata_fp = protein_metadata_fp)
}
