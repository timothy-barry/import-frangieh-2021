###############################################################################
#
# Import Frangieh et al (2020) data.
#
#
# Notes:
###############################################################################

### retrieve top-level data directory ###
frangieh_dir <- .get_config_path("LOCAL_FRANGIEH_2021_DATA_DIR")

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


for (dir in c(processed_gene_dir, processed_gRNA_dir, processed_protein_dir)) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}


### import gene data ###
cat("Reading gene expression matrix from file...\n")
gene_expr_filename <- sprintf("%sraw/single-cell-portal/other/RNA_expression.csv", frangieh_dir)
gene_expr_metadata_filename <- sprintf("%sraw/single-cell-portal/metadata/RNA_metadata.csv", frangieh_dir)
gene_expr_cluster_filename <- sprintf("%sraw/single-cell-portal/cluster/5fd0e449771a5b0db7207711/RNA_UMAP_cluster.csv", frangieh_dir)

# read gene expression data
gene_expr_data <- readr::read_csv(gene_expr_filename,
  n_max = 3,
  col_types = list(
    GENE = readr::col_character(),
    .default = readr::col_double()
  )
)
# read gene expression metadata
gene_expr_metadata <- readr::read_csv(gene_expr_metadata_filename,
  col_types = list(
    MOI = readr::col_integer(),
    UMI_count = readr::col_double(),
    .default = readr::col_character()
  ),
  comment = "TYPE"
)
# read gene expression clusters
gene_expr_clusters <- readr::read_csv(gene_expr_cluster_filename,
  col_types = list(
    NAME = readr::col_character(),
    X = readr::col_double(),
    Y = readr::col_double()
  ),
  comment = "TYPE"
)

# extract cell barcodes and gene names from gene expression data
cell_barcodes_gene <- colnames(gene_expr_data)[-1]
gene_names <- gene_expr_data |>
  dplyr::select(GENE) |>
  dplyr::rename(gene_name = GENE)

# convert to matrix
cat("Converting gene expression data to matrix...\n")
gene_expr_data_norm <- gene_expr_data |>
  dplyr::select(-GENE) |>
  as.matrix()

# remove gene_expr_data from workspace to save memory
rm(gene_expr_data)

# undo normalization
cat("Undoing normalization...\n")
gene_expr_data_unnorm <- gene_expr_data_norm |>
  exp() |>
  sweep(2, gene_expr_metadata$UMI_count / 1e6, "*") |>
  round()

# save to ondisc matrix
cat("Creating ODM for gene expression matrix...\n")
odm_fp <- sprintf("%s/gene_expression_matrix.odm", processed_gene_dir)
metadata_fp <- sprintf("%s/gene_expression_metadata.rds", processed_gene_dir)
ondisc::create_ondisc_matrix_from_R_matrix(
  r_matrix = gene_expr_data_raw,
  barcodes = cell_barcodes_gene,
  features_df = gene_names,
  odm_fp = odm_fp,
  metadata_fp = metadata_fp
) |>
  # add condition information to cell covariates
  ondisc::mutate_cell_covariates(
    condition = gene_expr_metadata$condition[match(cell_barcodes_gene, gene_expr_metadata$NAME)],
    cluster_x = gene_expr_clusters$X[match(cell_barcodes_gene, gene_expr_clusters$NAME)],
    cluster_y = gene_expr_clusters$Y[match(cell_barcodes_gene, gene_expr_clusters$NAME)]
  ) |>
  # save to disk
  ondisc::save_odm(metadata_fp = metadata_fp)

### import protein data ###

cat("Reading protein expression matrix from file...\n")
prot_expr_filename <- sprintf("%sraw/single-cell-portal/other/raw_CITE_expression.csv", frangieh_dir)
prot_expr_data <- readr::read_csv(prot_expr_filename,
  n_max = 3,
  col_types = list(
    `...1` = readr::col_character(),
    .default = readr::col_integer()
  )
) |>
  dplyr::rename(Protein = `...1`)

# extract cell barcodes and protein names
cell_barcodes_protein <- colnames(prot_expr_data)[-1]
protein_names <- prot_expr_data |>
  dplyr::select(Protein) |>
  dplyr::rename(protein_name = Protein)

# convert to matrix
cat("Converting protein expression data to matrix...\n")
prot_expr_data <- prot_expr_data |>
  dplyr::select(-Protein) |>
  as.matrix()

# save to ondisc matrix
cat("Creating ODM for protein expression matrix...\n")
odm_fp <- sprintf("%s/protein_expression_matrix.odm", processed_protein_dir)
metadata_fp <- sprintf("%s/protein_expression_metadata.rds", processed_protein_dir)
ondisc::create_ondisc_matrix_from_R_matrix(
  r_matrix = prot_expr_data,
  barcodes = cell_barcodes_protein,
  features_df = protein_names,
  odm_fp = odm_fp,
  metadata_fp = metadata_fp
) |>
  # add condition information to cell covariates
  ondisc::mutate_cell_covariates(condition = gene_expr_metadata$condition[match(cell_barcodes_protein, gene_expr_metadata$NAME)]) |>
  # save to disk
  ondisc::save_odm(metadata_fp = metadata_fp)

### import gRNA data ###

cat("Reading gRNA assignments from file...\n")
gRNA_assignments_filename <- sprintf("%sraw/single-cell-portal/documentation/all_sgRNA_assignments.txt", frangieh_dir)
gRNA_list_filename <- sprintf("%sraw/supp_tables/41588_2021_779_MOESM3_ESM.xlsx", frangieh_dir)


gRNA_assignments <- readr::read_csv(gRNA_assignments_filename)
gRNA_list <- readxl::read_excel(gRNA_list_filename,
  sheet = 1,  # first sheet corresponds to gRNA list
  skip = 2,   # two first lines are header
  n_max = 818 # only the first 818 gRNAs used for perturb-CITE-seq
) 

# convert the gRNA assignments to indices that can be passed to sparseMatrix()
gRNA_mapping <- gRNA_assignments |>
  dplyr::mutate(cell_idx = dplyr::row_number()) |>
  dplyr::rowwise() |>
  dplyr::mutate(sgRNAs = ifelse(is.na(sgRNAs), NA, strsplit(sgRNAs, split = ","))) |>
  dplyr::mutate(sgRNA_idx = ifelse(any(is.na(sgRNAs)),
    list(),
    list(match(sgRNAs, gRNA_list$`Guide Name`))
  ), ) |>
  dplyr::select(-sgRNAs) |>
  tidyr::unnest(sgRNA_idx)

# create a sparse matrix of gRNA assignments
gRNA_assignment_matrix <- Matrix::sparseMatrix(
  i = gRNA_mapping$sgRNA_idx,
  j = gRNA_mapping$cell_idx,
  # specifying x is unnecessary because this is a
  # sparse logical matrix, but ondisc does not currently
  # support sparse logical matrix input
  x = nrow(gRNA_mapping),
  dims = c(nrow(gRNA_list), nrow(gRNA_assignments))
)

# extract the experimental design information
experimental_design <- gRNA_list |>
  dplyr::rowwise() |>
  dplyr::mutate(
    # those gRNAs with "SITE" in their names are non-targeting
    target_type = ifelse(grepl("SITE", `Guide Name`),
      "non-targeting",
      "gene"
    ),
    target = ifelse(target_type == "non-targeting",
      "non-targeting",
      strsplit(`Guide Name`, split = "_")[[1]][1]
    )
  ) |>
  dplyr::ungroup()

cat("Creating ODM for gRNA expression matrix...\n")
cell_barcodes <- gRNA_assignments$Cell
gRNA_barcodes <- gRNA_list |>
  dplyr::rename(
    gRNA_barcode = `sgRNA Sequence`,
    gRNA_name = `Guide Name`
  ) |>
  dplyr::select(gRNA_barcode, gRNA_name)
odm_fp <- sprintf("%s/gRNA_assignments_ungrouped.odm", processed_gRNA_dir)
metadata_fp <- sprintf("%s/gRNA_assignments_ungrouped_metadata.rds", processed_gRNA_dir)
ondisc::create_ondisc_matrix_from_R_matrix(
  r_matrix = gRNA_assignment_matrix,
  barcodes = cell_barcodes,
  features_df = gRNA_barcodes,
  odm_fp = odm_fp,
  metadata_fp = metadata_fp
) |>
  # add gRNA metadata to feature covariates
  ondisc::mutate_feature_covariates(
    target = experimental_design$target,
    target_type = experimental_design$target_type
  ) |>
  # save to disk
  ondisc::save_odm(metadata_fp = metadata_fp)

cat("Done.\n")
