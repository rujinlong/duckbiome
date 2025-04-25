load_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), here::here("pc028e1e2.sqlite"))
  return(con)
}




add_binID <- function(df_binID2alias, df_annotation) {
  df_annotation %>%
    dplyr::mutate(bin_alias = str_replace(prot_id, "_.*", "")) %>%
    dplyr::left_join(df_binID2alias, by = "bin_alias") %>%
    dplyr::select(-bin_alias) %>%
    # put the column bin_id in the first column
    dplyr::select(binID, everything())
}



read_coverm_mag <- function(fin_abundance) {
  df <- fread(fin_abundance) %>%
    tibble::column_to_rownames("Genome") %>%
    as.matrix()

  return(df)
}



extract_ancombc_rst <- function(tse) {
  res <- tse$res
  res_summary <- res %>%
    dplyr::select(taxon,
           lfc = lfc_groupZF,
           se = se_groupZF,
           p_val = p_groupZF,
           q_val = q_groupZF,
           diff_abundant = diff_groupZF)
}
