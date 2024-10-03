#' Fill Missing Celltypes with Zero and Retain Sample Metadata
#'
#' This function fills missing cell type rates in the provided dataset by adding any missing
#' combinations of samples and cell types, setting their rate to 0. Additionally, it retains
#' and matches all other metadata related to the samples.
#'
#' @param data A dataframe that contains at least sample, celltype, and rate columns.
#' @param sample_col A string specifying the column name for sample identifiers. Default is "Sample".
#' @param celltype_col A string specifying the column name for cell types. Default is "Celltype".
#' @param rate_col A string specifying the column name for the rate values. Default is "rate".
#'
#' @return A dataframe where missing sample-celltype combinations have been filled in with
#'         a rate of 0, and all metadata related to the samples has been retained and matched.
#' @export
#'
cwas_fill_missing_celltype <- function(data, sample_col = "Sample", celltype_col = "Celltype", rate_col = "rate") {
  # 提取与Sample相关的元数据列（除去Celltype和rate）
  meta_cols <- setdiff(names(data), c(celltype_col, rate_col))
  meta_data <- data[, meta_cols] %>% distinct()  # 保证元数据不重复

  # 找到所有的 sample 和 celltype 组合
  complete_samples <- expand.grid(Sample = unique(data[[sample_col]]),
                                  Celltype = unique(data[[celltype_col]]))

  # 合并实际数据，缺失值补0
  df_filled <- merge(complete_samples, data[, c(sample_col, celltype_col, rate_col)],
                     by = c(sample_col, celltype_col), all.x = TRUE)
  df_filled[[rate_col]][is.na(df_filled[[rate_col]])] <- 0

  # 根据Sample匹配元数据信息，保留所有列
  df_filled <- df_filled %>% left_join(meta_data, by = sample_col)
  df_filled <- unique(df_filled)
  # 返回补全后的数据
  return(df_filled)
}


