rm(list=ls())
gc()

icd_ver <- 10

eg_data <-
  dplyr::data_frame(id = c('A', 'B', 'C'),
                    dx1 = c(NA, NA, sample(get_codes(icd_ver)[["hemato_immu", "dx"]], 1)),
                    dx2 = c("acode", sample(get_codes(icd_ver)[["gi", "dx"]], 2)),
                    pc1 = c("bcode", sample(get_codes(icd_ver)[["cvd", "pc"]], 1), 'ccode')
                    )

ccc(eg_data,
    id,
    dx_cols = dplyr::starts_with("dx"),
    pc_cols = dplyr::starts_with("pc"),
    icdv = icd_ver)

