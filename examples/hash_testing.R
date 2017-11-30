rm(list=ls())
gc()


l <- c("1", "22", "333", "4444", "55555", "666666")
Filter(function(x) stri_length(x) > 3, l)

save(death_small_data_example, file = 'test_compress.rdata', compress = "bzip2")
saveRDS(death_small_data_example, file = 'test.rds', compress = 'bzip2')
feather::write_feather(death_small_data_example, 'test.feather')
save(city, country, file="1.RData")

trim <- death_small_data_example[1:1000000, ]

############################################################################################
# setup
############################################################################################
rm(list=ls())
gc()
library(dplyr)
library(pccc)
dataset <- feather::read_feather('sample_icd9_data.feather')
patient <- dataset[1, ]
dx <- as.matrix(patient %>% dplyr:: select(starts_with("dx")))
dx[1, 'dx14'] <- '14'
dx[1, 'dx16'] <- '042'
dx[1, 'dx20'] <- '0050'
dx[1, 'dx25'] <- '25513'
pc <- as.matrix(patient %>% dplyr:: select(starts_with("pc")))



dx <- c('1234', '4321', '0000', '9999', '4444')
c <- '1'
library(microbenchmark)
library(stringi)
microbenchmark(
any(stringi::stri_startswith_fixed(dx, c),na.rm = TRUE),
any(stri_startswith_fixed(dx, c),na.rm = TRUE),
times = 100000)


###############################################################################
rm(list=ls())
gc()
library(pccc)
dataset <- feather::read_feather('sample_icd9_data.feather')
small <- dataset[1:100000, ]
rm(dataset)
gc()

###############################################################################
###############################################################################
###############################################################################
icd_ver <- 9
eg_data <-
  dplyr::data_frame(id = c('A'),
                    dx1 = c(sample(pccc::get_codes(icd_ver)[["hemato_immu", "dx"]], 1)),
                    dx2 = c(NA),
                    pc1 = c(sample(pccc::get_codes(icd_ver)[["cvd", "pc"]], 1))
  )
print(eg_data)
ccc(eg_data,
    id,
    dx_cols = dplyr::starts_with("dx"),
    pc_cols = dplyr::starts_with("pc"),
    icdv = icd_ver)

icd_ver <- 9
ccc(dataset,
      id,
      dx_cols = dplyr::starts_with("dx"),
      pc_cols = dplyr::starts_with("pc"),
      icdv = icd_ver)

###############################################################################
codes <- small[1, ]


library(microbenchmark)
microbenchmark(
  substr(codes, 1L, 2L),
  stri_sub(codes, from=1L, 2L),
  stringi::stri_sub(codes, from=1L, 2L),
  stri_sub(codes, from=1L, to=2L),
  stringr::str_sub(codes, start = 1L, end = 2L),
  times = 50000
)


isLongEnough <- function(X)
{
  ifelse(stri_length(X)>=3, TRUE, FALSE)
}
v<-c(1,2,3,4,5,5,5,5)
v [ isGoodNumber(v) == TRUE ]


microbenchmark(
nchar("string"),
stri_length("string"),
times = 10000)