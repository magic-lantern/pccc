#' Complex Chronic Conditions (CCC)
#'
#' Generate CCC and CCC subcategory flags and the number of categories.
#'
#' It is recommended that you view the codes defining the CCC via
#' \code{\link{get_codes}} and make sure that the ICD codes in your data set are
#' formatted in the same way.  The ICD codes used for CCC are character strings
#' must be formatted as follows:
#' \itemize{
#' \item *Do not* use decimal points or other separators
#' \item ICD 9 codes: Codes less than 10 should be left padded with 2 zeros. Codes
#' less than 100 should be left padded with 1 zero.
#' }
#'
#' See See `vignette("pccc-overview")` for more details.
#'
#' @references
#' See \code{\link{pccc-package}} for published paper on the topic of identifying
#' Complex Chronix Conditions
#'
#' @param data a \code{data.frame} containing a patient id and all the ICD-9-CM
#' or ICD-10-CM codes.  The \code{data.frame} passed to the function should be
#' in wide format.
#' @param id bare name of the column containing the patient id
#' @param dx_cols,pc_cols column names with the diagnostic codes and procedure
#' codes respectively.  These argument are passed to \code{\link[dplyr]{select}}.
#' @param icdv ICD version 9 or 10
#'
#' @seealso \code{\link{get_codes}} to view the ICD codes used to define the
#' CCC.  \code{\link[dplyr]{select}} for more examples and details on how to
#' identify and select the diagnostic and procedure code columns.
#'
#' @return A \code{data.frame} with a column for the subject id and integer (0
#' or 1) columns for each each of the categories.
#'
#' @example examples/ccc.R
#'
#' @export
ccc <- function(data, id, dx_cols = NULL, pc_cols = NULL, icdv) {
  UseMethod("ccc")
}

#' @method ccc data.frame
#' @export
ccc.data.frame <- function(data, id, dx_cols, pc_cols, icdv) {

  if (missing(dx_cols) & missing(pc_cols)) {
    stop("dx_cols and pc_cols are both missing.  At least one must not be.",
         call. = FALSE)
  }

  if (!missing(dx_cols)) {
    dxmat <- as.matrix(dplyr::mutate_all(dplyr::select(data, !!dplyr::enquo(dx_cols)), as.character))
  } else {
    dxmat <- matrix("", nrow = nrow(data))
  }

  if (!missing(pc_cols)) {
    pcmat <- as.matrix(dplyr::mutate_all(dplyr::select(data, !!dplyr::enquo(pc_cols)), as.character))
  } else {
    pcmat <- matrix("", nrow = nrow(data))
  }

  ids <- dplyr::select(data, !!dplyr::enquo(id))

  # # test with 100k rows of data
  # timing: r version: 103.123 sec elapsed - uses similar logic to C++ version
  # timing: r hash version: 17.759 sec elapsed - env to store ICD codes for comparison
  # timing: r flat hash version: 20.058 sec elapsed - first version of using env to store ICD codes for comparison
  # timing: c++ version: 6.593 sec elapsed - original
  # #
  library(tictoc)
  tic("timing: r version")
  dplyr::bind_cols(ids, ccc_mat_r(dxmat, pcmat, icdv))
  toc()

  tic("timing: r hash version")
  dplyr::bind_cols(ids, ccc_hash_r(dxmat, pcmat, icdv))
  toc()

  tic("timing: r flat hash version")
  dplyr::bind_cols(ids, ccc_hash_flat_r(dxmat, pcmat, icdv))
  toc()

  tic("timing: c++ version")
  dplyr::bind_cols(ids, ccc_mat_rcpp(dxmat, pcmat, icdv))
  toc()
}

ccc_mat_r <- function(dx, pc, version = 9L) {
  set_code_set(version)

  dx_neuromusc <- c("3180","3181","3182","330","33111","33119","3314","33189","3319","3320","3321",
                    "3330","3332","3334","3335","3337","3339","334","335","343","34501","34581","3590","3591",
                    "3592","3593","3361","3368","3379","3418","34290","3440","34481","3449","34511",
                    "3453","34541","34561","34571","34591","3481","3484","3491","43401","43491","359","740",
                    "741","742","7595","78003","9962","99663","V452","V5301","V5302")

  dx_cvd <- c("4161","4168","4169","4240","4258","4242","4243","4250","4251","4252","4253",
                      "4254","426","4270","4271","4272","4273","4274","4276","4277","4278","4279","4280","42883",
                      "4291","4293","43311","7450","7451","7452","7453","7456","746","7471","7472","7473","7474",
                      "74781","74789","9960","9961","99661","99662","V421","V422","V432","V433","V450","V4581",
                      "V533")

  dx_respiratory <- c("32725","4162","5163","51631","51884","5190","2770","748","7704","V426",
                              "V440","V4576","V460","V461","V550")

  dx_renal <- c("34461","585","5964","59653","59654","753","99668","V420","V445","V446",
                        "V451","V4573","V4574","V536","V555","V556","V56")

  dx_gi <- c("4530","5364","555","556","5571","5602","5647","5714","5715","5716","5717",
                     "5718","5719","7503","751","V427","V4283","V4284","V441","V442","V443","V444","V5350",
                     "V5351","V5359","V551","V552","V553","V554")

  dx_hemato_immu <- c("042","043","044","135","279","2820","2821","2822","2823","2824","2825",
                              "2826","2881","2882","284","2860","2863","28732","28733","28739","28801","28802",
                              "2884","4460","4461","44621","4464","4465","4466","4467","7100","7101","7103","V08")

  dx_metabolic <- c("243","2532","2535","2536","2539","2550","25513","2552","270","271","272",
                            "2750","2751","2752","2753","2772","2773","2774","2775","2776","2778","2779",
                            "V4585","V5391","V6546")

  dx_congeni_genetic <- c("2594","5533","7373","7560","7561","7562","7563","7564","7565","7566",
                                  "7567","758","7597","7598","7599")

  dx_malignancy <- c("14","15","16","17","18","19","20","23","V4281","V4282")

  dx_neonatal <- c("76401","76402","76411","76412","76421","76422","76491","76492","76501",
                           "76502","76511","76512","76521","76522","76523","7670","7674","7685","7687","7689",
                           "7702","7707","7710","7711","77213","77214","7733","7734","7747","7765","77753","7780",
                           "7797")

  dx_tech_dep <- c("5190","5364","9960","9961","9962","9964","99661","99662","99663","99666","99667",
                           "99668","9969","V433","V440","V441","V442","V443","V444","V445","V446","V450","V451",
                           "V452","V460","V461","V462","V5301","V5302","V5331","V5332","V5339","V5350","V5351",
                           "V5359","V536","V550","V551","V552","V553","V554","V555","V556")

  dx_transplant <- c("99680","99681","99682","99683","99684","99685","99686","99687","99688","99689",
                             "V421","V422","V426","V427","V4281","V4282","V4283","V4284","V432","V4585","V5391")

  pc_neuromusc <- c("0152","0153","0221","0222","0231","0232","0233","0234","0235","0239","0241",
                            "0242","0293","0371","0372","0379","0393","0397","0492")

  pc_cvd <- c("0050","0051","0053","0054","0055","0057","1751","1752","3581","3582","3583",
                      "3584","3741","3751","3752","3753","3754","3755","3760","3761","3763","3765","3766",
                      "3767","3768","3771","3772","3774","3776","3779","3780","3781","3782","3783","3785",
                      "3786","3787","3789","3794","3795","3796","3797","3798","3981","3982","3983","3984",
                      "3985","8945","8946","8947","8948","8949")

  pc_respiratory <- c("303","304","3121","3129","3141","3174","3241","3249","3250","3259","3321",
                              "3350","3351","3352","336", "3485","9655","9723")

  pc_renal <- c("3895","3927","3942","3993","3994","3995","5498","5502","5503","5504","5512",
                        "5551","5552","5553","5554","5561","5569","5593","5594","5597","5641","5642","5651","5652",
                        "5661","5662","5671","5672","5673","5674","5675","5679","5721","5722","5771","5779","5993",
                        "5994","8607","9645","9646","9647")

  pc_gi <- c("253","254","4210","4211","4242","4281","4311","4319","4391","4399","4412",
                     "4432","4438","4439","4563","4581","4582","4583","4613","4622","4623","4632","4640","4641",
                     "4643","4697","504", "5051","5059","526","527","5280","5282","5283","5284","5285","5286",
                     "5471","9624","9636","9702")

  pc_hemato_immu <- c("4100","4101","4102","4103","4104","4105","4106","4107","4108","4109","415",
                              "4194")

  pc_metabolic <- c("064","0652","0681","073","0764","0765","0768","0769","624","6241","645","6551",
                            "6553","6561","6563","6841","6849","6851","6859","6861","6869","6871","6879","8606")

  pc_malignancy <- c("0010","9925")

  pc_tech_dep <- c("0221","0222","0231","0232","0233","0234","0235","0239","0241","0242","0293","0371",
                           "0372","0379","0393","0397","0492","0050","0051","0053","0054","0055","0057","1751","1752",
                           "3741","3752","3753","3754","3755","3760","3761","3763","3765","3766","3767","3768","3771",
                           "3772","3774","3776","3779","3780","3781","3782","3783","3785","3786","3787","3789","3794",
                           "3795","3796","3797","3798","3981","3982","3983","3984","3985","8945","8946","8947","8948",
                           "8949","3121","3129","3141","3174","3321","3485","9655","9723","3895","3927","3942","3993",
                           "3994","3995","5498","5502","5503","5504","5512","5593","5594","5597","5651","5652","5661",
                           "5662","5672","5673","5674","5675","5721","5722","5993","5994","8607","9645","9646","9647",
                           "4210","4211","4281","4311","4319","4412","4432","4438","4439","4613","4622","4623","4632",
                           "4640","4641","4643","9624","9636","9702","8100","8101","8102","8103","8104","8105","8106",
                           "8107","8108","8109","8130","8131","8132","8133","8134","8135","8136","8137","8138","8139",
                           "8451")

  pc_transplant <- c("3751","3350","3351","3352","336", "5561","5569","4697","5051","5059","5280","5282",
                             "5283","5284","5285","5286","4100","4101","4102","4103","4104","4105","4106","4107","4108",
                             "4109","4194","0091","0092","0093")

  # using a matrix here instead of a data frame results in significant performance improvements
  out <- matrix(0L,
               nrow = nrow(dx),
               ncol = 13,
               dimnames = list(c(),                 # row names
                               c('neuromusc',       # column names - 1
                                 'cvd',             # 2
                                 'respiratory',     # 3
                                 'renal',           # 4
                                 'gi',              # 5
                                 'hemato_immu',     # 6
                                 'metabolic',       # 7
                                 'congeni_genetic', # 8
                                 'malignancy',      # 9
                                 'neonatal',        # 10
                                 'tech_dep',        # 11
                                 'transplant',      # 12
                                 'ccc_flag')))      # 13

  for(i in 1:nrow(dx)) {
    dx_row <- dx[i,]
    pc_row <- pc[i,]
    ccc_flag <- 0L

    if(find_match(dx_row, pc_row, dx_neuromusc, pc_neuromusc)) {
      out[i, 1] <- 1L
      ccc_flag <- 1L
    }
    else if(find_match(dx_row, pc_row, dx_cvd, pc_cvd)) {
      out[i,2] <- 1L
      ccc_flag <- 1L
    }
    else if(find_match(dx_row, pc_row, dx_respiratory, pc_respiratory)) {
      out[i, 3] <- 1L
      ccc_flag <- 1L
    }
    else if(find_match(dx_row, pc_row, dx_renal, pc_renal)) {
      out[i, 4] <- 1L
      ccc_flag <- 1L
    }
    else if(find_match(dx_row, pc_row, dx_gi, pc_gi)) {
      out[i, 5] <- 1L
      ccc_flag <- 1L
    }
    else if (find_match(dx_row, pc_row, dx_hemato_immu, pc_hemato_immu)) {
      out[i, 6] <- 1L
      ccc_flag <- 1L
    }
    else if (find_match(dx_row, pc_row, dx_metabolic, pc_metabolic)) {
      out[i, 7] <- 1L
      ccc_flag <- 1L
    }
    else if (find_match(dx_row, NA, dx_congeni_genetic)) {
      out[i, 8] <- 1L
      ccc_flag <- 1L
    }
    else if (find_match(dx_row, pc_row, dx_malignancy, pc_malignancy)) {
      out[i, 9] <- 1L
      ccc_flag <- 1L
    }
    else if (find_match(dx_row, NA, dx_neonatal)) {
      out[i, 10] <- 1L
      ccc_flag <- 1L
    }

    if (find_match(dx_row, pc_row, dx_tech_dep, pc_tech_dep)) {
      out[i, 11] <- 1L
      ccc_flag <- 1L
    }

    if (find_match(dx_row, pc_row, dx_transplant, pc_transplant)) {
      out[i, 12] <- 1L
      ccc_flag <- 1L
    }

    # Through various tests - recorded at https://gist.github.com/magic-lantern/648afc6963b1bc89dd1a6efc74eac2c0
    # this method was found to be fastest in general address cells in a matrix by number/integer
    # rather than name; also logic based on a flag is quicker than doing a row sum
    if(ccc_flag)
      out[i, 13] <- 1L
  }

  as.data.frame(out)
}

# other than this string comparison function, the pure R version is identical in performance
# to the C++ version. Just need to fix this portion...
find_match <- function(dx,
                       pc,
                       dx_codes,
                       pc_codes = NULL){

  for (c in dx_codes) {
    # imported stringi in NAMESPACE - faster than calling stringi:: stri_startswith...
    if(any(stri_startswith_fixed(dx, c),na.rm = TRUE))
      return(1L)
  }

  for (c in pc_codes) {
    # imported stringi in NAMESPACE - faster than calling stringi:: stri_startswith...
    if(any(stri_startswith_fixed(pc, c),na.rm = TRUE))
      return(1L)
  }
  #return 0 if no match
  0L
}

# using new list of envs
ccc_hash_r <- function(dx, pc, version = 9L) {
  out <- matrix(0L,
                nrow = nrow(dx),
                ncol = 13,
                dimnames = list(c(),                 # row names
                                c('neuromusc',       # column names - 1
                                  'cvd',             # 2
                                  'respiratory',     # 3
                                  'renal',           # 4
                                  'gi',              # 5
                                  'hemato_immu',     # 6
                                  'metabolic',       # 7
                                  'congeni_genetic', # 8
                                  'malignancy',      # 9
                                  'neonatal',        # 10
                                  'tech_dep',        # 11
                                  'transplant',      # 12
                                  'ccc_flag')))      # 13

  primary_codes <- get_primary_codes(version)
  tech_codes <- get_codes_subset(icdv = version,  ccc_subset_name = 'tech_dep')
  transplant_codes <- get_codes_subset(icdv = version, ccc_subset_name = 'transplant')

  for(i in 1:nrow(dx)) {
    dx_row <- dx[i, ]
    pc_row <- pc[i, ]

    for (l in pkg.env[["min_length"]]:pkg.env[["max_length"]]) {
      # tested several methods to trim the string to the correct length, and the base substr()
      # is fastest
      dx_trimmed <- substr(dx_row, 1, l)
      pc_trimmed <- substr(pc_row, 1, l)

      p_dx <- primary_codes[[l]][['dx']]
      tc_dx <- tech_codes[[l]][['dx']]
      tp_dx <- transplant_codes[[l]][['dx']]

      p_pc <- primary_codes[[l]][['pc']]
      tc_pc <- tech_codes[[l]][['pc']]
      tp_pc <- transplant_codes[[l]][['pc']]

      # check diagnosis codes for mutually exclusive cccs
      for (c in dx_trimmed) {
        match <- p_dx[[c]]
        if (!is.null(match)) {
          out[i, match] <- 1L
          out[i, 13] <- 1L
        }

        # look at tech_dependency codes
        match <- tc_dx[[c]]
        if (!is.null(match)) {
          out[i, match] <- 1L
          out[i, 13] <- 1L
        }

        # look at transplant related codes
        match <- tp_dx[[c]]
        if (!is.null(match)) {
          out[i, match] <- 1L
          out[i, 13] <- 1L
        }
      }

      # # check procedure codes
      for (c in pc_trimmed) {
        match <- p_pc[[c]]
        if (!is.null(match)) {
          out[i, match] <- 1L
          out[i, 13] <- 1L
        }

        # look at tech_dependency codes
        match <- tc_pc[[c]]
        if (!is.null(match)) {
          out[i, match] <- 1L
          out[i, 13] <- 1L
        }

        # look at transplant related codes
        match <- tp_pc[[c]]
        if (!is.null(match)) {
          out[i, match] <- 1L
          out[i, 13] <- 1L
        }
      }
    }
  }

  as.data.frame(out)
}

# using new list of envs
ccc_hash_flat_r <- function(dx, pc, version = 9L) {
  out <- matrix(0L,
                nrow = nrow(dx),
                ncol = 13,
                dimnames = list(c(),                 # row names
                                c('neuromusc',       # column names - 1
                                  'cvd',             # 2
                                  'respiratory',     # 3
                                  'renal',           # 4
                                  'gi',              # 5
                                  'hemato_immu',     # 6
                                  'metabolic',       # 7
                                  'congeni_genetic', # 8
                                  'malignancy',      # 9
                                  'neonatal',        # 10
                                  'tech_dep',        # 11
                                  'transplant',      # 12
                                  'ccc_flag')))      # 13

  primary_codes <- get_collapsed_codes(version)
  tech_codes <- get_collapsed_codes_subset(icdv = version,  ccc_subset_name = 'tech_dep')
  transplant_codes <- get_collapsed_codes_subset(icdv = version, ccc_subset_name = 'transplant')


  for(i in 1:nrow(dx)) {
    pt_codes <- c(dx[i, ], pc[i, ])

    # look at the 'mutually exclusive' cccs
    for (l in pkg.env[["min_length"]]:pkg.env[["max_length"]]) {
      trimmed <- substr(pt_codes, 1, l)
      for (c in trimmed) {
        if (!is.null(primary_codes[[l]][[c]])) {
          out[i, primary_codes[[l]][[c]]] <- 1L
          out[i, 13] <- 1L
        }
      }
    }

    # look at tech_dependency codes
    for (l in pkg.env[["min_length"]]:pkg.env[["max_length"]]) {
      trimmed <- substr(pt_codes, 1, l)
      for (c in trimmed) {
        if (!is.null(tech_codes[[l]][[c]])) {
          out[i, tech_codes[[l]][[c]]] <- 1L
          out[i, 13] <- 1L
        }
      }
    }

    # look at transplant related codes
    for (l in pkg.env[["min_length"]]:pkg.env[["max_length"]]) {
      trimmed <- substr(pt_codes, 1, l)
      for (c in trimmed) {
        if (!is.null(transplant_codes[[l]][[c]])) {
          out[i, transplant_codes[[l]][[c]]] <- 1L
          out[i, 13] <- 1L
        }
      }
    }
  }

  as.data.frame(out)
}


