// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

#include <string>
#include <Rcpp.h>
#include "pccc.h"

// [[Rcpp::export]]
Rcpp::DataFrame ccc_mat_rcpp(Rcpp::CharacterMatrix& dx, Rcpp::CharacterMatrix& pc, int version = 9)
{
  codes cdv(version);

  Rcpp::CharacterVector ccc_mat_rcpp_col_names(codes::col_names);
  ccc_mat_rcpp_col_names.push_back("ccc_flag");

  Rcpp::IntegerMatrix outmat(dx.nrow(), 13);

  Rcpp::CharacterVector dx_row;
  Rcpp::CharacterVector pc_row;
  std::vector<std::string> dx_str;
  std::vector<std::string> pc_str;

  for (size_t i=0; i < dx.nrow(); ++i) {
    dx_row = dx.row(i);
    pc_row = pc.row(i);
    dx_str = Rcpp::as<std::vector<std::string>>(dx_row);
    pc_str = Rcpp::as<std::vector<std::string>>(pc_row);

    if(cdv.neuromusc(dx_str, pc_str))
      outmat(i,  0) = 1;
    else if(cdv.cvd(dx_str, pc_str))
      outmat(i,  1) = 1;
    else if(cdv.respiratory(dx_str, pc_str))
      outmat(i,  2) = 1;
    else if(cdv.renal(dx_str, pc_str))
      outmat(i,  3) = 1;
    else if(cdv.gi(dx_str, pc_str))
      outmat(i,  4) = 1;
    else if (cdv.hemato_immu(dx_str, pc_str))
      outmat(i,  5) = 1;
    else if (cdv.metabolic(dx_str, pc_str))
      outmat(i,  6) = 1;
    else if (cdv.congeni_genetic(dx_str))
      outmat(i,  7) = 1;
    else if (cdv.malignancy(dx_str, pc_str))
      outmat(i,  8) = 1;
    else if (cdv.neonatal(dx_str))
      outmat(i,  9) = 1;

    outmat(i, 10) = cdv.tech_dep(dx_str, pc_str);
    outmat(i, 11) = cdv.transplant(dx_str, pc_str);

    if (sum(outmat.row(i))) {
      outmat(i, 12) = 1;
    }
    Rcpp::checkUserInterrupt();
  }

  outmat.attr("dimnames") = Rcpp::List::create(Rcpp::CharacterVector::create(),
              ccc_mat_rcpp_col_names
  );

//  return outmat;
 Rcpp::DataFrame out = Rcpp::internal::convert_using_rfunction(outmat, "as.data.frame");
 return out;

}
