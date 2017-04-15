// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <Rcpp.h>
#include <string.h>
#include "pccc.h"

struct cccpar : public RcppParallel::Worker
{
  // input
  const RcppParallel::RMatrix<std::string> input;
  int version;

  // output
  RcppParallel::RMatrix<double> outmat;

  // initialization
  cccpar(const Rcpp::CharacterMatrix input, int version, Rcpp::NumericMatrix outmat) 
    : input(input), version(version), outmat(outmat) {}

  // the operator
  void operator()(std::size_t begin, std::size_t end) {
    ccc_codes cds(version);

    std::vector<std::string> these_codes(input.ncol());

    size_t i, j;


    for (size_t i = begin; i < end; ++i) {
      for(j = 0; j < input.ncol(); ++j) {
        these_codes[j] = input(i, j);
      }

      outmat(i,  0) = cds.neuromusc(these_codes);
      outmat(i,  1) = cds.cvd(these_codes);
      outmat(i,  2) = cds.respiratory(these_codes);
      outmat(i,  3) = cds.renal(these_codes);
      outmat(i,  4) = cds.gi(these_codes);
      outmat(i,  5) = cds.hemato_immu(these_codes);
      outmat(i,  6) = cds.metabolic(these_codes);
      outmat(i,  7) = cds.congeni_genetic(these_codes);
      outmat(i,  8) = cds.malignancy(these_codes);
      outmat(i,  9) = cds.neonatal(these_codes);
      outmat(i, 10) = cds.tech_dep(these_codes);
      outmat(i, 11) = cds.transplant(these_codes);

      if (outmat(i,  0) + outmat(i,  1) + outmat(i,  2) + outmat(i,  3) + 
          outmat(i,  4) + outmat(i,  5) + outmat(i,  6) + outmat(i,  7) +
          outmat(i,  8) + outmat(i,  9) + outmat(i, 10) + outmat(i, 11)) {
        outmat(i, 12) = 1;
      } else {
        outmat(i, 12) = 0;
      }
    } 
  }
};


// [[Rcpp::export]]
Rcpp::DataFrame ccc_rcpp_par(Rcpp::CharacterMatrix& MAT, int version) {

  Rcpp::NumericMatrix OUT(MAT.nrow(), 13);

  cccpar CCCPAR(MAT, version, OUT);

  RcppParallel::parallelFor(0, MAT.nrow(), CCCPAR);

  Rcpp::DataFrame OUTD = Rcpp::internal::convert_using_rfunction(OUT, "as.data.frame");

  OUTD.attr("names") = Rcpp::CharacterVector::create("neuromusc",
      "cvd", "respiratory", "renal", "gi", "hemato_immu", "metabolic",
      "congeni_genetic", "malignancy", "neonatal", "tech_dep", "transplant",
      "ccc_flag");

  return OUTD; 
}
