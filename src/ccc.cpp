// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

#include <string>
#include <Rcpp.h>
#include "pccc.h"
#include <sstream>
#include <fstream>
#include <vector>
#include <array>
#include <chrono>

std::string trim(std::string const& str)
{
  // from https://codereview.stackexchange.com/questions/40124/trim-white-space-from-string
  
  if(str.empty())
    return str;
  
  std::size_t firstScan = str.find_first_not_of(' ');
  std::size_t first     = firstScan == std::string::npos ? str.length() : firstScan;
  std::size_t last      = str.find_last_not_of(' ');
  return str.substr(first, last-first+1);
}

// [[Rcpp::export]]
void test_ccc() {
  // Run this function with this line:
  // .Call('_pccc_test_ccc', PACKAGE = 'pccc')
  std::cout << " calling test_ccc\n";
  
  const int MAX_ROWS = 10;
  std::string line;
  
  auto start = std::chrono::system_clock::now();

  std::ifstream infile("C:/HCUPData/KID/KID_2009_Core.ASC");

  std::array<std::string, 121> row;
  std::vector<std::array<std::string, 121> > matrix;
  size_t mitr, ritr;
  
  Rcpp::CharacterMatrix dx(MAX_ROWS, 44);
  Rcpp::CharacterMatrix pc(MAX_ROWS, 44);
  Rcpp::CharacterVector dx_row(29);
  Rcpp::CharacterVector pc_row(15);
  
  int widths [2][121] = {
    {0,5,13,87,92,97,102,107,112,117,122,127,132,137,142,147,152,157,162,167,172,177,182,187,192,197,202,207,212,287,292,297,302,307,405,409,413,417,421,425,429,433,437,441,445,449,453,457,461,465},
    {4,12,86,91,96,101,106,111,116,121,126,131,136,141,146,151,156,161,166,171,176,181,186,191,196,201,206,211,286,291,296,301,306,404,408,412,416,420,424,428,432,436,440,444,448,452,456,460,464,600}
  };
  
  int row_count = 0;
  while (std::getline(infile, line) && row_count <= MAX_ROWS)
  {
    std::cout << row_count << "-----------------------\n" << line << "------------------------\n";
    for (ritr = 0; ritr < row.size(); ritr++) {
      row[ritr] = trim(line.substr(widths[0][ritr], widths[1][ritr] - widths[0][ritr] + 1));
    }
    
    //build the diagnostic code vector
    for (ritr = 3; ritr < 28; ritr++){
      dx[row_count, 0].push_back
      dx_row.push_back(row[ritr]);
    }
    for (ritr = 29; ritr < 33; ritr++){
      dx_row.push_back(row[ritr]);
    }
    
    // build the procedure code vector
    for (ritr = 34; ritr < 49; ritr++){
      pc_row.push_back(row[ritr]);
    }

    matrix.push_back(row);
    dx.
    dx[row_count] = dx_row;
    //pc.push_back(pc_row);
    
    row_count++;
  }
  
  std::cout << "matrix.size(): " << matrix.size() << "\n";
  
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
  std::cout << "time to load file: " << duration.count() << "ms\n";
  
  std::array<std::string, 121> r;
  for (mitr = 0; mitr < matrix.size(); ++mitr) {
    r = matrix[mitr];
    std::cout << "r[0]: ]" << r[0] << "[ ";
    std::cout << "r[1]: ]" << r[1] << "[ ";
    std::cout << "r[2]: ]" << r[2] << "[\n";
    for (ritr = 0; ritr < row.size(); ritr++) {
      //std::cout << r[ritr];
    }
  } 
  

  std::cout << " finsihed test_ccc\n";
}

// [[Rcpp::export]]
Rcpp::DataFrame ccc_mat_rcpp(Rcpp::CharacterMatrix& dx, Rcpp::CharacterMatrix& pc, int version = 9)
{ 
  codes cdv(version);

  Rcpp::IntegerMatrix outmat(dx.nrow(), 13);

  Rcpp::CharacterVector dx_row;
  Rcpp::CharacterVector pc_row;
  std::vector<std::string> dx_str;
  std::vector<std::string> pc_str;

  for (int i=0; i < dx.nrow(); ++i) { 

    dx_row = dx.row(i);
    pc_row = pc.row(i);
    dx_str = Rcpp::as<std::vector<std::string>>(dx_row);
    pc_str = Rcpp::as<std::vector<std::string>>(pc_row);
    outmat(i,  0) = cdv.neuromusc(dx_str, pc_str);
    outmat(i,  1) = cdv.cvd(dx_str, pc_str);
    outmat(i,  2) = cdv.respiratory(dx_str, pc_str);
    outmat(i,  3) = cdv.renal(dx_str, pc_str);
    outmat(i,  4) = cdv.gi(dx_str, pc_str);
    outmat(i,  5) = cdv.hemato_immu(dx_str, pc_str);
    outmat(i,  6) = cdv.metabolic(dx_str, pc_str);
    outmat(i,  7) = cdv.congeni_genetic(dx_str);
    outmat(i,  8) = cdv.malignancy(dx_str, pc_str);
    outmat(i,  9) = cdv.neonatal(dx_str);
    outmat(i, 10) = cdv.tech_dep(dx_str, pc_str);
    outmat(i, 11) = cdv.transplant(dx_str, pc_str);
    outmat(i, 12) = 0; 

    if (sum(outmat.row(i))) {
      outmat(i, 12) = 1; 
    } 
  }

  outmat.attr("dimnames") = Rcpp::List::create(Rcpp::CharacterVector::create(),
      Rcpp::CharacterVector::create("Neuromuscular", "CVD", "Respiratory", "Renal", "GI", "Hemato_immu", "Metabolic", "Congeni_genetic", "Malignancy", "Neonatal", "Tech_dep", "Transplant", "ccc_flag")
      );

//  return outmat; 
 Rcpp::DataFrame out = Rcpp::internal::convert_using_rfunction(outmat, "as.data.frame");
 return out;

}
