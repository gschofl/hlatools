#include <Rcpp.h>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>
using namespace Rcpp;

template <typename... T>
auto zip(const T&... containers) -> boost::iterator_range<boost::zip_iterator<decltype(boost::make_tuple(std::begin(containers)...))>>
{
  auto zip_begin = boost::make_zip_iterator(boost::make_tuple(std::begin(containers)...));
  auto zip_end = boost::make_zip_iterator(boost::make_tuple(std::end(containers)...));
  return boost::make_iterator_range(zip_begin, zip_end);
}

std::vector<std::string> string_split(const std::string& str, const std::string& sep = ":") {
  std::vector<std::string> tokens;
  // Find first non-separating character
  std::size_t first_pos = str.find_first_not_of(sep, 0);
  // Find first separator
  std::size_t sep_pos = str.find_first_of(sep, first_pos);
  while (sep_pos != std::string::npos || first_pos != std::string::npos) {
    // Push field
    tokens.push_back( str.substr(first_pos, sep_pos - first_pos) );
    // Update first position
    first_pos = str.find_first_not_of(sep, sep_pos);
    // Update separator position
    sep_pos = str.find_first_of(sep, first_pos);
  }
  return tokens;
}

std::vector<int> fields_to_ints(const std::vector<std::string>& fields) {
  std::vector<int> out(fields.size());
  for (int i = 0; i < fields.size(); i++) {
    try {
      out[i] = std::stoi(fields[i]);
    } catch(std::invalid_argument&) {
      // Catches NMDP codes and sorts them to the front.
      // NMDP codes themselves are sorted by their first letter only!
      int tmp = fields[i][0];
      out[i] = -100 + tmp;
    }
  }
  return out;
}

bool compare_hla_alleles(const std::string& a1, const std::string& a2) {
  // HLA allele fields
  std::vector<int> f1;
  std::vector<int> f2;

  // Separate HLA prefix and fields
  std::vector<std::string> v1 = string_split(a1, "*");

  if (v1.size() == 2 /* prefix exists */ ) {
    f1 = fields_to_ints( string_split(v1[1], ":") );
  } else {
    f1 = fields_to_ints( string_split(v1[0], ":") );
  }

  std::vector<std::string> v2 = string_split(a2, "*");
  if (v2.size() == 2 /* prefix exists */ ) {
    f2 = fields_to_ints( string_split(v2[1], ":") );
  } else {
    f2 = fields_to_ints( string_split(v2[0], ":") );
  }

  auto f1_it = f1.begin();
  auto f2_it = f2.begin();
  while (f1_it != f1.end() && f2_it != f2.end()) {
    if (*f1_it < *f2_it) {
      return true;
    }
    if (*f1_it > *f2_it) {
      return false;
    }
    f1_it++;
    f2_it++;
    // If a1 is at the last field but a2 goes on then a1 < a2
    if (f1_it == f1.end() && f2_it != f2.end()) {
      return true;
    }
    // If a1 goes on but a2 is at the last field then a1 > a2
    if (f1_it != f1.end() && f2_it == f2.end()) {
      return false;
    }
  }
  return false;
}

//' Sort HLA alleles by field
//'
//' Sort HLA alleles based on the numeric values of each successive field.
//' NMDP codes are sorted to the front, alphabetically based on their first letter.
//'
//' @param alleles A vector of HLA alleles either prefixed (e.g., \emph{"HLA-A*01:01:01:02"})
//' or without prefix (e.g.: \emph{"01:01:01:02"}).
//' @return A sorted vector of HLA alleles.
//' @export
// [[Rcpp::export]]
std::vector<std::string> hla_sort(std::vector<std::string> alleles) {
  std::sort(alleles.begin(), alleles.end(), compare_hla_alleles);
  return alleles;
}

//' Get the HLA allele code prefix
//'
//' @param alleles A vector of HLA allele codes (e.g.,  \emph{HLA-A*01:01:01:02}).
//' @return A vector of HLA allele code prefixes (e.g., \emph{"HLA-A"} or \emph{""} if
//' the allele code was not prefixed).
//' @export
// [[Rcpp::export]]
std::vector<std::string> hla_prefix(std::vector<std::string> a) {
  int i;
  std::vector<std::string> tmp;
  std::vector<std::string> out(a.size());
  for (i = 0; i < a.size(); i++) {
    tmp = string_split(a[i], "*");
    if (tmp.size() == 2 /* e.g. HLA-A*01:01:01 */ ) {
      out[i] = tmp[0];
    } else /* e.g. 01:01:01 */ {
      out[i] = "";
    }
  }

  return out;
}

//' Get the first field from a HLA allele code
//'
//' @param alleles A vector of HLA allele codes (e.g.,  \emph{HLA-A*01:01:01:02}).
//' @return A vector of first fields.
//' @export
// [[Rcpp::export]]
std::vector<std::string> hla_field1(std::vector<std::string> a) {
  int i;
  std::vector<std::string> tmp;
  std::vector<std::string> out(a.size());
  for (i = 0; i < a.size(); i++) {
    tmp = string_split(a[i], "*"); // Split off prefix if exists
    if (tmp.size() == 2 /* e.g. HLA-A*01:01:01 */ ) {
      out[i] = string_split(tmp[1], ":")[0];
    } else /* e.g. 01:01:01 */ {
      out[i] = string_split(tmp[0], ":")[0];
    }
  }

  return out;
}

//' Get the second field from a HLA allele code
//'
//' @param alleles A vector of HLA allele codes (e.g.,\emph{HLA-A*01:01:01:02}).
//' @return A vector of second fields or or \emph{""} if no second field exists).
//' @export
// [[Rcpp::export]]
std::vector<std::string> hla_field2(std::vector<std::string> a) {
  int i;
  std::vector<std::string> tmp;
  std::vector<std::string> tmp2;
  std::vector<std::string> out(a.size());
  for (i = 0; i < a.size(); i++) {
    tmp = string_split(a[i], "*");
    if (tmp.size() == 2 /* e.g. HLA-A*01:01:01 */ ) {
      tmp2 = string_split(tmp[1], ":");
      if (tmp2.size() > 1)
        out[i] = tmp2[1];
      else
        out[i] = "";
    } else /* e.g. 01:01:01 */ {
      tmp2 = string_split(tmp[0], ":");
      if (tmp2.size() > 1)
        out[i] = tmp2[1];
      else
        out[i] = "";
    }
  }

  return out;
}

// [[Rcpp::export]]
std::vector<std::string> hla_allele_to_genotype(std::vector<std::string> a1, std::vector<std::string> a2) {
  std::vector<std::string> genotype;
  for (auto x : zip(a1, a2)) {
    std::string a, b;
    boost::tie(a, b) = x;
    if (compare_hla_alleles(a, b)) {
      genotype.push_back( a + '/' + b );
    } else {
      genotype.push_back( b + '/' + a );
    }
  }
  return genotype;
}


