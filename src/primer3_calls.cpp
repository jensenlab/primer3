
#include <Rcpp.h>
using namespace Rcpp;

#include "oligotm.h"
#include "thal.h"


// [[Rcpp::export]]
NumericVector call_oligo_tm(CharacterVector oligos,
                            double mv,
                            double dv,
                            double dntp,
                            double dna,
                            int tp,
                            int sc) {
  int n = oligos.size();
  NumericVector tm(n);
  for (int i=0; i<n; i++) {
    tm[i] = oligotm(oligos[i], dna, mv, dv, dntp, (tm_method_type) tp, (salt_correction_type) sc);
  }
  
  return tm;
}

// [[Rcpp::export]]
NumericVector call_seq_tm(CharacterVector oligos,
                            double mv,
                            double dv,
                            double dntp,
                            double dna,
                            int nn_max_len,
                            int tp,
                            int sc) {
  int n = oligos.size();
  NumericVector tm(n);
  for (int i=0; i<n; i++) {
    tm[i] = seqtm(oligos[i], dna, mv, dv, dntp, nn_max_len, (tm_method_type) tp, (salt_correction_type) sc);
  }
  
  return tm;
}

// [[Rcpp::export]]
int is_thal_init() {
  return config_loaded;
}

// [[Rcpp::export]]
int call_thal_init(CharacterVector config_path) {
  thal_results output;
  int error = get_thermodynamic_values(config_path[0], &output);
  if (error) {
    fprintf(stderr, "%s\n", output.msg);
  };
  return error;
}

// [[Rcpp::export]]
void call_thal_free() {
  destroy_thal_structures();
}


// [[Rcpp::export]]
List call_thal(CharacterVector oligo1,
               CharacterVector oligo2,
               int debug,
               int alignment_type,
               int maxloop,
               double mv,
               double dv,
               double dntp,
               double dna,
               double temp,
               int temp_only,
               int dimer,
               int print_output) {
  thal_args targs = {
    debug,
    (thal_alignment_type) alignment_type,
    maxloop,
    mv,
    dv,
    dntp,
    dna,
    temp,
    temp_only,
    dimer
  };
  
  thal_results results;
  double NULL_REAL = 123456789.123456789;
  
  int n = oligo1.size();
  LogicalVector structure_found(n);
  NumericVector temps(n);
  NumericVector ds(n);
  NumericVector dh(n);
  NumericVector dg(n);
  IntegerVector align_end_1(n);
  IntegerVector align_end_2(n);
  
  for (int i=0; i<n; i++) {
    results.temp = NULL_REAL;
    results.ds = NULL_REAL;
    results.dh = NULL_REAL;
    results.dg = NULL_REAL;
    thal(oligo1[i], oligo2[i], &targs, &results, print_output);
    structure_found[i] = results.no_structure != 1;
    temps[i] = results.temp != NULL_REAL ? results.temp : NA_REAL;
    ds[i] = results.ds != NULL_REAL ? results.ds : NA_REAL;
    dh[i] = results.dh != NULL_REAL ? results.dh : NA_REAL;
    dg[i] = results.dg != NULL_REAL ? results.dg : NA_REAL;
    align_end_1[i] = results.align_end_1;
    align_end_2[i] = results.align_end_2;
  }
  
  
  return List::create(Named("structure_found", structure_found),
                      Named("temp", temps),
                      Named("ds", ds),
                      Named("dh", dh),
                      Named("dg", dg),
                      Named("align_end_1", align_end_1),
                      Named("align_end_2", align_end_2));
}
