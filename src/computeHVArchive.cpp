#include <Rcpp.h>
#include <stdlib.h>
using namespace Rcpp;

bool dominates(const NumericVector& a, const NumericVector& b) {
  return
    is_true(all(b - a >= 0)) && 
    is_true(any(b - a > 0));
}

struct NumericVectorComparator {
  bool operator() (const NumericVector& a,  
                const NumericVector& b) { 
    return (a[0] < b[0]);
  }
};

double hv_contribution(const NumericVector& ind,
                       const NumericVector& prev,
                       const NumericVector& next) {
  NumericVector delta = ind - pmax(prev, next);
  return delta[0] * delta[1];
}

// [[Rcpp::export]]
NumericVector computeHVArchiveCPP(NumericMatrix archive, NumericVector ref_point) {
  int n_individuals = archive.ncol();
  
  NumericVector hv_values(archive.ncol());
  double current_hv = 0;
  
  NumericVector prev_ind;
  NumericVector next_ind;
  
  std::set<NumericVector, NumericVectorComparator> nondominated;
  
  for (int idx = 0; idx < n_individuals; ++idx) {
    NumericVector individual = archive.column(idx);
    
    if (dominates(individual, ref_point)) {
        std::vector<NumericVector> to_delete;
      
        auto it = nondominated.upper_bound(individual);
        
        bool is_dominated = (it != nondominated.begin()) &&
                            (dominates(*prev(it), individual) || is_true(all(*prev(it) - individual == 0)));
        
        if (!is_dominated) {
          NumericVector prev_ind;
          
          if (it == nondominated.begin()) {
            prev_ind = {(*it)[0], ref_point[1]};
          } else {
            prev_ind = *prev(it);
          }
          
          while (it != nondominated.end() &&
                 dominates(individual, (*it))) {
            to_delete.push_back((*it));
  
            if (next(it) == nondominated.end()) {
              next_ind = {ref_point[0], (*it)[1]};
            } else {
              next_ind = *next(it);
            }
  
            current_hv -= hv_contribution((*it), prev_ind, next_ind);
  
            ++it;
          }
  
          for (auto v : to_delete) {
            nondominated.erase(v);
          }
        
          nondominated.insert(individual);
          
          it = nondominated.find(individual);
          
          if (it == nondominated.begin()) {
            prev_ind = {(*it)[0], ref_point[1]};
          } else {
            prev_ind = *(prev(it));
          }
          
          if (next(it) == nondominated.end()) {
            next_ind = {ref_point[0], (*it)[1]};
          } else {
            next_ind = *(next(it));
          }
          
          current_hv += hv_contribution((*it), prev_ind, next_ind);
        }
    }
    
    // Rcout << nondominated.size() << std::endl;
    hv_values[idx] = current_hv;
  }
  
  return hv_values;
}
