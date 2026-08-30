#include <RcppArmadillo.h>
#include <algorithm>
#include <map>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// distance 1
// Function to calculate Euclidean distance
// [[Rcpp::export]]
double euclideanDistance(
    const arma::vec& x,
    const arma::vec& y
) {
  return arma::norm(x - y, 2);
}

// distance 2
// Function to calculate manhattan distance
double manhattanDistance(
    const arma::vec& x,
    const arma::vec& y
) {
  return arma::sum(arma::abs(x - y));
}

// [[Rcpp::export]]
Rcpp::NumericVector rowMeansC(Rcpp::NumericMatrix mat) {
  int nRows = mat.nrow();
  int nCols = mat.ncol();
  
  // If the matrix has only one column, return the input as a numeric vector
  if (nCols == 1) {
    return mat(_, 0); // return the first column of the matrix
  }
  
  Rcpp::NumericVector means(nRows);
  
  for (int i = 0; i < nRows; i++) {
    double sum = 0;
    for (int j = 0; j < nCols; j++) {
      sum += mat(i, j);
    }
    means[i] = sum / nCols;
  }
  
  return means;
}



// [[Rcpp::export]]
List list_nn(const arma::mat& mat, int k = 10) {
  int n = mat.n_cols; 
  List neighbors(n); // List to store nearest neighbors for each column
  
  // Compute Euclidean distances between columns
  arma::mat dists(n, n, fill::zeros); // Matrix to store distances
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double dist = arma::norm(mat.col(i) - mat.col(j), 2); // Compute distance between columns i and j
      dists(i, j) = dist;
      dists(j, i) = dist; // Distance matrix is symmetric
    }
  }
  
  // Find and sort nearest neighbors for each column
  for (int i = 0; i < n; ++i) {
    uvec indices = sort_index(dists.row(i).t()); // Use Armadillo to sort distances and obtain sorted indices
    IntegerVector idx = wrap(indices); // Convert to Rcpp::IntegerVector and adjust for 1-based indexing
    idx = idx + 1; // Adjust for 1-based indexing
    idx = idx[idx != i + 1];// Remove the column itself from its neighbors
    
    // Only keep the top 10 closest neighbors
    if (idx.size() > k) {
      idx = idx[Range(0, k-1)]; // Keep only the first 10 indices (closest neighbors)
    }
    
    neighbors[i] = idx;
  }
  
  return neighbors;
}


// [[Rcpp::export]]
arma::mat dist_mat(const arma::mat& mat) {
  int n = mat.n_cols; 
  
  // Compute Euclidean distances between columns
  arma::mat dists(n, n, fill::zeros); // Matrix to store distances
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double dist = norm(mat.col(i) - mat.col(j), 2); // Compute distance between columns i and j
      dists(i, j) = dist;
      dists(j, i) = dist; // Distance matrix is symmetric
    }
  }
  return dists;
}


// [[Rcpp::export]]
double SilhouetteDistancecpp(const Rcpp::NumericVector& ocp, const Rcpp::NumericMatrix& distmat) {
  int n = ocp.size();
  Rcpp::NumericVector a(n, 0.0);
  Rcpp::NumericVector b(n, 0.0);
  Rcpp::NumericVector sw(n, 0.0);
  
  for (int i = 0; i < n; ++i) {
    int cluster = ocp[i];
    std::vector<int> ai_ind;
    
    // Get indices where ocp equals the current cluster
    for (int j = 0; j < n; ++j) {
      if (ocp[j] == cluster) {
        ai_ind.push_back(j);
      }
    }
    
    // Compute 'a' value for the i-th observation
    if (ai_ind.size() < 2) {
      sw[i] = 0.0;
    } else {
      double sum_a = 0.0;
      for (int k : ai_ind) {
        sum_a += distmat(i, k);
      }
      a[i] = sum_a / (ai_ind.size() - 1);
      
      // Find the nearest cluster
      std::vector<int> other_clusters;
      for (int j = 0; j < n; ++j) {
        if (ocp[j] != cluster && std::find(other_clusters.begin(), other_clusters.end(), ocp[j]) == other_clusters.end()) {
          other_clusters.push_back(ocp[j]);
        }
      }
      
      double min_b = Rcpp::NumericVector::get_na();
      for (int nearest : other_clusters) {
        std::vector<int> nearest_ind;
        for (int j = 0; j < n; ++j) {
          if (ocp[j] == nearest) {
            nearest_ind.push_back(j);
          }
        }
        
        double dist_to_nearest = 0.0;
        for (int k : nearest_ind) {
          dist_to_nearest += distmat(i, k);
        }
        dist_to_nearest /= nearest_ind.size();
        
        if (Rcpp::NumericVector::is_na(min_b) || dist_to_nearest < min_b) {
          min_b = dist_to_nearest;
        }
      }
      b[i] = min_b;
      
      // Compute silhouette width for point i
      sw[i] = (b[i] - a[i]) / std::max(a[i], b[i]);
    }
  }
  
  return Rcpp::mean(sw);
}



// Helper function to create contingency table manually
// [[Rcpp::export]]
arma::mat create_contingency_table(
    Rcpp::NumericVector rcjk,
    Rcpp::NumericVector ocjk
) {
  std::map<std::pair<int, int>, int> counts;
  int n = rcjk.size();
  
  for (int i = 0; i < n; ++i) {
    std::pair<int, int> key = std::make_pair(rcjk[i], ocjk[i]);
    counts[key]++;
  }
  
  // Create a matrix to hold the result
  arma::mat result(counts.size(), 3);
  int idx = 0;
  for (auto const& entry : counts) {
    result(idx, 0) = entry.first.first;
    result(idx, 1) = entry.first.second;
    result(idx, 2) = entry.second;
    idx++;
  }
  
  return result;
}


// [[Rcpp::export]]
double DunnIndexCpp(Rcpp::NumericVector ocp, arma::mat distmat) {
  // Access the base R environment and the 'which' function
  Environment base = Environment::base_env();
  Function which = base["which"];
  
  Rcpp::NumericVector ucp = unique(ocp); // unique cluster identifiers
  int m = ucp.size(); // number of unique clusters
  if (m == 1) return 0; // if only one cluster, return 0
  
  arma::vec num_dist; // vector to hold minimum distances between clusters
  
  // Loop to calculate minimum distance between different clusters
  for (int i = 0; i < m - 1; i++) {
    for (int j = i + 1; j < m; j++) {
      // Use R's which function to find indices of points belonging to clusters i and j
      Rcpp::IntegerVector sc_index = as<Rcpp::IntegerVector>(which(ocp == ucp[i]));
      Rcpp::IntegerVector tc_index = as<Rcpp::IntegerVector>(which(ocp == ucp[j]));
      
      // Convert to zero-based indexing for use with Armadillo (subtract 1 from R's 1-based indexing)
      uvec sc_arma_index = as<uvec>(sc_index) - 1;
      uvec tc_arma_index = as<uvec>(tc_index) - 1;
      
      // Get the minimum distance between points in cluster i and cluster j
      double min_dist = distmat(sc_arma_index, tc_arma_index).min();
      num_dist.insert_rows(num_dist.n_rows, 1); // expand the vector
      num_dist(num_dist.n_rows - 1) = min_dist;
    }
  }
  
  double num = num_dist.min(); // Get the minimum inter-cluster distance
  
  arma::vec denom_dist; // vector to hold maximum distances within clusters
  
  // Loop to calculate maximum distance within each cluster
  for (int i = 0; i < m; i++) {
    // Use R's which function to find indices of points in the cluster
    Rcpp::IntegerVector rc_index = as<Rcpp::IntegerVector>(which(ocp == ucp[i]));
    
    // Convert to zero-based indexing
    uvec rc_arma_index = as<uvec>(rc_index) - 1;
    
    // Get the maximum distance within the cluster
    double max_dist = distmat(rc_arma_index, rc_arma_index).max();
    denom_dist.insert_rows(denom_dist.n_rows, 1); // expand the vector
    denom_dist(denom_dist.n_rows - 1) = max_dist;
  }
  
  double denom = denom_dist.max(); // Get the maximum intra-cluster distance
  
  // Calculate Dunn Index
  double score = num / denom;
  return score;
}




// [[Rcpp::export]]
double BioSIndexCpp(const arma::vec& oc, const arma::vec& bd, const Rcpp::List& rc) {
  int l = oc.n_elem; // number of cells
  double res = 0.0; // overall result
  
  for(int j = 0; j < l; ++j) { // remove cell j
    arma::vec rcj = as<arma::vec>(rc[j]); // get removed cluster groups
    arma::vec ocj = oc;
    ocj.shed_row(j); // get original cluster groups without j
    
    arma::vec bdj = bd;
    bdj.shed_row(j); // bio data when j removed
    arma::vec ubj = unique(bdj); // unique groups when j removed
    int Nj = ubj.n_elem; // number of unique groups when j removed
    
    double jsum = 0.0; // internal sum;
    for(int i = 0; i < Nj; ++i) {
      arma::uvec bdi_index = find(bdj == ubj(i)); // get all cells in i'th unique group
      int Ni = bdi_index.n_elem; // number of cells in i'th group
      
      double isum = 0.0; // internal sum2;
      for(int x = 0; x < Ni; ++x) {
        for(int y = 0; y < Ni; ++y) {
          if(x != y) { // loop all cells in i'th group except self
            arma::uvec ocj_x = find(ocj == ocj(bdi_index(x)));
            arma::uvec rcj_y = find(rcj == rcj(bdi_index(y)));
            arma::uvec inter = intersect(ocj_x, rcj_y);
            isum += static_cast<double>(inter.n_elem) / ocj_x.n_elem;
          }
        }
      }
      
      if(Ni > 1) { // Avoid division by zero
        jsum += isum / (Ni * (Ni - 1));
      }
    }
    res += jsum / Nj;
  }
  return res / l;
}


