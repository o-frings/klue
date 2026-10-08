// RcppArmadillo kernel for the LCMNL hot path.
//
// This is a line-for-line port of the pure-R reference kernel in
// R/03_estimate.R (.mnl_eval, .mnl_score, .panel_sum, the C>=2 neg_ll/grad_ll
// and the EM E-step). The R version stays the reference implementation; a
// testthat case asserts the two agree to 1e-8 on fixed inputs. The design
// arrays are copied once into an XPtr<LcmnlData> per estimation call, so the
// optimiser's ~100 nll/grad evaluations touch no R<->C++ data marshalling.
//
// Parameter layout per class: [beta_1..beta_{n_beta}, asc_1..asc_{n_asc}],
// the reference-alternative ASC fixed at 0. For C>=2 the (C-1) free class
// deltas are appended (last class delta = 0).

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

// Design data, built once per estimation and reused across optim iterations.
struct LcmnlData {
  arma::cube X;       // (T_total, n_beta, J): slice j = X[[j]]
  arma::mat  ch_ind;  // (T_total, J): 1 where alt j chosen
  arma::mat  Xc;      // (T_total, n_beta): chosen-alternative attributes
  int T_per_n, N, n_beta, n_asc, J;
};

// Per-task log-likelihood (and optionally choice probabilities) for one class.
// betas: length n_beta; asc_full: length J with the reference entry 0.
// Mirrors .mnl_eval: V[,j] = X[[j]] %*% betas + asc[j]; row-max Vm;
// tll = rowSums(ch_ind * V) - Vm - log(rowSums(exp(V - Vm))).
static inline void mnl_core(const LcmnlData& d,
                            const arma::vec& betas, const arma::vec& asc_full,
                            arma::vec& tll, arma::mat& probs, bool want_probs) {
  const uword T = d.X.n_rows;
  const uword J = d.X.n_slices;
  arma::mat V(T, J);
  for (uword j = 0; j < J; ++j) V.col(j) = d.X.slice(j) * betas + asc_full(j);
  arma::vec Vm = max(V, 1);                 // rowwise max -> length T
  arma::mat eV = exp(V.each_col() - Vm);    // subtract Vm from each column
  arma::vec rse = sum(eV, 1);               // rowsum -> length T
  arma::vec chosenV = sum(d.ch_ind % V, 1); // rowSums(ch_ind * V)
  tll = chosenV - Vm - log(rse);
  if (want_probs) { probs = eV; probs.each_col() /= rse; }
}

// colSums(matrix(x, T_per_n, N)): sum each respondent's T_per_n-row block.
// No-copy view over x's memory (column-major, matches R matrix()).
static inline arma::vec panel_sum(const arma::vec& x, int T_per_n, int N) {
  arma::mat M(const_cast<double*>(x.memptr()), T_per_n, N, false, true);
  return vectorise(sum(M, 0));              // colSums -> length N
}

// Split a class parameter block into betas and the full-length ASC vector.
static inline void split_class_par(const LcmnlData& d, const arma::vec& par,
                                    int off, arma::vec& betas, arma::vec& asc_full) {
  betas = par.subvec(off, off + d.n_beta - 1);
  asc_full.zeros(d.J);
  for (int a = 0; a < d.n_asc; ++a) asc_full(a) = par(off + d.n_beta + a);
}

// Panel log-likelihood (length N) for one class given its npc-vector.
static inline arma::vec class_panel_loglik(const LcmnlData& d, const arma::vec& par,
                                           int off) {
  arma::vec betas, asc_full, tll; arma::mat probs;
  split_class_par(d, par, off, betas, asc_full);
  mnl_core(d, betas, asc_full, tll, probs, false);
  return panel_sum(tll, d.T_per_n, d.N);
}

// Row-wise softmax of (log_panel + log_pi): returns weights and sample LL.
static inline void posterior_weights(const arma::mat& log_panel,
                                     const arma::vec& log_pi,
                                     arma::mat& w, double& LL) {
  arma::mat log_joint = log_panel.each_row() + log_pi.t();
  arma::vec lm = max(log_joint, 1);
  arma::vec denom = lm + log(sum(exp(log_joint.each_col() - lm), 1));
  w = exp(log_joint.each_col() - denom);
  LL = accu(denom);
}

// ---------------------------------------------------------------------------
// Exported entry points
// ---------------------------------------------------------------------------

//' @keywords internal
// [[Rcpp::export]]
SEXP cpp_lcmnl_data(const arma::cube& X, const arma::mat& ch_ind,
                    const arma::mat& Xc, int T_per_n, int N,
                    int n_beta, int n_asc, int J) {
  LcmnlData* d = new LcmnlData{X, ch_ind, Xc, T_per_n, N, n_beta, n_asc, J};
  Rcpp::XPtr<LcmnlData> ptr(d, true);
  return ptr;
}

// Single-class (optionally row-weighted) MNL negative log-likelihood.
// [[Rcpp::export]]
double cpp_mnl_nll(SEXP data, const arma::vec& par,
                   Rcpp::Nullable<Rcpp::NumericVector> rw_ = R_NilValue) {
  Rcpp::XPtr<LcmnlData> d(data);
  arma::vec betas, asc_full, tll; arma::mat probs;
  split_class_par(*d, par, 0, betas, asc_full);
  mnl_core(*d, betas, asc_full, tll, probs, false);
  if (rw_.isNull()) return -accu(tll);
  arma::vec rw = Rcpp::as<arma::vec>(rw_.get());
  return -accu(rw % tll);
}

// Gradient of the single-class (optionally weighted) MNL: length npc.
// [[Rcpp::export]]
arma::vec cpp_mnl_grad(SEXP data, const arma::vec& par,
                       Rcpp::Nullable<Rcpp::NumericVector> rw_ = R_NilValue) {
  Rcpp::XPtr<LcmnlData> d(data);
  arma::vec betas, asc_full, tll; arma::mat probs;
  split_class_par(*d, par, 0, betas, asc_full);
  mnl_core(*d, betas, asc_full, tll, probs, true);

  arma::mat EX(d->X.n_rows, d->n_beta, fill::zeros);
  for (int j = 0; j < d->J; ++j) EX += d->X.slice(j).each_col() % probs.col(j);
  arma::mat sbeta = d->Xc - EX;             // T x n_beta

  arma::vec g(d->n_beta + d->n_asc);
  bool has_w = !rw_.isNull();
  arma::vec rw;
  if (has_w) rw = Rcpp::as<arma::vec>(rw_.get());
  for (int k = 0; k < d->n_beta; ++k)
    g(k) = has_w ? -accu(rw % sbeta.col(k)) : -accu(sbeta.col(k));
  for (int k = 0; k < d->n_asc; ++k) {
    arma::vec sasc = d->ch_ind.col(k) - probs.col(k);
    g(d->n_beta + k) = has_w ? -accu(rw % sasc) : -accu(sasc);
  }
  return g;
}

// Panel log-likelihood (length N) for one class npc-vector (EM E-step).
// [[Rcpp::export]]
arma::vec cpp_mnl_panel_loglik(SEXP data, const arma::vec& par) {
  Rcpp::XPtr<LcmnlData> d(data);
  return class_panel_loglik(*d, par, 0);
}

// Negative LCMNL log-likelihood (C >= 2).
// [[Rcpp::export]]
double cpp_lcmnl_nll(SEXP data, const arma::vec& par, int C) {
  Rcpp::XPtr<LcmnlData> d(data);
  int npc = d->n_beta + d->n_asc;
  arma::vec deltas(C, fill::zeros);
  for (int c = 0; c < C - 1; ++c) deltas(c) = par(C * npc + c);
  double dm = deltas.max();
  arma::vec log_pi = deltas - dm - std::log(accu(exp(deltas - dm)));

  arma::mat log_panel(d->N, C);
  for (int c = 0; c < C; ++c)
    log_panel.col(c) = class_panel_loglik(*d, par, c * npc);

  arma::mat log_joint = log_panel.each_row() + log_pi.t();
  arma::vec lm = max(log_joint, 1);
  arma::vec denom = lm + log(sum(exp(log_joint.each_col() - lm), 1));
  return -accu(denom);
}

// Gradient of the LCMNL log-likelihood (C >= 2): length C*npc + (C-1).
// [[Rcpp::export]]
arma::vec cpp_lcmnl_grad(SEXP data, const arma::vec& par, int C) {
  Rcpp::XPtr<LcmnlData> d(data);
  int npc = d->n_beta + d->n_asc;
  int n_free = C * npc + (C - 1);

  arma::vec deltas(C, fill::zeros);
  for (int c = 0; c < C - 1; ++c) deltas(c) = par(C * npc + c);
  double dm = deltas.max();
  arma::vec ed = exp(deltas - dm);
  arma::vec pi_c = ed / accu(ed);

  arma::mat log_panel(d->N, C);
  std::vector<arma::mat> pscores(C);        // each N x npc
  for (int c = 0; c < C; ++c) {
    arma::vec betas, asc_full, tll; arma::mat probs;
    split_class_par(*d, par, c * npc, betas, asc_full);
    mnl_core(*d, betas, asc_full, tll, probs, true);
    log_panel.col(c) = panel_sum(tll, d->T_per_n, d->N);

    arma::mat EX(d->X.n_rows, d->n_beta, fill::zeros);
    for (int j = 0; j < d->J; ++j) EX += d->X.slice(j).each_col() % probs.col(j);
    arma::mat sbeta = d->Xc - EX;           // T x n_beta

    arma::mat ps(d->N, npc);
    for (int k = 0; k < d->n_beta; ++k)
      ps.col(k) = panel_sum(sbeta.col(k), d->T_per_n, d->N);
    for (int k = 0; k < d->n_asc; ++k) {
      arma::vec sasc = d->ch_ind.col(k) - probs.col(k);
      ps.col(d->n_beta + k) = panel_sum(sasc, d->T_per_n, d->N);
    }
    pscores[c] = ps;
  }

  arma::mat w; double LL;
  posterior_weights(log_panel, log(pi_c), w, LL);

  arma::vec g(n_free, fill::zeros);
  for (int c = 0; c < C; ++c)               // -colSums(w[,c] * pscores[[c]])
    g.subvec(c * npc, c * npc + npc - 1) = -(pscores[c].t() * w.col(c));
  for (int c = 0; c < C - 1; ++c)
    g(C * npc + c) = -(accu(w.col(c)) - (double) d->N * pi_c(c));
  return g;
}
