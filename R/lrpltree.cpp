#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector lrpltree(NumericMatrix scores11, NumericMatrix scores10, NumericMatrix scores01, NumericMatrix scores00,
                       NumericMatrix targetX, int depth = 2, NumericMatrix s, NumericVector ns) {

  // 0) treat everyone vs treat no one
int n = targetX.nrow();
NumericVector rl0(n, 0.0);
NumericMatrix WW0(2, 2);
WW0(0, 0) = -9999;
WW0(0, 1) = 0;
double W0 = (2 / (n * (n - 1))) * sum(scores11(_, 0) * rl0[scores11(_, 1)] * rl0[scores11(_, 2)] +
                                        scores10(_, 0) * rl0[scores10(_, 1)] * (1 - rl0[scores10(_, 2)]) +
                                        scores01(_, 0) * (1 - rl0[scores01(_, 1)]) * rl0[scores01(_, 2)] +
                                        scores00(_, 0) * (1 - rl0[scores00(_, 1)]) * (1 - rl0[scores00(_, 2)]));

NumericVector rl1(n, 1.0);
double W1 = (2 / (n * (n - 1))) * sum(scores11(_, 0) * rl1[scores11(_, 1)] * rl1[scores11(_, 2)] +
                                        scores10(_, 0) * rl1[scores10(_, 1)] * (1 - rl1[scores10(_, 2)]) +
                                        scores01(_, 0) * (1 - rl1[scores01(_, 1)]) * rl1[scores01(_, 2)] +
                                        scores00(_, 0) * (1 - rl1[scores00(_, 1)]) * (1 - rl1[scores00(_, 2)]));
WW0(0, 0) = W0;
WW0(1, 0) = W1;
WW0(1, 1) = 1;
WW0 = WW0(Rcpp::which_max(Rcpp::c(W0, W1, WW0(0, 0))), _);

if (depth == 0) {
  return WW0;
}

// 1) Do single split (rl12 means treat node1 and not node 2, rl21 treat node 2 and not node 1)
NumericVector WW1(4, 0.0);
WW1(2) = -9999;
WW1(3) = 12;
for (int ii = 0; ii < targetX.ncol(); ii++) {
  for (int jj = 0; jj < ns[ii]; jj++) {
    NumericVector rl12 = (targetX(_, ii) <= s(jj, ii));
    double W12 = (2 / (n * (n - 1))) * sum(scores11(_, 0) * rl12[scores11(_, 1)] * rl12[scores11(_, 2)] +
                                             scores10(_, 0) * rl12[scores10(_, 1)] * (1 - rl12[scores10(_, 2)]) +
                                             scores01(_, 0) * (1 - rl12[scores01(_, 1)]) * rl12[scores01(_, 2)] +
                                             scores00(_, 0) * (1 - rl12[scores00(_, 1)]) * (1 - rl12[scores00(_, 2)]));

    NumericVector rl21 = (targetX(_, ii) > s(jj, ii));
    double W21 = (2 / (n * (n - 1))) * sum(scores11(_, 0) * rl21[scores11(_, 1)] * rl21[scores11(_, 2)] +
                                             scores10(_, 0) * rl21[scores10(_, 1)] * (1 - rl21[scores10(_, 2)]) +
                                             scores01(_, 0) * (1 - rl21[scores01(_, 1)]) * rl21[scores01(_, 2)] +
                                             scores00(_, 0) * (1 - rl21[scores00(_, 1)]) * (1 - rl21[scores00(_, 2)]));
    WW1(0) = ii + 1;
    WW1(1) = jj + 1;
    WW1(2) = W12;
    WW1(3) = 12;
    WW1 = WW1(Rcpp::which_max(Rcpp::c(W12, W21, WW1(2))), _);

    // 2 Depth 2 tree
    if (depth == 2) {
      //sub-trees
      NumericVector WW2(8, 0.0);
      WW2(6) = -9999;
      WW2(7) = 1;
      for (int kk = 0; kk < targetX.ncol(); kk++) {
        for (int ll = 0; ll < ns[kk]; ll++) {
          for (int pp = 0; pp < targetX.ncol(); pp++) {
            for (int tt = 0; tt < ns[pp]; tt++) {
              NumericVector tnode = rl12 * (targetX(_, kk) <= s(ll, kk)) +
                2 * rl12 * (targetX(_, kk) > s(ll, kk)) +
                3 * rl21 * (targetX(_, pp) <= s(tt, pp)) +
                4 * rl21 * (targetX(_, pp) > s(tt, pp));

              NumericVector rl_ttnt = (tnode == 1 | tnode == 2 | tnode == 4);
              NumericVector rl_tttn = (tnode == 1 | tnode == 2 | tnode == 3);

              NumericVector rl_tntn = (tnode == 1 | tnode == 3);
              NumericVector rl_tnnt = (tnode == 1 | tnode == 4);
              NumericVector rl_tntt = (tnode == 1 | tnode == 3 | tnode == 4);
              NumericVector rl_tnnn = (tnode == 1);

              NumericVector rl_ntnt = (tnode == 2 | tnode == 4);
              NumericVector rl_nttn = (tnode == 2 | tnode == 3);
              NumericVector rl_nttt = (tnode == 2 | tnode == 3 | tnode == 4);
              NumericVector rl_ntnn = (tnode == 2);

              NumericVector rl_nnnt = (tnode == 4);
              NumericVector rl_nntn = (tnode == 3);

              NumericMatrix rls(n, 12);
              rls(_, 0) = rl_ttnt;
              rls(_, 1) = rl_tttn;
              rls(_, 2) = rl_tntn;
              rls(_, 3) = rl_tnnt;
              rls(_, 4) = rl_tntt;
              rls(_, 5) = rl_tnnn;
              rls(_, 6) = rl_ntnt;
              rls(_, 7) = rl_nttn;
              rls(_, 8) = rl_nttt;
              rls(_, 9) = rl_ntnn;
              rls(_, 10) = rl_nnnt;
              rls(_, 11) = rl_nntn;

              for (int rr = 0; rr < rls.ncol(); rr++) {
                double W1234 = (2 / (n * (n - 1))) * sum(scores11(_, 0) * rls(_, rr)[scores11(_, 1)] * rls(_, rr)[scores11(_, 2)] +
                                                           scores10(_, 0) * rls(_, rr)[scores10(_, 1)] * (1 - rls(_, rr)[scores10(_, 2)]) +
                                                           scores01(_, 0) * (1 - rls(_, rr)[scores01(_, 1)]) * rls(_, rr)[scores01(_, 2)] +
                                                           scores00(_, 0) * (1 - rls(_, rr)[scores00(_, 1)]) * (1 - rls(_, rr)[scores00(_, 2)]));
                WW2(0) = ii + 1;
                WW2(1) = jj + 1;
                WW2(2) = kk + 1;
                WW2(3) = ll + 1;
                WW2(4) = pp + 1;
                WW2(5) = tt + 1;
                WW2(6) = W1234;
                WW2(7) = rr + 1;
                WW2 = WW2(Rcpp::which_max(Rcpp::c(W1234, WW2(6))), _);
              }
            }
          }
        }
      }
      return Rcpp::cbind(Rcpp::rbind(WW2, Rcpp::NumericMatrix(n - 1, 9)), tnodemax);
    }
  }
}
return Rcpp::cbind(WW1, W0);
}
