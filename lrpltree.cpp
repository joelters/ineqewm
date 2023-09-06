#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector lrpltree(NumericMatrix scores11, NumericMatrix scores10, NumericMatrix scores01, NumericMatrix scores00,
                       NumericMatrix targetX, int depth = 2, NumericMatrix s = 1, NumericVector ns = 1) {

  // 0) treat everyone vs treat no one
  int n = targetX.nrow();
  NumericVector rl0(n, 0.0);
  NumericVector rl1(n, 1.0);
  NumericVector WW0(2, 0.0);
  WW0[0] = -9999;
  WW0[1] = 0;
  NumericVector WW1(4, 1.0);
  WW1[2] = -9999;
  WW1[3] = 12;
  NumericVector WW2(8, 1.0);
  WW2(6) = -9999;
  NumericVector tnodemax(n,0.0);

  int jt = scores11(_,1).size();
  NumericVector rl0_11_1(jt);
  NumericVector rl0_11_2(jt);
  NumericVector rl0_10_1(jt);
  NumericVector rl0_10_2(jt);
  NumericVector rl0_01_1(jt);
  NumericVector rl0_01_2(jt);
  NumericVector rl0_00_1(jt);
  NumericVector rl0_00_2(jt);

  NumericVector rl1_11_1(jt);
  NumericVector rl1_11_2(jt);
  NumericVector rl1_10_1(jt);
  NumericVector rl1_10_2(jt);
  NumericVector rl1_01_1(jt);
  NumericVector rl1_01_2(jt);
  NumericVector rl1_00_1(jt);
  NumericVector rl1_00_2(jt);

  NumericVector rl12_11_1(jt);
  NumericVector rl12_11_2(jt);
  NumericVector rl12_10_1(jt);
  NumericVector rl12_10_2(jt);
  NumericVector rl12_01_1(jt);
  NumericVector rl12_01_2(jt);
  NumericVector rl12_00_1(jt);
  NumericVector rl12_00_2(jt);

  NumericVector rl21_11_1(jt);
  NumericVector rl21_11_2(jt);
  NumericVector rl21_10_1(jt);
  NumericVector rl21_10_2(jt);
  NumericVector rl21_01_1(jt);
  NumericVector rl21_01_2(jt);
  NumericVector rl21_00_1(jt);
  NumericVector rl21_00_2(jt);

  NumericVector rls_11_1(jt);
  NumericVector rls_11_2(jt);
  NumericVector rls_10_1(jt);
  NumericVector rls_10_2(jt);
  NumericVector rls_01_1(jt);
  NumericVector rls_01_2(jt);
  NumericVector rls_00_1(jt);
  NumericVector rls_00_2(jt);

  for (int i = 0; i < jt; i++) {
    rl0_11_1[i] = rl0[scores11(_,1)[i] - 1];
    rl0_11_2[i] = rl0[scores11(_,2)[i] - 1];
    rl0_10_1[i] = rl0[scores10(_,1)[i] - 1];
    rl0_10_2[i] = rl0[scores10(_,2)[i] - 1];
    rl0_01_1[i] = rl0[scores01(_,1)[i] - 1];
    rl0_01_2[i] = rl0[scores01(_,2)[i] - 1];
    rl0_00_1[i] = rl0[scores00(_,1)[i] - 1];
    rl0_00_2[i] = rl0[scores00(_,2)[i] - 1];

    rl1_11_1[i] = rl1[scores11(_,1)[i] - 1];
    rl1_11_2[i] = rl1[scores11(_,2)[i] - 1];
    rl1_10_1[i] = rl1[scores10(_,1)[i] - 1];
    rl1_10_2[i] = rl1[scores10(_,2)[i] - 1];
    rl1_01_1[i] = rl1[scores01(_,1)[i] - 1];
    rl1_01_2[i] = rl1[scores01(_,2)[i] - 1];
    rl1_00_1[i] = rl1[scores00(_,1)[i] - 1];
    rl1_00_2[i] = rl1[scores00(_,2)[i] - 1];
  }

  double W0 = (2 / (n * (n - 1))) * sum(scores11(_, 0) * rl0_11_1 * rl0_11_2 +
               scores10(_, 0) * rl0_10_1 * (1 - rl0_10_2) +
               scores01(_, 0) * (1 - rl0_01_1) * rl0_01_2 +
               scores00(_, 0) * (1 - rl0_00_1) * (1 - rl0_00_2));
  NumericVector a0(2,0.0);
  a0[0] = W0;
  a0[1] = 0;



  double W1 = (2 / (n * (n - 1))) * sum(scores11(_, 0) * rl1_11_1 * rl1_11_2 +
               scores10(_, 0) * rl1_10_1 * (1 - rl1_10_2) +
               scores01(_, 0) * (1 - rl1_01_1) * rl1_01_2 +
               scores00(_, 0) * (1 - rl1_00_1) * (1 - rl1_00_2));
  NumericVector a1(2,0.0);
  a1[0] = W1;
  a1[1] = 1;
  NumericMatrix WW0aux(3,2);
  WW0aux(0,_) = a0;
  WW0aux(1,_) = a1;
  WW0aux(2,_) = WW0;

  NumericVector WW0s(3,0.0);
  WW0s[0] = W0;
  WW0s[1] = W1;
  WW0s[2] = WW0[0];

  int index = 0;
  double maxVal = WW0s[0];
  for (int i = 1; i < WW0s.size(); i++) {
    if (WW0s[i] > maxVal) {
      maxVal = WW0s[i];
      index = i;
    }
  }

  WW0 = WW0aux(index,_);

    if (depth == 0) {
      return WW0;
    }

    // 1) Do single split (rl12 means treat node1 and not node 2, rl21 treat node 2 and not node 1)
    for (int ii = 0; ii < targetX.ncol(); ii++) {
      for (int jj = 0; jj < ns[ii]; jj++) {
        int nn = targetX(_, ii).size();
        NumericVector rl12(nn);
        for (int i = 0; i < nn; i++) {
          if (targetX(_, ii)[i] <= s(jj,ii)) {
            rl12[i] = 1;
          } else {
            rl12[i] = 0;
          }
        }

        for (int i = 0; i < jt; i++) {
          rl12_11_1[i] = rl12[scores11(_,1)[i] - 1];
          rl12_11_2[i] = rl12[scores11(_,2)[i] - 1];
          rl12_10_1[i] = rl12[scores10(_,1)[i] - 1];
          rl12_10_2[i] = rl12[scores10(_,2)[i] - 1];
          rl12_01_1[i] = rl12[scores01(_,1)[i] - 1];
          rl12_01_2[i] = rl12[scores01(_,2)[i] - 1];
          rl12_00_1[i] = rl12[scores00(_,1)[i] - 1];
          rl12_00_2[i] = rl12[scores00(_,2)[i] - 1];
        }

        double W12 = (2 / (n * (n - 1))) * sum(scores11(_, 0) * rl12_11_1 * rl12_11_2 +
                      scores10(_, 0) * rl12_10_1 * (1 - rl12_10_2) +
                      scores01(_, 0) * (1 - rl12_01_1) * rl12_01_2 +
                      scores00(_, 0) * (1 - rl12_00_1) * (1 - rl12_00_2));

        int nn1 = targetX(_, ii).size();
        NumericVector rl21(nn1);
        for (int i = 0; i < nn1; i++) {
          if (targetX(_, ii)[i] > s(jj, ii)) {
            rl21[i] = 1;
          } else {
            rl21[i] = 0;
          }
        }


        for (int i = 0; i < jt; i++) {
          rl21_11_1[i] = rl21[scores11(_,1)[i] - 1];
          rl21_11_2[i] = rl21[scores11(_,2)[i] - 1];
          rl21_10_1[i] = rl21[scores10(_,1)[i] - 1];
          rl21_10_2[i] = rl21[scores10(_,2)[i] - 1];
          rl21_01_1[i] = rl21[scores01(_,1)[i] - 1];
          rl21_01_2[i] = rl21[scores01(_,2)[i] - 1];
          rl21_00_1[i] = rl21[scores00(_,1)[i] - 1];
          rl21_00_2[i] = rl21[scores00(_,2)[i] - 1];
        }

        double W21 = (2 / (n * (n - 1))) * sum(scores11(_, 0) * rl21_11_1 * rl21_11_2 +
                      scores10(_, 0) * rl21_10_1 * (1 - rl21_10_2) +
                      scores01(_, 0) * (1 - rl21_01_1) * rl21_01_2 +
                      scores00(_, 0) * (1 - rl21_00_1) * (1 - rl21_00_2));

        NumericVector a12(4,0.0);
        a12[0] = ii + 1;
        a12[1] = jj + 1;
        a12[2] = W12;
        a12[3] = 12;

        NumericVector a21(4,0.0);
        a21[0] = ii + 1;
        a21[1] = jj + 1;
        a21[2] = W21;
        a21[3] = 21;

        NumericMatrix WW1221aux(3,4);
        WW1221aux(0,_) = a12;
        WW1221aux(1,_) = a21;
        WW1221aux(2,_) = WW1;

        NumericVector WW1221s(3,0.0);
          WW1221s[0] = W12;
        WW1221s[1] = W21;
        WW1221s[2] = WW1[2];

        int index = 0;
        double maxVal = WW1221s[0];
        for (int i = 1; i < WW1221s.size(); i++) {
          if (WW1221s[i] > maxVal) {
            maxVal = WW1221s[i];
            index = i;
          }
        }

        WW1 = WW1221aux(index,_);

        // 2 Depth 2 tree
        if (depth == 2) {
          //sub-trees
          for (int kk = 0; kk < targetX.ncol(); kk++) {
            for (int ll = 0; ll < ns[kk]; ll++) {
              for (int pp = 0; pp < targetX.ncol(); pp++) {
                for (int tt = 0; tt < ns[pp]; tt++) {

                  int nn2 = targetX(_, kk).size();
                  NumericVector rl12_1(nn2);
                  NumericVector rl12_2(nn2);
                  NumericVector rl12_3(nn2);
                  NumericVector rl12_4(nn2);
                  for (int i = 0; i < nn2; i++) {
                    if (targetX(_, kk)[i] <= s(ll, kk)) {
                      rl12_1[i] = 1;
                    } else {
                      rl12_1[i] = 0;
                    }

                    if (targetX(_, kk)[i] > s(ll, kk)) {
                      rl12_2[i] = 1;
                    } else {
                      rl12_2[i] = 0;
                    }

                    if (targetX(_, pp)[i] <= s(tt, pp)) {
                      rl12_3[i] = 1;
                    } else {
                      rl12_3[i] = 0;
                    }

                    if (targetX(_, pp)[i] > s(tt, pp)) {
                      rl12_4[i] = 1;
                    } else {
                      rl12_4[i] = 0;
                    }
                  }
                  NumericVector tnode = rl12 * rl12_1 +
                    2 * rl12 * rl12_2 +
                    3 * rl21 * rl12_3 +
                    4 * rl21 * rl12_4;

                  NumericVector rl_ttnt(n,0.0);
                  NumericVector rl_tttn(n,0.0);
                  NumericVector rl_tntn(n,0.0);
                  NumericVector rl_tnnt(n,0.0);
                  NumericVector rl_tntt(n,0.0);
                  NumericVector rl_tnnn(n,0.0);
                  NumericVector rl_ntnt(n,0.0);
                  NumericVector rl_nttn(n,0.0);
                  NumericVector rl_nttt(n,0.0);
                  NumericVector rl_ntnn(n,0.0);
                  NumericVector rl_nnnt(n,0.0);
                  NumericVector rl_nntn(n,0.0);

                  rl_ttnt = (tnode == 1) | (tnode == 2) | (tnode == 4);
                  rl_tttn = (tnode == 1) | (tnode == 2) | (tnode == 3);

                  rl_tntn = (tnode == 1) | (tnode == 3);
                  rl_tnnt = (tnode == 1) | (tnode == 4);
                  rl_tntt = (tnode == 1) | (tnode == 3) | (tnode == 4);
                  rl_tnnn = (tnode == 1);

                  rl_ntnt = (tnode == 2) | (tnode == 4);
                  rl_nttn = (tnode == 2) | (tnode == 3);
                  rl_nttt = (tnode == 2) | (tnode == 3) | (tnode == 4);
                  rl_ntnn = (tnode == 2);

                  rl_nnnt = (tnode == 4);
                  rl_nntn = (tnode == 3);

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

                    for (int i = 0; i < jt; i++) {
                      rls_11_1[i] = rls(_,rr)[scores11(_,1)[i] - 1];
                      rls_11_2[i] = rls(_,rr)[scores11(_,2)[i] - 1];
                      rls_10_1[i] = rls(_,rr)[scores10(_,1)[i] - 1];
                      rls_10_2[i] = rls(_,rr)[scores10(_,2)[i] - 1];
                      rls_01_1[i] = rls(_,rr)[scores01(_,1)[i] - 1];
                      rls_01_2[i] = rls(_,rr)[scores01(_,2)[i] - 1];
                      rls_00_1[i] = rls(_,rr)[scores00(_,1)[i] - 1];
                      rls_00_2[i] = rls(_,rr)[scores00(_,2)[i] - 1];
                    }

                    double W1234 = (2 / (n * (n - 1))) * sum(scores11(_, 0) * rls_11_1 * rls_11_2 +
                                    scores10(_, 0) * rls_10_1 * (1 - rls_10_2) +
                                    scores01(_, 0) * (1 - rls_01_1) * rls_01_2 +
                                    scores00(_, 0) * (1 - rls_00_1) * (1 - rls_00_2));
                    NumericVector aa2(8,0.0);
                    aa2[0] = ii + 1;
                    aa2[1] = jj + 1;
                    aa2[2] = kk + 1;
                    aa2[3] = ll + 1;
                    aa2[4] = pp + 1;
                    aa2[5] = tt + 1;
                    aa2[6] = W1234;
                    aa2[7] = rr + 1;

                    if (W1234 >= WW2[6]){
                      WW2 = aa2;
                      tnodemax = tnode;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    if (depth == 2){
      NumericMatrix aux2(n,10);
      aux2.fill(0);
      WW2.push_back(W0);
      aux2(0,_) = WW2;
      aux2(_,9) = tnodemax;
      return aux2;
    }
    else{
      WW1.push_back(W0);
      return WW1;
    }
}
