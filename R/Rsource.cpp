// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <RcppParallel.h>
#include <unsupported/Eigen/SpecialFunctions>
#include <unsupported/Eigen/KroneckerProduct>
#include <iostream>
#include <vector>
using namespace Rcpp;
using namespace Eigen;
using namespace RcppParallel;
using namespace std;

// [[Rcpp::export]]
MatrixXd test(MatrixXd mat1, MatrixXd mat2){
  VectorXd newCol(2);
  newCol[0] = 4;
  newCol[1] = 6;
  MatrixXd out = MatrixXd::Zero(2, 2);
  out.col(0) = newCol;
  return out;
}

VectorXd ifElse(const VectorXd &a, const double b, const VectorXd &return1, const double return2){
  int len = a.size();

  VectorXd res(len);
  for(int i = 0 ; i < len ; i++){
    if(a[i] > b){
      res[i] = return1[i];
    } else{
      res[i] = return2;
    }
  }
  return res;
}

MatrixXd ifElseMat(const MatrixXd&a, const double b, const MatrixXd&return1, const double return2){
  int nCol = a.cols();
  int nRow = a.rows();

  MatrixXd res(nRow, nCol);
  for(int i = 0 ; i < nRow ; i++){
    for(int j = 0 ; j < nCol ; j++){
      if(a(i, j) > b){
        res(i, j) = return1(i, j);
      } else{
        res(i, j) = return2;
      }
    }
  }
  return res;
}

// [[Rcpp::export]]
VectorXd EMCpp(MatrixXd X, VectorXd y, int maxiter=50, double maxerr=0.01, double lim=0.01){
  int XCol = X.cols();
  VectorXd beta0 = VectorXd::Constant(XCol, y.sum()/XCol);
  VectorXd beta;

  VectorXd cSum = X.colwise().sum();
  for(int i = 0 ; i < maxiter ; i++){
    MatrixXd X1 = (X.transpose().array().colwise() * beta0.array()).transpose();
    VectorXd rSum = X1.rowwise().sum();
    rSum = ifElse(rSum, 0, rSum, 1);
    MatrixXd X2 = X1.array().colwise()/rSum.array();
    beta = ((X2.transpose() * y).array()) / cSum.array();
    VectorXd norm = ifElse(beta0, lim, beta0, lim);
    VectorXd err = ((beta - beta0).array().abs())/norm.array();
    if(err.maxCoeff() < maxerr) break;
    beta0 = beta;
  }
  return beta;
}


// [[Rcpp::export]]
MatrixXd fn1C(MatrixXd X, VectorXd bt, double logNormI, bool modify=true){
  VectorXd gamma0 = bt;
  if(modify){
    VectorXd smallValue = VectorXd::Constant(gamma0.size(), 0.00000001);
    VectorXd logNorm = VectorXd::Constant(gamma0.size(), logNormI);
    gamma0 = Eigen::digamma((bt+smallValue).array());
    gamma0 = gamma0 - logNorm;
    gamma0 = gamma0.array().exp();
  }
  MatrixXd tmp = (X.transpose().array().colwise() * gamma0.array()).transpose();
  VectorXd rSum = tmp.rowwise().sum();
  rSum = ifElse(rSum, 0, rSum, 1);
  MatrixXd res = tmp.array().colwise()/rSum.array();
  return(res);
}

// [[Rcpp::export]]
VectorXd fn2C(MatrixXd Xi, VectorXd Yi, VectorXd cSum){
  MatrixXd bt = Xi.transpose() * Yi;
  bt = bt.array()/cSum.array();
  return(bt);
}

// [[Rcpp::export]]
MatrixXd estimateBETA(const MatrixXd X, const MatrixXd Ymat, VectorXd bet,
                      int maxiter=100, double maxerr=0.01, double lim=0.01, bool modify=true){
  Eigen::Map<VectorXd> beta0(bet.data(), bet.size());
  int q = X.rows();
  int nSample = Ymat.cols();
  VectorXd cSum = X.colwise().sum();
  VectorXd totCount = Ymat.colwise().sum();
  VectorXd smallValue = VectorXd::Constant(totCount.size(), 0.00000001);
  VectorXd logNorm = Eigen::digamma((totCount + smallValue).array());

  Rcpp::List YList(nSample);
  for(int i = 0 ; i < nSample ; i++){
    VectorXd temp = Ymat.col(i);
    YList[i] = temp;
  }

  Rcpp::List BTList(nSample);
  Rcpp::List X2List(nSample);
  int nIsoform = X.cols();
  MatrixXd beta(nIsoform, nSample);

  for(int l = 0 ; l < maxiter ; l++){
    for(int i = 0 ; i < nSample ; i++){
      VectorXd temp1 = beta0.segment(i*nIsoform, nIsoform);
      BTList[i] = temp1;
    }

    for(int i = 0 ; i < nSample; i++){
      MatrixXd temp = fn1C(X, BTList[i], logNorm[i]);
      X2List[i] = temp;
    }

    for(int i = 0 ; i < nSample ; i++){
      beta.col(i) = fn2C(X2List[i], YList[i], cSum);
    }
    Eigen::Map<VectorXd> vecBeta(beta.data(), beta.size());
    VectorXd norm = ifElse(beta0, lim, beta0, lim);
    VectorXd err = ((vecBeta.array() - beta0.array()).abs())/norm.array();

    if(err.maxCoeff() < maxerr) break;
    beta0 = vecBeta;
  }
  return(beta.transpose());
}

// [[Rcpp::export]]
MatrixXd estimateX(MatrixXd xEM, MatrixXd beta, MatrixXd yMat){
  MatrixXd x0 = xEM;
  Eigen::Map<VectorXd> y(yMat.data(), yMat.size());
  int q = yMat.rows();
  int n = yMat.cols();
  int p = beta.cols();

  VectorXd betaColSum = beta.colwise().sum();

  VectorXd BTSum(p*q);
  for(int i = 0 ; i < p ; i++){
    for(int j = 0 ; j < q ; j++){
      BTSum(i*q+j) = betaColSum(i);
    }
  }
  MatrixXd x1(q*beta.rows(), p);
  for(int i = 0 ; i < p ; i++){
    x1.col(i) = Eigen::kroneckerProduct(beta.col(i), x0.col(i));
  }
  VectorXd x1RowSum = x1.rowwise().sum();
  VectorXd x1Norm = ifElse(x1RowSum, 0, x1RowSum, 0.1);
  MatrixXd x2 = x1.array().colwise() / x1Norm.array();
  MatrixXd yCell = x2.array().colwise() * y.array();

  VectorXd tmp(q*p);
  for(int i = 0 ; i < p ; i++){
    Eigen::Map<MatrixXd> temp1(yCell.col(i).data(), q, n);
    MatrixXd temp = temp1.transpose();
    VectorXd tempColSum = temp.colwise().sum();
    tmp.segment(i*q, q) = tempColSum;
  }

  VectorXd tmpNorm = ifElse(BTSum, 0, BTSum, 0.1);
  tmp = tmp.array()/tmpNorm.array();
  Eigen::Map<MatrixXd> res(tmp.data(), q, p);
  VectorXd cSum = res.colwise().sum();
  VectorXd resNorm = ifElse(cSum, 0, cSum, 0.1);
  MatrixXd resT = res.transpose();
  resT = resT.array().colwise()/resNorm.array();
  res = resT.transpose();
  
  for(int i = 0 ; i < p ; i++){
    if(cSum(i) < 0.00001){
      res.col(i) = xEM.col(i);
    }
  }

  return(res);
}

// [[Rcpp::export]]
List AEMCpp(const MatrixXd x0, const MatrixXd yMat, int maxiterX=10, double maxErr=0.01,
                double lim=0.01, int maxiterBT=100, bool modify=true){
  MatrixXd xEST0 = x0;
  VectorXd ySum = yMat.colwise().sum();
  int x0Col = x0.cols();
  ySum = ySum.array()/x0Col;
  int lenYSum = ySum.size();

  VectorXd beta0(x0Col * lenYSum);
  for(int i = 0 ; i < lenYSum ; i++){
    for(int j = 0 ; j < x0Col ; j++){
      beta0(i*x0Col+j) = ySum(i);
    }
  }
  MatrixXd beta;
  MatrixXd xEST;
  for(int i = 0 ; i < maxiterX ; i++){
    beta = estimateBETA(xEST0, yMat, beta0, maxiterBT, maxErr, lim, modify=modify);
    xEST = estimateX(xEST0, beta, yMat);
    
    MatrixXd normErr = ifElseMat(xEST0.array().abs(), lim, xEST0, lim);
    MatrixXd err = (xEST.array() - xEST0.array()).abs();
    err = err.array()/normErr.array();
    if(err.maxCoeff() < maxErr) break;
    xEST0 = xEST;
    MatrixXd betaT = beta.transpose();
    Eigen::Map<VectorXd> tempBeta(betaT.data(), betaT.size());
    beta0 = tempBeta;
  }
  List res = List::create(Named("X") = xEST , _["BETA"] = beta);
  return(res);
}

// [[Rcpp::export]]
MatrixXd getSECpp(MatrixXd X, MatrixXd b, MatrixXd y){
  MatrixXd yHat = X * b.transpose();
  MatrixXd err = y.array() - yHat.array();
  int yHatCol = yHat.cols();
  int yHatRow = yHat.rows();

  MatrixXd W(yHatRow, yHatCol);
  for(int i = 0 ; i < yHatRow ; i++){
    for(int j = 0 ; j < yHatCol ; j++){
      W(i, j) = 1/(yHat(i, j) + 1);
    }
  }

  int wCol = W.cols();
  MatrixXd SE(wCol, X.cols());
  for(int i = 0 ; i < wCol ; i ++){
    MatrixXd WX = X.array().colwise() * W.col(i).array();
    MatrixXd XWX = X.transpose() * WX;
    MatrixXd XWXi = XWX.inverse();
    VectorXd err2 = err.col(i).array() * err.col(i).array();
    MatrixXd temp = WX.array().colwise() * err2.array();
    MatrixXd XWeWX = WX.transpose() * temp ;
    MatrixXd mat = XWXi * XWeWX * XWXi;
    VectorXd se = mat.diagonal().array().sqrt();
    SE.row(i) = se;
  }
  return(SE);
}

//###############################################################################

struct EMForMatrix : public Worker
{
   // source matrix
  const RMatrix<double> X;
  const RMatrix<double> y;
   
  // destination matrix
  RMatrix<double> out;
   
  EMForMatrix(const NumericMatrix X, const NumericMatrix y, NumericMatrix out) 
      : X(X), y(y), out(out) {}

  void operator()(std::size_t begin, std::size_t end) {
    // x: 13 5
    // y: 13 20
    // out: 5 20
    // i : 0 - 20
    int outRow = out.nrow();
    int yRow = y.nrow();
    int xRow = X.nrow();
    int xCol = X.ncol();

    MatrixXd Xs(xRow, xCol);
    for(unsigned i = 0 ; i < xRow ; i++){
      for(unsigned j = 0 ; j < xCol ; j++){
        Xs(i, j) = X(i, j);
      }
    }
    for(unsigned i = begin ; i < end ; i++){
      VectorXd ys(yRow);
      for(int j = 0 ; j < yRow ; j++){
        ys(j) = y(j, i);
      }

      VectorXd res = EMCpp(Xs, ys);
      for(unsigned j = 0 ; j < outRow ; j++){
        out(j, i) = res(j);
      }
     }
   }
};

// [[Rcpp::export]]
NumericMatrix parallelEMCpp(NumericMatrix X, NumericMatrix y) {

  NumericMatrix out(X.cols(), y.cols());

  EMForMatrix em(X, y, out);

  parallelFor(0, y.cols(), em);
  
  return out;
}

//#############################################################################
// currently using mcaplly, have not solve thread safety problem for List
// // [[Rcpp::export]]
// Rcpp::List parallelAEM(Rcpp::List x, Rcpp::List y, int nThread=2){
//   Rcpp::List res(x.size());
//   int n = res.size();
//   // omp_set_num_threads(nThread) ;
//   // omp_set_dynamic(0);         // Explicitly disable dynamic teams
//   omp_set_num_threads(2);
//   #pragma omp parallel for 
//   for(int i = 0 ; i < x.size() ; i++){
//     // auto cpu = sched_getcpu();
//     // cout << cpu << endl;
//     MatrixXd xs = Rcpp::as<MatrixXd>(x[i]);
//     MatrixXd ys = Rcpp::as<MatrixXd>(y[i]);
//     cout << xs << endl;
//     Rcpp::List r = AEMCpp(xs, ys);
//     MatrixXd temp = r[0];
//     cout << temp << endl;
//     // res(i) = r;
//   }
//   return(res);
// }