
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "SAIGE.hpp"
#include "DenseGRM.hpp"
#include "UTIL.hpp"

#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds

namespace POLMM {

POLMMClass::POLMMClass(arma::mat t_muMat,
                       arma::mat t_iRMat,
                       arma::mat t_Cova,
                       arma::uvec t_yVec,
                       arma::sp_mat t_SparseGRM,
                       double t_tau,
                       bool t_printPCGInfo,
                       double t_tolPCG,
                       int t_maxiterPCG,
                       double t_varRatio, 
                       double t_SPA_Cutoff,
                       bool t_flagSparseGRM)     // In region-based analysis, we use SparseGRM, in marker-based analysis, we do not use SparseGRM
{
  m_muMat = t_muMat;
  m_iRMat = t_iRMat;
  m_Cova = t_Cova;
  m_varRatio = t_varRatio;
  m_SPA_Cutoff = t_SPA_Cutoff;
  m_flagSparseGRM = t_flagSparseGRM;
  
  m_SparseGRM = t_SparseGRM;
  m_printPCGInfo = t_printPCGInfo;
  m_tolPCG = t_tolPCG;
  m_maxiterPCG = t_maxiterPCG;
  
  m_n = m_muMat.n_rows;
  m_J = m_muMat.n_cols;
  m_p = m_Cova.n_cols;
  
  m_CovaMat = getCovaMat(m_Cova, m_J);       // n(J-1) x p
  
  m_tau = t_tau;
  
  // std::cout << "test0" << std::endl;
  // std::this_thread::sleep_for (std::chrono::seconds(1));
  
  if(m_flagSparseGRM == true){
    
    // std::cout << "test01" << std::endl;
    // std::this_thread::sleep_for (std::chrono::seconds(1));
    
    m_InvBlockDiagSigma = getInvBlockDiagSigma();
  }
  
  // std::cout << "test1" << std::endl;
  // std::this_thread::sleep_for (std::chrono::seconds(1));
  
  // output for Step 2
  arma::mat XR_Psi_R(m_p, m_n * (m_J-1));                // p x n(J-1)
  for(int k = 0; k < m_p; k++){
    arma::mat xMat = Vec2Mat(m_CovaMat.col(k), m_n, m_J);
    arma::vec temp = Mat2Vec(getPsixMat(xMat / m_iRMat) / m_iRMat, m_n, m_J);
    XR_Psi_R.row(k) = temp.t();
  }
  
  m_XXR_Psi_RX = m_Cova * inv(XR_Psi_R * m_CovaMat);               // (n x p) * (p x p) = n x p
  
  // sum each (J-1) rows to 1 row: p x n(J-1) -> p x n
  m_XR_Psi_R = sumCols(XR_Psi_R, m_J);      // p x n
  
  m_yVec = t_yVec;
  arma::mat yMat(m_n, m_J, arma::fill::zeros);
  for(int i = 0; i < m_n; i++)
    yMat(i, m_yVec(i)) = 1;
  
  arma::mat ymuMat = yMat - m_muMat;                      // n x J
  arma::mat RymuMat = ymuMat.cols(0, m_J-2) / t_iRMat;    // n x (J-1): R %*% (y - mu)
  m_RymuVec = sumCols(RymuMat, m_J);                      // n x 1
  
  // std::cout << "test2" << std::endl;
  // std::this_thread::sleep_for (std::chrono::seconds(1));
  
  if(m_flagSparseGRM == true){
    arma::mat iSigma_CovaMat(m_n * (m_J-1), m_p);
    getPCGofSigmaAndCovaMat(m_CovaMat, iSigma_CovaMat);
    arma::mat XSigmaX = inv(m_CovaMat.t() * iSigma_CovaMat);
    m_iSigmaX_XSigmaX = iSigma_CovaMat * XSigmaX;
  }
  
  setRPsiR();
}

// http://thecoatlessprofessor.com/programming/set_rs_seed_in_rcpp_sequential_case/
void set_seed(unsigned int seed){
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}

void POLMMClass::setPOLMMObj(bool t_flagSparseGRM,       // if 1, then use SparseGRM, otherwise, use DenseGRM
                             DenseGRM::DenseGRMClass* t_ptrDenseGRMObj,
                             PLINK::PlinkClass* t_ptrPlinkObj,
                             arma::mat t_Cova,
                             arma::uvec t_yVec,     // should be from 0 to J-1
                             arma::vec t_beta,
                             arma::vec t_bVec,
                             arma::vec t_eps,           // 
                             double t_tau,
                             arma::sp_mat t_SparseGRM,    // results of function getKinMatList()
                             Rcpp::List t_controlList)
{
  setControlList(t_controlList);
  setPOLMMInner(t_Cova, t_yVec, t_beta,  t_bVec,  t_eps,  t_tau);
  
  m_ptrPlinkObj = t_ptrPlinkObj;
  
  // if t_flagSparseGRM = 1, then use "SparseGRM" methods, otherwise, use "DenseGRM" methods
  m_flagSparseGRM = t_flagSparseGRM;
  if(m_flagSparseGRM){
    m_SparseGRM = t_SparseGRM; // edited on 03-27-2021, update later
    m_ZMat_sp = setZMat_sp();
    m_M = 0;
  }else{
    m_M = t_ptrDenseGRMObj->getM();
    m_ptrDenseGRMObj = t_ptrDenseGRMObj;
  }

  getTraceRandMat();
}

void POLMMClass::getTraceRandMat()
{
  arma::vec t1  = getTime();
  for(unsigned int itrace = 0; itrace < m_tracenrun; itrace++)
  {
    arma::vec uVec = nb(m_n * (m_J-1));
    uVec = uVec * 2 - 1;
    m_TraceRandMat.col(itrace) = uVec;
    arma::vec ZuVec = ZMat(uVec);
    // m_V_TRM.col(itrace) = tZMat(getKinbVecPOLMM(ZuVec, "none"));
    
    arma::vec tempVec = getKinbVecPOLMM(ZuVec, "none");
    m_V_TRM.col(itrace) = tZMat(tempVec);
  }
  
  arma::vec t2  = getTime();
  std::string info = "calculate " + std::to_string(m_tracenrun) + " genKinbVec()";
  printTime(t1, t2, info);
}

void POLMMClass::setPOLMMInner(arma::mat t_Cova,
                               arma::uvec t_yVec,     // should be from 0 to J-1
                               arma::vec t_beta,
                               arma::vec t_bVec,
                               arma::vec t_eps,           // 
                               double t_tau)
{
  m_n = t_Cova.n_rows;
  m_p = t_Cova.n_cols;
  m_J = arma::max(t_yVec) + 1;
  
  m_CovaMat = getCovaMat(t_Cova, m_J);
  m_yVec = t_yVec;
  m_yMat = getyMat(t_yVec);
  
  m_Cova = t_Cova;
  m_beta = t_beta;
  m_bVec = t_bVec;
  m_eps = t_eps;
  m_tau = t_tau;
  
  setArray();
  set_seed(m_seed);
}

// This function only uses variance ratio and does not use sparse GRM
void POLMMClass::getMarkerPval(arma::vec t_GVec, 
                               double& t_Beta, 
                               double& t_seBeta, 
                               double& t_pval, 
                               double t_altFreq,
                               double& t_zScore)
{
  arma::vec adjGVec = getadjGFast(t_GVec);
  
  double Stat = getStatFast(adjGVec);
  arma::vec VarWVec = getVarWVec(adjGVec);
  double VarW = sum(VarWVec);
  double VarS = VarW * m_varRatio;
  
  double StdStat = std::abs(Stat) / sqrt(VarS);
  double pvalNorm = 2 * arma::normcdf(-1*StdStat);
  double pval = pvalNorm;
  
  arma::vec K1roots = {3, -3};
  if(StdStat > m_SPA_Cutoff){
    
    arma::uvec posG1 = arma::find(t_GVec != 0);
    double VarW1 = sum(VarWVec(posG1));
    double VarW0 = VarW - VarW1;
    double Ratio0 = VarW0 / VarW;
    
    Rcpp::List resSPA = MAIN_SPA(Stat, adjGVec, K1roots, VarS, VarW, Ratio0, posG1);
    pval = resSPA["pval"];
  }
  
  t_pval = pval;
  t_Beta = Stat / VarS;
  t_seBeta = t_Beta / StdStat;
  t_zScore = Stat / sqrt(VarS);
}

// This function should use sparse GRM since in region-based analysis
// since most of the variants are low-frequency variants or rare variants. 
void POLMMClass::getRegionPVec(arma::vec t_GVec, 
                               double& t_Stat,
                               double& t_Beta, 
                               double& t_seBeta, 
                               double& t_pval0, 
                               double& t_pval1,
                               arma::vec& t_P1Vec, 
                               arma::vec& t_P2Vec)
{
  arma::vec adjGVec = getadjGFast(t_GVec);
  double Stat = getStatFast(adjGVec);
  
  arma::vec ZPZ_adjGVec = get_ZPZ_adjGVec(adjGVec);
  double VarS = as_scalar(adjGVec.t() * ZPZ_adjGVec);
  double StdStat = std::abs(Stat) / sqrt(VarS);
  double pvalNorm = 2 * arma::normcdf(-1*StdStat);
  double pval = pvalNorm;
  
  arma::vec K1roots = {3, -3};
  if(StdStat > m_SPA_Cutoff){
    
    arma::vec VarWVec = getVarWVec(adjGVec);
    double VarW = sum(VarWVec);
    double VarS = VarW * m_varRatio;
    
    arma::uvec posG1 = arma::find(t_GVec != 0);
    double VarW1 = sum(VarWVec(posG1));
    double VarW0 = VarW - VarW1;
    double Ratio0 = VarW0 / VarW;
    
    Rcpp::List resSPA = MAIN_SPA(Stat, adjGVec, K1roots, VarS, VarW, Ratio0, posG1);
    pval = resSPA["pval"];
  }
  
  t_pval0 = pvalNorm;
  t_pval1 = pval;
  t_Stat = Stat;
  t_Beta = Stat / VarS;
  t_seBeta = t_Beta / StdStat;
  
  t_P1Vec = adjGVec;
  t_P2Vec = ZPZ_adjGVec;
  
  // getMarkerPval(t_GVec, t_Beta, t_seBeta, t_pval0, altFreq);
}

arma::vec POLMMClass::getadjGFast(arma::vec t_GVec)
{
  // To increase computational efficiency when lots of GVec elements are 0
  arma::vec XR_Psi_RG(m_p, arma::fill::zeros);
  for(int i = 0; i < m_n; i++){
    if(t_GVec(i) != 0){
      XR_Psi_RG += m_XR_Psi_R.col(i) * t_GVec(i);
    }
  }
  arma::vec adjGVec = t_GVec - m_XXR_Psi_RX * XR_Psi_RG;
  return adjGVec;
}

double POLMMClass::getStatFast(arma::vec t_adjGVec)         // n x 1
{
  double Stat = 0;
  for(int i = 0; i < m_n; i++){
    if(t_adjGVec(i) != 0){
      Stat += t_adjGVec(i) * m_RymuVec(i);
    }
  }
  return Stat;
}

arma::vec POLMMClass::get_ZPZ_adjGVec(arma::vec t_adjGVec)
{
  arma::vec adjGVecLong = Vec2LongVec(t_adjGVec, m_n, m_J);  // that is Z %*% adjGVec
  
  arma::vec iSigmaGVec(m_n * (m_J-1), arma::fill::zeros);
  
  getPCGofSigmaAndVector(adjGVecLong, iSigmaGVec);  // iSigmaGVec = Sigma^-1 %*% Z %*% adjGVec
  
  arma::vec PZ_adjGVec = iSigmaGVec - m_iSigmaX_XSigmaX * (m_CovaMat.t() * iSigmaGVec);
  arma::vec ZPZ_adjGVec = LongVec2Vec(PZ_adjGVec, m_n, m_J);
  
  return ZPZ_adjGVec;
}

// use PCG to calculate iSigma_xMat = Sigma^-1 %*% xMat
void POLMMClass::getPCGofSigmaAndCovaMat(arma::mat t_xMat,              // matrix with dim of n(J-1) x p
                                         arma::mat& t_iSigma_xMat)      // matrix with dim of n(J-1) x p
{
  int p1 = t_xMat.n_cols;
  for(int i = 0; i < p1; i++){
    arma::vec y1Vec = t_xMat.col(i);
    arma::vec iSigma_y1Vec = t_iSigma_xMat.col(i);
    getPCGofSigmaAndVector(y1Vec, iSigma_y1Vec);
    t_iSigma_xMat.col(i) = iSigma_y1Vec;
  }
}


// use PCG to calculate xVec = Sigma^-1 %*% yVec
void POLMMClass::getPCGofSigmaAndVector(arma::vec t_y1Vec,    // vector with length of n(J-1)
                                        arma::vec& t_xVec,    // vector with length of n(J-1)
                                        std::string t_excludechr)
{
  // if(m_flagSparseGRM){
  // setSigmaMat_sp();
  // t_xVec = spsolve(SigmaMat_sp, t_y1Vec);
  // }else{
  arma::mat xMat = convert2(t_xVec, m_n, m_J);
  arma::mat y1Mat = convert2(t_y1Vec, m_n, m_J);
  // r2Vec and z2Vec are for the current step; r1Vec and z1Vec are for the previous step
  unsigned int iter = 0;
  arma::mat r2Mat = y1Mat - getSigmaxMat(xMat, t_excludechr);  // n x (J-1): r0 = y1Mat- Sigma %*% xMat
  double meanL2 = sqrt(getInnerProd(r2Mat, r2Mat)) / sqrt(m_n * (m_J-1));
  if(meanL2 <= m_tolPCG){
    // do nothing, xMat is already close to (Sigma)^-1 %*% y1Mat
  }else{
    iter++;
    arma::cube InvBlockDiagSigma = getInvBlockDiagSigma();
    arma::mat z2Mat = solverBlockDiagSigma(InvBlockDiagSigma, r2Mat);
    //
    arma::mat z1Mat, r1Mat;
    double beta1 = 0;
    arma::mat pMat = z2Mat;
    arma::mat ApMat = getSigmaxMat(pMat, t_excludechr);
    double alpha = getInnerProd(z2Mat, r2Mat) / getInnerProd(pMat, ApMat);
    xMat = xMat + alpha * pMat;
    r1Mat = r2Mat;
    z1Mat = z2Mat;
    r2Mat = r1Mat - alpha * ApMat;
    
    meanL2 = sqrt(getInnerProd(r2Mat, r2Mat)) / sqrt(m_n * (m_J-1));
    
    while (meanL2 > m_tolPCG && iter < m_maxiterPCG){
      iter++;
      
      //  z2Mat = minvMat % r2Mat;
      z2Mat = solverBlockDiagSigma(InvBlockDiagSigma, r2Mat);
      //
      beta1 = getInnerProd(z2Mat, r2Mat) / getInnerProd(z1Mat, r1Mat);
      pMat = z2Mat + beta1 * pMat;
      ApMat = getSigmaxMat(pMat, t_excludechr);
      alpha = getInnerProd(z2Mat, r2Mat) / getInnerProd(pMat, ApMat);
      
      xMat = xMat + alpha * pMat;
      r1Mat = r2Mat;
      z1Mat = z2Mat;
      r2Mat = r1Mat - alpha * ApMat;
      meanL2 = sqrt(getInnerProd(r2Mat, r2Mat)) / sqrt(m_n * (m_J-1));
    }
  }
  
  t_xVec = convert1(xMat, m_n, m_J);
  if (iter >= m_maxiterPCG){
    std::cout << "pcg did not converge. You may increase maxiter number." << std::endl;
  }
  if(m_showInfo)
    std::cout << "iter from getPCG1ofSigmaAndVector " << iter << std::endl; 
  // }
}

// use PCG to calculate xVec = Sigma^-1 %*% yVec
void POLMMClass::getPCGofSigmaAndVector(arma::vec t_y1Vec,    // vector with length of n(J-1)
                                        arma::vec& t_xVec)    // vector with length of n(J-1)
{
  arma::mat xMat = Vec2Mat(t_xVec, m_n, m_J);
  arma::mat y1Mat = Vec2Mat(t_y1Vec, m_n, m_J);
  
  // r2Vec and z2Vec are for the current step; r1Vec and z1Vec are for the previous step
  unsigned int iter = 0;
  arma::mat r2Mat = y1Mat - getSigmaxMat(xMat);  // n x (J-1): r0 = y1Mat- Sigma %*% xMat
  double meanL2 = sqrt(getInnerProd(r2Mat, r2Mat)) / sqrt(m_n * (m_J-1));
  if(meanL2 <= m_tolPCG){
    // do nothing, xMat is already close to (Sigma)^-1 %*% y1Mat
  }else{
    
    iter++;
    arma::mat z2Mat = solverBlockDiagSigma(r2Mat);
    //
    arma::mat z1Mat, r1Mat;
    double beta1 = 0;
    arma::mat pMat = z2Mat;
    arma::mat ApMat = getSigmaxMat(pMat);
    double alpha = getInnerProd(z2Mat, r2Mat) / getInnerProd(pMat, ApMat);
    xMat = xMat + alpha * pMat;
    r1Mat = r2Mat;
    z1Mat = z2Mat;
    r2Mat = r1Mat - alpha * ApMat;
    
    meanL2 = sqrt(getInnerProd(r2Mat, r2Mat)) / sqrt(m_n * (m_J-1));
    
    while (meanL2 > m_tolPCG && iter < m_maxiterPCG){
      
      iter++;
      //  z2Mat = minvMat % r2Mat;
      z2Mat = solverBlockDiagSigma(r2Mat);
      //
      beta1 = getInnerProd(z2Mat, r2Mat) / getInnerProd(z1Mat, r1Mat);
      pMat = z2Mat + beta1 * pMat;
      ApMat = getSigmaxMat(pMat);
      alpha = getInnerProd(z2Mat, r2Mat) / getInnerProd(pMat, ApMat);
      xMat = xMat + alpha * pMat;
      r1Mat = r2Mat;
      z1Mat = z2Mat;
      r2Mat = r1Mat - alpha * ApMat;
      meanL2 = sqrt(getInnerProd(r2Mat, r2Mat)) / sqrt(m_n * (m_J-1));
    }
  }
  
  t_xVec = Mat2Vec(xMat, m_n, m_J);
  if (iter >= m_maxiterPCG){
    std::cout << "pcg did not converge. You may increase maxiter number." << std::endl;
  }
  if(m_printPCGInfo)
    std::cout << "iter from getPCG1ofSigmaAndVector " << iter << std::endl; 
  // }
}

// yMat = Sigma %*% xMat
arma::mat POLMMClass::getSigmaxMat(arma::mat& t_xMat)   // matrix: n x (J-1) 
{
  arma::mat iR_xMat = m_iRMat % t_xMat;
  arma::mat iPsi_iR_xMat = getiPsixMat(iR_xMat);
  arma::mat yMat = m_iRMat % iPsi_iR_xMat;
  if(m_tau == 0){}
  else{
    arma::vec tZ_xMat = arma::sum(t_xMat, 1);  // rowSums(xMat): n x 1
    arma::vec V_tZ_xMat = m_SparseGRM * tZ_xMat;
    yMat.each_col() += m_tau * V_tZ_xMat;
  }
  return yMat;
}

// outMat = iPsiMat %*% xMat, iPsiMat is determined by muMat
arma::mat POLMMClass::getiPsixMat(arma::mat t_xMat)   // matrix with dim of n x (J-1)
{
  arma::mat iPsi_xMat(m_n, m_J-1);
  for(int i = 0; i < m_n; i ++){   // loop for samples
    double sumx = arma::sum(t_xMat.row(i));
    for(int j = 0; j < m_J-1; j++){
      iPsi_xMat(i,j) = sumx / m_muMat(i, m_J-1) + t_xMat(i,j) / m_muMat(i,j);
    }
  }
  return iPsi_xMat;
}

// outMat = PsiMat %*% xMat, PsiMat is determined by muMat
arma::mat POLMMClass::getPsixMat(arma::mat t_xMat)   // matrix: n x (J-1)
{
  arma::mat Psi_xMat(m_n, m_J-1);
  // loop for samples
  for(int i = 0; i < m_n; i++){
    arma::rowvec muVec(m_J-1);
    for(int j = 0; j < m_J-1; j++){
      Psi_xMat(i,j) = m_muMat(i,j) * t_xMat(i,j);
      muVec(j) = m_muMat(i,j);
    }
    double sum_mu_x = sum(Psi_xMat.row(i));
    Psi_xMat.row(i) -= muVec * sum_mu_x; 
  }
  return Psi_xMat;
}

// // used in getPCGofSigmaAndVector()
// arma::cube POLMMClass::getInvBlockDiagSigma()
// {
//   // get diagonal elements of GRM
//   arma::vec DiagGRM;
//   DiagGRM = m_tau * m_SparseGRM.diag();
//   
//   arma::cube InvBlockDiagSigma(m_J-1, m_J-1, m_n, arma::fill::zeros);
//   for(int i = 0; i < m_n; i++){
//     for(int j2 = 0; j2 < m_J-1; j2++){
//       for(int j1 = 0; j1 < m_J-1; j1++){
//         double temp = m_iRMat(i,j2) * (1 / m_muMat(i, m_J-1)) * m_iRMat(i,j1) + DiagGRM(i);
//         if(j2 == j1){
//           temp += m_iRMat(i,j2) * (1 / m_muMat(i,j2)) * m_iRMat(i,j1); 
//         }
//         InvBlockDiagSigma(j2, j1, i) = temp;
//       }
//     }
//     InvBlockDiagSigma.slice(i) = inv(InvBlockDiagSigma.slice(i));
//   }
//   return InvBlockDiagSigma;
// }

// used in getPCGofSigmaAndVector()
arma::cube POLMMClass::getInvBlockDiagSigma()
{
  
  // std::cout << "test011" << std::endl;
  // std::cout << m_flagSparseGRM << std::endl;
  // std::this_thread::sleep_for (std::chrono::seconds(1));
  
  arma::vec DiagGRM;
  if(m_flagSparseGRM){
    // get diagonal elements of GRM
    
    // std::cout << "test02" << std::endl;
    // std::this_thread::sleep_for (std::chrono::seconds(1));
    
    DiagGRM = m_tau * m_SparseGRM.diag();
    
    // std::cout << "test03" << std::endl;
    // std::this_thread::sleep_for (std::chrono::seconds(1));
    
  }else{
    arma::vec* pDiagStdGeno = m_ptrDenseGRMObj->getDiagStdGeno();
    DiagGRM = m_tau * (*pDiagStdGeno);   // n x 1
  }
  
  double temp;
  //
  arma::cube InvBlockDiagSigma(m_J-1, m_J-1, m_n, arma::fill::zeros);
  for(int i = 0; i < m_n; i++){
    for(int j2 = 0; j2 < m_J-1; j2++){
      for(int j1 = 0; j1 < m_J-1; j1++){
        temp = m_iRMat(i,j2) * (1 / m_muMat(i, m_J-1)) * m_iRMat(i,j1) + DiagGRM(i);
        if(j2 == j1){
          temp += m_iRMat(i,j2) * (1 / m_muMat(i,j2)) * m_iRMat(i,j1); 
        }
        InvBlockDiagSigma(j2, j1, i) = temp;
      }
    }
    InvBlockDiagSigma.slice(i) = inv(InvBlockDiagSigma.slice(i));
  }
  return InvBlockDiagSigma;
}

arma::mat POLMMClass::solverBlockDiagSigma(arma::mat& t_xMat)     // n x (J-1)
{
  arma::mat outMat(m_n, m_J-1);
  for(int i = 0; i < m_n; i++){
    outMat.row(i) = t_xMat.row(i) * m_InvBlockDiagSigma.slice(i); // could invert matrix?? be careful!
  }
  return outMat;
}

Rcpp::List POLMMClass::MAIN_SPA(double t_Stat,
                                arma::vec t_adjGVec,
                                arma::vec t_K1roots,
                                double t_VarP,
                                double t_VarW,
                                double t_Ratio0,
                                arma::uvec t_posG1)
{
  // std::cout << "t_VarP:\t" << t_VarP << std::endl;
  // std::cout << "t_VarW:\t" << t_VarW << std::endl;
  
  Rcpp::List resSPA = fastSaddle_Prob(t_Stat, t_VarP, t_VarW, t_Ratio0, t_K1roots,
                                      t_adjGVec.elem(t_posG1), m_muMat.rows(t_posG1), m_iRMat.rows(t_posG1));
  return resSPA;
}

void POLMMClass::setRPsiR()   
{
  // arma::cube RPsiR(J-1, J-1, n, arma::fill::zeros);
  arma::vec RPsiRVec(m_n, arma::fill::zeros);
  arma::mat muRMat = m_muMat.cols(0, m_J-2) / m_iRMat;
  for(int i = 0; i < m_n; i++){
    for(int j1 = 0; j1 < m_J-1; j1++){
      // RPsiR(j1,j1,i) += muRMat(i,j1) / iRMat(i,j1) - muRMat(i,j1) * muRMat(i,j1);
      RPsiRVec(i) += muRMat(i,j1) / m_iRMat(i,j1) - muRMat(i,j1) * muRMat(i,j1);
      for(int j2 = j1+1; j2 < m_J-1; j2++){
        // RPsiR(j1,j2,i) -= muRMat(i,j1) * muRMat(i,j2);
        RPsiRVec(i) -= 2* muRMat(i,j1) * muRMat(i,j2);
      }
    }
  }
  m_RPsiR = RPsiRVec;
}

arma::vec POLMMClass::getVarWVec(arma::vec adjGVec)
{
  arma::vec VarWVec = m_RPsiR % pow(adjGVec, 2);
  return VarWVec;
}

double K0(double t_x,
          arma::mat t_muMat,     // N x (J-1)
          arma::mat t_cMat,      // N x (J-1)
          double t_m1)           // sum(muMat * cMat)
{
  arma::mat temp1Mat = - t_muMat + t_muMat % exp(t_cMat * t_x);
  arma::vec temp1Vec = log(1 + arma::sum(temp1Mat, 1));   // arma::sum(Mat, 1) is rowSums()
  double y = sum(temp1Vec) - t_m1 * t_x;
  
  // std::cout << "sum(temp1Vec):\t" << sum(temp1Vec) << std::endl;
  // std::cout << "t_m1:\t" << t_m1 << std::endl;
  // std::cout << "t_x:\t" << t_x << std::endl;
  
  return y;
}

arma::vec K12(double t_x,
              arma::mat t_muMat,
              arma::mat t_cMat,
              double t_m1)
{
  arma::mat temp0Mat = t_muMat % exp(t_cMat * t_x);
  arma::mat temp1Mat = - t_muMat + temp0Mat;
  arma::mat temp2Mat = temp0Mat % t_cMat;
  arma::mat temp3Mat = temp2Mat % t_cMat;
  
  arma::vec temp1Vec = 1 + arma::sum(temp1Mat, 1);
  arma::vec temp2Vec = arma::sum(temp2Mat, 1);
  arma::vec temp3Vec = arma::sum(temp3Mat, 1);
  
  arma::vec yVec(2);
  yVec(0) = sum(temp2Vec / temp1Vec) - t_m1;
  // yMat[i,2] = sum((temp3Vec*temp1Vec-temp2Vec^2)/temp1Vec^2, na.rm=TRUE);
  yVec(1) = sum((temp3Vec % temp1Vec - pow(temp2Vec, 2)) / pow(temp1Vec, 2));
  
  return yVec;
}

Rcpp::List fastgetroot_K1(double t_Stat,
                          double t_initX,
                          double t_Ratio0,
                          arma::mat t_muMat,
                          arma::mat t_cMat,
                          double t_m1)
{
  double x = t_initX;
  double K1 = 0;
  double K2 = 0;
  double diffX = arma::datum::inf;
  bool converge = true;
  double tol = 0.0001;
  int maxiter = 100;
  int iter = 0;
  
  for(iter = 0; iter < maxiter; iter ++){
    double oldX = x;
    double oldDiffX = diffX;
    double oldK1 = K1;
    
    arma::vec K12Vec = K12(x, t_muMat, t_cMat, t_m1);
    
    K1 = K12Vec(0) - t_Stat + t_Ratio0 * x;
    K2 = K12Vec(1) + t_Ratio0;
    
    diffX = -1 * K1 / K2;
    
    // std::cout << "iter:\t" << iter << std::endl;
    // std::cout << "x:\t" << x << std::endl;
    // std::cout << "K12Vec(1):\t" << K12Vec(1) << std::endl;
    // std::cout << "t_Ratio0:\t" << t_Ratio0 << std::endl;
    // std::cout << "K1:\t" << K1 << std::endl;
    // std::cout << "K2:\t" << K2 << std::endl;
    // std::cout << "diffX:\t" << diffX << std::endl;
    
    if(!std::isfinite(K1)){
      // checked it on 07/05:
      // if the solution 'x' tends to infinity, 'K2' tends to 0, and 'K1' tends to 0 very slowly.
      // then we can set the one sided p value as 0 (instead of setting converge = F)
      x = arma::sign(t_Stat) * arma::datum::inf;
      K2 = 0;
      break;
    }
    
    if(arma::sign(K1) != arma::sign(oldK1)){
      while(std::abs(diffX) > std::abs(oldDiffX) - tol){
        diffX = diffX / 2;
      }
    }
    if(std::abs(diffX) < tol) break;
    
    x = oldX + diffX;
    
  }
  
  if(iter == maxiter - 1) 
    converge = false;
  
  Rcpp::List yList = Rcpp::List::create(Rcpp::Named("root") = x,
                                        Rcpp::Named("iter") = iter,
                                        Rcpp::Named("converge") = converge,
                                        Rcpp::Named("K2") = K2);
  return yList;
}

double fastGet_Saddle_Prob(double t_Stat,
                           double t_zeta,
                           double t_K2,
                           double t_Ratio0,
                           arma::mat t_muMat,
                           arma::mat t_cMat,
                           double t_m1,          // sum(muMat * cMat)
                           bool t_lowerTail)
{
  double k1 = K0(t_zeta, t_muMat, t_cMat, t_m1) + 0.5 * pow(t_zeta, 2) * t_Ratio0;
  double k2 = t_K2;
  double pval = 0;
  if(std::isfinite(k1) && std::isfinite(k2))
  {
    double w = arma::sign(t_zeta) * sqrt(2 * (t_zeta * t_Stat - k1));
    double v = t_zeta * sqrt(t_K2);
    
    double Z = w + 1/w * log(v/w);
    pval = arma::normcdf(arma::sign(t_lowerTail-0.5) * Z);
  }
  
  return pval;
}

// add partial normal approximation to speed up the SPA
Rcpp::List fastSaddle_Prob(double t_Stat,
                           double t_VarP,
                           double t_VarW,
                           double t_Ratio0,      // Ratio of variance (G==0)
                           arma::vec t_K1roots,  // 2 x 1
                           arma::vec t_adjGVec1, // N1 x 1, where N1 is length(G!=0)
                           arma::mat t_muMat1,   // N1 x J
                           arma::mat t_iRMat1)   // N1 x (J-1)
{
  // std::cout<< "Start fastSaddle_Prob()....." << std::endl;
  
  int J = t_muMat1.n_cols;
  int N1 = t_muMat1.n_rows;
  
  t_muMat1 = t_muMat1.cols(0, J-2);
  double adjStat = t_Stat / sqrt(t_VarP);
  
  double sqrtVarW = sqrt(t_VarW);
  
  arma::mat cMat(N1, J-1);
  for(int i = 0; i < N1; i ++){
    for(int j = 0; j < J-1; j ++){
      cMat(i,j) = t_adjGVec1(i) / t_iRMat1(i,j) / sqrtVarW;
    }
  }
  
  double m1 = arma::accu(t_muMat1 % cMat);
  
  Rcpp::List outUni1 = fastgetroot_K1(std::abs(adjStat), std::min(t_K1roots(0), 5.0), 
                                      t_Ratio0, t_muMat1, cMat, m1);
  Rcpp::List outUni2 = fastgetroot_K1(-1 * std::abs(adjStat), std::max(t_K1roots(1), -5.0), 
                                      t_Ratio0, t_muMat1, cMat, m1);
  
  bool converge = false;
  double pval = 0; 
  arma::vec K1roots;
  
  bool outUnit1Converge = outUni1["converge"];
  bool outUnit2Converge = outUni2["converge"];
  
  // std::cout << "adjStat:\t" << adjStat << std::endl;
  // std::cout << "t_Ratio0:\t" << t_Ratio0 << std::endl;
  // double root1 = outUni1["root"];
  // double root2 = outUni2["root"];
  // std::cout << "outUni1.root:\t" << root1 << std::endl;
  // std::cout << "outUni2.root:\t" << root2 << std::endl;
    
  if(outUnit1Converge == true && outUnit2Converge == true){
    
    double p1 = fastGet_Saddle_Prob(std::abs(adjStat), outUni1["root"], 
                                    outUni1["K2"], t_Ratio0, t_muMat1, cMat, m1, false);
    
    // double root = outUni1["root"];
    // double K2 = outUni1["K2"];
    
    // std::cout << "outUni1:\t" << root << "\t" << K2 << std::endl;
    // std::cout << "p1:\t" << p1 << std::endl;
    
    double p2 = fastGet_Saddle_Prob(-1 * std::abs(adjStat), outUni2["root"], 
                                    outUni2["K2"], t_Ratio0, t_muMat1, cMat, m1, true);
    
    // root = outUni2["root"];
    // K2 = outUni2["K2"];
    // std::cout << "outUni2:\t" << root << "\t" << K2 << std::endl;
    // std::cout << "p2:\t" << p2 << std::endl;
    
    pval = p1 + p2;
    
    converge = true;
    K1roots = {outUni1["root"], outUni2["root"]};
  }else{
    std::cout << "SPA does not converge, use normal approximation p value." << std::endl;
    pval = 2 * arma::normcdf(-1 * std::abs(adjStat));
    K1roots = t_K1roots;
  }
  
  if((!std::isfinite(pval)) || (pval == 0)){
    std::cout << "SPA does not give a valid p value, use normal approximation p value." << std::endl;
    pval = 2 * arma::normcdf(-1 * std::abs(adjStat));
    K1roots = t_K1roots;
  }
  
  Rcpp::List yList = Rcpp::List::create(Rcpp::Named("pval") = pval,
                                        Rcpp::Named("converge") = converge,
                                        Rcpp::Named("K1roots") = K1roots);
  return yList;
}

void POLMMClass::setSeqMat(int t_NonZero_cutoff)         // number of subjects, if the number of subjects is less than or equal to this value, we use ER
{
  int n = t_NonZero_cutoff;
  std::cout << "Setting m_SeqMat for Efficient Resampling (ER)...." << std::endl;
  uint32_t nER = pow(m_J, n);       // J^n
  arma::Col<uint32_t> y = arma::linspace<arma::Col<uint32_t>>(0, nER-1, nER);  // nER x 1 matrix: seq(0, nER-1, 1)
  m_SeqMat.resize(n, nER);
  uint32_t powJ = nER / m_J;         // J^(n-1)
  arma::Col<uint8_t> SeqVec(nER);
  for(int i = 0; i < n; i++){
    int pos_row = n - 1 - i;
    for(uint32_t j = 0; j < nER; j++){
      SeqVec(j) = y(j) / powJ;
      y(j) = y(j) - SeqVec(j) * powJ;
    }
    powJ = powJ / m_J;           // J^(n-2), ..., J^0
    m_SeqMat.row(pos_row) = SeqVec.t();
  }
}

// arma::umat updateSeqMat(arma::umat t_SeqMat, // n x J^n matrix
//                         int t_n1,            // number of subjects, should be < n
//                         int t_J)             // number of levels
// {
//   int nER = pow(t_J, t_n1);
//   arma::umat PartSeqMat = t_SeqMat.submat(0, 0, t_n1-1, nER-1);
//   return PartSeqMat;
// }

double POLMMClass::MAIN_ER(arma::vec t_GVec,
                           arma::uvec t_posG1)
{
  int N1 = t_posG1.size();
  uint32_t nER = pow(m_J, N1);
  arma::Mat<uint8_t> SeqMat = m_SeqMat.submat(0, 0, N1-1, nER-1);
  double pvalER = getPvalER(m_yVec.elem(t_posG1), t_GVec.elem(t_posG1), m_muMat.rows(t_posG1), m_iRMat.rows(t_posG1), SeqMat);
  
  return pvalER;
}

// Main function: note that n is the number of subjects with Geno != 0
double getPvalER(arma::uvec t_yVec,     // N1 x 1 vector, from 0 to J-1
                 arma::vec t_GVec,      // N1 x 1 vector,
                 arma::mat t_muMat,     // N1 x J matrix,
                 arma::mat t_iRMat,     // N1 x (J-1) matrix
                 arma::Mat<uint8_t> t_SeqMat) // N1 x nER
{
  uint32_t nER = t_SeqMat.n_cols;
  arma::vec StatVec = getStatVec(t_SeqMat, t_GVec, t_muMat, t_iRMat);
  
  arma::Col<uint8_t> yVec = arma::conv_to<arma::Col<uint8_t>>::from(t_yVec);
  double StatObs = arma::as_scalar(getStatVec(yVec, t_GVec, t_muMat, t_iRMat));
  
  double eps = 1e-10;
  
  double pvalER_pos = 0;
  double pvalER_neg = 0;
  double absStatObs = std::abs(StatObs);
  for(uint32_t i = 0; i < nER; i++){
    // double absStatTmp = std::abs(StatVec(i));
    // if(absStatObs < absStatTmp - eps){
    //   pvalER += getProb(t_SeqMat.col(i), t_muMat);
    // }else if(absStatObs < absStatTmp + eps){
    //   pvalER += 0.5 * getProb(t_SeqMat.col(i), t_muMat);
    // }
    double StatTmp = StatVec(i);
    
    if(StatTmp > absStatObs + eps){
      pvalER_pos += getProb(t_SeqMat.col(i), t_muMat);
    }else if(StatTmp > absStatObs - eps){
      pvalER_pos += 0.5 * getProb(t_SeqMat.col(i), t_muMat);
    }
    
    if(StatTmp < -1 * absStatObs - eps){
      pvalER_neg += getProb(t_SeqMat.col(i), t_muMat);
    }else if(StatTmp < -1 * absStatObs + eps){
      pvalER_neg += 0.5 * getProb(t_SeqMat.col(i), t_muMat);
    }
  }
  
  // std::cout << "pvalER_pos:\t" << pvalER_pos << std::endl;
  // std::cout << "pvalER_neg:\t" << pvalER_neg << std::endl;
  double pvalER = pvalER_pos + pvalER_neg;
  return pvalER;
}

arma::vec getStatVec(arma::Mat<uint8_t> t_SeqMat,   // n x J^n matrix
                     arma::vec t_GVec,      // n x 1 vector, where n is number of subjects with Geno != 0
                     arma::mat t_muMat,     // n x J matrix, where n is number of subjects with Geno != 0
                     arma::mat t_iRMat)     // n x (J-1) matrix
{
  int n = t_muMat.n_rows;
  int J = t_muMat.n_cols;
  int nER = t_SeqMat.n_cols;
  
  arma::vec StatVec(nER);
  
  arma::mat A(n, J-1);
  for(int i = 0; i < J-1; i++){
    A.col(i) = t_GVec / t_iRMat.col(i);
  }
  
  double a1 = arma::accu(A % t_muMat.cols(0, J-2));
  
  for(int i = 0; i < nER; i++){
    double a2 = 0;
    for(int j = 0; j < n; j++){
      int idxL = t_SeqMat(j, i); // from 0 to J-1, level index
      if(idxL != J-1){
        a2 += A(j,idxL);
      }
    }
    
    StatVec(i) = a2 - a1;
  }
  
  return StatVec;
}

double getProbOne(arma::Col<uint8_t> t_SeqVec,  // n x 1
                  arma::mat t_muMat)    // n x J
{
  int n = t_muMat.n_rows;
  double tempProb = 1;
  for(int j = 0; j < n; j++){
    tempProb *= t_muMat(j, t_SeqVec(j));
  }
  return tempProb;
}


double getProb(arma::Mat<uint8_t> t_SeqMat,  // n x m matrix, where m \leq J^n is the number of resampling with abs(stat) > stat_obs
               arma::mat t_muMat)            // n x J matrix
{
  int nER = t_SeqMat.n_cols;
  
  double prob = 0;
  
  for(int i = 0; i < nER; i++){
    arma::Col<uint8_t> SeqVec = t_SeqMat.col(i);
    double tempProb = getProbOne(SeqVec, t_muMat);
    prob += tempProb;
  }
  
  return prob;
}

// update (m_eta, m_WMat, m_muMat, m_mMat, m_nuMat, m_iRMat, m_YMat) based on (m_beta, m_bVec, m_eps)
void POLMMClass::updateMats()
{
  // update (m_eta)
  m_eta = m_Cova * m_beta + m_bVec;
  
  // update (m_WMat, m_muMat, m_mMat, m_nuMat)
  double tmpExp, tmpnu0, tmpnu1;
  for(int i = 0; i < m_n; i ++){  // loop for samples
    tmpnu0 = 0;  // eps_0 = -Inf
    for(int j = 0; j < m_J-1; j ++){  // loop from eps_1 to eps_{J-1}
      tmpExp = exp(m_eps(j) - m_eta(i));
      tmpnu1 = tmpExp / (1 + tmpExp);
      m_muMat(i,j) = tmpnu1 - tmpnu0;
      m_WMat(i,j) = tmpnu1 * (1 - tmpnu1);
      m_mMat(i,j) = tmpnu1 + tmpnu0 - 1;
      m_nuMat(i,j) = tmpnu1;
      tmpnu0 = tmpnu1;
    }
    int j = m_J-1;      // eps_J = Inf
    tmpnu1 = 1;
    m_muMat(i,j) = tmpnu1 - tmpnu0;
    m_WMat(i,j) = tmpnu1 * (1 - tmpnu1);
    m_mMat(i,j) = tmpnu1 + tmpnu0 - 1;
    m_nuMat(i,j) = tmpnu1;
  }
  
  // update (iRMat)
  for(int i = 0; i < m_n; i ++){
    for(int j = 0; j < m_J-1; j ++){
      m_iRMat(i,j) = 1 / (m_mMat(i,j) - m_mMat(i, m_J-1));
    }
  }
  
  // update (YMat)
  arma::mat xMat = m_yMat.cols(0, m_J-2) - m_muMat.cols(0, m_J-2);
  arma::mat iPsi_xMat = getiPsixMat(xMat);
  for(int i = 0; i < m_n; i++){  // loop for samples
    for(int j = 0; j < m_J-1; j++)
      m_YMat(i,j) = m_eta(i) + (m_iRMat(i,j) * iPsi_xMat(i,j));
  }
}
// need update in case that eps(k+1) < eps(k)
void POLMMClass::updateEpsOneStep()
{
  // the first eps is fixed at 0
  arma::vec d1eps(m_J-2, arma::fill::zeros);
  arma::mat d2eps(m_J-2, m_J-2, arma::fill::zeros);
  double temp1, temp2, temp3;
  
  for(int k = 1; k < m_J-1; k++){
    for(int i = 0; i < m_n; i++){
      temp1 = m_yMat(i, k) / m_muMat(i, k) - m_yMat(i, k+1) / m_muMat(i, k+1);
      temp2 = - m_yMat(i, k) / m_muMat(i, k) / m_muMat(i, k) - m_yMat(i, k+1) / m_muMat(i, k+1) / m_muMat(i, k+1);
      
      d1eps(k - 1) += m_WMat(i, k) * temp1;
      d2eps(k - 1, k - 1) += m_WMat(i, k) * (1 - 2 * m_nuMat(i, k)) * temp1 + m_WMat(i, k) * m_WMat(i, k) * temp2;
      if(k < m_J-2){
        temp3 = m_WMat(i,k) * m_WMat(i, k+1) * m_yMat(i, k+1) / m_muMat(i, k+1) / m_muMat(i, k+1);
        d2eps(k-1, k) += temp3;
        d2eps(k, k-1) += temp3;
      }
    }
  }
  
  std::cout << "d2eps:\t" << d2eps << std::endl;
  std::cout << "d1eps:\t" << d1eps << std::endl;
  
  arma::vec deps = -1 * inv(d2eps) * d1eps;
  
  std::cout << "deps:\t" << deps << std::endl;
  
  for(int k = 1; k < m_J-1; k ++){
    m_eps(k) += deps(k-1);
  }
}

void POLMMClass::updateEps()
{
  for(unsigned int iter = 0; iter < m_maxiterEps; iter ++){
    
    arma::vec eps0 = m_eps;
    
    updateEpsOneStep();
    updateMats();
    
    std::cout << "iter:\t" << iter << std::endl;
    std::cout << "m_eps:\t" << m_eps << std::endl;
    std::cout << "eps0:\t" << eps0 << std::endl;
    std::cout << "m_tolEps:\t" << m_tolEps << std::endl;
    
    double diffeps = arma::max(arma::abs(m_eps - eps0)/(arma::abs(m_eps) + arma::abs(eps0) + m_tolEps));
    
    if(diffeps < m_tolEps){
      
      std::cout << "UpdateEps iter: " << iter << std::endl;
      std::cout << "eps: " << std::endl << m_eps << std::endl;
      break;
    }
  }
}

void POLMMClass::updatePara(std::string t_excludechr)
{
  getPCGofSigmaAndCovaMat(m_CovaMat, m_iSigma_CovaMat, t_excludechr);
  arma::vec YVec = convert1(m_YMat, m_n, m_J);
  getPCGofSigmaAndVector(YVec, m_iSigma_YVec, t_excludechr); 
  
  // update beta
  arma::mat XSigmaX = inv(m_CovaMat.t() * m_iSigma_CovaMat);
  arma::vec Cova_iSigma_YVec = m_CovaMat.t() * m_iSigma_YVec;
  m_beta = XSigmaX * Cova_iSigma_YVec;
  m_iSigmaX_XSigmaX = m_iSigma_CovaMat * XSigmaX;
  
  // update bVec
  arma::vec Z_iSigma_YVec = ZMat(m_iSigma_YVec);
  arma::vec Z_iSigma_Xbeta = ZMat(m_iSigma_CovaMat * m_beta);
  arma::vec tempVec = Z_iSigma_YVec - Z_iSigma_Xbeta;
  m_bVec = m_tau * getKinbVecPOLMM(tempVec, t_excludechr);
}

arma::mat POLMMClass::getVarRatio(arma::mat t_GMatRatio, std::string t_excludechr)
{
  std::cout << "Start estimating variance ratio...." << std::endl;
  Rcpp::List objP = getobjP(m_Cova, m_yMat, m_muMat, m_iRMat);
  
  arma::vec GVec(m_n);
  arma::rowvec VarOneSNP(5);
  
  arma::mat VarRatioMat(m_nSNPsVarRatio, 5);
  arma::mat newVarRatio(10, 5);
  
  unsigned int index = 0;
  int indexTot = 0;
  while(index < m_nSNPsVarRatio){
    GVec = t_GMatRatio.col(index);
    VarOneSNP = getVarOneSNP(GVec, t_excludechr, objP);
    VarRatioMat.row(index) = VarOneSNP;
    index++;
    indexTot++;
  }
  
  arma::vec VarRatio = VarRatioMat.col(4);
  double CV = calCV(VarRatio);
  std::cout << "nSNPs for CV: " << index << std::endl;
  std::cout << "CV: " << CV << std::endl;
  
  while(CV > m_CVcutoff && VarRatioMat.n_rows <= 100){
    int indexTemp = 0;
    while(indexTemp < 10){
      indexTot++;
      GVec = t_GMatRatio.col(indexTot);
      VarOneSNP = getVarOneSNP(GVec, t_excludechr, objP);
      newVarRatio.row(indexTemp) = VarOneSNP;
      index++;
      indexTemp++;
    }
    VarRatioMat.insert_rows(0, newVarRatio);
    arma::vec VarRatio = VarRatioMat.col(4);
    CV = calCV(VarRatio);
    std::cout << "nSNPs for CV: " << index << std::endl;
    std::cout << "CV: " << CV << std::endl;
  }
  return(VarRatioMat);
}

arma::rowvec POLMMClass::getVarOneSNP(arma::vec GVec,
                                      std::string excludechr,
                                      Rcpp::List objP)
{
  arma::rowvec VarOut(5);
  
  double AF = arma::sum(GVec) / GVec.size() / 2;
  
  if(AF > 0.5)
    AF = 1 - AF;
  
  Rcpp::List adjGList = outputadjGFast(GVec, objP);
  arma::vec adjGVec = adjGList["adjGVec"];
  double Stat = adjGList["Stat"];
  double VarW = adjGList["VarW"];
  double VarP = getVarP(adjGVec, excludechr);
  
  VarOut(0) = AF;
  VarOut(1) = Stat;
  VarOut(2) = VarW;
  VarOut(3) = VarP;
  VarOut(4) = VarP/VarW;
  return(VarOut);
}

// update parameters (except tau) until converge
void POLMMClass::updateParaConv(std::string t_excludechr)
{
  for(unsigned int iter = 0; iter < m_maxiter; iter ++){
    
    arma::vec beta0 = m_beta;
    
    // update beta and bVec
    updatePara(t_excludechr);
    updateMats();
    
    // update eps (cutpoints)
    updateEps();
    updateMats();
    
    std::cout << "beta: " << std::endl << m_beta << std::endl;
    
    double diffBeta = arma::max(arma::abs(m_beta - beta0)/(arma::abs(m_beta) + arma::abs(beta0) + m_tolBeta));
    std::cout << "diffBeta:\t" << diffBeta << std::endl << std::endl;
    if(diffBeta < m_tolBeta)
      break;
  }
}

void POLMMClass::updateTau()
{
  std::cout << "Start updating tau..." << std::endl;
  
  arma::vec YVec = convert1(m_YMat, m_n, m_J);
  getPCGofSigmaAndCovaMat(m_CovaMat, m_iSigma_CovaMat, "none");
  getPCGofSigmaAndVector(YVec, m_iSigma_YVec, "none"); 
  m_iSigmaX_XSigmaX = m_iSigma_CovaMat * inv(m_CovaMat.t() * m_iSigma_CovaMat);
  arma::vec PYVec = m_iSigma_YVec - m_iSigmaX_XSigmaX * (m_CovaMat.t() * m_iSigma_YVec);
  arma::vec ZPYVec = ZMat(PYVec);
  arma::vec VPYVec = tZMat(getKinbVecPOLMM(ZPYVec, "none"));
  
  getPCGofSigmaAndVector(VPYVec, m_iSigma_VPYVec, "none");
  
  arma::vec PVPYVec = m_iSigma_VPYVec - m_iSigmaX_XSigmaX * (m_CovaMat.t() * m_iSigma_VPYVec);
  double YPVPY = as_scalar(YVec.t() * PVPYVec);
  double YPVPVPY = as_scalar(VPYVec.t() * PVPYVec);
  // The below is to calculate trace
  
  getPCGofSigmaAndCovaMat(m_V_TRM, m_iSigma_V_TRM, "none");
  
  double tracePV = 0;
  int m = m_TraceRandMat.n_cols;
  
  for(int i = 0; i < m; i++){
    arma::vec iSigma_V_TRM_col = m_iSigma_V_TRM.col(i);
    arma::vec P_V_TRM_col = iSigma_V_TRM_col - m_iSigmaX_XSigmaX * (m_CovaMat.t() * iSigma_V_TRM_col);
    tracePV += as_scalar(m_TraceRandMat.col(i).t() * P_V_TRM_col);
  }
  tracePV /= m;
  // final step
  double deriv = 0.5 * YPVPY - 0.5 * tracePV;
  double AI = 0.5 * YPVPVPY;
  double dtau = deriv / AI;
  double tau0 = m_tau;
  m_tau = tau0 + dtau;
  while(m_tau < 0){
    dtau = dtau / 2;
    m_tau = tau0 + dtau;
  }
  if(m_tau < 1e-4){
    m_tau = 0;
  }
}

void POLMMClass::fitPOLMM()
{
  // initial vector
  arma::vec t1  = getTime();
  updateMats();
  
  // start iteration
  std::cout << "Start iteration ....." << std::endl;
  
  for(m_iter = 0; m_iter < m_maxiter; m_iter ++){
    
    // update fixed effect coefficients
    updateParaConv("none");
    
    // update tau
    double tau0 = m_tau;
    updateTau();
    
    if(std::isnan(m_tau))
      Rcpp::stop("Parameter tau is NA.");
    
    std::cout << "iter: " << m_iter << std::endl;
    std::cout << "beta: " << std::endl << m_beta << std::endl;
    std::cout << "tau: " << m_tau << std::endl << std::endl;
    
    double diffTau = std::abs(m_tau - tau0) / (std::abs(m_tau) + std::abs(tau0) + m_tolTau);
    
    if(diffTau < m_tolTau)
      break;
  }
  
  if(m_LOCO){
    
    // turn on LOCO option
    Rcpp::StringVector chrVec = m_ptrPlinkObj->getChrVec();
    Rcpp::StringVector uniqchr = unique(chrVec);
    
    std::cout << "uniqchr is " << uniqchr << std::endl;
    
    for(int i = 0; i < uniqchr.size(); i ++){
      
      std::string excludechr = std::string(uniqchr(i));
      std::cout << std::endl << "Leave One Chromosome Out: Chr " << excludechr << std::endl;
      
      updateParaConv(excludechr);
      
      arma::mat GMatRatio = m_ptrPlinkObj->getGMat(100, excludechr, m_minMafVarRatio, m_maxMissingVarRatio);
      arma::mat VarRatioMat = getVarRatio(GMatRatio, excludechr);
      double VarRatio = arma::mean(VarRatioMat.col(4));
      
      Rcpp::List temp = Rcpp::List::create(Rcpp::Named("muMat") = m_muMat,
                                           Rcpp::Named("iRMat") = m_iRMat,
                                           Rcpp::Named("VarRatioMat") = VarRatioMat,
                                           Rcpp::Named("VarRatio") = VarRatio);
      
      m_LOCOList[excludechr] = temp;
    }
    
  }else{
    
    // turn off LOCO option
    // if(!m_flagGMatRatio){
    
    arma::mat GMatRatio = m_ptrPlinkObj->getGMat(100, "none", m_minMafVarRatio, m_maxMissingVarRatio);
    arma::mat VarRatioMat = getVarRatio(GMatRatio, "none");
    double VarRatio = arma::mean(VarRatioMat.col(4));
    
    Rcpp::List temp = Rcpp::List::create(Rcpp::Named("muMat") = m_muMat,
                                         Rcpp::Named("iRMat") = m_iRMat,
                                         Rcpp::Named("VarRatioMat") = VarRatioMat,
                                         Rcpp::Named("VarRatio") = VarRatio);
    m_LOCOList["LOCO=F"] = temp;
  }
  
  // complete null POLMM fitting 
  arma::vec t2  = getTime();
  printTime(t1, t2, "fit the null POLMM.");
}



arma::mat POLMMClass::solverBlockDiagSigma(arma::cube& InvBlockDiagSigma,   // (J-1) x (J-1) x n
                                           arma::mat& xMat)                 // n x (J-1)
{
  arma::mat outMat(m_n, m_J-1);
  for(int i = 0; i < m_n; i++){
    outMat.row(i) = xMat.row(i) * InvBlockDiagSigma.slice(i); // could invert matrix?? be careful!
  }
  return(outMat);
}

// yMat = Sigma %*% xMat
arma::mat POLMMClass::getSigmaxMat(arma::mat t_xMat,   // matrix: n x (J-1) 
                                   std::string t_excludechr)
{
  arma::mat iR_xMat = m_iRMat % t_xMat;
  arma::mat iPsi_iR_xMat = getiPsixMat(iR_xMat);
  arma::mat yMat = m_iRMat % iPsi_iR_xMat;
  if(m_tau == 0){}
  else{
    arma::vec tZ_xMat = getRowSums(t_xMat);  // rowSums(xMat): n x 1
    arma::vec V_tZ_xMat = getKinbVecPOLMM(tZ_xMat, t_excludechr);
    yMat.each_col() += m_tau * V_tZ_xMat;
  }
  return(yMat);
}

arma::vec POLMMClass::getKinbVecPOLMM(arma::vec t_bVec, 
                                      std::string t_excludeChr)
{
  arma::vec KinbVec;
  
  if(m_flagSparseGRM){
    // arma::sp_mat temp = m_SparseGRM[t_excludeChr];
    // arma::sp_mat temp = m_SparseGRM;
    KinbVec = m_SparseGRM * t_bVec;
  }else{
    KinbVec = getKinbVec(t_bVec, m_ptrDenseGRMObj, t_excludeChr, m_grainSize);
  }
  Rcpp::checkUserInterrupt();
  return KinbVec;
}

// use PCG to calculate iSigma_xMat = Sigma^-1 %*% xMat
void POLMMClass::getPCGofSigmaAndCovaMat(arma::mat t_xMat,              // matrix with dim of n(J-1) x p
                                         arma::mat& t_iSigma_xMat,      // matrix with dim of n(J-1) x p
                                         std::string t_excludechr)
{
  int p1 = t_xMat.n_cols;
  for(int i = 0; i < p1; i++){
    
    arma::vec y1Vec = t_xMat.col(i);
    arma::vec iSigma_y1Vec = t_iSigma_xMat.col(i);
    getPCGofSigmaAndVector(y1Vec, iSigma_y1Vec, t_excludechr);
    
    t_iSigma_xMat.col(i) = iSigma_y1Vec;
  }
}

double POLMMClass::getVarP(arma::vec t_adjGVec,
                           std::string t_excludechr)
{
  arma::vec adjGVecLong = tZMat(t_adjGVec);
  arma::vec iSigmaGVec(m_n * (m_J-1), arma::fill::zeros);
  getPCGofSigmaAndVector(adjGVecLong, iSigmaGVec, t_excludechr);
  double VarP = as_scalar(adjGVecLong.t() * (iSigmaGVec - m_iSigmaX_XSigmaX * (m_CovaMat.t() * iSigmaGVec)));
  return(VarP);
}

arma::vec convert1(arma::mat xMat, // matrix: n x (J-1)
                   int n, int J) 
{
  arma::vec xVec(n*(J-1));
  int index = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      xVec(index) = xMat(i,j);
      index++;
    }
  }
  return(xVec);
}

arma::mat convert2(arma::vec xVec, // n(J-1) x 1 
                   int n, int J)
{
  arma::mat xMat(n,(J-1));
  int index = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      xMat(i,j) = xVec(index);
      index++;
    }
  }
  return(xMat);
}

double getVarWFast(arma::vec adjGVec,  // n x 1
                   arma::vec RPsiRVec) // n x 1
{
  int n = adjGVec.size();
  double VarW = 0;
  for(int i = 0; i < n; i++){
    VarW += RPsiRVec(i) * adjGVec(i) * adjGVec(i);
  }
  return(VarW);
}

// get a list for p value calculation in step 2
Rcpp::List getobjP(arma::mat t_Cova,     // matrix: n x p
                   arma::mat t_yMat,
                   arma::mat t_muMat,    // matrix: n x J
                   arma::mat t_iRMat)    // matrix: n x (J-1)
{
  int n = t_muMat.n_rows;
  int J = t_muMat.n_cols;
  int p = t_Cova.n_cols;
  
  // output for Step 2
  arma::mat XR_Psi_R(p, n*(J-1));                // p x n(J-1)
  arma::mat CovaMat = getCovaMat(t_Cova, J); // n(J-1) x p
  for(int k = 0; k < p; k++){
    arma::mat xMat = convert2(CovaMat.col(k), n, J);
    arma::vec temp = convert1(getPsixMat(xMat / t_iRMat, t_muMat) / t_iRMat, n, J);
    XR_Psi_R.row(k) = temp.t();
  }
  // arma::mat XXR_Psi_RX = CovaMat * inv(XR_Psi_R * CovaMat);             // (n(J-1) x p) * (p x p) = n(J-1) x p
  arma::mat XXR_Psi_RX_new = t_Cova * inv(XR_Psi_R * CovaMat);               // (n x p) * (p x p) = n x p
  
  // sum each (J-1) rows to 1 row: p x n(J-1) -> p x n
  arma::mat XR_Psi_R_new = sumCols(XR_Psi_R, J);      // p x n
  arma::mat ymuMat = t_yMat - t_muMat;                    // n x J
  arma::mat RymuMat = ymuMat.cols(0, J-2) / t_iRMat;    // n x (J-1): R %*% (y - mu)
  arma::mat RymuVec = sumCols(RymuMat, J);            // n x 1
  // arma::cube RPsiR = getRPsiR(muMat, iRMat, n, J, p); // (J-1) x (J-1) x n 
  arma::vec RPsiR = getRPsiR(t_muMat, t_iRMat, n, J, p); // (J-1) x (J-1) x n 
  
  Rcpp::List objP = Rcpp::List::create(Rcpp::Named("n")=n,
                                       Rcpp::Named("J")=J,
                                       Rcpp::Named("p")=p,
                                       Rcpp::Named("XXR_Psi_RX_new") = XXR_Psi_RX_new,
                                       Rcpp::Named("XR_Psi_R_new") = XR_Psi_R_new,           
                                       Rcpp::Named("RymuVec") = RymuVec,
                                       Rcpp::Named("RPsiR") = RPsiR,
                                       Rcpp::Named("muMat") = t_muMat,
                                       Rcpp::Named("iRMat") = t_iRMat);
  return(objP);
}

arma::vec getadjGFast(arma::vec GVec,
                      arma::mat XXR_Psi_RX_new,   // XXR_Psi_RX_new ( n x p )
                      arma::mat XR_Psi_R_new,     // XR_Psi_R_new ( p x n ), sum up XR_Psi_R ( p x n(J-1) ) for each subject 
                      int n, int p)
{
  // To increase computational efficiency when lots of GVec elements are 0
  arma::vec XR_Psi_RG1(p, arma::fill::zeros);
  for(int i = 0; i < n; i++){
    if(GVec(i) != 0){
      XR_Psi_RG1 += XR_Psi_R_new.col(i) * GVec(i);
    }
  }
  
  arma::vec adjGVec = GVec - XXR_Psi_RX_new * XR_Psi_RG1;
  return(adjGVec);
}

double getStatFast(arma::vec GVec,         // n x 1
                   arma::vec RymuVec)      // n x 1: row sum of the n x (J-1) matrix R %*% (yMat - muMat)
{
  int n = GVec.size();
  double Stat = 0;
  for(int i = 0; i < n; i++){
    if(GVec(i) != 0){
      Stat += GVec(i) * RymuVec(i);
    }
  }
  return(Stat);
}

Rcpp::List outputadjGFast(arma::vec GVec,
                          Rcpp::List objP)
{
  arma::vec adjGVec = getadjGFast(GVec, objP["XXR_Psi_RX_new"], objP["XR_Psi_R_new"], objP["n"], objP["p"]);
  double Stat = getStatFast(adjGVec, objP["RymuVec"]);
  double VarW = getVarWFast(adjGVec, objP["RPsiR"]);
  Rcpp::List outList = Rcpp::List::create(Rcpp::Named("adjGVec")=adjGVec,
                                          Rcpp::Named("Stat")=Stat,           
                                          Rcpp::Named("VarW")=VarW);
  
  return(outList);
}

// sum up each row: n1 x n2 matrix -> n1 x 1 vector
arma::vec getRowSums(arma::mat t_xMat)
{
  int n1 = t_xMat.n_rows;
  int n2 = t_xMat.n_cols;
  arma::vec y1Vec(n1, arma::fill::zeros);
  for(int i = 0; i < n1; i++){
    for(int j = 0; j < n2; j++){
      y1Vec(i) += t_xMat(i,j);
    }
  }
  return(y1Vec);
}

// get RPsiP: (J-1) x (J-1) x n 
// Only used in getVarWFast(): 
arma::vec getRPsiR(arma::mat t_muMat,
                   arma::mat t_iRMat,
                   int t_n, int t_J, int t_p)   
{
  // arma::cube RPsiR(J-1, J-1, n, arma::fill::zeros);
  arma::vec RPsiRVec(t_n, arma::fill::zeros);
  arma::mat muRMat = t_muMat.cols(0, t_J-2) / t_iRMat;
  for(int i = 0; i < t_n; i++){
    for(int j1 = 0; j1 < t_J-1; j1++){
      // RPsiR(j1,j1,i) += muRMat(i,j1) / iRMat(i,j1) - muRMat(i,j1) * muRMat(i,j1);
      RPsiRVec(i) += muRMat(i,j1) / t_iRMat(i,j1) - muRMat(i,j1) * muRMat(i,j1);
      for(int j2 = j1+1; j2 < t_J-1; j2++){
        // RPsiR(j1,j2,i) -= muRMat(i,j1) * muRMat(i,j2);
        RPsiRVec(i) -= 2* muRMat(i,j1) * muRMat(i,j2);
      }
    }
  }
  // return(RPsiR);
  return(RPsiRVec);
}

double calCV(arma::vec t_xVec){
  unsigned int n = t_xVec.size();
  double Mean = arma::mean(t_xVec);
  double Sd = arma::stddev(t_xVec);
  double CV = (Sd/Mean)/n;
  return CV;
}

// sum each (J-1) cols to 1 col: p x n(J-1) -> p x n (OR) p x (J-1) -> p x 1
arma::mat sumCols(arma::mat t_xMat,
                  int J)
{
  int n = t_xMat.n_cols / (J-1);
  int p = t_xMat.n_rows;
  arma::mat outMat(p, n, arma::fill::zeros);
  int index = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < J-1; j++){
      outMat.col(i) += t_xMat.col(index);
      index++;
    }
  }
  return(outMat);
}

// outMat = PsiMat %*% xMat, PsiMat is determined by muMat
arma::mat getPsixMat(arma::mat t_xMat,    // matrix: n x (J-1)
                     arma::mat t_muMat)   // matrix: n x J
{
  int n = t_muMat.n_rows;
  int J = t_muMat.n_cols;
  
  arma::mat Psi_xMat(n, J-1);
  // loop for samples
  for(int i = 0; i < n; i++){
    arma::rowvec muVec(J-1);
    for(int j = 0; j < J-1; j++){
      Psi_xMat(i,j) = t_muMat(i,j) * t_xMat(i,j);
      muVec(j) = t_muMat(i,j);
    }
    double sum_mu_x = sum(Psi_xMat.row(i));
    Psi_xMat.row(i) -= muVec * sum_mu_x; 
  }
  return(Psi_xMat);
}

// duplicate each row for (J-1) times: n x p -> n(J-1) x p
arma::mat getCovaMat(arma::mat t_Cova, unsigned int t_J)      
{
  unsigned int n = t_Cova.n_rows;
  unsigned int p = t_Cova.n_cols;
  
  arma::mat CovaMat(n * (t_J-1), p);
  int index = 0;
  for(unsigned int i = 0; i < n; i++){
    for(unsigned int j = 0; j < t_J-1; j++){
      CovaMat.row(index) = t_Cova.row(i);
      index++;
    }
  }
  return CovaMat;
}

Rcpp::List POLMMClass::getPOLMM()
{
  Rcpp::List outList = Rcpp::List::create(Rcpp::Named("N") = m_n,              // number of samples
                                          Rcpp::Named("M") = m_M,              // number of SNPs in Plink file
                                          Rcpp::Named("iter") = m_iter,
                                          Rcpp::Named("eta") = m_eta,          // X %*% beta + bVec
                                          Rcpp::Named("yVec") = m_yVec,        // matrix with dim of n x 1: observation
                                          Rcpp::Named("Cova") = m_Cova,        // matrix with dim of n(J-1) x p: covariates
                                          Rcpp::Named("muMat") = m_muMat,      // matrix with dim of n x J: probability
                                          Rcpp::Named("YMat") = m_YMat,        // matrix with dim of n x (J-1): working variables
                                          Rcpp::Named("beta") = m_beta,        // parameter for covariates
                                          Rcpp::Named("bVec") = m_bVec,        // terms of random effect 
                                          Rcpp::Named("tau") = m_tau,          // variance component
                                          Rcpp::Named("eps") = m_eps,          // cutpoints
                                          Rcpp::Named("LOCOList") = m_LOCOList);         
  
  return(outList);
}

}
// make a global variable for future usage
// static POLMM::POLMMClass* ptr_gPOLMMobj = NULL;

// ptr_gPOLMMobj = new POLMM::POLMMClass(t_muMat,
//                                       t_iRMat,
//                                       t_Cova,
//                                       t_yVec,
//                                       t_SparseGRM,
//                                       t_tau,
//                                       t_printPCGInfo,
//                                       t_tolPCG,
//                                       t_maxiterPCG);

// // [[Rcpp::export]]
// void setPOLMMobjInR(arma::mat t_muMat,
//                     arma::mat t_iRMat,
//                     arma::mat t_Cova,
//                     arma::vec t_yVec,
//                     Rcpp::List t_SPmatR,    // output of makeSPmatR()
//                     double t_tau,
//                     bool t_printPCGInfo,
//                     double t_tolPCG,
//                     int t_maxiterPCG)
// {
//   arma::umat locations = t_SPmatR["locations"];
//   arma::vec values = t_SPmatR["values"];
//   arma::sp_mat SparseGRM = arma::sp_mat(locations, values);
//   ptr_gPOLMMobj = new POLMM::POLMMClass(t_muMat,
//                                         t_iRMat,
//                                         t_Cova,
//                                         t_yVec,
//                                         SparseGRM,
//                                         t_tau,
//                                         t_printPCGInfo,
//                                         t_tolPCG,
//                                         t_maxiterPCG);
// }

// get RPsiP: (J-1) x (J-1) x n 
// Only used in getVarWFast(): 

